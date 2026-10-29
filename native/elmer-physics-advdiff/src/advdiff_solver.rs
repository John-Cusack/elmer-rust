//! `AdvectionDiffusionSolver`: the `Procedure = "AdvectionDiffusion"
//! "AdvectionDiffusionSolver"` solver, driving the element kernels through
//! [`SolverCtx`] the way AdvectionDiffusion.F90:67 does.

use elmer_core::{
    normalize_key, ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF,
};
use elmer_fe::{stab_param, ElmerQuadrature, LagrangeBasis};
use elmer_solver::backend::TimeTerms;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::compose::{advdiff_boundary_local, advdiff_bulk_local, AdvDiffConfig, AdvDiffNodal, Stabilization};
use crate::error::AdvDiffError;
use crate::fe::{condensate, BasisHessian, Fe};
use crate::preflight::{
    check_advdiff_supported, concentration_units, AdvDiffSections, ADVDIFF_PROCEDURE,
};

impl From<AdvDiffError> for SolverError {
    fn from(e: AdvDiffError) -> Self {
        match e {
            AdvDiffError::Unsupported { .. } => SolverError::Unsupported(e.to_string()),
            e => SolverError::fatal("AdvectionDiffusion", e.to_string()),
        }
    }
}

/// Registers `("AdvectionDiffusion", "AdvectionDiffusionSolver")`.
pub fn register_advdiff(registry: &mut SolverRegistry) {
    registry.register(ADVDIFF_PROCEDURE.0, ADVDIFF_PROCEDURE.1, advdiff_factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`AdvDiffSolver`].
pub fn advdiff_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(AdvDiffSolver::default()))
}

/// Second derivatives of the linear simplex bases (202, 303, 504), which are
/// identically zero in physical space. The other element types need the
/// `GlobalSecondDerivatives` of `elmer-fe`, which it does not provide yet.
struct SimplexHessian;

impl BasisHessian for SimplexHessian {
    fn eval_hessian(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        _uvw: [f64; 3],
        out: &mut Vec<[[f64; 3]; 3]>,
    ) -> Result<(), String> {
        if !matches!(etype.0, 202 | 303 | 504) {
            return Err(format!(
                "SUPG: second derivatives of element type {} are not available",
                etype.0
            ));
        }
        out.clear();
        out.resize(xyz.len(), [[0.0; 3]; 3]);
        Ok(())
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum Method {
    Galerkin,
    Bubbles,
    Stabilized,
}

/// `GetString(Solver % Values, 'Stabilization Method')` or, without it,
/// the `Stabilize` and `Bubbles` flags (bubbles being the default).
fn stabilization_method(solver: &ParamMap) -> Result<Method, SolverError> {
    if let Some(s) = solver.string("Stabilization Method") {
        return match s.to_ascii_lowercase().as_str() {
            "stabilized" => Ok(Method::Stabilized),
            "bubbles" => Ok(Method::Bubbles),
            other => Err(SolverError::fatal(
                "AdvectionDiffusion",
                format!("Unknown stabilization method: {other}"),
            )),
        };
    }
    match (solver.logical("Stabilize"), solver.logical("Bubbles")) {
        (None, None) => Ok(Method::Bubbles),
        (Some(true), Some(true)) => Err(SolverError::fatal(
            "AdvectionDiffusion",
            "Choose either Bubbles or Stabilize!",
        )),
        (Some(true), _) => Ok(Method::Stabilized),
        (_, Some(true)) => Ok(Method::Bubbles),
        _ => Ok(Method::Galerkin),
    }
}

fn is_cartesian(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(()),
        Some(s) => {
            let s = normalize_key(s);
            if s.starts_with("cartesian") {
                Ok(())
            } else {
                Err(SolverError::Unsupported(format!(
                    "Coordinate System = {s}: DiffuseConvectiveGenCompose is not ported"
                )))
            }
        }
    }
}

pub struct AdvDiffSolver {
    method: Option<Method>,
    local: LocalSystem,
    tmp: LocalSystem,
}

impl Default for AdvDiffSolver {
    fn default() -> Self {
        Self {
            method: None,
            local: LocalSystem::zeros(0),
            tmp: LocalSystem::zeros(0),
        }
    }
}

fn real_at_nodes(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    key: &str,
    e: ElemId,
) -> Result<Option<Vec<f64>>, SolverError> {
    match p {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
}

/// A scalar variable's values at `nodes`; `None` if there is no such
/// variable, an `Err` where it has no dof at a node.
fn field_at_nodes(
    ctx: &SolverCtx<'_>,
    name: &str,
    nodes: &[u32],
) -> Result<Option<Vec<f64>>, SolverError> {
    let Some(v) = ctx.sim.vars.get(name) else {
        return Ok(None);
    };
    nodes
        .iter()
        .map(|&n| {
            v.slot(n).map(|s| v.values[s * v.dofs]).ok_or_else(|| {
                SolverError::fatal(
                    "AdvectionDiffusion",
                    format!("variable {name} has no value at node {n}"),
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// `<name> Diffusivity` by shape (`ListGetRealArray`): 1x1 isotropic, n x 1
/// diagonal, n x n full. Scalars may depend on fields; tensors are constant.
fn diffusivity(
    ctx: &SolverCtx<'_>,
    mat: &ParamMap,
    key: &str,
    e: ElemId,
    n: usize,
) -> Result<Vec<[[f64; 3]; 3]>, SolverError> {
    match mat.get(key) {
        None => Err(SolverError::fatal(
            "AdvectionDiffusion",
            format!("Failed to find material keyword \"{key}\""),
        )),
        Some(Value::Real { rows, cols, data }) if rows * cols > 1 => {
            let (r, c) = (*rows, *cols);
            let mut k = [[0.0; 3]; 3];
            if c == 1 {
                for i in 0..r.min(3) {
                    k[i][i] = data[i];
                }
            } else {
                for i in 0..r.min(3) {
                    for j in 0..c.min(3) {
                        k[i][j] = data[i * c + j];
                    }
                }
            }
            Ok(vec![k; n])
        }
        Some(Value::Dependent(d)) if d.rows * d.cols > 1 => Err(SolverError::Unsupported(
            format!("\"{key}\": a non-constant tensor"),
        )),
        Some(_) => {
            let d = ctx
                .real_at_nodes(mat, key, e)?
                .ok_or_else(|| SolverError::fatal("AdvectionDiffusion", format!("{key} missing")))?;
            Ok(d.iter()
                .map(|&x| {
                    let mut k = [[0.0; 3]; 3];
                    for (i, row) in k.iter_mut().enumerate() {
                        row[i] = x;
                    }
                    k
                })
                .collect())
        }
    }
}

impl Solver for AdvDiffSolver {
    /// `AdvectionDiffusionSolver_init` plus the refusals: Time Derivative
    /// Order 1, then every keyword the kernels would otherwise ignore.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|f| f.dofs != 1) {
            return Err(SolverError::Unsupported(
                "AdvectionDiffusion solves a scalar variable (Variable DOFs = 1)".into(),
            ));
        }
        is_cartesian(ctx)?;
        self.method = Some(stabilization_method(ctx.params())?);
        ctx.set_time_derivative_order(1)?;
        let comp = normalize_key(ctx.variable_name().unwrap_or_default());

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut bodies: Vec<usize> = Vec::new();
        for e in elems {
            if let Some(b) = ctx.body_params(e) {
                if !bodies.contains(&b.body_id) {
                    bodies.push(b.body_id);
                }
            }
        }
        let model = &ctx.sim.model;
        let bcs: Vec<&dyn Params> = model
            .boundary_conditions
            .iter()
            .map(|p| p as &dyn Params)
            .collect();
        for b in bodies {
            let Some(bp) = model.body_params(b) else {
                continue;
            };
            let sections = AdvDiffSections {
                solver: ctx.params(),
                equation: bp.equation.map(|p| p as &dyn Params),
                material: bp.material.map(|p| p as &dyn Params),
                boundaries: bcs.clone(),
            };
            check_advdiff_supported(&sections, &comp)?;
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        const CALLER: &str = "AdvectionDiffusion";
        let method = match self.method {
            Some(m) => m,
            None => stabilization_method(ctx.params())?,
        };
        let dim = ctx.mesh.dim();
        let basis = LagrangeBasis::new(dim);
        let quad = ElmerQuadrature;
        let hessian = SimplexHessian;
        let fe = Fe {
            basis: &basis,
            quad: &quad,
            hessian: Some(&hessian),
        };
        let comp = normalize_key(ctx.variable_name().unwrap_or_default());
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let constants = ctx.constants().clone();

        ctx.info(CALLER, "-------------------------------------", 4);
        ctx.info(CALLER, "Solving for species transport", 4);
        ctx.info(CALLER, "-------------------------------------", 4);

        let mut dofs = Vec::new();
        let mut error_written = false;
        for iter in 1..=maxiter {
            ctx.info(CALLER, &format!("Iteration step: {iter}"), 4);
            ctx.default_initialize()?;

            // `AbsoluteMass` keeps the value of the last body visited, which
            // the boundary loop reads.
            let mut absolute_mass = false;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let nodes: Vec<u32> = el.nodes.to_vec();
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let (nodal, abs, convective) =
                    self.element_nodal(ctx, e, &comp, &nodes, n, &constants)?;
                absolute_mass = abs;

                let stabilization = match method {
                    Method::Galerkin => Stabilization::None,
                    Method::Bubbles => Stabilization::Bubbles,
                    Method::Stabilized => {
                        let (m_k, h_k) = stab_param(etype, &xyz, dim, false)
                            .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
                        Stabilization::Supg { h_k, m_k }
                    }
                };
                let cfg = AdvDiffConfig { dim, stabilization };
                advdiff_bulk_local(&fe, etype, &xyz, &nodal, &cfg, &mut self.local)?;

                if transient {
                    if self.local.n > n {
                        bubble_first_order_time(ctx, &mut self.local, e, dt)?;
                    } else {
                        ctx.default_first_order_time(&mut self.local, e, dt)?;
                    }
                }

                if self.local.n > n && method == Method::Bubbles && convective {
                    let nb = self.local.n;
                    condensate(n, n, &mut self.local.stiff, &mut self.local.force, None)?;
                    self.tmp = LocalSystem::zeros(n);
                    for p in 0..n {
                        for q in 0..n {
                            self.tmp.stiff[p * n + q] = self.local.stiff[p * nb + q];
                        }
                        self.tmp.force[p] = self.local.force[p];
                    }
                    ctx.default_update(&self.tmp, e)?;
                } else {
                    ctx.default_update(&self.local, e)?;
                }
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let bc = ctx.bc_params(e);
                let mut alpha = vec![0.0; n];
                let mut s_ext = vec![0.0; n];
                if let Some(a) = real_at_nodes(ctx, bc, "Mass Transfer Coefficient", e)? {
                    alpha = a;
                    s_ext = real_at_nodes(ctx, bc, "External Concentration", e)?
                        .unwrap_or_else(|| vec![0.0; n]);
                    if !absolute_mass {
                        if !error_written {
                            ctx.warn(
                                CALLER,
                                "Mass transfer coefficient possible to use only with \
                                 absolute mass concentrations; ignoring mass transfer BC",
                            );
                            error_written = true;
                        }
                        alpha = vec![0.0; n];
                        s_ext = vec![0.0; n];
                    }
                }
                let flux = real_at_nodes(ctx, bc, &format!("{comp} Flux"), e)?;
                let load: Vec<f64> = (0..n)
                    .map(|j| alpha[j] * s_ext[j] + flux.as_ref().map_or(0.0, |f| f[j]))
                    .collect();
                advdiff_boundary_local(&fe, etype, &xyz, &alpha, &load, &mut self.local)?;
                if transient {
                    self.local.mass = Some(vec![0.0; n * n]);
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}

/// `Default1stOrderTime` on an element with bubble rows. Elmer's local
/// time schemes loop over the whole allocated local matrix, so the mass
/// matrix enters the bubble rows and columns too, while the previous
/// solution has only the `n` nodal entries (the bubbles' previous values
/// read as zero).
fn bubble_first_order_time(
    ctx: &SolverCtx<'_>,
    local: &mut LocalSystem,
    e: ElemId,
    dt: f64,
) -> Result<(), SolverError> {
    const CALLER: &str = "Default1stOrderTime";
    let nb = local.n;
    let mut idx = Vec::new();
    ctx.element_dofs(e, &mut idx);
    let name = ctx
        .variable_name()
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let v = ctx
        .sim
        .vars
        .get(name)
        .ok_or_else(|| SolverError::NoVariable(name.into()))?;
    let pad = |h: &[f64]| -> Vec<f64> {
        let mut out: Vec<f64> = idx
            .iter()
            .map(|&i| if i == NO_DOF { 0.0 } else { h[i as usize] })
            .collect();
        out.resize(nb, 0.0);
        out
    };
    let prev: Vec<Vec<f64>> = v.prev_values.iter().map(|h| pad(h)).collect();
    let cur = pad(&v.values);
    let slot = ctx
        .sim
        .slots
        .iter()
        .find(|s| s.id == ctx.solver_id())
        .ok_or_else(|| SolverError::fatal(CALLER, "solver slot not found"))?;
    let backend = ctx
        .sim
        .backends
        .time
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("time backend (elmer-timeint)".into()))?;
    let terms = TimeTerms {
        dt,
        params: &slot.params,
        simulation: &ctx.sim.model.simulation,
        prev: &prev,
        cur: &cur,
        done_time: slot.done_time,
        order: slot.order,
        dts: &ctx.sim.prev_sizes,
    };
    backend.add_time_terms(local, &terms)
}

impl AdvDiffSolver {
    /// The nodal coefficients of element `e`, with its body's
    /// `AbsoluteMass` and whether the convection flag is on, following the
    /// per-body setup of AdvectionDiffusion.F90:338-736.
    fn element_nodal(
        &self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        comp: &str,
        nodes: &[u32],
        n: usize,
        constants: &ParamMap,
    ) -> Result<(AdvDiffNodal, bool, bool), SolverError> {
        const CALLER: &str = "AdvectionDiffusion";
        let bp = ctx.body_params(e).ok_or_else(|| {
            SolverError::fatal(CALLER, format!("element {} has no body", e.0))
        })?;
        let empty = ParamMap::new();
        let eq = bp.equation.unwrap_or(&empty);
        let mat = bp.material.unwrap_or(&empty);

        let flag = eq
            .string(&format!("{comp} Convection"))
            .or_else(|| eq.string("Convection"))
            .map(str::to_ascii_lowercase)
            .unwrap_or_default();
        let absolute = matches!(concentration_units(eq).as_deref(), Some("absolute mass"));
        let convective = flag == "constant" || flag == "computed";

        // Temperature: the field named by the material, else `Temperature`
        // if it exists, else `Reference Temperature`.
        let temp_name = mat
            .string("Temperature Field Variable")
            .map(str::to_string)
            .unwrap_or_else(|| "temperature".into());
        let temperature = match field_at_nodes(ctx, &temp_name, nodes)? {
            Some(t) => t,
            None if mat.string("Temperature Field Variable").is_some() => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!("No temperature variable {temp_name} available"),
                ))
            }
            None => real_at_nodes(ctx, Some(mat), "Reference Temperature", e)?
                .unwrap_or_else(|| vec![0.0; n]),
        };
        let phi_name = mat
            .string(&format!("{comp} EPotential Field Variable"))
            .or_else(|| mat.string("EPotential Field Variable"))
            .map(str::to_string);
        let potential = match &phi_name {
            Some(name) => field_at_nodes(ctx, name, nodes)?.ok_or_else(|| {
                SolverError::fatal(
                    CALLER,
                    format!("No electrical potential variable {name} available"),
                )
            })?,
            None => match field_at_nodes(ctx, "potential", nodes)? {
                Some(p) => p,
                None => real_at_nodes(ctx, Some(mat), "Reference EPotential", e)?
                    .unwrap_or_else(|| vec![0.0; n]),
            },
        };

        let diff = diffusivity(ctx, mat, &format!("{comp} Diffusivity"), e, n)?;
        let soret = real_at_nodes(ctx, Some(mat), &format!("{comp} Soret Diffusivity"), e)?;
        let ion_charge = ctx.creal(mat, &format!("{comp} Ion Charge"))?;
        let mut ne_const = match ion_charge {
            Some(q) => {
                let unit = ctx.creal(constants, "Unit Charge")?.ok_or_else(|| {
                    SolverError::fatal(CALLER, "Constants: Unit Charge not found")
                })?;
                let kb = ctx.creal(constants, "Boltzmann Constant")?.ok_or_else(|| {
                    SolverError::fatal(CALLER, "Constants: Boltzmann Constant not found")
                })?;
                temperature.iter().map(|t| q / t * unit / kb).collect()
            }
            None => Vec::new(),
        };
        let migration = ion_charge.is_some_and(|q| q != 0.0);

        let density = real_at_nodes(ctx, Some(mat), "Density", e)?.ok_or_else(|| {
            SolverError::fatal(CALLER, "Failed to find material keyword \"Density\"")
        })?;

        let mut velocity = vec![[0.0; 3]; n];
        if flag == "constant" {
            for i in 0..3 {
                let key = format!("Convection Velocity {}", i + 1);
                let v = match real_at_nodes(ctx, Some(eq), &format!("{comp} {key}"), e)? {
                    Some(v) => Some(v),
                    None => real_at_nodes(ctx, Some(mat), &key, e)?,
                };
                if let Some(v) = v {
                    for (vel, x) in velocity.iter_mut().zip(v) {
                        vel[i] = x;
                    }
                }
            }
        } else if flag == "computed" {
            let name = eq
                .string(&format!("{comp} Convection Field Variable"))
                .or_else(|| mat.string("Convection Field Variable"))
                .unwrap_or("flow solution");
            let v = ctx.sim.vars.get(name).ok_or_else(|| {
                SolverError::fatal(CALLER, "Give > Convection Field Variable <")
            })?;
            if !matches!(v.dofs, 3 | 4) {
                return Err(SolverError::Unsupported(format!(
                    "Convection Field Variable {name} with {} dofs",
                    v.dofs
                )));
            }
            for (vel, &node) in velocity.iter_mut().zip(nodes) {
                if let Some(s) = v.slot(node) {
                    let b = s * v.dofs;
                    vel[0] = v.values[b];
                    vel[1] = v.values[b + 1];
                    if v.dofs == 4 {
                        vel[2] = v.values[b + 2];
                    }
                }
            }
        }

        let mut mesh_velocity = Vec::new();
        if let Some(mv) = ctx.sim.vars.get("mesh velocity") {
            mesh_velocity = nodes
                .iter()
                .map(|&node| {
                    let mut m = [0.0; 3];
                    if let Some(s) = mv.slot(node) {
                        for (k, x) in m.iter_mut().enumerate().take(mv.dofs.min(3)) {
                            *x = mv.values[s * mv.dofs + k];
                        }
                    }
                    m
                })
                .collect();
        }

        let c1: Vec<f64> = if convective || migration {
            if absolute {
                vec![1.0; n]
            } else {
                density.clone()
            }
        } else {
            vec![0.0; n]
        };
        let (ct, c2): (Vec<f64>, Vec<[[f64; 3]; 3]>) = if absolute {
            (vec![1.0; n], diff)
        } else {
            for (ne, &d) in ne_const.iter_mut().zip(&density) {
                *ne = if d != 0.0 { *ne / d } else { 0.0 };
            }
            let c2 = diff
                .iter()
                .zip(&density)
                .map(|(k, &d)| {
                    let mut r = *k;
                    for row in r.iter_mut() {
                        for x in row.iter_mut() {
                            *x *= d;
                        }
                    }
                    r
                })
                .collect();
            (density.clone(), c2)
        };
        let c0 = if absolute { vec![1.0; n] } else { density };

        let mut load = vec![0.0; n];
        if let Some(bf) = bp.body_force {
            if let Some(l) = real_at_nodes(ctx, Some(bf), &format!("{comp} Diffusion Source"), e)? {
                load = l;
            }
        }

        let nodal = AdvDiffNodal {
            ct,
            c1,
            c2,
            velocity,
            mesh_velocity,
            load,
            c0,
            soret_diffusivity: soret.unwrap_or_default(),
            temperature,
            ne_const,
            potential,
        };
        Ok((nodal, absolute, convective))
    }
}
