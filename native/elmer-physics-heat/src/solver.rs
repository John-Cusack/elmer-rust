//! `HeatSolver`: the `Procedure = "HeatSolve" "HeatSolver"` solver, driving the
//! element kernels through [`SolverCtx`] the way `HeatSolver` in
//! HeatSolve.F90:254 drives `LocalMatrixVec` and `LocalMatrixBC`.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::material::{
    Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, Stabilization,
};
use crate::preflight::stabilized;
use crate::radiation::{is_diffuse_gray, temps4, DiffuseGray};
use crate::radiators::{is_radiator_bc, Radiators};
use crate::{
    check_supported, heat_boundary_local, heat_bulk_local, HeatError, HeatMaterialAtNodes,
    HeatSections, PROCEDURE,
};

/// Builds the element basis and quadrature the solver integrates with. The
/// registry's factories are plain `fn` pointers, so the front end (`elmer-fe`)
/// hands its implementations over through [`install_fe`] instead.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs the basis and quadrature provider used by [`factory`]. Returns
/// `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Registers `("HeatSolve", "HeatSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`HeatSolver`]. Fails
/// with `NotWired` until [`install_fe`] has been called.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let provider = FE
        .get()
        .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))?;
    let (basis, quad) = provider();
    Ok(Box::new(HeatSolver::new(basis, quad)))
}

impl From<HeatError> for SolverError {
    fn from(e: HeatError) -> Self {
        match e {
            HeatError::Unsupported { .. } => SolverError::Unsupported(e.to_string()),
            e => SolverError::fatal("HeatSolver", e.to_string()),
        }
    }
}

pub struct HeatSolver {
    pub(crate) basis: Box<dyn ElementBasis + Send>,
    pub(crate) quad: Box<dyn Quadrature + Send>,
    pub(crate) local: LocalSystem,
    /// The diffuse gray radiation, found and factored on the first solve.
    diffuse: Option<Option<DiffuseGray>>,
    /// The point radiators, likewise.
    radiators: Option<Option<Radiators>>,
}

impl HeatSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: LocalSystem::zeros(0),
            diffuse: None,
            radiators: None,
        }
    }
}

/// `ListGetElementReal` of a scalar keyword at the nodes of element `e`:
/// constants, tables, MATC and procedures through `elmer-expr`, evaluated at
/// the current iterate of the fields they depend on.
pub(crate) fn real_at_nodes(
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

/// `Heat Conductivity` by shape: 1x1 scalar, n x 1 diagonal, n x n full.
/// A scalar may depend on fields; a tensor must be constant.
pub(crate) fn conductivity(
    ctx: &SolverCtx<'_>,
    p: &ParamMap,
    e: ElemId,
    n: usize,
) -> Result<Option<Conductivity>, SolverError> {
    let key = "Heat Conductivity";
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) if rows * cols > 1 => {
            let (r, c) = (*rows, *cols);
            if (2..=3).contains(&r) && c == 1 {
                let mut d = [0.0; 3];
                d[..r].copy_from_slice(&data[..r]);
                Ok(Some(Conductivity::Diagonal(vec![d; n])))
            } else if (2..=3).contains(&r) && r == c {
                let mut k = [[0.0; 3]; 3];
                for i in 0..r {
                    for j in 0..r {
                        k[i][j] = data[i * c + j];
                    }
                }
                Ok(Some(Conductivity::Full(vec![k; n])))
            } else {
                Err(SolverError::Unsupported(format!(
                    "\"{key}\" with shape {r}x{c}"
                )))
            }
        }
        Some(Value::Dependent(d)) if d.rows * d.cols > 1 => Err(SolverError::Unsupported(format!(
            "\"{key}\": a non-constant tensor"
        ))),
        Some(_) => Ok(ctx.real_at_nodes(p, key, e)?.map(Conductivity::Scalar)),
    }
}

/// `Convection = Constant`: `Convection Velocity i` of the material at the
/// element nodes; a component without its own keyword is zero.
pub(crate) fn convection_velocity(
    ctx: &SolverCtx<'_>,
    m: Option<&ParamMap>,
    e: ElemId,
    n: usize,
) -> Result<Vec<[f64; 3]>, SolverError> {
    let mut v = vec![[0.0; 3]; n];
    for i in 0..3 {
        if let Some(c) = real_at_nodes(ctx, m, &format!("Convection Velocity {}", i + 1), e)? {
            for (vn, ci) in v.iter_mut().zip(c) {
                vn[i] = ci;
            }
        }
    }
    Ok(v)
}

/// The kernel inputs of bulk element `e` with `n` nodes: material and body
/// force evaluated at the element nodes. `legacy` forms the heat source as
/// HeatSolveLegacy does, the nodal product `Density * Heat Source`, instead
/// of `HeatSolveVec`'s product of interpolants.
pub(crate) fn bulk_inputs(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    n: usize,
    legacy: bool,
) -> Result<HeatMaterialAtNodes, SolverError> {
    const CALLER: &str = "HeatSolver";
    let bp = ctx
        .body_params(e)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
    let (m, bf) = (bp.material, bp.body_force);
    // Elmer's `init` accepts these; only a solver that runs needs the flow
    // solution they read, so the refusal waits for the solve (a solver with
    // `Exec Solver = Never` never gets here).
    if bp
        .equation
        .and_then(|q| q.string("Convection"))
        .is_some_and(|c| c.eq_ignore_ascii_case("computed"))
    {
        return Err(SolverError::Unsupported(
            "Convection = Computed: reading the flow solution as the convection field is \
             not implemented"
                .into(),
        ));
    }
    if m.is_some_and(|m| m.get("Compressibility Model").is_some()) {
        return Err(SolverError::Unsupported(
            "material keyword \"Compressibility Model\": compressible heat is not \
             implemented"
                .into(),
        ));
    }
    let convected = bp
        .equation
        .and_then(|q| q.string("Convection"))
        .is_some_and(|c| c.eq_ignore_ascii_case("constant"));
    let density = real_at_nodes(ctx, m, "Density", e)?;
    let mut heat_source = real_at_nodes(ctx, bf, "Heat Source", e)?;
    let mut volumetric_heat_source = real_at_nodes(ctx, bf, "Volumetric Heat Source", e)?;
    if legacy && volumetric_heat_source.is_none() {
        if let Some(src) = heat_source.take() {
            let rho = density.as_deref();
            volumetric_heat_source = Some(
                src.iter()
                    .enumerate()
                    .map(|(i, s)| s * rho.map_or(0.0, |r| r[i]))
                    .collect(),
            );
        }
    }
    Ok(HeatMaterialAtNodes {
        conductivity: match m {
            Some(m) => conductivity(ctx, m, e, n)?,
            None => None,
        },
        density,
        heat_capacity: real_at_nodes(ctx, m, "Heat Capacity", e)?,
        heat_transfer_multiplier: real_at_nodes(ctx, m, "Heat Transfer Multiplier", e)?,
        heat_source,
        volumetric_heat_source,
        convection_velocity: if convected {
            Some(convection_velocity(ctx, m, e, n)?)
        } else {
            None
        },
    })
}

pub(crate) fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) => {
            let s = elmer_core::normalize_key(s);
            if s.starts_with("cartesian") {
                Ok(false)
            } else if s == "axi symmetric" || s == "cylindric symmetric" {
                Ok(true)
            } else {
                Err(SolverError::Unsupported(format!("Coordinate System = {s}")))
            }
        }
    }
}

impl Solver for HeatSolver {
    /// Refuses every keyword in the sections this solver reads that the
    /// kernels would otherwise ignore (see [`check_supported`]), and
    /// `Convection = Computed`, which needs the flow solution field.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        axisymmetric(ctx)?;
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
            let sections = HeatSections {
                solver: ctx.params(),
                equation: bp.equation.map(|p| p as &dyn Params),
                material: bp.material.map(|p| p as &dyn Params),
                body_force: bp.body_force.map(|p| p as &dyn Params),
                boundaries: bcs.clone(),
            };
            check_supported(&sections)?;
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        const CALLER: &str = "HeatSolver";
        if ctx.params().logical("Legacy Assembly") == Some(true) {
            self.solve_legacy(ctx, transient)?;
            // `DefaultFinish` (HeatSolve.F90): `Apply Explicit Control`.
            return ctx.default_finish();
        }
        let axisym = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();

        ctx.info(
            CALLER,
            "------------------------------------------------",
            4,
        );
        ctx.info(CALLER, "Solving energy equation for temperature", 4);

        let mut bulk_cfg = HeatBulkConfig::new(dim);
        bulk_cfg.axisymmetric = axisym;
        bulk_cfg.transient = transient;
        let bnd_cfg = HeatBoundaryConfig {
            axisymmetric: axisym,
            quad_points: None,
        };
        let mut dofs = Vec::new();
        // SUPG (`Stabilize`) takes `Element % hK` and `Element %
        // StabilizationMK`, which `MeshStabParams` fills at mesh load with
        // `StabParam` in the model dimension.
        let supg = stabilized(ctx.params());
        let long_edge = ctx
            .simulation_params()
            .logical("Stabilization Use Longest Element Edge")
            .unwrap_or(false);

        if self.diffuse.is_none() {
            self.diffuse = Some(DiffuseGray::setup(ctx, &bnd, axisym)?);
        }
        if self.radiators.is_none() {
            self.radiators = Some(Radiators::setup(ctx, &bnd)?);
        }
        let stefan = ctx.constants().real("Stefan Boltzmann");
        if self.diffuse.as_ref().is_some_and(Option::is_some) && stefan.is_none() {
            return Err(SolverError::fatal(
                CALLER,
                "Stefan Boltzmann is not defined".to_string(),
            ));
        }

        for iter in 1..=maxiter {
            ctx.info(CALLER, &format!("Heat solver iteration: {iter}"), 4);
            ctx.default_initialize()?;
            if let Some(Some(dg)) = self.diffuse.as_mut() {
                dg.extend_topology(ctx)?;
            }
            let newton = {
                let id = ctx.solver_id();
                ctx.sim
                    .slots
                    .iter()
                    .find(|s| s.id == id)
                    .is_some_and(|s| s.state.newton_active)
            };
            let t4 = match self.diffuse {
                Some(Some(_)) => temps4(ctx, &bnd),
                _ => Default::default(),
            };

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = bulk_inputs(ctx, e, n, false)?;
                bulk_cfg.stabilization = if supg && mat.convection_velocity.is_some() {
                    let (m_k, h_k) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                        .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
                    Stabilization::Supg { h_k, m_k }
                } else {
                    Stabilization::None
                };
                heat_bulk_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &bulk_cfg,
                    &mut self.local,
                )?;
                if transient {
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_finish_bulk_assembly()?;

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let mut bc = {
                    let p = ctx.bc_params(e);
                    HeatBcAtNodes {
                        heat_flux: real_at_nodes(ctx, p, "Heat Flux", e)?,
                        heat_transfer_coefficient: real_at_nodes(
                            ctx,
                            p,
                            "Heat Transfer Coefficient",
                            e,
                        )?,
                        external_temperature: real_at_nodes(ctx, p, "External Temperature", e)?,
                    }
                };
                // Radiator BC: the radiators' power adds to the flux, unless
                // the element is a diffuse gray one.
                if let Some(Some(r)) = &self.radiators {
                    if is_radiator_bc(ctx, e) && !is_diffuse_gray(ctx, e) {
                        if let Some(q) = r.flux(e) {
                            let n = ctx.mesh.element(e).nodes.len();
                            let flux = bc.heat_flux.get_or_insert_with(|| vec![0.0; n]);
                            flux.iter_mut().for_each(|x| *x += q);
                        }
                    }
                }
                if bc.heat_flux.is_some() || bc.heat_transfer_coefficient.is_some() {
                    let el = ctx.mesh.element(e);
                    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                    heat_boundary_local(
                        &*self.basis,
                        &*self.quad,
                        el.etype,
                        &xyz,
                        &bc,
                        &bnd_cfg,
                        &mut self.local,
                    )?;
                    ctx.default_update(&self.local, e)?;
                }
                if let Some(Some(dg)) = &self.diffuse {
                    if is_diffuse_gray(ctx, e) {
                        let sigma = stefan.expect("checked above");
                        dg.assemble(ctx, &*self.basis, &*self.quad, e, &t4, newton, sigma)?;
                    }
                }
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        // `DefaultFinish` (HeatSolve.F90): `Apply Explicit Control`.
        ctx.default_finish()
    }
}
