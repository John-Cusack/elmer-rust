//! `StressSolve` / `StressSolver` as an `elmer_solver::Solver`.
//!
//! One steady linear pass: assemble the bulk and boundary systems with the
//! kernels of this crate, apply Dirichlet conditions, solve. Everything this
//! crate does not implement is refused up front instead of being dropped.

use std::sync::OnceLock;

use elmer_core::{
    ElemId, ElementBasis, ElementType, FieldKind, FieldRef, LocalSystem, MeshView, ParamMap, Params, Quadrature, Value,
};
use elmer_solver::sim::VarSpec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

use crate::boundary::{stress_boundary_local, stress_boundary_local_p, Spring, StressBoundaryInput};
use crate::bulk::{stress_bulk_local_p, stress_bulk_local, CoordSystem, IpProps, StressBulkInput};
use crate::material::{diagonal_c, ElasticLaw, Tensor3};
use crate::post::{self, stress_strain_local, StressPostInput};
use crate::props;

const CALLER: &str = "StressSolver";

/// Finite-element services, until `elmer-fe` provides them to every solver.
pub struct FeServices {
    pub basis: Box<dyn ElementBasis + Send + Sync>,
    pub quad: Box<dyn Quadrature + Send + Sync>,
}

static FE: OnceLock<FeServices> = OnceLock::new();

/// Install the process-wide basis and quadrature used by registry-built
/// `StressSolver`s. Returns the services back if already installed.
pub fn install_fe(fe: FeServices) -> Result<(), FeServices> {
    FE.set(fe)
}

/// Register `("StressSolve", "StressSolver")`.
pub fn register(reg: &mut SolverRegistry) {
    let (lib, entry) = crate::REGISTRY_ENTRY;
    reg.register(lib, entry, factory);
}

fn factory(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if FE.get().is_none() {
        return Err(SolverError::NotWired("element basis and quadrature (elmer-fe)".into()));
    }
    Ok(Box::new(StressSolver { fe: None }))
}

/// The linear elasticity solver.
pub struct StressSolver {
    fe: Option<FeServices>,
}

impl StressSolver {
    /// A solver that owns its basis and quadrature.
    pub fn new(fe: FeServices) -> Self {
        Self { fe: Some(fe) }
    }

    fn fe(&self) -> &FeServices {
        self.fe.as_ref().or_else(|| FE.get()).expect("checked by factory")
    }
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn flag(p: Option<&ParamMap>, key: &str) -> bool {
    p.and_then(|p| p.logical(key)).unwrap_or(false)
}

fn present_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

/// Model-wide refusals, from the solver section.
/// The basis, rule and dof count of element `e` from the solver's p dofs
/// (`ElementInfo`, `GaussPoints`, `GetElementNOFDOFs`).
struct PElement {
    basis: Box<dyn ElementBasis + Send + Sync>,
    rule: Vec<elmer_core::QuadPoint>,
    nd: usize,
}

fn p_element(d: &dyn elmer_solver::SolverDofs, e: ElemId) -> Result<PElement, SolverError> {
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    Ok(PElement {
        basis: d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
        rule: d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?,
        nd: slots.len(),
    })
}

/// `values` (at `mesh_nodes`) at `nodes` instead.
fn reorder(values: Vec<f64>, mesh_nodes: &[u32], nodes: &[u32]) -> Vec<f64> {
    if mesh_nodes == nodes {
        return values;
    }
    nodes
        .iter()
        .map(|n| values[mesh_nodes.iter().position(|m| m == n).expect("same node set")])
        .collect()
}

fn check_solver_params(p: &ParamMap) -> Result<(), SolverError> {
    for (key, what) in [
        ("large deflection", "large deflection"),
        ("incompressible", "incompressible formulation"),
        ("stability analysis", "stability analysis"),
        ("eigen analysis", "eigen analysis"),
        ("calculate pangle", "principal angle postprocessing"),
        ("constant bulk system", "constant bulk system"),
        ("model lumping", "model lumping"),
        ("update reference displacement", "staged construction"),
        ("legacy assembly", "legacy assembly"),
        ("geometric stiffness", "geometric stiffness"),
        ("quasi stationary", "quasi stationary"),
        ("calculate loads", "load postprocessing"),
    ] {
        if flag(Some(p), key) {
            return Err(unsupported(what));
        }
    }
    Ok(())
}

/// Which postprocessing fields the model asks for.
#[derive(Copy, Clone, Debug, Default)]
struct PostFlags {
    stress: bool,
    strain: bool,
    principal: bool,
    /// Some equation names `Calculate Stresses`, so only those bodies project.
    masked: bool,
}

impl PostFlags {
    fn any(&self) -> bool {
        self.stress || self.strain
    }
}

fn post_flags(ctx: &SolverCtx<'_>) -> PostFlags {
    let p = ctx.params();
    let masked = ctx.sim.model.equations.iter().any(|e| flag(Some(e), "calculate stresses"));
    PostFlags {
        stress: flag(Some(p), "calculate stresses") || masked,
        strain: flag(Some(p), "calculate strains"),
        principal: flag(Some(p), "calculate principal"),
        masked,
    }
}

/// The postprocessing fields `ElasticSolver_Init` exports
/// (`Stress.F90:1656` `StressFieldDefinition`): name, dofs and component names.
fn field_specs(cdim: usize, f: PostFlags) -> Vec<VarSpec> {
    let comp = ["xx", "yy", "zz", "xy", "yz", "xz"];
    let n = crate::post::output_components(cdim);
    let spec = |name: &str, components: Vec<String>| VarSpec {
        name: name.to_string(),
        dofs: components.len().max(1),
        components,
        kind: FieldKind::Nodal,
        output: true,
    };
    let tensor = |name: &str| spec(name, comp[..n].iter().map(|c| format!("{name}_{c}")).collect());
    let triple = |name: &str| spec(name, (1..=3).map(|k| format!("{name} {k}")).collect());
    let mut out = Vec::new();
    if f.stress {
        out.push(tensor("stress"));
        out.push(spec("vonmises", Vec::new()));
        if f.principal {
            out.push(triple("principal stress"));
            out.push(spec("tresca", Vec::new()));
        }
    }
    if f.strain {
        out.push(tensor("strain"));
        if f.principal {
            out.push(triple("principal strain"));
        }
    }
    out
}

struct ElementData {
    law: ElasticLaw,
    density: Vec<f64>,
    damping: Option<Vec<f64>>,
    rayleigh_alpha: Option<f64>,
    rayleigh_beta: Option<f64>,
    load: Vec<[f64; 3]>,
    pressure: Option<Vec<f64>>,
    temperature: Option<Vec<f64>>,
    heat_expansion: Option<Vec<[f64; 3]>>,
    ip: Option<IpProps>,
}

fn nodal_temperature(temp: &Option<FieldRef<'_>>, nodes: &[u32], reference: &[f64]) -> Option<Vec<f64>> {
    let t = temp.as_ref()?;
    let perm = t.perm?;
    Some(
        nodes
            .iter()
            .zip(reference)
            .map(|(&n, r)| match perm.get(n as usize) {
                Some(&s) if s != elmer_core::NO_DOF => t.values[s as usize] - r,
                _ => 0.0,
            })
            .collect(),
    )
}

fn parent_center(basis: &dyn ElementBasis, mesh: &dyn MeshView, parent: ElemId) -> Result<[f64; 3], SolverError> {
    let pe = mesh.element(parent);
    let uvw = match pe.etype.family() {
        3 | 7 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        5 => [0.25, 0.25, 0.25],
        6 => [0.0, 0.0, 1.0 / 3.0],
        _ => [0.0; 3],
    };
    let xyz: Vec<[f64; 3]> = pe.nodes.iter().map(|&n| mesh.coord(n)).collect();
    let mut be = elmer_core::BasisEval::default();
    basis
        .eval(pe.etype, &xyz, uvw, &mut be)
        .map_err(|e| SolverError::fatal(CALLER, e))?;
    Ok(std::array::from_fn(|a| (0..xyz.len()).map(|k| be.basis[k] * xyz[k][a]).sum()))
}

/// The 6x6 `C` of one node from a `Youngs Modulus` block, as `InputTensor`
/// expands a vector or a matrix.
fn c_of_block(rows: usize, cols: usize, d: &[f64]) -> [[f64; 6]; 6] {
    if rows == 1 || cols == 1 {
        return diagonal_c(d);
    }
    let mut c = diagonal_c(&[]);
    for i in 0..rows.min(6) {
        for j in 0..cols.min(6) {
            c[i][j] = d[i * cols + j];
        }
    }
    c
}

/// Diagonal of a `Heat Expansion Coefficient` block.
fn beta_of_block(r: usize, c: usize, d: &[f64]) -> [f64; 3] {
    std::array::from_fn(|i| {
        if r == 1 && c == 1 {
            d[0]
        } else if r == 1 || c == 1 {
            d.get(i).copied().unwrap_or(0.0)
        } else if i < r && i < c {
            d[i * c + i]
        } else {
            0.0
        }
    })
}

impl StressSolver {
    #[allow(clippy::too_many_arguments)]
    fn element_data(
        &self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        cdim: usize,
        axi: bool,
        temp: &Option<FieldRef<'_>>,
    ) -> Result<ElementData, SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let bp = ctx
            .body_params(e)
            .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
        let mat = bp
            .material
            .ok_or_else(|| SolverError::fatal(CALLER, format!("body {} has no material", bp.body_id)))?;
        for (key, what) in [("maxwell material", "Maxwell material"), ("umat subroutine", "UMAT")] {
            if mat.get(key).is_some() && (flag(Some(mat), key) || !matches!(mat.get(key), Some(Value::Logical(_)))) {
                return Err(unsupported(what));
            }
        }
        if let Some(bf) = bp.body_force {
            for prefix in [
                "stress load",
                "strain load",
                "inertial bodyforce",
                "stress volume source",
                "gravitational prestress",
            ] {
                if present_prefix(bf, prefix) {
                    return Err(unsupported(&format!("body force '{prefix}'")));
                }
            }
        }

        let young_ip = flag(Some(mat), "youngs modulus at ip");
        let poisson_ip = flag(Some(mat), "poisson ratio at ip");
        let beta_ip = flag(Some(mat), "heat expansion coefficient ip");
        let load_ip = flag(bp.body_force, "stress bodyforce at ip");
        let any_ip = young_ip || poisson_ip || beta_ip || load_ip;
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
        let shape = if any_ip {
            let fe = self.fe();
            crate::bulk::shape_at_points(&*fe.basis, &*fe.quad, el.etype, &xyz, None)
                .map_err(|x| SolverError::fatal(CALLER, x.to_string()))?
        } else {
            Vec::new()
        };
        let mut ip = IpProps::default();

        let law = if young_ip {
            // Isotropy is assumed with it, as in ElasticSolve.
            ip.young = props::scalar_ips(ctx, mat, "youngs modulus", e, &shape)?;
            let poisson = if poisson_ip {
                ip.poisson = props::scalar_ips(ctx, mat, "poisson ratio", e, &shape)?;
                vec![0.0; n]
            } else {
                props::scalar_nodes(ctx, mat, "poisson ratio", e)?
                    .ok_or_else(|| SolverError::fatal(CALLER, "material has no 'Poisson Ratio'"))?
            };
            ElasticLaw::Isotropic { young: vec![0.0; n], poisson }
        } else {
            let ym = props::node_array(ctx, mat, "youngs modulus", e)?
                .ok_or_else(|| SolverError::fatal(CALLER, "material has no 'Youngs Modulus'"))?;
            if ym.rows == 1 && ym.cols == 1 {
                let poisson = if poisson_ip {
                    ip.poisson = props::scalar_ips(ctx, mat, "poisson ratio", e, &shape)?;
                    vec![0.0; n]
                } else {
                    props::scalar_nodes(ctx, mat, "poisson ratio", e)?
                        .ok_or_else(|| SolverError::fatal(CALLER, "material has no 'Poisson Ratio'"))?
                };
                ElasticLaw::Isotropic { young: ym.data, poisson }
            } else {
                if cdim != 3 || axi {
                    return Err(unsupported("anisotropic material outside 3D (ElasticSolve refuses it too)"));
                }
                let c: Vec<[[f64; 6]; 6]> = (0..n).map(|k| c_of_block(ym.rows, ym.cols, ym.node(k))).collect();
                let rotation = if flag(Some(mat), "rotate elasticity tensor") {
                    let mut t: Tensor3 = [[0.0; 3]; 3];
                    for (i, row) in t.iter_mut().enumerate() {
                        let v = props::node_array(ctx, mat, &format!("material coordinates unit vector {}", i + 1), e)?
                            .ok_or_else(|| {
                                SolverError::fatal(CALLER, "'Rotate Elasticity Tensor' without unit vectors")
                            })?;
                        let v = v.node(0);
                        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                        if norm < f64::EPSILON {
                            return Err(SolverError::fatal(CALLER, "material coordinate unit vector too short"));
                        }
                        for (a, x) in row.iter_mut().enumerate() {
                            *x = v[a] / norm;
                        }
                    }
                    Some(t)
                } else {
                    None
                };
                ElasticLaw::Anisotropic { c, rotation }
            }
        };
        if (young_ip || poisson_ip) && !matches!(law, ElasticLaw::Isotropic { .. }) {
            return Err(unsupported("properties at integration points need an isotropic material"));
        }

        let density = props::scalar_nodes(ctx, mat, "density", e)?.unwrap_or_else(|| vec![0.0; n]);
        let damping = props::scalar_nodes(ctx, mat, "damping", e)?;
        let rayleigh_alpha = props::const_real(ctx, mat, "rayleigh damping alpha")?;
        let rayleigh_beta = props::const_real(ctx, mat, "rayleigh damping beta")?;

        let mut load = vec![[0.0; 3]; n];
        let mut pressure = None;
        if let Some(bf) = bp.body_force {
            if load_ip {
                let mut l = vec![[0.0; 3]; shape.len() / n];
                for i in 0..cdim {
                    if let Some(v) = props::scalar_ips(ctx, bf, &format!("stress bodyforce {}", i + 1), e, &shape)? {
                        for (li, x) in l.iter_mut().zip(v) {
                            li[i] = x;
                        }
                    }
                }
                ip.load = Some(l);
                ip.pressure = Some(
                    props::scalar_ips(ctx, bf, "stress pressure", e, &shape)?
                        .unwrap_or_else(|| vec![0.0; shape.len() / n]),
                );
            } else {
                for i in 0..cdim {
                    if let Some(v) = props::scalar_nodes(ctx, bf, &format!("stress bodyforce {}", i + 1), e)? {
                        for (l, x) in load.iter_mut().zip(v) {
                            l[i] = x;
                        }
                    }
                }
                pressure = props::scalar_nodes(ctx, bf, "stress pressure", e)?;
            }
        }

        let reference =
            props::scalar_nodes(ctx, mat, "reference temperature", e)?.unwrap_or_else(|| vec![0.0; n]);
        let temperature = nodal_temperature(temp, el.nodes, &reference);
        let mut heat_expansion = None;
        if beta_ip {
            if temperature.is_some() {
                ip.beta = props::scalar_ips(ctx, mat, "heat expansion coefficient", e, &shape)?;
            }
        } else if let (Some(b), true) = (props::node_array(ctx, mat, "heat expansion coefficient", e)?, temperature.is_some()) {
            heat_expansion = Some((0..n).map(|k| beta_of_block(b.rows, b.cols, b.node(k))).collect());
        }
        let ip = (ip.young.is_some()
            || ip.poisson.is_some()
            || ip.beta.is_some()
            || ip.load.is_some())
        .then_some(ip);
        Ok(ElementData {
            law,
            density,
            damping,
            rayleigh_alpha,
            rayleigh_beta,
            load,
            pressure,
            temperature,
            heat_expansion,
            ip,
        })
    }

    /// Temperature snapshot of the primary-variable-independent field, taken
    /// before the mutable assembly.
    fn temperature_field(ctx: &SolverCtx<'_>) -> Option<(Vec<u32>, Vec<f64>)> {
        ctx.variable("temperature")
            .and_then(|f| f.perm.map(|p| (p.to_vec(), f.values.to_vec())))
    }

    /// `ComputeStressAndStrain`: L2 projection of the elastic strain and stress
    /// onto the displacement's nodes, then von Mises and the principal values.
    fn postprocess(
        &self,
        ctx: &mut SolverCtx<'_>,
        flags: PostFlags,
        cdim: usize,
        axi: bool,
        temp_ref: &Option<FieldRef<'_>>,
    ) -> Result<(), SolverError> {
        let coord = if axi { CoordSystem::AxiSymmetric } else { CoordSystem::Cartesian };
        let ncomp = post::output_components(cdim);
        let (disp_perm, disp_vals, dofs) = {
            let f = ctx
                .primary()
                .ok_or_else(|| SolverError::fatal(CALLER, "solver has no displacement variable"))?;
            (
                f.perm
                    .ok_or_else(|| SolverError::fatal(CALLER, "displacement has no permutation"))?
                    .to_vec(),
                f.values.to_vec(),
                f.dofs,
            )
        };
        let nslots = disp_vals.len() / dofs;
        let mut mass = post::MassSystem::new(nslots);
        let mut fstrain = vec![0.0; nslots * ncomp];
        let mut fstress = vec![0.0; nslots * ncomp];

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let fe = self.fe();
        for e in elems {
            let body = ctx.body_params(e);
            if flags.masked && !flag(body.as_ref().and_then(|b| b.equation), "calculate stresses") {
                continue;
            }
            let plane_stress = flag(body.as_ref().and_then(|b| b.equation), "plane stress");
            let el = ctx.mesh.element(e);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let d = self.element_data(ctx, e, cdim, axi, temp_ref)?;
            let disp: Vec<[f64; 3]> = el
                .nodes
                .iter()
                .map(|&n| match disp_perm.get(n as usize) {
                    Some(&s) if s != elmer_core::NO_DOF => {
                        std::array::from_fn(|a| if a < dofs.min(3) { disp_vals[s as usize * dofs + a] } else { 0.0 })
                    }
                    _ => [0.0; 3],
                })
                .collect();
            let inp = StressPostInput {
                etype: el.etype,
                xyz: &xyz,
                cdim,
                coord,
                plane_stress,
                law: &d.law,
                disp: &disp,
                temperature: d.temperature.as_deref(),
                heat_expansion: d.heat_expansion.as_deref(),
                ip: d.ip.as_ref(),
                want_strain: flags.strain,
                want_stress: flags.stress,
                np: None,
            };
            let loc = stress_strain_local(&*fe.basis, &*fe.quad, &inp).map_err(|x| SolverError::fatal(CALLER, x.to_string()))?;
            let slots: Vec<Option<usize>> = el
                .nodes
                .iter()
                .map(|&n| disp_perm.get(n as usize).copied().filter(|&s| s != elmer_core::NO_DOF).map(|s| s as usize))
                .collect();
            mass.add(&slots, &loc.mass);
            for (p, s) in slots.iter().enumerate() {
                let Some(s) = s else { continue };
                for k in 0..ncomp {
                    if flags.strain {
                        fstrain[ncomp * s + k] += loc.strain[ncomp * p + k];
                    }
                    if flags.stress {
                        fstress[ncomp * s + k] += loc.stress[ncomp * p + k];
                    }
                }
            }
        }

        let project = |f: &[f64]| -> Result<Vec<f64>, SolverError> {
            let mut out = vec![0.0; nslots * ncomp];
            let mut b = vec![0.0; nslots];
            let mut x = vec![0.0; nslots];
            for k in 0..ncomp {
                for s in 0..nslots {
                    b[s] = f[ncomp * s + k];
                }
                mass.solve(&b, &mut x).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
                for s in 0..nslots {
                    out[ncomp * s + k] = x[s];
                }
            }
            Ok(out)
        };
        let set = |ctx: &mut SolverCtx<'_>, name: &str, v: Vec<f64>| -> Result<(), SolverError> {
            let var = ctx
                .sim
                .vars
                .get_mut(name)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("variable '{name}' does not exist")))?;
            if var.values.len() != v.len() {
                return Err(SolverError::fatal(CALLER, format!("variable '{name}' has the wrong size")));
            }
            var.values = v;
            var.values_changed = true;
            Ok(())
        };

        if flags.stress {
            ctx.info(CALLER, "Calculating stress components", 7);
            let stress = project(&fstress)?;
            let mut vm = vec![0.0; nslots];
            let (mut pri, mut tre) = (vec![0.0; nslots * 3], vec![0.0; nslots]);
            for s in 0..nslots {
                let t = post::vector_to_tensor(&stress[ncomp * s..ncomp * (s + 1)], ncomp);
                vm[s] = post::von_mises(&t);
                if flags.principal {
                    let w = post::principal_values(&t);
                    pri[3 * s..3 * s + 3].copy_from_slice(&w);
                    tre[s] = post::tresca(&w);
                }
            }
            set(ctx, "stress", stress)?;
            set(ctx, "vonmises", vm)?;
            if flags.principal {
                set(ctx, "principal stress", pri)?;
                set(ctx, "tresca", tre)?;
            }
        }
        if flags.strain {
            ctx.info(CALLER, "Calculating strain components", 7);
            let strain = project(&fstrain)?;
            if flags.principal {
                let mut pri = vec![0.0; nslots * 3];
                for s in 0..nslots {
                    let t = post::vector_to_tensor(&strain[ncomp * s..ncomp * (s + 1)], ncomp);
                    pri[3 * s..3 * s + 3].copy_from_slice(&post::principal_values(&t));
                }
                set(ctx, "principal strain", pri)?;
            }
            set(ctx, "strain", strain)?;
        }
        Ok(())
    }
}

impl Solver for StressSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_solver_params(ctx.params())?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(dim) {
            return Err(SolverError::fatal(
                CALLER,
                format!("'Variable DOFs' must equal the coordinate dimension {dim}, got {dofs:?}"),
            ));
        }
        // ElasticSolver_Init: the time derivative order is two unless the sif says otherwise.
        if ctx.params().get("time derivative order").is_none() {
            ctx.set_time_derivative_order(2)?;
        }
        // The postprocessing fields share the displacement's permutation.
        let flags = post_flags(ctx);
        if flags.any() {
            let name = ctx.variable_name().unwrap_or_default().to_string();
            let (perm, nslots) = {
                let v = ctx
                    .sim
                    .vars
                    .get(&name)
                    .ok_or_else(|| SolverError::fatal(CALLER, "solver has no displacement variable"))?;
                (v.perm.clone(), v.values.len() / v.dofs)
            };
            for spec in field_specs(dim, flags) {
                if ctx.sim.vars.id(&spec.name).is_some() {
                    continue;
                }
                let mut var = Variable::new(&spec.name, FieldKind::Nodal, spec.dofs, perm.clone(), vec![0.0; nslots * spec.dofs]);
                var.output = spec.output;
                ctx.sim.vars.add_with_components(var, &spec.components)?;
            }
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        check_solver_params(ctx.params())?;
        let second_order = ctx.params().integer("time derivative order") == Some(2);
        let cdim = ctx.mesh.dim();
        let coord_name = elmer_core::normalize_key(ctx.simulation_params().string("coordinate system").unwrap_or("cartesian"));
        let axi = match coord_name.as_str() {
            "axi symmetric" => true,
            s if s.starts_with("cartesian") => false,
            other => return Err(unsupported(&format!("coordinate system '{other}'"))),
        };
        let coord = if axi { CoordSystem::AxiSymmetric } else { CoordSystem::Cartesian };

        ctx.default_initialize()?;

        // Snapshot the temperature field (read-only) before the mutable assembly.
        let temp_owned = Self::temperature_field(ctx);
        let temp_ref = temp_owned.as_ref().map(|(p, v)| FieldRef {
            name: "temperature",
            kind: elmer_core::FieldKind::Nodal,
            dofs: 1,
            perm: Some(p),
            values: v,
        });

        let flags = post_flags(ctx);
        if ctx.solver_dofs().is_some() && flags.any() {
            return Err(unsupported("stress postprocessing (Calculate Stresses) on p elements"));
        }
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for e in elems {
            let el = ctx.mesh.element(e);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let pe = match ctx.solver_dofs() {
                Some(d) => Some(p_element(d, e)?),
                None => None,
            };
            let d = self.element_data(ctx, e, cdim, axi, &temp_ref)?;
            let plane_stress = flag(ctx.body_params(e).and_then(|b| b.equation), "plane stress");
            let inp = StressBulkInput {
                etype: el.etype,
                xyz: &xyz,
                cdim,
                coord,
                plane_stress,
                law: &d.law,
                density: &d.density,
                damping: d.damping.as_deref(),
                load: &d.load,
                pressure_load: d.pressure.as_deref(),
                stress_load: None,
                strain_load: None,
                temperature: d.temperature.as_deref(),
                heat_expansion: d.heat_expansion.as_deref(),
                want_mass: transient,
                np: None,
                ip: d.ip.as_ref(),
            };
            let mut local = LocalSystem::zeros(0);
            let fe = self.fe();
            match &pe {
                // StressCompose(..., n, ntot = GetElementNOFDOFs(), ...); the
                // displacement-only formulation neither adds nor condenses
                // local bubbles.
                Some(p) => stress_bulk_local_p(&*p.basis, &p.rule, p.nd, &inp, &mut local),
                None => stress_bulk_local(&*fe.basis, &*fe.quad, &inp, &mut local),
            }
            .map_err(|x| SolverError::fatal(CALLER, x.to_string()))?;
            if transient {
                let nn = local.n * local.n;
                let mass = local.mass.take().unwrap_or_else(|| vec![0.0; nn]);
                let mut damp = local.damp.take().unwrap_or_else(|| vec![0.0; nn]);
                if let Some(a) = d.rayleigh_alpha {
                    damp.iter_mut().zip(&mass).for_each(|(c, m)| *c += a * m);
                }
                if let Some(b) = d.rayleigh_beta {
                    damp.iter_mut().zip(&local.stiff).for_each(|(c, k)| *c += b * k);
                }
                if second_order {
                    // `Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )`, Rayleigh terms included.
                    local.mass = Some(mass);
                    local.damp = Some(damp);
                    ctx.default_second_order_time(&mut local, e, dt)?;
                } else {
                    // First order in time: the damping matrix carries d/dt
                    // (`Default1stOrderTime( DAMP, STIFF, FORCE )`), Rayleigh terms included.
                    local.mass = Some(damp);
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
            }
            ctx.default_update(&local, e)?;
        }

        let belems: Vec<ElemId> = ctx.boundary_elements().collect();
        for e in belems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let bc = ctx.bc_params(e).expect("boundary_elements filters on BC").clone();
            // A p boundary element lists its nodes in the parent's edge or
            // face order (AssignLocalNumber); its values follow that order.
            let p_nodes: Option<Vec<u32>> = ctx.solver_dofs().map(|d| d.element_nodes(e));
            if let Some(v) = bc.logical("fsi bc") {
                if v {
                    return Err(unsupported("FSI boundary condition"));
                }
            }
            for (prefix, what) in [
                ("stress load", "boundary 'Stress Load'"),
                ("model lumping", "model lumping boundary"),
                ("normal-tangential", "normal-tangential displacement"),
            ] {
                if present_prefix(&bc, prefix) {
                    return Err(unsupported(what));
                }
            }
            let nodal = |key: &str| -> Result<Option<Vec<f64>>, SolverError> {
                let v = props::scalar_nodes(ctx, &bc, key, e)?;
                Ok(match (&p_nodes, v) {
                    (Some(pn), Some(v)) => Some(reorder(v, el.nodes, pn)),
                    (_, v) => v,
                })
            };
            let mut traction = vec![[0.0; 3]; n];
            let mut got = false;
            for prefix in ["surface traction", "force"] {
                if present_prefix(&bc, prefix) {
                    for i in 0..cdim {
                        if let Some(v) = nodal(&format!("{prefix} {}", i + 1))? {
                            for (t, x) in traction.iter_mut().zip(v) {
                                t[i] = x;
                            }
                        }
                    }
                    got = true;
                    break;
                }
            }
            let mut beta = nodal("normal surface traction")?;
            if beta.is_none() {
                beta = nodal("normal force")?;
            }
            got |= beta.is_some();

            let mut spring_normal = None;
            let mut spring_tensor = None;
            if present_prefix(&bc, "spring") {
                if let Some(k) = nodal("spring")? {
                    spring_normal = Some(k);
                } else {
                    let mut t = vec![[[0.0; 3]; 3]; n];
                    for i in 0..cdim {
                        if let Some(v) = nodal(&format!("spring {}", i + 1))? {
                            for (tk, x) in t.iter_mut().zip(v) {
                                tk[i][i] = x;
                            }
                        }
                        for j in 0..cdim {
                            if let Some(v) = nodal(&format!("spring {}{}", i + 1, j + 1))? {
                                for (tk, x) in t.iter_mut().zip(v) {
                                    tk[i][j] = x;
                                }
                            }
                        }
                    }
                    spring_tensor = Some(t);
                }
                got = true;
            }
            if !got {
                continue;
            }

            let parent = el.parents[0]
                .or(el.parents[1])
                .ok_or_else(|| SolverError::fatal(CALLER, format!("boundary element {} has no parent", e.0)))?;
            let fe = self.fe();
            let center = parent_center(&*fe.basis, ctx.mesh, parent)?;
            let nodes: &[u32] = p_nodes.as_deref().unwrap_or(el.nodes);
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
            let nt = beta.unwrap_or_else(|| vec![0.0; n]);
            let spring = if let Some(sn) = &spring_normal {
                Spring::Normal(sn)
            } else if let Some(st) = &spring_tensor {
                Spring::Tensor(st)
            } else {
                Spring::None
            };
            let inp = StressBoundaryInput {
                etype: ElementType(el.etype.0),
                xyz: &xyz,
                cdim,
                coord,
                traction: &traction,
                normal_traction: &nt,
                spring,
                parent_center: center,
                np: None,
            };
            let mut local = LocalSystem::zeros(0);
            match ctx.solver_dofs() {
                Some(d) => {
                    let p = p_element(d, e)?;
                    stress_boundary_local_p(&*p.basis, &p.rule, p.nd, &inp, &mut local)
                }
                None => stress_boundary_local(&*fe.basis, &*fe.quad, &inp, &mut local),
            }
            .map_err(|x| SolverError::fatal(CALLER, x.to_string()))?;
            if transient {
                local.mass = Some(vec![0.0; local.n * local.n]);
                if second_order {
                    local.damp = Some(vec![0.0; local.n * local.n]);
                    ctx.default_second_order_time(&mut local, e, dt)?;
                } else {
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
            }
            ctx.default_update(&local, e)?;
        }

        ctx.default_dirichlet()?;
        ctx.default_solve()?;

        if flags.any() {
            self.postprocess(ctx, flags, cdim, axi, &temp_ref)?;
        }
        Ok(())
    }
}
