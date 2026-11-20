//! `VectorHelmholtzNodal.F90` (fem/src/modules), the segregated nodal
//! solver of the time-harmonic electric field: `Procedure =
//! "VectorHelmholtzNodal" "VectorHelmholtzNodal"` without `Monolithic
//! Solver`.
//!
//! [`init_section`] is `VectorHelmholtzNodal_Init`: the complex scalar
//! solver variable `Etmp` (not saved), the exported field `Elfield` with
//! the six components `EX re ... EZ im`, `Linear System Complex`.
//!
//! The solver solves one complex scalar system per component `i` of the
//! field. Bulk: `(1/mu) grad u . grad v - kappa omega^2 eps u v`. Boundary:
//! `Electric Robin Coefficient` (+ `im`) or, with `Absorbing BC`,
//! `i omega sqrt(eps0 mu0)`, scaled by `sqrt(1 - n_i^2)`; and `Magnetic
//! Boundary Load i` (+ `im`) on the right-hand side. The starting guess is
//! that component of `Elfield`, and the solution is written back to it.
//! The variable norm is `SQRT(SUM(Norm(1:dim)**2)/dim)` of the component
//! solves.
//!
//! Refused: the monolithic and curl-curl forms, the use as a
//! preconditioner, ports, good-conductor / layer / TEM / normal-tangential
//! conditions, conductivity, current densities, `Relative Reluctivity`,
//! Dirichlet conditions on the field and `Relative Integration Order`.

use num_complex::Complex64 as C;

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::metric::MetricStatus;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::case::has_prefix;

const CALLER: &str = "VectorHelmholtzNodal";

pub const PROC: (&str, &str) = ("VectorHelmholtzNodal", "VectorHelmholtzNodal");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROC.0, PROC.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(VhNodal))
}

pub fn is_section(p: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(p)
        .is_some_and(|(l, e)| l.trim_start_matches("./").eq_ignore_ascii_case(PROC.0) && e.eq_ignore_ascii_case(PROC.1))
}

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

/// `VectorHelmholtzNodal_Init` for the segregated solver (the solve
/// refuses the other forms, whose sections are left alone).
pub fn init_section(p: &mut ParamMap) {
    let on = |k: &str| p.logical(k) == Some(true);
    if on("monolithic solver") || on("curl-curl form") || on("preconditioning solver") {
        return;
    }
    if p.get("variable output").is_none() {
        p.insert("Variable Output", Value::Logical(false));
    }
    if p.get("variable").is_none() {
        p.insert("Variable", Value::Str("Etmp[Etmp re:1 Etmp im:1]".into()));
    }
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    p.insert(&format!("Exported Variable {k}"), Value::Str("Elfield[EX re:1 EX im:1 EY re:1 EY im:1 EZ re:1 EZ im:1]".into()));
    if p.get("linear system complex").is_none() {
        p.insert("Linear System Complex", Value::Logical(true));
    }
    p.insert("Time Derivative Order", Value::Integer(vec![0]));
}

struct VhNodal;

fn check(ctx: &SolverCtx<'_>, sp: &ParamMap) -> Result<(), SolverError> {
    for k in ["monolithic solver", "curl-curl form", "preconditioning solver", "follow p curvature"] {
        if sp.logical(k) == Some(true) {
            return Err(unsupported(format!("'{k}'")));
        }
    }
    if sp.get("relative integration order").is_some() {
        return Err(unsupported("'Relative Integration Order'"));
    }
    let model = &ctx.sim.model;
    let cs = model.simulation.string("coordinate system").unwrap_or("cartesian").to_ascii_lowercase();
    if !cs.trim().starts_with("cartesian") {
        return Err(SolverError::fatal(CALLER, "Only implemented for Cartesian problems!"));
    }
    for m in &model.materials {
        for k in ["electric conductivity", "relative reluctivity"] {
            if has_prefix(m, k) {
                return Err(unsupported(format!("material '{k}'")));
            }
        }
    }
    for bf in &model.body_forces {
        if has_prefix(bf, "current density") {
            return Err(unsupported("body force 'Current Density'"));
        }
    }
    const BC_KEYS: &[&str] = &[
        "port type",
        "good conductor bc",
        "tem potential",
        "layer electric conductivity",
        "layer relative reluctivity",
        "normal-tangential",
        "eigenfunction",
        "ex re",
        "ex im",
        "ey re",
        "ey im",
        "ez re",
        "ez im",
        "elfield",
        "etmp",
    ];
    for bc in &model.boundary_conditions {
        for k in BC_KEYS {
            if has_prefix(bc, k) {
                return Err(unsupported(format!("boundary condition '{k}'")));
            }
        }
        for key in bc.0.keys().filter(|k| k.starts_with("magnetic boundary load")) {
            let ok = (1..=3).any(|c| *key == format!("magnetic boundary load {c}") || *key == format!("magnetic boundary load im {c}"));
            if !ok {
                return Err(unsupported(format!("boundary condition '{key}'")));
            }
        }
    }
    Ok(())
}

/// `SUM(Basis(1:n) * v(1:n))`.
fn interp(basis: &[f64], v: &[f64]) -> f64 {
    let mut s = 0.0;
    for (b, x) in basis.iter().zip(v) {
        s += b * x;
    }
    s
}

/// `ListGetElementComplex(h, Basis, Element)` of `key` (and `key im`):
/// `None` when neither is given.
fn complex_at_ip(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str, im_key: &str, e: ElemId, basis: &[f64]) -> Result<Option<C>, SolverError> {
    let re = ctx.real_at_nodes(p, key, e)?;
    let im = ctx.real_at_nodes(p, im_key, e)?;
    if re.is_none() && im.is_none() {
        return Ok(None);
    }
    let r = re.map_or(0.0, |v| interp(basis, &v));
    let i = im.map_or(0.0, |v| interp(basis, &v));
    Ok(Some(C::new(r, i)))
}

/// `NormalVector(Element, Nodes)` of a linear triangle or quadrilateral
/// (unturned; only its squared components are used).
fn normal(xyz: &[[f64; 3]], code: u16) -> Result<[f64; 3], SolverError> {
    let (du, dv): (Vec<f64>, Vec<f64>) = match code {
        303 => (vec![-1.0, 1.0, 0.0], vec![-1.0, 0.0, 1.0]),
        404 => (vec![-0.25, 0.25, 0.25, -0.25], vec![-0.25, -0.25, 0.25, 0.25]),
        c => return Err(unsupported(format!("boundary elements of type {c}"))),
    };
    let dot = |d: &[f64], a: usize| d.iter().zip(xyz).fold(0.0, |s, (di, xi)| s + di * xi[a]);
    let (xu, yu, zu) = (dot(&du, 0), dot(&du, 1), dot(&du, 2));
    let (xv, yv, zv) = (dot(&dv, 0), dot(&dv, 1), dot(&dv, 2));
    let auu = xu * xu + yu * yu + zu * zu;
    let auv = xu * xv + yu * yv + zu * zv;
    let avv = xv * xv + yv * yv + zv * zv;
    let det = 1.0 / (auu * avv - auv * auv).sqrt();
    Ok([(yu * zv - yv * zu) * det, (xv * zu - xu * zv) * det, (xu * yv - xv * yu) * det])
}

/// The complex local system as the interleaved real one (`[re -im; im re]`
/// blocks).
fn to_real(n: usize, stiff: &[C], force: &[C]) -> LocalSystem {
    let m = 2 * n;
    let mut ls = LocalSystem { n: m, stiff: vec![0.0; m * m], mass: None, damp: None, force: vec![0.0; m] };
    for p in 0..n {
        for q in 0..n {
            let a = stiff[p * n + q];
            ls.stiff[2 * p * m + 2 * q] = a.re;
            ls.stiff[2 * p * m + 2 * q + 1] = -a.im;
            ls.stiff[(2 * p + 1) * m + 2 * q] = a.im;
            ls.stiff[(2 * p + 1) * m + 2 * q + 1] = a.re;
        }
        ls.force[2 * p] = force[p].re;
        ls.force[2 * p + 1] = force[p].im;
    }
    ls
}

fn eval(code: u16, xyz: &[[f64; 3]], uvw: [f64; 3], status: &mut MetricStatus) -> Result<BasisEval, SolverError> {
    let def = elmer_fe::element_type(code).map_err(|e| SolverError::fatal("ElementInfo", e.to_string()))?;
    let mut be = BasisEval::default();
    elmer_fe::info::element_info(def, xyz, uvw, 3, status, &mut be).map_err(|e| SolverError::fatal("ElementInfo", e.to_string()))?;
    Ok(be)
}

impl Solver for VhNodal {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("a transient run"));
        }
        let sp = ctx.params().clone();
        check(ctx, &sp)?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map_or(0, |v| v.dofs) / 2;
        if dofs != 1 {
            return Err(SolverError::fatal(CALLER, format!("Invalid number of dofs for the solver variable: {dofs}")));
        }
        let ef_name = "elfield";
        let ef_id = ctx.sim.vars.id(ef_name).ok_or_else(|| SolverError::fatal(CALLER, "Variable for Electric field not found!"))?;
        if ctx.sim.vars.by_id(ef_id).dofs != 2 * dim {
            return Err(SolverError::fatal(CALLER, "the electric field has the wrong number of components"));
        }

        // InitStuff.
        let mu0 = ctx.constants().real("permeability of vacuum");
        let mu0inv = 1.0 / mu0.unwrap_or(std::f64::consts::PI * 4.0e-7);
        let eps0 = ctx.constants().real("permittivity of vacuum").unwrap_or(8.854187817e-12);
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
        let first = *active.first().ok_or_else(|| SolverError::fatal(CALLER, "no active elements"))?;
        let omega = crate::calcfields2d::angular_frequency(ctx, &sp, first)?
            .ok_or_else(|| SolverError::fatal(CALLER, "We should have Omega!"))?;
        let rob0 = omega * (eps0 / mu0inv).sqrt();
        let kappa = C::new(1.0, 0.0);
        let mu = C::new(mu0inv, 0.0);

        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut norms = vec![0.0f64; dim];
        let mut slots = Vec::new();
        let mut status = MetricStatus::default();
        for (compi, norm_i) in norms.iter_mut().enumerate() {
            ctx.default_initialize()?;
            for &e in &active {
                ctx.element_slots(e, &mut slots);
                if slots.contains(&NO_DOF) {
                    return Err(SolverError::fatal(CALLER, format!("element {} lacks DOFs", e.0)));
                }
                let el = ctx.mesh.element(e);
                let (code, n) = (el.etype.0, el.nodes.len());
                let nd = slots.len();
                if nd != n {
                    return Err(unsupported("elements with other than nodal DOFs"));
                }
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let material = ctx.body_params(e).and_then(|b| b.material).cloned();
                let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                let mut stiff = vec![C::new(0.0, 0.0); nd * nd];
                let force = vec![C::new(0.0, 0.0); nd];
                for t in 0..rule.len() {
                    let be = eval(code, &xyz, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
                    let (basis, db) = (&be.basis[..n], &be.dbasis_dx[..n]);
                    let weight = rule.s[t] * be.det_j;
                    let eps = match &material {
                        Some(m) => complex_at_ip(ctx, m, "relative permittivity", "relative permittivity im", e, basis)?,
                        None => None,
                    };
                    let eps = match eps {
                        Some(x) => x * eps0,
                        None => C::new(eps0, 0.0),
                    };
                    // Weight * MuinvAtIp * MATMUL(dBasisdx, TRANSPOSE(dBasisdx)).
                    let wm = mu * weight;
                    for p in 0..nd {
                        for q in 0..nd {
                            let mut g = 0.0;
                            for k in 0..3 {
                                g += db[p][k] * db[q][k];
                            }
                            stiff[p * nd + q] += wm * g;
                        }
                    }
                    // - Weight * kappa * Omega**2 * epsAtIP * Basis(p) * Basis(1:nd).
                    let c = kappa * weight * (omega * omega) * eps;
                    for p in 0..nd {
                        let cp = c * basis[p];
                        for q in 0..nd {
                            stiff[p * nd + q] -= cp * basis[q];
                        }
                    }
                }
                let ls = to_real(nd, &stiff, &force);
                ctx.default_update(&ls, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            for &e in &boundary {
                let Some(bc) = ctx.bc_params(e).cloned() else { continue };
                ctx.element_slots(e, &mut slots);
                if slots.is_empty() || slots.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let (code, n) = (el.etype.0, el.nodes.len());
                let nd = slots.len();
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let nrm = normal(&xyz, code)?;
                let absorb = bc.logical("absorbing bc") == Some(true);
                let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                let mut stiff = vec![C::new(0.0, 0.0); nd * nd];
                let mut force = vec![C::new(0.0, 0.0); nd];
                for t in 0..rule.len() {
                    let be = eval(code, &xyz, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
                    let basis = &be.basis[..n];
                    let weight = rule.s[t] * be.det_j;
                    let mut l = [C::new(0.0, 0.0); 3];
                    let b = if absorb {
                        Some(C::new(0.0, rob0))
                    } else {
                        for (c, lc) in l.iter_mut().enumerate() {
                            let k = format!("magnetic boundary load {}", c + 1);
                            let ki = format!("magnetic boundary load im {}", c + 1);
                            *lc = complex_at_ip(ctx, &bc, &k, &ki, e, basis)?.unwrap_or(C::new(0.0, 0.0));
                        }
                        complex_at_ip(ctx, &bc, "electric robin coefficient", "electric robin coefficient im", e, basis)?
                    };
                    let Some(b) = b else { continue };
                    // FORCE(1:nd,i) - muinvAtIp * L(i) * Basis(1:nd) * Weight.
                    let ml = mu * l[compi];
                    for p in 0..nd {
                        force[p] -= ml * basis[p] * weight;
                    }
                    let c2 = 1.0 - nrm[compi] * nrm[compi];
                    let coeff = if c2 <= 0.0 { 0.0 } else { c2.sqrt() };
                    // STIFF(p,1:nd,i) - coeff * muinvAtIp * B * Basis(p) * Basis(1:nd) * Weight.
                    let cb = mu * coeff * b;
                    for p in 0..nd {
                        let cp = cb * basis[p];
                        for q in 0..nd {
                            stiff[p * nd + q] -= cp * basis[q] * weight;
                        }
                    }
                }
                let ls = to_real(nd, &stiff, &force);
                ctx.default_update(&ls, e)?;
            }
            ctx.default_finish_boundary_assembly()?;

            // EiVar % Values = this component of the field.
            let ef = ctx.sim.vars.by_id(ef_id);
            let (ef_vals, ef_perm) = (ef.values.clone(), ef.perm.clone());
            let ei = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            let ei_perm: Option<Vec<u32>> = ei.perm.map(|p| p.to_vec());
            copy_component(&ef_vals, ef_perm.as_deref(), 2 * dim, 2 * compi, ei.values, ei_perm.as_deref(), 2, 0)?;
            *norm_i = ctx.default_solve()?;
            let ei = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            let (ei_vals, ei_perm) = (ei.values.to_vec(), ei.perm.map(|p| p.to_vec()));
            let ef = ctx.sim.vars.by_id_mut(ef_id);
            let ef_perm = ef.perm.clone();
            copy_component(&ei_vals, ei_perm.as_deref(), 2, 0, &mut ef.values, ef_perm.as_deref(), 2 * dim, 2 * compi)?;
            ef.values_changed = true;
        }
        let mut s = 0.0;
        for x in &norms {
            s += x * x;
        }
        let id = ctx.sim.slots.iter().find(|s| s.id == ctx.solver_id()).and_then(|s| s.variable);
        if let Some(id) = id {
            ctx.sim.vars.by_id_mut(id).conv.norm = (s / dim as f64).sqrt();
        }
        Ok(())
    }
}

/// Copies the complex component starting at `from_off` of `from` (`from_dofs`
/// values per node) into the one at `to_off` of `to`, node by node through
/// the two permutations.
#[allow(clippy::too_many_arguments)]
fn copy_component(
    from: &[f64],
    from_perm: Option<&[u32]>,
    from_dofs: usize,
    from_off: usize,
    to: &mut [f64],
    to_perm: Option<&[u32]>,
    to_dofs: usize,
    to_off: usize,
) -> Result<(), SolverError> {
    match (from_perm, to_perm) {
        (Some(fp), Some(tp)) if fp.len() == tp.len() => {
            for (&a, &b) in fp.iter().zip(tp) {
                if a == NO_DOF || b == NO_DOF {
                    continue;
                }
                for k in 0..2 {
                    to[to_dofs * b as usize + to_off + k] = from[from_dofs * a as usize + from_off + k];
                }
            }
            Ok(())
        }
        (None, None) => {
            let n = from.len() / from_dofs;
            for i in 0..n {
                for k in 0..2 {
                    to[to_dofs * i + to_off + k] = from[from_dofs * i + from_off + k];
                }
            }
            Ok(())
        }
        _ => Err(SolverError::fatal(CALLER, "the electric field and the solver variable have different permutations")),
    }
}
