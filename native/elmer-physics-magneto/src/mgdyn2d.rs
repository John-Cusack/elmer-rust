//! `MagnetoDynamics2DHarmonic` (MagnetoDynamics2D.F90:1723): the
//! time-harmonic equation for the out-of-plane (or azimuthal, in axial
//! symmetry) magnetic vector potential on nodal elements,
//!
//! ```text
//!   -div(nu grad A) + i omega sigma A = J      (+ A/r terms in axial symmetry)
//! ```
//!
//! assembled through `elmer-solver`'s [`SolverCtx`] as an interleaved
//! complex system (`Potential[Potential re:1 Potential im:1]`, `Linear
//! System Complex`). Ported: `LocalMatrix` for linear isotropic materials
//! (reluctivity from `Relative Permeability`, `Permeability` or
//! `Reluctivity`, inverted as `Get2x2CMPLXTensorInverse` does), conductivity,
//! current density and magnetization loads; `LocalMatrixBC` (air gap,
//! surface current) and `LocalMatrixSkinBC`; Dirichlet conditions on the
//! two components by name. Refused: H-B curves, coils, London equations,
//! velocities, electrodynamics, the infinity BC, flux-density BCs,
//! imaginary material values and transient restarts.

use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, ElementType, FieldKind, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::metric::MetricStatus;
use elmer_solver::sim::node_order_perm;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};
use num_complex::Complex64;

const CALLER: &str = "MagnetoDynamics2DHarmonic";

/// `(library, entry)`.
pub const PROCEDURE: (&str, &str) = ("MagnetoDynamics2D", "MagnetoDynamics2DHarmonic");

/// The component names of the default variable
/// (`MagnetoDynamics2DHarmonic_Init`).
const COMPONENTS: [&str; 2] = ["potential re", "potential im"];

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MagnetoDynamics2DHarmonic))
}

pub struct MagnetoDynamics2DHarmonic;

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", what.into()))
}

fn slot_index(ctx: &SolverCtx<'_>) -> usize {
    let id = ctx.solver_id();
    ctx.sim.slots.iter().position(|s| s.id == id).expect("the running solver has a slot")
}

fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    let prefix = elmer_core::normalize_key(prefix);
    p.0.keys().any(|k| k.starts_with(&prefix))
}

fn check_supported(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for k in ["electrodynamics model", "transient restart", "apply integral bcs", "calculate loads"] {
        if p.logical(k) == Some(true) {
            return Err(unsupported(format!("solver option '{k}'")));
        }
    }
    let m = &ctx.sim.model;
    for mat in &m.materials {
        for k in ["h-b curve", "london lambda", "permittivity"] {
            if mat.get(k).is_some() {
                return Err(unsupported(format!("material '{k}'")));
            }
        }
        if mat.0.keys().any(|k| k.ends_with(" im")) {
            return Err(unsupported("imaginary material values"));
        }
    }
    for bf in &m.body_forces {
        for k in ["lorentz velocity", "angular velocity"] {
            if has_prefix(bf, k) {
                return Err(unsupported(format!("body force '{k}'")));
            }
        }
    }
    for bc in &m.boundary_conditions {
        if bc.logical("infinity bc") == Some(true) {
            return Err(unsupported("Infinity BC (needs NormalVector)"));
        }
        if has_prefix(bc, "magnetic flux density") {
            return Err(unsupported("Magnetic Flux Density BC"));
        }
    }
    Ok(())
}

/// `GetAngularFrequency(UElement=Element)`: `Frequency` (times 2 pi) or
/// `Angular Frequency` from the element's equation, material, then the
/// simulation and the solver section (Lists.F90:10543).
fn angular_frequency(ctx: &SolverCtx<'_>, e: ElemId) -> Result<f64, SolverError> {
    let bp = ctx.body_params(e);
    let lists = [
        bp.and_then(|b| b.equation),
        bp.and_then(|b| b.material),
        Some(ctx.simulation_params()),
        Some(ctx.params()),
    ];
    for p in lists.into_iter().flatten() {
        if let Some(f) = ctx.creal(p, "frequency")? {
            return Ok(2.0 * std::f64::consts::PI * f);
        }
        if let Some(w) = ctx.creal(p, "angular frequency")? {
            return Ok(w);
        }
    }
    Err(SolverError::fatal(CALLER, "Harmonic solution requires frequency!"))
}

/// `GetReluctivity` (MagnetoDynamics2D.F90:2790) for an isotropic scalar:
/// `Get2x2CMPLXTensorInverse` of `mu0 mu_r` (or `Permeability`) on the
/// diagonal, else `Reluctivity`. Returns the diagonal entry per node.
fn reluctivity(ctx: &SolverCtx<'_>, material: &ParamMap, e: ElemId, mu0: f64) -> Result<Vec<f64>, SolverError> {
    let inv = |a: f64| {
        // det = a*d - b*c with b = c = 0, then 1/det * d.
        let det = a * a - 0.0 * 0.0;
        1.0 / det * a
    };
    if let Some(mr) = ctx.real_at_nodes(material, "relative permeability", e)? {
        return Ok(mr.iter().map(|m| inv(mu0 * m)).collect());
    }
    if let Some(mu) = ctx.real_at_nodes(material, "permeability", e)? {
        return Ok(mu.iter().map(|&m| inv(m)).collect());
    }
    if let Some(nu) = ctx.real_at_nodes(material, "reluctivity", e)? {
        return Ok(nu);
    }
    Err(SolverError::fatal(
        "GetReluctivity",
        "Could not get either > Reluctivity > or > Relative Permeability < !",
    ))
}

fn interp(basis: &[f64], v: &[f64]) -> f64 {
    let mut s = 0.0;
    for (b, x) in basis.iter().zip(v) {
        s += b * x;
    }
    s
}

/// Complex `nd x nd` matrix and vector into Elmer's interleaved real local
/// system (`(2p, 2q) = re`, `(2p, 2q+1) = -im`, `(2p+1, 2q) = im`,
/// `(2p+1, 2q+1) = re`).
fn interleave(stiff: &[Complex64], force: &[Complex64], nd: usize, out: &mut LocalSystem) {
    let n = 2 * nd;
    *out = LocalSystem::zeros(n);
    for p in 0..nd {
        for q in 0..nd {
            let z = stiff[p * nd + q];
            out.stiff[(2 * p) * n + 2 * q] = z.re;
            out.stiff[(2 * p) * n + 2 * q + 1] = -z.im;
            out.stiff[(2 * p + 1) * n + 2 * q] = z.im;
            out.stiff[(2 * p + 1) * n + 2 * q + 1] = z.re;
        }
        out.force[2 * p] = force[p].re;
        out.force[2 * p + 1] = force[p].im;
    }
}

struct Point {
    s_detj: f64,
    basis: Vec<f64>,
    dbasis: Vec<[f64; 3]>,
    x: f64,
}

/// The Gauss points of `GaussPoints(Element)` with `ElementInfo` and the
/// axisymmetric `detJ * x`.
fn points(etype: ElementType, xyz: &[[f64; 3]], axi: bool) -> Result<Vec<Point>, SolverError> {
    let def = elmer_fe::element_type(etype.0).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
    let rule = elmer_fe::quadrature::gauss_points(etype, None).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
    let n = xyz.len();
    let mut status = MetricStatus::default();
    let mut out = Vec::with_capacity(rule.len());
    for t in 0..rule.len() {
        let mut be = BasisEval::default();
        elmer_fe::info::element_info(def, xyz, [rule.u[t], rule.v[t], rule.w[t]], 2, &mut status, &mut be)
            .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
        let mut det_j = be.det_j;
        let mut x = 0.0;
        if axi {
            x = interp(&be.basis[..n], &xyz.iter().map(|c| c[0]).collect::<Vec<_>>());
            det_j *= x;
        }
        out.push(Point {
            s_detj: rule.s[t] * det_j,
            basis: be.basis,
            dbasis: be.dbasis_dx,
            x,
        });
    }
    Ok(out)
}

/// Nodal coefficients of one bulk element.
struct Bulk {
    nu: Vec<f64>,
    sigma: Vec<f64>,
    m1: Vec<f64>,
    m2: Vec<f64>,
    load: Vec<f64>,
    omega: f64,
}

/// `LocalMatrix` (MagnetoDynamics2D.F90:2230), linear isotropic case.
fn bulk_local(etype: ElementType, xyz: &[[f64; 3]], axi: bool, c: &Bulk) -> Result<(Vec<Complex64>, Vec<Complex64>), SolverError> {
    let nd = xyz.len();
    let mut stiff = vec![Complex64::new(0.0, 0.0); nd * nd];
    let mut force = vec![Complex64::new(0.0, 0.0); nd];
    let im = Complex64::new(0.0, 1.0);
    for pt in points(etype, xyz, axi)? {
        let (b, db, x) = (&pt.basis[..nd], &pt.dbasis[..nd], pt.x);
        let load = interp(b, &c.load);
        let nu = interp(b, &c.nu);
        let sigma = interp(b, &c.sigma);
        let m = [interp(b, &c.m1), interp(b, &c.m2)];
        // i omega sigma <A, v>: IP%s*detJ*im*omega*C_ip*Basis(q)*Basis(p).
        let w = pt.s_detj * im * c.omega * sigma;
        for p in 0..nd {
            for q in 0..nd {
                stiff[p * nd + q] += w * b[q] * b[p];
            }
        }
        // Bt = curl of the basis (with the 1/r term), Ht = nu Bt.
        let bt: Vec<[f64; 2]> = (0..nd)
            .map(|p| {
                let mut v = [-db[p][1], db[p][0]];
                if axi {
                    v[1] += b[p] / x;
                }
                v
            })
            .collect();
        for p in 0..nd {
            let ht = [nu * bt[p][0], nu * bt[p][1]];
            for q in 0..nd {
                let k = ht[0] * bt[q][0] + ht[1] * bt[q][1];
                stiff[p * nd + q] += Complex64::new(pt.s_detj * k, 0.0);
            }
        }
        for p in 0..nd {
            let f = if axi {
                load * b[p] - m[0] * db[p][1] + m[1] * (db[p][0] + b[p] / x)
            } else {
                load * b[p] + m[0] * db[p][1] - m[1] * db[p][0]
            };
            force[p] += pt.s_detj * f;
        }
    }
    Ok((stiff, force))
}

impl Solver for MagnetoDynamics2DHarmonic {
    /// `MagnetoDynamics2DHarmonic_Init0/_Init`: a complex two-component
    /// variable, `Potential[Potential re:1 Potential im:1]` unless the
    /// section names one, and `Linear System Complex`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let ix = slot_index(ctx);
        if ctx.sim.slots[ix].variable.is_none() {
            let mesh = ctx.mesh;
            let sim = &mut *ctx.sim;
            let slot = &sim.slots[ix];
            let perm = match &sim.backends.assembly {
                Some(b) => b.build_perm_for(mesh, &slot.active, &slot.params)?,
                None => node_order_perm(mesh, &slot.active),
            };
            let perm: Arc<[u32]> = perm.perm.into();
            let n = perm.iter().filter(|&&s| s != NO_DOF).count();
            let v = Variable::new("potential", FieldKind::Nodal, 2, Some(perm.clone()), vec![0.0; 2 * n]);
            let comps: Vec<String> = COMPONENTS.iter().map(|s| s.to_string()).collect();
            let id = sim.vars.add_with_components(v, &comps)?;
            let slot = &mut sim.slots[ix];
            slot.variable = Some(id);
            slot.perm = Some(perm);
        } else if ctx.primary().is_some_and(|v| v.dofs != 2) {
            return Err(SolverError::fatal(
                CALLER,
                "Variable is not properly defined for time harmonic solver, use Potential[Potential re:1 Potential im:1]",
            ));
        }
        let params = &mut ctx.sim.slots[ix].params;
        if params.get("linear system complex").is_none() {
            params.insert("Linear System Complex", Value::Logical(true));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("a transient simulation"));
        }
        check_supported(ctx)?;
        let axi = match ctx.simulation_params().string("coordinate system").map(elmer_core::normalize_key) {
            None => false,
            Some(s) if s == "cartesian" || s == "cartesian 2d" => false,
            Some(s) if s == "axi symmetric" || s == "cylindric symmetric" => true,
            Some(s) => return Err(unsupported(format!("Coordinate System = {s}"))),
        };
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("a mesh that is not 2D"));
        }
        let mu0 = ctx.constants().real("permeability of vacuum").unwrap_or(std::f64::consts::PI * 4.0e-7);
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let maxiter = ctx.params().integer("nonlinear system max iterations").map_or(1, |n| n.max(1));
        let mut local = LocalSystem::zeros(0);
        let mut dofs = Vec::new();
        for _ in 0..maxiter {
            ctx.default_initialize()?;
            let mut omega = 0.0;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let bp = ctx
                    .body_params(e)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let (material, body_force) = (
                    bp.material
                        .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no material", e.0)))?,
                    bp.body_force,
                );
                omega = angular_frequency(ctx, e)?;
                let zero = || vec![0.0; n];
                let at = |p: Option<&ParamMap>, k: &str| -> Result<Vec<f64>, SolverError> {
                    Ok(match p {
                        Some(p) => ctx.real_at_nodes(p, k, e)?.unwrap_or_else(zero),
                        None => zero(),
                    })
                };
                let c = Bulk {
                    nu: reluctivity(ctx, material, e, mu0)?,
                    sigma: at(Some(material), "electric conductivity")?,
                    m1: at(Some(material), "magnetization 1")?,
                    m2: at(Some(material), "magnetization 2")?,
                    load: at(body_force, "current density")?,
                    omega,
                };
                if body_force.is_some_and(|b| b.get("current density im").is_some()) {
                    return Err(unsupported("Current Density im"));
                }
                let (s, f) = bulk_local(etype, &xyz, axi, &c)?;
                interleave(&s, &f, n, &mut local);
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let Some(bc) = ctx.bc_params(e) else { continue };
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let nd = xyz.len();
                let mut stiff = vec![Complex64::new(0.0, 0.0); nd * nd];
                let mut force = vec![Complex64::new(0.0, 0.0); nd];
                if let Some(cond) = ctx.real_at_nodes(bc, "layer electric conductivity", e)? {
                    // LocalMatrixSkinBC (:2675), with the last bulk element's omega.
                    let mu = ctx.real_at_nodes(bc, "layer relative permeability", e)?.unwrap_or_else(|| vec![0.0; nd]);
                    let mu_vacuum = 4.0 * std::f64::consts::PI * 1e-7;
                    let imu = Complex64::new(0.0, 1.0);
                    for pt in points(etype, &xyz, axi)? {
                        let b = &pt.basis[..nd];
                        let mu_ip = mu_vacuum * interp(b, &mu);
                        let cond_ip = interp(b, &cond);
                        let delta = (2.0 / (cond_ip * omega * mu_ip)).sqrt();
                        let inv_zs = Complex64::new(cond_ip * delta, 0.0) / (1.0 + imu);
                        for p in 0..nd {
                            for q in 0..nd {
                                stiff[p * nd + q] += pt.s_detj * (imu * omega * inv_zs) * b[p] * b[q];
                            }
                        }
                    }
                } else {
                    // LocalMatrixBC (:2607): air gap and surface current.
                    let air_gap = bc.logical("air gap") == Some(true);
                    let sc = ctx.real_at_nodes(bc, "surface current", e)?;
                    let sc_im = ctx.real_at_nodes(bc, "surface current im", e)?;
                    if !(air_gap || sc.is_some() || sc_im.is_some()) {
                        continue;
                    }
                    let len = if air_gap {
                        ctx.real_at_nodes(bc, "air gap length", e)?
                            .ok_or_else(|| SolverError::fatal("LocalMatrixBC", "\"Air Gap Length\" not found!"))?
                    } else {
                        vec![0.0; nd]
                    };
                    let gap_mu = ctx.real_at_nodes(bc, "air gap relative permeability", e)?.unwrap_or_else(|| vec![1.0; nd]);
                    let sc = sc.unwrap_or_else(|| vec![0.0; nd]);
                    let sc_im = sc_im.unwrap_or_else(|| vec![0.0; nd]);
                    for pt in points(etype, &xyz, axi)? {
                        let (b, db) = (&pt.basis[..nd], &pt.dbasis[..nd]);
                        if air_gap {
                            let mu = 4.0 * std::f64::consts::PI * 1e-7 * interp(b, &gap_mu);
                            let l = interp(b, &len);
                            for p in 0..nd {
                                for q in 0..nd {
                                    let g = db[p][0] * db[q][0] + db[p][1] * db[q][1] + db[p][2] * db[q][2];
                                    stiff[p * nd + q] += Complex64::new(pt.s_detj * (l / mu * g), 0.0);
                                }
                            }
                        }
                        let surf = Complex64::new(interp(b, &sc), interp(b, &sc_im));
                        for p in 0..nd {
                            force[p] += pt.s_detj * surf * b[p];
                        }
                    }
                }
                interleave(&stiff, &force, nd, &mut local);
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;
            apply_dirichlet(ctx)?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}

/// `DefaultDirichletBCs` for the two components by their names
/// (`ComponentName(name, DOF)`): `Potential re` on the even rows,
/// `Potential im` on the odd ones.
fn apply_dirichlet(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    let ix = slot_index(ctx);
    let var = ctx.sim.slots[ix].variable.ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
    let names: Vec<String> = {
        let v = ctx.sim.vars.by_id(var);
        let mut comps: Vec<(usize, String)> = ctx
            .sim
            .vars
            .names_in_order()
            .into_iter()
            .filter_map(|(name, h)| match h {
                elmer_solver::VarHandle::Component { parent, component } if parent == var => Some((component, name)),
                _ => None,
            })
            .collect();
        comps.sort();
        let comps: Vec<String> = comps.into_iter().map(|(_, n)| n).collect();
        if comps.len() == v.dofs {
            comps
        } else {
            (1..=v.dofs).map(|k| format!("{} {k}", v.name)).collect()
        }
    };
    let perm = ctx.sim.slots[ix].perm.clone().ok_or_else(|| SolverError::fatal(CALLER, "no perm"))?;
    let mut fixed = std::collections::BTreeMap::new();
    let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
    for e in bnd {
        let Some(bc) = ctx.bc_params(e) else { continue };
        for (k, name) in names.iter().enumerate() {
            if bc.get(&format!("{name} condition")).is_some() {
                return Err(unsupported(format!("Dirichlet mask '{name} Condition'")));
            }
            let Some(vals) = ctx.real_at_nodes(bc, name, e)? else { continue };
            for (j, &node) in ctx.mesh.element(e).nodes.iter().enumerate() {
                let s = perm[node as usize];
                if s != NO_DOF {
                    fixed.insert(2 * s + k as u32, vals[j]);
                }
            }
        }
    }
    let fixed: Vec<(u32, f64)> = fixed.into_iter().collect();
    let sim = &mut *ctx.sim;
    let backend = sim
        .backends
        .assembly
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
    let slot = &mut sim.slots[ix];
    let sys = slot
        .system
        .as_mut()
        .ok_or_else(|| SolverError::fatal("DefaultDirichletBCs", "call default_initialize first"))?;
    backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, &fixed, &slot.params)?;
    Ok(())
}
