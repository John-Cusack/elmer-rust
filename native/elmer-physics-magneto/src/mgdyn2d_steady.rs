//! `MagnetoDynamics2D` (MagnetoDynamics2D.F90:122), steady state: the real
//! out-of-plane (or azimuthal) vector potential on nodal elements,
//!
//! ```text
//!   -div(nu grad A) = J + curl M      (+ A/r terms in axial symmetry)
//! ```
//!
//! with `LocalMatrix` for linear materials (`GetReluctivity`: `Relative
//! Permeability`, `Permeability` or `Reluctivity`), the current density and
//! magnetisation loads, `LocalMatrixInfinityBC` (`Infinity BC`, with the
//! oriented `NormalVector`), `LocalMatrixBC` (air gap, surface current),
//! `DefaultDirichletBCs` and the nonlinear iteration loop, and the
//! `Coordinates` exported variable filled at the end.
//!
//! [`init_section`] is `MagnetoDynamics2D_Init` (one DOF, default variable
//! `Potential`). Refused: transient runs, H-B curves and Zirka materials,
//! London equations, velocities, coils, handle assembly, mortar BCs,
//! `Magnetic Flux Density` BCs, torque control and the lumped quantities of
//! `CalculateLumpedTransient` (`R Inner`, `Rotor Radius`, `Calculate
//! Potential`).

use elmer_core::{ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::mgdyn2d::{interp, points, reluctivity};

const CALLER: &str = "MagnetoDynamics2D";

pub const PROCEDURE: (&str, &str) = ("MagnetoDynamics2D", "MagnetoDynamics2D");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MagnetoDynamics2DSteady))
}

/// Whether `params` runs `MagnetoDynamics2D` (the steady/transient solver).
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `MagnetoDynamics2D_Init`: `Variable Dofs = 1`, `Variable = Potential`
/// unless given, `Apply Mortar BCs` and `Use Global Mass Matrix` unless given.
pub fn init_section(p: &mut ParamMap) {
    p.insert("Variable Dofs", Value::Integer(vec![1]));
    for (k, v) in [
        ("variable", Value::Str("Potential".into())),
        ("apply mortar bcs", Value::Logical(true)),
        ("use global mass matrix", Value::Logical(true)),
    ] {
        if p.get(k).is_none() {
            p.insert(k, v);
        }
    }
}

struct MagnetoDynamics2DSteady;

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", what.into()))
}

fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

fn check_supported(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for k in [
        "electrodynamics model",
        "handle assembly",
        "skip degenerate elements",
        "newton-raphson iteration",
        "constant mass matrix",
    ] {
        if p.logical(k) == Some(true) {
            return Err(unsupported(format!("solver option '{k}'")));
        }
    }
    if p.get("nonlinear system torque tolerance").is_some() {
        return Err(unsupported("'Nonlinear System Torque Tolerance'"));
    }
    let m = &ctx.sim.model;
    for mat in &m.materials {
        for k in ["h-b curve", "london lambda", "zirka material", "reluctivity function"] {
            if mat.get(k).is_some() {
                return Err(unsupported(format!("material '{k}'")));
            }
        }
    }
    for bf in &m.body_forces {
        for k in ["lorentz velocity", "angular velocity"] {
            if has_prefix(bf, k) {
                return Err(unsupported(format!("body force '{k}'")));
            }
        }
        if bf.logical("calculate potential") == Some(true) {
            return Err(unsupported("'Calculate Potential' (CalculateLumpedTransient)"));
        }
    }
    if m.bodies.iter().any(|b| b.get("r inner").is_some()) || m.simulation.get("rotor radius").is_some() {
        return Err(unsupported("torque (CalculateLumpedTransient)"));
    }
    if m.components.iter().any(|c| c.get("coil type").is_some()) {
        return Err(unsupported("coil components"));
    }
    for bc in &m.boundary_conditions {
        if has_prefix(bc, "magnetic flux density") {
            return Err(unsupported("Magnetic Flux Density BC (SetMagneticFluxDensityBC)"));
        }
        if bc.get("mortar bc").is_some() || bc.get("periodic bc").is_some() || bc.get("anti periodic bc").is_some() {
            return Err(unsupported("mortar / periodic BCs"));
        }
    }
    Ok(())
}

/// `NormalVector(Element, Nodes, u, v, .TRUE.)` of a linear boundary line
/// (elmer-physics-flow `boundary::outward_normal`): the unit normal of the
/// segment, turned away from the parent element's node centre.
fn outward_normal_line(xyz: &[[f64; 3]], parent: &[[f64; 3]]) -> [f64; 3] {
    let du: [f64; 3] = std::array::from_fn(|a| (xyz[1][a] - xyz[0][a]) / 2.0);
    let nrm = [-du[1], du[0], 0.0];
    let len = nrm.iter().map(|x| x * x).sum::<f64>().sqrt();
    if len <= 0.0 {
        return [0.0; 3];
    }
    let nrm = nrm.map(|x| x / len);
    let center = |p: &[[f64; 3]]| -> [f64; 3] {
        let n = p.len() as f64;
        std::array::from_fn(|a| p.iter().map(|x| x[a]).sum::<f64>() / n)
    };
    let (c, pc) = (center(xyz), center(parent));
    let d: f64 = (0..3).map(|a| (pc[a] - c[a]) * nrm[a]).sum();
    if d > 0.0 {
        nrm.map(|x| -x)
    } else {
        nrm
    }
}

impl Solver for MagnetoDynamics2DSteady {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("a transient simulation"));
        }
        check_supported(ctx)?;
        if ctx.primary().map(|v| v.dofs) != Some(1) {
            return Err(SolverError::fatal(CALLER, "the solver variable must have one component"));
        }
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
        let nb = ctx.mesh.n_bulk();
        let elems: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nb).collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let max_iter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let mut dofs = Vec::new();
        for _ in 0..max_iter {
            ctx.default_initialize()?;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let n = xyz.len();
                let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let material = bp.material.ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no material", e.0)))?;
                let body_force = bp.body_force;
                let zero = || vec![0.0; n];
                let at = |p: Option<&ParamMap>, k: &str| -> Result<Vec<f64>, SolverError> {
                    Ok(match p {
                        Some(p) => ctx.real_at_nodes(p, k, e)?.unwrap_or_else(zero),
                        None => zero(),
                    })
                };
                let nu = reluctivity(ctx, material, e, mu0)?;
                let m1 = at(Some(material), "magnetization 1")?;
                let m2 = at(Some(material), "magnetization 2")?;
                let load = at(body_force, "current density")?;
                let mut stiff = vec![0.0; n * n];
                let mut force = vec![0.0; n];
                for pt in points(etype, &xyz, axi)? {
                    let (b, db, x) = (&pt.basis[..n], &pt.dbasis[..n], pt.x);
                    let load_ip = interp(b, &load);
                    let nu_ip = interp(b, &nu);
                    let m_ip = [interp(b, &m1), interp(b, &m2)];
                    // Bt = curl of the basis, Ht = nu_tensor Bt (diagonal).
                    let bt: Vec<[f64; 2]> = (0..n)
                        .map(|p| {
                            let mut v = [db[p][1], -db[p][0]];
                            if axi {
                                v = [-v[0], -v[1]];
                                v[1] += b[p] / x;
                            }
                            v
                        })
                        .collect();
                    for p in 0..n {
                        let ht = [nu_ip * bt[p][0] + 0.0 * bt[p][1], 0.0 * bt[p][0] + nu_ip * bt[p][1]];
                        for q in 0..n {
                            // MATMUL(Ht, TRANSPOSE(Bt))(p, q).
                            let k = ht[0] * bt[q][0] + ht[1] * bt[q][1];
                            stiff[p * n + q] += pt.s_detj * k;
                        }
                    }
                    for p in 0..n {
                        let f = if axi {
                            load_ip * b[p] - m_ip[0] * db[p][1] + m_ip[1] * (db[p][0] + b[p] / x)
                        } else {
                            load_ip * b[p] + m_ip[0] * db[p][1] - m_ip[1] * db[p][0]
                        };
                        force[p] += pt.s_detj * f;
                    }
                }
                let ls = LocalSystem { n, stiff, mass: None, damp: None, force };
                ctx.default_update(&ls, e)?;
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
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let n = xyz.len();
                let mut stiff = vec![0.0; n * n];
                let mut force = vec![0.0; n];
                if bc.logical("infinity bc") == Some(true) {
                    // LocalMatrixInfinityBC.
                    let parent = el.parents.iter().flatten().find(|p| (p.0 as usize) < nb).copied();
                    let Some(parent) = parent else { continue };
                    if etype != ElementType(202) {
                        return Err(unsupported(format!("Infinity BC on element type {}", etype.0)));
                    }
                    let pmat = ctx
                        .body_params(parent)
                        .and_then(|b| b.material)
                        .ok_or_else(|| SolverError::fatal(CALLER, "Infinity BC parent without a material"))?;
                    // GetReluctivity(Material, R, n, Parent): the first n of
                    // the parent's nodal values.
                    let r = reluctivity(ctx, pmat, parent, mu0)?;
                    let pxyz: Vec<[f64; 3]> = ctx.mesh.element(parent).nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                    let normal = outward_normal_line(&xyz, &pxyz);
                    for pt in points(etype, &xyz, axi)? {
                        let b = &pt.basis[..n];
                        let mu = interp(b, &r[..n]);
                        let coord: [f64; 3] = std::array::from_fn(|a| interp(b, &xyz.iter().map(|c| c[a]).collect::<Vec<_>>()));
                        let mut cn = 0.0;
                        let mut cc = 0.0;
                        for a in 0..3 {
                            cn += coord[a] * normal[a];
                        }
                        for c in &coord {
                            cc += c * c;
                        }
                        let inf_ip = mu * cn / cc;
                        for p in 0..n {
                            for q in 0..n {
                                stiff[p * n + q] += pt.s_detj * inf_ip * b[q] * b[p];
                            }
                        }
                    }
                } else {
                    // LocalMatrixBC: air gap and surface current.
                    let air_gap = bc.logical("air gap") == Some(true);
                    let sc = ctx.real_at_nodes(bc, "surface current", e)?;
                    if !(air_gap || sc.is_some()) {
                        continue;
                    }
                    let (len, gap_mu) = if air_gap {
                        let len = ctx
                            .real_at_nodes(bc, "air gap length", e)?
                            .ok_or_else(|| SolverError::fatal("LocalMatrixBC", "\"Air Gap Length\" not found!"))?;
                        (len, ctx.real_at_nodes(bc, "air gap relative permeability", e)?.unwrap_or_else(|| vec![1.0; n]))
                    } else {
                        (vec![0.0; n], vec![1.0; n])
                    };
                    for pt in points(etype, &xyz, axi)? {
                        let (b, db) = (&pt.basis[..n], &pt.dbasis[..n]);
                        if air_gap {
                            let mu = 4.0 * std::f64::consts::PI * 1e-7 * interp(b, &gap_mu);
                            let l = interp(b, &len);
                            for p in 0..n {
                                for q in 0..n {
                                    let g = db[p][0] * db[q][0] + db[p][1] * db[q][1] + db[p][2] * db[q][2];
                                    stiff[p * n + q] += pt.s_detj * l / mu * g;
                                }
                            }
                        }
                        if let Some(sc) = &sc {
                            let s = interp(b, sc);
                            for p in 0..n {
                                force[p] += pt.s_detj * s * b[p];
                            }
                        }
                    }
                }
                let ls = LocalSystem { n, stiff, mass: None, damp: None, force };
                ctx.default_update(&ls, e)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        // The `Coordinates` exported variable, if any.
        if let Some(id) = ctx.sim.vars.id("coordinates") {
            let coords: Vec<[f64; 3]> = (0..ctx.mesh.n_nodes() as u32).map(|k| ctx.mesh.coord(k)).collect();
            let v = ctx.sim.vars.by_id_mut(id);
            if v.dofs == 3 {
                if let Some(perm) = v.perm.clone() {
                    for (i, c) in coords.iter().enumerate() {
                        let s = perm[i];
                        if s != NO_DOF {
                            v.values[3 * s as usize..3 * s as usize + 3].copy_from_slice(c);
                        }
                    }
                    v.values_changed = true;
                }
            }
        }
        Ok(())
    }
}
