//! `StatMagSolver` (fem/src/modules/StatMagSolve.F90) in 2D: the magnetic
//! vector potential with the axisymmetric operator Elmer composes
//! (`StatMagAxisCompose`, `HarmMagAxisCompose`: the `r` weight and the
//! `1/r` terms whatever the coordinate system), steady (one DOF) or
//! time-harmonic (two, `Harmonic Simulation`), and the Joule heating fields
//! (`Joule Field = omega^2/2 |A|^2`, `Joule Heating = max sigma * Joule
//! Field` over the elements at a node), scaled to `Desired Heating Power`
//! (or `Power Control`) by the integrated `2 pi r sigma Joule Field`.
//!
//! [`init_section`] is `StatMagSolver_Init` (`Variable DOFs` and the
//! exported fields). Refused: 3D (`StatMagCartesianCompose`) and `Calculate
//! Magnetic Flux` (`AxiSCurl`).

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value};
use elmer_fe::metric::MetricStatus;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

const CALLER: &str = "StatMagSolve";

pub const PROCEDURE: (&str, &str) = ("StatMagSolve", "StatMagSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StatMag))
}

/// Whether `params` runs `StatMagSolver`.
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

fn next_exported(p: &ParamMap) -> String {
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    format!("Exported Variable {k}")
}

/// `StatMagSolver_Init`.
pub fn init_section(p: &mut ParamMap) {
    if p.logical("harmonic simulation") == Some(true) {
        p.insert("Variable Dofs", Value::Integer(vec![2]));
    }
    let flux = p.logical("calculate magnetic flux") == Some(true);
    let flux_abs = p.logical("calculate magnetic flux abs") == Some(true);
    if flux {
        let k = next_exported(p);
        p.insert(&k, Value::Str("-dofs 3 Magnetic Flux Density".into()));
        if flux_abs {
            let k = next_exported(p);
            p.insert(&k, Value::Str("Magnetic Flux Density_abs".into()));
        }
    }
    if p.logical("calculate joule heating") == Some(true) {
        let k = next_exported(p);
        p.insert(&k, Value::Str("Joule Heating".into()));
        let k = next_exported(p);
        p.insert(&k, Value::Str("Joule Field".into()));
    }
}

struct StatMag;

impl Solver for StatMag {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        // `Solver % Matrix % COMPLEX = .TRUE.` for a harmonic simulation:
        // complex scaling and solution of the interleaved system.
        let dofs = ctx.primary().map_or(0, |v| v.dofs);
        let harmonic = ctx.params().logical("harmonic simulation").unwrap_or(dofs == 2);
        if harmonic && ctx.mesh.dim() < 3 {
            let id = ctx.solver_id();
            let ix = ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot");
            ctx.sim.slots[ix].params.insert("Linear System Complex", Value::Logical(true));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        self.run(ctx, dt, transient)
    }
}

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("StatMagSolver: {}", what.into()))
}

/// `SUM(a(1:n) * b(1:n))`.
fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `ListGetReal(section, key, n, NodeIndexes)` (zero when absent).
fn at_nodes(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Option<Vec<f64>>, SolverError> {
    match p {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
    .map(|v| v.map(|mut v| {
        v.resize(n, 0.0);
        v
    }))
}

/// One Gauss point: `s = SqrtElementMetric * S_Integ(t)`, `r = SUM(Basis *
/// x)`, the basis and its x/y derivatives.
struct Point {
    s: f64,
    r: f64,
    basis: Vec<f64>,
    dx: Vec<f64>,
    dy: Vec<f64>,
}

fn points(ctx: &SolverCtx<'_>, e: ElemId, status: &mut MetricStatus) -> Result<Vec<Point>, SolverError> {
    let el = ctx.mesh.element(e);
    let n = el.nodes.len();
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let xs: Vec<f64> = xyz.iter().map(|c| c[0]).collect();
    let def = elmer_fe::element_type(el.etype.0).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let mut out = Vec::with_capacity(rule.len());
    for t in 0..rule.len() {
        let mut be = elmer_core::BasisEval::default();
        elmer_fe::info::element_info(def, &xyz, [rule.u[t], rule.v[t], rule.w[t]], 2, status, &mut be)
            .map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
        let basis = be.basis[..n].to_vec();
        let r = sum(&basis, &xs);
        out.push(Point {
            s: be.det_j * rule.s[t],
            r,
            dx: be.dbasis_dx[..n].iter().map(|d| d[0]).collect(),
            dy: be.dbasis_dx[..n].iter().map(|d| d[1]).collect(),
            basis,
        });
    }
    Ok(out)
}

/// `ListGetAngularFrequency()` without a value list or element: the
/// simulation, the solver, then the equation of the first active element.
fn angular_frequency(ctx: &SolverCtx<'_>) -> Result<f64, SolverError> {
    let two_pi = 2.0 * std::f64::consts::PI;
    let first = ctx.active_elements().next();
    let eq = first.and_then(|e| ctx.body_params(e)).and_then(|b| b.equation);
    for p in [Some(ctx.simulation_params()), Some(ctx.params()), eq].into_iter().flatten() {
        if let Some(f) = ctx.creal(p, "frequency")? {
            return Ok(two_pi * f);
        }
        if let Some(w) = ctx.creal(p, "angular frequency")? {
            return Ok(w);
        }
    }
    Ok(0.0)
}

impl StatMag {
    fn run(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("a transient simulation"));
        }
        ctx.info(CALLER, " ", 4);
        ctx.info(CALLER, "------------------------------------------------", 4);
        ctx.info(CALLER, "Starting simulation", 4);
        let dim = ctx.mesh.dim();
        let mu0 = ctx.constants().real("permeability of vacuum").unwrap_or(std::f64::consts::PI * 4.0e-7);
        let dofs = ctx.primary().map(|v| v.dofs).ok_or_else(|| SolverError::fatal(CALLER, "the solver has no variable"))?;
        let p = ctx.params().clone();
        let harmonic = p.logical("harmonic simulation").unwrap_or(dofs == 2);
        if harmonic {
            ctx.info("StatMagSolver", "Assuming harmonic simulation", 4);
        } else {
            ctx.info("StatMagSolver", "Assuming steady-state simulation", 4);
        }
        if dim == 3 {
            return Err(unsupported("3D (StatMagCartesianCompose)"));
        }
        let mut omega = 0.0;
        if harmonic {
            if dofs != 2 {
                return Err(SolverError::fatal(
                    "StatMagSolver",
                    "For harmonic simulation there must be two components for the vector potential!",
                ));
            }
            omega = angular_frequency(ctx)?;
        } else if dofs != 1 {
            return Err(SolverError::fatal(
                "StatMagSolver",
                "For 2d steady-state cases there should be just one component for vector potential!",
            ));
        }
        if p.logical("calculate magnetic flux") == Some(true) {
            return Err(unsupported("'Calculate Magnetic Flux' (AxiSCurl)"));
        }
        let mut joule = p.logical("calculate joule heating") == Some(true);
        if joule && !harmonic {
            ctx.warn("StatMagSolver", "Computation of Joule heating is only relevan for harmonic systems!");
            joule = false;
        }
        let nb = ctx.mesh.n_bulk();
        let elems: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nb).collect();
        let max_iter = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let mut status = MetricStatus::default();
        for iter in 1..=max_iter {
            ctx.info(CALLER, &format!("Magnetic Field Iteration: {iter}"), 4);
            ctx.default_initialize()?;
            for &e in &elems {
                let n = ctx.mesh.element(e).nodes.len();
                let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal(CALLER, "element without a body"))?;
                let material = bp.material.ok_or_else(|| SolverError::fatal(CALLER, "body without a material"))?;
                let body_force = bp.body_force;
                let perm = match at_nodes(ctx, Some(material), "relative permeability", e, n)? {
                    Some(mr) => mr.iter().map(|m| mu0 * m).collect(),
                    None => at_nodes(ctx, Some(material), "magnetic permeability", e, n)?.unwrap_or_else(|| vec![mu0; n]),
                };
                let reluct: Vec<f64> = perm.iter().map(|m| 1.0 / m).collect();
                let mut cond = vec![0.0; n];
                if harmonic {
                    cond = match at_nodes(ctx, Some(material), "electrical conductivity", e, n)? {
                        Some(c) => {
                            ctx.warn(CALLER, "Use > Electric Conductivity < instead of electrial");
                            c
                        }
                        None => at_nodes(ctx, Some(material), "electric conductivity", e, n)?.unwrap_or_else(|| vec![0.0; n]),
                    };
                }
                let load = at_nodes(ctx, body_force, "current density", e, n)?.unwrap_or_else(|| vec![0.0; n]);
                let angle = if harmonic {
                    at_nodes(ctx, body_force, "current phase angle", e, n)?.unwrap_or_else(|| vec![0.0; n])
                } else {
                    vec![0.0; n]
                };
                let pts = points(ctx, e, &mut status)?;
                let nd = dofs * n;
                let mut stiff = vec![0.0; nd * nd];
                let mut force = vec![0.0; nd];
                for pt in &pts {
                    let (b, dx, dy, r, s) = (&pt.basis, &pt.dx, &pt.dy, pt.r, pt.s);
                    let f = sum(&load, b);
                    let rel = sum(&reluct, b);
                    if harmonic {
                        // HarmMagAxisCompose.
                        let c = sum(&cond, b);
                        let ang = (std::f64::consts::PI / 180.0) * sum(&angle, b);
                        for pp in 0..n {
                            for q in 0..n {
                                let mut a11 =
                                    dx[pp] * dx[q] * r + dy[pp] * dy[q] * r + b[pp] * dx[q] + b[q] * dx[pp] + b[pp] * b[q] / r;
                                a11 *= rel * s; // Reluctivity * s * a11
                                let a22 = a11;
                                let a21 = -c * omega * s * r * b[q] * b[pp];
                                let a12 = -a21;
                                stiff[(2 * pp) * nd + 2 * q] += a11;
                                stiff[(2 * pp + 1) * nd + 2 * q + 1] += a22;
                                stiff[(2 * pp + 1) * nd + 2 * q] += a21;
                                stiff[(2 * pp) * nd + 2 * q + 1] += a12;
                            }
                        }
                        for pp in 0..n {
                            force[2 * pp] += f * ang.cos() * b[pp] * r * s;
                            force[2 * pp + 1] += f * ang.sin() * b[pp] * r * s;
                        }
                    } else {
                        // StatMagAxisCompose.
                        for pp in 0..n {
                            for q in 0..n {
                                let mut mat =
                                    r * dx[pp] * dx[q] + r * dy[pp] * dy[q] + b[pp] * dx[q] + b[q] * dx[pp] + b[pp] * b[q] / r;
                                mat = mat * rel * s;
                                stiff[pp * nd + q] += mat;
                            }
                        }
                        for pp in 0..n {
                            force[pp] += f * b[pp] * r * s;
                        }
                    }
                }
                let ls = LocalSystem { n: nd, stiff, mass: None, damp: None, force };
                ctx.default_update(&ls, e)?;
            }
            ctx.info(CALLER, "Assembly done", 4);
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.info(CALLER, "Set boundaries done", 4);
            ctx.default_solve()?;
            if ctx.default_converged() {
                ctx.info(CALLER, &format!("Convergence after {iter} iterations"), 4);
                break;
            }
        }

        if harmonic && joule {
            self.joule_heating(ctx, &elems, omega, &p, &mut status)?;
        }
        ctx.info(CALLER, "All done, exiting", 4);
        ctx.info(CALLER, "------------------------------------------------", 4);
        Ok(())
    }
}

impl StatMag {
    fn joule_heating(
        &self,
        ctx: &mut SolverCtx<'_>,
        elems: &[ElemId],
        omega: f64,
        p: &ParamMap,
        status: &mut MetricStatus,
    ) -> Result<(), SolverError> {
        let (mvp, perm) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "the solver has no variable"))?;
            (v.values.to_vec(), v.perm.map(|p| p.to_vec()))
        };
        let perm = perm.ok_or_else(|| unsupported("a potential without a permutation"))?;
        let slot = |node: u32| perm.get(node as usize).copied().filter(|&s| s != elmer_core::NO_DOF).map(|s| s as usize);
        let mut abs_joule = ctx
            .variable("joule field")
            .map(|f| f.values.to_vec())
            .ok_or_else(|| SolverError::fatal(CALLER, "no Joule Field variable"))?;
        let mut heat = ctx
            .variable("joule heating")
            .map(|f| f.values.to_vec())
            .ok_or_else(|| SolverError::fatal(CALLER, "no Joule Heating variable"))?;
        let jc = 0.5 * (omega * omega);
        for node in 0..ctx.mesh.n_nodes() as u32 {
            if let Some(j) = slot(node) {
                let (jre, jim) = (mvp[2 * j], mvp[2 * j + 1]);
                abs_joule[j] = jc * (jre * jre + jim * jim);
            }
        }
        let mut total_heating = 0.0;
        let mut total_volume = 0.0;
        let two_pi = 2.0 * std::f64::consts::PI;
        for &e in elems {
            let nodes = ctx.mesh.element(e).nodes.to_vec();
            let n = nodes.len();
            let material = ctx
                .body_params(e)
                .and_then(|b| b.material)
                .ok_or_else(|| SolverError::fatal(CALLER, "body without a material"))?;
            let cond = match at_nodes(ctx, Some(material), "electrical conductivity", e, n)? {
                Some(c) => {
                    ctx.warn(CALLER, "Use electric conductivity instead of electrical");
                    c
                }
                None => at_nodes(ctx, Some(material), "electric conductivity", e, n)?.unwrap_or_else(|| vec![0.0; n]),
            };
            let ae: Vec<f64> = nodes.iter().map(|&k| slot(k).map_or(0.0, |j| abs_joule[j])).collect();
            // JouleIntegrate.
            for pt in points(ctx, e, status)? {
                let field = sum(&ae, &pt.basis);
                let c = sum(&cond, &pt.basis);
                for pp in 0..n {
                    total_volume += two_pi * pt.r * pt.s * pt.basis[pp];
                    total_heating += two_pi * pt.r * pt.s * pt.basis[pp] * field * c;
                }
            }
            for (i, &k) in nodes.iter().enumerate() {
                if let Some(j) = slot(k) {
                    let v = cond[i] * abs_joule[j];
                    if v > heat[j] {
                        heat[j] = v;
                    }
                }
            }
        }
        let _ = total_volume;
        let desired = match ctx.creal(p, "desired heating power")? {
            Some(d) => Some(d),
            None => ctx.creal(p, "power control")?,
        };
        if let Some(d) = desired.filter(|_| total_heating > 0.0) {
            let f = d / total_heating;
            abs_joule.iter_mut().for_each(|x| *x *= f);
            heat.iter_mut().for_each(|x| *x *= f);
        }
        ctx.info(
            "StatMagSolve",
            &format!("Joule Heating (W): {}", elmer_solver::fortran_fmt::es(total_heating, 15, 4)),
            4,
        );
        for (name, vals) in [("joule field", abs_joule), ("joule heating", heat)] {
            if let Some(f) = ctx.variable_mut(name) {
                f.values.copy_from_slice(&vals);
            }
            if let Some(id) = ctx.sim.vars.id(name) {
                ctx.sim.vars.by_id_mut(id).values_changed = true;
            }
        }
        Ok(())
    }
}
