//! `ComputeChange` (SolverBasics.F90:3571): relaxation, the norm of the new
//! iterate, the change measure (`RELC`) and the convergence flags, for the
//! nonlinear (`NS`) and steady-state/coupled (`SS`) levels.
//!
//! Not ported yet (they return [`SolverError::Unsupported`] when asked
//! for): soft-limiter-aware relaxation, a weight variable for
//! `Set Average Solution`, and parallel residual measures. The convergence
//! monitor file is not written.

use elmer_core::{Communicator, CsrMatrix, LinearOperator, Logger, ParamMap, Params};

use crate::error::SolverError;
use crate::fortran_fmt::{g, i2s};
use crate::norm::{compute_norm, select_norm_entries, NormSpec};
use crate::variable::Variable;

const CALLER: &str = "ComputeChange";

/// Everything `ComputeChange` reads or writes besides the variable list.
pub struct ChangeEnv<'a> {
    /// The solver section.
    pub params: &'a ParamMap,
    /// The `Simulation` section (`Global Abort Not Converged`, slices).
    pub simulation: &'a ParamMap,
    /// `true` for the steady-state (coupled) level, `false` for nonlinear.
    pub steady: bool,
    /// Current `coupled iter` (steady) or `nonlin iter` (nonlinear) value,
    /// read before any advance.
    pub iter_no: i32,
    /// Entries to consider; `None` means all values.
    pub nsize: Option<usize>,
    pub matrix: Option<&'a CsrMatrix<f64>>,
    pub rhs: Option<&'a [f64]>,
    /// Mesh node count (nodal norm).
    pub n_nodes: usize,
    pub log: &'a dyn Logger,
    pub comm: &'a dyn Communicator,
}

/// Solver-level state `ComputeChange` updates (`Solver % NewtonActive`,
/// `Solver % AitkenRelax`).
#[derive(Clone, Debug, PartialEq)]
pub struct SolverState {
    pub newton_active: bool,
    pub aitken_relax: f64,
}

impl Default for SolverState {
    fn default() -> Self {
        Self {
            newton_active: false,
            aitken_relax: 1.0,
        }
    }
}

/// What one call computed.
#[derive(Clone, Debug, PartialEq)]
pub struct ChangeReport {
    pub norm: f64,
    pub change: f64,
    /// The `NS`/`SS` line as Elmer prints it (Info level 3).
    pub message: String,
    /// The section asked to skip this level's change computation.
    pub skipped: bool,
}

/// Elmer `NumericalError`: fatal unless `Global Abort Not Converged = False`.
pub fn numerical_error(
    simulation: &ParamMap,
    caller: &str,
    msg: &str,
    is_fatal: Option<bool>,
    log: &dyn Logger,
) -> Result<(), SolverError> {
    let fatal = is_fatal.unwrap_or_else(|| simulation.logical("global abort not converged").unwrap_or(true));
    let msg = format!("Numerical Error: {msg}");
    if fatal {
        Err(SolverError::numerical(caller, msg))
    } else {
        log.warn(caller, &msg);
        Ok(())
    }
}

/// One `ComputeChange` call on `var`. The caller advances `nonlin iter`
/// (unless `Skip Advance Nonlinear iter`) and sets `var.conv.nonlin_iter`,
/// as the Fortran does before reaching this point; see
/// [`crate::Simulation::compute_change`].
pub fn compute_change(env: &ChangeEnv<'_>, var: &mut Variable, st: &mut SolverState) -> Result<ChangeReport, SolverError> {
    let p = env.params;
    let skipped = || ChangeReport {
        norm: 0.0,
        change: 0.0,
        message: String::new(),
        skipped: true,
    };

    let (measure, absolute, mut relaxation, mut relax, aitken, relax_before, relative_p, skip_constraints, nodal_norm);
    // No residual mode for steady state analysis.
    let mut residual = false;
    let iter_no = env.iter_no;
    if env.steady {
        if p.logical("skip compute steady state change") == Some(true) {
            env.log.info(CALLER, "Skipping the computation of steady state change", 15);
            return Ok(skipped());
        }
        measure = p.string("steady state convergence measure").unwrap_or("norm").to_string();
        absolute = p
            .logical("steady state convergence absolute")
            .or_else(|| p.logical("use absolute norm for convergence"))
            .unwrap_or(false);
        relaxation = p.real("steady state relaxation factor").unwrap_or(1.0);
        relax = (relaxation - 1.0).abs() > f64::EPSILON;
        aitken = p.logical("steady state aitken relaxation") == Some(true);
        if aitken {
            relax = true;
        }
        if relax {
            if let Some(after) = p.integer("steady state relaxation after") {
                if after >= iter_no {
                    relax = false;
                }
            }
        }
        nodal_norm = p.logical("steady state nodal norm") == Some(true);
        relax_before = if relax {
            p.logical("steady state relaxation before").unwrap_or(true)
        } else {
            true
        };
        relative_p = false;
        skip_constraints = false;
    } else {
        if !st.newton_active {
            if let Some(i) = p.integer("nonlinear system newton after iterations") {
                if i <= iter_no {
                    st.newton_active = true;
                }
            }
        }
        if p.logical("skip compute nonlinear change") == Some(true) {
            env.log.info(CALLER, "Skipping the computation of nonlinear change", 15);
            return Ok(skipped());
        }
        residual = p.logical("linear system residual mode") == Some(true);
        measure = p.string("nonlinear system convergence measure").unwrap_or("norm").to_string();
        absolute = p
            .logical("nonlinear system convergence absolute")
            .or_else(|| p.logical("use absolute norm for convergence"))
            .unwrap_or(false);
        relaxation = p.real("nonlinear system relaxation factor").unwrap_or(1.0);
        relax = (relaxation - 1.0).abs() > f64::EPSILON;
        aitken = p.logical("nonlinear system aitken relaxation") == Some(true);
        if aitken {
            relax = true;
        }
        let mut rp = false;
        if relax {
            if let Some(after) = p.integer("nonlinear system relaxation after") {
                if after >= var.conv.nonlin_iter {
                    relax = false;
                }
            }
            rp = p.logical("relative pressure relaxation") == Some(true);
            if rp {
                env.log.info(CALLER, "Using relative pressure relaxation", 10);
            }
        }
        relative_p = rp;
        skip_constraints = p.logical("nonlinear system convergence without constraints") == Some(true);
        nodal_norm = p.logical("nonlinear system nodal norm") == Some(true);
        relax_before = if relax {
            p.logical("nonlinear system relaxation before").unwrap_or(true)
        } else {
            true
        };
    }
    let mut n = env.nsize.unwrap_or(var.values.len());
    if skip_constraints {
        if let Some(a) = env.matrix {
            n = n.min(a.nrows);
        }
    }
    if nodal_norm {
        if let Some(perm) = var.perm.as_deref() {
            let nn = env.n_nodes.min(perm.len());
            let i = perm[..nn]
                .iter()
                .map(|&s| if s == elmer_core::NO_DOF { 0 } else { s as usize + 1 })
                .max()
                .unwrap_or(0);
            n = n.min(i * var.dofs);
        }
    }

    // x0: the steady or nonlinear previous iterate; without one the Fortran
    // points x0 at x itself, which the elementwise formulas below honour.
    let mut x0: Option<Vec<f64>> = if env.steady {
        var.steady_values.clone()
    } else {
        var.nonlin_values.clone()
    };

    if aitken {
        let a = env
            .matrix
            .ok_or_else(|| SolverError::fatal(CALLER, "Aitken relaxation needs the system matrix"))?;
        let x0v = x0.clone().unwrap_or_else(|| var.values.clone());
        let mut dx = match var.delta_values.take() {
            Some(d) if d.len() == n => d,
            _ => vec![0.0; n],
        };
        aitken_relaxation(p, iter_no, &mut relaxation, st.aitken_relax, a, &var.values, &x0v, &dx, var.dofs, env.log);
        st.aitken_relax = relaxation;
        for i in 0..n {
            dx[i] = var.values[i] - x0v[i];
        }
        var.delta_values = Some(dx);
    }

    if !skip_constraints {
        if let Some(v0) = &x0 {
            if v0.len() != var.values.len() {
                env.log.info(
                    CALLER,
                    &format!(
                        "WARNING: Possible mismatch in length of vectors ({} vs. {})!",
                        i2s(var.values.len() as i64),
                        i2s(v0.len() as i64)
                    ),
                    10,
                );
                if n > v0.len() {
                    let mut y = v0.clone();
                    y.extend_from_slice(&var.values[v0.len()..n]);
                    x0 = Some(y);
                }
            }
        }
    }

    let dofs = var.dofs;
    let mut poffset = 0.0;
    if relative_p {
        let x = &var.values;
        let x0r = x0.as_deref().unwrap_or(x);
        let s1 = seq_sum((dofs - 1..n).step_by(dofs).map(|i| x[i]));
        let s0 = seq_sum((dofs - 1..n).step_by(dofs).map(|i| x0r[i]));
        let dpsum = reduce_sum(env.comm, s1 - s0);
        let nsum = reduce_sum(env.comm, n as f64 / dofs as f64);
        poffset = (1.0 - relaxation) * dpsum / nsum;
    }

    let relax_now = |x: &mut [f64], x0: Option<&[f64]>| {
        for i in 0..n {
            let a = x0.map_or(x[i], |v| v[i]);
            x[i] = (1.0 - relaxation) * a + relaxation * x[i];
        }
        if relative_p {
            for i in (dofs - 1..n).step_by(dofs) {
                x[i] += poffset;
            }
        }
    };
    if residual {
        // The solve gave the update dx of A dx = b - A x0 (SolveSystem).
        let v0 = x0
            .as_deref()
            .ok_or_else(|| SolverError::fatal(CALLER, "residual mode needs the previous iterate"))?;
        let r = if relax && relax_before { relaxation } else { 1.0 };
        for i in 0..n {
            var.values[i] = v0[i] + r * var.values[i];
        }
    } else if relax && relax_before {
        relax_now(&mut var.values, x0.as_deref());
    }

    let prev_norm = if env.steady { var.conv.prev_norm } else { var.conv.norm };

    let spec = NormSpec::from_params(p, dofs);
    let norm_of = |x: &[f64], n: usize| {
        let (vals, m) = select_norm_entries(&spec, x, n, var.perm.as_deref(), env.n_nodes, env.log);
        compute_norm(&spec, &vals, m, env.comm)
    };
    let norm = norm_of(&var.values, n);
    var.conv.norm = norm;

    if norm.is_nan() {
        env.log.info(
            CALLER,
            &format!("Norm and previous norm: {}{} size: {}", g(norm, 15, 8), g(prev_norm, 15, 8), n),
            3,
        );
        numerical_error(env.simulation, CALLER, "Norm of solution appears to be NaN", None, env.log)?;
    }

    let max_norm = if env.steady {
        p.real("steady state max norm")
    } else {
        p.real("nonlinear system max norm")
    };
    if let Some(mx) = max_norm {
        if norm > mx {
            env.log.info(
                CALLER,
                &format!("Computed norm and given max norm: {}{}", g(norm, 15, 8), g(mx, 15, 8)),
                4,
            );
            numerical_error(env.simulation, CALLER, "Norm of solution exceeded given bounds", None, env.log)?;
        }
    }

    let x0_or_x = |i: usize, x: &[f64]| x0.as_ref().map_or(x[i], |v| v[i]);
    let mut change = 0.0;
    match measure.as_str() {
        "residual" => {
            if env.comm.size() > 1 {
                return Err(SolverError::Unsupported("parallel residual convergence measure".into()));
            }
            let a = env
                .matrix
                .ok_or_else(|| SolverError::fatal(CALLER, "residual measure needs the system matrix"))?;
            let b = env
                .rhs
                .ok_or_else(|| SolverError::fatal(CALLER, "residual measure needs the system rhs"))?;
            let x0v: Vec<f64> = (0..var.values.len()).map(|i| x0_or_x(i, &var.values)).collect();
            let mut r = vec![0.0; a.nrows.max(n)];
            a.apply(&x0v, &mut r[..a.nrows]);
            for i in 0..n {
                r[i] -= b[i];
            }
            // r is ALLOCATE(r(n)): component strides end at n.
            change = norm_of(&r[..n], n);
            let bnorm = norm_of(b, n);
            if !absolute && bnorm > 0.0 {
                change /= bnorm;
            }
        }
        "linear system residual" => {
            if env.comm.size() > 1 {
                return Err(SolverError::Unsupported("parallel residual convergence measure".into()));
            }
            let a = env
                .matrix
                .ok_or_else(|| SolverError::fatal(CALLER, "residual measure needs the system matrix"))?;
            let b = env
                .rhs
                .ok_or_else(|| SolverError::fatal(CALLER, "residual measure needs the system rhs"))?;
            let mut r = vec![0.0; a.nrows.max(n)];
            a.apply(&var.values, &mut r[..a.nrows]);
            for i in 0..n {
                r[i] -= b[i];
            }
            change = dot(&r[..n], &r[..n]).sqrt();
            if !absolute {
                let bnorm = dot(&b[..n], &b[..n]).sqrt();
                if bnorm > 0.0 {
                    change /= bnorm;
                }
            }
        }
        "solution" => {
            let r: Vec<f64> = (0..n).map(|i| var.values[i] - x0_or_x(i, &var.values)).collect();
            change = norm_of(&r, n);
            if !absolute && norm + prev_norm > 0.0 {
                change = change * 2.0 / (norm + prev_norm);
            }
        }
        "norm" => {
            change = (norm - prev_norm).abs();
            if !absolute && norm + prev_norm > 0.0 {
                change = change * 2.0 / (norm + prev_norm);
            }
        }
        other => env.log.warn(CALLER, &format!("Unknown convergence measure: {other}")),
    }

    if env.simulation.integer("number of slices").unwrap_or(0) > 1 {
        env.log.info(CALLER, "Communicating maximum norm in single mesh operation", 10);
        change = reduce_max(env.comm, change);
    }

    if env.steady {
        let prev_change = var.conv.steady_change;
        var.conv.steady_change = change;
        if let Some(tol) = p.real("steady state convergence tolerance") {
            var.conv.steady_converged = i32::from(change <= tol);
        }
        if let Some(tol) = p.real("steady state divergence limit") {
            if change > tol && iter_no > 1 && change > prev_change {
                env.log.info(CALLER, "Steady state iteration diverged over tolerance", 5);
                var.conv.steady_converged = 2;
            }
        }
        if let Some(tol) = p.real("steady state exit condition") {
            if tol > 0.0 {
                env.log.info(CALLER, "Nonlinear iteration condition enforced by exit condition", 6);
                var.conv.steady_converged = 3;
            }
        }
    } else {
        let prev_change = var.conv.nonlin_change;
        var.conv.nonlin_change = change;
        var.conv.nonlin_converged = 0;
        let max_iter = p.integer("nonlinear system max iterations").unwrap_or(0);
        if let Some(tol) = p.real("nonlinear system convergence tolerance") {
            if change <= tol {
                var.conv.nonlin_converged = 1;
            } else if iter_no >= max_iter {
                if p.logical("nonlinear system abort not converged") == Some(true) {
                    return Err(SolverError::fatal(CALLER, "Nonlinear iteration did not converge to tolerance"));
                }
                env.log.info(CALLER, "Nonlinear iteration did not converge to tolerance", 6);
            }
        }
        if let Some(tol) = p.real("nonlinear system divergence limit") {
            if change > tol && ((iter_no > 1 && change > prev_change) || iter_no >= max_iter) {
                if p.logical("nonlinear system abort diverged") == Some(true) {
                    return Err(SolverError::fatal(CALLER, "Nonlinear iteration diverged over limit"));
                }
                env.log.info(CALLER, "Nonlinear iteration diverged over limit", 6);
                var.conv.nonlin_converged = 2;
            }
        }
        if let Some(tol) = p.real("nonlinear system exit condition") {
            if tol > 0.0 {
                env.log.info(CALLER, "Nonlinear iteration condition enforced by exit condition", 6);
                var.conv.nonlin_converged = 3;
            }
        }
        if var.conv.nonlin_converged > 1 {
            if let Some(min_iter) = p.integer("nonlinear system min iterations") {
                if iter_no < min_iter {
                    env.log.info(CALLER, "Enforcing continuation of iteration", 6);
                    var.conv.nonlin_converged = 0;
                }
            }
        }
        if !st.newton_active {
            if let Some(tol) = p.real("nonlinear system newton after tolerance") {
                if change < tol {
                    st.newton_active = true;
                }
            }
        }
    }

    if relax && !relax_before {
        relax_now(&mut var.values, x0.as_deref());
        var.conv.norm = norm_of(&var.values, n);
    }

    let name = p.string("equation").map(str::to_string).unwrap_or_else(|| var.name.clone());
    let message = format!(
        "{} (ITER={}) (NRM,RELC): ({}{} ) :: {}",
        if env.steady { "SS" } else { "NS" },
        i2s(iter_no as i64),
        g(norm, 15, 8),
        g(change, 15, 8),
        name.trim_end()
    );
    env.log.info(CALLER, &message, 3);

    let set_avg = if env.steady {
        p.logical("nonlinear system set average solution")
    } else {
        p.logical("linear system set average solution")
    };
    if set_avg == Some(true) {
        if env.comm.size() > 1 {
            return Err(SolverError::fatal(CALLER, "Setting average value not implemented in parallel!"));
        }
        if p.string("average solution weight variable").is_some() {
            return Err(SolverError::Unsupported("Average Solution Weight Variable".into()));
        }
        let x = &mut var.values;
        let ctarget = p.real("average solution value").unwrap_or(0.0) - seq_sum(x.iter().copied()) / x.len() as f64;
        for v in x.iter_mut() {
            *v += ctarget;
        }
    }

    Ok(ChangeReport {
        norm,
        change,
        message,
        skipped: false,
    })
}

/// `CalculateAitkenRelaxation` (contained in ComputeChange), serial: the
/// owned permutation is the identity.
#[allow(clippy::too_many_arguments)]
fn aitken_relaxation(
    p: &ParamMap,
    iter: i32,
    r: &mut f64,
    r0: f64,
    a: &CsrMatrix<f64>,
    x0: &[f64],
    x1: &[f64],
    dx: &[f64],
    dofs: usize,
    log: &dyn Logger,
) {
    if iter <= 1 {
        return;
    }
    let citer = p.integer("aitken relaxation factor after iterations").unwrap_or(0).max(1);
    let k = if dofs > 1 {
        p.integer("aitken relaxation component").unwrap_or(0)
    } else {
        0
    };
    let n = a.nrows.min(dx.len()).min(x0.len()).min(x1.len());
    let (mut aa, mut dd) = (0.0f64, 0.0f64);
    for i in 0..n {
        if k > 0 && (i % dofs) as i32 + 1 != k {
            continue;
        }
        let ds = x0[i] - x1[i] - dx[i];
        aa += dx[i] * ds;
        dd += ds * ds;
    }
    if dd < f64::MIN_POSITIVE {
        log.info(
            "CalculateAitkenRelaxation",
            "Vanishing denominator, keeping the previous relaxation factor",
            8,
        );
        return;
    }
    let mut r1 = -r0 * aa / dd;
    log.info(
        "CalculateAitkenRelaxation",
        &format!("Aitken relaxation factor suggested: {}", crate::fortran_fmt::es(r1, 12, 3)),
        10,
    );
    if iter > citer {
        if let Some(s) = p.real("aitken relaxation factor relaxation") {
            r1 = s * r1 + (1.0 - s) * r0;
        }
        if let Some(s) = p.real("aitken relaxation factor minimum") {
            r1 = s.max(r1);
        }
        if let Some(s) = p.real("aitken relaxation factor maximum") {
            r1 = s.min(r1);
        }
        if let Some(s) = p.real("aitken relaxation factor increase") {
            if r1.abs() > s * r0.abs() {
                r1 = (s * r0).abs().copysign(r1);
            }
        }
        if let Some(s) = p.real("aitken relaxation factor decrease") {
            if r1.abs() < s * r0.abs() {
                r1 = (s * r0).abs().copysign(r1);
            }
        }
        log.info(
            "CalculateAitkenRelaxation",
            &format!("Aitken relaxation factor regulated: {}", crate::fortran_fmt::es(r1, 12, 3)),
            5,
        );
        *r = r1;
    }
}

fn seq_sum(it: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for v in it {
        s += v;
    }
    s
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    seq_sum(a.iter().zip(b).map(|(x, y)| x * y))
}

fn reduce_sum(comm: &dyn Communicator, v: f64) -> f64 {
    let mut b = [v];
    comm.sum_f64(&mut b);
    b[0]
}

fn reduce_max(comm: &dyn Communicator, v: f64) -> f64 {
    let mut b = [v];
    comm.max_f64(&mut b);
    b[0]
}
