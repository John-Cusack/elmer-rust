//! `FindOptimum` (fem/src/modules/FindOptimum.F90), the secant method: the
//! solver's (global) variable holds the parameters; each visit after the
//! first reads the cost (`Cost Function` of the Simulation, else of the
//! section, else the Simulation constant named by `Cost Function Name`),
//! less `Cost Function Target`, made absolute or negated on request, keeps
//! the best parameters and sets the global `Found Better` to 1 on an
//! improvement (-1 otherwise), and moves the one free parameter:
//! first to its start value, then by `Step Size`, then by the relaxed
//! secant step `relax * f1 (x1 - x0) / (f1 - f0)` limited to `Max Step
//! Size` (nothing once `|f1| < Optimization Tolerance`). With `Optimal
//! Finish` the last visit restores the best parameters.
//!
//! Refused: methods other than `secant`, `Optimal Restart`, `Internal
//! History`, and the `Best File` and `History File` outputs.

use elmer_core::{ParamMap, Params};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

const CALLER: &str = "FindOptimum";
const PROCEDURE: (&str, &str) = ("FindOptimum", "FindOptimum");

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |p| {
        check(p)?;
        Ok(Box::new(FindOptimum::default()))
    });
}

fn check(p: &ParamMap) -> Result<(), SolverError> {
    let method = p.string("optimization method").map(elmer_core::normalize_key);
    if method.as_deref() != Some("secant") {
        return Err(SolverError::Unsupported(format!(
            "{CALLER}: Optimization Method {} (only secant is ported)",
            method.as_deref().unwrap_or("(none)")
        )));
    }
    for key in ["optimal restart", "internal history"] {
        if p.logical(key) == Some(true) {
            return Err(SolverError::Unsupported(format!("{CALLER}: {key}")));
        }
    }
    for key in ["best file", "history file"] {
        if p.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("{CALLER}: {key}")));
        }
    }
    Ok(())
}

/// `ES15.6E3`: `ES15.6` with a three-digit exponent.
fn es_e3(x: f64) -> String {
    let s = elmer_solver::fortran_fmt::es(x, 15, 6);
    match s.rfind(['+', '-']) {
        Some(i) if s.len() - i == 3 && s.starts_with(' ') => format!("{}{}0{}", &s[1..i], &s[i..=i], &s[i + 1..]),
        _ => s,
    }
}

#[derive(Default)]
struct FindOptimum {
    /// `OptimizationsDone`.
    done: usize,
    fixed: Vec<bool>,
    min_param: Vec<f64>,
    max_param: Vec<f64>,
    best: Vec<f64>,
    min_cost: f64,
    optimal_finish: bool,
    no_values: i64,
    opt_tol: f64,
    /// SecantSearch's SAVEd state.
    secant_no: usize,
    step: f64,
    maxstep: f64,
    relax: f64,
    x: [f64; 3],
    f: [f64; 2],
}

impl FindOptimum {
    fn secant(&mut self, param: &mut [f64], cost: f64, p: &ParamMap) -> Result<(), SolverError> {
        let free: Vec<usize> = (0..param.len()).filter(|&i| !self.fixed[i]).collect();
        if free.len() != 1 {
            return Err(SolverError::fatal(CALLER, "Secant search implemented only for one parameter"));
        }
        let j = free[0];
        self.secant_no += 1;
        let no = self.secant_no;
        if no == 1 {
            self.maxstep = p.real("max step size").unwrap_or(0.0);
            self.step = p.real("step size").unwrap_or(1.0e-3 * (self.max_param[j] - self.min_param[j]));
            self.relax = p.real("relaxation factor").unwrap_or(1.0);
        }
        let [_, x1, x2] = self.x;
        self.x = [x1, x2, 0.0];
        self.f = [self.f[1], cost];
        let [x0, x1, _] = self.x;
        let [f0, f1] = self.f;
        let x2 = if no <= 2 {
            param[j] + (no - 1) as f64 * self.step
        } else if f1.abs() < self.opt_tol {
            x1
        } else {
            let mut dx = self.relax * f1 * (x1 - x0) / (f1 - f0);
            if dx.abs() > self.maxstep {
                dx = self.maxstep.copysign(dx);
            }
            x1 - dx
        };
        self.x[2] = x2;
        param[j] = x2;
        Ok(())
    }
}

impl Solver for FindOptimum {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        let mut param: Vec<f64> = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "There are no parameters to optimize!"))?.values.to_vec();
        let n = param.len();
        if n == 0 {
            return Err(SolverError::fatal(CALLER, "There are no parameters to optimize!"));
        }
        let mut cost = 0.0;
        let mut finish = false;
        if self.done == 0 {
            ctx.warn(CALLER, "The solver will soon be removed!");
            if ctx.sim.vars.get("found better").is_none() {
                ctx.sim.vars.add(Variable::global("found better", -1.0))?;
            }
            self.optimal_finish = p.logical("optimal finish") == Some(true);
            self.no_values = ctx.simulation_params().integer("timestep intervals").map_or(0, i64::from);
            self.opt_tol = p.real("optimization tolerance").unwrap_or(0.0);
            self.min_param = vec![-f64::MAX; n];
            self.max_param = vec![f64::MAX; n];
            self.best = vec![0.0; n];
            self.fixed = vec![false; n];
            self.min_cost = f64::MAX;
            for i in 0..n {
                let k = i + 1;
                self.fixed[i] = p.logical(&format!("fixed parameter {k}")) == Some(true);
                let init = p.real(&format!("initial parameter {k}"));
                param[i] = init.unwrap_or(0.0);
                if !self.fixed[i] {
                    let minv = p.real(&format!("min parameter {k}"));
                    let maxv = p.real(&format!("max parameter {k}"));
                    if let Some(v) = minv {
                        self.min_param[i] = v;
                    }
                    if let Some(v) = maxv {
                        self.max_param[i] = v;
                    }
                    if init.is_none() {
                        match (minv, maxv) {
                            (Some(a), Some(b)) => param[i] = 0.5 * (a + b),
                            (Some(a), None) => param[i] = a,
                            (None, Some(b)) => param[i] = b,
                            (None, None) => {}
                        }
                    }
                }
            }
            if self.fixed.iter().all(|&f| f) {
                ctx.warn(CALLER, "All parameters are fixed, no optimization!");
                ctx.primary_mut().expect("checked above").values.copy_from_slice(&param);
                return Ok(());
            }
        } else {
            ctx.sim.vars.set_scalar("found better", -1.0);
            let sim = ctx.simulation_params().clone();
            cost = match ctx.creal(&sim, "cost function")? {
                Some(c) => c,
                None => match ctx.creal(&p, "cost function")? {
                    Some(c) => c,
                    None => {
                        let name = p
                            .string("cost function name")
                            .ok_or_else(|| SolverError::fatal(CALLER, "Give Cost Function or its name"))?;
                        sim.real(name).ok_or_else(|| SolverError::fatal(CALLER, "Cost with the given name was not found"))?
                    }
                },
            };
            if let Some(t) = p.real("cost function target") {
                cost -= t;
            }
            if p.logical("cost function absolute") == Some(true) {
                cost = cost.abs();
            } else if p.logical("cost function maximize") == Some(true) {
                cost = -cost;
            }
            if cost < self.min_cost {
                ctx.sim.vars.set_scalar("found better", 1.0);
                self.min_cost = cost;
                self.best.copy_from_slice(&param);
                ctx.info(CALLER, &format!("Found New Minimum Cost:{}", es_e3(self.min_cost)), 4);
            }
            if self.optimal_finish && self.done as i64 == self.no_values - 1 {
                ctx.info(CALLER, "Peforming the last step with the best so far", 4);
                param.copy_from_slice(&self.best);
                finish = true;
            }
        }
        if !finish {
            self.secant(&mut param, cost, &p)?;
        }
        ctx.primary_mut().expect("checked above").values.copy_from_slice(&param);
        self.done += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn secant_steps() {
        let mut o = FindOptimum { fixed: vec![false], min_param: vec![-f64::MAX], max_param: vec![f64::MAX], ..Default::default() };
        let mut p = ParamMap::new();
        let r = |v: f64| elmer_core::Value::Real { rows: 1, cols: 1, data: vec![v] };
        p.insert("Step Size", r(0.1));
        p.insert("Max Step Size", r(1.0));
        p.insert("Relaxation Factor", r(0.5));
        let mut x = [0.0];
        // Start, then one step, then secant on f(x) = x - 0.5.
        o.secant(&mut x, 99.0, &p).unwrap();
        assert_eq!(x, [0.0]);
        o.secant(&mut x, -0.5, &p).unwrap();
        assert_eq!(x, [0.1]);
        o.secant(&mut x, 0.1 - 0.5, &p).unwrap();
        // x1 = 0.1, x0 = 0, f1 = -0.4, f0 = -0.5: dx = 0.5 * -0.4 * 0.1 / 0.1.
        assert!((x[0] - (0.1 + 0.2)).abs() < 1e-12);
        assert!(check(&p).is_err());
        assert_eq!(es_e3(8.316270e-4), "  8.316270E-004");
        assert_eq!(es_e3(-3.042330e-4), " -3.042330E-004");
    }
}
