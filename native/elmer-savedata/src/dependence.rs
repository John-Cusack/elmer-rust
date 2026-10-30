//! `SaveDependence` (fem/src/modules/SaveData/SaveDependence.F90): the
//! functions `Expression 1..` of the solver section tabulated on
//! `Number of Points` equally spaced arguments between `Lower Limit` and
//! `Upper Limit`, written as `I6,ES15.6` then one `ES15.6` per expression
//! (two with `Take Derivative`). With `Show Norm Index = j` the solver's
//! norm is the RMS of column `j`.

use elmer_core::{normalize_key, DependentKind, ParamMap, Params, Value};
use elmer_expr::curve::{cubic_spline, derivate_curve, interpolate_curve};
use elmer_saveio::fortran::{fortran_es, fortran_i};
use elmer_solver::{Solver, SolverCtx, SolverError, Variable};

use crate::scalars::split_output_path;

const CALLER: &str = "SaveDependence";

pub struct SaveDependence {
    section: ParamMap,
}

impl SaveDependence {
    pub fn from_params(params: &ParamMap) -> Result<Self, SolverError> {
        if params.logical("filename numbering") == Some(true) {
            return Err(SolverError::Unsupported("SaveDependence: Filename Numbering".into()));
        }
        Ok(SaveDependence { section: params.clone() })
    }

    fn norm_index(&self) -> usize {
        self.section.integer("show norm index").unwrap_or(0).max(0) as usize
    }
}

/// `ListGetFun(Params, name, x [, dfdx])`.
fn list_get_fun(ctx: &SolverCtx<'_>, v: &Value, x: f64, der: bool) -> Result<(f64, f64), SolverError> {
    match v {
        Value::Real { data, .. } if data.len() == 1 => Ok((data[0], 0.0)),
        Value::Dependent(d) => match &d.kind {
            DependentKind::Table { x: tv, y: fv, cubic } => {
                let spline = cubic.then(|| cubic_spline(tv, fv, false));
                Ok((interpolate_curve(tv, fv, x, spline.as_deref()), derivate_curve(tv, fv, x, spline.as_deref())))
            }
            DependentKind::Matc(expr) if !der => {
                let matc = ctx.sim.backends.matc.as_deref().ok_or_else(|| SolverError::NotWired("MATC engine".into()))?;
                let f = matc.eval(expr, &[x]).map_err(|e| SolverError::fatal(CALLER, e))?;
                Ok((f.first().copied().unwrap_or(0.0), 0.0))
            }
            _ => Err(SolverError::Unsupported(format!("{CALLER}: this kind of Expression{}", if der { " with Take Derivative" } else { "" }))),
        },
        _ => Err(SolverError::Unsupported(format!("{CALLER}: Expression that is not a function"))),
    }
}

impl Solver for SaveDependence {
    /// `SaveDependence_init`: with a norm index, the global variable
    /// `<Equation>_var` unless the section names a `Variable`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if self.norm_index() == 0 || self.section.get("variable").is_some() {
            return Ok(());
        }
        let eq = self.section.string("equation").unwrap_or(CALLER).to_string();
        let name = normalize_key(&format!("{eq}_var"));
        let id = match ctx.sim.vars.id(&name) {
            Some(id) => id,
            None => {
                let mut v = Variable::global(&name, 0.0);
                v.output = false;
                ctx.sim.vars.add(v)?
            }
        };
        let sid = ctx.solver_id();
        let ix = ctx.sim.slots.iter().position(|s| s.id == sid).expect("running solver has a slot");
        ctx.sim.slots[ix].variable = Some(id);
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 && ctx.comm.rank() > 0 {
            return Ok(());
        }
        ctx.info(CALLER, "--------------------------------", 4);
        ctx.info(CALLER, "Saving dependencies in a table", 4);
        let p = &self.section;
        let filename = p.string("filename").unwrap_or("dep.dat").to_string();
        let norm_ind = self.norm_index();
        let der = p.logical("take derivative") == Some(true);
        let n = p.integer("number of points").ok_or_else(|| SolverError::fatal(CALLER, "Number of points"))?;
        if n < 2 {
            return Err(SolverError::fatal(CALLER, "Number of points must be at least 2"));
        }
        let x0 = p.real("lower limit").ok_or_else(|| SolverError::fatal(CALLER, "Lower Limit"))?;
        let x1 = p.real("upper limit").ok_or_else(|| SolverError::fatal(CALLER, "Upper Limit"))?;
        let exprs: Vec<&Value> = (1..=100).map_while(|j| p.get(&format!("expression {j}"))).collect();
        if exprs.is_empty() {
            ctx.warn(CALLER, "No parameter given!");
            return Ok(());
        }
        let (_, path) = split_output_path(ctx, p.string("output directory"), &filename, true)?;
        let mut text = String::new();
        let mut norm = 0.0;
        for i in 1..=n {
            let w = (i - 1) as f64 / (n - 1) as f64;
            let x = x0 + w * (x1 - x0);
            text.push_str(&fortran_i(i as i64, 6));
            text.push_str(&fortran_es(x, 15, 6));
            for (j, e) in exprs.iter().enumerate() {
                let (f, dfdx) = list_get_fun(ctx, e, x, der)?;
                text.push_str(&fortran_es(f, 15, 6));
                if der {
                    text.push_str(&fortran_es(dfdx, 15, 6));
                }
                if norm_ind == j + 1 {
                    norm += f * f;
                }
            }
            text.push_str(" \n");
        }
        std::fs::write(&path, text).map_err(|e| SolverError::fatal(CALLER, format!("{}: {e}", path.display())))?;
        if norm_ind > 0 {
            let norm = (norm / n as f64).sqrt();
            let sid = ctx.solver_id();
            let ix = ctx.sim.slots.iter().position(|s| s.id == sid).expect("running solver has a slot");
            if let Some(id) = ctx.sim.slots[ix].variable {
                let v = ctx.sim.vars.by_id_mut(id);
                v.values.iter_mut().for_each(|x| *x = norm);
                v.conv.norm = norm;
            }
        }
        Ok(())
    }
}
