//! Keyword values that are not plain constants (tables, MATC, procedures),
//! evaluated with `elmer-expr` against the simulation's fields.
//!
//! `Procedure` values resolve against the simulation's own
//! [`UdfRegistry`] ([`crate::Backends::udfs`]) first, then the process-wide
//! one. MATC values go to the installed [`MatcEngine`]
//! ([`crate::Backends::matc`]); without one they fail.

use std::rc::Rc;

use elmer_core::{MeshView, ParamMap, Params, Value};
use elmer_expr::{compile, compile_with, CompileOptions, Compiled, FieldLookup, MatcEngine, UdfRegistry};

use crate::error::SolverError;
use crate::variable::VariableStore;

/// The fields of a simulation as `elmer-expr` reads them: a variable or
/// component alias at a node through its perm, a global's first value,
/// mesh coordinates and the `time` global.
pub struct SimFields<'a> {
    pub vars: &'a VariableStore,
    pub mesh: &'a dyn MeshView,
}

/// `VariableGet` + `Perm` at one node: `None` where the variable has no dof.
pub(crate) fn lookup(vars: &VariableStore, var: &str, node: u32) -> Option<f64> {
    let v = vars.strided(var)?;
    match v.perm {
        Some(perm) => {
            let s = *perm.get(node as usize)?;
            if s == elmer_core::NO_DOF {
                None
            } else {
                v.values.get(s as usize * v.stride + v.offset).copied()
            }
        }
        // Globals and perm-less fields: one value (globals) or node-indexed.
        None if v.values.len() == 1 => Some(v.values[0]),
        None => v.values.get(node as usize * v.stride + v.offset).copied(),
    }
}

impl FieldLookup for SimFields<'_> {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        lookup(self.vars, var, node)
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }

    fn time(&self) -> f64 {
        self.vars.scalar("time").unwrap_or(0.0)
    }
}

/// Stand-in until a MATC binding is installed in [`crate::Backends::matc`]:
/// every MATC value fails.
pub struct NoMatc;

impl MatcEngine for NoMatc {
    fn eval(&self, expr: &str, _args: &[f64]) -> Result<Vec<f64>, String> {
        Err(format!("no MATC engine installed for: {expr}"))
    }
}

/// A MATC engine shared with its other users. MATC keeps one evaluator per
/// process (variables defined while the `.sif` is read stay visible to the
/// keyword values evaluated later), so the front end and the solver hold the
/// same non-`Send` handle; this forwards to it.
#[derive(Clone)]
pub struct SharedMatc(pub Rc<dyn MatcEngine>);

impl MatcEngine for SharedMatc {
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.0.eval(expr, args)
    }

    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.0.eval_with(var, expr, args)
    }
}

/// A [`MatcEngine`] over a text MATC evaluator such as
/// `elmer_sif::matc::matc_eval` (Elmer's `matc_c`, Load.c). Arguments are
/// bound the way `SetGetMatcParams` (Lists.F90:128) does through
/// `MatcSetParams`: the array `tx` (or `st`) becomes the 1 x n row of
/// argument values, here by an assignment whose numbers round-trip
/// exactly; no arguments bind nothing (`GetMatcReal`). The printed result
/// is read back as `READ(res,*)` does, so values carry MATC's 12 printed
/// digits, as in Elmer.
///
/// Binding and evaluating are two calls into one process-wide evaluator:
/// install it in one simulation at a time ([`crate::Simulation::set_matc`]
/// takes a non-`Send` handle for that reason).
pub struct TextMatc<F>(pub F);

impl<F: Fn(&str) -> Result<String, String>> MatcEngine for TextMatc<F> {
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.eval_with("tx", expr, args)
    }

    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        if !args.is_empty() {
            // `{:e}` is the shortest exact form and never has a space after a
            // sign, so MATC reads one row element per number.
            let row: Vec<String> = args.iter().map(|a| format!("{a:e}")).collect();
            (self.0)(&format!("{var} = {}", row.join(" ")))?;
        }
        let out = (self.0)(expr)?;
        out.split_whitespace()
            .map(|t| t.parse::<f64>().map_err(|_| format!("MATC result {out:?} of {expr:?} is not numeric")))
            .collect()
    }
}

/// [`compile_with`] against `udfs`, falling back to the process-wide
/// registry for procedures `udfs` does not know. Any other compile error
/// comes back from the property's evaluation, as with [`compile`].
pub fn compile_value(v: &Value, udfs: &UdfRegistry) -> Compiled {
    compile_with(v, udfs, &CompileOptions::default()).unwrap_or_else(|_| compile(v))
}

/// `ListGetReal` of a scalar keyword at `nodes`, procedures resolved in the
/// process-wide registry only. See [`eval_at_nodes_with`].
pub fn eval_at_nodes(
    v: &Value,
    nodes: &[u32],
    fields: &SimFields<'_>,
    matc: Option<&dyn MatcEngine>,
) -> Result<Vec<f64>, SolverError> {
    eval_at_nodes_with(v, nodes, fields, matc, UdfRegistry::global())
}

/// `ListGetReal` of a scalar keyword at `nodes`.
pub fn eval_at_nodes_with(
    v: &Value,
    nodes: &[u32],
    fields: &dyn FieldLookup,
    matc: Option<&dyn MatcEngine>,
    udfs: &UdfRegistry,
) -> Result<Vec<f64>, SolverError> {
    if let Value::Real { data, .. } = v {
        if data.len() == 1 {
            return Ok(vec![data[0]; nodes.len()]);
        }
    }
    let c = compile_value(v, udfs);
    if c.ncomp() != 1 {
        return Err(SolverError::Unsupported(format!("{}x{} value where a scalar is expected", c.rows(), c.cols())));
    }
    let mut out = vec![0.0; nodes.len()];
    c.eval_nodes(nodes, fields, matc.unwrap_or(&NoMatc), &mut out)?;
    Ok(out)
}

/// `ListGetConstReal` of a scalar keyword, procedures resolved in the
/// process-wide registry only.
pub fn eval_const(v: &Value, matc: Option<&dyn MatcEngine>) -> Result<f64, SolverError> {
    eval_const_with(v, matc, UdfRegistry::global())
}

/// `ListGetConstReal` of a scalar keyword: a constant, a MATC expression
/// evaluated once with no arguments, or a constant procedure at the origin.
pub fn eval_const_with(v: &Value, matc: Option<&dyn MatcEngine>, udfs: &UdfRegistry) -> Result<f64, SolverError> {
    if let Value::Real { data, .. } = v {
        if data.len() == 1 {
            return Ok(data[0]);
        }
    }
    Ok(compile_value(v, udfs).eval_const(matc.unwrap_or(&NoMatc), [0.0; 3])?)
}

/// `ListGetCReal` (Lists.F90:4591): `ListGetReal` at node 1, so a value may
/// depend on `time` (a MATC expression with no `Variable` binds it as `st`)
/// or on a field's first node. `None` when the section lacks the keyword.
pub fn eval_creal(
    section: &ParamMap,
    key: &str,
    fields: &dyn FieldLookup,
    matc: Option<&dyn MatcEngine>,
    udfs: &UdfRegistry,
) -> Result<Option<f64>, SolverError> {
    let Some(v) = section.get(key) else { return Ok(None) };
    if let Value::Real { data, .. } = v {
        if let [x] = data.as_slice() {
            return Ok(Some(*x));
        }
    }
    Ok(Some(eval_at_nodes_with(v, &[0], fields, matc, udfs)?[0]))
}

/// `section` with each of `keys` that holds a non-constant value replaced
/// by its [`eval_creal`] value, for code that reads those keywords with
/// [`Params::real`]. Borrows `section` when nothing needs evaluating.
pub fn resolve_creals<'p>(
    section: &'p ParamMap,
    keys: &[&str],
    fields: &dyn FieldLookup,
    matc: Option<&dyn MatcEngine>,
    udfs: &UdfRegistry,
) -> Result<std::borrow::Cow<'p, ParamMap>, SolverError> {
    let dependent = |k: &str| matches!(section.get(k), Some(Value::Dependent(_)));
    if !keys.iter().any(|k| dependent(k)) {
        return Ok(std::borrow::Cow::Borrowed(section));
    }
    let mut p = section.clone();
    for k in keys.iter().filter(|k| dependent(k)) {
        let x = eval_creal(section, k, fields, matc, udfs)?.expect("present: dependent");
        p.insert(k, Value::real(x));
    }
    Ok(std::borrow::Cow::Owned(p))
}
