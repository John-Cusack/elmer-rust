//! Property evaluation: constants, tensors, `Variable` tables (linear and
//! cubic), MATC expressions and `Procedure` user functions, evaluated at mesh
//! nodes or integration points.
//!
//! Reference: `ListGetConstReal` / `ListGetReal` / `ListGetRealArray` /
//! `ListGetElementReal` in `fem/src/Lists.F90`, the curve routines in
//! `fem/src/GeneralUtils.F90` (see [`curve`]) and symbol lookup in
//! `fem/src/Load.c` (see [`udf`]).
//!
//! A `elmer_core::Value` is turned into a [`Compiled`] property once; the
//! per-node loop then does no string work beyond one batched field fetch per
//! dependency variable ([`FieldLookup::values`]).
#![forbid(unsafe_code)]

mod compiled;
pub mod curve;
pub mod udf;

pub use compiled::{compile, compile_with, CompileOptions, Compiled, ExprError};
pub use udf::{ArrayFn, ConstFn, NodeFn, Udf, UdfRegistry};

/// Marker for "no value at this node", the `HUGE(1.0_dp)` sentinel of
/// `VarsToValuesOnNodes`. A property whose dependency is missing at a node
/// evaluates to zero there, as in `ListGetReal`.
pub const MISSING: f64 = f64::MAX;

/// Read access to the fields a property may depend on.
pub trait FieldLookup {
    /// Value of the named variable at a node (`VariableGet` + `Perm`), or `None`
    /// where the variable has no value there. Global variables ignore `node`.
    /// Names reach the lookup lower-cased and trimmed; `time` and
    /// `coordinate [1-3]` are answered by [`time`](Self::time) and
    /// [`coord`](Self::coord) and never get here. `prev <name>` is passed
    /// through verbatim for lookups that keep previous time levels.
    fn value(&self, var: &str, node: u32) -> Option<f64>;

    fn coord(&self, node: u32) -> [f64; 3];

    fn time(&self) -> f64;

    /// Batched [`value`](Self::value): `out[i]` for `nodes[i]`, [`MISSING`]
    /// where absent. Evaluation calls this once per variable per batch, so an
    /// implementation can resolve the name once. The default loops over
    /// `value`.
    fn values(&self, var: &str, nodes: &[u32], out: &mut [f64]) {
        for (o, &n) in out.iter_mut().zip(nodes) {
            *o = self.value(var, n).unwrap_or(MISSING);
        }
    }
}

/// MATC evaluation, implemented by the MATC binding (S1). Results are
/// flattened row-major.
pub trait MatcEngine {
    /// Evaluates `expr` with `args` bound to the MATC array `tx` (`tx(0)`,
    /// `tx(1)`, ...), as `SetGetMatcParams` does, and returns the result
    /// values.
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String>;

    /// As [`eval`](Self::eval) but names the array the arguments are bound
    /// to: `tx` for variable-dependent properties, `st` for the time-only form
    /// of `ListGetReal` (`Real MATC` with no `Variable`). Engines that bind a
    /// fixed name can ignore `var`.
    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        let _ = var;
        self.eval(expr, args)
    }
}
