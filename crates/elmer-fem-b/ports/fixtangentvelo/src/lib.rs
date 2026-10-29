//! Port of `fem/src/modules/FixTangentVelo.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/FixTangentVelo.F90` — two F77-global
//! entries (`fixtangentvelo_`, `fixtangentvelo_init_`) with one
//! `CONTAINS` internal. The emitted dynamic symbols are exactly those
//! two (proven by `nm -D` on the reference `.so`); the internal emits
//! only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/fixtangentvelo_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (edge vectors, the 3x3 solve, the `MAX`/`MIN`
//! accumulation, `ComputeNorm`) and the order-sensitive reductions
//! (`SUM`) stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FixTangentVelo.F90 statement for
//! statement):
//! - The `_init_` entry maps to the same `ListAddNewLogical` /
//!   `ListAddNewString` calls (`No Matrix`, `-nooutput -global
//!   FixTangent_var`).
//! - `SAVE` semantics replicated: first-visit sizing of `FixPerm`,
//!   `VariableAddVector('FixVelo')`, and the `Visited`/`Nodes`/
//!   `FixVar`/`VeloVar`/`FixVals`/`ActiveCoord`/`DownstreamOnly` state
//!   live as wrapper-side SAVE state with identical lifetimes; the
//!   `Params`/`Mesh` entry-scope pointers are re-fetched per call and
//!   resolve to the same targets.
//! - The entry prints the same `Info` line on every invocation, zeroes
//!   `FixVals` per invocation, sweeps the same boundary-element range,
//!   guards the `VectorValuesRange` report with the same
//!   `InfoActive(20)` check, applies the fix only when the `Fix Tangent
//!   Velocity` keyword is present, and assigns the `ComputeNorm` result
//!   to `Solver % Variable` exactly as the original.
//! - The `ApplyTangentFix` internal procedure is kept verbatim as a
//!   module procedure and invoked with `.FALSE.`, exactly as the
//!   original calls it.
//! - `Model`, `dt`, and the transient flag are unused by the entries,
//!   exactly as in the original (received and ignored).
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ftv_begin(h_solver: *mut c_void);
    fn ftv_need_first() -> c_int;
    fn ftv_first();
    fn ftv_zero();
    fn ftv_num_bulk() -> c_int;
    fn ftv_num_bnd() -> c_int;
    fn ftv_elem(t: c_int);
    fn ftv_range();
    fn ftv_fix_flag() -> c_int;
    fn ftv_apply();
    fn ftv_tail() -> c_double;
    fn ftv_init_begin(h_solver: *mut c_void);
    fn ftv_init_body();
}

/// `fixtangentvelo_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `FixTangentVelo_Init` (`No Matrix` flag, default
/// global `FixTangent_var` variable).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "fixtangentvelo_init_"]
pub unsafe extern "C" fn fixtangentvelo_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ftv_init_begin(solver);
        ftv_init_body();
    }
}

/// `fixtangentvelo_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`fixtangentvelo_init_entry`].
#[export_name = "fixtangentvelo_"]
pub unsafe extern "C" fn fixtangentvelo_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ftv_begin(solver);
        if ftv_need_first() != 0 {
            ftv_first();
        }
        ftv_zero();
        let nb = ftv_num_bulk();
        let nn = ftv_num_bnd();
        for t in nb + 1..=nb + nn {
            ftv_elem(t);
        }
        ftv_range();
        if ftv_fix_flag() != 0 {
            ftv_apply();
        }
        let _nrm: c_double = ftv_tail();
    }
}
