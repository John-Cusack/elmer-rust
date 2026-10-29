//! Port of `fem/src/modules/VelocityPrecond.F90` (Elmer fork) for lane
//! T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib
//! with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/VelocityPrecond.F90` — two F77-global
//! entries, `velocityprecond_init0_` (twelve unconditional `ListAdd*`
//! solver-param defaults) and `velocityprecond_` (body is a single
//! comment: a no-op). The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/velocityprecond_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset. The `CHARACTER` literals
//! and the `OPTIONAL` dummy of `GetSolverParams` stay inside Fortran, so
//! the hidden-length contract is never crossed.
//!
//! Faithfulness notes (mirrors VelocityPrecond.F90 statement for
//! statement):
//! - `VelocityPrecond_Init0`: `SolverParams => GetSolverParams()` then
//!   the twelve `ListAdd*` calls in source order, with identical
//!   `ListAdd*` vs `ListAddNew*` choice per statement (the `New` variants
//!   keep a user-preset value, which the config-B differential observes).
//! - `VelocityPrecond`: the body is one comment line, so the entry is a
//!   no-op that receives and ignores all four dummies, exactly as the
//!   original.
//! - No `SAVE` state anywhere in the original; the accessor's pointer is
//!   re-acquired on every `Init0` call like the original's local. No
//!   Rust threads, no Rust statics.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn vpa_begin();
    fn vpa_add_exec_solver();
    fn vpa_add_variable_output();
    fn vpa_add_bubbles();
    fn vpa_add_skip_change();
    fn vpa_add_back_rotate();
    fn vpa_add_variable();
    fn vpa_add_linsys_solver();
    fn vpa_add_linsys_method();
    fn vpa_add_linsys_maxiter();
    fn vpa_add_linsys_precond();
    fn vpa_add_linsys_tol();
    fn vpa_add_linsys_abort();
}

/// `velocityprecond_init0_` — the solver init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`; all four dummies are unused by the
/// body except through the ambient `GetSolverParams()` context, exactly
/// as in the original.
///
/// # Safety
///
/// Called by the Elmer loader with its own context current; takes no
/// readable arguments (all dummies ignored).
#[export_name = "velocityprecond_init0_"]
pub unsafe extern "C" fn velocity_precond_init0_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        vpa_begin();
        vpa_add_exec_solver();
        vpa_add_variable_output();
        vpa_add_bubbles();
        vpa_add_skip_change();
        vpa_add_back_rotate();
        vpa_add_variable();
        vpa_add_linsys_solver();
        vpa_add_linsys_method();
        vpa_add_linsys_maxiter();
        vpa_add_linsys_precond();
        vpa_add_linsys_tol();
        vpa_add_linsys_abort();
    }
}

/// `velocityprecond_` — the solver body. The original body is a single
/// comment line, so this is a no-op receiving and ignoring all four
/// dummies.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "velocityprecond_"]
pub unsafe extern "C" fn velocity_precond_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}
