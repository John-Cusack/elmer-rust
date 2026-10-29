//! Port of `fem/src/modules/PressurePrecond.F90` (Elmer @ 9f6af2f85, the
//! lane fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/PressurePrecond.F90` — two F77-global
//! entries, `pressureprecond_` (a dummy no-op body: one commented-out
//! print) and `pressureprecond_init0_` (twelve `ListAdd*` default
//! initializations for block preconditioning with the ParStokes solver).
//! The emitted dynamic symbols are exactly `pressureprecond_` and
//! `pressureprecond_init0_` (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/pressureprecond_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset. The `ListAdd*` string and
//! real literals live verbatim in the Fortran wrappers (the same
//! compiler evaluates them, in the same statement order as the
//! original), so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirror PressurePrecond.F90 statement for statement):
//! - `pressureprecond_` is a no-op, exactly as in the original (the only
//!   body line is a comment); all four dummies are received and ignored.
//! - `pressureprecond_init0_` reproduces the twelve `ListAdd*` calls in
//!   order: `GetSolverParams()` first (the original ignores its `Solver`
//!   dummy and reads the current-solver context instead), then the five
//!   unconditional adds, then the seven `ListAddNew*` add-if-absent
//!   defaults — so a `.sif` that pre-sets e.g. `Variable` keeps its value
//!   in both arms identically.
//! - No Rust threads: the solver entries run serially under ElmerSolver,
//!   and the port holds no shared Rust state at all (the wrapper-side
//!   `acc_params` pointer is set-before-use on every entry call, exactly
//!   like the original's `SolverParams` local).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ppa_begin();
    fn ppa_init_exec_solver();
    fn ppa_init_variable_output();
    fn ppa_init_bubbles();
    fn ppa_init_skip_change();
    fn ppa_init_back_rotate();
    fn ppa_init_variable();
    fn ppa_init_lin_solver();
    fn ppa_init_lin_method();
    fn ppa_init_lin_maxiter();
    fn ppa_init_lin_precond();
    fn ppa_init_lin_tol();
    fn ppa_init_lin_abort();
}

/// `pressureprecond_init0_` — the `_Init0` default-install entry.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`; every dummy is unused by the body
/// beyond the `GetSolverParams()` context lookup, exactly as in the
/// original.
///
/// # Safety
///
/// The loader passes the Elmer `Model`/`Solver` addresses; `dt` and the
/// transient flag are received and ignored, as in the original.
#[export_name = "pressureprecond_init0_"]
pub unsafe extern "C" fn pressure_precond_init0_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ppa_begin();
        ppa_init_exec_solver();
        ppa_init_variable_output();
        ppa_init_bubbles();
        ppa_init_skip_change();
        ppa_init_back_rotate();
        ppa_init_variable();
        ppa_init_lin_solver();
        ppa_init_lin_method();
        ppa_init_lin_maxiter();
        ppa_init_lin_precond();
        ppa_init_lin_tol();
        ppa_init_lin_abort();
    }
}

/// `pressureprecond_` — the dummy solver entry. The original body is a
/// single commented-out print, so the port is a no-op with the same
/// four ignored dummies.
///
/// # Safety
///
/// See [`pressure_precond_init0_entry`]; all arguments are ignored here.
#[export_name = "pressureprecond_"]
pub unsafe extern "C" fn pressure_precond_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}
