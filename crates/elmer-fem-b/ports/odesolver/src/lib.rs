//! Port of `fem/src/modules/OdeSolver.F90` (the campaign fork) for lane
//! T021-B, at the solver-plugin boundary (whole `.so` -> Rust cdylib
//! with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/OdeSolver.F90` — two F77-global entries
//! `odesolver_` and `odesolver_init_` (proven by `nm -D` on the reference
//! `.so`); the `OdeMatrix` `CONTAINS` internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/odesolver_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset.
//! Order-sensitive fills (`CRS_SetMatrixElement` sequences, whole-array
//! `= 0` fills) and all formatted I/O (`Info`/`Fatal`, `I2S` string
//! building) stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors OdeSolver.F90 statement for statement):
//! - The original has no `SAVE` state (all locals are per-call), so the
//!   port keeps no `static` state at all; mass/damping allocation
//!   lifetimes live in the Elmer matrix, exactly as in the original.
//! - `Model`, `dt`, `TransientSimulation` are unused by both bodies,
//!   exactly as in the original (received and ignored).
//! - `MaxIter` defaults to 1 when `Nonlinear System Max Iterations` is
//!   absent; `TimeOrder` reads `Time Derivative Order` with no default,
//!   exactly as in the original.
//! - The `IF(.FALSE.)` `PRINT` block is dead code in the original and is
//!   not transcribed (it never executes).
//! - No Rust threads: the solver calls entries serially, and with no
//!   shared state there is nothing to guard.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn oda_begin(h_solver: *const c_void);
    fn oda_init_params();
    fn oda_info_start();
    fn oda_params_begin();
    fn oda_maxiter() -> c_int;
    fn oda_timeorder() -> c_int;
    fn oda_info_iter(iter: c_int);
    fn oda_default_initialize();
    fn oda_matrix(timeorder: c_int);
    fn oda_finish_bulk();
    fn oda_finish_bnd();
    fn oda_finish_assembly();
    fn oda_solve() -> c_double;
    fn oda_converged() -> c_int;
    fn oda_info_done();
}

/// `odesolver_init_` — the initializer entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; only the solver context
/// matters (via `GetSolverParams`), the rest are unused, exactly as in
/// the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "odesolver_init_"]
pub unsafe extern "C" fn odesolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { oda_begin(solver) };
    unsafe { oda_init_params() };
}

/// `odesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "odesolver_"]
pub unsafe extern "C" fn odesolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { oda_begin(solver) };
    unsafe { oda_info_start() };

    unsafe { oda_params_begin() };
    let max_iter = unsafe { oda_maxiter() };
    let time_order = unsafe { oda_timeorder() };

    for iter in 1..=max_iter {
        if max_iter > 1 {
            unsafe { oda_info_iter(iter) };
        }

        unsafe { oda_default_initialize() };

        unsafe { oda_matrix(time_order) };

        unsafe {
            oda_finish_bulk();
            oda_finish_bnd();
        }

        // This sets the time-integration and is therefore imperative.
        unsafe { oda_finish_assembly() };

        // Norm is local and unused, as in the original.
        let _norm: c_double = unsafe { oda_solve() };

        if unsafe { oda_converged() } == 1 {
            break;
        }
    }

    unsafe { oda_info_done() };
}
