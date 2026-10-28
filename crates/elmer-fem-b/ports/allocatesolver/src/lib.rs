//! Port of `fem/src/modules/AllocateSolver.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/AllocateSolver.F90` — two F77-global
//! entries (`allocatesolver_`, `allocatesolver_init_`), no `CONTAINS`
//! internals. The emitted dynamic symbols are exactly those two (proven by
//! `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/allocate_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset.
//!
//! Faithfulness notes (mirrors AllocateSolver.F90 statement for statement):
//! - `AllocateSolver_init` calls `ListAddNewString(Solver % Values,
//!   'Exec Solver', 'never')` and nothing else; the wrapper issues the
//!   identical call against the same list object, so first-match-wins
//!   lookup order against `.sif`-supplied `Exec Solver` values is
//!   preserved by construction (the append lands at the list tail in both
//!   arms).
//! - `AllocateSolver` calls `Info('AllocateSolver', ...)` with the exact
//!   message text; the `.sif` `Exec Solver` value decides whether the
//!   entry runs at all (Elmer-side dispatch, identical in both arms).
//! - `Model`, `dt`, `TransientSimulation` are unused by both bodies,
//!   exactly as in the original (received and ignored).
//! - No `SAVE` state anywhere in the original; none here either.
//! - No Rust threads: the solver entries run serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ala_begin(h_solver: *mut c_void);
    fn ala_init_exec_never();
    fn ala_info();
}

/// `allocatesolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Argument order and by-reference travel match the Fortran dummy
/// list `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "allocatesolver_init_"]
pub unsafe extern "C" fn allocatesolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ala_begin(solver);
        ala_init_exec_never();
    }
}

/// `allocatesolver_` — the solver entry. Runs only when the Elmer-side
/// `Exec Solver` dispatch selects it (the `_init` above defaults it to
/// `never`), exactly as in the original.
///
/// # Safety
///
/// Same contract as [`allocatesolver_init_entry`].
#[export_name = "allocatesolver_"]
pub unsafe extern "C" fn allocatesolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ala_begin(solver);
        ala_info();
    }
}
