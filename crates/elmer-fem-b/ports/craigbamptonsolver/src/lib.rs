//! Port of `fem/src/modules/CraigBamptonSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/CraigBamptonSolver.F90` — two F77-global
//! entries (`craigbamptonsolver_`, `craigbamptonsolver_init_`) with four
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/craigbamptonsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM(x*Ax)`, the Frobenius norm) and every formatted/file WRITE stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors CraigBamptonSolver.F90 statement for
//! statement):
//! - The `_Init` entry's conditional `Variable` default
//!   (`-nooutput -global <Equation>_var` when `Calculate Matrix Norm` is
//!   set and no `Variable` is present yet) lives in the wrapper, keeping
//!   the string handling in Fortran; the entry only binds the handle.
//! - `Model`, `dt`, and the transient flag are unused by the main body
//!   itself (received and ignored), except that `Model` provides the
//!   solver table the elasticity-solver lookup scans — exactly as in the
//!   original.
//! - The dead `IF (SolverId == 0)` re-check before the constraint-mode
//!   report is not replicated: the lookup above already Fatals on zero,
//!   so it is unreachable (no observable state either way).
//! - The reduction work arrays are allocated once per invocation at the
//!   original extents (`Ahat(NoComponentModes, NoComponentModes)` zeroed,
//!   `Ax(A % NumberOfRows)`) and dropped at return, matching the
//!   original's per-call allocate/deallocate pair. Allocation failure
//!   calls the same `Fatal('CraigBamptonSolver', ...)` via `try_reserve`.
//!   The original only deallocates `Ahat` (its `Ax` leaks); freeing both
//!   here is unobservable — and the `Save Reduction Matrix = False`
//!   input, where the reference would hit `DEALLOCATE` of a
//!   never-allocated `Ahat`, is not exercised by the proof.
//! - `A % Values` is swapped to the bulk (stiffness) / bulk-mass values
//!   per reduction pass and restored to the saved pointer on every path,
//!   exactly as in the original; the bulk-mass-missing `Warn` + `CYCLE`
//!   path is kept.
//! - No Rust threads: the solver entries run serially; all entry state
//!   lives wrapper-side with identical lifetimes, so no `Mutex` guard is
//!   needed.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn cba_init_begin(h_solver: *mut c_void);
    fn cba_init_body();
    fn cba_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn cba_info_start();
    fn cba_resolve() -> c_int;
    fn cba_fatal_noname();
    fn cba_fatal_nosolver();
    fn cba_fatal_nomatrix();
    fn cba_fatal_nobulk();
    fn cba_fatal_noeigen();
    fn cba_fatal_alloc();
    fn cba_save_dofs();
    fn cba_want_basis() -> c_int;
    fn cba_save_basis();
    fn cba_want_matrix() -> c_int;
    fn cba_ncomp() -> c_int;
    fn cba_nrows() -> c_int;
    fn cba_matrix_begin(mno: c_int) -> c_int;
    fn cba_mode_pair(
        k: c_int,
        l: c_int,
        ncomp: c_int,
        n: c_int,
        ahat: *mut c_double,
        ax: *mut c_double,
    );
    fn cba_save_matrix(ncomp: c_int, ahat: *mut c_double);
    fn cba_want_norm() -> c_int;
    fn cba_norm(mno: c_int, ncomp: c_int, ahat: *mut c_double);
    fn cba_restore();
    fn cba_info_done();
}

/// `craigbamptonsolver_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `CraigBamptonSolver_Init` (conditional default
/// `Variable` for the matrix-norm global).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "craigbamptonsolver_init_"]
pub unsafe extern "C" fn craigbamptonsolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        cba_init_begin(solver);
        cba_init_body();
    }
}

/// `craigbamptonsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`craigbamptonsolver_init_entry`]; additionally
/// `model` must be the Elmer `Model_t` address the loader passes.
#[export_name = "craigbamptonsolver_"]
pub unsafe extern "C" fn craigbamptonsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        cba_begin(model, solver);
        cba_info_start();
        match cba_resolve() {
            0 => {}
            1 => cba_fatal_noname(),
            2 => cba_fatal_nosolver(),
            3 => cba_fatal_nomatrix(),
            4 => cba_fatal_nobulk(),
            _ => cba_fatal_noeigen(),
        }
        cba_save_dofs();
        if cba_want_basis() != 0 {
            cba_save_basis();
        }
        if cba_want_matrix() != 0 {
            let ncomp = cba_ncomp() as usize;
            let n = cba_nrows() as usize;
            let mut ahat = Vec::<c_double>::new();
            let mut ax = Vec::<c_double>::new();
            if ahat.try_reserve(ncomp * ncomp).is_err() || ax.try_reserve(n).is_err() {
                cba_fatal_alloc();
            }
            ahat.resize(ncomp * ncomp, 0.0);
            ax.resize(n, 0.0);
            for mno in 1..=2 {
                if cba_matrix_begin(mno) != 0 {
                    continue;
                }
                for k in 1..=ncomp {
                    for l in 1..=ncomp {
                        cba_mode_pair(
                            k as c_int,
                            l as c_int,
                            ncomp as c_int,
                            n as c_int,
                            ahat.as_mut_ptr(),
                            ax.as_mut_ptr(),
                        );
                    }
                }
                cba_save_matrix(ncomp as c_int, ahat.as_mut_ptr());
                if cba_want_norm() != 0 {
                    cba_norm(mno, ncomp as c_int, ahat.as_mut_ptr());
                }
            }
        }
        cba_restore();
        cba_info_done();
    }
}
