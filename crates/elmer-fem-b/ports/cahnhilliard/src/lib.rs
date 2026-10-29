//! Port of `fem/src/modules/CahnHilliard.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/CahnHilliard.F90` — one F77-global entry
//! `cahnhilliardsolver_` with three `CONTAINS` internals (`BulkAssembly`,
//! `LocalMatrix`, `BoundaryAssembly`, `BoundaryLocalMatrix`; the last is
//! an empty stub in the original). The emitted dynamic symbol is exactly
//! `cahnhilliardsolver_` (proven by `nm -D` on the reference `.so`); the
//! internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/cahnhilliard_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the 2x2-block Galerkin assembly, `SUM`
//! reductions, `SQRT`, the list-directed `WRITE`+`Info` reports) stays in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors CahnHilliard.F90 statement for statement):
//! - The nonlinear iteration loop (`DO iter=1,NonlinIter` with the `EXIT`
//!   on `NonlinTol`) is driven from Rust with identical bounds and the
//!   same convergence queries; the `!$omp parallel do` drivers serialize
//!   to index order `1..=n`, which matches the serial reference build
//!   (the reference `.so` carries no OpenMP symbols, so its `!$omp`
//!   directives are comments and its assembly order is already serial).
//! - `BoundaryLocalMatrix` is an empty stub in the original, so the
//!   boundary sweep replicates only its driver calls
//!   (`GetBoundaryElement` + `ActiveBoundaryElement` + the `n`/`nd`
//!   getters) with no further action, exactly as the original.
//! - `NewtonLinearization` is set-but-never-read in the original; the
//!   wrapper replicates the flag and its update (`ch_newton_update`)
//!   with identical inputs (no observable effect, kept for fidelity).
//! - `dt` is received and ignored: the original's `Default1stOrderTime`
//!   call takes no timestep argument, and nothing else reads `dt`.
//!   `Model` is likewise unused by the body, exactly as in the original.
//! - `Norm = DefaultSolve()` is kept and reported via the same two
//!   `WRITE(Message,*)` + `Info(..., Level=4)` lines, inside Fortran so
//!   the list-directed formatting is byte-identical.
//! - No Rust threads: the solver entry runs serially, so no `static`
//!   sharing guard is needed — the port holds no Rust-side state at all
//!   (all SAVE-equivalent state lives in the wrappers with the original's
//!   lifetimes); the entry body performs no allocation, hence no
//!   `Fatal`-on-alloc path exists to replicate.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ch_begin(h_solver: *mut c_void, transient: c_int);
    fn ch_read_params();
    fn ch_nonlin_iter() -> c_int;
    fn ch_default_initialize();
    fn ch_nof_active() -> c_int;
    fn ch_local_elem(t: c_int);
    fn ch_default_finish_bulk();
    fn ch_nof_boundary() -> c_int;
    fn ch_boundary_elem(t: c_int);
    fn ch_default_finish();
    fn ch_default_dirichlet();
    fn ch_solve() -> c_double;
    fn ch_newton_update(iter: c_int);
    fn ch_converged() -> c_int;
}

/// `cahnhilliardsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are
/// received and ignored, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read (the Elmer `LOGICAL` flag).
#[export_name = "cahnhilliardsolver_"]
pub unsafe extern "C" fn cahnhilliard_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        ch_begin(solver, *transient);
        ch_read_params();
        let nonlin = ch_nonlin_iter();
        // Mirrors `DO iter=1,NonlinIter`: `1..=nonlin` is empty exactly
        // when the Fortran loop trips zero times (NonlinIter < 1).
        for iter in 1..=nonlin {
            // BulkAssembly: serial index order matches the serial
            // reference build (no OpenMP in the reference `.so`).
            ch_default_initialize();
            let nact = ch_nof_active();
            for t in 1..=nact {
                ch_local_elem(t);
            }
            ch_default_finish_bulk();
            // BoundaryAssembly over the empty BoundaryLocalMatrix stub.
            let nbnd = ch_nof_boundary();
            for t in 1..=nbnd {
                ch_boundary_elem(t);
            }
            ch_default_finish();
            ch_default_dirichlet();
            // `Norm = DefaultSolve()` + the two Level=4 reports.
            let _norm: c_double = ch_solve();
            ch_newton_update(iter);
            // `IF ( RelativeChange < NonlinTol ) EXIT`.
            if ch_converged() != 0 {
                break;
            }
        }
    }
}
