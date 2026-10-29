//! Port of `fem/src/modules/HelmholtzSolve.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/HelmholtzSolve.F90` — two F77-global
//! entries (`helmholtzsolver_`, `helmholtzsolver_init_`) with three
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those
//! two (proven by `nm -D` on the reference `.so`); the internals emit
//! only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/helmholtzsolve_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the Galerkin assembly, complex
//! admittance/convection algebra, `MATMUL`/`SUM` reductions, the
//! `CPUTime`-based reports) stays in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors HelmholtzSolve.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the entry SAVE block is allocated once
//!   at first-mesh sizing (`MaxElementDOFs`) and reallocated only when
//!   `Solver % MeshChanged`; the internal procedures' `SAVE Nodes` live
//!   as wrapper-side SAVE state with identical lifetimes (consumed
//!   before the next fill, as in the original).
//! - Allocation failure calls the same `Fatal('HelmholtzSolve', ...)`
//!   via the `STAT=istat` check, instead of aborting. The Rust side
//!   allocates nothing, so no `try_reserve` gate applies.
//! - Pure-query touches with no observable effect are replicated as
//!   touches where they could trap (the `PresPerm` association and the
//!   `Solver % Variable % Norm` read) and documented as skipped where
//!   they cannot (`GetSimulation()`, `GetEquation()`, the pre-loop
//!   `GetElementNOFNodes()` sizing read whose value is overwritten).
//! - `Model`, `dt`, and the transient flag are unused by the entry body
//!   itself, exactly as in the original (received and ignored), except
//!   for `Model % Simulation` in the `res: Frequency` list update.
//! - The `iter:` `Assembly: (s)` / `Solve: (s)` reports print raw
//!   `CPUTime()` values; the proof's norm class strips those
//!   wall-clock values while comparing the `iter:` structure.
//! - No Rust threads: the solver entries run serially, so the
//!   wrapper-side `SAVE` state needs no guard.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn hhz_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn hhz_ensure_alloc();
    fn hhz_pre() -> c_int;
    fn hhz_need_iface() -> c_int;
    fn hhz_iter_begin(iter: c_int);
    fn hhz_nactive() -> c_int;
    fn hhz_bulk(t: c_int);
    fn hhz_finish_bulk();
    fn hhz_nbnd() -> c_int;
    fn hhz_bnd(t: c_int);
    fn hhz_iface(t: c_int);
    fn hhz_assemble_done();
    fn hhz_dirichlet();
    fn hhz_solve(iter: c_int) -> c_double;
    fn hhz_converged() -> c_int;
    fn hhz_finish();
    fn hhz_init_begin(h_solver: *mut c_void);
    fn hhz_init_body();
}

/// `helmholtzsolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `HelmholtzSolver_init` (force `Linear System Complex`,
/// allocate the preconditioning matrix when a damp coefficient is given).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "helmholtzsolver_init_"]
pub unsafe extern "C" fn helmholtz_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        hhz_init_begin(solver);
        hhz_init_body();
    }
}

/// `helmholtzsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`helmholtz_solver_init_entry`]; additionally `model`
/// must be the Elmer `Model_t` address (used for the `res: Frequency`
/// list update, as in the original).
#[export_name = "helmholtzsolver_"]
pub unsafe extern "C" fn helmholtz_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        hhz_begin(model, solver);
        hhz_ensure_alloc();
        let niter = hhz_pre();
        let need_iface = hhz_need_iface();
        for iter in 1..=niter {
            hhz_iter_begin(iter);
            let na = hhz_nactive();
            for t in 1..=na {
                hhz_bulk(t);
            }
            hhz_finish_bulk();
            let nb = hhz_nbnd();
            for t in 1..=nb {
                hhz_bnd(t);
            }
            if need_iface != 0 {
                for t in 1..=nb {
                    hhz_iface(t);
                }
            }
            hhz_assemble_done();
            hhz_dirichlet();
            let _norm: c_double = hhz_solve(iter);
            if hhz_converged() != 0 {
                break;
            }
        }
        hhz_finish();
    }
}
