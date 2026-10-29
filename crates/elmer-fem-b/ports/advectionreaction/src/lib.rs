//! Port of `fem/src/modules/AdvectionReaction.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/AdvectionReaction.F90` — three F77-global
//! entries (`advectionreactionsolver_`, `advectionreactionsolver_init_`,
//! `advectionreactionsolver_post_`) with four `CONTAINS` internals. The
//! emitted dynamic symbols are exactly those three (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/advectionreaction_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), formatted I/O (`WRITE`/`Info`), whole-array updates, and the
//! DG jump/boundary algebra stay in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors AdvectionReaction.F90 statement for
//! statement):
//! - `SAVE` replicated exactly: `MASS`/`STIFF`/`LOAD`/`FORCE`/`Velo`/
//!   `MeshVelo`/`Gamma`/`UpperLimit`/`LowerLimit`, `AllocationsDone`,
//!   `DIM`, `VariableName` live as wrapper-side SAVE with identical
//!   lifetimes; the internal procedures' `SAVE Nodes` states (one per
//!   internal, two for the boundary one) likewise live wrapper-side.
//! - `Constants => GetConstants()` is a pure read whose value the
//!   original never uses; documented as skipped (no side effects).
//! - `RelativeChange = Solver % Variable % NonlinChange` is assigned but
//!   never read in the original; the read is replicated and discarded.
//! - Pre-loop `CorrectedUpperLimit`/`CorrectedLowerLimit` increments are
//!   unobservable (the original resets both to zero before the only
//!   prints); the wrapper keeps module counters reset at the same point.
//! - `CPUTIme()` is the same external as `CPUTime()` (Fortran is
//!   case-insensitive); the wrapper calls the same spelling.
//! - `FlowDOFs` is only set on the `computed` convection path; the
//!   wrapper zero-initialises it so the `none`/`constant` + mesh-velocity
//!   path (undefined in the original) subtracts nothing. The defined
//!   paths are statement-identical. Proof configs carry no mesh velocity,
//!   so the `ANY(MeshVelo /= 0)` gate is `.FALSE.` in both arms.
//! - `Model`, `dt` are received and ignored, exactly as in the original
//!   (`Default1stOrderTime` reads `dt` from the simulation itself).
//! - Allocation failure calls the same
//!   `Fatal('AdvectionReaction', ...)` via the `STAT=istat` check.
//! - No Rust threads: the solver entries run serially; the Rust side holds
//!   no shared state at all (all SAVE state is wrapper-side).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ara_begin(h_solver: *const c_void);
    fn ara_info_header();
    fn ara_faces() -> c_int;
    fn ara_num_boundary_elems() -> c_int;
    fn ara_alloc();
    fn ara_read_params(
        nonlin_min: *mut c_int,
        nonlin_max: *mut c_int,
        limit: *mut c_int,
    );
    fn ara_iter_start(iter: c_int, maxiter: c_int);
    fn ara_default_initialize();
    fn ara_nof_active() -> c_int;
    fn ara_bulk_assemble(t: c_int, active: c_int, limit: c_int);
    fn ara_first_order_time();
    fn ara_default_update();
    fn ara_face_force_zero();
    fn ara_face_assemble(t: c_int);
    fn ara_finish_bulk();
    fn ara_boundary_assemble(t: c_int) -> c_int;
    fn ara_finish_boundary();
    fn ara_finish_assembly();
    fn ara_solve_tail(iter: c_int, maxiter: c_int) -> c_int;
    fn ara_limit_pass(do_limit: c_int);
    fn ara_init_begin(h_solver: *const c_void);
    fn ara_init_body();
    fn ara_post_begin(h_solver: *const c_void);
    fn ara_post_body() -> c_int;
}

/// `advectionreactionsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are unused
/// by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read if non-null.
#[export_name = "advectionreactionsolver_"]
pub unsafe extern "C" fn advectionreaction_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let is_transient = !transient.is_null() && *transient != 0;
        ara_begin(solver);
        ara_info_header();
        let nfaces = ara_faces();
        ara_alloc();
        let mut nonlin_min: c_int = 0;
        let mut nonlin_max: c_int = 0;
        let mut limit: c_int = 0;
        ara_read_params(&mut nonlin_min, &mut nonlin_max, &mut limit);
        for iter in 1..=nonlin_max {
            ara_iter_start(iter, nonlin_max);
            ara_default_initialize();
            let active = ara_nof_active();
            for t in 1..=active {
                ara_bulk_assemble(t, active, limit);
                if is_transient {
                    ara_first_order_time();
                }
                ara_default_update();
            }
            ara_face_force_zero();
            for t in 1..=nfaces {
                ara_face_assemble(t);
            }
            ara_finish_bulk();
            let nbound = ara_num_boundary_elems();
            for t in 1..=nbound {
                if ara_boundary_assemble(t) == 0 {
                    continue;
                }
                if is_transient {
                    ara_first_order_time();
                }
                ara_default_update();
            }
            ara_finish_boundary();
            ara_finish_assembly();
            if ara_solve_tail(iter, nonlin_max) == 1 {
                break;
            }
        }
        ara_limit_pass(limit);
    }
}

/// `advectionreactionsolver_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `AdvectionReactionSolver_init` (DG flag, default
/// `Variable`, nodal-average export branch).
///
/// # Safety
///
/// Same contract as [`advectionreaction_solver_entry`].
#[export_name = "advectionreactionsolver_init_"]
pub unsafe extern "C" fn advectionreaction_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ara_init_begin(solver);
        ara_init_body();
    }
}

/// `advectionreactionsolver_post_` — the `_post` entry Elmer calls after
/// the solve when `PostSolver Active` is set (nodal averaging of the DG
/// field). Mirrors `AdvectionReactionSolver_post`.
///
/// # Safety
///
/// Same contract as [`advectionreaction_solver_entry`].
#[export_name = "advectionreactionsolver_post_"]
pub unsafe extern "C" fn advectionreaction_solver_post_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ara_post_begin(solver);
        let _code: c_int = ara_post_body();
    }
}
