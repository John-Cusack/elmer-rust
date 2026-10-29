//! Port of `fem/src/modules/TurbulenceModels/KomegaLegacy.F90` (Elmer @
//! c0b4bc6, the campaign fork) for lane T022-A, at the solver-plugin
//! boundary (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90` —
//! two F77-global entries (`komegalegacy_`, `komegalegacy_init_`) with
//! two `CONTAINS` internals. The emitted dynamic symbols are exactly
//! those two (proven by `nm -D` on the reference `.so`); the internals
//! emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries'
//! derived dummies (`Model`, `Solver`) cross as opaque addresses.
//! Every component access and every Elmer library call goes through
//! the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/komegalegacy_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. All floating-point work (the
//! LocalMatrix Galerkin assembly, the wall-distance MINVAL, the K/Omega
//! clip floors, the `CondensatePTransient` bubble history) and all
//! formatted I/O (`WRITE(Message,*)`/`Info`) stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors KomegaLegacy.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the first-mesh SAVE
//!   arrays (sized by `Mesh % MaxElementDOFs` and `Variable % DOFs`)
//!   and the per-element bubble history `bx`/`bxprev` live as
//!   wrapper-side SAVE state with identical lifetimes (allocated once,
//!   never resized; the timestep shift runs once per timestep, guarded
//!   by the same `GetTimestep() /= BubbleTimestep` check).
//! - Allocation failure calls the same `Fatal('KOmega', ...)` via the
//!   `STAT=istat` check, instead of aborting.
//! - Dead stores are documented skips, never executed: `at`/`at0`
//!   (CPUTime/RealTime assigned, never read), `NewtonTol`/`NewtonIter`
//!   (read, never used — still read, preserving the call sequence),
//!   `body_id` (set to -1, never read), `FlowSol`, `StiffMatrix`,
//!   `Stabilize`, `NewtonLinearization`, `eq_id`, `bf_id`, and the
//!   `KMax`/`EMax` MAXVAL reductions (pure, results never used) — like
//!   DivergenceSolver's `Visited`.
//! - `dt`/`TransientSimulation` cross by value; `dt` feeds only the two
//!   `CondensatePTransient` calls (`Default1stOrderTime` reads its own
//!   state), exactly as in the original. `Model` is received and stored
//!   (its components are read in wrappers, as in the original).
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn ko_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, transient: c_int);
    fn ko_prologue() -> c_int;
    fn ko_params();
    fn ko_noslip();
    fn ko_alloc();
    fn ko_niter() -> c_int;
    fn ko_iter_head(iter: c_int) -> c_int;
    fn ko_bulk_elem(t: c_int);
    fn ko_bulk_tail();
    fn ko_finish_assembly();
    fn ko_boundary();
    fn ko_dirichlet();
    fn ko_solve_tail() -> c_int;
    fn ko_init_begin(h_solver: *mut c_void, transient: c_int);
    fn ko_init_body();
}

static STATE: Mutex<()> = Mutex::new(());

/// `komegalegacy_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `KOmegaLegacy_Init` (p-bubble global-system default
/// and transient two-solve minimums).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read if non-null.
#[export_name = "komegalegacy_init_"]
pub unsafe extern "C" fn komegalegacy_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let _g = STATE.lock().unwrap();
        let tr = if transient.is_null() { 0 } else { *transient };
        ko_init_begin(solver, tr);
        ko_init_body();
    }
}

/// `komegalegacy_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`komegalegacy_init_entry`]; `dt` and `transient`
/// are read (the body consumes them).
#[export_name = "komegalegacy_"]
pub unsafe extern "C" fn komegalegacy_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let _g = STATE.lock().unwrap();
        let dtv = *dt;
        let tr = *transient;
        ko_begin(model, solver, dtv, tr);
        if ko_prologue() != 0 {
            return;
        }
        ko_params();
        ko_noslip();
        ko_alloc();
        let niter = ko_niter();
        for iter in 1..=niter {
            let noactive = ko_iter_head(iter);
            for t in 1..=noactive {
                ko_bulk_elem(t);
            }
            ko_bulk_tail();
            ko_finish_assembly();
            ko_boundary();
            ko_dirichlet();
            if ko_solve_tail() != 0 {
                break;
            }
        }
    }
}
