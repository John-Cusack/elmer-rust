//! Port of `fem/src/modules/StructuredProjectToPlane.F90` for lane T021-B,
//! at the solver-plugin boundary (whole `.so` → Rust cdylib with the same
//! dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/StructuredProjectToPlane.F90` — two F77-global
//! entries `structuredprojecttoplane_` and `structuredprojecttoplane_init_`
//! (the trailing underscore is part of the emitted dynamic symbol, proven by
//! `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/structuredprojecttoplane_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. All computation, order-sensitive
//! reductions (`SUM`, whole-array updates), formatted I/O (`WRITE(Message)`,
//! `Info`/`Warn`/`Fatal`) and the `SELECT CASE` operator dispatch stay
//! verbatim in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors the original statement for statement):
//! - `SAVE` semantics replicated explicitly: the whole SAVE list of the
//!   original (`Initialized`, `TopPointer`, `BotPointer`, `MidPointer`,
//!   `UpPointer`, `DownPointer`, `MaskPerm`, `Coord`, `TopPerm`, `MidPerm`,
//!   `BotPerm`, `UnitPerm`, `TopField`, `FieldOut`, `FieldIn`, `TopNodes`,
//!   `MidNodes`, `BotNodes`, `nsize`, `nnodes`, `NoLayers`, `MaskExist`,
//!   `MidLayerExists`, `Nodes`, `UnitVector`, `Visited`) lives as
//!   module-level SAVE state in the accessor with identical lifetimes; the
//!   `InvDGPerm` allocatable stays a per-call local (never SAVEd), exactly
//!   as in the original.
//! - Allocations (`TopPerm`, `BotPerm`, `MidPerm`, `TopField`, `UnitPerm`,
//!   `InvDGPerm`) stay in Fortran with the original bare `ALLOCATE`
//!   semantics (no `STAT=` in the original), so failure behavior is
//!   identical by construction; the Rust side allocates nothing.
//! - `Model`, `dt`, `Transient` are received and forwarded opaquely; the
//!   `GetSolverParams()` (no-arg, `CurrentModel % Solver`) resolution,
//!   `CPUTime`, `Info`, and `Message` traffic are verbatim in the wrappers.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn spa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn spa_init();
    fn spa_setup();
    fn spa_project();
}

/// Serial-entry guard: the Fortran side holds all (SAVE) state; the Rust
/// side only serializes entry execution.
static GUARD: Mutex<()> = Mutex::new(());

/// `structuredprojecttoplane_init_` — the solver init entry. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "structuredprojecttoplane_init_"]
pub unsafe extern "C" fn structuredprojecttoplane_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _g = GUARD.lock().unwrap();
    unsafe {
        spa_begin(model as *const c_void, solver as *const c_void);
        spa_init();
    }
}

/// `structuredprojecttoplane_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`.
///
/// # Safety
///
/// Same contract as [`structuredprojecttoplane_init_entry`].
#[export_name = "structuredprojecttoplane_"]
pub unsafe extern "C" fn structuredprojecttoplane_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _g = GUARD.lock().unwrap();
    unsafe {
        spa_begin(model as *const c_void, solver as *const c_void);
        spa_setup();
        spa_project();
    }
}
