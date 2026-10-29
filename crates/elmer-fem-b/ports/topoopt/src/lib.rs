//! Port of `fem/src/modules/TopoOpt.F90` (Elmer @ 9f6af2f85, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/TopoOpt.F90` — three F77-global entries
//! `topoopt_`, `topoopt_init_`, `topoopt_init0_` with `CONTAINS`
//! internals. The emitted dynamic symbols are exactly those three (proven
//! by `nm -D` on the reference `.so`); the internals emit only local
//! symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses; `dt` and
//! `TransientSimulation` cross by reference and are ignored, exactly as
//! in the original. Every component access and every Elmer library call
//! goes through the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/topoopt_acc.F90`, compiled against the built `.mod` files
//! of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`MATMUL`,
//! `SUM`), formatted I/O (`PRINT`/`WRITE`), whole-array updates, and the
//! UpdateDensities bisection loop stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors TopoOpt.F90 statement for statement):
//! - `SAVE` semantics replicated explicitly: the call-count latch
//!   `TimesVisited` lives in the guarded process-wide `STATE` (zeroed at
//!   load, incremented once per entry return, read for the first-touch
//!   dispatch and the `GetCycleMode` count, exactly like the original),
//!   while every other member of the original's `SAVE` list (filter
//!   matrices, elemental fields, scratch arrays, mesh/physics handles,
//!   filter strings, PDE-latch `DoMatrix`) lives as wrapper-side SAVE
//!   state with identical lifetimes. `GetCycleMode`'s period logic is
//!   transcribed verbatim, fed by the Rust latch.
//! - The mid-cycle early exit (`cMode == 1 .OR. cMode == 2` → `GOTO 1`,
//!   skipping reduction, report, filters, update, and norm assignment but
//!   still applying the stiffness multiplier and the visit count) is
//!   mirrored by the same predicate in the entry.
//! - Per-call `ValueList` scalars (`penal`, `volFrac`, `emin`, `efrac`)
//!   are recomputed from the same keywords on every call, exactly as the
//!   original's per-call locals; promoting them to module state is
//!   unobservable across wrapper boundaries within one call.
//! - `UpdateDensities`' `SAVE move, tol` latch, `ApplyPDEFilter`'s first-
//!   call matrix-assembly latch, and the internal procedures' `SAVE
//!   Nodes` live in Fortran with the original lifetimes.
//! - No Rust threads: the solver calls entries serially; the `Mutex` is
//!   only a `static`-sharing guard held for the counter read/increment,
//!   never across the solve (and nothing re-enters Rust mid-call).
//! - The Rust side allocates nothing, so there is no allocation-failure
//!   path to diverge; the Fortran side keeps the original `ALLOCATE`
//!   statements verbatim.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn toa_begin(model: *mut c_void, solver: *mut c_void);
    fn toa_init0();
    fn toa_init();
    fn toa_info_open();
    fn toa_first_touch();
    fn toa_info_kmat();
    fn toa_params();
    fn toa_recopy_xtopo();
    fn toa_cmode(nt: c_int) -> c_int;
    fn toa_zero_obj_dc();
    fn toa_objective_gradients();
    fn toa_mid_early();
    fn toa_parallel_obj();
    fn toa_report_objective();
    fn toa_pre_filter();
    fn toa_update();
    fn toa_post_filter();
    fn toa_post_phys_info();
    fn toa_assign_norm();
    fn toa_bw();
    fn toa_zero_level();
    fn toa_xmult();
}

/// First-call SAVE latch: mirrors `TimesVisited` (zero-initialized, and
/// zero-initialization implies `SAVE` in the original).
struct TopoState {
    times_visited: c_int,
}

impl TopoState {
    const fn new() -> Self {
        Self { times_visited: 0 }
    }
}

static STATE: Mutex<TopoState> = Mutex::new(TopoState::new());

/// `topoopt_init0_` — the phase-0 init entry: registers the primary
/// solver index and the adjoint export through the wrapper, which owns
/// the `ListAdd`/`ListGet` traffic verbatim.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are received and ignored, as in the original).
#[export_name = "topoopt_init0_"]
pub unsafe extern "C" fn topoopt_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    toa_begin(model, solver);
    toa_init0();
}

/// `topoopt_init_` — the init entry: declares the elemental `topo`
/// fields (plus `bw`/`levelset` on their branches) and the solver
/// variable through the wrapper, which owns the whole body verbatim.
///
/// # Safety
///
/// Same contract as [`topoopt_init0_entry`].
#[export_name = "topoopt_init_"]
pub unsafe extern "C" fn topoopt_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    toa_begin(model, solver);
    toa_init();
}

/// `topoopt_` — the solver entry. Argument order and by-reference travel
/// match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// Control flow mirrors the original body statement for statement:
/// first-touch sizing on `TimesVisited == 0`, per-call parameter reads,
/// cycle-mode dispatch (reset + gradients always; reduction, report,
/// pre/post filtering, update, norm, and extractions only off the
/// mid-cycle path), then the unconditional stiffness-multiplier tail and
/// the visit-count increment.
///
/// # Safety
///
/// Same contract as [`topoopt_init0_entry`].
#[export_name = "topoopt_"]
pub unsafe extern "C" fn topoopt_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    toa_begin(model, solver);
    // Short-lived latch read: the guard is released before any Elmer
    // call, so nothing is held across the solve.
    let visited = STATE.lock().unwrap().times_visited;
    toa_info_open();
    if visited == 0 {
        toa_first_touch();
    }
    toa_info_kmat();
    toa_params();
    if visited == 0 {
        toa_recopy_xtopo();
    }
    let mode = toa_cmode(visited);
    if mode == 0 || mode == 1 {
        toa_zero_obj_dc();
    }
    toa_objective_gradients();
    if mode == 1 || mode == 2 {
        // Mid of cycle, finishing early: skip straight to the label-1
        // tail (stiffness multiplier), exactly like the original GOTO.
        toa_mid_early();
    } else {
        toa_parallel_obj();
        toa_report_objective();
        toa_pre_filter();
        toa_update();
        toa_post_filter();
        toa_post_phys_info();
        toa_assign_norm();
        toa_bw();
        toa_zero_level();
    }
    toa_xmult();
    STATE.lock().unwrap().times_visited += 1;
}
