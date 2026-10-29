//! Port of `fem/src/modules/ReloadData.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ReloadData.F90` — one F77-global entry
//! `reloadsolution_` with no `CONTAINS` internals. The emitted dynamic
//! symbol is exactly `reloadsolution_` (proven by `nm -D` on the reference
//! `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/reloaddata_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. Order-
//! sensitive formatted I/O (`WRITE(Message, *)`) and the `LoadRestartFile`
//! / `Info` calls stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ReloadData.F90 statement for statement):
//! - `SAVE` semantics replicated: `Round` / `FirstStep` / `LastStep` /
//!   `Visit` / `StepsBetween` / `ContReading` live in the `static`
//!   `ReloadState`, initialized exactly once (`Round = 0`, `Visit = 0` at
//!   declaration; the rest on the first `Visit == 0` pass, as in the
//!   original). The per-invocation `RestartFile` local is staged in
//!   wrapper SAVE (filled before use, consumed before the next fill, as
//!   the original re-reads it every invocation).
//! - `GetInteger` returns `0` when unfound, so the `StartingStep /= 1`
//!   comparison sees `0` exactly as the original does (replicated by
//!   returning the real `GetInteger` value together with its `Found` flag).
//! - `MOD` with a positive span is Rust `%` (identical truncation); the
//!   `k == 0` / `k + FirstStep - 1` wrap logic is transcribed verbatim.
//! - `LoadRestartFile` runs with the `.TRUE.` optional only on the
//!   continuous arm (two wrapper spellings, matching the original's
//!   present/absent optional), never approximated with an explicit false.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn rda_begin(h_solver: *const c_void);
    fn rda_get_int(which: c_int, val: *mut c_int, found: *mut c_int);
    fn rda_get_cont(val: *mut c_int, found: *mut c_int);
    fn rda_info_cont_disabled();
    fn rda_has_restart(found: *mut c_int);
    fn rda_info_reading();
    fn rda_load(step: c_int, cont: c_int);
}

/// `which` selects the keyword: 1 `Reload Range Minimum`,
/// 2 `Reload Range Maximum`, 3 `Reload Starting Position`,
/// 4 `Reload Reading Intervals`. Returns the real `GetInteger` value
/// (`0` when unfound) plus its `Found` flag.
fn get_int(which: c_int) -> (c_int, bool) {
    let (mut val, mut found) = (0, 0);
    unsafe { rda_get_int(which, &mut val, &mut found) };
    (val, found != 0)
}

/// First-invocation SAVE state: `FirstStep`, `LastStep`, `StepsBetween`,
/// `ContReading`, plus the running `Round` / `Visit` counters.
struct ReloadState {
    round: c_int,
    first: c_int,
    last: c_int,
    visit: c_int,
    between: c_int,
    cont: bool,
}

impl ReloadState {
    const fn new() -> Self {
        Self {
            round: 0,
            first: 0,
            last: 0,
            visit: 0,
            between: 0,
            cont: false,
        }
    }
}

static STATE: Mutex<ReloadState> = Mutex::new(ReloadState::new());

/// `reloadsolution_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "reloadsolution_"]
pub unsafe extern "C" fn reload_solution_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { rda_begin(solver) };

    // IF ( Visit == 0 ) THEN — first-invocation defaults.
    if state.visit == 0 {
        let (v, got) = get_int(1);
        state.first = if got { v } else { 1 };

        let (v, got) = get_int(2);
        state.last = if got { v } else { 100000 };

        let (starting, got) = get_int(3);
        if got {
            state.round += starting - state.first;
        }

        let (v, got) = get_int(4);
        state.between = if got { v } else { 1 };

        let (mut val, mut found) = (0, 0);
        unsafe { rda_get_cont(&mut val, &mut found) };
        state.cont = if found != 0 { val != 0 } else { true };

        // Disabled unless the reload starts at the very first step.
        // `starting` is the raw GetInteger value (0 when unfound), so an
        // absent Reload Starting Position disables continuous reading,
        // exactly as in the original.
        if state.cont && (starting != 1 || state.first != 1) {
            state.cont = false;
            unsafe { rda_info_cont_disabled() };
        }
    }

    // IF ( MOD( Visit, StepsBetween ) == 0 ) THEN
    if state.visit % state.between == 0 {
        let mut found = 0;
        unsafe { rda_has_restart(&mut found) };

        if found != 0 {
            if state.visit == 0 {
                unsafe { rda_info_reading() };
            }

            state.round += 1;
            let span = 1 + state.last - state.first;
            let mut k = state.round % span;
            if k == 0 {
                k = state.last;
            } else {
                k += state.first - 1;
            }

            unsafe { rda_load(k, i32::from(state.cont)) };
        }
    }

    state.visit += 1;
}
