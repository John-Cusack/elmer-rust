//! Port of `fem/src/modules/ThermoElectricSolver.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ThermoElectricSolver.F90` — two F77-global
//! entries `thermoelectricsolver_` and `thermoelectricsolver_init_` with
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/thermoelectricsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`), formatted I/O (`WRITE`/`Message`) and whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ThermoElectricSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the per-element bubble history
//!   `bx` / `bxprev` (sized `2*MaxBDOFs*(bulk+bnd)` on first transient touch
//!   with `MaxBDOFs > 0`, zeroed) plus the `BubbleHistoryAllocated` / `stimestep`
//!   latch live in the guarded process-wide `STATE`, with identical
//!   lifetimes (`bxprev = bx` exactly when `GetTimestep()` moves).
//!   `NewtonLinearization` is SAVE too (its declaration carries `= .FALSE.`,
//!   which implies `SAVE`): once set it latches process-wide and is never
//!   reset, so the Newton branch stays on across timesteps exactly like the
//!   original (this is what makes timestep-2+ Newton continuations match).
//! - `LOAD` is intentionally left uninitialized when no body force is
//!   associated, exactly like the original's automatic array; every proof
//!   config provides a body force, so the proof runs on defined inputs.
//! - The Newton/Picard branch, the `Stabilize` flag (read once per solve),
//!   and the `GlobalBubbles` condensation dispatch match the original's
//!   predicates.
//! - `Model`, `dt`, `TransientSimulation` cross by reference exactly as in
//!   the original dummy list; `dt`'s value reaches `CondensatePTransient`
//!   unchanged, and the transient flag is only ever tested, never set.
//! - Allocation failure calls the same `Fatal('ThermoElectricSolver', ...)`
//!   via `try_reserve`, instead of aborting.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn tea_begin(h_model: *const c_void, h_solver: *const c_void);
    fn tea_init_fixup(transient: c_int);
    fn tea_mesh_maxbdofs() -> c_int;
    fn tea_mesh_bulk() -> c_int;
    fn tea_mesh_bnd() -> c_int;
    fn tea_get_timestep() -> c_int;
    fn tea_stabilize() -> c_int;
    fn tea_nl_params(
        ni: *mut c_int,
        ntol: *mut c_double,
        mi: *mut c_int,
        mtol: *mut c_double,
        mi_found: *mut c_int,
        mtol_found: *mut c_int,
    );
    fn tea_default_initialize();
    fn tea_nof_active() -> c_int;
    fn tea_default_finish_bulk();
    fn tea_nof_bnd() -> c_int;
    fn tea_default_finish();
    fn tea_default_dirichlet();
    fn tea_default_solve() -> c_double;
    fn tea_nonlin_change() -> c_double;
    fn tea_report(norm: c_double, rel: c_double);
    fn tea_fatal_alloc();
    fn tea_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn tea_elem_index() -> c_int;
    fn tea_global_bubbles() -> c_int;
    fn tea_lm_std(n: c_int, nd: c_int, nb: c_int, stab: c_int, newton: c_int, transient: c_int);
    #[allow(clippy::too_many_arguments)]
    fn tea_lm_bubble(
        n: c_int,
        nd: c_int,
        nb: c_int,
        stab: c_int,
        newton: c_int,
        dt: c_double,
        bprev: *const c_double,
        bcur: *mut c_double,
    );
    fn tea_begin_bnd(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn tea_bnd_assemble(n: c_int, nd: c_int, transient: c_int);
}

/// Bubble-history SAVE state: `bx` / `bxprev` sized
/// `2*MaxBDOFs*(NumberOfBulkElements+NumberOfBoundaryElements)` on first
/// transient touch, plus the allocation latch and last-seen timestep.
/// `stimestep = -1` is the pristine value, exactly like the original.
struct ThermoState {
    bx: Vec<c_double>,
    bxprev: Vec<c_double>,
    stride: usize,
    stimestep: c_int,
    bubble_allocated: bool,
    newton_linearization: bool,
}

impl ThermoState {
    const fn new() -> Self {
        Self {
            bx: Vec::new(),
            bxprev: Vec::new(),
            stride: 0,
            stimestep: -1,
            bubble_allocated: false,
            newton_linearization: false,
        }
    }
}

static STATE: Mutex<ThermoState> = Mutex::new(ThermoState::new());

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (it would overflow silently); growth only triggers where the original
/// is already out of bounds, so valid inputs never observe it.
fn ensure_bx(state: &mut ThermoState, need: usize) {
    if need <= state.bx.len() {
        return;
    }
    if state.bx.try_reserve(need - state.bx.len()).is_err()
        || state.bxprev.try_reserve(need - state.bxprev.len()).is_err()
    {
        unsafe { tea_fatal_alloc() };
        std::process::abort();
    }
    state.bx.resize(need, 0.0);
    state.bxprev.resize(need, 0.0);
}

/// `thermoelectricsolver_init_` — the solver init entry. Runs the
/// condensed-bubble min-iteration fixup through the wrapper, which owns
/// the whole `BLOCK` verbatim (including the `ListGetString` / hidden
/// length traffic that cannot cross the boundary).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "thermoelectricsolver_init_"]
pub unsafe extern "C" fn thermoelectric_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { tea_begin(model, solver) };
    let tflag: c_int = if transient.is_null() {
        0
    } else {
        unsafe { *transient }
    };
    unsafe { tea_init_fixup(tflag) };
}

/// `thermoelectricsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null.
#[export_name = "thermoelectricsolver_"]
pub unsafe extern "C" fn thermoelectric_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { tea_begin(model, solver) };

    let transient_flag: c_int = if transient.is_null() {
        0
    } else {
        unsafe { *transient }
    };
    let is_transient = transient_flag != 0;
    let dt_val: c_double = if dt.is_null() { 0.0 } else { unsafe { *dt } };

    let stabilize = unsafe { tea_stabilize() };

    // Per-element bubble history: allocated once, on first transient touch
    // with bubbles present (mirrors the SAVE allocation block).
    if is_transient && unsafe { tea_mesh_maxbdofs() } > 0 && !state.bubble_allocated {
        let stride = (2 * unsafe { tea_mesh_maxbdofs() }) as usize;
        let total =
            stride * ((unsafe { tea_mesh_bulk() } + unsafe { tea_mesh_bnd() }) as usize);
        if state.bx.try_reserve(total).is_err()
            || state.bxprev.try_reserve(total).is_err()
        {
            unsafe { tea_fatal_alloc() };
            std::process::abort();
        }
        state.bx.resize(total, 0.0);
        state.bxprev.resize(total, 0.0);
        state.stride = stride;
        state.bubble_allocated = true;
    }

    // Roll the previous-timestep history exactly on timestep change.
    if is_transient && unsafe { tea_mesh_maxbdofs() } > 0 {
        let ts = unsafe { tea_get_timestep() };
        if ts != state.stimestep {
            let len = state.bx.len();
            std::ptr::copy_nonoverlapping(state.bx.as_ptr(), state.bxprev.as_mut_ptr(), len);
            state.stimestep = ts;
        }
    }

    let (mut newton_iter, mut newton_tol) = (0 as c_int, 0.0 as c_double);

    let (mut nonlin_iter, mut nonlin_tol) = (0 as c_int, 0.0 as c_double);
    let (mut mi_found, mut mtol_found) = (0 as c_int, 0 as c_int);
    unsafe {
        tea_nl_params(
            &mut newton_iter,
            &mut newton_tol,
            &mut nonlin_iter,
            &mut nonlin_tol,
            &mut mi_found,
            &mut mtol_found,
        )
    };
    if mi_found == 0 {
        nonlin_iter = 1;
    }
    let _ = mtol_found;

    // `NewtonLinearization` is SAVE (declaration init): latched in state,
    // never reset here.
    for iter in 1..=nonlin_iter {
        unsafe { tea_default_initialize() };
        // BulkAssembly: serial element loop (the original's `!$omp`
        // pragmas compile to serial under the proof flags).
        let active = unsafe { tea_nof_active() };
        for t in 1..=active {
            let (mut n, mut nd0, mut nb) = (0, 0, 0);
            unsafe { tea_begin_active(t, &mut n, &mut nd0, &mut nb) };
            let nd = nd0 + nb;
            let newton: c_int = if state.newton_linearization { 1 } else { 0 };
            if is_transient && nb > 0 && unsafe { tea_global_bubbles() } == 0 {
                // Transient condensed-bubble path with per-element history.
                let stride = state.stride;
                let elem_index = unsafe { tea_elem_index() } as usize;
                let boff = stride * (elem_index - 1);
                ensure_bx(&mut state, boff + (2 * nb) as usize);
                let (bprev_ptr, bcur_ptr) = {
                    let s = boff;
                    let e = s + (2 * nb) as usize;
                    (
                        state.bxprev[s..e].as_ptr(),
                        state.bx[s..e].as_mut_ptr(),
                    )
                };
                unsafe {
                    tea_lm_bubble(n, nd, nb, stabilize, newton, dt_val, bprev_ptr, bcur_ptr)
                };
            } else {
                unsafe { tea_lm_std(n, nd, nb, stabilize, newton, transient_flag) };
            }
        }
        unsafe { tea_default_finish_bulk() };
        // BoundaryAssembly.
        let nbnd = unsafe { tea_nof_bnd() };
        for t in 1..=nbnd {
            let (mut n, mut nd, mut is_active) = (0, 0, 0);
            unsafe { tea_begin_bnd(t, &mut n, &mut nd, &mut is_active) };
            if is_active == 0 {
                continue;
            }
            unsafe { tea_bnd_assemble(n, nd, transient_flag) };
        }
        unsafe { tea_default_finish() };
        unsafe { tea_default_dirichlet() };

        let norm: c_double = unsafe { tea_default_solve() };
        let rel: c_double = unsafe { tea_nonlin_change() };

        unsafe { tea_report(norm, rel) };

        if rel < newton_tol || iter >= newton_iter {
            state.newton_linearization = true;
        }
        if rel < nonlin_tol {
            break;
        }
    }
}
