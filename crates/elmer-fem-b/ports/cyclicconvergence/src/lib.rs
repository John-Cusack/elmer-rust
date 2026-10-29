//! Port of `fem/src/modules/CyclicConvergence.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/CyclicConvergence.F90` — two F77-global
//! entries (`cyclicconvergence_`, `cyclicconvergence_init_`) with one
//! `CONTAINS` internal. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); the internal emits only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/cyclicconvergence_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. All floating-point work (the
//! `FLOOR` cycle index, `SUM`/`SQRT`/`MINVAL`/`MAXVAL` reductions, the
//! error ratios), all formatted I/O (`Message`, the `ES12.3` entry
//! report, the `7ES15.6` convergence file), and the whole-array updates
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors CyclicConvergence.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the entry's `SAVE Visited, valTable,
//!   tsize, prevCycle, prevVal, prevDev, prevminv, prevmaxv, prevrangev,
//!   iCycle, ConvCycle, Converged` state and the internal procedure's
//!   `SAVE ConvVisited` live as wrapper-side SAVE with identical
//!   lifetimes; the growth doubling (`2*tsize`, zero-filled tail) and the
//!   reset tail (`iCycle = 0`, `prevCycle = nCycle`) are transcribed
//!   exactly. The internal `ALLOCATE`s carry no `STAT=` check in the
//!   original, so neither do the wrappers.
//! - `Model`, `dt` (past the `1.0e-3*dt` cycle-index offset), and the
//!   transient flag cross the boundary and are used exactly as in the
//!   original; nothing else is added or skipped.
//! - No Rust threads: the solver entries run serially, so no `Mutex` is
//!   needed — all shared state lives in the Fortran wrapper module.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn cca_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn cca_info_start();
    fn cca_percycle(dt: c_double);
    fn cca_first_touch();
    fn cca_resolve(val: *mut c_double);
    fn cca_push(val: c_double);
    fn cca_same_cycle() -> c_int;
    fn cca_stats();
    fn cca_min_pending() -> c_int;
    fn cca_is_converged() -> c_int;
    fn cca_evaluate();
    fn cca_produce();
    fn cca_reset();
    fn cci_init_begin(h_solver: *mut c_void);
    fn cci_init_body();
}

/// `cyclicconvergence_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `CyclicConvergence_Init` (the `Equation` and
/// `Variable` list entries).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "cyclicconvergence_init_"]
pub unsafe extern "C" fn cyclicconvergence_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        cci_init_begin(solver);
        cci_init_body();
    }
}

/// `cyclicconvergence_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`cyclicconvergence_init_entry`]; additionally `dt`
/// must be a valid read (only the `1.0e-3*dt` offset uses it).
#[export_name = "cyclicconvergence_"]
pub unsafe extern "C" fn cyclicconvergence_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        let dtv: c_double = *dt;
        cca_begin(model, solver);
        cca_info_start();
        cca_percycle(dtv);
        cca_first_touch();
        let mut val: c_double = 0.0;
        cca_resolve(&mut val);
        cca_push(val);
        // Same cycle, nothing to do!
        if cca_same_cycle() != 0 {
            return;
        }
        cca_stats();
        if cca_min_pending() != 0 {
            cca_reset();
            return;
        }
        if cca_is_converged() == 0 {
            cca_evaluate();
        } else {
            cca_produce();
        }
        // Restart the counter
        cca_reset();
    }
}
