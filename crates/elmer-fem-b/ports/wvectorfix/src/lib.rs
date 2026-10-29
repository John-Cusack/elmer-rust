//! Port of `fem/src/modules/WVectorFix.F90` (Elmer @ c0b4bc6) for lane
//! T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib
//! with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/WVectorFix.F90` — two F77-global
//! entries `wvectorfix_` and `wvectorfix_init_` with two `CONTAINS`
//! internals (`LocalCorrMatrix`, `LocalCorrCurrent`). The emitted dynamic
//! symbols are exactly `wvectorfix_` + `wvectorfix_init_` (proven by `nm
//! -D` on the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/wvectorfix_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. The
//! two internal procedures are transcribed whole in the wrappers, so
//! order-sensitive reductions (`MATMUL`, `SUM`), the LU factor/solve,
//! whole-array updates and formatted I/O hold bit-identity by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors WVectorFix.F90 statement for statement):
//! - `SAVE` semantics replicated explicitly: the internals' `SAVE Nodes`
//!   live as wrapper-side SAVE state with identical lifetimes (consumed
//!   before the next fill, as in the original); the entry's `Dnode_index`
//!   dead read (ParEnv touch + last-element `NodeIndexes(1)`) is kept as
//!   a wrapper touch call with its value discarded, exactly as the
//!   original discards it.
//! - The per-flux loop keeps the original's branch structure: a missing
//!   `J Vector E` (`varindex == 2`) `CYCLE`s; a missing `W Vector E`
//!   fatals inside the wrapper with the original message. `Active =
//!   GetNOFActive()` is read once and reused by both element loops.
//! - `Norm = DefaultSolve()` is kept and discarded; `nsize`/`StiffMatrix`/
//!   `ForceVector`/`Mesh`/`SolVar` locals are pure handle reads with no
//!   downstream use, so they need no wrapper traffic.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies,
//!   exactly as in the original (received and ignored).
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn wva_begin(h_solver: *const c_void);
    fn wva_init_variable();
    fn wva_info_banner();
    fn wva_info_apply(v: c_int);
    fn wva_info_fixed(v: c_int);
    fn wva_find_flux(v: c_int) -> c_int;
    fn wva_default_initialize();
    fn wva_nof_active() -> c_int;
    fn wva_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn wva_corr_matrix(n: c_int, nd: c_int);
    fn wva_dead_node();
    fn wva_default_solve() -> c_double;
    fn wva_corr_current(n: c_int, nd: c_int);
}

/// Sharing guard for the wrapper-side SAVE state (entries run serially).
static GUARD: Mutex<()> = Mutex::new(());

/// `wvectorfix_init_` — solver init: sizes nothing, just registers the
/// `-nooutput WFixVarTmp` variable name on the solver's ValueList.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "wvectorfix_init_"]
pub unsafe extern "C" fn wvectorfix_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        wva_begin(solver);
        wva_init_variable();
    }
}

/// `wvectorfix_` — the solver entry. For each flux variable (`W Vector
/// E`, then `J Vector E` unless absent): assemble the elemental
/// divergence-correction potential system, solve it, and correct the
/// elemental field to be divergence free.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "wvectorfix_"]
pub unsafe extern "C" fn wvectorfix_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe { wva_begin(solver) };

    unsafe { wva_info_banner() };

    for v in 1..=2 as c_int {
        // Missing J Vector E (v == 2) CYCLEs; missing W Vector E fatals
        // inside the wrapper with the original message.
        if unsafe { wva_find_flux(v) } == 0 {
            continue;
        }
        unsafe { wva_info_apply(v) };

        unsafe { wva_default_initialize() };
        // Read once, reused by both element loops, as in the original.
        let active = unsafe { wva_nof_active() };
        for t in 1..=active {
            let (mut n, mut nd) = (0, 0);
            unsafe { wva_begin_active(t, &mut n, &mut nd) };
            unsafe { wva_corr_matrix(n, nd) };
        }

        // Dead Dirichlet-node line (value discarded, as original).
        unsafe { wva_dead_node() };

        // Kept and discarded, as in the original.
        let _norm: c_double = unsafe { wva_default_solve() };

        unsafe { wva_info_fixed(v) };
        for t in 1..=active {
            let (mut n, mut nd) = (0, 0);
            unsafe { wva_begin_active(t, &mut n, &mut nd) };
            unsafe { wva_corr_current(n, nd) };
        }
    }
}
