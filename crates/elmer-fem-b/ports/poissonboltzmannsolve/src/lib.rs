//! Port of `fem/src/modules/PoissonBoltzmannSolve.F90` (Elmer @ the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/PoissonBoltzmannSolve.F90` — one
//! F77-global entry `poissonboltzmannsolve_` with four `CONTAINS`
//! internals. The emitted dynamic symbol is exactly
//! `poissonboltzmannsolve_` (proven by `nm -D` on the reference `.so`);
//! the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/poissonboltzmannsolve_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), transcendental evaluations (`SINH`/`COSH`/`EXP`), whole-array
//! updates, and formatted I/O (`WRITE`/`Info`) stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors PoissonBoltzmannSolve.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: first-mesh sizing (`MaxElementNodes`)
//!   lives wrapper-side with identical lifetimes; `NewtonLinearization`
//!   keeps its declaration-init (`.FALSE.`) latch — set `TRUE` on the
//!   same conditions, never reset per call. Rust owns only the
//!   `AllocationsDone` flag and the `DO iter` counter with the same
//!   `EXIT` shape.
//! - Allocation failure calls the same `Fatal('StatElecSolve', ...)` in
//!   Fortran (the `STAT=istat` checks), instead of aborting.
//! - Dead reads (`ForceVector`, entry `Norm`) are absent: they feed
//!   nothing downstream in the original either.
//! - `dt`/`TransientSimulation` are unused by the body, exactly as in
//!   the original (received and ignored).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.
//! - Nondimensional-mode `Permittivity` is left uninitialized exactly as
//!   in the original; the energy path that would read it stays off
//!   unless the inputs request it (same defined-inputs contract).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pbb_begin(h_model: *const c_void, h_solver: *const c_void);
    fn pbb_local_nodes() -> c_int;
    fn pbb_pre(nonlin: *mut c_int, newtontol: *mut c_double);
    fn pbb_first_alloc();
    fn pbb_refresh();
    fn pbb_iter_info(iter: c_int);
    fn pbb_newton_update(iter: c_int);
    fn pbb_default_initialize();
    fn pbb_bulk_assemble();
    fn pbb_neumann();
    fn pbb_finish();
    fn pbb_dirichlet();
    fn pbb_cpu() -> c_double;
    fn pbb_info_assembly(at: c_double);
    fn pbb_info_solve(st: c_double);
    fn pbb_solve_core();
    fn pbb_info_normblock(change: *mut c_double, conv: *mut c_int);
    fn pbb_flux_tail();
}

/// First-call flag only (`AllocationsDone`); all first-mesh buffers live
/// wrapper-side with the original's SAVE lifetimes.
struct PoissonBoltzmannState {
    allocations_done: bool,
}

impl PoissonBoltzmannState {
    const fn new() -> Self {
        Self {
            allocations_done: false,
        }
    }
}

static STATE: Mutex<PoissonBoltzmannState> = Mutex::new(PoissonBoltzmannState::new());

/// `poissonboltzmannsolve_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "poissonboltzmannsolve_"]
pub unsafe extern "C" fn poissonboltzmannsolve_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe {
        pbb_begin(model, solver);

        // IF (LocalNodes <= 0) RETURN
        if pbb_local_nodes() <= 0 {
            return;
        }

        let (mut nonlin, mut newtontol) = (0 as c_int, 0.0 as c_double);
        pbb_pre(&mut nonlin, &mut newtontol);

        // First-mesh permanent storage (bug-compatible: sized once).
        if !state.allocations_done {
            pbb_first_alloc();
            state.allocations_done = true;
        }

        pbb_refresh();

        // DO iter = 1, NonlinearIter ... IF (NonlinConverged == 1) EXIT
        let mut iter: c_int = 1;
        while iter <= nonlin {
            let at0 = pbb_cpu();
            pbb_default_initialize();
            pbb_iter_info(iter);
            // IF (iter > NewtonIter) NewtonLinearization = .TRUE.
            pbb_newton_update(iter);
            pbb_bulk_assemble();
            pbb_neumann();
            pbb_finish();
            pbb_dirichlet();
            let at1 = pbb_cpu();
            pbb_info_assembly(at1 - at0);

            let st0 = pbb_cpu();
            pbb_solve_core();
            let st1 = pbb_cpu();
            pbb_info_solve(st1 - st0);

            let (mut change, mut conv) = (0.0 as c_double, 0 as c_int);
            pbb_info_normblock(&mut change, &mut conv);
            let _ = change;
            if conv == 1 {
                break;
            }
            iter += 1;
        }

        pbb_flux_tail();
    }
}
