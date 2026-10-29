//! Port of `fem/src/modules/FilmFlowSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FilmFlowSolver.F90` — three F77-global
//! entries (`filmflowsolver_`, `filmflowsolver_init_`,
//! `filmflowsolver_init0_`); the `CONTAINS` internals
//! (`FrictionLawPraks`, `LocalBulkMatrix`, `LocalBoundaryMatrix`) emit
//! only local symbols. The emitted dynamic symbols are exactly those
//! three (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses; `dt` and the
//! transient flag cross by value after one by-reference read each,
//! matching the Fortran dummy travel. Every component access and every
//! Elmer library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/filmflowsolver_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself
//! contains no layout assumption and no pinned offset. All
//! floating-point work (the bulk/boundary Galerkin assembly, `MATMUL`/
//! `SUM` reductions, `WHERE` clamps) and all formatted/whole-array I/O
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FilmFlowSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: all entry-level SAVE arrays and flags
//!   plus the internals' `SAVE Nodes` / `SAVE Visited, MinPres, MinSpeed`
//!   live as wrapper-side SAVE state with identical lifetimes and the
//!   same first-mesh sizing (`(mdim+1)*(MaxElementDOFs+4*BDOFs)`).
//!   `FirstRound` starts `.TRUE.` and latches `.FALSE.` exactly as in
//!   the original; `itime0` starts at -1; `DoneWeight` latches once.
//! - Allocation failure calls the same `Fatal('FilmFlowSolver', ...)`
//!   via the `STAT=istat` check, instead of aborting.
//! - The nonlinear-sweep structure (`100 CONTINUE` + `DO iter`,
//!   `NSCondensate` vs `Default1stOrderTime` vs `DefaultUpdateEquations`
//!   sequencing, boundary assembly, `DefaultSensitivity` re-entry with
//!   `maxiter = 1`, heating pre/post blocks, `InfoActive`-gated
//!   diagnostics) is transcribed in the entry below; every per-element
//!   and per-phase body is a wrapper call in the original order.
//! - `FrictionLawPraks`'s `SAVE Visited` flag is set-but-never-read in
//!   the original and is not replicated (no observable state).
//! - The original's non-`SAVE` locals that are only set on the first
//!   visit (e.g. `SurfAc`) keep their first-visit meaning across
//!   re-entry; they feed only commented-out code, so the defined path
//!   is unaffected.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn ffw_init0_begin(h_model: *const c_void, h_solver: *const c_void);
    fn ffw_init0_body();
    fn ffw_init_begin(h_model: *const c_void, h_solver: *const c_void);
    fn ffw_init_body();
    fn ffw_begin(
        h_model: *const c_void,
        h_solver: *const c_void,
        dt: c_double,
        transient: c_int,
        code: *mut c_int,
    );

    fn ffw_heating_pre();
    fn ffw_ac_begin();
    fn ffw_maxiter() -> c_int;
    fn ffw_iter_begin() -> c_int;
    fn ffw_bulk_elem(t: c_int);
    fn ffw_finish_bulk();
    fn ffw_ac_report();
    fn ffw_nbnd() -> c_int;
    fn ffw_boundary_elem(t: c_int);
    fn ffw_iter_end() -> c_double;
    fn ffw_nonlin_converged() -> c_int;
    fn ffw_firstround_done();
    fn ffw_sensitivity() -> c_int;
    fn ffw_finish();
    fn ffw_heating_post();
    fn ffw_diagnostics();
}

static GUARD: Mutex<()> = Mutex::new(());

/// Run one nonlinear sweep: `DO iter = 1, maxiter` with bulk assembly,
/// FSI-flux report, boundary assembly, solve, and the convergence /
/// `FirstRound` latch. Mirrors the `100 CONTINUE` target region.
unsafe fn nonlinear_sweep(maxiter: c_int) {
    for _iter in 1..=maxiter {
        let active = unsafe { ffw_iter_begin() };
        for t in 1..=active {
            unsafe { ffw_bulk_elem(t) };
        }
        unsafe { ffw_finish_bulk() };
        unsafe { ffw_ac_report() };
        let nbnd = unsafe { ffw_nbnd() };
        for t in 1..=nbnd {
            unsafe { ffw_boundary_elem(t) };
        }
        let _norm: c_double = unsafe { ffw_iter_end() };
        if unsafe { ffw_nonlin_converged() } != 0 {
            break;
        }
        unsafe { ffw_firstround_done() };
    }
}

/// `filmflowsolver_init0_` — the `_Init0` entry Elmer calls at load
/// time. Mirrors `FilmFlowSolver_Init0` (time derivative order, `p`
/// element string honoring `Serendipity P Elements`).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "filmflowsolver_init0_"]
pub unsafe extern "C" fn filmflowsolver_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        ffw_init0_begin(_model as *const c_void, solver as *const c_void);
        ffw_init0_body();
    }
}

/// `filmflowsolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `FilmFlowSolver_Init` (`Model Dimension`-driven
/// `Variable`, norm DOFs, bubble/global-mass defaults, heating and
/// sensitivity exports).
///
/// # Safety
///
/// Same contract as [`filmflowsolver_init0_entry`].
#[export_name = "filmflowsolver_init_"]
pub unsafe extern "C" fn filmflowsolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        ffw_init_begin(_model as *const c_void, solver as *const c_void);
        ffw_init_body();
    }
}

/// `filmflowsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads (their values drive the
/// transient term and the heating correction, as in the original).
#[export_name = "filmflowsolver_"]
pub unsafe extern "C" fn filmflowsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    let dt_val: c_double = unsafe { *dt };
    let transient_val: c_int = unsafe { *transient };
    // Skip-First-Solution early RETURN path.
    let mut code: c_int = 0;
    unsafe {
        ffw_begin(
            _model as *const c_void,
            solver as *const c_void,
            dt_val,
            transient_val,
            &mut code,
        );
    }
    if code != 0 {
        return;
    }
    unsafe {
        ffw_heating_pre();
        ffw_ac_begin();
    }
    let mut maxiter = unsafe { ffw_maxiter() };
    // `IF (DefaultSensitivity()) THEN maxiter = 1; GOTO 100` — the goto
    // targets the sweep head, so a second pass re-checks the flag.
    loop {
        unsafe { nonlinear_sweep(maxiter) };
        if unsafe { ffw_sensitivity() } == 0 {
            break;
        }
        maxiter = 1;
    }
    unsafe {
        ffw_finish();
        ffw_heating_post();
        ffw_diagnostics();
    }
}
