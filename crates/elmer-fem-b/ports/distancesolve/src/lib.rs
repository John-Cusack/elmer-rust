//! Port of `fem/src/modules/DistanceSolve.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/DistanceSolve.F90` — three F77-global
//! entries (`distancesolver_`, `distancesolver1_`, `distancesolver_init_`)
//! with `CONTAINS` internals. The emitted dynamic symbols are exactly
//! those three (proven by `nm -D` on the reference `.so`); the internals
//! emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/distancesolve_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`, `MAXVAL`, `MINVAL`, `SQRT`), the
//! `RANDOM_NUMBER` draw, `SortR`, formatted `WRITE`s, and whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors DistanceSolve.F90 statement for statement):
//! - `SAVE` semantics replicated: `STIFF`/`FORCE` allocated once at
//!   first-mesh sizing (`MaxElementDOFs`), never resized; the
//!   `LocalMatrix` `SAVE Nodes` lives as wrapper-side SAVE state with
//!   identical lifetime (consumed before the next fill, as in the
//!   original). The `DistanceSolver1` per-call arrays are allocated at
//!   entry and released at return with the same lifetimes as the
//!   original's locals.
//! - Allocation failure in the entry-1 init calls the same
//!   `Fatal('PoissonSolve', ...)` (the original's literal callee name,
//!   kept verbatim) via the `STAT=istat` check.
//! - Two deliberate output-neutral deviations: `DistanceSolver1` leaks
//!   its `bdist`/`bd`/`work` arrays per call in the original; the
//!   wrappers release every per-call array at return (and guard the
//!   allocation site) so repeated invocation cannot fail on
//!   already-allocated arrays. The Rust side passes the `DummyDistance`
//!   flag as 0/1; the wrapper branches exactly as the original `IF`.
//! - `Model`, `dt`, and the transient flag of `DistanceSolver1` are
//!   unused by the body, exactly as in the original (received and
//!   ignored); `Model` of `DistanceSolver` feeds the Noslip-wall loop.
//! - The parallel (`PEs > 1`) gather branch is transcribed into a
//!   wrapper (linking against the same serial MPI stubs as the reference
//!   `.so`); it never executes in the serial proof.
//! - No Rust threads: the solver entries run serially; the port holds no
//!   `static` state at all (all persistent state is wrapper-side SAVE).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    // -- DistanceSolver --
    fn ds_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn ds_info_using();
    fn ds_first_init();
    fn ds_read_params() -> c_int;
    fn ds_noslip();
    fn ds_active() -> c_int;
    fn ds_default_initialize();
    fn ds_elem(t: c_int);
    fn ds_finish_bulk();
    fn ds_finish();
    fn ds_dirichlet();
    fn ds_solve() -> c_double;
    fn ds_converged() -> c_int;
    fn ds_info_done();
    // -- DistanceSolver_init --
    fn ds_init_begin(h_solver: *mut c_void);
    fn ds_init_body();
    // -- DistanceSolver1 --
    fn d1_begin(h_solver: *mut c_void);
    fn d1_info_head();
    fn d1_mesh_var_info();
    fn d1_alloc();
    fn d1_read_dummy() -> c_int;
    fn d1_bulk_loop();
    fn d1_bc_loop();
    fn d1_pe_check() -> c_int;
    fn d1_serial_gather();
    fn d1_parallel_gather();
    fn d1_compute(dummy: c_int);
    fn d1_maxinfo();
    fn d1_finish();
}

/// `distancesolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `DistanceSolver_Init` (default nonlinear convergence
/// tolerance).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "distancesolver_init_"]
pub unsafe extern "C" fn distancesolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ds_init_begin(solver);
        ds_init_body();
    }
}

/// `distancesolver_` — the PDE-based distance solver entry. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "distancesolver_"]
pub unsafe extern "C" fn distancesolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ds_begin(model, solver);
        ds_info_using();
        ds_first_init();
        let maxiter = ds_read_params();
        ds_noslip();
        for _ in 0..maxiter {
            let active = ds_active();
            ds_default_initialize();
            for t in 1..=active {
                ds_elem(t);
            }
            ds_finish_bulk();
            ds_finish();
            ds_dirichlet();
            let _norm: c_double = ds_solve();
            if ds_converged() != 0 {
                break;
            }
        }
        ds_info_done();
    }
}

/// `distancesolver1_` — the geometric distance solver entry. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// Same contract as [`distancesolver_init_entry`].
#[export_name = "distancesolver1_"]
pub unsafe extern "C" fn distancesolver1_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        d1_begin(solver);
        d1_info_head();
        d1_mesh_var_info();
        d1_alloc();
        let dummy = d1_read_dummy();
        d1_bulk_loop();
        d1_bc_loop();
        match d1_pe_check() {
            1 => return,
            2 => d1_parallel_gather(),
            _ => d1_serial_gather(),
        }
        d1_compute(dummy);
        d1_maxinfo();
        d1_finish();
    }
}
