//! Port of `fem/src/modules/FilmCompressibility.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FilmCompressibility.F90` — two
//! F77-global entries (`filmcompressibility_`,
//! `filmcompressibility_init_`) with four `CONTAINS` internals. The
//! emitted dynamic symbols are exactly those two (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/filmcompressibility_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive FP
//! (`MATMUL`/`SUM`, the `WHERE` relaxation update, `MAX`), formatted
//! I/O (`WRITE`/`Info`), and whole-array updates stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors FilmCompressibility.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: all SAVE arrays and the
//!   `Nodes` state live as wrapper-side SAVE state with identical
//!   lifetimes (consumed before the next fill, as in the original);
//!   `nVisited`/`AllocationsDone` persist across invocations. The Rust
//!   side holds no state, so no `Mutex` is needed (the solver entry
//!   runs serially).
//! - Allocation failure calls the same
//!   `Fatal('FilmCompressibility', ...)` via the `STAT=istat` check,
//!   instead of aborting.
//! - The author's stale-`Found` read after the `Hydrostatic Pressure`
//!   lookup replicates automatically: one wrapper-side `Found` cell is
//!   threaded in the exact original call order.
//! - `m = SIZE(Var % Values)` is a dead read (never used) and is
//!   skipped; the entry-`dim` touch call is kept. `Norm =
//!   DefaultSolve()` is kept and discarded, exactly as in the original.
//! - `Model`, `dt`, and the transient flag travel to the wrappers
//!   (`Model` feeds the node/material/constant lookups); the transient
//!   flag only gates the author's transient-only `Fatal`.
//! - The `Levelset Initialize` path (`LocalLS` + `SetACInit`, early
//!   `RETURN`) is fully transcribed but never taken by the proof
//!   configs (the keyword is absent, so `LSInit` is false).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fcm_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn fcm_info_header();
    fn fcm_default_start();
    fn fcm_check_transient(tr: c_int);
    fn fcm_dim_touch();
    fn fcm_assoc_acinst();
    fn fcm_read_modes();
    fn fcm_sensmode() -> c_int;
    fn fcm_find_pressure();
    fn fcm_find_height();
    fn fcm_ranges20();
    fn fcm_sens_block();
    fn fcm_gravity();
    fn fcm_fixcheck() -> c_int;
    fn fcm_alloc();
    fn fcm_volreset();
    fn fcm_lsinit() -> c_int;
    fn fcm_nof_active() -> c_int;
    fn fcm_ls_begin();
    fn fcm_localls(t: c_int);
    fn fcm_ls_reduce();
    fn fcm_setacinit();
    fn fcm_dpmax_reset();
    fn fcm_localmaxima(t: c_int);
    fn fcm_max_reduce();
    fn fcm_max_report();
    fn fcm_read_diffusion();
    fn fcm_read_coeffs();
    fn fcm_default_initialize();
    fn fcm_assembly_elem(t: c_int);
    fn fcm_finish_boundary();
    fn fcm_finish_assembly();
    fn fcm_dirichlet();
    fn fcm_solve() -> c_double;
    fn fcm_default_finish();
    fn fcm_final_reduce();
    fn fcm_acinst_range();
    fn fcm_report();
    fn fcm_update();
    fn fcm_alldone_info();
    fn fcm_init_begin(h_solver: *mut c_void);
    fn fcm_init_body();
}

/// `filmcompressibility_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `FilmCompressibility_Init` (time derivative
/// order, `Variable = 'AC inst'`, exported `AC` variable).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "filmcompressibility_init_"]
pub unsafe extern "C" fn filmcompressibility_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fcm_init_begin(solver);
        fcm_init_body();
    }
}

/// `filmcompressibility_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `transient` must be a valid read (the Fortran `LOGICAL` flag).
#[export_name = "filmcompressibility_"]
pub unsafe extern "C" fn filmcompressibility_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        fcm_begin(model, solver);
        fcm_info_header();
        fcm_default_start();
        fcm_check_transient(i32::from(*transient != 0));
        fcm_dim_touch();
        fcm_assoc_acinst();
        fcm_read_modes();
        fcm_find_pressure();
        fcm_find_height();
        fcm_ranges20();
        if fcm_sensmode() != 0 {
            fcm_sens_block();
        }
        fcm_gravity();
        if fcm_fixcheck() != 0 {
            return;
        }
        fcm_alloc();
        fcm_volreset();
        if fcm_lsinit() != 0 {
            fcm_ls_begin();
            let n = fcm_nof_active();
            for t in 1..=n {
                fcm_localls(t);
            }
            fcm_ls_reduce();
            fcm_setacinit();
            return;
        }
        fcm_dpmax_reset();
        let n = fcm_nof_active();
        for t in 1..=n {
            fcm_localmaxima(t);
        }
        fcm_max_reduce();
        fcm_max_report();
        fcm_read_diffusion();
        fcm_read_coeffs();
        fcm_default_initialize();
        let n = fcm_nof_active();
        for t in 1..=n {
            fcm_assembly_elem(t);
        }
        fcm_finish_boundary();
        fcm_finish_assembly();
        fcm_dirichlet();
        let _norm: c_double = fcm_solve();
        fcm_default_finish();
        fcm_final_reduce();
        fcm_acinst_range();
        fcm_report();
        fcm_update();
        fcm_alldone_info();
    }
}
