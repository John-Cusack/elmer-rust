//! Port of `fem/src/modules/VectorHelmholtzNodal.F90` (Elmer @ 9f6af2f,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/VectorHelmholtzNodal.F90` — two
//! F77-global entries, `vectorhelmholtznodal_init_` and
//! `vectorhelmholtznodal_`, with two `CONTAINS` internals. The emitted
//! dynamic symbols are exactly `vectorhelmholtznodal_init_` and
//! `vectorhelmholtznodal_` (proven by `nm -D` on the reference `.so`);
//! the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/vectorhelmholtznodal_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`), formatted I/O (`PRINT`), and whole-array/strided
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors VectorHelmholtzNodal.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the entry `Proj` matrix pointer, the
//!   internal-procedure SAVE (`Nodes`, `ValueHandle`s, `Basis`/
//!   `dBasisdx`/`STIFF`/`FORCE` first-mesh sizing, `GotPort`) live as
//!   wrapper-side SAVE state with identical lifetimes (consumed before
//!   the next fill, as in the original). The entry automatic `Norm` is
//!   zeroed at entry; every read follows the per-component
//!   `DefaultSolve` write, as in the original. The one piece of
//!   procedure-level state the BIND(C) split must carry across calls
//!   (the boundary `Found`, including the Parent-guarded stale read) is
//!   a wrapper-side mirror primed by the head wrapper, exactly the
//!   original single-procedure lifetime.
//! - `Model`, `dt`, `TransientSimulation` are received and ignored,
//!   exactly as in the original (the main body never touches `dt` or
//!   the transient flag; the entry dummy `Solver` travels via the
//!   wrapper-captured handle).
//! - The dead `Absorb` head fetch (re-fetched per IP in the original)
//!   is kept verbatim, and the unused `F`, `C`, `Ext` BC declarations
//!   are kept, so the call sequence matches the original statement
//!   for statement.
//! - No Rust threads: the solver entries run serially, so no `static`
//!   sharing guard is needed — the Rust side holds no cross-invocation
//!   state at all (every loop counter and branch predicate is a call
//!   local); all persistent state lives in the Fortran wrappers.
//! - The `im` unit-imaginary factor is the `VectorHelmholtzUtils`
//!   parameter, exactly as in the original unit.
//! - Allocation failure calls the same `Fatal('VectorHelmholtzNodal',
//!   ...)` in the wrappers (the `STAT=allocstat` check), never
//!   aborting directly.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn vhn_init_begin(model: *mut c_void, solver: *mut c_void);
    fn vhn_init_precuse() -> c_int;
    fn vhn_init_monolithic() -> c_int;
    fn vhn_init_curlcurl() -> c_int;
    fn vhn_init_set_monolithic();
    fn vhn_init_var_prec();
    fn vhn_init_var_field();
    fn vhn_init_var_etmp();
    fn vhn_init_export_prec();
    fn vhn_init_export_field();
    fn vhn_init_complex_time();
    fn vhn_init_findeigen() -> c_int;
    fn vhn_init_eigensolver();
    fn vhn_begin(model: *mut c_void, solver: *mut c_void);
    fn vhn_info_header();
    fn vhn_dim() -> c_int;
    fn vhn_mesh_params();
    fn vhn_dofs_branch(compn: *mut c_int, segregated: *mut c_int);
    fn vhn_field_setup() -> c_int;
    fn vhn_follow() -> c_int;
    fn vhn_follow_apply();
    fn vhn_default_start();
    fn vhn_maxiter() -> c_int;
    fn vhn_relorder() -> c_int;
    fn vhn_init_stuff();
    fn vhn_proj();
    fn vhn_info_compi(compi: c_int);
    fn vhn_default_initialize();
    fn vhn_nof_active() -> c_int;
    fn vhn_nof_boundary() -> c_int;
    fn vhn_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn vhn_begin_boundary(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn vhn_lm_setup(h: *mut c_int);
    fn vhn_lm_nip() -> c_int;
    fn vhn_lm_ip(t: c_int, compi: c_int, nd: c_int);
    fn vhn_lm_update(compi: c_int, nd: c_int);
    fn vhn_finish_bulk();
    fn vhn_apply_residual();
    fn vhn_info_range();
    fn vhn_lbc_head(h: *mut c_int, proceed: *mut c_int, nip: *mut c_int);
    fn vhn_lbc_ip(t: c_int, nd: c_int);
    fn vhn_lbc_update(compi: c_int, nd: c_int);
    fn vhn_finish_bnd();
    fn vhn_finish_assem();
    fn vhn_dirichlet(compi: c_int);
    fn vhn_gather(compi: c_int);
    fn vhn_solve(compi: c_int);
    fn vhn_scatter(compi: c_int);
    fn vhn_project();
    fn vhn_finish();
    fn vhn_final_norm();
    fn vhn_info_done();
}

/// `vectorhelmholtznodal_init_` — initialization of the primary solver.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt`/transient are
/// received and ignored, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtznodal_init_"]
pub unsafe extern "C" fn vectorhelmholtz_nodal_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { vhn_init_begin(model, solver) };
    let precuse = unsafe { vhn_init_precuse() };
    let mut monolithic = unsafe { vhn_init_monolithic() };
    let curlcurl = unsafe { vhn_init_curlcurl() };
    // curl-curl form forces the monolithic branch, and the flag with it.
    if curlcurl != 0 && monolithic == 0 {
        unsafe { vhn_init_set_monolithic() };
        monolithic = 1;
    }
    if monolithic != 0 {
        if precuse != 0 {
            unsafe { vhn_init_var_prec() };
        } else {
            unsafe { vhn_init_var_field() };
        }
    } else {
        unsafe { vhn_init_var_etmp() };
        if precuse != 0 {
            unsafe { vhn_init_export_prec() };
        } else {
            unsafe { vhn_init_export_field() };
        }
    }
    unsafe { vhn_init_complex_time() };
    if unsafe { vhn_init_findeigen() } != 0 {
        unsafe { vhn_init_eigensolver() };
    }
}

/// `vectorhelmholtznodal_` — the nodal vector-Helmholtz solver entry.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`/`dt`/transient
/// are received and ignored by the body, exactly as in the original
/// (the model handle is still captured for `FollowCurvedBoundary` and
/// the constants lookup).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "vectorhelmholtznodal_"]
pub unsafe extern "C" fn vectorhelmholtz_nodal_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { vhn_begin(model, solver) };
    unsafe { vhn_info_header() };
    let _dim: c_int = unsafe { vhn_dim() };
    unsafe { vhn_mesh_params() };
    let (mut compn, mut segregated) = (0, 0);
    unsafe { vhn_dofs_branch(&mut compn, &mut segregated) };
    let precuse = unsafe { vhn_field_setup() };
    if unsafe { vhn_follow() } != 0 {
        unsafe { vhn_follow_apply() };
    }
    unsafe { vhn_default_start() };
    // Assigned and never read in the original either (there is no
    // nonlinear iteration loop over maxiter); the query is kept so the
    // call sequence matches.
    let _maxiter = unsafe { vhn_maxiter() };
    let _relorder: c_int = unsafe { vhn_relorder() };
    unsafe { vhn_init_stuff() };
    unsafe { vhn_proj() };
    for compi in 1..=compn {
        if segregated != 0 {
            unsafe { vhn_info_compi(compi) };
        }
        unsafe { vhn_default_initialize() };
        let active = unsafe { vhn_nof_active() };
        let mut init_handles: c_int = 1;
        for t in 1..=active {
            let (mut _n, mut nd, mut nb) = (0, 0, 0);
            unsafe { vhn_begin_active(t, &mut _n, &mut nd, &mut nb) };
            let ndnb = nd + nb;
            unsafe { vhn_lm_setup(&mut init_handles) };
            let nip = unsafe { vhn_lm_nip() };
            for ip in 1..=nip {
                unsafe { vhn_lm_ip(ip, compi, ndnb) };
            }
            unsafe { vhn_lm_update(compi, ndnb) };
        }
        unsafe { vhn_finish_bulk() };
        // Residual back-transformation when used as a preconditioner
        // (a no-op for the component-wise form, as in the original).
        if precuse != 0 {
            unsafe { vhn_apply_residual() };
        }
        unsafe { vhn_info_range() };
        let nbnd = unsafe { vhn_nof_boundary() };
        let mut init_handles_bc: c_int = 1;
        for t in 1..=nbnd {
            let (mut _n, mut nd, mut active_el) = (0, 0, 0);
            unsafe { vhn_begin_boundary(t, &mut _n, &mut nd, &mut active_el) };
            if active_el == 0 {
                continue;
            }
            let (mut proceed, mut nip) = (0, 0);
            unsafe { vhn_lbc_head(&mut init_handles_bc, &mut proceed, &mut nip) };
            if proceed == 0 {
                continue;
            }
            for ip in 1..=nip {
                unsafe { vhn_lbc_ip(ip, nd) };
            }
            unsafe { vhn_lbc_update(compi, nd) };
        }
        unsafe {
            vhn_finish_bnd();
            vhn_finish_assem();
        }
        unsafe { vhn_dirichlet(compi) };
        // Segregated component shuffle around the linear solve.
        if segregated != 0 {
            unsafe { vhn_gather(compi) };
        }
        unsafe { vhn_solve(compi) };
        if segregated != 0 {
            unsafe { vhn_scatter(compi) };
        }
    }
    if precuse != 0 {
        unsafe { vhn_project() };
    }
    unsafe { vhn_finish() };
    if segregated != 0 {
        unsafe { vhn_final_norm() };
    }
    unsafe { vhn_info_done() };
}
