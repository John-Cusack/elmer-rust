//! Port of `fem/src/modules/FreeSurfaceSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FreeSurfaceSolver.F90` — six F77-global
//! entries (`freesurfacesolver_`, `freesurfacesolver_rotinit_`,
//! `freesurfacetomeshupdate_`, `freesurfacetomeshupdate1_`,
//! `freesurfacetomeshupdate2_`, `freesurfacetomeshupdate3_`; proven by
//! `nm -D` on the reference `.so`) with one `CONTAINS` internal
//! (`LocalMatrix`). The emitted dynamic symbols are exactly those six;
//! the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/freesurfacesolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`, `MAXVAL`/`MINVAL`), formatted `WRITE` + `Info` /
//! `Warn` / `Fatal` reports, the element assembly, and the whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FreeSurfaceSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: variable-name/`PrevSize` reallocation
//!   tracking, allocation-done flag, all STIFF/MASS/... work arrays, the
//!   `SubstantialSurface`/`NormalFlux`/`Bubbles`/rotation state, limiter
//!   arrays, `OldValues`/`OldRHS`, and the mesh-update function's
//!   `FirstTime`/`t0` all live as wrapper-side SAVE state with identical
//!   lifetimes. The Rust side holds no cross-call state (no `Mutex`
//!   needed); the nonlinear-iteration and element loops re-query wrapper
//!   flags exactly where the original re-tests its locals.
//! - Dead locals (`LocalNodes`, `smallestpossiblenumber`, `firstTime`,
//!   `stat`, `UseBodyForce`, `maxdh_comm`, `LocalBottom`, `DIM` touch in
//!   the init routine) are assigned-never-read in the original and are
//!   documented as skipped; the `Pure` coordinate query they came from is
//!   kept wherever it feeds a live branch.
//! - The four `FreeSurfaceToMeshUpdate` functions share one wrapper core,
//!   exactly as 1/2/3 call the common internal function in the original;
//!   the core keeps the original's uninitialized-`SAVE` `FirstTime`/`t0`
//!   (zero-fill in the serial build, so the `FirstTime` branch is never
//!   taken) and the stale-`Found` read inside that never-taken branch.
//!   The Rust trampolines only forward the by-reference dummies and the
//!   axis selector and return the `REAL(KIND=dp)` result, matching the
//!   F77 call convention of the user-function loader.
//! - `Model`, `dt`, and the transient flag are unused by the init-routine
//!   and mesh-update bodies, exactly as in the original (received and
//!   ignored); `dt` is likewise unused by the main body, while the
//!   transient flag drives the time-discretization branches.
//! - Allocation failure calls the same `Fatal(..., 'Memory allocation
//!   error ...')` via the `STAT=istat` check, never aborting directly.
//! - No Rust threads: the solver entries run serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fsr_begin(h_model: *const c_void, h_solver: *const c_void);
    fn fsr_params();
    fn fsr_bind();
    fn fsr_num_nodes() -> c_int;
    fn fsr_ref_perm(i: c_int) -> c_int;
    fn fsr_rot_node(i: c_int);
    fn fss_mu(
        h_model: *const c_void,
        h_node: *const c_int,
        arr: *const c_double,
        axis: c_int,
    ) -> c_double;
    fn fss_begin(h_model: *const c_void, h_solver: *const c_void);
    fn fss_info_header();
    fn fss_realloc() -> c_int;
    fn fss_free_count() -> c_int;
    fn fss_warn_nofree();
    fn fss_read();
    fn fss_apply() -> c_int;
    fn fss_need_old() -> c_int;
    fn fss_nliter() -> c_int;
    fn fss_rot_check() -> c_int;
    fn fss_info_dims();
    fn fss_allocate();
    fn fss_bind_residual();
    fn fss_maybe_reset();
    fn fss_iter_begin(it: c_int, nlit: c_int);
    fn fss_nactive() -> c_int;
    fn fss_elem_begin(t: c_int, n: *mut c_int, skip: *mut c_int);
    fn fss_elem_limits(n: c_int);
    fn fss_use_local() -> c_int;
    fn fss_elem_maxdisp(n: c_int);
    fn fss_elem_velo(n: c_int, transient: c_int);
    fn fss_elem_meshvelo();
    fn fss_elem_source(n: c_int);
    fn fss_elem_old(n: c_int);
    fn fss_local_matrix(ncoord: c_int, transient: c_int);
    fn fss_elem_time(transient: c_int);
    fn fss_elem_update(n: c_int);
    fn fss_finish();
    fn fss_limits_pre();
    fn fss_info_assdone();
    fn fss_solve(relchange: *mut c_double) -> c_double;
    fn fss_residual(it: c_int);
    fn fss_relax();
    fn fss_iter_tail(it: c_int) -> c_int;
}

/// `freesurfacesolver_rotinit_` — the rotated-free-surface init routine.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored, exactly as in the original).
#[export_name = "freesurfacesolver_rotinit_"]
pub unsafe extern "C" fn freesurface_solver_rotinit_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fsr_begin(model, solver);
        fsr_params();
        fsr_bind();
        let nn = fsr_num_nodes();
        for i in 1..=nn {
            if fsr_ref_perm(i) == 0 {
                continue;
            }
            fsr_rot_node(i);
        }
    }
}

/// Shared trampoline core for the `FreeSurfaceToMeshUpdate` user
/// functions. `axis` selects the mesh-update component (1, 2, or 3);
/// the general entry forwards the caller's axis.
///
/// # Safety
///
/// Same contract as [`freesurface_solver_rotinit_entry`]; `nodenumber`
/// and `inarray` must be valid reads if non-null, `inarray` holding at
/// least two values.
unsafe fn mesh_update_trampoline(
    model: *mut c_void,
    nodenumber: *const c_int,
    inarray: *const c_double,
    axis: c_int,
) -> c_double {
    unsafe { fss_mu(model, nodenumber, inarray, axis) }
}

/// `freesurfacetomeshupdate1_` — x-component mesh update user function.
///
/// # Safety
///
/// See [`mesh_update_trampoline`].
#[export_name = "freesurfacetomeshupdate1_"]
pub unsafe extern "C" fn freesurface_to_mesh_update1_entry(
    model: *mut c_void,
    nodenumber: *const c_int,
    inarray: *const c_double,
) -> c_double {
    unsafe { mesh_update_trampoline(model, nodenumber, inarray, 1) }
}

/// `freesurfacetomeshupdate2_` — y-component mesh update user function.
///
/// # Safety
///
/// See [`mesh_update_trampoline`].
#[export_name = "freesurfacetomeshupdate2_"]
pub unsafe extern "C" fn freesurface_to_mesh_update2_entry(
    model: *mut c_void,
    nodenumber: *const c_int,
    inarray: *const c_double,
) -> c_double {
    unsafe { mesh_update_trampoline(model, nodenumber, inarray, 2) }
}

/// `freesurfacetomeshupdate3_` — z-component mesh update user function.
///
/// # Safety
///
/// See [`mesh_update_trampoline`].
#[export_name = "freesurfacetomeshupdate3_"]
pub unsafe extern "C" fn freesurface_to_mesh_update3_entry(
    model: *mut c_void,
    nodenumber: *const c_int,
    inarray: *const c_double,
) -> c_double {
    unsafe { mesh_update_trampoline(model, nodenumber, inarray, 3) }
}

/// `freesurfacetomeshupdate_` — general-axis mesh update user function.
///
/// # Safety
///
/// See [`mesh_update_trampoline`]; `axis` must be a valid read if
/// non-null.
#[export_name = "freesurfacetomeshupdate_"]
pub unsafe extern "C" fn freesurface_to_mesh_update_entry(
    model: *mut c_void,
    nodenumber: *const c_int,
    inarray: *const c_double,
    axis: *const c_int,
) -> c_double {
    unsafe { mesh_update_trampoline(model, nodenumber, inarray, *axis) }
}

/// `freesurfacesolver_` — the free-surface solver entry. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` is unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt` must be a valid read if non-null (it is ignored);
/// `transient` must be a valid read if non-null (it drives the
/// time-discretization branches).
#[export_name = "freesurfacesolver_"]
pub unsafe extern "C" fn freesurface_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let transient_flag: c_int = if transient.is_null() {
        0
    } else {
        (unsafe { *transient != 0 }) as c_int
    };
    unsafe {
        fss_begin(model, solver);
        fss_info_header();
        let _realloc = fss_realloc();
        if fss_free_count() == 0 {
            fss_warn_nofree();
            return;
        }
        fss_read();
        let apply = fss_apply() != 0;
        let need_old = fss_need_old() != 0;
        let nlit = fss_nliter();
        let _rotate = fss_rot_check();
        fss_info_dims();
        fss_allocate();
        if apply {
            fss_bind_residual();
        }
        fss_maybe_reset();
        for iter in 1..=nlit {
            fss_iter_begin(iter, nlit);
            let na = fss_nactive();
            for t in 1..=na {
                let (mut n, mut skip) = (0, 0);
                fss_elem_begin(t, &mut n, &mut skip);
                if skip != 0 {
                    continue;
                }
                if apply {
                    fss_elem_limits(n);
                }
                if fss_use_local() != 0 {
                    fss_elem_maxdisp(n);
                }
                fss_elem_velo(n, transient_flag);
                fss_elem_meshvelo();
                fss_elem_source(n);
                if transient_flag != 0 {
                    fss_elem_old(n);
                }
                fss_local_matrix(n, transient_flag);
                fss_elem_time(transient_flag);
                fss_elem_update(n);
            }
            fss_finish();
            if apply {
                fss_limits_pre();
            }
            fss_info_assdone();
            let mut _relchange = 0.0;
            let _norm = fss_solve(&mut _relchange);
            if apply {
                fss_residual(iter);
            }
            if need_old {
                fss_relax();
            }
            if fss_iter_tail(iter) != 0 {
                break;
            }
        }
    }
}
