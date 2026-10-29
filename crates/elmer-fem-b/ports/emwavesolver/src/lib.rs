//! Port of `fem/src/modules/EMWaveSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/EMWaveSolver.F90` — five F77-global
//! entries (`emwavesolver_`, `emwavesolver_init0_`,
//! `emwavecalcfields_`, `emwavecalcfields_init_`,
//! `emwavecalcfields_init0_`) with `CONTAINS` internals. The emitted
//! dynamic symbols are exactly those five (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/emwavesolver_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the edge-element Galerkin assembly with its
//! `SUM`/`MATMUL` reductions, the whole-array fills, the strided eigen
//! copies, the LU solves) stays in Fortran wrappers transcribed statement
//! for statement, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors EMWaveSolver.F90 statement for statement):
//! - `SAVE` semantics replicated: the entry work arrays
//!   (`STIFF`/`DAMP`/`MASS`/`FORCE`) are allocated once at first-mesh
//!   sizing (`MaxElementDOFs`) behind a `static` guard and never resized;
//!   the internal procedures' `SAVE` state (`Nodes`, `ValueHandle`s, the
//!   `SAVE`-by-initialization `ConstantBulkInUse`, `Visited`, and the
//!   `LocalAssembly` work arrays) lives wrapper-side with identical
//!   lifetimes (consumed before the next fill, as in the original).
//! - The main-entry nonlinear loop (`NoIterationsMax`, `Edge Basis`
//!   refactorize add/remove, `DefaultConverged` exit) is driven from Rust;
//!   the bulk/boundary assembly sweeps (including the `InitHandles`
//!   first-element detection and the `ConstantBulkMatrix` short-circuit)
//!   are single wrapper calls, so the library call sequence matches.
//! - `EMWaveCalcFields` keeps the author's structure: the per-eigen
//!   iteration (`GOTO 10` loop) and the active-element sweep are driven
//!   from Rust; `LocalAssembly`, `GlobalSol`, `LocalSol`,
//!   `AddLocalFaceTerms`, `LocalJumps` are wrapper-internal procedures
//!   with the same host-variable sharing as the original `CONTAINS`
//!   bodies (including the `Fsave` RHS-pointer juggling and the `dofs`
//!   column counter).
//! - Allocation failure of the entry SAVE arrays calls the same
//!   `Fatal('EMWaveSolver', ...)` via `try_reserve` (the `STAT=istat`
//!   check), instead of aborting; the `CalcFields` per-invocation
//!   workspaces use the same treatment with the sibling message.
//! - `Model`, `dt`, `TransientSimulation` are unused by the main body
//!   (received and ignored), exactly as in the original; the `Perm`
//!   touch and the `Mesh => GetMesh()` call are kept.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn emw_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn emw_info_start();
    fn emw_prologue();
    fn emw_maxdofs() -> c_int;
    fn emw_noitmax() -> c_int;
    fn emw_edgebasis() -> c_int;
    fn emw_default_start();
    fn emw_bulk(
        stiff: *mut c_double,
        damp: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        n: c_int,
    );
    fn emw_bnd(
        stiff: *mut c_double,
        damp: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        n: c_int,
    );
    fn emw_solve() -> c_double;
    fn emw_converged() -> c_int;
    fn emw_norefactor();
    fn emw_rmrefactor();
    fn emw_finish();
    fn emw_done();
    fn emw_fatal_alloc();
    fn emw_i0_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn emw_i0_body();
    fn emw_cfi0_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn emw_cfi0_body();
    fn emw_cfi_begin(h_solver: *mut c_void);
    fn emw_cfi_body();
    fn emw_cf_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn emw_cf_setup();
    fn emw_cf_counts(
        nofeigen: *mut c_int,
        dofs: *mut c_int,
        nmax: *mut c_int,
        nrhs: *mut c_int,
        ncols: *mut c_int,
    );
    fn emw_cf_nactive() -> c_int;
    fn emw_cf_eigenhead(ieigen: c_int);
    fn emw_cf_elem(
        t: c_int,
        ieigen: c_int,
        mass: *mut c_double,
        force: *mut c_double,
        gf: *mut c_double,
        pivot: *mut c_int,
        n0: c_int,
        d0: c_int,
        nr: c_int,
        nc: c_int,
    );
    fn emw_cf_post(
        ieigen: c_int,
        mass: *mut c_double,
        force: *mut c_double,
        gf: *mut c_double,
        n0: c_int,
        d0: c_int,
        nr: c_int,
        nc: c_int,
    );
    fn emw_cf_done();
    fn emw_cf_fatal_alloc();
}

/// First-mesh SAVE state: `STIFF(N,N)`, `DAMP(N,N)`, `MASS(N,N)`,
/// `FORCE(N)` plus the `AllocationsDone` flag. `N` is the first call's
/// `MaxElementDOFs`.
struct EmwState {
    n: usize,
    stiff: Vec<c_double>,
    damp: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    done: bool,
}

impl EmwState {
    const fn new() -> Self {
        Self {
            n: 0,
            stiff: Vec::new(),
            damp: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            done: false,
        }
    }
}

static STATE: Mutex<EmwState> = Mutex::new(EmwState::new());

/// `emwavesolver_init0_` — the `_Init0` entry Elmer calls once at solver
/// setup. Mirrors `EMWaveSolver_Init0` (edge-element `Element` string,
/// `Hcurl Basis` / `Variable` / `Time derivative order` defaults, vacuum
/// constants check).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "emwavesolver_init0_"]
pub unsafe extern "C" fn emwavesolver_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emw_i0_begin(model, solver);
        emw_i0_body();
    }
}

/// `emwavesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// Same contract as [`emwavesolver_init0_entry`].
#[export_name = "emwavesolver_"]
pub unsafe extern "C" fn emwavesolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emw_begin(model, solver);
        emw_info_start();
        emw_prologue();
        let n = emw_maxdofs();
        let mut state = STATE.lock().unwrap();
        if !state.done {
            let nn = n as usize;
            if state.stiff.try_reserve(nn * nn).is_err()
                || state.damp.try_reserve(nn * nn).is_err()
                || state.mass.try_reserve(nn * nn).is_err()
                || state.force.try_reserve(nn).is_err()
            {
                emw_fatal_alloc();
            }
            state.stiff.resize(nn * nn, 0.0);
            state.damp.resize(nn * nn, 0.0);
            state.mass.resize(nn * nn, 0.0);
            state.force.resize(nn, 0.0);
            state.n = nn;
            state.done = true;
        }
        let n = state.n as c_int;
        let (stiff, damp, mass, force) = (
            state.stiff.as_mut_ptr(),
            state.damp.as_mut_ptr(),
            state.mass.as_mut_ptr(),
            state.force.as_mut_ptr(),
        );
        let noit = emw_noitmax();
        let edgebasis = emw_edgebasis() != 0;
        emw_default_start();
        for _ in 0..noit {
            emw_bulk(stiff, damp, mass, force, n);
            emw_bnd(stiff, damp, mass, force, n);
            let _norm: c_double = emw_solve();
            if emw_converged() != 0 {
                break;
            }
            if edgebasis {
                emw_norefactor();
            }
        }
        if edgebasis {
            emw_rmrefactor();
        }
        emw_finish();
        emw_done();
    }
}

/// `emwavecalcfields_init0_` — the `_Init0` entry Elmer calls once at
/// solver setup. Mirrors `EMWaveCalcFields_Init0` (primary-solver lookup,
/// eigen flags, DG solver-instance injection for elemental fields).
///
/// # Safety
///
/// Same contract as [`emwavesolver_init0_entry`].
#[export_name = "emwavecalcfields_init0_"]
pub unsafe extern "C" fn emwavecalcfields_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emw_cfi0_begin(model, solver);
        emw_cfi0_body();
    }
}

/// `emwavecalcfields_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `EMWaveCalcFields_Init` (dummy `Variable`,
/// no-refactorize default, `Exported Variable` requests).
///
/// # Safety
///
/// Same contract as [`emwavesolver_init0_entry`].
#[export_name = "emwavecalcfields_init_"]
pub unsafe extern "C" fn emwavecalcfields_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _ = model;
    unsafe {
        emw_cfi_begin(solver);
        emw_cfi_body();
    }
}

/// `emwavecalcfields_` — the postprocessing entry. Mirrors
/// `EMWaveCalcFields` (nodal/elemental field assembly, per-eigen
/// iteration, face-term averaging). Per-invocation workspaces
/// (`MASS`, `FORCE`, `Pivot`, `GForce`) are sized from the setup query
/// and dropped at return, matching the original's per-call
/// allocate/deallocate pair.
///
/// # Safety
///
/// Same contract as [`emwavesolver_init0_entry`].
#[export_name = "emwavecalcfields_"]
pub unsafe extern "C" fn emwavecalcfields_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emw_cf_begin(model, solver);
        emw_cf_setup();
        let (mut nofeigen, mut dofs, mut nmax, mut nrhs, mut ncols): (
            c_int,
            c_int,
            c_int,
            c_int,
            c_int,
        ) = (0, 0, 0, 0, 0);
        emw_cf_counts(&mut nofeigen, &mut dofs, &mut nmax, &mut nrhs, &mut ncols);
        // Zero-field configurations leave a dimension at 0 (the original
        // leaves the corresponding array unallocated); keep a 1-wide
        // dummy so the C pointers stay valid on paths the original never
        // touches. Defined paths (dofs/ncols > 0) are unaffected.
        let n0 = (nmax as usize).max(1);
        let d0 = (dofs as usize).max(1);
        let nr = (nrhs as usize).max(1);
        let nc = (ncols as usize).max(1);
        let mut mass = Vec::<c_double>::new();
        let mut force = Vec::<c_double>::new();
        let mut gf = Vec::<c_double>::new();
        let mut pivot = Vec::<c_int>::new();
        if mass.try_reserve(n0 * n0).is_err()
            || force.try_reserve(n0 * d0).is_err()
            || gf.try_reserve(nr * nc).is_err()
            || pivot.try_reserve(n0).is_err()
        {
            emw_cf_fatal_alloc();
        }
        mass.resize(n0 * n0, 0.0);
        force.resize(n0 * d0, 0.0);
        gf.resize(nr * nc, 0.0);
        pivot.resize(n0, 0);
        for ie in 1..=nofeigen {
            emw_cf_eigenhead(ie);
            let na = emw_cf_nactive();
            for t in 1..=na {
                emw_cf_elem(
                    t,
                    ie,
                    mass.as_mut_ptr(),
                    force.as_mut_ptr(),
                    gf.as_mut_ptr(),
                    pivot.as_mut_ptr(),
                    n0 as c_int,
                    d0 as c_int,
                    nr as c_int,
                    nc as c_int,
                );
            }
            emw_cf_post(
                ie,
                mass.as_mut_ptr(),
                force.as_mut_ptr(),
                gf.as_mut_ptr(),
                n0 as c_int,
                d0 as c_int,
                nr as c_int,
                nc as c_int,
            );
        }
        emw_cf_done();
    }
}
