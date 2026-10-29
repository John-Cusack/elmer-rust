//! Port of `fem/src/modules/ModelMixedPoisson.F90` (Elmer @ 9f6af2f85,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ModelMixedPoisson.F90` — three F77-global
//! entries `mixedpoisson_`, `mixedpoisson_init0_`, `mixedpoisson_post_`
//! with `CONTAINS` internals (`LocalMatrix`, `LocalMatrixBC`, `GetFlux`).
//! The emitted dynamic symbols are exactly those three (proven by `nm -D`
//! on the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/modelmixedpoisson_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`) and whole-array updates stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ModelMixedPoisson.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `STIFF`/`MASS`/`FORCE` allocated once at
//!   first-mesh sizing (`MaxElementDOFs`), never resized; the internal
//!   procedures' `SAVE Nodes` / `SAVE ValueHandle` state and first-call
//!   `InitHandles` flags live as wrapper-side SAVE state with identical
//!   lifetimes (reset at exactly the points the original assigns
//!   `InitHandles = .TRUE.`).
//! - Allocation failure calls the same `Fatal('MixedPoisson', ...)` via
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the main body
//!   except for the transient flag, exactly as in the original;
//!   `DefaultStart()` ordering, `SecondFamily` sampling, bulk/boundary
//!   finish sequencing, and `Norm = DefaultSolve()` (kept and discarded)
//!   match the original statement order.
//! - The `Cond >= 0` gate around `LocalMatrixBC` is evaluated in Rust
//!   from the wrapper-computed `SUM(...)` value, matching the original
//!   branch exactly.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn mmp_begin(h_solver: *const c_void);
    fn mmp_init0_begin(h_model: *const c_void);
    fn mmp_init0_go();
    fn mmp_default_start();
    fn mmp_coord_dim() -> c_int;
    fn mmp_mesh_maxdofs() -> c_int;
    fn mmp_nof_active() -> c_int;
    fn mmp_default_initialize();
    fn mmp_second_family() -> c_int;
    fn mmp_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn mmp_lm_reset();
    #[allow(clippy::too_many_arguments)]
    fn mmp_lm_go(
        ld: c_int,
        stiff: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        n: c_int,
        nd: c_int,
        nb: c_int,
        dim: c_int,
        second_family: c_int,
        need_mass: c_int,
    );
    fn mmp_first_order_time(ld: c_int, stiff: *mut c_double, mass: *mut c_double, force: *mut c_double);
    fn mmp_default_update(ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn mmp_default_finish_bulk();
    fn mmp_nof_bnd() -> c_int;
    fn mmp_begin_bnd(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn mmp_bnd_cond(n: c_int, cond: *mut c_double);
    fn mmp_lbc_reset();
    fn mmp_lbc_go(n: c_int, nd: c_int, second_family: c_int);
    fn mmp_default_finish_bnd();
    fn mmp_default_finish();
    fn mmp_default_dirichlet();
    fn mmp_default_solve() -> c_double;
    fn mmp_fatal_alloc();
    fn mmp_post_go();
}

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut MixedPoissonState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.mass.resize(need * need, 0.0);
    state.force.resize(need, 0.0);
    state.n = need;
}

/// First-mesh SAVE state: `STIFF(N,N)`, `MASS(N,N)`, `FORCE(N)` plus the
/// `AllocationsDone` flag. `N` is the first call's `MaxElementDOFs`.
struct MixedPoissonState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl MixedPoissonState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<MixedPoissonState> = Mutex::new(MixedPoissonState::new());

/// `mixedpoisson_init0_` — solver initialization entry. All dummies except
/// `Model` are unused by the body, exactly as in the original; the whole
/// body runs in one Fortran wrapper so list surgery is bit-identical.
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the loader passes.
#[export_name = "mixedpoisson_init0_"]
pub unsafe extern "C" fn mixedpoisson_init0_entry(
    model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        mmp_init0_begin(model);
        mmp_init0_go();
    }
}

/// `mixedpoisson_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are unused
/// by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read if non-null (Fortran LOGICAL).
#[export_name = "mixedpoisson_"]
pub unsafe extern "C" fn mixedpoisson_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { mmp_begin(solver) };

    unsafe { mmp_default_start() };

    let dim = unsafe { mmp_coord_dim() };

    // Mesh => GetMesh(); N = Mesh % MaxElementDOFs — just big enough.
    // First-mesh permanent storage, done first time only.
    if !state.allocations_done {
        let n = unsafe { mmp_mesh_maxdofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
        {
            unsafe { mmp_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.allocations_done = true;
    }

    // System assembly.
    let second_family = unsafe { mmp_second_family() };
    let need_mass: c_int = if !transient.is_null() && unsafe { *transient } != 0 {
        1
    } else {
        0
    };

    let active = unsafe { mmp_nof_active() };
    unsafe { mmp_default_initialize() };

    unsafe { mmp_lm_reset() };
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { mmp_begin_active(t, &mut n, &mut nd, &mut nb) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        ensure_state(&mut state, (n as usize).max(nd as usize));
        let nmax = state.n as c_int;

        unsafe {
            mmp_lm_go(
                nmax,
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.force.as_mut_ptr(),
                n,
                nd,
                nb,
                dim,
                second_family,
                need_mass,
            )
        };

        if need_mass != 0 {
            unsafe {
                mmp_first_order_time(
                    nmax,
                    state.stiff.as_mut_ptr(),
                    state.mass.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }

        unsafe {
            mmp_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe { mmp_default_finish_bulk() };

    unsafe { mmp_lbc_reset() };
    let nbnd = unsafe { mmp_nof_bnd() };
    for t in 1..=nbnd {
        let (mut n, mut nd, mut active_el) = (0, 0, 0);
        unsafe { mmp_begin_bnd(t, &mut n, &mut nd, &mut active_el) };
        if active_el == 0 {
            continue;
        }
        ensure_state(&mut state, (n as usize).max(nd as usize));
        let mut cond: c_double = 0.0;
        unsafe { mmp_bnd_cond(n, &mut cond) };
        if cond >= 0.0 {
            unsafe { mmp_lbc_go(n, nd, second_family) };
        }
    }

    unsafe {
        mmp_default_finish_bnd();
        mmp_default_finish();
        mmp_default_dirichlet();
    }

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { mmp_default_solve() };
}

/// `mixedpoisson_post_` — postprocessing entry for the main solver. The
/// whole body runs in one Fortran wrapper (including the `GetFlux`
/// internal with its `MATMUL`/`InvertMatrix` reductions), so the
/// postprocessing output is bit-identical by construction.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "mixedpoisson_post_"]
pub unsafe extern "C" fn mixedpoisson_post_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        mmp_begin(solver);
        mmp_post_go();
    }
}
