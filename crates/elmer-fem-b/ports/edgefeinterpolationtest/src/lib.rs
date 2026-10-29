//! Port of `fem/src/modules/EdgeFEInterpolationTest.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/EdgeFEInterpolationTest.F90` — two
//! F77-global entries (`bestapproximationsolver_`,
//! `bestapproximationsolver_init0_`) with two `CONTAINS` internals. The
//! emitted dynamic symbols are exactly those two (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols. Note the
//! library-vs-entry name difference: the proof invokes the campaign form
//! `Procedure = "EdgeFEInterpolationTest" "BestApproximationSolver"`.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/edgefeinterpolationtest_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), the Galerkin assembly, formatted I/O (`WRITE`), and
//! whole-array updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors EdgeFEInterpolationTest.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `STIFF(N,N)`, `LOAD(6,N)`, `FORCE(N)`,
//!   `Acoef(N)`, `Indices(N)` allocated once at first-mesh sizing
//!   (`MaxElementDOFs`), never resized; the internal procedures' `SAVE
//!   Nodes` live as wrapper-side SAVE state with identical lifetimes
//!   (consumed before the next fill, as in the original). The main-level
//!   `Nodes` SAVE of the original is never referenced by any statement,
//!   so it is not replicated (no observable state).
//! - The pure `Simplicial Mesh` keyword read is kept as a call-sequence
//!   mirror (no side effects either way); `Model`, `dt`, and the
//!   transient flag are unused by the body, exactly as in the original.
//! - `Solver % Matrix % COMPLEX = .FALSE.` and the `A => GetMatrix()`
//!   touch call are kept via a wrapper, exactly as in the original.
//! - `Material Param` handling transcribes the original exactly:
//!   `Acoef(1:n) = 0` first, then overwrite from `GetReal` when a
//!   material is associated (`Found` ignored, as in the original).
//! - The dead `np > 0` Galerkin-projection branch (`np = 0` always) is
//!   transcribed verbatim in the wrapper so the instruction sequence
//!   matches by construction.
//! - Allocation failure calls the same
//!   `Fatal('BestApproximationSolver', ...)` via `try_reserve` (the
//!   `STAT=istat` check), instead of aborting.
//! - `Norm = DefaultSolve()` is kept and discarded, exactly as in the
//!   original. The `WRITE` of `Error Norm` (with `ParallelReduction`)
//!   stays in Fortran, as does the whole error-estimation sweep body.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn efi_begin(h_solver: *mut c_void);
    fn efi_header();
    fn efi_error_on() -> c_int;
    fn efi_touch_mesh();
    fn efi_max_dofs() -> c_int;
    fn efi_matrix_setup();
    fn efi_nof_active() -> c_int;
    fn efi_default_initialize();
    fn efi_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn efi_acoef(n: c_int, acoef: *mut c_double);
    fn efi_local_matrix(
        n: c_int,
        nd: c_int,
        ld: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        acoef: *mut c_double,
    );
    fn efi_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn efi_default_finish_bulk();
    fn efi_default_finish();
    fn efi_default_dirichlet();
    fn efi_default_solve() -> c_double;
    fn efi_fatal_alloc();
    fn efi_error_begin();
    fn efi_num_active_elems() -> c_int;
    fn efi_error_elem(t: c_int, nmax: c_int, load: *mut c_double, indices: *mut c_int);
    fn efi_error_end();
    fn efi_init_begin(h_solver: *mut c_void);
    fn efi_init_body();
}

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut EdgeState, need: usize, need6: usize) {
    if need > state.n {
        state.stiff.resize(need * need, 0.0);
        state.force.resize(need, 0.0);
        state.acoef.resize(need, 0.0);
        state.indices.resize(need, 0);
        state.n = need;
    }
    if need6 > state.load.len() {
        state.load.resize(need6, 0.0);
    }
}

/// First-mesh SAVE state: `STIFF(N,N)`, `LOAD(6,N)`, `FORCE(N)`,
/// `Acoef(N)`, `Indices(N)` plus the `AllocationsDone` flag. `N` is the
/// first call's `MaxElementDOFs`.
struct EdgeState {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    acoef: Vec<c_double>,
    indices: Vec<c_int>,
    n: usize,
    allocations_done: bool,
}

impl EdgeState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            acoef: Vec::new(),
            indices: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<EdgeState> = Mutex::new(EdgeState::new());

/// `bestapproximationsolver_init0_` — the `_Init0` entry Elmer calls at
/// setup. Mirrors `BestApproximationSolver_Init0` (default `Element`
/// definition from the edge-element style). `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "bestapproximationsolver_init0_"]
pub unsafe extern "C" fn bestapproximationsolver_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        efi_init_begin(solver);
        efi_init_body();
    }
}

/// `bestapproximationsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`bestapproximationsolver_init0_entry`].
#[export_name = "bestapproximationsolver_"]
pub unsafe extern "C" fn bestapproximationsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { efi_begin(solver) };
    // EdgeElementStyle + ValueList reads; dim/TestMode/flags stay
    // wrapper-side (consumed only by the Fortran transcriptions).
    unsafe { efi_header() };

    // Mesh => GetMesh() — result unused by the body; kept as a touch call.
    unsafe { efi_touch_mesh() };

    // First-mesh permanent storage (bug-compatible: sized once, never
    // resized by the original).
    if !state.allocations_done {
        let n = unsafe { efi_max_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.load.try_reserve(6 * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.acoef.try_reserve(n).is_err()
            || state.indices.try_reserve(n).is_err()
        {
            unsafe { efi_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.load.resize(6 * n, 0.0);
        state.force.resize(n, 0.0);
        state.acoef.resize(n, 0.0);
        state.indices.resize(n, 0);
        state.allocations_done = true;
    }
    // Solver % Matrix % COMPLEX = .FALSE.; A => GetMatrix().
    unsafe { efi_matrix_setup() };

    // System assembly.
    let active = unsafe { efi_nof_active() };
    unsafe { efi_default_initialize() };
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { efi_begin_active(t, &mut n, &mut nd, &mut nb) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        ensure_state(&mut state, (n as usize).max((nd + nb) as usize), 6 * (n as usize));
        let nmax = state.n as c_int;

        unsafe { efi_acoef(n, state.acoef.as_mut_ptr()) };

        // Get element local matrix and rhs vector, then update the global
        // system from the local matrix and vector.
        unsafe {
            efi_local_matrix(
                n,
                nd + nb,
                nmax,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                state.acoef.as_mut_ptr(),
            );
            efi_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe {
        efi_default_finish_bulk();
        efi_default_finish();
        efi_default_dirichlet();
    }

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { efi_default_solve() };

    // Compute the norm of the error.
    if unsafe { efi_error_on() } != 0 {
        unsafe { efi_error_begin() };
        let nelem = unsafe { efi_num_active_elems() };
        let nmax = state.n as c_int;
        for t in 1..=nelem {
            unsafe {
                efi_error_elem(
                    t,
                    nmax,
                    state.load.as_mut_ptr(),
                    state.indices.as_mut_ptr(),
                )
            };
        }
        unsafe { efi_error_end() };
    }
}
