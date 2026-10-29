//! Port of `fem/src/modules/WPotentialSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/WPotentialSolver.F90` — three F77-global
//! entries `wsolve_`, `wsolve_dummy_`, `wsolve_init0_` (proven by `nm -D`
//! on the reference `.so`); the `CONTAINS` internals and the commented-out
//! `SaveWPotSolution`/`BoundaryCondition` bodies emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/wpotentialsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`/`MATMUL`), the `Info` call, whole-array updates, string
//! comparisons, and the solver-list surgery of `Wsolve_Init0` stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors WPotentialSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `STIFF`/`LOAD`/`FORCE`/`Tcoef`/`RotM`
//!   allocated once at first-mesh sizing (`MaxElementDOFs`), never
//!   resized; the `NoRotM` host SAVE (initialized `.FALSE.` once,
//!   rewritten only in the stranded branch, read by the foil path) is
//!   carried in the Rust state with the same lifetime; `Init0`'s
//!   `visited` counter and the internal procedures' `SAVE Nodes` /
//!   `SAVE First` / variable-handle state live as wrapper-side SAVE
//!   state with identical lifetimes.
//! - Allocation failure calls the same `Fatal('Wsolve', ...)` via
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the `Wsolve` body
//!   except for `Model % NumberOfBodies`, exactly as in the original;
//!   `Mesh => GetMesh()` is kept as a touch call (it also publishes the
//!   mesh handle the save kernel's `VariableGet` needs); `Norm =
//!   DefaultSolve()` is kept and discarded.
//! - `Load(1:n) = GetReal(...)` runs only when the body force is
//!   associated; the boundary loop keeps the electrode queries and the
//!   `LOAD` fills even though the commented-out `BoundaryCondition` call
//!   means the assembled zeros are the only contribution, exactly as in
//!   the original.
//! - `Tcoef = GetElectricConductivityTensor(...)` is a whole-array
//!   assignment of a `(3,3,n)` result into the `(3,3,N)` SAVE array, so
//!   valid inputs have `N == n` per element, as in the original (the
//!   `SIZE(Tcoef,3) /= n` reallocation check is a same-size no-op there
//!   and is skipped: nothing observable can depend on it).
//! - The commented-out `SaveWPotSolution` call and `BoundaryCondition`
//!   internal procedure are dead in the original and stay absent here.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn wpo_init0_begin(h_model: *const c_void, h_solver: *const c_void);
    fn wpo_init0_enter() -> c_int;
    fn wpo_init0_extend();
    fn wpo_begin(h_model: *const c_void, h_solver: *const c_void);
    fn wpo_assoc_matrix() -> c_int;
    fn wpo_add_components();
    fn wpo_max_elem_dofs() -> c_int;
    fn wpo_getmesh_touch();
    fn wpo_numberofbodies() -> c_int;
    fn wpo_nof_active() -> c_int;
    fn wpo_default_initialize();
    fn wpo_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn wpo_body_id() -> c_int;
    fn wpo_bodyforce() -> c_int;
    fn wpo_load_source(n: c_int, load: *mut c_double, found: *mut c_int);
    /// 0 = no component, 1 = stranded, 2 = massive, 3 = foil winding,
    /// -1 = unknown (the original Fatals on these).
    fn wpo_comp_class() -> c_int;
    fn wpo_comp_norotm() -> c_int;
    fn wpo_fatal_coil();
    fn wpo_rotm_fill(class: c_int, ld: c_int, rotm: *mut c_double, n: c_int, norotm: c_int);
    fn wpo_tcoef_fill(class: c_int, ld: c_int, tcoef: *mut c_double, n: c_int);
    #[allow(clippy::too_many_arguments)]
    fn wpo_localmatrix(
        ld: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        load: *mut c_double,
        tcoef: *mut c_double,
        rotm: *mut c_double,
        n: c_int,
        nd: c_int,
        class: c_int,
        norotm: c_int,
    );
    fn wpo_condensate(nd: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn wpo_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn wpo_num_bnd_elems() -> c_int;
    fn wpo_begin_boundary(t: c_int, n: *mut c_int, active: *mut c_int, family: *mut c_int);
    fn wpo_bc();
    fn wpo_bc_electrode(which: c_int, val: *mut c_int, found: *mut c_int);
    fn wpo_default_finish();
    fn wpo_default_dirichlet();
    fn wpo_default_solve() -> c_double;
    fn wpo_wnorms(nofbodies: c_int, wnorms: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn wpo_save_element(
        ld: c_int,
        tcoef: *mut c_double,
        rotm: *mut c_double,
        n: c_int,
        wnorm: c_double,
        class: c_int,
        norotm: c_int,
    );
    fn wpo_fatal_alloc();
}

/// Grow the assembly SAVE arrays when an element needs more than the
/// first-mesh sizing. The original never grows (its sizing comment
/// assumes the first mesh is big enough); growth only triggers where
/// the original would overflow silently, and the leading sections used
/// by the wrappers are unaffected, so valid inputs never observe it.
/// `Tcoef`/`RotM` stay fixed at the first-mesh size, exactly like the
/// original's `(3,3,N)` SAVE arrays.
fn ensure_state(state: &mut WpoState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.load.resize(need, 0.0);
    state.force.resize(need, 0.0);
    state.n = need;
}

/// First-mesh SAVE state: `STIFF(N,N)`, `LOAD(N)`, `FORCE(N)`,
/// `Tcoef(3,3,N)`, `RotM(3,3,N)`, the `AllocationsDone` flag, and the
/// `NoRotM` host SAVE (initialized `.FALSE.` once). `N` is the first
/// call's `MaxElementDOFs`.
struct WpoState {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    tcoef: Vec<c_double>,
    rotm: Vec<c_double>,
    n: usize,
    allocations_done: bool,
    norotm: bool,
}

impl WpoState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            tcoef: Vec::new(),
            rotm: Vec::new(),
            n: 0,
            allocations_done: false,
            norotm: false,
        }
    }
}

static STATE: Mutex<WpoState> = Mutex::new(WpoState::new());

/// `wsolve_init0_` — solver-list surgery creating the on-the-fly DG
/// solver. The `visited` counter increment and the Discontinuous
/// Galerkin early return live in the wrapper, exactly like the
/// original's `SAVE visited` and first statement.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "wsolve_init0_"]
pub unsafe extern "C" fn wsolve_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { wpo_init0_begin(model, solver) };
    if unsafe { wpo_init0_enter() } != 0 {
        return;
    }
    unsafe { wpo_init0_extend() };
}

/// `wsolve_dummy_` — the never-executed DG placeholder solver. Empty,
/// exactly as in the original.
///
/// # Safety
///
/// Takes the standard entry addresses and ignores them.
#[export_name = "wsolve_dummy_"]
pub unsafe extern "C" fn wsolve_dummy_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}

/// Coil class for the current element, mirroring the original's
/// `SELECT CASE (CoilType)` including the `Fatal` on unknown types.
unsafe fn coil_class() -> c_int {
    let class = unsafe { wpo_comp_class() };
    if class < 0 {
        unsafe { wpo_fatal_coil() };
        std::process::abort();
    }
    class
}

/// `wsolve_` — the solver entry. Argument order and by-reference travel
/// match the Fortran dummy list `(Model, Solver, dt,
/// TransientSimulation)`; `dt` and the transient flag are unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "wsolve_"]
pub unsafe extern "C" fn wsolve_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { wpo_begin(model, solver) };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { wpo_assoc_matrix() } == 0 {
        return;
    }
    unsafe { wpo_add_components() };

    // First-mesh permanent storage (sized once, never resized on valid
    // inputs; bug-compatible with the original's AllocationsDone gate).
    if !state.allocations_done {
        let n = unsafe { wpo_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.tcoef.try_reserve(9 * n).is_err()
            || state.rotm.try_reserve(9 * n).is_err()
        {
            unsafe { wpo_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.tcoef.resize(9 * n, 0.0);
        state.rotm.resize(9 * n, 0.0);
        state.allocations_done = true;
    }

    // Mesh => GetMesh() — kept as a touch call; it also publishes the
    // mesh handle the save kernel's VariableGet needs.
    unsafe { wpo_getmesh_touch() };

    let active = unsafe { wpo_nof_active() };
    unsafe { wpo_default_initialize() };

    // System assembly.
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { wpo_begin_active(t, &mut n, &mut nd, &mut nb) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        ensure_state(&mut state, (n as usize).max((nd + nb) as usize));
        let nmax = state.n as c_int;

        state.load.fill(0.0);
        let bf_associated = unsafe { wpo_bodyforce() } != 0;
        let mut found = 0;
        if bf_associated {
            unsafe { wpo_load_source(n, state.load.as_mut_ptr(), &mut found) };
        }
        let _ = found;

        let class = unsafe { coil_class() };
        if class == 1 {
            state.norotm = unsafe { wpo_comp_norotm() } != 0;
        }
        let norotm = c_int::from(state.norotm);
        unsafe {
            wpo_rotm_fill(class, nmax, state.rotm.as_mut_ptr(), n, norotm);
            wpo_tcoef_fill(class, nmax, state.tcoef.as_mut_ptr(), n);
        }

        state.stiff.fill(0.0);
        state.force.fill(0.0);
        unsafe {
            wpo_localmatrix(
                nmax,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                state.load.as_mut_ptr(),
                state.tcoef.as_mut_ptr(),
                state.rotm.as_mut_ptr(),
                n,
                nd + nb,
                class,
                norotm,
            );
            wpo_condensate(nd, nb, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            wpo_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    // Boundary assembly (the BoundaryCondition call is commented out in
    // the original, so the zeroed STIFF/FORCE are the only contribution;
    // the electrode queries and LOAD fills are kept for fidelity).
    let nbnd = unsafe { wpo_num_bnd_elems() };
    for t in 1..=nbnd {
        let (mut n, mut active_el, mut family) = (0, 0, 0);
        unsafe { wpo_begin_boundary(t, &mut n, &mut active_el, &mut family) };
        ensure_state(&mut state, n as usize);
        let nmax = state.n as c_int;
        if active_el == 0 {
            continue;
        }
        // No evaluation of Neumann BC's on points.
        if family == 1 {
            continue;
        }
        unsafe { wpo_bc() };
        state.force.fill(0.0);
        state.stiff.fill(0.0);

        let (mut pos_el, mut neg_el) = (0, 0);
        let mut found = 0;
        unsafe {
            wpo_bc_electrode(1, &mut pos_el, &mut found);
            wpo_bc_electrode(2, &mut neg_el, &mut found);
        }
        let _ = found;

        state.load.fill(0.0);
        if pos_el != 0 {
            state.load.fill(1.0);
        }
        if neg_el != 0 {
            state.load.fill(-1.0);
        }

        unsafe {
            wpo_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe {
        wpo_default_finish();
        wpo_default_dirichlet();
    }

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { wpo_default_solve() };

    // Wnorms = GetWNormsForBodies(Model % NumberOfBodies)
    let nbf = unsafe { wpo_numberofbodies() };
    let mut wnorms = vec![0.0 as c_double; nbf.max(0) as usize];
    unsafe { wpo_wnorms(nbf, wnorms.as_mut_ptr()) };

    // Per-element DG save pass.
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { wpo_begin_active(t, &mut n, &mut nd, &mut nb) };
        let _ = (nd, nb);
        ensure_state(&mut state, n as usize);
        let nmax = state.n as c_int;

        let class = unsafe { coil_class() };
        if class == 1 {
            state.norotm = unsafe { wpo_comp_norotm() } != 0;
        }
        let norotm = c_int::from(state.norotm);
        unsafe {
            wpo_rotm_fill(class, nmax, state.rotm.as_mut_ptr(), n, norotm);
            wpo_tcoef_fill(class, nmax, state.tcoef.as_mut_ptr(), n);
        }

        let body = unsafe { wpo_body_id() } as usize;
        let wnorm = wnorms[body.saturating_sub(1)];
        unsafe {
            wpo_save_element(
                nmax,
                state.tcoef.as_mut_ptr(),
                state.rotm.as_mut_ptr(),
                n,
                wnorm,
                class,
                norotm,
            );
        }
    }
}
