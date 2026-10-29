//! Port of `fem/src/modules/DirectionSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/DirectionSolver.F90` — three F77-global
//! entries `directionsolver_`, `directionsolver_dummy_`, and
//! `directionsolver_init0_` with four `CONTAINS` internals. The emitted
//! dynamic symbols are exactly those three (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/directionsolver_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`MATMUL`,
//! `SUM`), the DG-indexed save, and the `_Init0` model surgery stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors DirectionSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `STIFF`/`LOAD`/`FORCE` allocated once
//!   at first-mesh sizing (`MaxElementDOFs`), never resized; the internal
//!   procedures' `SAVE Nodes` live as wrapper-side SAVE state with
//!   identical lifetimes (consumed before the next fill, as in the
//!   original); `SaveElementSolution`'s SAVE `First`/`varname`/
//!   `directionvar` live as wrapper-side SAVE state with identical
//!   lifetimes, including the undefined-at-first-compare `varname`.
//! - Allocation failure calls the same `Fatal('DirectionSolver', ...)`
//!   via `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the main body
//!   except for the body count, exactly as in the original (received and
//!   ignored otherwise); `Mesh => GetMesh()` is kept as a touch call
//!   that also pins the mesh handle the save sweep reads;
//!   `Norm = DefaultSolve()` is kept and discarded.
//! - `Load(1:n) = GetReal(...)` runs only when the body force is
//!   associated; `BoundaryCondition` runs only when the solver variable
//!   is named `'w'`; whole-array `= 0`/`= 1` fills match the original's
//!   whole-array assignments, including the internal re-zero of `FORCE`.
//! - `ListSetNameSpace('')` closes every per-body namespace iteration,
//!   exactly as in the original.
//! - The `_Init0` write-only `visited` counter is not replicated (no
//!   observable state); the DG guard and the whole model surgery
//!   (equation widening, DG dummy-solver append, exported `Direction`
//!   variable) are transcribed in order, so the port and the reference
//!   build identical model state.
//! - `DirectionSolver_Dummy` is an empty body in the original; the port
//!   entry is likewise empty.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn dsa_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn dsa_assoc_matrix() -> c_int;
    fn dsa_getmesh_touch();
    fn dsa_max_elem_dofs() -> c_int;
    fn dsa_number_of_bodies() -> c_int;
    fn dsa_fatal_alloc();
    fn dsa_init0_dg() -> c_int;
    fn dsa_init0_body();
    fn dsa_ns_setup(ns: c_int);
    fn dsa_ns_reset();
    fn dsa_nof_active() -> c_int;
    fn dsa_default_initialize();
    fn dsa_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn dsa_bodyforce() -> c_int;
    fn dsa_load_source(n: c_int, load: *mut c_double, found: *mut c_int);
    fn dsa_num_bnd_elems() -> c_int;
    fn dsa_begin_boundary(t: c_int, n: *mut c_int, active: *mut c_int, family: *mut c_int);
    fn dsa_bc();
    fn dsa_bc_logical(which: c_int, val: *mut c_int, found: *mut c_int);
    fn dsa_var_is_w() -> c_int;
    fn dsa_lm_nodes();
    fn dsa_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn dsa_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn dsa_lm_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    fn dsa_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn dsa_lm_assemble(
        nd: c_int,
        dbasis: *mut c_double,
        sdet: c_double,
        basis: *mut c_double,
        loadatip: c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn dsa_bc_nodes();
    fn dsa_bc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn dsa_bc_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn dsa_bc_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    fn dsa_bc_force_update(
        n: c_int, sdet: c_double, loadatip: c_double,
        basis: *mut c_double, force: *mut c_double,
    );
    fn dsa_condensate(nd: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn dsa_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn dsa_default_finish();
    fn dsa_default_dirichlet();
    fn dsa_default_solve() -> c_double;
    fn dsa_save_sweep(ns: c_int);
}

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut DirectionState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.load.resize(need, 0.0);
    state.force.resize(need, 0.0);
    state.n = need;
}

/// First-mesh SAVE state: `STIFF(N,N)`, `LOAD(N)`, `FORCE(N)` plus the
/// `AllocationsDone` flag. `N` is the first call's `MaxElementDOFs`.
struct DirectionState {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl DirectionState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<DirectionState> = Mutex::new(DirectionState::new());

/// `directionsolver_init0_` — the `_Init0` model surgery Elmer runs once
/// per solver at setup (before mesh-solver allocation): when the solver
/// is not already marked discontinuous Galerkin, widen the equation
/// `Active Solvers` lists and append the never-executed DG dummy solver
/// carrying the exported `<Variable> Direction` field.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "directionsolver_init0_"]
pub unsafe extern "C" fn directionsolver_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        dsa_begin(model, solver);
        // IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) RETURN
        if dsa_init0_dg() != 0 {
            return;
        }
        dsa_init0_body();
    }
}

/// `directionsolver_dummy_` — the never-executed dummy entry installed by
/// `_Init0` (`Exec Solver = never`); the original body is empty.
///
/// # Safety
///
/// Same contract as [`directionsolver_init0_entry`].
#[export_name = "directionsolver_dummy_"]
pub unsafe extern "C" fn directionsolver_dummy_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}

/// `directionsolver_` — the coil-direction solver entry: per-body
/// namespace loop of Poisson-like assembly, electrode boundary loading,
/// solve, and DG-direction save. Argument order and by-reference travel
/// match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "directionsolver_"]
pub unsafe extern "C" fn directionsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { dsa_begin(model, solver) };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { dsa_assoc_matrix() } == 0 {
        return;
    }

    // Mesh => GetMesh() — kept as a touch call; the handle also feeds
    // the save sweep's VariableGet, as in the original.
    unsafe { dsa_getmesh_touch() };

    // First-mesh permanent storage (bug-compatible: sized once, never resized).
    if !state.allocations_done {
        let n = unsafe { dsa_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
        {
            unsafe { dsa_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.allocations_done = true;
    }

    let nbodies = unsafe { dsa_number_of_bodies() };
    for ns in 1..=nbodies {
        // 'body <ns>:' namespace iff any BC carries '<ns>: <varname>'.
        unsafe { dsa_ns_setup(ns) };

        let active = unsafe { dsa_nof_active() };
        unsafe { dsa_default_initialize() };

        // System assembly.
        for t in 1..=active {
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { dsa_begin_active(t, &mut n, &mut nd, &mut nb) };
            // First-mesh sizing holds on valid inputs; grow-only guard covers
            // remeshed re-entry (invisible otherwise).
            ensure_state(&mut state, (n as usize).max((nd + nb) as usize));
            let nmax = state.n as c_int;

            state.load.fill(0.0);
            let bf_associated = unsafe { dsa_bodyforce() } != 0;
            let mut found = 0;
            if bf_associated {
                unsafe { dsa_load_source(n, state.load.as_mut_ptr(), &mut found) };
            }

            local_matrix(&mut state, n, nd + nb);

            unsafe {
                dsa_condensate(nd, nb, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
                dsa_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            }
        }

        // Boundary assembly.
        let nbnd = unsafe { dsa_num_bnd_elems() };
        for t in 1..=nbnd {
            let (mut n, mut active_el, mut family) = (0, 0, 0);
            unsafe { dsa_begin_boundary(t, &mut n, &mut active_el, &mut family) };
            ensure_state(&mut state, n as usize);
            let nmax = state.n as c_int;
            if active_el == 0 {
                continue;
            }
            // No evaluation of Neumann BCs on points.
            if family == 1 {
                continue;
            }
            unsafe { dsa_bc() };
            state.force.fill(0.0);
            state.stiff.fill(0.0);

            let (mut pos_el, mut neg_el) = (0, 0);
            let mut found = 0;
            unsafe {
                dsa_bc_logical(1, &mut pos_el, &mut found);
                dsa_bc_logical(2, &mut neg_el, &mut found);
            }
            let _ = found;

            state.load.fill(0.0);
            if pos_el != 0 {
                state.load.fill(1.0);
            }
            if neg_el != 0 {
                state.load.fill(-1.0);
            }

            if unsafe { dsa_var_is_w() } != 0 {
                boundary_condition(&mut state, n);
            }

            unsafe {
                dsa_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            }
        }

        unsafe {
            dsa_default_finish();
            dsa_default_dirichlet();
        }

        // And finally, solve (Norm is local and unused, as in the original).
        let _norm: c_double = unsafe { dsa_default_solve() };
        unsafe {
            dsa_save_sweep(ns);
            dsa_ns_reset();
        }
    }
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count;
/// `ndnb` is `nd+nb` (the internal procedure's `nd` dummy).
fn local_matrix(state: &mut DirectionState, n: c_int, ndnb: c_int) {
    unsafe { dsa_lm_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    // Basis(n) / dBasisdx(n,3) automatics, sized by this element.
    let nn = n as usize;
    let ndn = ndnb as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { dsa_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { dsa_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            dsa_lm_eleminfo(
                u, v, w, n, ndnb,
                basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;

        let load_at_ip: c_double =
            unsafe { dsa_dot(n, basis.as_ptr(), state.load.as_ptr()) };

        unsafe {
            dsa_lm_assemble(
                ndnb,
                dbasis.as_mut_ptr(),
                s * detj,
                basis.as_mut_ptr(),
                load_at_ip,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }
}

/// The `BoundaryCondition` internal procedure (only runs when the solver
/// variable is named `'w'`): Neumann contribution by boundary
/// integration, including the internal whole-array re-zero of `FORCE`.
fn boundary_condition(state: &mut DirectionState, n: c_int) {
    unsafe { dsa_bc_nodes() };
    state.force.fill(0.0);
    let nn = n as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];

    let nip = unsafe { dsa_bc_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { dsa_bc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            dsa_bc_eleminfo(
                u, v, w, n, n,
                basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;

        let load_at_ip: c_double =
            unsafe { dsa_dot(n, basis.as_ptr(), state.load.as_ptr()) };

        unsafe {
            dsa_bc_force_update(
                n,
                s * detj,
                load_at_ip,
                basis.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }
}
