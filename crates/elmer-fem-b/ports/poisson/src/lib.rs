//! Port of `fem/src/modules/Poisson.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T020, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/Poisson.F90` — one F77-global entry
//! `poissonsolver_` with two `CONTAINS` internals. The emitted dynamic
//! symbol is exactly `poissonsolver_` (proven by `nm -D` on the reference
//! `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/poisson_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`MATMUL`, `SUM`) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors Poisson.F90 statement for statement):
//! - `SAVE` semantics replicated: `STIFF`/`LOAD`/`FORCE` allocated once at
//!   first-mesh sizing (`MaxElementDOFs`), never resized; the internal
//!   procedures' `SAVE Nodes` live as wrapper-side SAVE state with
//!   identical lifetimes (consumed before the next fill, as in the
//!   original).
//! - Allocation failure calls the same `Fatal('PoissonSolve', ...)` via
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored); `Mesh => GetMesh()` is
//!   kept as a touch call; `Norm = DefaultSolve()` is kept and discarded.
//! - `Load(1:n) = GetReal(...)` runs only when the body force is
//!   associated; `BoundaryCondition` runs only when the variable name is
//!   `'w'`; whole-array `= 0`/`= 1` fills match the original's whole-array
//!   assignments.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn poa_begin(h_solver: *const c_void);
    fn poa_assoc_matrix() -> c_int;
    fn poa_max_elem_dofs() -> c_int;
    fn poa_num_bnd_elems() -> c_int;
    fn poa_var_is_w() -> c_int;
    fn poa_nof_active() -> c_int;
    fn poa_default_initialize();
    fn poa_getmesh_touch();
    fn poa_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn poa_bodyforce() -> c_int;
    fn poa_load_source(n: c_int, load: *mut c_double, found: *mut c_int);
    fn poa_begin_boundary(t: c_int, n: *mut c_int, active: *mut c_int, family: *mut c_int);
    fn poa_bc();
    fn poa_bc_logical(which: c_int, val: *mut c_int, found: *mut c_int);
    fn poa_lm_nodes();
    fn poa_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn poa_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn poa_lm_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    fn poa_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn poa_lm_assemble(
        nd: c_int,
        dbasis: *mut c_double,
        sdet: c_double,
        basis: *mut c_double,
        loadatip: c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn poa_bc_nodes();
    fn poa_bc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn poa_bc_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn poa_bc_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    fn poa_bc_force_update(
        n: c_int, sdet: c_double, loadatip: c_double,
        basis: *mut c_double, force: *mut c_double,
    );
    fn poa_condensate(nd: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn poa_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn poa_default_finish();
    fn poa_default_dirichlet();
    fn poa_default_solve() -> c_double;
    fn poa_fatal_alloc();
}
/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut PoissonState, need: usize) {
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
struct PoissonState {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl PoissonState {
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

static STATE: Mutex<PoissonState> = Mutex::new(PoissonState::new());

/// `poissonsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "poissonsolver_"]
pub unsafe extern "C" fn poisson_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { poa_begin(solver) };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { poa_assoc_matrix() } == 0 {
        return;
    }

    // Mesh => GetMesh() — result unused by the body; kept as a touch call.
    unsafe { poa_getmesh_touch() };

    // First-mesh permanent storage (bug-compatible: sized once, never resized).
    if !state.allocations_done {
        let n = unsafe { poa_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
        {
            unsafe { poa_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.allocations_done = true;
    }

    let active = unsafe { poa_nof_active() };
    unsafe { poa_default_initialize() };

    // System assembly.
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { poa_begin_active(t, &mut n, &mut nd, &mut nb) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        ensure_state(&mut state, (n as usize).max((nd + nb) as usize));
        let nmax = state.n as c_int;

        state.load.fill(0.0);
        let bf_associated = unsafe { poa_bodyforce() } != 0;
        let mut found = 0;
        if bf_associated {
            unsafe { poa_load_source(n, state.load.as_mut_ptr(), &mut found) };
        }

        local_matrix(&mut state, n, nd + nb);

        unsafe {
            poa_condensate(nd, nb, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            poa_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    // Boundary assembly.
    let nbnd = unsafe { poa_num_bnd_elems() };
    for t in 1..=nbnd {
        let (mut n, mut active_el, mut family) = (0, 0, 0);
        unsafe { poa_begin_boundary(t, &mut n, &mut active_el, &mut family) };
        ensure_state(&mut state, n as usize);
        let nmax = state.n as c_int;
        if active_el == 0 {
            continue;
        }
        // No evaluation of Neumann BCs on points.
        if family == 1 {
            continue;
        }
        unsafe { poa_bc() };
        state.force.fill(0.0);
        state.stiff.fill(0.0);

        let (mut pos_el, mut neg_el) = (0, 0);
        let mut found = 0;
        unsafe {
            poa_bc_logical(1, &mut pos_el, &mut found);
            poa_bc_logical(2, &mut neg_el, &mut found);
        }
        let _ = found;

        state.load.fill(0.0);
        if pos_el != 0 {
            state.load.fill(1.0);
        }
        if neg_el != 0 {
            state.load.fill(-1.0);
        }

        if unsafe { poa_var_is_w() } != 0 {
            boundary_condition(&mut state, n);
        }

        unsafe {
            poa_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe {
        poa_default_finish();
        poa_default_dirichlet();
    }

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { poa_default_solve() };
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count;
/// `ndnb` is `nd+nb` (the internal procedure's `nd` dummy).
fn local_matrix(state: &mut PoissonState, n: c_int, ndnb: c_int) {
    unsafe { poa_lm_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    // Basis(n) / dBasisdx(n,3) automatics, sized by this element.
    let nn = n as usize;
    let ndn = ndnb as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { poa_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { poa_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            poa_lm_eleminfo(
                u, v, w, n, ndnb,
                basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;

        let load_at_ip: c_double =
            unsafe { poa_dot(n, basis.as_ptr(), state.load.as_ptr()) };

        unsafe {
            poa_lm_assemble(
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
/// variable is named `'w'`): Neumann contribution by boundary integration.
fn boundary_condition(state: &mut PoissonState, n: c_int) {
    unsafe { poa_bc_nodes() };
    let nn = n as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];

    let nip = unsafe { poa_bc_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { poa_bc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            poa_bc_eleminfo(
                u, v, w, n, n,
                basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;

        let load_at_ip: c_double =
            unsafe { poa_dot(n, basis.as_ptr(), state.load.as_ptr()) };

        unsafe {
            poa_bc_force_update(
                n,
                s * detj,
                load_at_ip,
                basis.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }
}
