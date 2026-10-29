//! Port of `fem/src/modules/PoissonDG.F90` (Elmer @ the campaign fork)
//! for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/PoissonDG.F90` — one F77-global entry
//! `poissonsolver_` with four `CONTAINS` internals. The emitted dynamic
//! symbol is exactly `poissonsolver_` (proven by `nm -D` on the reference
//! `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/poissondg_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`SUM`), `MATMUL` updates, and whole-array fills
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors PoissonDG.F90 statement for statement):
//! - `SAVE` semantics replicated: `STIFF`/`LOAD`/`FORCE`/`EpsilonBoundary`
//!   allocated once at first-mesh sizing (`N = 2*MaxElementDOFs`, with
//!   `STIFF` `(2N,2N)`), never resized; the internal procedures' `SAVE
//!   Nodes` (+ IP state, + per-variant `LeftOut`) live as wrapper-side
//!   SAVE state with identical lifetimes (consumed before the next fill,
//!   as in the original).
//! - Allocation failure calls the same `Fatal('PoissonDG', ...)` via
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored); `Mesh => GetMesh()` is
//!   kept as a touch call; `Norm = DefaultSolve()` is kept and discarded.
//! - `LOAD(1:n) = GetReal(...)` runs unconditionally (also when the body
//!   force is missing — `GetReal` on a null list is the original's own
//!   statement, kept in the wrapper); `FORCE` is zeroed before the edge
//!   loop exactly as in the original.
//! - `LocalJumps` and `LocalJumpsDiscontBC` are textually identical in the
//!   original (only the names differ — verified by diff); each gets its
//!   own wrapper set with its own SAVE state, and the entry dispatches on
//!   the `'Discont BC'` logical exactly as the original does.
//! - `LocalMatrixBoundary` is dead code in the original (never called from
//!   the entry) and has no wrappers; `EpsilonBoundary` is allocated (part
//!   of the `STAT` check) but never read on any live path, so it is kept
//!   as a zeroed SAVE mirror only. `LocalSolution` (a local `POINTER`
//!   allocated but never used, not even `SAVE`d) has no observable effect
//!   and is omitted.
//! - `gamma2` is read and discarded exactly as in the original (it only
//!   feeds the dead `LocalMatrixBoundary`); `beta` defaults to `gamma`
//!   when the list has no `'beta'` entry.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pda_begin(h_solver: *const c_void);
    fn pda_touch_mesh();
    fn pda_first_n() -> c_int;
    fn pda_edge_counts(dim: *mut c_int, np: *mut c_int);
    fn pda_nof_active() -> c_int;
    fn pda_default_initialize();
    fn pda_gamma(which: c_int) -> c_double;
    fn pda_beta(gamma: c_double) -> c_double;
    fn pda_begin_active(t: c_int, n: *mut c_int);
    fn pda_bodyforce();
    fn pda_load_source(n: c_int, load: *mut c_double, found: *mut c_int);
    fn pda_lm_nodes();
    fn pda_lm_ip_begin() -> c_int;
    fn pda_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    fn pda_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn pda_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn pda_lm_assemble(
        n: c_int,
        dbasis: *mut c_double,
        sdet: c_double,
        basis: *mut c_double,
        loadatip: c_double,
        ld: c_int,
        stiff: *mut c_double,
        lf: c_int,
        force: *mut c_double,
    );
    fn pda_begin_edge(
        t: c_int,
        n: *mut c_int,
        nl: *mut c_int,
        nr: *mut c_int,
        active: *mut c_int,
        hasright: *mut c_int,
    );
    fn pda_edge_discont() -> c_int;
    fn pda_lj_setup(n: c_int, nl: c_int, nr: c_int, he: *mut c_double);
    fn pda_lj_ip_begin() -> c_int;
    fn pda_lj_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn pda_lj_point(
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
        n: c_int,
        nl: c_int,
        nr: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        leftbasis: *mut c_double,
        leftdbasis: *mut c_double,
        rightbasis: *mut c_double,
        rightdbasis: *mut c_double,
        normal: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn pda_lj_assemble(
        nl: c_int,
        nr: c_int,
        gamma: c_double,
        he: c_double,
        leftbasis: *mut c_double,
        leftdbasis: *mut c_double,
        rightbasis: *mut c_double,
        rightdbasis: *mut c_double,
        normal: *mut c_double,
        s: c_double,
        ld: c_int,
        stiff: *mut c_double,
    );
    fn pda_ld_setup(n: c_int, nl: c_int, nr: c_int, he: *mut c_double);
    fn pda_ld_ip_begin() -> c_int;
    fn pda_ld_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn pda_ld_point(
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
        n: c_int,
        nl: c_int,
        nr: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        leftbasis: *mut c_double,
        leftdbasis: *mut c_double,
        rightbasis: *mut c_double,
        rightdbasis: *mut c_double,
        normal: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn pda_ld_assemble(
        nl: c_int,
        nr: c_int,
        gamma: c_double,
        he: c_double,
        leftbasis: *mut c_double,
        leftdbasis: *mut c_double,
        rightbasis: *mut c_double,
        rightdbasis: *mut c_double,
        normal: *mut c_double,
        s: c_double,
        ld: c_int,
        stiff: *mut c_double,
    );
    fn pda_default_update(ld: c_int, stiff: *mut c_double, lf: c_int, force: *mut c_double);
    fn pda_default_update_edge(ld: c_int, stiff: *mut c_double, lf: c_int, force: *mut c_double);
    fn pda_default_finish();
    fn pda_default_dirichlet();
    fn pda_default_solve() -> c_double;
    fn pda_fatal_alloc();
}

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut PoissonDgState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(4 * need * need, 0.0);
    state.load.resize(need, 0.0);
    state.force.resize(need, 0.0);
    state.eps.resize(need, 0.0);
    state.n = need;
}

/// First-mesh SAVE state: `STIFF(2N,2N)`, `LOAD(N)`, `FORCE(N)`,
/// `EpsilonBoundary(N)` plus the `AllocationsDone` flag. `N` is the first
/// call's `2*MaxElementDOFs`.
struct PoissonDgState {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    eps: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl PoissonDgState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            eps: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<PoissonDgState> = Mutex::new(PoissonDgState::new());

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
    unsafe { pda_begin(solver) };

    // Mesh => GetMesh() — kept as a touch call.
    unsafe { pda_touch_mesh() };

    // First-mesh permanent storage (sized once, never resized).
    if !state.allocations_done {
        let n = unsafe { pda_first_n() } as usize;
        state.n = n;
        if state.stiff.try_reserve(4 * n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.eps.try_reserve(n).is_err()
        {
            unsafe { pda_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(4 * n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.eps.resize(n, 0.0);
        state.allocations_done = true;
    }

    // Bulk assembly.
    let active = unsafe { pda_nof_active() };
    unsafe { pda_default_initialize() };
    for t in 1..=active {
        let mut n = 0;
        unsafe { pda_begin_active(t, &mut n) };
        ensure_state(&mut state, n as usize);
        let (ld, lf) = (2 * state.n as c_int, state.n as c_int);

        unsafe { pda_bodyforce() };
        let mut found = 0;
        unsafe { pda_load_source(n, state.load.as_mut_ptr(), &mut found) };
        let _ = found;

        local_matrix(&mut state, n);

        unsafe {
            pda_default_update(ld, state.stiff.as_mut_ptr(), lf, state.force.as_mut_ptr());
        }
    }

    // Edge (DG jump) assembly.
    let gamma1: c_double = unsafe { pda_gamma(1) };
    let beta: c_double = unsafe { pda_beta(gamma1) };

    let (mut dim, mut np) = (0, 0);
    unsafe { pda_edge_counts(&mut dim, &mut np) };
    let _ = dim;

    // FORCE = 0.0d0 (whole array, exactly as in the original).
    state.force.fill(0.0);

    for t in 1..=np {
        let (mut n, mut nl, mut nr, mut active_el, mut hasright) = (0, 0, 0, 0, 0);
        unsafe { pda_begin_edge(t, &mut n, &mut nl, &mut nr, &mut active_el, &mut hasright) };
        if active_el == 0 {
            continue;
        }
        if hasright == 0 {
            continue;
        }
        ensure_state(&mut state, (n as usize).max((nl + nr) as usize));
        let (ld, lf) = (2 * state.n as c_int, state.n as c_int);

        let discont = unsafe { pda_edge_discont() } != 0;
        let gamma = if discont { beta } else { gamma1 };

        // STIFF = 0.0d0 (first statement of either LocalJumps variant).
        state.stiff.fill(0.0);

        local_jumps(&mut state, n, nl, nr, gamma, discont);

        unsafe {
            pda_default_update_edge(ld, state.stiff.as_mut_ptr(), lf, state.force.as_mut_ptr());
        }
    }

    // gamma2 is read and discarded in the original (it only feeds the
    // dead LocalMatrixBoundary); kept for statement fidelity.
    let _gamma2: c_double = unsafe { pda_gamma(2) };

    unsafe {
        pda_default_finish();
        pda_default_dirichlet();
    }

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { pda_default_solve() };
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count.
fn local_matrix(state: &mut PoissonDgState, n: c_int) {
    unsafe { pda_lm_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    // Basis(n) / dBasisdx(n,3) automatics, sized by this element.
    let nn = n as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];

    let nip = unsafe { pda_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { pda_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            pda_lm_eleminfo(
                u,
                v,
                w,
                n,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            )
        };
        let _ = stat;

        let load_at_ip: c_double =
            unsafe { pda_dot(n, basis.as_ptr(), state.load.as_ptr()) };

        unsafe {
            pda_lm_assemble(
                n,
                dbasis.as_mut_ptr(),
                s * detj,
                basis.as_mut_ptr(),
                load_at_ip,
                2 * state.n as c_int,
                state.stiff.as_mut_ptr(),
                state.n as c_int,
                state.force.as_mut_ptr(),
            )
        };
    }
}

/// The `LocalJumps` / `LocalJumpsDiscontBC` internal procedures (textually
/// identical bodies, separate SAVE state): DG jump stabilization by edge
/// integration. `discont` selects the wrapper set, mirroring the entry's
/// `'Discont BC'` dispatch.
#[allow(clippy::too_many_arguments)]
fn local_jumps(
    state: &mut PoissonDgState,
    n: c_int,
    nl: c_int,
    nr: c_int,
    gamma: c_double,
    discont: bool,
) {
    let mut he = 0.0;
    if discont {
        unsafe { pda_ld_setup(n, nl, nr, &mut he) };
    } else {
        unsafe { pda_lj_setup(n, nl, nr, &mut he) };
    }

    // Basis(n) / parent-basis automatics, sized by this edge.
    let (nn, nln, nrn) = (n as usize, nl as usize, nr as usize);
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];
    let mut leftbasis = vec![0.0 as c_double; nln];
    let mut leftdbasis = vec![0.0 as c_double; nln * 3];
    let mut rightbasis = vec![0.0 as c_double; nrn];
    let mut rightdbasis = vec![0.0 as c_double; nrn * 3];
    let mut normal = [0.0 as c_double; 3];

    let nip = if discont {
        unsafe { pda_ld_ip_begin() }
    } else {
        unsafe { pda_lj_ip_begin() }
    };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        let (mut detj, mut stat) = (0.0, 0);
        if discont {
            unsafe { pda_ld_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
            unsafe {
                pda_ld_point(
                    &mut u,
                    &mut v,
                    &mut w,
                    &mut s,
                    n,
                    nl,
                    nr,
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    leftbasis.as_mut_ptr(),
                    leftdbasis.as_mut_ptr(),
                    rightbasis.as_mut_ptr(),
                    rightdbasis.as_mut_ptr(),
                    normal.as_mut_ptr(),
                    &mut detj,
                    &mut stat,
                )
            };
        } else {
            unsafe { pda_lj_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
            unsafe {
                pda_lj_point(
                    &mut u,
                    &mut v,
                    &mut w,
                    &mut s,
                    n,
                    nl,
                    nr,
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    leftbasis.as_mut_ptr(),
                    leftdbasis.as_mut_ptr(),
                    rightbasis.as_mut_ptr(),
                    rightdbasis.as_mut_ptr(),
                    normal.as_mut_ptr(),
                    &mut detj,
                    &mut stat,
                )
            };
        }
        let _ = stat;

        if discont {
            unsafe {
                pda_ld_assemble(
                    nl,
                    nr,
                    gamma,
                    he,
                    leftbasis.as_mut_ptr(),
                    leftdbasis.as_mut_ptr(),
                    rightbasis.as_mut_ptr(),
                    rightdbasis.as_mut_ptr(),
                    normal.as_mut_ptr(),
                    s,
                    2 * state.n as c_int,
                    state.stiff.as_mut_ptr(),
                )
            };
        } else {
            unsafe {
                pda_lj_assemble(
                    nl,
                    nr,
                    gamma,
                    he,
                    leftbasis.as_mut_ptr(),
                    leftdbasis.as_mut_ptr(),
                    rightbasis.as_mut_ptr(),
                    rightdbasis.as_mut_ptr(),
                    normal.as_mut_ptr(),
                    s,
                    2 * state.n as c_int,
                    state.stiff.as_mut_ptr(),
                )
            };
        }
    }
}
