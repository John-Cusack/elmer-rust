//! Port of `fem/src/modules/HelmholtzProjection.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/HelmholtzProjection.F90` — four
//! F77-global entries (`helmholtzprojector_`,
//! `helmholtzprojector_init_`, `removekernelcomponent_`,
//! `removekernelcomponent_init0_`) each with one `CONTAINS` internal
//! (`LocalMatrix`). The emitted dynamic symbols are exactly those four
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/helmholtzprojection_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`), whole-array updates, and the final permutation
//! redefinition stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors HelmholtzProjection.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: first-mesh sizing (`MaxElementDOFs`)
//!   with one `AllocationsDone` flag per entry (`Stiff`/`Force`/`PotSol`
//!   for the projector; plus `Sol`/`PhiSol` for the kernel remover).
//!   The internal procedures' `SAVE Nodes` live as wrapper-side SAVE
//!   state with identical lifetimes (consumed before the next fill, as in
//!   the original). `Sol` is allocated but never read, exactly as in the
//!   original (kept so the allocation sequence matches).
//! - Allocation failure calls the same `Fatal(...)` via `try_reserve`
//!   (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies except
//!   through the lookup loops over `Model % NumberOfSolvers`, exactly as
//!   in the original (received and passed opaquely); `Mesh => GetMesh()`
//!   is kept as a touch call; `Norm = DefaultSolve()` is kept and
//!   discarded.
//! - `HelmholtzProjector` deliberately omits `DefaultFinishAssembly`
//!   (only bulk finish + Dirichlet), while `RemoveKernelComponent` calls
//!   both, exactly as in the original.
//! - The projector's `ElementInfo` failure flag (`Stat`) is returned but
//!   unused, exactly as in the original (no check on it).
//! - No Rust threads: the solver entries run serially, so each `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn hp_begin(h_model: *const c_void, h_solver: *const c_void);
    fn hp_default_start();
    fn hp_dim() -> c_int;
    fn hp_params_begin() -> c_int;
    fn hp_constbulk() -> c_int;
    fn hp_ready() -> c_int;
    fn hp_lookup();
    fn hp_nof_active() -> c_int;
    fn hp_default_initialize(ready: c_int);
    fn hp_elem_begin(t: c_int, n: *mut c_int, npot: *mut c_int, ndpot: *mut c_int);
    fn hp_get_potsol(buf: *mut c_double, nb: c_int);
    fn hp_nodes_touch();
    fn hp_ip_begin() -> c_int;
    fn hp_check_2d();
    fn hp_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn hp_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int, ne: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        wbasis: *mut c_double, curlw: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn hp_proj_point(
        ne: c_int, dim: c_int, npot: c_int, ndpot: c_int,
        potsol: *const c_double, sdet: c_double,
        n: c_int, ld: c_int,
        dbasis: *mut c_double, wbasis: *mut c_double,
        stiff: *mut c_double, force: *mut c_double,
        ready: c_int,
    );
    fn hp_default_update_equations(ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn hp_default_update_force(ld: c_int, force: *mut c_double);
    fn hp_finish_bulk(cb: c_int, ready: c_int);
    fn hp_default_dirichlet();
    fn hp_zero_values();
    fn hp_default_solve() -> c_double;
    fn hp_fatal_alloc();
    fn hp_init_begin(h_solver: *const c_void);
    fn hp_init_body();
    fn rk_begin(h_model: *const c_void, h_solver: *const c_void);
    fn rk_default_start();
    fn rk_dim() -> c_int;
    fn rk_params_begin() -> c_int;
    fn rk_constbulk() -> c_int;
    fn rk_ready() -> c_int;
    fn rk_lookup();
    fn rk_nof_active() -> c_int;
    fn rk_default_initialize(ready: c_int);
    fn rk_elem_begin(
        t: c_int, n: *mut c_int, nd: *mut c_int,
        npot: *mut c_int, ndpot: *mut c_int,
    );
    fn rk_fatal_dofcount();
    fn rk_get_potsol(buf: *mut c_double, nb: c_int);
    fn rk_get_phisol(buf: *mut c_double, nb: c_int);
    fn rk_nodes_touch();
    fn rk_ip_begin() -> c_int;
    fn rk_check_2d();
    fn rk_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn rk_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        wbasis: *mut c_double, curlw: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn rk_point(
        nd: c_int, dim: c_int, n: c_int,
        phisol: *const c_double, sdet: c_double,
        dbasis: *mut c_double, wbasis: *mut c_double,
        stiff: *mut c_double, force: *mut c_double,
        ready: c_int, ld: c_int,
    );
    fn rk_default_update_equations(ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn rk_default_update_force(ld: c_int, force: *mut c_double);
    fn rk_finish_bulk(cb: c_int, ready: c_int);
    fn rk_default_finish();
    fn rk_default_dirichlet();
    fn rk_zero_values();
    fn rk_default_solve() -> c_double;
    fn rk_perm_size(which: c_int) -> c_int;
    fn rk_perm_val(which: c_int, i: c_int) -> c_int;
    fn rk_values_get(which: c_int, j: c_int) -> c_double;
    fn rk_values_set(which: c_int, j: c_int, v: c_double);
    fn rk_fatal_nonmatch();
    fn rk_fatal_differ();
    fn rk_fatal_alloc();
    fn rk_init0_begin(h_solver: *const c_void);
    fn rk_init0_body();
}

/// First-mesh SAVE state of `HelmholtzProjector`: `Stiff(N,N)`,
/// `Force(N)`, `PotSol(1,N)` plus the `AllocationsDone` flag.
struct HpState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    potsol: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl HpState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            potsol: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

/// First-mesh SAVE state of `RemoveKernelComponent`: `Stiff(N,N)`,
/// `Force(N)`, `Sol(1,N)`, `PotSol(1,N)`, `PhiSol(N)` plus the
/// `AllocationsDone` flag.
struct RkState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    sol: Vec<c_double>,
    potsol: Vec<c_double>,
    phisol: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl RkState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            sol: Vec::new(),
            potsol: Vec::new(),
            phisol: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static HP_STATE: Mutex<HpState> = Mutex::new(HpState::new());
static RK_STATE: Mutex<RkState> = Mutex::new(RkState::new());

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn hp_ensure(state: &mut HpState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.force.resize(need, 0.0);
    state.potsol.resize(need, 0.0);
    state.n = need;
}

/// Same grow-once guard for the kernel remover's SAVE arrays.
fn rk_ensure(state: &mut RkState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.force.resize(need, 0.0);
    state.sol.resize(need, 0.0);
    state.potsol.resize(need, 0.0);
    state.phisol.resize(need, 0.0);
    state.n = need;
}

/// `helmholtzprojector_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `HelmholtzProjector_Init` (linear-system
/// refactorize default).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "helmholtzprojector_init_"]
pub unsafe extern "C" fn helmholtzprojector_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        hp_init_begin(solver);
        hp_init_body();
    }
}

/// `helmholtzprojector_` — compute the H1-regular scalar field whose
/// gradient is the kernel component of the projected curl-conforming
/// field. Argument order and by-reference travel match the Fortran dummy
/// list `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "helmholtzprojector_"]
pub unsafe extern "C" fn helmholtzprojector_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = HP_STATE.lock().unwrap();
    unsafe { hp_begin(model, solver) };

    let dim = unsafe { hp_dim() };
    let maxdofs = unsafe { hp_params_begin() } as usize;

    // Allocate some permanent storage, this is done first time only.
    if !state.allocations_done {
        state.n = maxdofs;
        if state.stiff.try_reserve(maxdofs * maxdofs).is_err()
            || state.force.try_reserve(maxdofs).is_err()
            || state.potsol.try_reserve(maxdofs).is_err()
        {
            unsafe { hp_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(maxdofs * maxdofs, 0.0);
        state.force.resize(maxdofs, 0.0);
        state.potsol.resize(maxdofs, 0.0);
        state.allocations_done = true;
    }

    unsafe { hp_default_start() };

    // Check if accelerated assembly is desired.
    let constbulk = unsafe { hp_constbulk() };
    let ready = unsafe { hp_ready() };

    // Find the variable which is projected (Fatal paths inside).
    unsafe { hp_lookup() };

    // System assembly.
    let active = unsafe { hp_nof_active() };
    unsafe { hp_default_initialize(ready) };

    for t in 1..=active {
        let (mut n, mut npot, mut ndpot) = (0, 0, 0);
        unsafe { hp_elem_begin(t, &mut n, &mut npot, &mut ndpot) };
        hp_ensure(&mut state, n as usize);
        let ld = state.n as c_int;
        let ne = ndpot - npot;

        unsafe { hp_get_potsol(state.potsol.as_mut_ptr(), ld) };

        // Get element local matrix and rhs vector (LocalMatrix).
        unsafe { hp_nodes_touch() };
        state.stiff.fill(0.0);
        state.force.fill(0.0);

        let nn = n as usize;
        let nen = ne as usize;
        let mut basis = vec![0.0 as c_double; nn];
        let mut dbasis = vec![0.0 as c_double; nn * 3];
        let mut wbasis = vec![0.0 as c_double; nen * 3];
        let mut curlw = vec![0.0 as c_double; nen * 3];

        let nip = unsafe { hp_ip_begin() };
        unsafe { hp_check_2d() };
        for t2 in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { hp_ip_point(t2, &mut u, &mut v, &mut w, &mut s) };

            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                hp_eleminfo(
                    u, v, w, n, n, ne,
                    basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                    wbasis.as_mut_ptr(), curlw.as_mut_ptr(),
                    &mut detj, &mut stat,
                )
            };
            let _ = stat;

            unsafe {
                hp_proj_point(
                    ne, dim, npot, ndpot,
                    state.potsol.as_ptr(), s * detj,
                    n, ld,
                    dbasis.as_mut_ptr(), wbasis.as_mut_ptr(),
                    state.stiff.as_mut_ptr(), state.force.as_mut_ptr(),
                    ready,
                )
            };
        }

        // Update global matrix and rhs vector from local matrix & vector.
        unsafe {
            if ready != 0 {
                hp_default_update_force(ld, state.force.as_mut_ptr());
            } else {
                hp_default_update_equations(
                    ld,
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                );
            }
        }
    }

    unsafe { hp_finish_bulk(constbulk, ready) };

    // In this case no need to call DefaultFinishAssembly.
    unsafe {
        hp_default_dirichlet();
    }

    // Zero values may be a better initial guess than the previous solution.
    unsafe { hp_zero_values() };

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { hp_default_solve() };
}

/// `removekernelcomponent_init0_` — the `_Init0` entry Elmer calls once at
/// solver setup. Mirrors `RemoveKernelComponent_Init0` (linear-system
/// refactorize default plus the automatic `Element` definition from the
/// edge-element style check).
///
/// # Safety
///
/// Same contract as [`helmholtzprojector_init_entry`].
#[export_name = "removekernelcomponent_init0_"]
pub unsafe extern "C" fn removekernelcomponent_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        rk_init0_begin(solver);
        rk_init0_body();
    }
}

/// `removekernelcomponent_` — apply the Helmholtz projection by
/// subtracting the kernel component (grad phi in curl-conforming basis)
/// from the projected field. Argument order and by-reference travel match
/// the Fortran dummy list `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`helmholtzprojector_entry`].
#[export_name = "removekernelcomponent_"]
pub unsafe extern "C" fn removekernelcomponent_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = RK_STATE.lock().unwrap();
    unsafe { rk_begin(model, solver) };

    let dim = unsafe { rk_dim() };
    let maxdofs = unsafe { rk_params_begin() } as usize;

    // Allocate some permanent storage, this is done first time only.
    if !state.allocations_done {
        state.n = maxdofs;
        if state.stiff.try_reserve(maxdofs * maxdofs).is_err()
            || state.force.try_reserve(maxdofs).is_err()
            || state.sol.try_reserve(maxdofs).is_err()
            || state.potsol.try_reserve(maxdofs).is_err()
            || state.phisol.try_reserve(maxdofs).is_err()
        {
            unsafe { rk_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(maxdofs * maxdofs, 0.0);
        state.force.resize(maxdofs, 0.0);
        state.sol.resize(maxdofs, 0.0);
        state.potsol.resize(maxdofs, 0.0);
        state.phisol.resize(maxdofs, 0.0);
        state.allocations_done = true;
    }

    unsafe { rk_default_start() };

    // Check if accelerated assembly is desired.
    let constbulk = unsafe { rk_constbulk() };
    let ready = unsafe { rk_ready() };

    // Find the projected and kernel variables (Fatal paths inside).
    unsafe { rk_lookup() };

    // System assembly.
    let active = unsafe { rk_nof_active() };
    unsafe { rk_default_initialize(ready) };

    for t in 1..=active {
        let (mut n, mut nd, mut npot, mut ndpot) = (0, 0, 0, 0);
        unsafe { rk_elem_begin(t, &mut n, &mut nd, &mut npot, &mut ndpot) };

        // The DOF counts for the potential variable.
        if nd != ndpot - npot {
            unsafe { rk_fatal_dofcount() };
        }
        rk_ensure(&mut state, (n as usize).max(nd as usize));
        let ld = state.n as c_int;

        unsafe {
            rk_get_potsol(state.potsol.as_mut_ptr(), ld);
            rk_get_phisol(state.phisol.as_mut_ptr(), ld);
        }

        // Get element local matrix and rhs vector (LocalMatrix).
        unsafe { rk_nodes_touch() };
        state.stiff.fill(0.0);
        state.force.fill(0.0);

        let nn = n as usize;
        let ndn = nd as usize;
        let mut basis = vec![0.0 as c_double; nn];
        let mut dbasis = vec![0.0 as c_double; nn * 3];
        let mut wbasis = vec![0.0 as c_double; ndn * 3];
        let mut curlw = vec![0.0 as c_double; ndn * 3];

        let nip = unsafe { rk_ip_begin() };
        unsafe { rk_check_2d() };
        for t2 in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { rk_ip_point(t2, &mut u, &mut v, &mut w, &mut s) };

            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                rk_eleminfo(
                    u, v, w, n, nd,
                    basis.as_mut_ptr(), dbasis.as_mut_ptr(),
                    wbasis.as_mut_ptr(), curlw.as_mut_ptr(),
                    &mut detj, &mut stat,
                )
            };
            let _ = stat;

            unsafe {
                rk_point(
                    nd, dim, n,
                    state.phisol.as_ptr(), s * detj,
                    dbasis.as_mut_ptr(), wbasis.as_mut_ptr(),
                    state.stiff.as_mut_ptr(), state.force.as_mut_ptr(),
                    ready, ld,
                )
            };
        }

        // Update global matrix and rhs vector from local matrix & vector.
        unsafe {
            if ready != 0 {
                rk_default_update_force(ld, state.force.as_mut_ptr());
            } else {
                rk_default_update_equations(
                    ld,
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                );
            }
        }
    }

    unsafe {
        rk_finish_bulk(constbulk, ready);
        rk_default_finish();
        rk_default_dirichlet();
    }

    // Zero values may be a better initial guess than the previous solution.
    unsafe { rk_zero_values() };

    let _norm: c_double = unsafe { rk_default_solve() };

    // Finally, redefine the potential variable: A := A - W.
    let n1 = unsafe { rk_perm_size(1) };
    let n2 = unsafe { rk_perm_size(2) };
    if n1 != n2 {
        unsafe { rk_fatal_differ() };
    }
    for i in 1..=n1 {
        let j = unsafe { rk_perm_val(1, i) };
        if j < 1 {
            continue;
        }
        let k = unsafe { rk_perm_val(2, i) };
        if k < 1 {
            unsafe { rk_fatal_nonmatch() };
        }
        let v = unsafe { rk_values_get(2, k) - rk_values_get(1, j) };
        unsafe { rk_values_set(2, k, v) };
    }
}
