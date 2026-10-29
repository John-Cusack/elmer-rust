//! Port of `fem/src/modules/TransportEquation.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/TransportEquation.F90` — one F77-global
//! entry `transportequationsolver_` (proven by `nm -D` on the reference
//! `.so`); the `CONTAINS` internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape, T020 Poisson pattern):
//! the entry's derived-type dummies (`Model`, `Solver`) cross as opaque
//! addresses. Every component access and every Elmer library call goes
//! through the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/transportequation_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), formatted I/O (`WRITE`), and whole-array updates stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors TransportEquationSolver statement for
//! statement):
//! - `SAVE` semantics replicated: `ThisSolverCalls` is a Rust-side
//!   `Mutex`-guarded counter (even calls return immediately, exactly like
//!   the original's `MOD(ThisSolverCalls,2)` gate); the first-mesh sizing
//!   (`N = Model % MaxElementNodes`, `D = Solver % Variable % DOFs`)
//!   allocates the `LocalStiffMatrix/LocalForce/LocalMassMatrix/
//!   LocalDampMatrix/V1/V2/V3` scratch once (`AllocationsDone`), while
//!   `ElementNodes`, `U0Var`, and the per-internal-procedure Gauss-point
//!   states live as wrapper-side SAVE with identical lifetimes.
//! - `Add2ndOrderTime2`'s `NodeIndexes` dummy receives `UPerm(NodeIndexes)`
//!   (permuted equation numbers, not element nodes) at both call sites;
//!   the wrapper reads `Solver % Variable % Perm` at
//!   `acc_elem % NodeIndexes(i)` for the same `L`, never the raw node id.
//! - Allocation failure calls the same
//!   `Fatal('TransportEquationSolver', 'Memory allocation error,
//!   aborting.')` via `try_reserve` (the `STAT=istat` check), never
//!   aborting directly.
//! - `Add2ndOrderTime2`'s uninitialized `Method` `SELECT CASE` is
//!   transcribed verbatim in the wrapper: its only case is `DEFAULT`, so
//!   the indeterminate value is unobservable and both arms take the same
//!   path deterministically.
//! - `PrevValues(:,3/4/5)` (X/V/A) require the proof configs to set
//!   `Time Derivative Order = 2` (Elmer core then sizes `PrevValues` to 7
//!   columns); without it the reference itself reads out of bounds. The
//!   port stays on this defined path, like the reference.
//! - `Model` and `TransientSimulation` are unused by the body, exactly as
//!   in the original (received and ignored); `Norm = DefaultSolve()` is
//!   kept and discarded.
//! - The ` Assembly: (s)` / ` Solve:    (s)` `Info` lines print process
//!   CPU times, so the proof script norms them out as wall-clock lines.
//! - No Rust threads: the solver entry runs serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn teq_begin(h_model: *const c_void, h_solver: *const c_void);
    fn teq_info_init();
    fn teq_var_dofs() -> c_int;
    fn teq_max_elem_nodes() -> c_int;
    fn teq_init_cond(dt: c_double);
    fn teq_at_start();
    fn teq_at_done();
    fn teq_init_zero();
    fn teq_init_timestep();
    fn teq_eqname();
    fn teq_adv_setup(flag: *mut c_int, velcomps: *mut c_int);
    fn teq_num_bulk() -> c_int;
    fn teq_num_bnd() -> c_int;
    fn teq_begin_bulk(t: c_int, n: *mut c_int, proceed: *mut c_int);
    fn teq_begin_bnd(t: c_int, n: *mut c_int, proceed: *mut c_int);
    fn teq_elem_nodes(n: c_int);
    fn teq_fill_vel(n: c_int, nn: c_int, v1: *mut c_double, v2: *mut c_double, v3: *mut c_double);
    fn teq_lm_zero(
        ld: c_int,
        stiff: *mut c_double,
        mass: *mut c_double,
        damp: *mut c_double,
        force: *mut c_double,
    );
    fn teq_lm_ip_begin() -> c_int;
    fn teq_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn teq_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn teq_lm_ip(
        ld: c_int,
        n: c_int,
        vc: c_int,
        dim: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        sdet: c_double,
        v1: *mut c_double,
        v2: *mut c_double,
        v3: *mut c_double,
        stiff: *mut c_double,
        mass: *mut c_double,
    );
    fn teq_lmb_ip_begin() -> c_int;
    fn teq_lmb_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn teq_lmb_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn teq_lmb_ip(
        ld: c_int,
        n: c_int,
        vc: c_int,
        dim: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        sdet: c_double,
        u: c_double,
        v: c_double,
        w: c_double,
        v1: *mut c_double,
        v2: *mut c_double,
        v3: *mut c_double,
        damp: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn teq_add2nd(
        ld: c_int,
        mass: *mut c_double,
        damp: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        dt: c_double,
        n: c_int,
        dofs: c_int,
    );
    fn teq_update_global(ld: c_int, n: c_int, dofs: c_int, stiff: *mut c_double, force: *mut c_double);
    fn teq_finish();
    fn teq_dirichlet();
    fn teq_solve() -> c_double;
    fn teq_st_start();
    fn teq_print_times();
    fn teq_restore_dt(dt: c_double);
    fn teq_post();
    fn teq_fatal_alloc();
}

/// Grow the scratch buffers when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut TransState, ld: usize, nmax: usize) {
    if ld <= state.ld {
        return;
    }
    let n2 = ld * ld;
    state.stiff.resize(n2, 0.0);
    state.mass.resize(n2, 0.0);
    state.damp.resize(n2, 0.0);
    state.force.resize(ld, 0.0);
    state.v1.resize(nmax, 0.0);
    state.v2.resize(nmax, 0.0);
    state.v3.resize(nmax, 0.0);
    state.ld = ld;
    state.nmax = nmax;
}

/// First-mesh SAVE state: `LocalStiffMatrix/LocalMassMatrix/LocalDampMatrix((N*D)^2)`,
/// `LocalForce(N*D)`, `V1/V2/V3(N)` plus the `AllocationsDone` flag. `N` is the
/// first call's `MaxElementNodes`, `D` the solver variable's `DOFs`.
struct TransState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    damp: Vec<c_double>,
    force: Vec<c_double>,
    v1: Vec<c_double>,
    v2: Vec<c_double>,
    v3: Vec<c_double>,
    ld: usize,
    nmax: usize,
    dofs: usize,
    allocations_done: bool,
    calls: u32,
}

impl TransState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            damp: Vec::new(),
            force: Vec::new(),
            v1: Vec::new(),
            v2: Vec::new(),
            v3: Vec::new(),
            ld: 0,
            nmax: 0,
            dofs: 0,
            allocations_done: false,
            calls: 0,
        }
    }
}

static STATE: Mutex<TransState> = Mutex::new(TransState::new());

/// `transportequationsolver_` — the convective transport solver entry.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt` must be a valid read (the time step, used scaled).
#[export_name = "transportequationsolver_"]
pub unsafe extern "C" fn transport_equation_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { teq_begin(model, solver) };

    // IF (MOD(ThisSolverCalls,2)==1) ... END IF — even calls are a no-op.
    state.calls += 1;
    if state.calls.is_multiple_of(2) {
        return;
    }

    let dtv: c_double = unsafe { *dt };

    unsafe { teq_info_init() };

    // First-mesh permanent storage (bug-compatible: sized once, never resized).
    if !state.allocations_done {
        let nmax = unsafe { teq_max_elem_nodes() } as usize;
        let dofs = unsafe { teq_var_dofs() } as usize;
        let ld = nmax * dofs;
        if state.stiff.try_reserve(ld * ld).is_err()
            || state.mass.try_reserve(ld * ld).is_err()
            || state.damp.try_reserve(ld * ld).is_err()
            || state.force.try_reserve(ld).is_err()
            || state.v1.try_reserve(nmax).is_err()
            || state.v2.try_reserve(nmax).is_err()
            || state.v3.try_reserve(nmax).is_err()
        {
            unsafe { teq_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(ld * ld, 0.0);
        state.mass.resize(ld * ld, 0.0);
        state.damp.resize(ld * ld, 0.0);
        state.force.resize(ld, 0.0);
        state.v1.resize(nmax, 0.0);
        state.v2.resize(nmax, 0.0);
        state.v3.resize(nmax, 0.0);
        state.ld = ld;
        state.nmax = nmax;
        state.dofs = dofs;
        state.allocations_done = true;
    }
    let dofs = state.dofs as c_int;

    // Initial conditions + Solver % dt = 2*dt (original lines 154-236).
    unsafe { teq_init_cond(dtv) };

    unsafe { teq_at_start() };
    unsafe { teq_init_zero() };
    unsafe { teq_init_timestep() };
    unsafe { teq_eqname() };

    let (mut advflag, mut velcomps) = (0, 0);
    unsafe { teq_adv_setup(&mut advflag, &mut velcomps) };
    let _ = advflag;

    // Bulk assembly (original lines 249-323).
    let nb = unsafe { teq_num_bulk() };
    for t in 1..=nb {
        let (mut n, mut proceed) = (0, 0);
        unsafe { teq_begin_bulk(t, &mut n, &mut proceed) };
        if proceed == 0 {
            continue;
        }
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        let need = state.nmax.max(n as usize);
        ensure_state(&mut state, (n as usize) * (dofs as usize), need);
        let ld = state.ld as c_int;
        unsafe { teq_elem_nodes(n) };
        unsafe {
            teq_fill_vel(
                n,
                state.nmax as c_int,
                state.v1.as_mut_ptr(),
                state.v2.as_mut_ptr(),
                state.v3.as_mut_ptr(),
            )
        };
        unsafe {
            teq_lm_zero(
                ld,
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.damp.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
        local_matrix(&mut state, n, velcomps);
        unsafe {
            teq_add2nd(
                ld,
                state.mass.as_mut_ptr(),
                state.damp.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                2.0 * dtv,
                n,
                dofs,
            );
            teq_update_global(ld, n, dofs, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    // Boundary assembly (original lines 328-412).
    let nbn = unsafe { teq_num_bnd() };
    for t in nb + 1..=nb + nbn {
        let (mut n, mut proceed) = (0, 0);
        unsafe { teq_begin_bnd(t, &mut n, &mut proceed) };
        if proceed == 0 {
            continue;
        }
        let need = state.nmax.max(n as usize);
        ensure_state(&mut state, (n as usize) * (dofs as usize), need);
        let ld = state.ld as c_int;
        unsafe { teq_elem_nodes(n) };
        unsafe {
            teq_fill_vel(
                n,
                state.nmax as c_int,
                state.v1.as_mut_ptr(),
                state.v2.as_mut_ptr(),
                state.v3.as_mut_ptr(),
            )
        };
        unsafe {
            teq_lm_zero(
                ld,
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.damp.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
        local_matrix_boundary(&mut state, n, velcomps);
        unsafe {
            teq_add2nd(
                ld,
                state.mass.as_mut_ptr(),
                state.damp.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                2.0 * dtv,
                n,
                dofs,
            );
            teq_update_global(ld, n, dofs, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe { teq_finish() };
    unsafe { teq_dirichlet() };

    unsafe { teq_at_done() };
    unsafe { teq_st_start() };

    // And finally, solve (Norm is local and unused, as in the original).
    let _norm: c_double = unsafe { teq_solve() };
    unsafe { teq_print_times() };

    unsafe { teq_restore_dt(dtv) };

    // Outflow-BC retrieval + advection-variable update (original 450-473).
    unsafe { teq_post() };
}

/// The `LocalMatrix` internal procedure: element matrix by numerical
/// integration. `n` is the element node count.
fn local_matrix(state: &mut TransState, n: c_int, velcomps: c_int) {
    let nn = n as usize;
    // Basis(n) / dBasisdx(n,3) automatics, sized by this element.
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];
    let ld = state.ld as c_int;

    let nip = unsafe { teq_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { teq_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            teq_lm_eleminfo(
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

        unsafe {
            teq_lm_ip(
                ld,
                n,
                state.dofs as c_int,
                velcomps,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                s * detj,
                state.v1.as_mut_ptr(),
                state.v2.as_mut_ptr(),
                state.v3.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
            )
        };
    }
}

/// The `LocalMatrixBoundary` internal procedure: boundary damping matrix
/// by boundary integration. `n` is the element node count.
fn local_matrix_boundary(state: &mut TransState, n: c_int, velcomps: c_int) {
    let nn = n as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];
    let ld = state.ld as c_int;

    let nip = unsafe { teq_lmb_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { teq_lmb_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            teq_lmb_eleminfo(
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

        unsafe {
            teq_lmb_ip(
                ld,
                n,
                state.dofs as c_int,
                velcomps,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                s * detj,
                u,
                v,
                w,
                state.v1.as_mut_ptr(),
                state.v2.as_mut_ptr(),
                state.v3.as_mut_ptr(),
                state.damp.as_mut_ptr(),
            )
        };
    }
}

