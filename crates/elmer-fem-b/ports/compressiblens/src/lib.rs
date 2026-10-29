//! Port of `fem/src/modules/CompressibleNS.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/CompressibleNS.F90` — one F77-global
//! entry `compressiblens_` with one `CONTAINS` internal. The emitted
//! dynamic symbol is exactly `compressiblens_` (proven by `nm -D` on the
//! reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/compressiblens_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. All floating-point work (the
//! LocalMatrix Galerkin assembly, `MATMUL`/`SUM` reductions, the
//! `CPUTime`-based reports) stays in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors CompressibleNS.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: first-mesh sizing (`MaxElementDOFs`,
//!   `(dim+2)*p` blocks) allocated once via `try_reserve`; the internal
//!   procedure's `SAVE Nodes` lives as wrapper-side SAVE state with
//!   identical lifetime (consumed before the next fill, as in the
//!   original). `ALocal`/`SLocal` are allocated but never referenced,
//!   exactly as in the original; `Pressure` is allocated and never
//!   touched, exactly as in the original.
//! - The entry-only zeroing (`Velocity = 0`) runs on every call; the
//!   per-element `STIFF`/`MASS`/`FORCE` zeroing and the `1:nd`-bounded
//!   nodal fills live in the wrappers, so stale tails past `nd` persist
//!   across elements exactly as the original `SAVE` arrays do. `LOAD` is
//!   zeroed over the whole array for each element, then columns 1–3 are
//!   rewritten only when a body force is associated, exactly as in the
//!   original.
//! - `LocalMatrix` forces `PicardIteration = .TRUE.` and
//!   `Stabilization = .FALSE.` on entry, overriding the caller; the
//!   wrapper folds those forced values. The `IF (.FALSE.)` convection
//!   alternative is transcribed verbatim (its `ELSE` arm is live).
//! - `PrevNorm` (set but never read) is not replicated; `Norm` is kept
//!   and discarded. `Model` and `dt` are received and ignored, exactly
//!   as in the original; the transient branch calls the same
//!   `Default1stOrderTime`.
//! - Allocation failure calls the same
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - The `Level=4` timing reports use the same `(a,F8.2)` formats under
//!   the same `NavierStokesSolver` program name; the wall-clock values
//!   are masked by the proof's merged wall-clock norm class.
//! - No Rust threads: the solver entry runs serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn cns_begin(h_solver: *mut c_void);
    fn cns_max_dofs() -> c_int;
    fn cns_coord_dim() -> c_int;
    fn cns_done_time() -> c_int;
    fn cns_nl_params(nliter: *mut c_int, nltol: *mut c_double);
    fn cns_default_start();
    fn cns_default_initialize();
    fn cns_nof_active() -> c_int;
    fn cns_begin_elem(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        nb: *mut c_int,
        idx: *mut c_int,
        psz: c_int,
    );
    fn cns_bodyforce() -> c_int;
    fn cns_load_cols(n: c_int, pld: c_int, load: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn cns_material(
        n: c_int,
        pld: c_int,
        mu: *mut c_double,
        cv: *mut c_double,
        k: *mut c_double,
        gam: *mut c_double,
        lam: *mut c_double,
        t0: *mut c_double,
        rho0: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn cns_elem_values(
        first: c_int,
        newts: c_int,
        iter: c_int,
        d1: c_int,
        pld: c_int,
        nd: c_int,
        velo: *mut c_double,
        temp: *mut c_double,
        rho: *mut c_double,
        idx: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn cns_local_matrix(
        nloc: c_int,
        pld: c_int,
        d1: c_int,
        stiff: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        load: *mut c_double,
        tnod: *mut c_double,
        rnod: *mut c_double,
        mnod: *mut c_double,
        vnod: *mut c_double,
        n: c_int,
        nd: c_int,
        dd: c_int,
        gam: c_double,
        cv: c_double,
        k: c_double,
        la: c_double,
        t0: c_double,
        rho0: c_double,
    );
    fn cns_first_order(nloc: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn cns_default_update(nloc: c_int, stiff: *mut c_double, force: *mut c_double);
    fn cns_finish_bulk();
    fn cns_finish_assembly();
    fn cns_dirichlet();
    fn cns_iter_tail() -> c_double;
    fn cns_nonlin_converged() -> c_int;
    fn cns_default_finish();
    fn cns_fatal_alloc();
}

/// First-mesh SAVE state: `STIFF(nloc,nloc)`, `Mass(nloc,nloc)`,
/// `FORCE(nloc)`, `LOAD(p,5)`, nodal `rho`/`Temperature`/`mu`/`Pressure`,
/// `Velocity(d1,p)`, `ALocal(m+1,m+1)`, `SLocal(2p,2p)`, `Indexes(p)`,
/// plus the `AllocationsDone` flag and the `CurrentDoneTime` counter.
/// `p` is the first call's `MaxElementDOFs`; `nloc = (dim+2)*p`,
/// `m = dim*p`, `d1 = dim+1`.
struct CnsState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    rho: Vec<c_double>,
    temp: Vec<c_double>,
    mu: Vec<c_double>,
    press: Vec<c_double>,
    velo: Vec<c_double>,
    alocal: Vec<c_double>,
    slocal: Vec<c_double>,
    indexes: Vec<c_int>,
    p: usize,
    dim: usize,
    nloc: usize,
    current_done: c_int,
    allocations_done: bool,
}

impl CnsState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            rho: Vec::new(),
            temp: Vec::new(),
            mu: Vec::new(),
            press: Vec::new(),
            velo: Vec::new(),
            alocal: Vec::new(),
            slocal: Vec::new(),
            indexes: Vec::new(),
            p: 0,
            dim: 0,
            nloc: 0,
            current_done: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<CnsState> = Mutex::new(CnsState::new());

/// `compressiblens_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads (they travel by reference, as in
/// the original call sequence).
#[export_name = "compressiblens_"]
pub unsafe extern "C" fn compressiblens_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    let is_transient = unsafe { *transient } != 0;
    unsafe { cns_begin(solver) };

    // NewTimestep / FirstSolve from the SAVE CurrentDoneTime counter.
    let first = i32::from(state.current_done == 0);
    let mut newts = 0;
    if state.current_done != unsafe { cns_done_time() } {
        newts = 1;
        state.current_done += 1;
    }

    // First-mesh permanent storage (bug-compatible: sized once from the
    // first mesh, never resized).
    if !state.allocations_done {
        let p = unsafe { cns_max_dofs() } as usize;
        // dim was captured into the wrapper SAVE by cns_begin.
        let dim = unsafe { cns_coord_dim() } as usize;
        let nloc = (dim + 2) * p;
        let m = dim * p;
        let d1 = dim + 1;
        state.p = p;
        state.dim = dim;
        state.nloc = nloc;
        if state.stiff.try_reserve(nloc * nloc).is_err()
            || state.mass.try_reserve(nloc * nloc).is_err()
            || state.force.try_reserve(nloc).is_err()
            || state.load.try_reserve(p * 5).is_err()
            || state.rho.try_reserve(p).is_err()
            || state.temp.try_reserve(p).is_err()
            || state.mu.try_reserve(p).is_err()
            || state.press.try_reserve(p).is_err()
            || state.velo.try_reserve(d1 * p).is_err()
            || state.alocal.try_reserve((m + 1) * (m + 1)).is_err()
            || state.slocal.try_reserve(2 * p * 2 * p).is_err()
            || state.indexes.try_reserve(p).is_err()
        {
            unsafe { cns_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(nloc * nloc, 0.0);
        state.mass.resize(nloc * nloc, 0.0);
        state.force.resize(nloc, 0.0);
        state.load.resize(p * 5, 0.0);
        state.rho.resize(p, 0.0);
        state.temp.resize(p, 0.0);
        state.mu.resize(p, 0.0);
        state.press.resize(p, 0.0);
        state.velo.resize(d1 * p, 0.0);
        state.alocal.resize((m + 1) * (m + 1), 0.0);
        state.slocal.resize(2 * p * 2 * p, 0.0);
        state.indexes.resize(p, 0);
        state.allocations_done = true;
    }

    let (p, dim, nloc) = (state.p, state.dim, state.nloc);
    let d1 = (dim + 1) as c_int;

    // Velocity = 0 on every entry, exactly as in the original.
    state.velo.fill(0.0);

    let (mut nliter, mut nltol) = (0, 0.0);
    unsafe { cns_nl_params(&mut nliter, &mut nltol) };
    let _ = nltol;

    unsafe { cns_default_start() };

    for iter in 1..=nliter {
        unsafe { cns_default_initialize() };

        let nelem = unsafe { cns_nof_active() };
        for t in 1..=nelem {
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe {
                cns_begin_elem(
                    t,
                    &mut n,
                    &mut nd,
                    &mut nb,
                    state.indexes.as_mut_ptr(),
                    p as c_int,
                )
            };
            let _ = nb;

            // LOAD = 0 over the whole array each element, then columns 1-3
            // only when a body force is associated (mirrors the original).
            state.load.fill(0.0);
            if unsafe { cns_bodyforce() } != 0 {
                unsafe { cns_load_cols(n, p as c_int, state.load.as_mut_ptr()) };
            }

            let (mut cv, mut k, mut gam, mut lam, mut t0, mut rho0) =
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            unsafe {
                cns_material(
                    n,
                    p as c_int,
                    state.mu.as_mut_ptr(),
                    &mut cv,
                    &mut k,
                    &mut gam,
                    &mut lam,
                    &mut t0,
                    &mut rho0,
                )
            };

            let stiff = state.stiff.as_mut_ptr();
            let mass = state.mass.as_mut_ptr();
            let force = state.force.as_mut_ptr();
            let load = state.load.as_mut_ptr();
            let temp = state.temp.as_mut_ptr();
            let rho = state.rho.as_mut_ptr();
            let mu = state.mu.as_mut_ptr();
            let velo = state.velo.as_mut_ptr();
            let idx = state.indexes.as_mut_ptr();
            unsafe {
                cns_elem_values(
                    first, newts, iter, d1, p as c_int, nd, velo, temp, rho, idx,
                );
                cns_local_matrix(
                    nloc as c_int,
                    p as c_int,
                    d1,
                    stiff,
                    mass,
                    force,
                    load,
                    temp,
                    rho,
                    mu,
                    velo,
                    n,
                    nd,
                    dim as c_int,
                    gam,
                    cv,
                    k,
                    lam,
                    t0,
                    rho0,
                );
                if is_transient {
                    cns_first_order(nloc as c_int, mass, stiff, force);
                }
                cns_default_update(nloc as c_int, stiff, force);
            }
        }

        unsafe {
            cns_finish_bulk();
            cns_finish_assembly();
            cns_dirichlet();
        }

        let _norm: c_double = unsafe { cns_iter_tail() };

        if unsafe { cns_nonlin_converged() } != 0 {
            break;
        }
    }

    unsafe { cns_default_finish() };
}
