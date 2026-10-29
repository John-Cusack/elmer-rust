//! Port of `fem/src/modules/DFTSolver.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/DFTSolver.F90` — three F77-global
//! entries (`chargedensitysolver_`, `wavefunctionsolver_`,
//! `poissonsolver_`) with four `CONTAINS` internals (one `LocalMatrix`
//! per entry plus `GRPULAY3` inside `PoissonSolver`). The emitted dynamic
//! symbols are exactly those three entries (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/dftm_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`, `DOT_PRODUCT`, complex eigenvector
//! products), formatted I/O (`WRITE(Message,*)`/`PRINT*`/`INFO`/`WARN`
//! with identical caller strings, typos included), the `LUSolve`
//! relaxation solve, `EXP`-based mixing factors, and whole-array updates
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches. The three `LocalMatrix`
//! integrations and `GRPULAY3` are transcribed statement for statement
//! as Fortran module procedures over wrapper-side SAVE state with
//! identical lifetimes; the three entries' control flow (first-time
//! blocks, assembly loops, relaxation dispatch, GRPULAY3 gating) is
//! transcribed here in Rust.
//!
//! Faithfulness notes (mirror DFTSolver.F90 statement for statement):
//! - `SAVE` semantics replicated: first-mesh work arrays (`Indexes`,
//!   `STIFF`, `FORCE`, `MASS`, local potential/charge arrays) sized once
//!   at `MaxElementDOFs`, never resized on valid inputs (a grow-only
//!   guard covers remeshed re-entry, invisible otherwise); every other
//!   SAVE scalar/pointer/array (eigen counts and weights, nuclei tables,
//!   relaxation parameters, `OldSolution`, GRPULAY3's storage) lives
//!   wrapper-side with identical lifetimes.
//! - `ChargeDensitySolver`'s `ALLOCATE` has no `STAT=` check in the
//!   original, so its three arrays use plain `resize` (OOM aborts, as an
//!   unchecked `ALLOCATE` failure does); `WaveFunctionSolver` and
//!   `PoissonSolver` keep the originals' `STAT=istat` checks via
//!   `try_reserve` calling the same `Fatal(..., 'Memory allocation
//!   error.')`.
//! - `PoissonSolver`'s dead `EigsVar => VariableGet(...)` pointer query
//!   (never read afterwards) is skipped: a pure lookup with no side
//!   effects.
//! - `PoissonSolver`'s `RF` local is assigned in the `CASE(1)`/`CASE(2)`
//!   mixing arms only; it is read by the GRPULAY3 gate on every pass.
//!   The wrapper keeps it in SAVE storage (zero-initialized), while the
//!   original's stack slot is undefined on paths that never assign it
//!   (`RelaxType` 0 or 3 on entry). Proof configs always take a mixing
//!   arm (`constant`/`exponential`), so the gate only ever reads a
//!   freshly assigned value; the unassigned paths are not exercised.
//! - The commented-out MPI alternatives (`MPI_ALLREDUCE` etc.) are
//!   skipped exactly as in the original (serial proof tree); the serial
//!   `DOT_PRODUCT` branches are transcribed.
//! - `Norm = DefaultSolve()` is kept and discarded in all three entries,
//!   as in the original. `Model`, `dt`, and the transient flag are unused
//!   by the entry bodies themselves (received and ignored), except that
//!   `Model` travels into the wrappers as the cross-solver handle.
//! - No Rust threads: the solver entries run serially, so each `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn dfm_begin(h_model: *const c_void, h_solver: *const c_void);
    fn dfm_max_elem_dofs() -> c_int;
    fn dfm_nof_active() -> c_int;
    fn dfm_default_initialize();
    fn dfm_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn dfm_default_update_mass(ld: c_int, g: *mut c_double);
    fn dfm_default_finish();
    fn dfm_default_dirichlet();
    fn dfm_default_solve() -> c_double;
    // ChargeDensitySolver wrappers.
    fn dfm_cds_info_start();
    fn dfm_cds_firsttime();
    fn dfm_cds_zero_values();
    fn dfm_cds_begin_active(
        t: c_int,
        indexes: *mut c_int,
        maxn: c_int,
        n: *mut c_int,
        nd: *mut c_int,
    );
    fn dfm_cds_remap(nd: c_int, indexes: *mut c_int, maxn: c_int);
    #[allow(clippy::too_many_arguments)]
    fn dfm_cds_lm(
        stiff: *mut c_double,
        force: *mut c_double,
        indexes: *mut c_int,
        nmax: c_int,
        nd: c_int,
    );
    fn dfm_cds_msg_max();
    fn dfm_cds_msg_min();
    // WaveFunctionSolver wrappers.
    fn dfm_wfs_info_start();
    fn dfm_wfs_firsttime();
    fn dfm_wfs_fatal_alloc();
    fn dfm_wfs_getlocal(hart: *mut c_double, chg: *mut c_double, nmax: c_int);
    fn dfm_wfs_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn dfm_wfs_lm(
        stiff: *mut c_double,
        mass: *mut c_double,
        hart: *mut c_double,
        chg: *mut c_double,
        nmax: c_int,
        nd: c_int,
    );
    fn dfm_wfs_maybe_shift();
    fn dfm_wfs_eig_store();
    fn dfm_wfs_eig_head();
    fn dfm_wfs_eig_report();
    // PoissonSolver wrappers.
    fn dfm_ps_info_start();
    fn dfm_ps_firsttime();
    fn dfm_ps_fatal_alloc();
    fn dfm_ps_relaxtype() -> c_int;
    fn dfm_ps_calculating() -> c_int;
    fn dfm_ps_save_old();
    fn dfm_ps_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn dfm_ps_getlocal(chg: *mut c_double, nmax: c_int);
    #[allow(clippy::too_many_arguments)]
    fn dfm_ps_lm(
        stiff: *mut c_double,
        force: *mut c_double,
        chg: *mut c_double,
        nmax: c_int,
        nd: c_int,
        total: *mut c_double,
    );
    fn dfm_ps_normalize(total: c_double);
    fn dfm_ps_add_nuclei();
    fn dfm_ps_info_norelax();
    fn dfm_ps_info_grp();
    fn dfm_ps_mix_const();
    fn dfm_ps_mix_exp();
    fn dfm_ps_grpulay_gate();
    fn dfm_ps_relchange();
    fn dfm_ps_msg_total(total: c_double);
    fn dfm_ps_msg_nuctotal();
    fn dfm_ps_msg_min();
    fn dfm_ps_msg_max();
}

/// Grow one entry's first-mesh SAVE arrays when an element needs more
/// than the first-mesh sizing (adaptive remesh across invocations). The
/// original never grows (its comment assumes the first mesh is big
/// enough); growth only triggers where the original would overflow
/// silently, and the leading sections used by the wrappers are
/// unaffected, so valid inputs never observe it.
fn grow_n(state_n: &mut usize, need: usize) {
    if need > *state_n {
        *state_n = need;
    }
}

/// First-mesh SAVE state for `ChargeDensitySolver`: `Indexes(N)`,
/// `STIFF(N,N)`, `FORCE(N)` plus the `FirstTime` flag.
struct CdsState {
    indexes: Vec<c_int>,
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    first: bool,
    n: usize,
}

impl CdsState {
    const fn new() -> Self {
        Self {
            indexes: Vec::new(),
            stiff: Vec::new(),
            force: Vec::new(),
            first: true,
            n: 0,
        }
    }
}

/// First-mesh SAVE state for `WaveFunctionSolver`: `STIFF(N,N)`,
/// `MASS(N,N)`, `FORCE(N)`, `LocalHartreePotential(N)`,
/// `LocalChargeDensity(N)` plus the `AllocationsDone` flag.
struct WfsState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    hart: Vec<c_double>,
    chg: Vec<c_double>,
    alloc: bool,
    n: usize,
}

impl WfsState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            hart: Vec::new(),
            chg: Vec::new(),
            alloc: false,
            n: 0,
        }
    }
}

/// First-mesh SAVE state for `PoissonSolver`: `STIFF(N,N)`, `FORCE(N)`,
/// `LocalChargeDensity(N)` plus the `FirstTime` flag.
struct PsState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    chg: Vec<c_double>,
    first: bool,
    n: usize,
}

impl PsState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            chg: Vec::new(),
            first: true,
            n: 0,
        }
    }
}

static CDS_STATE: Mutex<CdsState> = Mutex::new(CdsState::new());
static WFS_STATE: Mutex<WfsState> = Mutex::new(WfsState::new());
static PS_STATE: Mutex<PsState> = Mutex::new(PsState::new());

/// `chargedensitysolver_` — the charge-density entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "chargedensitysolver_"]
pub unsafe extern "C" fn charge_density_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = CDS_STATE.lock().unwrap();
    unsafe { dfm_begin(model, solver) };

    unsafe { dfm_cds_info_start() };

    // FirstTimeOnly block: first-mesh sizing, eigen-count/weight reads,
    // eigen-solver search (all messages inside the wrapper).
    if state.first {
        let n = unsafe { dfm_max_elem_dofs() } as usize;
        state.n = n;
        // The original ALLOCATEs without STAT= here; plain resize aborts
        // identically on OOM.
        state.indexes.resize(n, 0);
        state.stiff.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        unsafe { dfm_cds_firsttime() };
        state.first = false;
    }

    // Solver % Variable % Values = 0; DefaultInitialize; assembly loop.
    unsafe { dfm_cds_zero_values() };
    unsafe { dfm_default_initialize() };

    let active = unsafe { dfm_nof_active() };
    for t in 1..=active {
        let (mut n, mut nd) = (0, 0);
        unsafe {
            dfm_cds_begin_active(
                t,
                state.indexes.as_mut_ptr(),
                state.n as c_int,
                &mut n,
                &mut nd,
            )
        };
        grow_n(&mut state.n, n as usize);
        let nnow = state.n;
        if state.indexes.len() < nnow {
            state.indexes.resize(nnow, 0);
        }
        if state.stiff.len() < nnow * nnow {
            state.stiff.resize(nnow * nnow, 0.0);
        }
        if state.force.len() < nnow {
            state.force.resize(nnow, 0.0);
        }
        let nmax = state.n as c_int;
        state.stiff.fill(0.0);
        state.force.fill(0.0);
        unsafe {
            dfm_cds_remap(nd, state.indexes.as_mut_ptr(), state.n as c_int);
            dfm_cds_lm(
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                state.indexes.as_mut_ptr(),
                nmax,
                nd,
            );
            dfm_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe {
        dfm_default_finish();
        dfm_default_dirichlet();
    }

    // Norm is local and discarded, as in the original.
    let _norm: c_double = unsafe { dfm_default_solve() };

    // Greatest/smallest reports (the node count is read inside).
    unsafe {
        dfm_cds_msg_max();
        dfm_cds_msg_min();
    }
}

/// `wavefunctionsolver_` — the Kohn-Sham eigenproblem entry. Same
/// argument contract as above.
///
/// # Safety
///
/// Same contract as [`charge_density_solver_entry`].
#[export_name = "wavefunctionsolver_"]
pub unsafe extern "C" fn wave_function_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = WFS_STATE.lock().unwrap();
    unsafe { dfm_begin(model, solver) };

    unsafe { dfm_wfs_info_start() };

    // FirstTimeOnly block: first-mesh sizing (STAT= check), XC-type
    // read, shift read, eigen-count read, EigStorage allocation.
    if !state.alloc {
        let n = unsafe { dfm_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.hart.try_reserve(n).is_err()
            || state.chg.try_reserve(n).is_err()
        {
            unsafe { dfm_wfs_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.hart.resize(n, 0.0);
        state.chg.resize(n, 0.0);
        unsafe { dfm_wfs_firsttime() };
        state.alloc = true;
    }

    unsafe { dfm_default_initialize() };

    state.hart.fill(0.0);
    state.chg.fill(0.0);
    state.force.fill(0.0);

    let active = unsafe { dfm_nof_active() };
    for t in 1..=active {
        let (mut n, mut nd) = (0, 0);
        unsafe { dfm_wfs_begin_active(t, &mut n, &mut nd) };
        grow_n(&mut state.n, n as usize);
        let nnow = state.n;
        if state.stiff.len() < nnow * nnow {
            state.stiff.resize(nnow * nnow, 0.0);
        }
        if state.mass.len() < nnow * nnow {
            state.mass.resize(nnow * nnow, 0.0);
        }
        if state.force.len() < nnow {
            state.force.resize(nnow, 0.0);
        }
        if state.hart.len() < nnow {
            state.hart.resize(nnow, 0.0);
        }
        if state.chg.len() < nnow {
            state.chg.resize(nnow, 0.0);
        }
        let nmax = state.n as c_int;
        unsafe {
            dfm_wfs_getlocal(state.hart.as_mut_ptr(), state.chg.as_mut_ptr(), nmax);
        }
        state.stiff.fill(0.0);
        state.mass.fill(0.0);
        unsafe {
            dfm_wfs_lm(
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.hart.as_mut_ptr(),
                state.chg.as_mut_ptr(),
                nmax,
                nd,
            );
            dfm_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            dfm_default_update_mass(nmax, state.mass.as_mut_ptr());
        }
    }

    unsafe {
        dfm_default_finish();
        dfm_default_dirichlet();
        // Eigenvalue shift, previous-eigenvalue store, eigensolve,
        // absolute-change reports.
        dfm_wfs_maybe_shift();
        dfm_wfs_eig_store();
    }
    let _norm: c_double = unsafe { dfm_default_solve() };
    unsafe {
        dfm_wfs_eig_head();
        dfm_wfs_eig_report();
    }
}

/// `poissonsolver_` — the Hartree-potential entry. Same argument
/// contract as above.
///
/// # Safety
///
/// Same contract as [`charge_density_solver_entry`].
#[export_name = "poissonsolver_"]
pub unsafe extern "C" fn poisson_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = PS_STATE.lock().unwrap();
    unsafe { dfm_begin(model, solver) };

    unsafe { dfm_ps_info_start() };

    // FirsTimeOnly block: first-mesh sizing (STAT= check), solver
    // search, nuclei tables, OldSolution, relaxation parameters.
    // `first` stays true through the assembly below and is cleared at
    // the end, exactly as in the original.
    if state.first {
        let n = unsafe { dfm_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.chg.try_reserve(n).is_err()
        {
            unsafe { dfm_ps_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.chg.resize(n, 0.0);
        unsafe { dfm_ps_firsttime() };
    }

    // Previous-solution store (skipped once GRPulay takes over).
    unsafe {
        if dfm_ps_relaxtype() != 3 {
            dfm_ps_save_old();
        }
    }

    // The original's `EigsVar => VariableGet(...)` result is never read
    // afterwards; the pure pointer query is skipped (no side effects).

    unsafe { dfm_default_initialize() };

    state.chg.fill(1.0);
    let mut total: c_double = 0.0;

    let active = unsafe { dfm_nof_active() };
    for t in 1..=active {
        let (mut n, mut nd) = (0, 0);
        unsafe { dfm_ps_begin_active(t, &mut n, &mut nd) };
        grow_n(&mut state.n, n as usize);
        let nnow = state.n;
        if state.stiff.len() < nnow * nnow {
            state.stiff.resize(nnow * nnow, 0.0);
        }
        if state.force.len() < nnow {
            state.force.resize(nnow, 0.0);
        }
        if state.chg.len() < nnow {
            state.chg.resize(nnow, 0.0);
        }
        let nmax = state.n as c_int;
        if !state.first {
            unsafe { dfm_ps_getlocal(state.chg.as_mut_ptr(), nmax) };
        }
        state.stiff.fill(0.0);
        state.force.fill(0.0);
        unsafe {
            dfm_ps_lm(
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                state.chg.as_mut_ptr(),
                nmax,
                nd,
                &mut total,
            );
            dfm_default_update(nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe {
        dfm_default_finish();
        dfm_ps_normalize(total);
        dfm_default_dirichlet();
        dfm_ps_add_nuclei();
    }
    let _norm: c_double = unsafe { dfm_default_solve() };

    // Relaxation dispatch, GRPULAY3 call + gate, relaxed-solution
    // bookkeeping, result reports.
    unsafe {
        match dfm_ps_relaxtype() {
            0 => dfm_ps_info_norelax(),
            1 => dfm_ps_mix_const(),
            2 => dfm_ps_mix_exp(),
            3 => dfm_ps_info_grp(),
            _ => {}
        }
        dfm_ps_grpulay_gate();
        if dfm_ps_relaxtype() == 3 && dfm_ps_calculating() != 0 {
            dfm_ps_relchange();
        }
        dfm_ps_msg_total(total);
        dfm_ps_msg_nuctotal();
        dfm_ps_msg_min();
        dfm_ps_msg_max();
    }

    state.first = false;
}
