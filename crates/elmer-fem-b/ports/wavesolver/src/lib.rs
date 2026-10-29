//! Port of `fem/src/modules/WaveSolver.F90` (Elmer, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/WaveSolver.F90` — two F77-global
//! entries `wavesolver_` and `wavesolver_init_` (proven by `nm -D` on
//! the reference `.so`); the `CONTAINS` internals emit only local
//! symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape, T020 Poisson pattern):
//! the entries' derived-type dummies (`Model`, `Solver`) cross as opaque
//! addresses. Every component access and every Elmer library call goes
//! through the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/wavesolver_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`) and whole-array updates stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors WaveSolver.F90 statement for statement):
//! - `SAVE` semantics replicated: the wrapper-side `SAVE` mirrors the
//!   original's internal-procedure `SAVE` (per-procedure `ValueHandle`
//!   sets and `Nodes`) with identical lifetimes; the re-init flag is a
//!   per-iteration local set `.TRUE.` before each element loop, exactly
//!   like `InitHandles`.
//! - `MASS`/`STIFF`/`DAMP`/`FORCE` are per-element automatics in the
//!   original (no `SAVE`); the port reuses zeroed scratch buffers of the
//!   first-mesh sizing, observationally identical.
//! - `Model`, `dt` are unused by the body, exactly as in the original
//!   (received and ignored); `Norm = DefaultSolve()` is kept and
//!   discarded; `dim` is queried and passed through, unused, as in the
//!   original.
//! - `TransientSimulation` travels as the Fortran `LOGICAL` by reference
//!   the loader passes; only its truth value is read.
//! - `IF (.NOT.ASSOCIATED(BC)) RETURN` and the
//!   `OutflowBC/NotPresentAnywhere` early return live in the head
//!   wrapper with the original statement order.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn wva_init();
    fn wva_begin(h_solver: *const c_void);
    fn wva_info_start();
    fn wva_max_iter() -> c_int;
    fn wva_coord_dim() -> c_int;
    fn wva_need_mass() -> c_int;
    fn wva_nonlin_converged() -> c_int;
    fn wva_avg_zero();
    fn wva_default_start();
    fn wva_default_initialize();
    fn wva_nof_active() -> c_int;
    fn wva_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn wva_nof_boundary() -> c_int;
    fn wva_begin_boundary(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn wva_default_finish_bulk();
    fn wva_default_finish_bnd();
    fn wva_default_finish_assem();
    fn wva_default_dirichlet();
    fn wva_default_solve() -> c_double;
    fn wva_default_finish();
    fn wva_lm_init(ph: *mut c_int);
    fn wva_lm_nodes();
    fn wva_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn wva_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn wva_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn wva_lm_ip(
        nd: c_int,
        ld: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        sdet: c_double,
        stiff: *mut c_double,
        damp: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn wva_lm_tail(
        nd: c_int,
        nb: c_int,
        ld: c_int,
        mass: *mut c_double,
        damp: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        transient: c_int,
        needmass: c_int,
    );
    fn wva_lbc_head(ph: *mut c_int, outflow: *mut c_int, proceed: *mut c_int);
    fn wva_lbc_nodes();
    fn wva_lbc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn wva_lbc_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn wva_lbc_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn wva_lbc_ip(
        nd: c_int,
        ld: c_int,
        basis: *mut c_double,
        sdet: c_double,
        damp: *mut c_double,
        force: *mut c_double,
        outflow: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn wva_lbc_tail(
        nd: c_int,
        ld: c_int,
        mass: *mut c_double,
        damp: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        transient: c_int,
        needmass: c_int,
    );
}

/// Grow the scratch buffers when an element needs more than the
/// first-mesh sizing (adaptive remesh across invocations). The original
/// uses per-call automatics (no overflow possible); growth only triggers
/// where the original would simply allocate bigger automatics, and the
/// leading sections used by the wrappers are unaffected, so valid inputs
/// never observe it.
fn ensure_state(state: &mut WaveState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.damp.resize(need * need, 0.0);
    state.mass.resize(need * need, 0.0);
    state.force.resize(need, 0.0);
    state.n = need;
}

/// Per-element scratch (`STIFF`/`DAMP`/`MASS`/`FORCE` automatics). `n`
/// is the first call's `MaxElementDOFs`-derived sizing.
struct WaveState {
    stiff: Vec<c_double>,
    damp: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    n: usize,
}

impl WaveState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            damp: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            n: 0,
        }
    }
}

static STATE: Mutex<WaveState> = Mutex::new(WaveState::new());

/// `wavesolver_init_` — `Time derivative order = 2` on the solver params.
///
/// # Safety
///
/// Called by the Elmer loader with the solver context current, exactly
/// like the original entry.
#[export_name = "wavesolver_init_"]
pub unsafe extern "C" fn wave_solver_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { wva_init() };
}

/// `wavesolver_` — the transient wave solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read for the Fortran `LOGICAL` flag.
#[export_name = "wavesolver_"]
pub unsafe extern "C" fn wave_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { wva_begin(solver) };

    unsafe { wva_info_start() };

    // SolverPars => GetSolverParams(); maxiter; dim; NeedMass.
    let maxiter = unsafe { wva_max_iter() };
    let _dim: c_int = unsafe { wva_coord_dim() };
    let needmass = unsafe { wva_need_mass() };
    let transient_flag: c_int = if transient.is_null() {
        0
    } else {
        unsafe { i32::from(*transient != 0) }
    };

    unsafe { wva_default_start() };

    for _ in 1..=maxiter.max(1) {
        // System assembly.
        unsafe { wva_default_initialize() };
        let mut init_handles: c_int = 1;
        let active = unsafe { wva_nof_active() };
        for t in 1..=active {
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { wva_begin_active(t, &mut n, &mut nd, &mut nb) };
            let ndnb = nd + nb;
            ensure_state(&mut state, (n as usize).max(ndnb as usize));
            let ld = state.n as c_int;
            unsafe { wva_lm_init(&mut init_handles) };
            local_matrix(&mut state, ndnb);
            unsafe {
                wva_lm_tail(
                    ndnb,
                    nb,
                    ld,
                    state.mass.as_mut_ptr(),
                    state.damp.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                    transient_flag,
                    needmass,
                );
            }
        }

        unsafe { wva_default_finish_bulk() };

        let mut init_handles_bc: c_int = 1;
        let nbnd = unsafe { wva_nof_boundary() };
        for t in 1..=nbnd {
            let (mut n, mut nd, mut active_el) = (0, 0, 0);
            unsafe { wva_begin_boundary(t, &mut n, &mut nd, &mut active_el) };
            if active_el == 0 {
                continue;
            }
            let (mut outflow, mut proceed) = (0, 0);
            unsafe { wva_lbc_head(&mut init_handles_bc, &mut outflow, &mut proceed) };
            if proceed == 0 {
                continue;
            }
            ensure_state(&mut state, nd as usize);
            let ld = state.n as c_int;
            local_matrix_bc(&mut state, nd, outflow);
            unsafe {
                wva_lbc_tail(
                    nd,
                    ld,
                    state.mass.as_mut_ptr(),
                    state.damp.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                    transient_flag,
                    needmass,
                );
            }
        }

        unsafe {
            wva_default_finish_bnd();
            wva_default_finish_assem();
            wva_default_dirichlet();
        }

        // And finally, solve.
        let _norm: c_double = unsafe { wva_default_solve() };

        if unsafe { wva_nonlin_converged() } > 0 {
            break;
        }

        unsafe { wva_avg_zero() };
    }

    unsafe { wva_default_finish() };
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `ndnb` is `nd+nb` (the internal
/// procedure's `nd` dummy).
fn local_matrix(state: &mut WaveState, ndnb: c_int) {
    let ld = state.n;
    // Zero the leading sections (the original's `MASS = 0` etc. on the
    // per-call automatics).
    state.stiff[..ld * ld].fill(0.0);
    state.damp[..ld * ld].fill(0.0);
    state.mass[..ld * ld].fill(0.0);
    state.force[..ld].fill(0.0);

    unsafe { wva_lm_nodes() };

    // Basis(ndnb) / dBasisdx(ndnb,3) automatics, sized by this element.
    let nn = ndnb as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];

    let nip = unsafe { wva_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { wva_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            wva_lm_eleminfo(
                u,
                v,
                w,
                ndnb,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;

        unsafe {
            wva_lm_ip(
                ndnb,
                state.n as c_int,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                s * detj,
                state.stiff.as_mut_ptr(),
                state.damp.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.force.as_mut_ptr(),
            );
        }
    }
}

/// The `LocalMatrixBC` internal procedure: boundary damping and flux by
/// boundary integration.
fn local_matrix_bc(state: &mut WaveState, nd: c_int, outflow: c_int) {
    let ld = state.n;
    // Zero the leading sections (the original's `STIFF = 0` etc.).
    state.stiff[..ld * ld].fill(0.0);
    state.damp[..ld * ld].fill(0.0);
    state.mass[..ld * ld].fill(0.0);
    state.force[..ld].fill(0.0);

    unsafe { wva_lbc_nodes() };

    let ndn = nd as usize;
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { wva_lbc_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { wva_lbc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            wva_lbc_eleminfo(
                u,
                v,
                w,
                nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;

        unsafe {
            wva_lbc_ip(
                nd,
                state.n as c_int,
                basis.as_mut_ptr(),
                s * detj,
                state.damp.as_mut_ptr(),
                state.force.as_mut_ptr(),
                outflow,
            );
        }
    }
}
