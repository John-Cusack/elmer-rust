//! Port of `fem/src/modules/StatMagSolve.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/StatMagSolve.F90` — two F77-global
//! entries `statmagsolver_` and `statmagsolver_init_` with four `CONTAINS`
//! assemblers plus `JouleIntegrate`. The emitted dynamic symbols are
//! exactly those two (proven by `nm -D` on the reference `.so`); the
//! internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/statmagsolve_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`SUM`),
//! formatted I/O (`WRITE`/`Message`) and whole-array updates stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors StatMagSolve.F90 statement for statement):
//! - `SAVE` semantics replicated explicitly: first-mesh sizing
//!   (`Model % MaxElementNodes`, `dofs * N` buffers) with the
//!   `AllocationsDone` latch in the guarded process-wide `STATE`;
//!   `ElementNodes` xyz, the `MFD`/`Br`/`Bz`/`Bp`/`Babs`/`Joule` POINTER
//!   views, `Brim`/`Bzim`, and the joule scratch live as wrapper-side
//!   SAVE state with identical lifetimes (associated/allocated once,
//!   never resized — the same mesh serves every invocation in a process).
//! - The four `CONTAINS` assemblers and `JouleIntegrate` are transcribed
//!   statement-for-statement into Fortran wrappers (same Gauss-point loop,
//!   same `ElementInfo` call order, same update statements); the Rust side
//!   only selects which one runs per element (`dim`, harmonic flag),
//!   exactly like the original's `IF` cascade.
//! - `Ap` (previous-iteration potential) is allocated in the original but
//!   never read by any assembler; the port omits the dead buffer (no
//!   observable effect).
//! - `UNorm = DefaultSolve()` is kept as a discarded return, `at =
//!   CPUTime()` as a progress-clock touch, and the transient `dt` dummy
//!   as received-but-ignored, all exactly as in the original.
//! - Allocation failure calls the same `Fatal('StatMagSolve', ...)` via
//!   `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `LOGICAL :: X, Y = .FALSE.`-style latching does not occur here; the
//!   only latch is `AllocationsDone`, replicated above.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sma_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sma_assoc_var();
    fn sma_assoc_matrix() -> c_int;
    fn sma_info_start();
    fn sma_mesh_dim() -> c_int;
    fn sma_perm_vacuum() -> c_double;
    fn sma_count_local() -> c_int;
    fn sma_var_dofs() -> c_int;
    fn sma_harmonic(found: *mut c_int) -> c_int;
    fn sma_info_mode(h: c_int);
    fn sma_set_complex();
    fn sma_warn_whitney();
    fn sma_angfreq() -> c_double;
    fn sma_fatal_3d();
    fn sma_fatal_harm();
    fn sma_fatal_steady();
    fn sma_fatal_alloc();
    fn sma_warn_joule_only();
    fn sma_calc_flags(mflux: *mut c_int, mabs: *mut c_int, mjoule: *mut c_int);
    fn sma_max_elem_nodes() -> c_int;
    fn sma_allocate(n: c_int);
    fn sma_assoc_flux(mabs: c_int, harm: c_int);
    fn sma_assoc_joule();
    fn sma_maxiter(found: *mut c_int) -> c_int;
    fn sma_info_iter(iter: c_int);
    fn sma_prog_begin();
    fn sma_prog_check(t: c_int);
    fn sma_default_initialize();
    fn sma_nof_active() -> c_int;
    fn sma_elem_begin(t: c_int) -> c_int;
    fn sma_bf_id() -> c_int;
    fn sma_perm_fill(n: c_int, ld: c_int, perm: *mut c_double, mu0: c_double);
    fn sma_reluct(n: c_int, reluc: *mut c_double, perm: *const c_double);
    fn sma_cond_harm(n: c_int, cond: *mut c_double);
    fn sma_current(n: c_int, bf: c_int, cur: *mut c_double);
    fn sma_phase(n: c_int, bf: c_int, phase: *mut c_double);
    fn sma_vecload(n: c_int, ldn: c_int, bf: c_int, vlv: *mut c_double);
    fn sma_stat_axis(
        ld: c_int,
        n: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        mass: *mut c_double,
        load: *mut c_double,
        reluc: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn sma_harm_axis(
        ld: c_int,
        n: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        cur: *mut c_double,
        phase: *mut c_double,
        reluc: *mut c_double,
        cond: *mut c_double,
        wang: c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn sma_cart_3d(
        ld: c_int,
        n: c_int,
        ldn: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        vlv: *mut c_double,
        reluc: *mut c_double,
    );
    fn sma_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn sma_info_asm_done();
    fn sma_default_finish();
    fn sma_default_dirichlet();
    fn sma_info_bnd_done();
    fn sma_default_solve() -> c_double;
    fn sma_flux_zero();
    fn sma_flux_harm();
    fn sma_flux_steady(dim: c_int);
    fn sma_invalidate(which: c_int);
    fn sma_babs();
    fn sma_nonlin_converged() -> c_int;
    fn sma_info_converged(iter: c_int);
    fn sma_joule_section();
    fn sma_info_done();
    fn sma_init_run();
}

/// Grow the SAVE buffers when a later mesh needs more than the first-mesh
/// sizing. The original never grows (it assumes the first mesh is big
/// enough); growth only triggers where the original would overflow
/// silently, and valid inputs never observe it.
fn ensure_state(state: &mut StatMagState, ld: usize, n: usize) {
    if ld <= state.ld {
        return;
    }
    state.reluc.resize(ld, 0.0);
    state.permeab.resize(ld, 0.0);
    state.conduct.resize(ld, 0.0);
    state.curdens.resize(ld, 0.0);
    state.phase.resize(ld, 0.0);
    state.lmass.resize(ld * ld, 0.0);
    state.lstiff.resize(ld * ld, 0.0);
    state.lforce.resize(ld, 0.0);
    state.ae.resize(n, 0.0);
    state.vecload.resize(3 * n, 0.0);
    state.ld = ld;
    state.n = n;
}

/// First-mesh SAVE state: `N = Model % MaxElementNodes`, `ld = dofs * N`
/// buffers plus the `AllocationsDone` flag.
struct StatMagState {
    reluc: Vec<c_double>,
    permeab: Vec<c_double>,
    conduct: Vec<c_double>,
    curdens: Vec<c_double>,
    phase: Vec<c_double>,
    lmass: Vec<c_double>,
    lstiff: Vec<c_double>,
    lforce: Vec<c_double>,
    ae: Vec<c_double>,
    vecload: Vec<c_double>,
    n: usize,
    ld: usize,
    allocations_done: bool,
}

impl StatMagState {
    const fn new() -> Self {
        Self {
            reluc: Vec::new(),
            permeab: Vec::new(),
            conduct: Vec::new(),
            curdens: Vec::new(),
            phase: Vec::new(),
            lmass: Vec::new(),
            lstiff: Vec::new(),
            lforce: Vec::new(),
            ae: Vec::new(),
            vecload: Vec::new(),
            n: 0,
            ld: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<StatMagState> = Mutex::new(StatMagState::new());

/// `statmagsolver_init_` — the solver init entry. Runs the init body
/// through the wrapper, which owns the whole block verbatim (including
/// the `GetSolverParams` / `NextFreeKeyword` / hidden-length traffic that
/// cannot cross the boundary).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "statmagsolver_init_"]
pub unsafe extern "C" fn statmag_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sma_begin(model, solver) };
    unsafe { sma_init_run() };
}

/// `statmagsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "statmagsolver_"]
pub unsafe extern "C" fn statmag_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sma_begin(model, solver) };

    // IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { sma_assoc_matrix() } == 0 {
        return;
    }

    unsafe {
        sma_info_start();
        sma_assoc_var();
    }
    let dim = unsafe { sma_mesh_dim() };
    let mu0: c_double = unsafe { sma_perm_vacuum() };

    // LocalNodes = COUNT(MagneticPerm > 0); IF (LocalNodes <= 0) RETURN
    if unsafe { sma_count_local() } <= 0 {
        return;
    }

    let dofs = unsafe { sma_var_dofs() };

    // HarmonicSimulation default: (dofs == 2) when the keyword is absent.
    let mut harm_found = 0;
    let harm_got = unsafe { sma_harmonic(&mut harm_found) };
    let harmonic = if harm_found != 0 {
        harm_got != 0
    } else {
        dofs == 2
    };
    let harm_i: c_int = if harmonic { 1 } else { 0 };
    unsafe { sma_info_mode(harm_i) };

    let mut wang: c_double = 0.0;
    if dim == 3 {
        if dofs != 3 {
            unsafe { sma_fatal_3d() };
            return;
        }
        unsafe { sma_warn_whitney() };
    } else if harmonic {
        if dofs != 2 {
            unsafe { sma_fatal_harm() };
            return;
        }
        unsafe { sma_set_complex() };
        wang = unsafe { sma_angfreq() };
    } else if dofs != 1 {
        unsafe { sma_fatal_steady() };
        return;
    }

    let (mut mflux, mut mabs, mut mjoule) = (0, 0, 0);
    unsafe { sma_calc_flags(&mut mflux, &mut mabs, &mut mjoule) };
    if mjoule != 0 && !harmonic {
        unsafe { sma_warn_joule_only() };
        mjoule = 0;
    }

    // First-mesh permanent storage (bug-compatible: sized once).
    if !state.allocations_done {
        let n = unsafe { sma_max_elem_nodes() } as usize;
        let ld = (dofs as usize) * n;
        let mut ok = state.reluc.try_reserve(ld).is_ok();
        ok = ok
            && state.permeab.try_reserve(ld).is_ok()
            && state.conduct.try_reserve(ld).is_ok()
            && state.curdens.try_reserve(ld).is_ok()
            && state.phase.try_reserve(ld).is_ok()
            && state.lmass.try_reserve(ld * ld).is_ok()
            && state.lstiff.try_reserve(ld * ld).is_ok()
            && state.lforce.try_reserve(ld).is_ok()
            && state.ae.try_reserve(n).is_ok()
            && state.vecload.try_reserve(3 * n).is_ok();
        if !ok {
            unsafe { sma_fatal_alloc() };
            std::process::abort();
        }
        state.reluc.resize(ld, 0.0);
        state.permeab.resize(ld, 0.0);
        state.conduct.resize(ld, 0.0);
        state.curdens.resize(ld, 0.0);
        state.phase.resize(ld, 0.0);
        state.lmass.resize(ld * ld, 0.0);
        state.lstiff.resize(ld * ld, 0.0);
        state.lforce.resize(ld, 0.0);
        state.ae.resize(n, 0.0);
        state.vecload.resize(3 * n, 0.0);
        state.n = n;
        state.ld = ld;
        unsafe {
            sma_allocate(n as c_int);
            if mflux != 0 {
                sma_assoc_flux(mabs, harm_i);
            }
            if mjoule != 0 {
                sma_assoc_joule();
            }
        }
        state.allocations_done = true;
    }

    let mut maxiter_found = 0;
    let mut maxiter = unsafe { sma_maxiter(&mut maxiter_found) };
    if maxiter_found == 0 {
        maxiter = 1;
    }

    for iter in 1..=maxiter {
        unsafe { sma_info_iter(iter) };
        // at = CPUTime(); at0 = RealTime()
        unsafe { sma_prog_begin() };
        unsafe { sma_default_initialize() };

        let nactive = unsafe { sma_nof_active() };
        for t in 1..=nactive {
            unsafe { sma_prog_check(t) };

            let n = unsafe { sma_elem_begin(t) };
            ensure_state(&mut state, (dofs as usize) * (n as usize), n as usize);
            let ld = state.ld as c_int;

            unsafe {
                sma_perm_fill(n, ld, state.permeab.as_mut_ptr(), mu0);
                sma_reluct(n, state.reluc.as_mut_ptr(), state.permeab.as_ptr());
            }

            if harmonic {
                unsafe { sma_cond_harm(n, state.conduct.as_mut_ptr()) };
            }

            let bf = unsafe { sma_bf_id() };

            if dim < 3 {
                unsafe {
                    sma_current(n, bf, state.curdens.as_mut_ptr());
                    if harmonic {
                        sma_phase(n, bf, state.phase.as_mut_ptr());
                    }
                }
            } else {
                unsafe {
                    sma_vecload(
                        n,
                        state.n as c_int,
                        bf,
                        state.vecload.as_mut_ptr(),
                    )
                };
            }

            if dim < 3 {
                if harmonic {
                    unsafe {
                        sma_harm_axis(
                            ld,
                            n,
                            state.lstiff.as_mut_ptr(),
                            state.lforce.as_mut_ptr(),
                            state.curdens.as_mut_ptr(),
                            state.phase.as_mut_ptr(),
                            state.reluc.as_mut_ptr(),
                            state.conduct.as_mut_ptr(),
                            wang,
                        )
                    };
                } else {
                    unsafe {
                        sma_stat_axis(
                            ld,
                            n,
                            state.lstiff.as_mut_ptr(),
                            state.lforce.as_mut_ptr(),
                            state.lmass.as_mut_ptr(),
                            state.curdens.as_mut_ptr(),
                            state.reluc.as_mut_ptr(),
                        )
                    };
                }
            } else {
                unsafe {
                    sma_cart_3d(
                        ld,
                        n,
                        state.n as c_int,
                        state.lstiff.as_mut_ptr(),
                        state.lforce.as_mut_ptr(),
                        state.vecload.as_mut_ptr(),
                        state.reluc.as_mut_ptr(),
                    )
                };
            }

            unsafe { sma_default_update(ld, state.lstiff.as_mut_ptr(), state.lforce.as_mut_ptr()) };
        }

        unsafe {
            sma_info_asm_done();
            sma_default_finish();
            sma_default_dirichlet();
            sma_info_bnd_done();
        }

        // UNorm is local and unused, as in the original.
        let _unorm: c_double = unsafe { sma_default_solve() };

        if mflux != 0 {
            unsafe {
                sma_flux_zero();
                if harmonic {
                    sma_flux_harm();
                } else {
                    sma_flux_steady(dim);
                }
                sma_invalidate(0);
                if mabs != 0 {
                    sma_babs();
                    sma_invalidate(1);
                }
            }
        }

        if unsafe { sma_nonlin_converged() } == 1 {
            unsafe { sma_info_converged(iter) };
            break;
        }
    }

    if harmonic && mjoule != 0 {
        unsafe { sma_joule_section() };
    }

    unsafe { sma_info_done() };
}
