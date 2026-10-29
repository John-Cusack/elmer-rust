//! Port of `fem/src/modules/Smitc.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/Smitc.F90` — two F77-global entries
//! `smitcsolver_` and `smitcsolver_init0_` with three `CONTAINS`
//! internals (`LocalMatrix`, `PerforatedElasticity`,
//! `CovariantInterpolation`). The emitted dynamic symbols are exactly
//! `smitcsolver_` and `smitcsolver_init0_` (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/smitc_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. The whole
//! `LocalMatrix` internal procedure (with `PerforatedElasticity` and
//! `CovariantInterpolation` inlined, including the `MATMUL` reductions
//! and the `SUM` load/material interpolations) stays in the Fortran
//! wrapper, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors Smitc.F90 statement for statement):
//! - `SAVE` semantics replicated: `STIFF`/`MASS`/`DAMP`/`FORCE` and the
//!   material/load vectors are allocated once at first-mesh sizing
//!   (`MaxElementDOFs`, `Variable % DOFs`), never resized (bug-
//!   compatible: a larger second mesh is out of contract, exactly as in
//!   the original, whose fixed SAVE arrays would overflow); the
//!   internal-procedure `SAVE ElementNodes` and IP state plus the
//!   `Indexes` SAVE pointer live as wrapper-side SAVE state with
//!   identical lifetimes.
//! - `LOGICAL :: AllocationsDone = .FALSE.` latches process-wide via the
//!   `allocations_done` flag, never reset per call (SAVE caution).
//! - `n = GetElementDOFs(Indexes)` is evaluated-then-discarded exactly as
//!   in the original (`GetElementNOFNodes` overwrites `n`); the `Indexes`
//!   fill itself is write-only, kept wrapper-side.
//! - The `Model % Materials(mat_id)` read in the BodyForce-index branch
//!   (not `BodyForces`) is transcribed verbatim.
//! - Allocation failure calls the same `Fatal('SmitcSolver', ...)` via
//!   `try_reserve` (the `STAT=istat` check), never abort directly; the
//!   gravity-missing path calls the same `Fatal`.
//! - `Model`, `dt` are unused by the body, exactly as in the original
//!   (received and ignored); `Norm = DefaultSolve()` is kept and
//!   discarded.
//! - The dead `FMIN4 IF(.FALSE.)` block is omitted (never executes).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sma_init0();
    fn sma_begin(h_solver: *const c_void, h_model: *const c_void);
    fn sma_dofs() -> c_int;
    fn sma_max_elem_dofs() -> c_int;
    fn sma_nof_active() -> c_int;
    fn sma_nonlin_converged() -> c_int;
    fn sma_alloc_indexes(nmax: c_int);
    fn sma_relax_alpha(alpha: *mut c_double, found: *mut c_int);
    fn sma_harmonic_mode() -> c_int;
    fn sma_kernel_version() -> c_int;
    fn sma_use_gravity() -> c_int;
    fn sma_gravity_value(grav: *mut c_double, found: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn sma_local_matrix(
        ld: c_int,
        nn: c_int,
        n: c_int,
        dofs: c_int,
        stiff: *mut c_double,
        damp: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        load: *mut c_double,
        poisson: *mut c_double,
        young: *mut c_double,
        thickness: *mut c_double,
        tension: *mut c_double,
        dampcoef: *mut c_double,
        springcoef: *mut c_double,
        holefraction: *mut c_double,
        holesize: *mut c_double,
        density: *mut c_double,
        relax: c_double,
        kernel: c_int,
        holecorr: c_int,
        gotholetype: c_int,
        grav: c_double,
    );
    fn sma_mat_index(idx: *mut c_int, got: *mut c_int);
    fn sma_bf_index(idx: *mut c_int, got: *mut c_int);
    fn sma_material_select(got: c_int, mat_id: c_int);
    fn sma_bodyforce_select(got: c_int, mat_id: c_int);
    fn sma_hole_correction() -> c_int;
    fn sma_max_iter(maxiter: *mut c_int, got: *mut c_int);
    fn sma_info_banner();
    fn sma_default_start();
    fn sma_default_initialize();
    fn sma_cpu() -> c_double;
    fn sma_info_assembly(at: c_double);
    fn sma_info_solve(st: c_double);
    fn sma_begin_active(t: c_int, n_dofs: *mut c_int, n_nodes: *mut c_int);
    fn sma_bodyforce() -> c_int;
    fn sma_material() -> c_int;
    fn sma_load_pressure(nn: c_int, n: c_int, load: *mut c_double, load2: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sma_material_vecs(
        nn: c_int,
        n: c_int,
        poisson: *mut c_double,
        density: *mut c_double,
        thickness: *mut c_double,
        young: *mut c_double,
        tension: *mut c_double,
    );
    fn sma_dampcoef(nn: c_int, n: c_int, dampcoef: *mut c_double);
    fn sma_springcoef(nn: c_int, n: c_int, springcoef: *mut c_double);
    fn sma_hole(nn: c_int, n: c_int, holesize: *mut c_double, holefraction: *mut c_double)
        -> c_int;
    fn sma_default_2nd_order(
        ld: c_int,
        mass: *mut c_double,
        damp: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn sma_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn sma_default_update_mass(ld: c_int, m: *mut c_double);
    fn sma_default_update_damp(ld: c_int, b: *mut c_double);
    fn sma_default_finish_bulk();
    fn sma_default_finish_boundary();
    fn sma_default_finish_assembly();
    fn sma_default_dirichlet();
    fn sma_default_solve() -> c_double;
    fn sma_default_finish();
    fn sma_fatal_alloc();
    fn sma_fatal_gravity();
}

/// First-mesh SAVE state: `STIFF`/`MASS`/`DAMP` (`ld` x `ld`), `FORCE`
/// (`ld`), the ten `N` material/load vectors, plus the `AllocationsDone`
/// latch. `N` is the first call's `MaxElementDOFs`, `ld = DOFs * N`.
struct SmitcState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    damp: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    load2: Vec<c_double>,
    poisson: Vec<c_double>,
    density: Vec<c_double>,
    thickness: Vec<c_double>,
    young: Vec<c_double>,
    tension: Vec<c_double>,
    dampcoef: Vec<c_double>,
    springcoef: Vec<c_double>,
    holefraction: Vec<c_double>,
    holesize: Vec<c_double>,
    n: usize,
    dofs: usize,
    allocations_done: bool,
}

impl SmitcState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            damp: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            load2: Vec::new(),
            poisson: Vec::new(),
            density: Vec::new(),
            thickness: Vec::new(),
            young: Vec::new(),
            tension: Vec::new(),
            dampcoef: Vec::new(),
            springcoef: Vec::new(),
            holefraction: Vec::new(),
            holesize: Vec::new(),
            n: 0,
            dofs: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<SmitcState> = Mutex::new(SmitcState::new());

/// `smitcsolver_init0_` — the init entry: declares the `-dofs 3
/// Deflection` variable and the plate-solver list flags.
///
/// # Safety
///
/// All arguments are unused by the body, exactly as in the original.
#[export_name = "smitcsolver_init0_"]
pub unsafe extern "C" fn smitc_solver_init0(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sma_init0() };
}

/// `smitcsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `transient` must be a valid read if non-null.
#[export_name = "smitcsolver_"]
pub unsafe extern "C" fn smitc_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sma_begin(solver, model) };

    // Allocate some permanent storage, this is done first time only.
    if !state.allocations_done {
        let dofs = unsafe { sma_dofs() } as usize;
        let n = unsafe { sma_max_elem_dofs() } as usize;
        let ld = dofs * n;
        state.dofs = dofs;
        state.n = n;
        let ok = state.stiff.try_reserve(ld * ld).is_ok()
            && state.mass.try_reserve(ld * ld).is_ok()
            && state.damp.try_reserve(ld * ld).is_ok()
            && state.force.try_reserve(ld).is_ok()
            && state.load.try_reserve(n).is_ok()
            && state.load2.try_reserve(n).is_ok()
            && state.poisson.try_reserve(n).is_ok()
            && state.density.try_reserve(n).is_ok()
            && state.thickness.try_reserve(n).is_ok()
            && state.young.try_reserve(n).is_ok()
            && state.tension.try_reserve(n).is_ok()
            && state.dampcoef.try_reserve(n).is_ok()
            && state.springcoef.try_reserve(n).is_ok()
            && state.holefraction.try_reserve(n).is_ok()
            && state.holesize.try_reserve(n).is_ok();
        if !ok {
            unsafe { sma_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(ld * ld, 0.0);
        state.mass.resize(ld * ld, 0.0);
        state.damp.resize(ld * ld, 0.0);
        state.force.resize(ld, 0.0);
        state.load.resize(n, 0.0);
        state.load2.resize(n, 0.0);
        state.poisson.resize(n, 0.0);
        state.density.resize(n, 0.0);
        state.thickness.resize(n, 0.0);
        state.young.resize(n, 0.0);
        state.tension.resize(n, 0.0);
        state.dampcoef.resize(n, 0.0);
        state.springcoef.resize(n, 0.0);
        state.holefraction.resize(n, 0.0);
        state.holesize.resize(n, 0.0);
        unsafe { sma_alloc_indexes(n as c_int) };
        state.allocations_done = true;
    }

    let dofs = state.dofs as c_int;
    let nmax = state.n as c_int;
    let ld = state.dofs * state.n;
    let ld_c = ld as c_int;

    // Do some additional initialization, and go for it.
    unsafe { sma_info_banner() };

    let mut relax = 0.0;
    let mut found = 0;
    unsafe { sma_relax_alpha(&mut relax, &mut found) };
    if found == 0 {
        relax = 1.0;
    }

    let eigen = unsafe { sma_harmonic_mode() };
    let kernel = unsafe { sma_kernel_version() };

    let mut grav = 0.0;
    if unsafe { sma_use_gravity() } != 0 {
        let mut gfound = 0;
        unsafe { sma_gravity_value(&mut grav, &mut gfound) };
        if gfound == 0 {
            unsafe { sma_fatal_gravity() };
            std::process::abort();
        }
    }

    unsafe { sma_default_start() };

    let mut maxiter = 0;
    let mut gotit = 0;
    unsafe { sma_max_iter(&mut maxiter, &mut gotit) };
    if gotit == 0 {
        maxiter = 1;
    }

    let transient_sim = unsafe { *transient } != 0;

    for _iter in 1..=maxiter {
        let at = unsafe { sma_cpu() };
        unsafe { sma_default_initialize() };

        let mut mat_id = 0;
        let mut got_mat_id = 0;
        unsafe { sma_mat_index(&mut mat_id, &mut got_mat_id) };
        unsafe { sma_material_select(got_mat_id, mat_id) };

        let mut bf_id = 0;
        let mut got_bf_id = 0;
        unsafe { sma_bf_index(&mut bf_id, &mut got_bf_id) };
        unsafe { sma_bodyforce_select(got_bf_id, mat_id) };

        let holecorr = unsafe { sma_hole_correction() };

        // Do the assembly.
        let active = unsafe { sma_nof_active() };
        for t in 1..=active {
            let (mut _ndofs_count, mut n) = (0, 0);
            unsafe { sma_begin_active(t, &mut _ndofs_count, &mut n) };

            if got_bf_id == 0 {
                unsafe { sma_bodyforce() };
            }

            state.load.fill(0.0);
            state.load2.fill(0.0);
            unsafe {
                sma_load_pressure(
                    nmax,
                    n,
                    state.load.as_mut_ptr(),
                    state.load2.as_mut_ptr(),
                )
            };

            if got_mat_id == 0 {
                unsafe { sma_material() };
            }

            unsafe {
                sma_material_vecs(
                    nmax,
                    n,
                    state.poisson.as_mut_ptr(),
                    state.density.as_mut_ptr(),
                    state.thickness.as_mut_ptr(),
                    state.young.as_mut_ptr(),
                    state.tension.as_mut_ptr(),
                );
                sma_dampcoef(nmax, n, state.dampcoef.as_mut_ptr());
                sma_springcoef(nmax, n, state.springcoef.as_mut_ptr());
            }

            let mut got_holetype = 0;
            if holecorr != 0 {
                got_holetype = unsafe {
                    sma_hole(
                        nmax,
                        n,
                        state.holesize.as_mut_ptr(),
                        state.holefraction.as_mut_ptr(),
                    )
                };
            }

            // Get element local matrix, and rhs vector.
            state.stiff.fill(0.0);
            state.mass.fill(0.0);
            state.damp.fill(0.0);
            state.force.fill(0.0);
            unsafe {
                sma_local_matrix(
                    ld_c,
                    nmax,
                    n,
                    dofs,
                    state.stiff.as_mut_ptr(),
                    state.damp.as_mut_ptr(),
                    state.mass.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                    state.load.as_mut_ptr(),
                    state.poisson.as_mut_ptr(),
                    state.young.as_mut_ptr(),
                    state.thickness.as_mut_ptr(),
                    state.tension.as_mut_ptr(),
                    state.dampcoef.as_mut_ptr(),
                    state.springcoef.as_mut_ptr(),
                    state.holefraction.as_mut_ptr(),
                    state.holesize.as_mut_ptr(),
                    state.density.as_mut_ptr(),
                    relax,
                    kernel,
                    holecorr,
                    got_holetype,
                    grav,
                )
            };

            if transient_sim {
                unsafe {
                    sma_default_2nd_order(
                        ld_c,
                        state.mass.as_mut_ptr(),
                        state.damp.as_mut_ptr(),
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                    )
                };
            }

            // Update global matrix and rhs vector from local matrix & vector.
            unsafe {
                sma_default_update(
                    ld_c,
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };

            if eigen != 0 {
                unsafe {
                    sma_default_update_mass(ld_c, state.mass.as_mut_ptr());
                    sma_default_update_damp(ld_c, state.damp.as_mut_ptr());
                }
            }
        }
        unsafe {
            sma_default_finish_bulk();
            sma_default_finish_boundary();
            sma_default_finish_assembly();
        }

        // Dirichlet boundary conditions.
        unsafe { sma_default_dirichlet() };

        let at_now = unsafe { sma_cpu() };
        unsafe { sma_info_assembly(at_now - at) };

        // Solve the system and we are done.
        let st = unsafe { sma_cpu() };
        let _norm: c_double = unsafe { sma_default_solve() };
        let st_now = unsafe { sma_cpu() };
        unsafe { sma_info_solve(st_now - st) };

        if unsafe { sma_nonlin_converged() } == 1 {
            break;
        }
    }

    unsafe { sma_default_finish() };
}
