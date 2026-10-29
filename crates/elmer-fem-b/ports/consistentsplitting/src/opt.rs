//! `optimalsolutionupdate_` — port of OptimalSolutionUpdate.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE` replicated: `STIFF`/`Mass` (n×n, n = (dim+1)·maxdofs),
//!   `FORCE` (n), `LOAD` (p×4), `rho`/`mu` (p), `Velocity` ((dim+1)×p
//!   column-major), `Indexes` (p, Fortran-side), `Snew`/`R` (rows),
//!   `S`/`V` (rows×maxiter), `Round`/`CurrentDoneTime` (timestep tracking).
//! - Dead declarations not replicated: `ALocal`/`MLocal`/`PLocal` (only
//!   touched inside the dead `IF (.FALSE.)` outflow block),
//!   `TractionBCIndeces`/`NumberOfNormalTractionNodes` (no live use),
//!   `PrevSol`/`NonLinRes` (commented-out diagnostics only), `Norm`
//!   (initialized 0, never updated — this entry has no solve of its own).
//! - Dead values skipped with identical observables: `AK`/`rotterm`/
//!   `StabTerms` accumulation (read only in dead code), `GradDivParam`
//!   (passed but never used), `CPUTime` (pure timing), `atime`.
//! - `NonlinearIterationMethod` is overwritten to `'picard'` upstream, so
//!   Newton/Hybrid are always false — hardcoded with the same branch
//!   structure transcribed around the flags.
//! - `ConvectionStabilization` is forced `.FALSE.` upstream; hardcoded.
//! - The `IF (.FALSE.)` outflow/bubble-stabilization block is omitted.
//! - Steady-state with convection is Fatal, exactly as upstream — staged
//!   differentials must be transient.
//! - `PRINT *,'Residual norm for initial guess'` (Round==1) and the
//!   `WRITE(Message,...)`/`Info` residual reports are replicated in the
//!   wrappers (deterministic given identical inputs).
//! - Serial tree: `Parallel` is queried and transcribed; the parallel path
//!   (GCRUpdate et al.) is transcribed and compile-checked.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn cso_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn cso_num_active_elems() -> c_int;
    fn cso_nof_active() -> c_int;
    fn cso_dirichlet();
    fn cso_new_timestep(newts: *mut c_int);
    fn cso_round() -> c_int;
    fn cso_const_system(cb_in_use: *mut c_int);
    fn cso_dim() -> c_int;
    fn cso_block_precond(val: *mut c_int, found: *mut c_int);
    fn cso_max_iter() -> c_int;
    fn cso_mesh_max_dofs() -> c_int;
    fn cso_matrix_rows() -> c_int;
    fn cso_conv_flags(conv: *mut c_int, found: *mut c_int);
    fn cso_alloc(n: c_int);
    fn cso_fatal_steady();
    fn cso_fatal_velo();
    fn cso_fatal_alloc();
    fn cso_default_init();
    fn cso_bulk_restore();
    fn cso_progress_due(at0: c_double) -> c_int;
    fn cso_progress(nactive: c_int, t: c_int, at0: *mut c_double);
    fn cso_realtime() -> c_double;
    fn cso_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn cso_bodyforce(n: c_int, load: *mut c_double);
    fn cso_material(n: c_int, rho: *mut c_double, mu: *mut c_double);
    fn cso_gather_vel(nd: c_int, vel: *mut c_double);
    fn cso_lm_nodes();
    fn cso_lm_ip_begin() -> c_int;
    fn cso_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn cso_lm_eleminfo(
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
    fn cso_lm_point(
        n: c_int,
        nd: c_int,
        dim: c_int,
        stab: c_int,
        convect: c_int,
        newton: c_int,
        hybrid: c_int,
        sdet: c_double,
        load: *mut c_double,
        rho_n: *mut c_double,
        mu_n: *mut c_double,
        vel_n: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        mass: *mut c_double,
    );
    fn cso_default_time(m_ld: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn cso_update_eq(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn cso_finish_bulk();
    fn cso_finish();
    fn cso_parallel() -> c_int;
    fn cso_info_assembly_done();
    fn cso_set_bc(dof: c_int, ndofs: c_int);
    fn cso_bind_flow(ok: *mut c_int);
    fn cso_gather_snew(active: c_int, nrows: c_int, snew: *mut c_double);
    fn cso_orthogonalize(
        n: c_int,
        maxrounds: c_int,
        round: c_int,
        snew: *mut c_double,
        s: *mut c_double,
        v: *mut c_double,
        r: *mut c_double,
    );
    fn cso_parallel_update(
        n: c_int,
        maxrounds: c_int,
        round: c_int,
        snew: *mut c_double,
        s: *mut c_double,
        v: *mut c_double,
        r: *mut c_double,
    );
}

struct State {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    mass: Vec<c_double>,
    load: Vec<c_double>,
    rho: Vec<c_double>,
    mu: Vec<c_double>,
    vel: Vec<c_double>,
    snew: Vec<c_double>,
    r: Vec<c_double>,
    s: Vec<c_double>,
    v: Vec<c_double>,
    p: usize,
    n: usize,
    dim: usize,
    maxiter: usize,
    allocations_done: bool,
}

impl State {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            mass: Vec::new(),
            load: Vec::new(),
            rho: Vec::new(),
            mu: Vec::new(),
            vel: Vec::new(),
            snew: Vec::new(),
            r: Vec::new(),
            s: Vec::new(),
            v: Vec::new(),
            p: 0,
            n: 0,
            dim: 0,
            maxiter: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<State> = Mutex::new(State::new());

fn ensure(state: &mut State, p: usize, dim: usize, rows: usize, maxiter: usize) {
    if p <= state.p && dim == state.dim && rows * maxiter <= state.s.len() {
        return;
    }
    let pp = p.max(state.p);
    let n = (dim + 1) * pp;
    state.stiff.resize(n * n, 0.0);
    state.force.resize(n, 0.0);
    state.mass.resize(n * n, 0.0);
    state.load.resize(pp * 4, 0.0);
    state.rho.resize(pp, 0.0);
    state.mu.resize(pp, 0.0);
    state.vel.resize((dim + 1) * pp, 0.0);
    state.snew.resize(rows, 0.0);
    state.r.resize(rows, 0.0);
    state.s.resize(rows * maxiter, 0.0);
    state.v.resize(rows * maxiter, 0.0);
    state.p = pp;
    state.n = n;
    state.dim = dim;
    state.maxiter = maxiter;
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { cso_begin(model, solver) };
    let transient_flag: bool = unsafe { *transient } != 0;

    let mut newts = 0;
    unsafe { cso_new_timestep(&mut newts) };
    let round = unsafe { cso_round() };

    let mut cb_in_use = 0;
    unsafe { cso_const_system(&mut cb_in_use) };

    let dim = unsafe { cso_dim() } as usize;
    // Newton/Hybrid/Picard: NonlinearIterationMethod is overwritten to
    // 'picard' upstream — hardcoded picard with the flags transcribed.
    let (newton, hybrid) = (0, 0);

    let (mut block, mut block_found) = (0, 0);
    unsafe { cso_block_precond(&mut block, &mut block_found) };
    let _ = block_found;

    let maxiter = unsafe { cso_max_iter() } as usize;
    let p0 = unsafe { cso_mesh_max_dofs() } as usize;
    let rows = unsafe { cso_matrix_rows() } as usize;

    if !state.allocations_done {
        let n = (dim + 1) * p0;
        state.p = p0;
        state.n = n;
        state.dim = dim;
        state.maxiter = maxiter;
        if state.stiff.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.load.try_reserve(p0 * 4).is_err()
            || state.rho.try_reserve(p0).is_err()
            || state.mu.try_reserve(p0).is_err()
            || state.vel.try_reserve((dim + 1) * p0).is_err()
            || state.snew.try_reserve(rows).is_err()
            || state.r.try_reserve(rows).is_err()
            || state.s.try_reserve(rows * maxiter).is_err()
            || state.v.try_reserve(rows * maxiter).is_err()
        {
            unsafe { cso_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.load.resize(p0 * 4, 0.0);
        state.rho.resize(p0, 0.0);
        state.mu.resize(p0, 0.0);
        state.vel.resize((dim + 1) * p0, 0.0);
        state.snew.resize(rows, 0.0);
        state.r.resize(rows, 0.0);
        state.s.resize(rows * maxiter, 0.0);
        state.v.resize(rows * maxiter, 0.0);
        unsafe { cso_alloc(p0 as c_int) };
        state.allocations_done = true;
    }

    let (mut convect, mut conv_found) = (0, 0);
    unsafe { cso_conv_flags(&mut convect, &mut conv_found) };
    if conv_found == 0 {
        convect = 1;
    }
    let _ = conv_found;
    // ConvectionStabilization is forced .FALSE. upstream.
    let stab = 0;

    let mut at0 = unsafe { cso_realtime() };

    if cb_in_use != 0 {
        unsafe { cso_bulk_restore() };
    } else {
        unsafe { cso_default_init() };
        let active_src = unsafe { cso_num_active_elems() };
        for t in 1..=active_src {
            if unsafe { cso_progress_due(at0) } != 0 {
                unsafe { cso_progress(active_src, t, &mut at0) };
            }
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { cso_begin_active(t, &mut n, &mut nd, &mut nb) };
            let _ = nb;
            ensure(&mut state, (n as usize).max(nd as usize), dim, rows, maxiter);
            let nmax = state.n;

            unsafe { cso_bodyforce(n, state.load.as_mut_ptr()) };
            unsafe { cso_material(n, state.rho.as_mut_ptr(), state.mu.as_mut_ptr()) };

            for v in state.vel.iter_mut() {
                *v = 0.0;
            }
            if convect == 0 {
                // Velocity stays zeroed.
            } else if transient_flag {
                unsafe { cso_gather_vel(nd, state.vel.as_mut_ptr()) };
            } else {
                unsafe { cso_fatal_steady() };
                std::process::abort();
            }

            for v in state.stiff.iter_mut() {
                *v = 0.0;
            }
            for v in state.force.iter_mut() {
                *v = 0.0;
            }
            for v in state.mass.iter_mut() {
                *v = 0.0;
            }
            unsafe { cso_lm_nodes() };
            let nip = unsafe { cso_lm_ip_begin() };
            let mut basis = vec![0.0 as c_double; nd as usize];
            let mut dbasis = vec![0.0 as c_double; nd as usize * 3];
            for k in 1..=nip {
                let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
                unsafe { cso_lm_ip_point(k, &mut u, &mut v, &mut w, &mut s) };
                let (mut detj, mut stat) = (0.0, 0);
                unsafe {
                    cso_lm_eleminfo(
                        u,
                        v,
                        w,
                        nd,
                        basis.as_mut_ptr(),
                        dbasis.as_mut_ptr(),
                        &mut detj,
                        &mut stat,
                    )
                };
                let _ = stat;
                unsafe {
                    cso_lm_point(
                        n,
                        nd,
                        dim as c_int,
                        stab,
                        convect,
                        newton,
                        hybrid,
                        s * detj,
                        state.load.as_mut_ptr(),
                        state.rho.as_mut_ptr(),
                        state.mu.as_mut_ptr(),
                        state.vel.as_mut_ptr(),
                        basis.as_mut_ptr(),
                        dbasis.as_mut_ptr(),
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                        state.mass.as_mut_ptr(),
                    )
                };
            }

            if transient_flag {
                unsafe {
                    cso_default_time(
                        nmax as c_int,
                        state.mass.as_mut_ptr(),
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                    )
                };
            }
            unsafe {
                cso_update_eq(
                    nmax as c_int,
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }

        unsafe { cso_finish_bulk() };
    }

    unsafe { cso_finish() };

    let parallel = unsafe { cso_parallel() };
    if parallel != 0 {
        unsafe { cso_dirichlet() };
    } else {
        unsafe { cso_set_bc(1, dim as c_int + 1) };
        unsafe { cso_set_bc(2, dim as c_int + 1) };
        unsafe { cso_set_bc(3, dim as c_int + 1) };
        if dim > 2 {
            unsafe { cso_set_bc(4, dim as c_int + 1) };
        }
    }

    unsafe { cso_info_assembly_done() };

    for v in state.snew.iter_mut() {
        *v = 0.0;
    }
    let mut ok_flow = 0;
    unsafe { cso_bind_flow(&mut ok_flow) };
    if ok_flow == 0 {
        unsafe { cso_fatal_velo() };
        std::process::abort();
    }
    let active = unsafe { cso_nof_active() };
    unsafe { cso_gather_snew(active, rows as c_int, state.snew.as_mut_ptr()) };

    if parallel != 0 {
        unsafe {
            cso_parallel_update(
                rows as c_int,
                maxiter as c_int,
                round,
                state.snew.as_mut_ptr(),
                state.s.as_mut_ptr(),
                state.v.as_mut_ptr(),
                state.r.as_mut_ptr(),
            )
        };
    } else {
        unsafe {
            cso_orthogonalize(
                rows as c_int,
                maxiter as c_int,
                round,
                state.snew.as_mut_ptr(),
                state.s.as_mut_ptr(),
                state.v.as_mut_ptr(),
                state.r.as_mut_ptr(),
            )
        };
    }
}

