//! `velocitysolver_` — port of VelocitySolver.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE` replicated: `STIFF`/`MASS` (n×n, n = dim·maxdofs),
//!   `FORCE` (n), `LOAD` (t×4, t = maxdofs), `rho`/`mu` (t),
//!   `Velocity` (dim×t column-major), `Pressure` (t),
//!   `CurrentDoneTime` (timestep tracking). `PrevSol`/`NonLinRes` are
//!   allocated by the source on serial trees but only read by
//!   commented-out diagnostics; they are documented, not replicated.
//! - `Stabilization` is forced `.FALSE.` by the source after reading the
//!   parameter (no side effects in the read); the port hardcodes false.
//! - `Pnorm` exists only for commented-out diagnostics; omitted.
//! - The `IF (.FALSE.)` outflow block in OptimalSolutionUpdate-style code
//!   does not exist here; the rotational outflow block is live and
//!   transcribed. The commented progress `WRITE`s/`Info` are omitted.
//! - `Norm = DefaultSolve()` kept and discarded in every arm.
//! - Serial tree: `Parallel` is queried and transcribed; the parallel arm
//!   (`DefaultDirichletBCs` + `DefaultSolve`) lives in wrappers and is
//!   compile-checked.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn csv_begin(h_model: *const c_void, h_solver: *const c_void);
    fn csv_parallel() -> c_int;
    fn csv_conv_flags(conv: *mut c_int, block: *mut c_int, found_c: *mut c_int, found_b: *mut c_int);
    fn csv_new_timestep() -> c_int;
    fn csv_const_bulk(cb: *mut c_int, cb_in_use: *mut c_int, newts: c_int);
    fn csv_ilu_recompute(newts: c_int);
    fn csv_convection(
        block: c_int,
        picard: *mut c_int,
        newton: *mut c_int,
        hybrid: *mut c_int,
        form: *mut c_int,
        nliter: *mut c_int,
    );
    fn csv_dim() -> c_int;
    fn csv_mesh_max_dofs() -> c_int;
    fn csv_matrix_rows() -> c_int;
    fn csv_bind_vars(block: c_int, ok_flow: *mut c_int, ok_pres: *mut c_int);
    fn csv_fatal_flow();
    fn csv_fatal_pres();
    fn csv_fatal_alloc();
    fn csv_alloc(n: c_int);
    fn csv_default_init();
    fn csv_bulk_restore();
    fn csv_nof_active() -> c_int;
    fn csv_begin_active(t: c_int, n: *mut c_int, nb: *mut c_int, nd: *mut c_int);
    fn csv_bodyforce(n: c_int, load: *mut c_double);
    fn csv_material(n: c_int, rho: *mut c_double, mu: *mut c_double);
    fn csv_gather(
        block: c_int,
        newts: c_int,
        iter: c_int,
        nd: c_int,
        vel: *mut c_double,
        pres: *mut c_double,
    );
    fn csv_lm_nodes();
    fn csv_lm_ip_begin() -> c_int;
    fn csv_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn csv_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        ntot: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn csv_lm_point(
        n: c_int,
        nd: c_int,
        ntot: c_int,
        dim: c_int,
        convect: c_int,
        form: c_int,
        newton: c_int,
        bulk: c_int,
        sdet: c_double,
        load: *mut c_double,
        rho_n: *mut c_double,
        mu_n: *mut c_double,
        vel_n: *mut c_double,
        pres_n: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        mass: *mut c_double,
    );
    fn csv_default_time(m_ld: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn csv_update_eq(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn csv_update_force(ld: c_int, f: *mut c_double);
    fn csv_finish_bulk_cond(in_use: c_int);
    fn csv_finish();
    fn csv_dirichlet();
    fn csv_solve() -> c_double;
    fn csv_nof_boundary() -> c_int;
    fn csv_begin_boundary(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        active: *mut c_int,
        family: *mut c_int,
    );
    fn csv_bc();
    fn csv_bc_associated() -> c_int;
    fn csv_bc_outflow(val: *mut c_int, found: *mut c_int);
    fn csv_gather_outflow(newts: c_int, iter: c_int, nd: c_int, vel: *mut c_double);
    fn csv_lb_nodes();
    fn csv_lb_ip_begin() -> c_int;
    fn csv_lb_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn csv_lb_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn csv_lb_normal(u: c_double, v: c_double, normal: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn csv_lb_point(
        nd: c_int,
        dim: c_int,
        newton: c_int,
        rho: c_double,
        sdet: c_double,
        vel_n: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        normal: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn csv_set_bc(dof: c_int, ndofs: c_int);
}

struct State {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    mass: Vec<c_double>,
    load: Vec<c_double>,
    rho: Vec<c_double>,
    mu: Vec<c_double>,
    vel: Vec<c_double>,
    pres: Vec<c_double>,
    m: usize,
    n: usize,
    dim: usize,
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
            pres: Vec::new(),
            m: 0,
            n: 0,
            dim: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<State> = Mutex::new(State::new());

fn ensure(state: &mut State, m: usize, dim: usize, _rows: usize) {
    let n = dim * m;
    if m <= state.m && dim == state.dim {
        return;
    }
    let mm = m.max(state.m);
    let nn = dim * mm;
    let _ = n;
    state.stiff.resize(nn * nn, 0.0);
    state.force.resize(nn, 0.0);
    state.mass.resize(nn * nn, 0.0);
    state.load.resize(mm * 4, 0.0);
    state.rho.resize(mm, 0.0);
    state.mu.resize(mm, 0.0);
    state.vel.resize(dim * mm, 0.0);
    state.pres.resize(mm, 0.0);
    state.m = mm;
    state.n = nn;
    state.dim = dim;
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { csv_begin(model, solver) };

    let parallel = unsafe { csv_parallel() };
    let (mut convect, mut block, mut found_c, mut found_b) = (0, 0, 0, 0);
    unsafe { csv_conv_flags(&mut convect, &mut block, &mut found_c, &mut found_b) };
    if found_c == 0 {
        convect = 1;
    }
    if found_b == 0 {
        block = 1;
    }

    let newts = unsafe { csv_new_timestep() };
    let (mut cb_matrix, mut cb_in_use) = (0, 0);
    unsafe { csv_const_bulk(&mut cb_matrix, &mut cb_in_use, newts) };
    if cb_matrix != 0 {
        unsafe { csv_ilu_recompute(newts) };
    }

    let (mut picard, mut newton, mut hybrid, mut form, mut nliter) = (0, 0, 0, 0, 1);
    if convect != 0 {
        unsafe {
            csv_convection(
                block,
                &mut picard,
                &mut newton,
                &mut hybrid,
                &mut form,
                &mut nliter,
            )
        };
    }
    let _ = (picard, hybrid);

    let dim = unsafe { csv_dim() } as usize;
    let m0 = unsafe { csv_mesh_max_dofs() } as usize;
    let rows = unsafe { csv_matrix_rows() } as usize;
    let _ = rows;

    if !state.allocations_done {
        let n = dim * m0;
        state.m = m0;
        state.n = n;
        state.dim = dim;
        if state.stiff.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.load.try_reserve(m0 * 4).is_err()
            || state.rho.try_reserve(m0).is_err()
            || state.mu.try_reserve(m0).is_err()
            || state.vel.try_reserve(dim * m0).is_err()
            || state.pres.try_reserve(m0).is_err()
        {
            unsafe { csv_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.load.resize(m0 * 4, 0.0);
        state.rho.resize(m0, 0.0);
        state.mu.resize(m0, 0.0);
        state.vel.resize(dim * m0, 0.0);
        state.pres.resize(m0, 0.0);
        unsafe { csv_alloc(m0 as c_int) };
        // Serial trees allocate PrevSol/NonLinRes in the source; they feed
        // only commented-out diagnostics, so no Rust mirror is kept.
        let _ = parallel;
        state.allocations_done = true;
    }

    // Stabilization is forced .FALSE. upstream right after the parameter
    // read (which has no side effects); the port hardcodes false.
    for v in state.vel.iter_mut() {
        *v = 0.0;
    }

    let (mut ok_flow, mut ok_pres) = (0, 0);
    unsafe { csv_bind_vars(block, &mut ok_flow, &mut ok_pres) };
    if block != 0 && ok_flow == 0 {
        unsafe { csv_fatal_flow() };
        std::process::abort();
    }
    if block == 0 && ok_pres == 0 {
        unsafe { csv_fatal_pres() };
        std::process::abort();
    }

    for iter in 1..=nliter {
        let active = unsafe { csv_nof_active() };
        if cb_in_use != 0 {
            unsafe { csv_bulk_restore() };
        } else {
            unsafe { csv_default_init() };
        }

        for t in 1..=active {
            let (mut n, mut nb, mut nd) = (0, 0, 0);
            unsafe { csv_begin_active(t, &mut n, &mut nb, &mut nd) };
            let m_need = (n as usize).max(nd as usize).max(nb as usize);
            ensure(&mut state, m_need, dim, rows);
            let nmax = state.n;

            unsafe { csv_bodyforce(n, state.load.as_mut_ptr()) };
            unsafe { csv_material(n, state.rho.as_mut_ptr(), state.mu.as_mut_ptr()) };
            unsafe {
                csv_gather(
                    block,
                    newts,
                    iter,
                    nd,
                    state.vel.as_mut_ptr(),
                    state.pres.as_mut_ptr(),
                )
            };

            for v in state.stiff.iter_mut() {
                *v = 0.0;
            }
            for v in state.force.iter_mut() {
                *v = 0.0;
            }
            for v in state.mass.iter_mut() {
                *v = 0.0;
            }
            unsafe { csv_lm_nodes() };
            let nip = unsafe { csv_lm_ip_begin() };
            let ntot = (nd + nb) as usize;
            let mut basis = vec![0.0 as c_double; ntot];
            let mut dbasis = vec![0.0 as c_double; ntot * 3];
            for k in 1..=nip {
                let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
                unsafe { csv_lm_ip_point(k, &mut u, &mut v, &mut w, &mut s) };
                let (mut detj, mut stat) = (0.0, 0);
                unsafe {
                    csv_lm_eleminfo(
                        u,
                        v,
                        w,
                        ntot as c_int,
                        basis.as_mut_ptr(),
                        dbasis.as_mut_ptr(),
                        &mut detj,
                        &mut stat,
                    )
                };
                let _ = stat;
                unsafe {
                    csv_lm_point(
                        n,
                        nd,
                        ntot as c_int,
                        dim as c_int,
                        convect,
                        form,
                        newton,
                        cb_in_use,
                        s * detj,
                        state.load.as_mut_ptr(),
                        state.rho.as_mut_ptr(),
                        state.mu.as_mut_ptr(),
                        state.vel.as_mut_ptr(),
                        state.pres.as_mut_ptr(),
                        basis.as_mut_ptr(),
                        dbasis.as_mut_ptr(),
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                        state.mass.as_mut_ptr(),
                    )
                };
            }

            unsafe {
                csv_default_time(
                    nmax as c_int,
                    state.mass.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
            if cb_in_use != 0 {
                unsafe { csv_update_force(nmax as c_int, state.force.as_mut_ptr()) };
            } else {
                unsafe {
                    csv_update_eq(
                        nmax as c_int,
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                    )
                };
            }
        }

        unsafe { csv_finish_bulk_cond(cb_in_use) };

        // Terms on the outflow boundary (rotational convection form only).
        if convect != 0 && form == 1 {
            let nbnd = unsafe { csv_nof_boundary() };
            for t in 1..=nbnd {
                let (mut n, mut nd, mut active_el, mut family) = (0, 0, 0, 0);
                unsafe { csv_begin_boundary(t, &mut n, &mut nd, &mut active_el, &mut family) };
                if active_el == 0 {
                    continue;
                }
                if family == 1 {
                    continue;
                }
                unsafe { csv_bc() };
                if unsafe { csv_bc_associated() } == 0 {
                    continue;
                }
                let (mut oval, mut ofound) = (0, 0);
                unsafe { csv_bc_outflow(&mut oval, &mut ofound) };
                let _ = ofound;
                if oval == 0 {
                    continue;
                }
                unsafe { csv_gather_outflow(newts, iter, nd, state.vel.as_mut_ptr()) };
                for v in state.stiff.iter_mut() {
                    *v = 0.0;
                }
                for v in state.force.iter_mut() {
                    *v = 0.0;
                }
                unsafe { csv_lb_nodes() };
                let nip = unsafe { csv_lb_ip_begin() };
                let mut basis = vec![0.0 as c_double; nd as usize];
                let mut dbasis = vec![0.0 as c_double; nd as usize * 3];
                let mut normal = vec![0.0 as c_double; dim];
                for k in 1..=nip {
                    let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
                    unsafe { csv_lb_ip_point(k, &mut u, &mut v, &mut w, &mut s) };
                    let (mut detj, mut stat) = (0.0, 0);
                    unsafe {
                        csv_lb_eleminfo(
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
                    unsafe { csv_lb_normal(u, v, normal.as_mut_ptr()) };
                    unsafe {
                        csv_lb_point(
                            nd,
                            dim as c_int,
                            newton,
                            state.rho[0],
                            s * detj,
                            state.vel.as_mut_ptr(),
                            basis.as_mut_ptr(),
                            dbasis.as_mut_ptr(),
                            normal.as_mut_ptr(),
                            state.stiff.as_mut_ptr(),
                            state.force.as_mut_ptr(),
                        )
                    };
                }
                unsafe {
                    csv_update_eq(
                        state.n as c_int,
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                    )
                };
            }
        }

        unsafe { csv_finish() };

        if parallel != 0 {
            unsafe { csv_dirichlet() };
            let _norm: c_double = unsafe { csv_solve() };
        } else {
            unsafe { csv_set_bc(1, dim as c_int) };
            unsafe { csv_set_bc(2, dim as c_int) };
            if dim > 2 {
                unsafe { csv_set_bc(3, dim as c_int) };
            }
            // Both arms end in DefaultSolve (the nonlinear-residual probe
            // is commented out upstream).
            let _norm: c_double = unsafe { csv_solve() };
        }
    }
}
