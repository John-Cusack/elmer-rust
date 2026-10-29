//! `pressuresolver_` — port of PressureSolver.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE` replicated: `STIFF` (n×n), `LOAD` (n), `FORCE` (n), `r`
//!   (rows), `CVelo`/`PVelo` (3×n column-major), `PPres` (n×2
//!   column-major). `TotalForce`/`TotalArea`/`ParentIndexes`/`mu` live only
//!   in the traction wrapper (same lifetimes, Fortran-side), exactly where
//!   the source uses them.
//! - `LOAD` is zeroed per element and never filled (the `Source` lines are
//!   commented out upstream); the zero fill is kept.
//! - `Norm = DefaultSolve()` is kept and discarded.
//! - The traction post-pass runs only when block preconditioning is off,
//!   guarded by the same runtime query as the source.
//! - Serial tree: `Parallel` is queried and transcribed; the parallel arms
//!   live in wrappers and compile-checked (never taken serially).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn csp_begin(h_model: *const c_void, h_solver: *const c_void);
    fn csp_parallel() -> c_int;
    fn csp_block_precond(val: *mut c_int, found: *mut c_int);
    fn csp_mesh_max_dofs() -> c_int;
    fn csp_matrix_rows() -> c_int;
    fn csp_num_bcs() -> c_int;
    fn csp_bind_vars(block: c_int, ok_velo: *mut c_int, ok_div: *mut c_int, ok_pres: *mut c_int);
    fn csp_fatal_flow();
    fn csp_fatal_velo();
    fn csp_fatal_div();
    fn csp_alloc(n: c_int, nbc: c_int);
    fn csp_default_init();
    fn csp_nof_active() -> c_int;
    fn csp_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn csp_fill_block(ntot: c_int, nd: c_int, cvelo: *mut c_double, pvelo: *mut c_double, ppres: *mut c_double);
    fn csp_fill_std(ntot: c_int, nd: c_int, cvelo: *mut c_double, pvelo: *mut c_double, ppres: *mut c_double);
    fn csp_lm_nodes();
    fn csp_lm_ip_begin() -> c_int;
    fn csp_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn csp_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn csp_lm_point(
        n: c_int,
        nd: c_int,
        ntot: c_int,
        dt: c_double,
        sdet: c_double,
        ppres: *mut c_double,
        cvelo: *mut c_double,
        pvelo: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn csp_lm_load_zero(nd: c_int, load: *mut c_double);
    fn csp_lm_correct(
        n: c_int,
        nd: c_int,
        ntot: c_int,
        ppres: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn csp_update_eq(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn csp_finish();
    fn csp_fatal_alloc();
    fn csp_solve() -> c_double;
    fn csp_nof_boundary() -> c_int;
    fn csp_begin_boundary(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        active: *mut c_int,
        family: *mut c_int,
    );
    fn csp_bc();
    fn csp_bc_associated() -> c_int;
    fn csp_bc_outflow(val: *mut c_int, found: *mut c_int);
    fn csp_outflow_surgery(block: c_int, nd: c_int, n: c_int);
    fn csp_traction();
}

struct State {
    stiff: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    r: Vec<c_double>,
    cvelo: Vec<c_double>,
    pvelo: Vec<c_double>,
    ppres: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl State {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            r: Vec::new(),
            cvelo: Vec::new(),
            pvelo: Vec::new(),
            ppres: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<State> = Mutex::new(State::new());

fn ensure(state: &mut State, need: usize, rows: usize) {
    if need <= state.n {
        return;
    }
    let n = need;
    state.stiff.resize(n * n, 0.0);
    state.load.resize(n, 0.0);
    state.force.resize(n, 0.0);
    state.cvelo.resize(3 * n, 0.0);
    state.pvelo.resize(3 * n, 0.0);
    state.ppres.resize(2 * n, 0.0);
    if state.r.len() < rows {
        state.r.resize(rows, 0.0);
    }
    state.n = n;
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { csp_begin(model, solver) };

    let _parallel = unsafe { csp_parallel() };
    let (mut block_val, mut block_found) = (0, 0);
    unsafe { csp_block_precond(&mut block_val, &mut block_found) };
    let _ = block_found;
    let block = block_val;

    let dt_v: c_double = unsafe { *dt };

    // Variable binding with the source's Fatal checks, in order.
    let (mut ok_velo, mut ok_div, mut ok_pres) = (0, 0, 0);
    unsafe { csp_bind_vars(block, &mut ok_velo, &mut ok_div, &mut ok_pres) };
    if block != 0 && ok_pres == 0 {
        unsafe { csp_fatal_flow() };
        std::process::abort();
    }
    if ok_velo == 0 {
        unsafe { csp_fatal_velo() };
        std::process::abort();
    }
    if ok_div == 0 {
        unsafe { csp_fatal_div() };
        std::process::abort();
    }

    if !state.allocations_done {
        let n = unsafe { csp_mesh_max_dofs() } as usize;
        let t = unsafe { csp_matrix_rows() } as usize;
        let nbc = unsafe { csp_num_bcs() };
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.r.try_reserve(t).is_err()
            || state.cvelo.try_reserve(3 * n).is_err()
            || state.pvelo.try_reserve(3 * n).is_err()
            || state.ppres.try_reserve(2 * n).is_err()
        {
            unsafe { csp_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.r.resize(t, 0.0);
        state.cvelo.resize(3 * n, 0.0);
        state.pvelo.resize(3 * n, 0.0);
        state.ppres.resize(2 * n, 0.0);
        unsafe { csp_alloc(n as c_int, nbc) };
        state.allocations_done = true;
    }

    // System assembly.
    let active = unsafe { csp_nof_active() };
    unsafe { csp_default_init() };
    for t in 1..=active {
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { csp_begin_active(t, &mut n, &mut nd, &mut nb) };
        let ntot = (nd + nb) as usize;
        let rows = state.r.len();
        ensure(&mut state, ntot.max(n as usize), rows);
        let nmax = state.n;

        if block != 0 {
            unsafe {
                csp_fill_block(
                    ntot as c_int,
                    nd,
                    state.cvelo.as_mut_ptr(),
                    state.pvelo.as_mut_ptr(),
                    state.ppres.as_mut_ptr(),
                )
            };
        } else {
            unsafe {
                csp_fill_std(
                    ntot as c_int,
                    nd,
                    state.cvelo.as_mut_ptr(),
                    state.pvelo.as_mut_ptr(),
                    state.ppres.as_mut_ptr(),
                )
            };
        }

        unsafe { csp_lm_load_zero(nd, state.load.as_mut_ptr()) };

        for v in state.stiff.iter_mut() {
            *v = 0.0;
        }
        for v in state.force.iter_mut() {
            *v = 0.0;
        }
        unsafe { csp_lm_nodes() };
        let nip = unsafe { csp_lm_ip_begin() };
        let mut basis = vec![0.0 as c_double; ntot];
        let mut dbasis = vec![0.0 as c_double; ntot * 3];
        for k in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { csp_lm_ip_point(k, &mut u, &mut v, &mut w, &mut s) };
            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                csp_lm_eleminfo(
                    u,
                    v,
                    w,
                    n,
                    nd,
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    &mut detj,
                    &mut stat,
                )
            };
            let _ = stat;
            unsafe {
                csp_lm_point(
                    n,
                    nd,
                    ntot as c_int,
                    dt_v,
                    s * detj,
                    state.ppres.as_mut_ptr(),
                    state.cvelo.as_mut_ptr(),
                    state.pvelo.as_mut_ptr(),
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }

        unsafe {
            csp_lm_correct(
                n,
                nd,
                ntot as c_int,
                state.ppres.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };

        unsafe {
            csp_update_eq(
                nmax as c_int,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }

    unsafe { csp_finish() };

    // Outflow-boundary matrix surgery.
    let nbnd = unsafe { csp_nof_boundary() };
    for t in 1..=nbnd {
        let (mut n, mut nd, mut active_el, mut family) = (0, 0, 0, 0);
        unsafe { csp_begin_boundary(t, &mut n, &mut nd, &mut active_el, &mut family) };
        if active_el == 0 {
            continue;
        }
        if family == 1 {
            continue;
        }
        unsafe { csp_bc() };
        if unsafe { csp_bc_associated() } == 0 {
            continue;
        }
        let (mut oval, mut ofound) = (0, 0);
        unsafe { csp_bc_outflow(&mut oval, &mut ofound) };
        let _ = ofound;
        if oval != 0 {
            unsafe { csp_outflow_surgery(block, nd, n) };
        }
    }

    // And finally, solve (Norm local and unused, as in the original).
    let _norm: c_double = unsafe { csp_solve() };

    // Traction post-pass only without block preconditioning.
    if block == 0 {
        unsafe { csp_traction() };
    }
}
