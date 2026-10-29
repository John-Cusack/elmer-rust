//! `compressibilitysolver_` — port of CompressibilitySolver.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE` replicated: `MASS`/`STIFF` (n×n), `FORCE` (n),
//!   `CVelo` (dim×n), `Indexes` (n), `AllocationsDone`, `CurrentDoneTime`
//!   (declared but never read — replicated as a counter like the source's
//!   unused increment site, which this entry lacks; kept for shape parity
//!   with sibling entries that test it).
//! - Dead declarations `Pvelo`, `ppres`, `LOAD` are never allocated or
//!   touched by the source; they are not replicated.
//! - `VeloVar` null is not tested by the source (no ASSOCIATED check), so
//!   the port adds none: a missing `VelocityTot` crashes identically.
//! - `LOAD`/`Pvelo`/`ppres` are absent; `MASS` is zeroed per element and
//!   never read back — passed as scratch, exactly like the source.
//! - `Norm = DefaultSolve()` is kept and discarded.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn csc_begin(h_model: *const c_void, h_solver: *const c_void);
    fn csc_dim() -> c_int;
    fn csc_max_elem_dofs() -> c_int;
    fn csc_const_bulk(val: *mut c_int, found: *mut c_int);
    fn csc_const_bulk_in_use(cb: c_int) -> c_int;
    fn csc_velovar();
    fn csc_alloc_indexes(n: c_int);
    fn csc_default_init(cb: c_int);
    fn csc_nof_active() -> c_int;
    fn csc_begin_active(t: c_int, n: *mut c_int, nb: *mut c_int, nd: *mut c_int);
    fn csc_fill_cvelo(dim: c_int, nd: c_int, cvelo: *mut c_double);
    fn csc_viscosity() -> c_double;
    fn csc_lm_nodes();
    fn csc_lm_ip_begin() -> c_int;
    fn csc_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn csc_lm_eleminfo(
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
    fn csc_lm_point(
        dim: c_int,
        n: c_int,
        nd: c_int,
        ntot: c_int,
        cb: c_int,
        visc: c_double,
        sdet: c_double,
        cvelo: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn csc_update_eq(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn csc_update_force(ld: c_int, f: *mut c_double);
    fn csc_finish_bulk(cb_matrix: c_int, cb_in_use: c_int);
    fn csc_finish();
    fn csc_dirichlet();
    fn csc_solve() -> c_double;
    fn csc_fatal_alloc();
}

struct State {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    mass: Vec<c_double>,
    cvelo: Vec<c_double>,
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
            cvelo: Vec::new(),
            n: 0,
            dim: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<State> = Mutex::new(State::new());

/// Grow-only guard (adaptive remesh across invocations); the source sizes
/// once from the first mesh and would overflow silently on a finer remesh.
fn ensure(state: &mut State, need: usize, dim: usize) {
    if need <= state.n && dim == state.dim {
        return;
    }
    let n = need.max(state.n);
    state.stiff.resize(n * n, 0.0);
    state.force.resize(n, 0.0);
    state.mass.resize(n * n, 0.0);
    state.cvelo.resize(dim * n, 0.0);
    state.n = n;
    state.dim = dim;
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { csc_begin(model, solver) };

    let dim = unsafe { csc_dim() } as usize;
    let (mut cb_val, mut cb_found) = (0, 0);
    unsafe { csc_const_bulk(&mut cb_val, &mut cb_found) };
    let _ = cb_found;
    let cb_in_use = unsafe { csc_const_bulk_in_use(cb_val) };

    unsafe { csc_velovar() };

    // First-mesh permanent storage (bug-compatible: sized once).
    if !state.allocations_done {
        let n = unsafe { csc_max_elem_dofs() } as usize;
        state.n = n;
        state.dim = dim;
        if state.stiff.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.cvelo.try_reserve(dim * n).is_err()
        {
            unsafe { csc_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.cvelo.resize(dim * n, 0.0);
        unsafe { csc_alloc_indexes(n as c_int) };
        state.allocations_done = true;
    }

    let active = unsafe { csc_nof_active() };
    unsafe { csc_default_init(cb_in_use) };

    for t in 1..=active {
        let (mut n, mut nb, mut nd) = (0, 0, 0);
        unsafe { csc_begin_active(t, &mut n, &mut nb, &mut nd) };
        let ntot = (nd + nb) as usize;
        ensure(&mut state, ntot.max(n as usize), dim);
        let nmax = state.n;
        unsafe { csc_fill_cvelo(dim as c_int, nd, state.cvelo.as_mut_ptr()) };
        let visc = unsafe { csc_viscosity() };
        // LocalMatrix: zero + integrate + bubble elimination (wrapper owns
        // the per-IP FLOPs; Rust owns the IP loop and the zero fills).
        for v in state.stiff.iter_mut() {
            *v = 0.0;
        }
        for v in state.force.iter_mut() {
            *v = 0.0;
        }
        for v in state.mass.iter_mut() {
            *v = 0.0;
        }
        unsafe { csc_lm_nodes() };
        let nip = unsafe { csc_lm_ip_begin() };
        let mut basis = vec![0.0 as c_double; ntot];
        let mut dbasis = vec![0.0 as c_double; ntot * 3];
        for k in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { csc_lm_ip_point(k, &mut u, &mut v, &mut w, &mut s) };
            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                csc_lm_eleminfo(
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
                csc_lm_point(
                    dim as c_int,
                    n,
                    nd,
                    ntot as c_int,
                    cb_in_use,
                    visc,
                    s * detj,
                    state.cvelo.as_mut_ptr(),
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }

        if cb_in_use != 0 {
            unsafe { csc_update_force(nmax as c_int, state.force.as_mut_ptr()) };
        } else {
            unsafe {
                csc_update_eq(
                    nmax as c_int,
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }
    }

    unsafe { csc_finish_bulk(cb_val, cb_in_use) };
    unsafe { csc_finish() };
    unsafe { csc_dirichlet() };

    // Solved norm is local and unused, as in the original.
    let _norm: c_double = unsafe { csc_solve() };
}
