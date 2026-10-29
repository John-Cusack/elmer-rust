//! Port of `fem/src/modules/Stokes.F90` (the Elmer fork) for lane T021-B,
//! at the solver-plugin boundary (whole `.so` → Rust cdylib with the same
//! dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/Stokes.F90` — one F77-global entry
//! `stokessolver_` with `CONTAINS` internals. The emitted dynamic symbol is
//! exactly `stokessolver_` (proven by `nm -D` on the reference `.so`); the
//! internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/stokes_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`MATMUL`, `SUM`, `DOT_PRODUCT`), formatted I/O
//! (`PRINT`/`WRITE`+`Info`), and whole-array updates stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors Stokes.F90 statement for statement):
//! - `SAVE` semantics replicated: first-mesh sizing (`MaxElementDOFs`)
//!   fixes `STIFF`/`Mass`/`LOAD`/`FORCE`/`rho`/`mu`/`Velocity`/`ALocal`/
//!   `MLocal`/`PLocal`/`PrevSol`/`NonLinRes` plus `AllocationsDone`,
//!   `CurrentDoneTime`, `Norm`, and `NumberOfNormalTractionNodes`; the
//!   wrapper side mirrors the internal-procedure `SAVE Nodes`/`SAVE IP`
//!   state and the `SAVE Indexes` / `SAVE TractionBCIndeces` arrays with
//!   identical lifetimes (consumed before the next fill, as in the
//!   original). Growth guards only trigger on remeshed re-entry, where
//!   the original would overflow silently.
//! - `Newton`/`Hybrid` carry the Fortran initialization-implies-`SAVE`
//!   latch but are reassigned from `Nonlinear Iteration Method` on every
//!   entry before any internal procedure reads them, exactly as in the
//!   original; `PicardIteration`/`NewTimeStep` are write-only in the
//!   original and are omitted (no observable effect).
//! - The `IF (.FALSE.)` convergence-study block (and its callees
//!   `ComputePressureError` / `ComputeVelocityError`), the `IF (.FALSE.)`
//!   standard-convection form, the `IF (.TRUE.) ... ELSE` alternative
//!   stabilization/boundary branches, and the `IF (.FALSE.)`-guarded
//!   `LCondensateStabilizationBubble` call + callee are statically dead
//!   and are omitted (documented, unobservable).
//! - `GradDivParam` is accepted (it travels into `LocalMatrix`) but unused
//!   in the body, exactly as in the original.
//! - `BiCGStab` is transcribed (it is only referenced by commented-out
//!   call sites) and marked `#[allow(dead_code)]`; it is unreachable.
//! - Allocation failure calls the same `Fatal('NavierStokesSolver', ...)`
//!   via `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `Model`, `dt`, and `TransientSimulation` are received and used
//!   exactly as in the original (`Model` for `SetBoundaryConditions` via
//!   wrappers and the block scan; `dt` for the GCR call; the transient
//!   flag for the `PrevValues`/`Default1stOrderTime` selection).
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads. (No
//!   Rust-defined callbacks are registered with the solver, so holding
//!   the guard across the solve cannot deadlock.)

use std::ffi::c_void;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sta_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn sta_done_time() -> c_int;
    fn sta_dim() -> c_int;
    fn sta_max_elem_dofs() -> c_int;
    fn sta_num_rows() -> c_int;
    fn sta_nof_active() -> c_int;
    fn sta_num_bnd_elems() -> c_int;
    fn sta_list_logical(which: c_int, key: *const c_char, keylen: c_int, found: *mut c_int) -> c_int;
    fn sta_list_const_real(
        key: *const c_char, keylen: c_int, found: *mut c_int, use_minv: c_int, minv: c_double,
        with_found: c_int,
    ) -> c_double;
    fn sta_list_integer(
        key: *const c_char, keylen: c_int, found: *mut c_int, use_minv: c_int, minv: c_int,
        with_found: c_int,
    ) -> c_int;
    fn sta_list_creal(key: *const c_char, keylen: c_int, found: *mut c_int) -> c_double;
    fn sta_get_string(
        which: c_int, idx: c_int, buf: *mut c_char, buflen: c_int, found: *mut c_int,
    );
    fn sta_fatal_alloc();
    fn sta_block_find();
    fn sta_block_alloc() -> c_int;
    fn sta_block_fill();
    fn sta_block_rows(which: c_int) -> c_int;
    fn sta_zero_crs(which: c_int);
    fn sta_cpu_time() -> c_double;
    fn sta_real_time() -> c_double;
    fn sta_default_initialize();
    fn sta_info_progress(t: c_int);
    fn sta_info_nonlin(res: c_double);
    fn sta_info_snorm(norm: c_double);
    fn sta_info_atime(atime: c_double);
    fn sta_info_stime(stime: c_double);
    fn sta_outer_timelevel(cur: c_int);
    fn sta_outer_nonlin_iter(iter: c_int);
    fn sta_print_msg(which: c_int);
    fn sta_gcr_info0(res: c_double);
    fn sta_gcr_info(k: c_int, res: c_double);
    fn sta_bicg_write(i: c_int, res: c_double);
    fn sta_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn sta_bodyforce() -> c_int;
    fn sta_load_force(k: c_int, n: c_int, load: *mut c_double, found: *mut c_int);
    fn sta_mat_rho_mu(n: c_int, rho: *mut c_double, mu: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sta_gather_vel(which: c_int, dim: c_int, ldv: c_int, nd: c_int, vel: *mut c_double);
    fn sta_ip_begin(nip: *mut c_int, ch: *mut c_double);
    fn sta_schur_begin(nip: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn sta_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sta_lm_eleminfo(
        u: c_double, v: c_double, w: c_double, n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double, detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn sta_lm_eleminfo2(
        u: c_double, v: c_double, w: c_double, n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double, ddb: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn sta_lm_velo(
        dim: c_int, ldv: c_int, nd: c_int, vel: *const c_double, basis: *const c_double,
        out: *mut c_double,
    );
    fn sta_lm_w(
        need12: c_int, ldv: c_int, nd: c_int, vel: *const c_double, dd: *const c_double,
        w: *mut c_double,
    );
    fn sta_lm_grad(ldv: c_int, nd: c_int, vel: *const c_double, dd: *const c_double, g: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sta_lm_murho(
        n: c_int, basis: *const c_double, nodalmu: *const c_double, nodalrho: *const c_double,
        mu: *mut c_double, rho: *mut c_double,
    );
    fn sta_lm_loadatip(
        n: c_int, dim: c_int, ldl: c_int, load: *const c_double, basis: *const c_double,
        out: *mut c_double,
    );
    fn sta_lm_hk_mk(hk: *mut c_double, mk: *mut c_double);
    fn sta_lb_normal(u: c_double, v: c_double, nl: *mut c_double);
    fn sta_default_1st_order(ld: c_int, m: *mut c_double, a: *mut c_double, f: *mut c_double);
    fn sta_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn sta_update_precond(which: c_int, ld: c_int, nd: c_int, ndofs: c_int, local: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sta_begin_boundary(
        t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int, family: *mut c_int,
    );
    fn sta_bc();
    fn sta_bc_associated() -> c_int;
    fn sta_default_finish_bulk();
    fn sta_default_finish_assembly();
    fn sta_set_bcs(mat: c_int, dof: c_int, ndofs: c_int, with_rhs: c_int);
    fn sta_prev_sol(which: c_int, n: c_int, prev: *mut c_double);
    fn sta_matvec_rhs(n: c_int, prev: *const c_double, res: *mut c_double);
    fn sta_pnorm_vec(n: c_int, x: *const c_double) -> c_double;
    fn sta_rhs_pnorm() -> c_double;
    fn sta_values_pnorm(n: c_int) -> c_double;
    fn sta_pdot(n: c_int, x: *const c_double, y: *const c_double) -> c_double;
    fn sta_relax(n: c_int, prev: *const c_double, factor: c_double);
    fn sta_default_solve() -> c_double;
    fn sta_values_copy(x: *mut c_double, n: c_int);
    fn sta_rhs_copy(b: *mut c_double, n: c_int);
    fn sta_values_set(x: *const c_double, n: c_int);
    fn sta_pmv_main(n: c_int, x: *const c_double, b: *mut c_double);
    fn sta_crs_matvec_which(which: c_int, n: c_int, x: *const c_double, b: *mut c_double);
    fn sta_crs_lusolve(which: c_int, n: c_int, b: *mut c_double);
    fn sta_crs_divmult(dim: c_int, u: *const c_double, v: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn sta_solve_lin(
        which: c_int, n: c_int, f: *const c_double, y: *mut c_double, dofs: c_int,
        norm: *mut c_double,
    );
    fn sta_prec_traction(fvec: *mut c_double, dvec: *const c_double);
    fn sta_precomp_toggle(which: c_int, val: c_int);
    fn sta_get_max_restarts(found: *mut c_int) -> c_int;
}

/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (remesh across invocations). The original never grows (its
/// allocation assumes the first mesh is big enough); growth only triggers
/// where the original would overflow silently, and the leading sections
/// used by the wrappers are unaffected, so valid inputs never observe it.
fn ensure_state(state: &mut StokesState, p: usize, nrows: usize) {
    if p <= state.p && nrows <= state.nrows {
        return;
    }
    let p = p.max(state.p);
    let nrows = nrows.max(state.nrows);
    let dim = state.dim;
    let n = (dim + 1) * p;
    let m = dim * p;
    state.stiff.resize(n * n, 0.0);
    state.mass.resize(n * n, 0.0);
    state.force.resize(n, 0.0);
    state.load.resize(p * 4, 0.0);
    state.rho.resize(p, 0.0);
    state.mu.resize(p, 0.0);
    state.velocity.resize((dim + 1) * p, 0.0);
    state.alocal.resize(m * m, 0.0);
    state.mlocal.resize(p * p, 0.0);
    state.plocal.resize(p * p, 0.0);
    state.prevsol.resize(nrows, 0.0);
    state.nonlinres.resize(nrows, 0.0);
    state.p = p;
    state.n = n;
    state.m = m;
    state.nrows = nrows;
}

/// First-mesh SAVE state. `p` is the first call's `MaxElementDOFs`,
/// `n = (dim+1)*p`, `m = dim*p`, plus the `AllocationsDone`,
/// `CurrentDoneTime`, `NumberOfNormalTractionNodes`, and `Norm` latches.
struct StokesState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    rho: Vec<c_double>,
    mu: Vec<c_double>,
    velocity: Vec<c_double>,
    alocal: Vec<c_double>,
    mlocal: Vec<c_double>,
    plocal: Vec<c_double>,
    prevsol: Vec<c_double>,
    nonlinres: Vec<c_double>,
    p: usize,
    n: usize,
    m: usize,
    nrows: usize,
    dim: usize,
    allocations_done: bool,
    current_done_time: c_int,
    num_traction: usize,
    norm: c_double,
}

impl StokesState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            rho: Vec::new(),
            mu: Vec::new(),
            velocity: Vec::new(),
            alocal: Vec::new(),
            mlocal: Vec::new(),
            plocal: Vec::new(),
            prevsol: Vec::new(),
            nonlinres: Vec::new(),
            p: 0,
            n: 0,
            m: 0,
            nrows: 0,
            dim: 0,
            allocations_done: false,
            current_done_time: 0,
            num_traction: 0,
            norm: 0.0,
        }
    }
}

static STATE: Mutex<StokesState> = Mutex::new(StokesState::new());

fn list_logical(which: c_int, key: &[u8]) -> (c_int, c_int) {
    let mut found = 0;
    let v = unsafe {
        sta_list_logical(
            which,
            key.as_ptr() as *const c_char,
            key.len() as c_int,
            &mut found,
        )
    };
    (v, found)
}

fn list_const_real(
    key: &[u8], use_minv: bool, minv: c_double, with_found: bool,
) -> (c_double, c_int) {
    let mut found = 0;
    let v = unsafe {
        sta_list_const_real(
            key.as_ptr() as *const c_char,
            key.len() as c_int,
            &mut found,
            i32::from(use_minv),
            minv,
            i32::from(with_found),
        )
    };
    (v, found)
}

fn list_integer(
    key: &[u8], use_minv: bool, minv: c_int, with_found: bool,
) -> (c_int, c_int) {
    let mut found = 0;
    let v = unsafe {
        sta_list_integer(
            key.as_ptr() as *const c_char,
            key.len() as c_int,
            &mut found,
            i32::from(use_minv),
            minv,
            i32::from(with_found),
        )
    };
    (v, found)
}

fn list_creal(key: &[u8]) -> (c_double, c_int) {
    let mut found = 0;
    let v = unsafe {
        sta_list_creal(key.as_ptr() as *const c_char, key.len() as c_int, &mut found)
    };
    (v, found)
}

fn get_string(which: c_int, idx: c_int) -> String {
    let mut buf = [0u8; 256];
    let mut found = 0;
    unsafe {
        sta_get_string(
            which,
            idx,
            buf.as_mut_ptr() as *mut c_char,
            buf.len() as c_int,
            &mut found,
        );
    }
    String::from_utf8_lossy(&buf).trim_end().to_string()
}

/// `stokessolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads (they select the transient and
/// GCR paths, as in the original).
#[export_name = "stokessolver_"]
pub unsafe extern "C" fn stokes_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sta_begin(model, solver) };
    let transient_simulation = unsafe { *transient } != 0;
    let dt_val: c_double = unsafe { *dt };

    // IF (CurrentDoneTime /= Solver % DoneTime) CurrentDoneTime += 1
    if state.current_done_time != unsafe { sta_done_time() } {
        state.current_done_time += 1;
    }

    let dim = unsafe { sta_dim() } as usize;
    // Convection stabilization is probed, then unconditionally disabled.
    let _ = list_logical(0, b"Stabilize");
    let convection_stabilization = false;
    // Grad-div stabilization (unavailable option; kept for fidelity).
    let (gds_val, gds_found) = list_logical(0, b"Grad-Div Stabilization");
    // The option is unavailable; both branches assign the same values.
    let grad_div_stabilization = gds_found != 0 && gds_val != 0;
    let grad_div_param = if grad_div_stabilization {
        let (v, f) = list_const_real(b"Grad-Div Stabilization Parameter", false, 0.0, true);
        if f != 0 { v } else { 0.0 }
    } else {
        0.0
    };
    let _ = (grad_div_stabilization, grad_div_param);

    let nim = get_string(0, 0);
    let (mut newton, mut hybrid);
    if nim == "hybrid" {
        newton = true;
        hybrid = true;
    } else if nim != "picard" {
        newton = true;
        hybrid = false;
    } else {
        newton = false;
        hybrid = false;
    }

    let (block_val, _) = list_logical(0, b"Block Preconditioning");
    let block_preconditioning = block_val != 0;
    if block_preconditioning {
        unsafe { sta_outer_timelevel(state.current_done_time) };
    }
    let (mut tolerance, mut max_iterations) = (0.0, 0);
    if block_preconditioning {
        (tolerance, _) = list_const_real(b"Linear System Convergence Tolerance", false, 0.0, false);
        (max_iterations, _) = list_integer(b"Linear System Max Iterations", false, 0, false);
    }
    let _ = max_iterations;

    // First-mesh permanent storage (bug-compatible: sized once).
    if !state.allocations_done {
        let p = unsafe { sta_max_elem_dofs() } as usize;
        let n = (dim + 1) * p;
        let m = dim * p;
        let nrows = unsafe { sta_num_rows() } as usize;
        state.dim = dim;
        if state.stiff.try_reserve(n * n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.load.try_reserve(p * 4).is_err()
            || state.rho.try_reserve(p).is_err()
            || state.mu.try_reserve(p).is_err()
            || state.velocity.try_reserve((dim + 1) * p).is_err()
            || state.alocal.try_reserve(m * m).is_err()
            || state.mlocal.try_reserve(p * p).is_err()
            || state.plocal.try_reserve(p * p).is_err()
            || state.prevsol.try_reserve(nrows).is_err()
            || state.nonlinres.try_reserve(nrows).is_err()
        {
            unsafe { sta_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.force.resize(n, 0.0);
        state.load.resize(p * 4, 0.0);
        state.rho.resize(p, 0.0);
        state.mu.resize(p, 0.0);
        state.velocity.resize((dim + 1) * p, 0.0);
        state.alocal.resize(m * m, 0.0);
        state.mlocal.resize(p * p, 0.0);
        state.plocal.resize(p * p, 0.0);
        state.prevsol.resize(nrows, 0.0);
        state.nonlinres.resize(nrows, 0.0);
        state.p = p;
        state.n = n;
        state.m = m;
        state.nrows = nrows;
        if block_preconditioning {
            state.num_traction = unsafe { sta_block_alloc() } as usize;
        }
        state.allocations_done = true;
    } else {
        let p = unsafe { sta_max_elem_dofs() } as usize;
        let nrows = unsafe { sta_num_rows() } as usize;
        ensure_state(&mut state, p, nrows);
    }

    if block_preconditioning {
        unsafe { sta_block_find() };
    }

    let (conv_val, conv_found) = list_logical(1, b"Convective");
    let mut convect = conv_val != 0;
    if conv_found == 0 {
        convect = true;
    }

    let mut atime = unsafe { sta_cpu_time() };
    let mut at0 = unsafe { sta_real_time() };

    let (nonlinear_iter, _) = list_integer(b"Nonlinear System Max Iterations", true, 0, false);
    let mut nonlinear_iter = nonlinear_iter;
    let (nonlinear_tol, _) =
        list_const_real(b"Nonlinear System Convergence Tolerance", true, 0.0, false);
    if !convect {
        nonlinear_iter = 0;
    }

    let mut nonlin_error = 0.0;
    let mut iter = 1;
    while iter <= nonlinear_iter + 1 {
        let (nt_val, nt_found) =
            list_const_real(b"Nonlinear System Newton After Tolerance", false, 0.0, true);
        if nt_found != 0 && iter > 1 && nonlin_error < nt_val {
            newton = true;
            hybrid = false;
        }
        let (mpi_val, mpi_found) =
            list_integer(b"Nonlinear System Newton After Iterations", false, 0, true);
        if mpi_found != 0 && iter > mpi_val {
            newton = true;
            hybrid = false;
        }

        unsafe { sta_default_initialize() };
        if block_preconditioning {
            unsafe {
                sta_zero_crs(0);
                sta_zero_crs(1);
                sta_zero_crs(2);
            }
        }

        let active = unsafe { sta_nof_active() };
        for t in 1..=active {
            if unsafe { sta_real_time() } - at0 > 1.0 {
                unsafe { sta_info_progress(t) };
                at0 = unsafe { sta_real_time() };
            }
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { sta_begin_active(t, &mut n, &mut nd, &mut nb) };
            let _ = nb;
            let (need_p, need_rows) = (state.p.max(nd as usize), state.nrows);
            ensure_state(&mut state, need_p, need_rows);
            let (ld, p_alloc) = (state.n as c_int, state.p);
            let dim_i = dim as c_int;

            state.load.fill(0.0);
            if unsafe { sta_bodyforce() } != 0 {
                for k in 1..=3 {
                    let mut found = 0;
                    let col = state.load.as_mut_ptr().wrapping_add((k as usize - 1) * p_alloc);
                    unsafe { sta_load_force(k, n, col, &mut found) };
                    let _ = found;
                }
            }

            unsafe {
                sta_mat_rho_mu(
                    n,
                    state.rho.as_mut_ptr(),
                    state.mu.as_mut_ptr(),
                );
            }

            state.velocity.fill(0.0);
            if convect {
                let which = if transient_simulation {
                    if iter == 1 { 0 } else { 1 }
                } else {
                    1
                };
                unsafe {
                    sta_gather_vel(
                        which,
                        dim_i,
                        (dim + 1) as c_int,
                        nd,
                        state.velocity.as_mut_ptr(),
                    );
                }
            }

            local_matrix(
                &mut state,
                n,
                nd,
                dim,
                convection_stabilization,
                convect,
                newton,
                hybrid,
                grad_div_param,
            );

            if transient_simulation {
                unsafe {
                    sta_default_1st_order(
                        ld,
                        state.mass.as_mut_ptr(),
                        state.stiff.as_mut_ptr(),
                        state.force.as_mut_ptr(),
                    );
                }
            }
            unsafe {
                sta_default_update(ld, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
            }

            if block_preconditioning {
                let m = state.m;
                let nld = state.n;
                state.alocal.fill(0.0);
                for p in 0..nd as usize {
                    for i in 0..dim {
                        for q in 0..nd as usize {
                            for j in 0..dim {
                                state.alocal[(p * dim + i) + (q * dim + j) * m] =
                                    state.stiff[(p * (dim + 1) + i)
                                        + (q * (dim + 1) + j) * nld];
                            }
                        }
                    }
                }
                unsafe {
                    sta_update_precond(0, m as c_int, nd, dim_i, state.alocal.as_mut_ptr());
                }
                schur_complement(&mut state, n, nd, dim);
                unsafe {
                    sta_update_precond(1, p_alloc as c_int, nd, 1, state.mlocal.as_mut_ptr());
                    sta_update_precond(2, p_alloc as c_int, nd, 1, state.plocal.as_mut_ptr());
                }
            }
        }
        unsafe { sta_default_finish_bulk() };

        if convect {
            let nbnd = unsafe { sta_num_bnd_elems() };
            for t in 1..=nbnd {
                let (mut n, mut nd, mut active_el, mut family) = (0, 0, 0, 0);
                unsafe { sta_begin_boundary(t, &mut n, &mut nd, &mut active_el, &mut family) };
                if active_el == 0 {
                    continue;
                }
                if family == 1 {
                    continue;
                }
                unsafe { sta_bc() };
                if unsafe { sta_bc_associated() } != 0 {
                    let (ob_val, _) = list_logical(2, b"Outflow Boundary");
                    if ob_val != 0 {
                        state.velocity.fill(0.0);
                        let which = if transient_simulation {
                            if iter == 1 { 0 } else { 1 }
                        } else {
                            1
                        };
                        unsafe {
                            sta_gather_vel(
                                which,
                                dim as c_int,
                                (dim + 1) as c_int,
                                nd,
                                state.velocity.as_mut_ptr(),
                            );
                        }
                        local_matrix_boundary(&mut state, n, nd, dim, newton, block_preconditioning);
                        let ld = state.n as c_int;
                        unsafe {
                            sta_default_update(
                                ld,
                                state.stiff.as_mut_ptr(),
                                state.force.as_mut_ptr(),
                            );
                        }
                        if block_preconditioning {
                            let m = state.m as c_int;
                            unsafe {
                                sta_update_precond(
                                    0,
                                    m,
                                    nd,
                                    dim as c_int,
                                    state.alocal.as_mut_ptr(),
                                );
                            }
                        }
                    }
                }
            }
        }

        unsafe { sta_default_finish_assembly() };

        let dim1 = (dim + 1) as c_int;
        let dim_i = dim as c_int;
        unsafe {
            sta_set_bcs(0, 1, dim1, 1);
            sta_set_bcs(0, 2, dim1, 1);
            sta_set_bcs(0, 3, dim1, 1);
        }
        if dim > 2 {
            unsafe { sta_set_bcs(0, 4, dim1, 1) };
        }

        if block_preconditioning {
            unsafe {
                sta_set_bcs(1, 1, dim_i, 0);
                sta_set_bcs(1, 2, dim_i, 0);
            }
            if dim > 2 {
                unsafe { sta_set_bcs(1, 3, dim_i, 0) };
            }
            if state.num_traction > 0 {
                unsafe { sta_block_fill() };
            }
        }

        atime = unsafe { sta_cpu_time() } - atime;

        let nrows = unsafe { sta_num_rows() };
        if iter == 1 && transient_simulation {
            unsafe { sta_prev_sol(0, nrows, state.prevsol.as_mut_ptr()) };
        } else {
            unsafe { sta_prev_sol(1, nrows, state.prevsol.as_mut_ptr()) };
        }

        unsafe { sta_matvec_rhs(nrows, state.prevsol.as_ptr(), state.nonlinres.as_mut_ptr()) };
        nonlin_error =
            unsafe { sta_pnorm_vec(nrows, state.nonlinres.as_ptr()) } / unsafe { sta_rhs_pnorm() };

        unsafe { sta_info_nonlin(nonlin_error) };

        if (nonlin_error < nonlinear_tol || iter == nonlinear_iter + 1) && convect {
            break;
        }
        let mut stime = unsafe { sta_cpu_time() };
        if block_preconditioning {
            unsafe { sta_outer_nonlin_iter(iter) };
            let l = unsafe { sta_block_rows(0) };
            let q = 2 * unsafe { sta_block_rows(1) };
            let mu1 = if state.mu.is_empty() { 0.0 } else { state.mu[0] };
            gcr_outer(
                &mut state,
                max_iterations,
                tolerance,
                dim,
                l,
                q,
                mu1,
                dt_val,
                transient_simulation,
            );
        } else {
            let _ = unsafe { sta_default_solve() };
        }
        let norm = unsafe { sta_values_pnorm(nrows) };
        state.norm = norm;
        unsafe { sta_info_snorm(norm) };

        stime = unsafe { sta_cpu_time() } - stime;
        unsafe {
            sta_info_atime(atime);
            sta_info_stime(stime);
        }

        let (rf_val, rf_found) = list_creal(b"Nonlinear System Relaxation Factor");
        if rf_found != 0 {
            unsafe { sta_relax(nrows, state.prevsol.as_ptr(), rf_val) };
        }
        iter += 1;
    }
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count.
#[allow(clippy::too_many_arguments)]
fn local_matrix(
    state: &mut StokesState,
    n: c_int,
    nd: c_int,
    dim: usize,
    stabilization: bool,
    convect: bool,
    newton: bool,
    hybrid: bool,
    _grad_div_param: c_double,
) {
    let ld = state.n;
    let nn = n as usize;
    let ndn = nd as usize;
    let d = dim;
    state.stiff.fill(0.0);
    state.force.fill(0.0);
    state.mass.fill(0.0);

    // Basis(nd) / dBasisdx(nd,3) automatics, sized by this element's DOFs
    // (bubble-inclusive, as in the original).
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let (mut nip, mut ch) = (0, 0.0);
    unsafe { sta_ip_begin(&mut nip, &mut ch) };

    let mut ak = 0.0;
    let mut rotterm = 0.0;
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { sta_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            sta_lm_eleminfo(
                u, v, w, n, nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;
        let sdet = s * detj;
        ak += sdet;

        let mut velo = [0.0 as c_double; 3];
        unsafe {
            sta_lm_velo(
                d as c_int,
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                basis.as_ptr(),
                velo.as_mut_ptr(),
            );
        }
        let mut wv = [0.0 as c_double; 3];
        let need12 = i32::from(d > 2 && (newton || stabilization));
        unsafe {
            sta_lm_w(
                need12,
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                dbasis.as_ptr(),
                wv.as_mut_ptr(),
            );
        }
        if stabilization {
            if d > 2 {
                rotterm += sdet
                    * ((wv[1] * velo[2] - wv[2] * velo[1]).powi(2)
                        + (wv[2] * velo[0] - wv[0] * velo[2]).powi(2)
                        + (wv[0] * velo[1] - wv[1] * velo[0]).powi(2));
            } else {
                rotterm += sdet * (wv[2] * wv[2] * (velo[0] * velo[0] + velo[1] * velo[1]));
            }
        }

        let mut grad = [[0.0 as c_double; 3]; 3];
        unsafe {
            sta_lm_grad(
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                dbasis.as_ptr(),
                grad.as_mut_ptr() as *mut c_double,
            );
        }

        let (mut mu, mut rho) = (0.0, 0.0);
        unsafe {
            sta_lm_murho(
                n,
                basis.as_ptr(),
                state.mu.as_ptr(),
                state.rho.as_ptr(),
                &mut mu,
                &mut rho,
            );
        }
        let mut load_at_ip = [0.0 as c_double; 4];
        unsafe {
            sta_lm_loadatip(
                n,
                d as c_int,
                state.p as c_int,
                state.load.as_ptr(),
                basis.as_ptr(),
                load_at_ip.as_mut_ptr(),
            );
        }

        // at(i,j) aliases STIFF(i0+a, j0+b); db(a,b) is dBasisdx(a+1,b+1).
        let db = |a: usize, b: usize| dbasis[a + b * ndn];
        for p in 0..ndn {
            for q in 0..ndn {
                let i0 = (d + 1) * p;
                let j0 = (d + 1) * q;
                for ii in 0..d {
                    for jj in 0..d {
                        state.stiff[(i0 + ii) + (j0 + ii) * ld] +=
                            sdet * mu * db(q, jj) * db(p, jj);
                        state.stiff[(i0 + ii) + (j0 + jj) * ld] +=
                            sdet * mu * db(q, ii) * db(p, jj);
                    }
                    state.mass[(i0 + ii) + (j0 + ii) * ld] +=
                        sdet * rho * basis[p] * basis[q];
                    if stabilization {
                        state.stiff[(i0 + ii) + (j0 + d) * ld] -=
                            sdet * basis[q] * db(p, ii);
                        if p < nn {
                            state.stiff[(i0 + d) + (j0 + ii) * ld] -=
                                sdet * db(q, ii) * basis[p];
                        }
                    } else {
                        if q < nn {
                            state.stiff[(i0 + ii) + (j0 + d) * ld] -=
                                sdet * basis[q] * db(p, ii);
                        }
                        if p < nn {
                            state.stiff[(i0 + d) + (j0 + ii) * ld] -=
                                sdet * db(q, ii) * basis[p];
                        }
                    }
                }

                if convect {
                    let r1 = (d + 1) * p;
                    let c1 = (d + 1) * q;
                    if d > 2 {
                        state.stiff[r1 + c1 * ld] += sdet * rho * velo[1] * db(q, 1) * basis[p]
                            + sdet * rho * velo[2] * db(q, 2) * basis[p];
                        state.stiff[r1 + (c1 + 1) * ld] -=
                            sdet * rho * velo[1] * db(q, 0) * basis[p];
                        state.stiff[r1 + (c1 + 2) * ld] -=
                            sdet * rho * velo[2] * db(q, 0) * basis[p];
                        state.stiff[(r1 + 1) + c1 * ld] -=
                            sdet * rho * velo[0] * db(q, 1) * basis[p];
                        state.stiff[(r1 + 1) + (c1 + 1) * ld] += sdet * rho * velo[0]
                            * db(q, 0)
                            * basis[p]
                            + sdet * rho * velo[2] * db(q, 2) * basis[p];
                        state.stiff[(r1 + 1) + (c1 + 2) * ld] -=
                            sdet * rho * velo[2] * db(q, 1) * basis[p];
                        state.stiff[(r1 + 2) + c1 * ld] -=
                            sdet * rho * velo[0] * db(q, 2) * basis[p];
                        state.stiff[(r1 + 2) + (c1 + 1) * ld] -=
                            sdet * rho * velo[1] * db(q, 2) * basis[p];
                        state.stiff[(r1 + 2) + (c1 + 2) * ld] += sdet * rho * velo[1]
                            * db(q, 1)
                            * basis[p]
                            + sdet * rho * velo[0] * db(q, 0) * basis[p];
                        if newton {
                            state.stiff[r1 + (c1 + 1) * ld] -=
                                sdet * rho * wv[2] * basis[q] * basis[p];
                            state.stiff[r1 + (c1 + 2) * ld] +=
                                sdet * rho * wv[1] * basis[q] * basis[p];
                            state.stiff[(r1 + 1) + c1 * ld] +=
                                sdet * rho * wv[2] * basis[q] * basis[p];
                            state.stiff[(r1 + 1) + (c1 + 2) * ld] -=
                                sdet * rho * wv[0] * basis[q] * basis[p];
                            state.stiff[(r1 + 2) + c1 * ld] -=
                                sdet * rho * wv[1] * basis[q] * basis[p];
                            state.stiff[(r1 + 2) + (c1 + 1) * ld] +=
                                sdet * rho * wv[0] * basis[q] * basis[p];
                        }
                    } else {
                        state.stiff[r1 + c1 * ld] +=
                            sdet * rho * velo[1] * db(q, 1) * basis[p];
                        state.stiff[r1 + (c1 + 1) * ld] -=
                            sdet * rho * velo[1] * db(q, 0) * basis[p];
                        state.stiff[(r1 + 1) + c1 * ld] -=
                            sdet * rho * velo[0] * db(q, 1) * basis[p];
                        state.stiff[(r1 + 1) + (c1 + 1) * ld] +=
                            sdet * rho * velo[0] * db(q, 0) * basis[p];
                        if newton {
                            state.stiff[r1 + (c1 + 1) * ld] -=
                                sdet * rho * wv[2] * basis[q] * basis[p];
                            state.stiff[(r1 + 1) + c1 * ld] +=
                                sdet * rho * wv[2] * basis[q] * basis[p];
                        }
                        if hybrid {
                            state.stiff[r1 + c1 * ld] -=
                                rho * sdet * grad[0][0] * basis[q] * basis[p];
                            state.stiff[r1 + (c1 + 1) * ld] -=
                                rho * sdet * grad[0][1] * basis[q] * basis[p];
                            state.stiff[(r1 + 1) + c1 * ld] -=
                                rho * sdet * grad[1][0] * basis[q] * basis[p];
                            state.stiff[(r1 + 1) + (c1 + 1) * ld] -=
                                rho * sdet * grad[1][1] * basis[q] * basis[p];
                        }
                    }
                }
            }

            let i0 = (d + 1) * p;
            for k in 0..=d {
                state.force[i0 + k] += sdet * load_at_ip[k] * basis[p];
            }
            if newton {
                state.force[i0] -= sdet * rho * basis[p] * velo[1] * wv[2];
                state.force[i0 + 1] += sdet * rho * basis[p] * velo[0] * wv[2];
                if hybrid && d == 2 {
                    state.force[i0] -= sdet
                        * rho
                        * basis[p]
                        * (grad[0][0] * velo[0] + grad[0][1] * velo[1]);
                    state.force[i0 + 1] -= sdet
                        * rho
                        * basis[p]
                        * (grad[1][0] * velo[0] + grad[1][1] * velo[1]);
                }
                if d > 2 {
                    state.force[i0] += sdet * rho * basis[p] * velo[2] * wv[1];
                    state.force[i0 + 1] -= sdet * rho * basis[p] * velo[2] * wv[0];
                    state.force[i0 + 2] -= sdet * rho * basis[p] * velo[0] * wv[1];
                    state.force[i0 + 2] += sdet * rho * basis[p] * velo[1] * wv[0];
                }
            }
        }
    }

    // Bubble nullification.
    for p in nn..ndn {
        let i1 = (d + 1) * (p + 1);
        state.force[i1 - 1] = 0.0;
        for c in 0..ld {
            state.stiff[(i1 - 1) + c * ld] = 0.0;
        }
        if !stabilization {
            state.stiff[(i1 - 1) + (i1 - 1) * ld] = 1.0;
        }
    }

    if stabilization {
        second_pass(state, n, nd, dim, ak, rotterm, ch);
    }
}
/// Convection-stabilization second sweep of `LocalMatrix` (only reachable
/// with stabilization enabled, which the entry unconditionally disables;
/// transcribed for fidelity).
fn second_pass(
    state: &mut StokesState,
    n: c_int,
    nd: c_int,
    dim: usize,
    ak: c_double,
    rotterm: c_double,
    ch0: c_double,
) {
    let ld = state.n;
    let nn = n as usize;
    let ndn = nd as usize;
    let d = dim;
    let rotterm = rotterm.sqrt() / ak.sqrt();
    let mut ch: c_double;
    let _ = ch0 / 3.0 * rotterm.sqrt();
    let (mut hk, mut mk) = (0.0, 0.0);
    unsafe { sta_lm_hk_mk(&mut hk, &mut mk) };
    let mut basis = vec![0.0 as c_double; ndn];
    let mut ddbasis = vec![0.0 as c_double; ndn * 9];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];
    // Second integration loop re-establishes the Gauss points and element
    // nodes, exactly as the original's reassigned `IP` local does (`ch`
    // is recomputed identically and discarded).
    let (mut nip, mut dummy_ch) = (0, 0.0);
    unsafe { sta_ip_begin(&mut nip, &mut dummy_ch) };
    let _ = dummy_ch;
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { sta_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            sta_lm_eleminfo2(
                u, v, w, n, nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                ddbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;
        let sdet = s * detj;
        let mut velo = [0.0 as c_double; 3];
        unsafe {
            sta_lm_velo(
                d as c_int,
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                basis.as_ptr(),
                velo.as_mut_ptr(),
            );
        }
        let mut wv = [0.0 as c_double; 3];
        unsafe {
            sta_lm_w(
                0,
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                dbasis.as_ptr(),
                wv.as_mut_ptr(),
            );
        }
        let db = |a: usize, b: usize| dbasis[a + b * ndn];
        // Exact MAXVALs (no rounding): a1/a2 over velocity rows 1..2.
        let mut a1 = 0.0 as c_double;
        let mut a2 = 0.0 as c_double;
        for k in 0..ndn {
            a1 = a1.max(state.velocity[k * (d + 1)].abs());
            a2 = a2.max(state.velocity[1 + k * (d + 1)].abs());
        }
        let vnorm = a1.max(a2);
        let (mut mu, mut rho) = (0.0, 0.0);
        unsafe {
            sta_lm_murho(
                n,
                basis.as_ptr(),
                state.mu.as_ptr(),
                state.rho.as_ptr(),
                &mut mu,
                &mut rho,
            );
        }
        let re = (rho * mk * hk * vnorm / (4.0 * mu)).min(1.0);
        ch = (3.0 as c_double / 100.0).sqrt() * (hk * re * rho / vnorm).sqrt();
        // Pressure-bubble construction (robust Picard branch).
        for p in nn..ndn {
            for q in nn..ndn {
                let mut sm = 0.0;
                for k in 0..d {
                    sm += db(p, k) * db(q, k);
                }
                state.stiff[(d + 1) * (p + 1) - 1 + ((d + 1) * (q + 1) - 1) * ld] +=
                    sdet * sm;
            }
            let c = -velo[1] * db(p, 0) + velo[0] * db(p, 1);
            for q in 0..ndn {
                state.stiff[(d + 1) * (p + 1) - 1 + ((d + 1) * q + 1) * ld] +=
                    rho * sdet * db(q, 0) * c;
                state.stiff[(d + 1) * (p + 1) - 1 + ((d + 1) * q) * ld] -=
                    rho * sdet * db(q, 1) * c;
            }
        }
        // Projection part of the equation of motion (robust branch).
        for p in 0..ndn {
            for q in nn..ndn {
                state.stiff[(d + 1) * p + ((d + 1) * (q + 1) - 1) * ld] +=
                    ch * ch * sdet * db(q, 0) * db(p, 1) * velo[1];
                state.stiff[(d + 1) * p + ((d + 1) * (q + 1) - 1) * ld] -=
                    ch * ch * sdet * db(q, 1) * db(p, 1) * velo[0];
                state.stiff[(d + 1) * p + 1 + ((d + 1) * (q + 1) - 1) * ld] -=
                    ch * ch * sdet * db(q, 0) * db(p, 0) * velo[1];
                state.stiff[(d + 1) * p + 1 + ((d + 1) * (q + 1) - 1) * ld] +=
                    ch * ch * sdet * db(q, 1) * db(p, 0) * velo[0];
            }
            let c2 = rho * ch * ch * db(p, 1) * (velo[0] * velo[0] + velo[1] * velo[1]);
            for q in 0..ndn {
                state.stiff[(d + 1) * p + ((d + 1) * q) * ld] += sdet * c2 * db(q, 1);
                state.stiff[(d + 1) * p + ((d + 1) * q + 1) * ld] -= sdet * c2 * db(q, 0);
            }
            let c2 = rho * ch * ch * db(p, 0) * (velo[0] * velo[0] + velo[1] * velo[1]);
            for q in 0..ndn {
                state.stiff[(d + 1) * p + 1 + ((d + 1) * q) * ld] -= sdet * c2 * db(q, 1);
                state.stiff[(d + 1) * p + 1 + ((d + 1) * q + 1) * ld] +=
                    sdet * c2 * db(q, 0);
            }
        }
    }
}
/// The `LocalMatrixBoundary` internal procedure: outflow-boundary assembly.
/// `n`/`nd` are the element node/dof counts.
fn local_matrix_boundary(
    state: &mut StokesState,
    n: c_int,
    nd: c_int,
    dim: usize,
    newton: bool,
    block: bool,
) {
    let ld = state.n;
    let m = state.m;
    let ndn = nd as usize;
    let d = dim;
    state.stiff.fill(0.0);
    state.force.fill(0.0);
    state.alocal.fill(0.0);
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];
    let mut nip = 0;
    unsafe { sta_schur_begin(&mut nip) };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { sta_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            sta_lm_eleminfo(
                u, v, w, n, nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;
        let sdet = s * detj;
        let mut nl = [0.0 as c_double; 3];
        unsafe { sta_lb_normal(u, v, nl.as_mut_ptr()) };
        let mut velo = [0.0 as c_double; 3];
        unsafe {
            sta_lm_velo(
                d as c_int,
                (d + 1) as c_int,
                nd,
                state.velocity.as_ptr(),
                basis.as_ptr(),
                velo.as_mut_ptr(),
            );
        }
        let mut sq: c_double = velo.iter().take(d).map(|v| v * v).sum();
        let rho1 = state.rho[0];
        sq *= rho1;
        let c = if newton { 1.0 } else { 0.5 };
        for p in 0..ndn {
            for (i, nl_i) in nl.iter().take(d).enumerate() {
                let fr = p * (d + 1) + i;
                if newton {
                    state.force[fr] += sdet * nl_i * basis[p] * 0.5 * sq;
                }
                for q in 0..ndn {
                    let c1 = q * (d + 1);
                    state.stiff[fr + c1 * ld] +=
                        c * sdet * basis[p] * nl_i * velo[0] * basis[q] * rho1;
                    state.stiff[fr + (c1 + 1) * ld] +=
                        c * sdet * basis[p] * nl_i * velo[1] * basis[q] * rho1;
                    if d > 2 {
                        state.stiff[fr + (c1 + 2) * ld] +=
                            c * sdet * basis[p] * nl_i * velo[2] * basis[q] * rho1;
                    }
                }
            }
        }
        if block {
            for p in 0..ndn {
                for q in 0..ndn {
                    for i in 0..d {
                        state.alocal[p * d + i + (q * d) * m] +=
                            c * sdet * basis[p] * nl[i] * velo[0] * basis[q] * rho1;
                        state.alocal[p * d + i + (q * d + 1) * m] +=
                            c * sdet * basis[p] * nl[i] * velo[1] * basis[q] * rho1;
                        if d > 2 {
                            state.alocal[p * d + i + (q * d + 2) * m] +=
                                c * sdet * basis[p] * nl[i] * velo[2] * basis[q] * rho1;
                        }
                    }
                }
            }
        }
    }
}

/// The `SchurComplementMatrix` internal procedure: pressure mass and
/// Laplacian blocks for the block preconditioner.
fn schur_complement(state: &mut StokesState, n: c_int, nd: c_int, dim: usize) {
    let ldp = state.p;
    let nn = n as usize;
    let ndn = nd as usize;
    let d = dim;
    state.mlocal.fill(0.0);
    state.plocal.fill(0.0);
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];
    let mut ddbasis = vec![0.0 as c_double; ndn * 9];
    let mut nip = 0;
    unsafe { sta_schur_begin(&mut nip) };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { sta_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            sta_lm_eleminfo2(
                u, v, w, n, nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                ddbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            );
        }
        let _ = stat;
        let sdet = s * detj;
        let db = |a: usize, b: usize| dbasis[a + b * ndn];
        let (mut mu, mut rho) = (0.0, 0.0);
        unsafe {
            sta_lm_murho(
                n,
                basis.as_ptr(),
                state.mu.as_ptr(),
                state.rho.as_ptr(),
                &mut mu,
                &mut rho,
            );
        }
        let _ = mu;
        for p in 0..nn {
            for q in 0..nn {
                state.mlocal[p + q * ldp] += sdet * basis[q] * basis[p];
                let mut sm = 0.0;
                for k in 0..d {
                    sm += db(q, k) * db(p, k);
                }
                state.plocal[p + q * ldp] += sdet * sm;
            }
        }
    }
    for p in nn..ndn {
        state.mlocal[p + p * ldp] = 1.0;
        state.plocal[p + p * ldp] = 1.0;
    }
}
/// The `GCROuterIteration` internal procedure (block path only): GCR(m)
/// iteration with block-factorization preconditioning.
#[allow(clippy::too_many_arguments)]
fn gcr_outer(
    state: &mut StokesState,
    rounds: c_int,
    tol: c_double,
    dim: usize,
    l: c_int,
    q: c_int,
    mu1: c_double,
    dt_val: c_double,
    transient_simulation: bool,
) {
    let n = unsafe { sta_num_rows() };
    let nn = n as usize;
    let nr = rounds as usize;
    let d = dim;
    let (inner_tol, _) = list_const_real(b"Linear System Convergence Tolerance", false, 0.0, false);
    let (inner_rounds, _) = list_integer(b"Linear System Max Iterations", false, 0, false);
    let mut found = 0;
    let _ = unsafe { sta_get_max_restarts(&mut found) };
    let mut x = vec![0.0 as c_double; nn];
    let mut b = vec![0.0 as c_double; nn];
    unsafe {
        sta_values_copy(x.as_mut_ptr(), n);
        sta_rhs_copy(b.as_mut_ptr(), n);
    }
    let mut r = vec![0.0 as c_double; nn];
    unsafe { sta_pmv_main(n, x.as_ptr(), r.as_mut_ptr()) };
    for i in 0..nn {
        r[i] = b[i] - r[i];
    }
    let rhs_norm = unsafe { sta_rhs_pnorm() };
    let mut res = unsafe { sta_pnorm_vec(n, r.as_ptr()) } / rhs_norm;
    unsafe { sta_gcr_info0(res) };
    let mut converged = res < tol;
    let mut v = vec![0.0 as c_double; nn * nr.max(1)];
    let mut s = vec![0.0 as c_double; nn * nr.max(1)];
    if !converged {
        let mut t1 = vec![0.0 as c_double; nn];
        let mut t2 = vec![0.0 as c_double; nn];
        for k in 1..=rounds as usize {
            t1.copy_from_slice(&r);
            preconditioning(
                state, n, &mut x, l, q, d, inner_tol, inner_rounds, mu1, dt_val,
                transient_simulation, &mut t1,
            );
            unsafe { sta_pmv_main(n, t1.as_ptr(), t2.as_mut_ptr()) };
            for i in 1..k {
                let beta = unsafe {
                    sta_pdot(
                        n,
                        v.as_ptr().wrapping_add((i - 1) * nn),
                        t2.as_ptr(),
                    )
                };
                for j in 0..nn {
                    t1[j] -= beta * s[(i - 1) * nn + j];
                    t2[j] -= beta * v[(i - 1) * nn + j];
                }
            }
            let alpha = unsafe { sta_pnorm_vec(n, t2.as_ptr()) };
            for j in 0..nn {
                t1[j] /= alpha;
                t2[j] /= alpha;
            }
            let beta = unsafe { sta_pdot(n, t2.as_ptr(), r.as_ptr()) };
            for j in 0..nn {
                x[j] += beta * t1[j];
                r[j] -= beta * t2[j];
            }
            s[(k - 1) * nn..k * nn].copy_from_slice(&t1);
            v[(k - 1) * nn..k * nn].copy_from_slice(&t2);
            res = unsafe { sta_pnorm_vec(n, r.as_ptr()) } / rhs_norm;
            unsafe { sta_gcr_info(k as c_int, res) };
            converged = res < tol;
            if converged {
                break;
            }
            if k == 1 {
                unsafe {
                    sta_precomp_toggle(0, 1);
                    sta_precomp_toggle(1, 1);
                    sta_precomp_toggle(2, 1);
                }
            }
        }
        unsafe {
            sta_precomp_toggle(0, 0);
            sta_precomp_toggle(1, 0);
            sta_precomp_toggle(2, 0);
        }
        let _ = converged;
    }
    unsafe { sta_values_set(x.as_ptr(), n) };
}

/// The `PreconditioningIteration` internal procedure (block path only).
#[allow(clippy::too_many_arguments)]
fn preconditioning(
    _state: &mut StokesState,
    n: c_int,
    x: &mut [c_double],
    l: c_int,
    q: c_int,
    dim: usize,
    tol: c_double,
    rounds: c_int,
    mu1: c_double,
    dt_val: c_double,
    transient_simulation: bool,
    r: &mut [c_double],
) {
    let nn = n as usize;
    let d = dim;
    let mut z = vec![0.0 as c_double; nn];
    z.copy_from_slice(x);
    let nl = l as usize;
    let nq = (q / 2) as usize;
    let mut f = vec![0.0 as c_double; nl];
    let mut y = vec![0.0 as c_double; nl];
    for j in 0..nn / (d + 1) {
        for m in 0..d {
            f[d * j + m] = r[(d + 1) * j + m];
        }
    }
    y.fill(0.0);
    unsafe { sta_print_msg(0) };
    let mut inorm = 0.0;
    unsafe { sta_solve_lin(0, l, f.as_ptr(), y.as_mut_ptr(), d as c_int, &mut inorm) };
    for j in 0..nn / (d + 1) {
        for m in 0..d {
            z[(d + 1) * j + m] += y[d * j + m];
        }
    }
    let mut g = vec![0.0 as c_double; nq];
    unsafe { sta_crs_divmult(d as c_int, z.as_ptr(), g.as_mut_ptr()) };
    let mut fvec = vec![0.0 as c_double; nq];
    let mut prevpres = vec![0.0 as c_double; nq];
    for j in 0..nq {
        fvec[j] = 2.0 * mu1 * g[j];
        prevpres[j] = x[(d + 1) * j + d];
    }
    let mut dvec = vec![0.0 as c_double; nq];
    unsafe { sta_print_msg(1) };
    unsafe {
        sta_solve_lin(1, q / 2, fvec.as_ptr(), dvec.as_mut_ptr(), 1, &mut inorm);
    }
    let mut tmpvec = vec![0.0 as c_double; nq];
    unsafe { sta_crs_matvec_which(3, q / 2, dvec.as_ptr(), tmpvec.as_mut_ptr()) };
    if transient_simulation {
        for j in 0..nq {
            fvec[j] = -tmpvec[j] - 1.0 / dt_val * g[j];
        }
    } else {
        for j in 0..nq {
            fvec[j] = -tmpvec[j];
        }
    }
    unsafe { sta_prec_traction(fvec.as_mut_ptr(), dvec.as_ptr()) };
    dvec.fill(0.0);
    unsafe { sta_print_msg(2) };
    unsafe {
        sta_solve_lin(2, q / 2, fvec.as_ptr(), dvec.as_mut_ptr(), 1, &mut inorm);
    }
    for j in 0..nn / (d + 1) {
        z[(d + 1) * j + d] = prevpres[j] + dvec[j];
    }
    for i in 0..nn {
        r[i] = z[i] - x[i];
    }
    let _ = (tol, rounds);
}

/// The internal ILU-preconditioned BiCGStab method. All of its call sites
/// in the original are commented out, so this is unreachable; it is
/// transcribed for fidelity.
#[allow(dead_code)]
fn bicgstab(n: c_int, which_mat: c_int, x: &mut [c_double], b: &[c_double], rounds: c_int, tol: c_double) {
    let nn = n as usize;
    let all_zero = x.iter().all(|&v| v == 0.0);
    if all_zero {
        x.copy_from_slice(b);
    }
    let mut r = vec![0.0 as c_double; nn];
    unsafe { sta_crs_matvec_which(which_mat, n, x.as_ptr(), r.as_mut_ptr()) };
    for i in 0..nn {
        r[i] = b[i] - r[i];
    }
    let ri = r.clone();
    let mut p = vec![0.0 as c_double; nn];
    let mut vv = vec![0.0 as c_double; nn];
    let mut omega = 1.0;
    let mut alpha = 0.0;
    let mut oldrho = 1.0;
    for i in 1..=rounds {
        let rho = unsafe { sta_pdot(n, r.as_ptr(), ri.as_ptr()) };
        let beta = alpha * rho / (oldrho * omega);
        for j in 0..nn {
            p[j] = r[j] + beta * (p[j] - omega * vv[j]);
        }
        vv.copy_from_slice(&p);
        unsafe { sta_crs_lusolve(which_mat, n, vv.as_mut_ptr()) };
        let t1 = vv.clone();
        unsafe { sta_crs_matvec_which(which_mat, n, t1.as_ptr(), vv.as_mut_ptr()) };
        let v1 = vv.clone();
        alpha = rho / unsafe { sta_pdot(n, ri.as_ptr(), vv.as_ptr()) };
        let mut s = vec![0.0 as c_double; nn];
        for j in 0..nn {
            s[j] = r[j] - alpha * vv[j];
        }
        let mut t = s.clone();
        unsafe { sta_crs_lusolve(which_mat, n, t.as_mut_ptr()) };
        let t2 = t.clone();
        unsafe { sta_crs_matvec_which(which_mat, n, t2.as_ptr(), t.as_mut_ptr()) };
        omega = unsafe { sta_pdot(n, t.as_ptr(), s.as_ptr()) }
            / unsafe { sta_pdot(n, t.as_ptr(), t.as_ptr()) };
        oldrho = rho;
        for j in 0..nn {
            r[j] = s[j] - omega * t[j];
            x[j] += alpha * t1[j] + omega * t2[j];
        }
        let res = unsafe { sta_pnorm_vec(n, r.as_ptr()) }
            / unsafe { sta_pnorm_vec(n, b.as_ptr()) };
        unsafe { sta_bicg_write(i, res) };
        if res < tol {
            break;
        }
        let _ = v1;
    }
}



