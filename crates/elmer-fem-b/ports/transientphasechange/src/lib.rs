//! Port of `fem/src/modules/TransientPhaseChange.F90` (Elmer @ 9f6af2f85,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/TransientPhaseChange.F90` — two F77-global
//! entries `transientphasechange_` (main) and `transientphasechange_init_`
//! (init); three `CONTAINS` internals. The emitted dynamic symbols are
//! exactly `transientphasechange_` + `transientphasechange_init_` (proven
//! by `nm -D` on the reference `.so`); the internals emit only local
//! symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/transientphasechange_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MINVAL`/`MAXVAL`), formatted I/O (`WRITE ... ES12.5`, `I2S`),
//! whole-array updates, and the three `CONTAINS` procedures (which mix
//! those with deep derived-type graph walks over `BoundaryInfo`,
//! parent elements and `Model` tables) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors TransientPhaseChange.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: latched scalars (`SubroutineVisited`,
//!   `Visited`, `prevpos0`, `pos0`, `UPull`, `Trip_node`,
//!   `PullVelocitySet`, `DoVelocityRelax`) live in the Rust `Mutex`
//!   state; latched arrays/pointers (`Nodes`, `NodalTemp`,
//!   `Conductivity`, `LatentHeat`, `Density`, local matrices,
//!   `IsBoundaryNode`, `PrevSurfaceVelo`, `CurrentLoads`, `PrevLoads`,
//!   `SurfaceVelo`, `NormalDir`, `AverageNormal`, `NrmSol`) live as
//!   wrapper-side SAVE state with identical lifetimes (consumed before
//!   the next fill, as in the original). Per-call associations
//!   (`Surface`, `PrevSurface`, `Temperature`, `ForceVector`,
//!   `StiffMatrix`, ...) are refreshed on every entry, exactly where the
//!   original re-associates its locals. `PullControl` is SAVE-listed but
//!   re-read before any read each call, so a per-call local is
//!   equivalent.
//! - Dead state omitted with zero observable effect (verified: assigned
//!   but never read, or declared but never referenced): `PSolver`,
//!   `AllocationsDone`, `RelativeChange`, `NonlinearTol`, `PrevUpull`,
//!   host scalars `d`/`s`/`u`/`v`/`w`/`yy`/`ds`, and the
//!   never-associated `Normals` pointer. The `NonlinearTol`
//!   `ListGetConstReal` call itself is kept (same call sequence).
//! - Allocation failure calls the same `Fatal(...)` via the same
//!   allocation sites (no `STAT` check in the original loads allocs, so
//!   the wrappers allocate the same way), instead of aborting.
//! - `GOTO 200` (Passive Steps) is a branch to the shared pull-position
//!   tail; `FirstTime = .FALSE.` runs on all paths, as in the original.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn tpa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn tpa_info_header();
    fn tpa_fatal_notransient();
    fn tpa_coord_dim() -> c_int;
    fn tpa_bind_primary();
    fn tpa_bind_temp();
    fn tpa_nonlin_iter() -> c_int;
    fn tpa_nonlintol_touch();
    fn tpa_pullcontrol(dim: c_int) -> c_int;
    fn tpa_triplefixed(dim: c_int) -> c_int;
    fn tpa_coupled_iter() -> c_int;
    fn tpa_timestep() -> c_int;
    fn tpa_nof_nodes() -> c_int;
    fn tpa_num_bnd() -> c_int;
    fn tpa_surfperm(node: c_int) -> c_int;
    fn tpa_nodecoord(node: c_int, j1: c_int) -> c_double;
    fn tpa_normaldir(jdet: c_int) -> c_int;
    fn tpa_alloc_bndflags();
    fn tpa_set_bndflag(idx: c_int);
    fn tpa_bnd_scan(t: c_int, n: *mut c_int, has: *mut c_int);
    fn tpa_bnd_nodeidx(i: c_int) -> c_int;
    fn tpa_report_bndcount();
    fn tpa_alloc_first();
    fn tpa_bind_normals(nliter: c_int);
    fn tpa_bind_velo();
    fn tpa_velrelax(val: *mut c_double, found: *mut c_int);
    fn tpa_alloc_prevvelo();
    fn tpa_passive_steps(val: *mut c_int, found: *mut c_int);
    fn tpa_pull_velo() -> c_double;
    fn tpa_use_loads() -> c_int;
    fn tpa_bind_loads(sign: *mut c_int, usefirst: *mut c_int);
    fn tpa_loads_snapshot();
    fn tpa_loadsperm(node: c_int) -> c_int;
    fn tpa_loadsval(j: c_int) -> c_double;
    fn tpa_curloads_set(i: c_int, v: c_double);
    fn tpa_report_loads();
    fn tpa_snap_prevvelo();
    fn tpa_assoc_mat();
    fn tpa_loadsrelax(val: *mut c_double, useavg: *mut c_int);
    fn tpa_fill_force(order: c_int, relax: c_double);
    fn tpa_default_init();
    fn tpa_nactive() -> c_int;
    fn tpa_ebegin(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn tpa_egather(n: c_int);
    fn tpa_evelo(n: c_int, nd: c_int, ul: c_int);
    fn tpa_eupdate();
    fn tpa_finish_asm();
    fn tpa_dirichlet();
    fn tpa_solve_velo() -> c_double;
    fn tpa_nonlin(change: *mut c_double, conv: *mut c_int);
    fn tpa_velrelax2() -> c_double;
    fn tpa_relax_velo(vr: c_double);
    fn tpa_report_velo();
    fn tpa_trip_velo(trip: c_int) -> c_double;
    fn tpa_info_pull(upull: c_double);
    fn tpa_speedup(val: *mut c_double, found: *mut c_int);
    fn tpa_disprelax(val: *mut c_double, stat: *mut c_int);
    fn tpa_solver_order() -> c_int;
    fn tpa_simple_update(spdup: c_double, dtv: c_double, upullv: c_double);
    fn tpa_pde_info();
    fn tpa_sbegin(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn tpa_sgather();
    #[allow(clippy::too_many_arguments)]
    fn tpa_surf_compute(nc: c_int, spdup: c_double, upullv: c_double, dtv: c_double);
    fn tpa_first_order();
    fn tpa_solve_surf() -> c_double;
    fn tpa_donetime() -> c_int;
    fn tpa_find_pull(trip: c_int);
    fn tpa_res_pull(pos0v: c_double, upullv: c_double, dopull: c_int);
    fn tpa_init_begin(h_model: *const c_void, h_solver: *const c_void);
    fn tpa_init_body();
}

/// Latched SAVE scalars (mirrors the original SAVE statement lifetimes).
/// Per-call values (`PullControl`, `TriplePointFixed`, `UseLoads`,
/// `DIM`, ...) are locals, exactly where the original re-reads them.
struct TpcState {
    first_time: bool,
    subvisited: c_int,
    visited: c_int,
    prevpos0: c_double,
    pos0: c_double,
    upull: c_double,
    trip_node: c_int,
    pull_vel_set: bool,
    do_velo_relax: bool,
}

impl TpcState {
    const fn new() -> Self {
        Self {
            first_time: true,
            subvisited: 0,
            visited: -1,
            prevpos0: 0.0,
            pos0: 0.0,
            upull: 0.0,
            trip_node: 0,
            pull_vel_set: false,
            do_velo_relax: false,
        }
    }
}

static STATE: Mutex<TpcState> = Mutex::new(TpcState::new());

/// `transientphasechange_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes; `dt`/`transient` must be valid reads.
#[export_name = "transientphasechange_"]
pub unsafe extern "C" fn transient_phase_change_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    let dtv: c_double = unsafe { *dt };
    unsafe {
        tpa_begin(model as *const c_void, solver as *const c_void);
        tpa_info_header();
        if *transient == 0 {
            tpa_fatal_notransient();
            std::process::abort();
        }
    }

    state.subvisited += 1;
    let dim = unsafe { tpa_coord_dim() };
    unsafe {
        tpa_bind_primary();
        tpa_bind_temp();
    }
    let nliter = unsafe { tpa_nonlin_iter() };
    unsafe { tpa_nonlintol_touch() };
    let pull_control = unsafe { tpa_pullcontrol(dim) } != 0;
    let triple_fixed = unsafe { tpa_triplefixed(dim) } != 0;
    let coupled = unsafe { tpa_coupled_iter() };
    let timestep = unsafe { tpa_timestep() };

    // First-time permanent setup (latched SAVE state).
    if state.first_time {
        let nnox = unsafe { tpa_nof_nodes() };
        let mut cmax = [f64::MIN; 3];
        let mut cmin = [f64::MAX; 3];
        let mut cmaxi = [0 as c_int; 3];
        let mut cmini = [0 as c_int; 3];
        for k in 1..=nnox {
            if unsafe { tpa_surfperm(k) } <= 0 {
                continue;
            }
            for j in 0..(dim as usize) {
                let xx = unsafe { tpa_nodecoord(k, j as c_int + 1) };
                if xx > cmax[j] {
                    cmax[j] = xx;
                    cmaxi[j] = k;
                }
                if xx < cmin[j] {
                    cmin[j] = xx;
                    cmini[j] = k;
                }
            }
        }
        // Direction of minimum change is the normal direction if not given.
        let mut jdet = 0_usize;
        for i in 0..(dim as usize) {
            if cmax[i] - cmin[i] < cmax[jdet] - cmin[jdet] {
                jdet = i;
            }
        }
        let normal_dir = unsafe { tpa_normaldir(jdet as c_int + 1) };
        unsafe { tpa_alloc_bndflags() };
        if dim == 2 {
            // 1-based Fortran index math: CoordMaxi(3 - NormalDir).
            let ti = (3 - normal_dir) as usize - 1;
            state.trip_node = cmaxi[ti];
            let axis_node = cmini[ti];
            unsafe {
                tpa_set_bndflag(tpa_surfperm(state.trip_node));
                tpa_set_bndflag(tpa_surfperm(axis_node));
            }
        }
        let nbel = unsafe { tpa_num_bnd() };
        for t in 1..=nbel {
            let (mut n, mut has) = (0, 0);
            unsafe { tpa_bnd_scan(t, &mut n, &mut has) };
            if has == 0 {
                continue;
            }
            for i in 1..=n {
                let node = unsafe { tpa_bnd_nodeidx(i) };
                let sp = unsafe { tpa_surfperm(node) };
                if sp > 0 {
                    unsafe { tpa_set_bndflag(sp) };
                }
            }
        }
        unsafe {
            tpa_report_bndcount();
            tpa_alloc_first();
            tpa_bind_normals(nliter);
            tpa_bind_velo();
        }
        let (mut vr, mut found) = (0.0, 0);
        unsafe { tpa_velrelax(&mut vr, &mut found) };
        if found != 0 {
            unsafe { tpa_alloc_prevvelo() };
            state.do_velo_relax = true;
        }
    }

    let (mut passive, mut pfound) = (0, 0);
    unsafe { tpa_passive_steps(&mut passive, &mut pfound) };
    let skip = pfound != 0 && passive >= state.subvisited;
    if !skip {
        if !state.pull_vel_set {
            state.upull = unsafe { tpa_pull_velo() };
            state.pull_vel_set = true;
        }

        let use_loads = unsafe { tpa_use_loads() } != 0;
        let (mut loads_sign, mut use_first) = (0, 0);
        if use_loads {
            unsafe { tpa_bind_loads(&mut loads_sign, &mut use_first) };
            if use_first == 0 && timestep > 1 && coupled == 1 {
                unsafe { tpa_loads_snapshot() };
            }
            let nnox = unsafe { tpa_nof_nodes() };
            for t in 1..=nnox {
                let i = unsafe { tpa_surfperm(t) };
                if i == 0 {
                    continue;
                }
                let j = unsafe { tpa_loadsperm(t) };
                if j == 0 {
                    continue;
                }
                let lv = unsafe { tpa_loadsval(j) };
                unsafe { tpa_curloads_set(i, loads_sign as c_double * lv) };
            }
            if use_first != 0 && coupled == 1 {
                unsafe { tpa_loads_snapshot() };
            }
            unsafe { tpa_report_loads() };
        }

        if state.do_velo_relax {
            unsafe { tpa_snap_prevvelo() };
        }
        unsafe { tpa_assoc_mat() };

        let (mut loads_relax, mut use_avg) = (0.0, 0);
        unsafe { tpa_loadsrelax(&mut loads_relax, &mut use_avg) };
        let mut loads_order: c_int = 0;
        if use_avg != 0 && (use_first == 0 || coupled > 1) {
            loads_order = std::cmp::min(2, timestep - 1);
        }

        // `DO iter = 1, NonlinearIter`: zero iterations when NonlinearIter < 1.
        for _iter in 1..=nliter {
            unsafe { tpa_default_init() };
            if use_loads {
                unsafe { tpa_fill_force(loads_order, loads_relax) };
            }
            let nactive = unsafe { tpa_nactive() };
            for t in 1..=nactive {
                let (mut n, mut nd) = (0, 0);
                unsafe {
                    tpa_ebegin(t, &mut n, &mut nd);
                    tpa_egather(n);
                    tpa_evelo(n, nd, use_loads as c_int);
                    tpa_eupdate();
                }
            }
            unsafe {
                tpa_finish_asm();
                tpa_dirichlet();
            }
            let _norm: c_double = unsafe { tpa_solve_velo() };
            let (mut _change, mut conv) = (0.0, 0);
            unsafe { tpa_nonlin(&mut _change, &mut conv) };
            if conv == 1 {
                break;
            }
        }

        if state.do_velo_relax {
            let vr = unsafe { tpa_velrelax2() };
            unsafe { tpa_relax_velo(vr) };
        }
        unsafe { tpa_report_velo() };

        if pull_control || triple_fixed {
            state.upull = unsafe { tpa_trip_velo(state.trip_node) };
            if pull_control {
                unsafe { tpa_info_pull(state.upull) };
            }
        }

        // Displacement update: direct or PDE-smoothed.
        let (mut speedup, mut sfound) = (0.0, 0);
        unsafe { tpa_speedup(&mut speedup, &mut sfound) };
        if sfound == 0 {
            speedup = 1.0;
        }
        let (mut disprelax, mut _dstat) = (0.0, 0);
        unsafe { tpa_disprelax(&mut disprelax, &mut _dstat) };
        let order = unsafe { tpa_solver_order() };
        if order <= 1 && disprelax.abs() < f64::EPSILON {
            unsafe { tpa_simple_update(speedup, dtv, state.upull) };
        } else {
            unsafe { tpa_pde_info() };
            unsafe { tpa_default_init() };
            let nactive = unsafe { tpa_nactive() };
            for t in 1..=nactive {
                let (mut n, mut nd) = (0, 0);
                unsafe {
                    tpa_sbegin(t, &mut n, &mut nd);
                    tpa_sgather();
                    tpa_surf_compute(n, speedup, state.upull, dtv);
                    tpa_first_order();
                    tpa_eupdate();
                }
            }
            unsafe {
                tpa_finish_asm();
                tpa_dirichlet();
            }
            let _norm: c_double = unsafe { tpa_solve_surf() };
        }

        // Pulled-distance integration.
        if pull_control || triple_fixed {
            let donetime = unsafe { tpa_donetime() };
            if state.visited != donetime {
                if donetime == 1 {
                    state.prevpos0 = 0.0;
                } else {
                    state.prevpos0 = state.pos0;
                }
                state.visited = donetime;
            }
            state.pos0 = state.prevpos0 + state.upull * dtv;
            if pull_control {
                unsafe { tpa_find_pull(state.trip_node) };
            }
        }
    }

    // 200 CONTINUE: pull position reporting (falls through on all paths).
    if pull_control || triple_fixed {
        unsafe { tpa_res_pull(state.pos0, state.upull, pull_control as c_int) };
    }

    state.first_time = false;
}

/// `transientphasechange_init_` — the solver init entry.
///
/// # Safety
///
/// Same handle contract as the main entry.
#[export_name = "transientphasechange_init_"]
pub unsafe extern "C" fn transient_phase_change_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        tpa_init_begin(model as *const c_void, solver as *const c_void);
        tpa_init_body();
    }
}
