//! Port of `fem/src/modules/PhaseChangeSolve.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/PhaseChangeSolve.F90` — two F77-global
//! entries, `phasechangesolve_` (the deprecated combined steady/transient
//! Lagrangian phase-change solver) and `phasechangesolve_init_`, plus the
//! `meltingheat_` co-entry (the `Heat Flux` real-procedure UDF called by
//! `HeatSolve`). The emitted dynamic symbols are exactly
//! `phasechangesolve_` + `phasechangesolve_init_` + `meltingheat_`
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/phasechangesolve_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MAXVAL`), formatted I/O (`Info`/`Warn` with `Message`,
//! `PRINT *`, `WRITE` to the isotherm report), whole-array updates, and
//! the three `CONTAINS` procedures (which mix those with deep
//! derived-type graph walks over `BoundaryInfo`, parent elements and
//! `Model` tables) stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors PhaseChangeSolve.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the Rust side owns the
//!   branch-latching SAVE scalars (`FirstTime`, `SubroutineVisited`,
//!   `Newton`); every other SAVE item (first-mesh heap: `PrevTemp`,
//!   `NodeDone`, `Nodes`, `x/y/z/Basis/NodalTemp/Conductivity/LatentHeat/
//!   TempDiff`, local matrices, `SurfaceMove`/`SurfaceMoveAve`; first-call
//!   scalars: `Trip_node`, `Axis_node`, `NoBNodes`, `NormalDir`,
//!   `TangentDirection`, `MeltPoint`, `ccum`, `prev*`, `pos0`/`prevpos0`,
//!   `Visited`, `PullVelocitySet`, `SurfaceVelocitySet`, `CoordMax/Min`,
//!   `UPull`; per-call config re-reads; and `MeltingHeat`'s SAVE `Nodes`
//!   buffer) lives as wrapper-side SAVE with identical lifetimes.
//! - `GOTO 200` (Passive Steps) is a branch to the shared
//!   pull-position tail; `FirstTime = .FALSE.` runs on all paths. The
//!   `dpos` slot is undefined on that path in the reference too
//!   (uninitialized host local); the wrapper holds a deterministic
//!   module slot instead. The proof configs never set `Passive Steps`,
//!   so that path is excluded from the differential (documented,
//!   transcribed path kept).
//! - The single-element-crossing `n == 1` edge and the empty-isotherm
//!   `Fatal` of the isotherm construction are transcribed verbatim; the
//!   proof mesh produces a multi-segment isotherm, staying on defined
//!   inputs.
//! - `MeltingHeat` crosses the C boundary as
//!   `double (void *, int *, double *)` (all dummies are interoperable:
//!   opaque `Model` handle + scalar `Node`/`t` by reference, matching the
//!   F77 calling convention; the function result travels in the floating
//!   register exactly as gfortran returns `REAL(KIND=dp)`). Its
//!   `Visited` flag is a fresh local every call in the original, so the
//!   `Normals` lookup re-runs on every call here too — no latch. The
//!   entry takes no lock (it owns no Rust state), so a `HeatSolve`
//!   re-entry into it can never deadlock against the solver entry's
//!   guard.
//! - Dead SAVE weight is dropped, never faked: `Weights`, `ForceVector`
//!   (host SAVE, shadowed per call), `FluxCorrect`, `AllocationsDone`
//!   (set, never read), `PrevUpull`, `PNodes` (host copy, never
//!   touched), and the unused host scalars/logicals. `NodalTemp` storage
//!   is allocated (for `Fatal` parity) but never filled, as in the
//!   original.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pca_begin(h_model: *const c_void, h_solver: *const c_void);
    fn pca_info_start();
    fn pca_algo(dtv: c_double, transient: c_int) -> c_int;
    fn pca_coord_dim() -> c_int;
    fn pca_fatal_2d();
    fn pca_bind_primary();
    fn pca_bind_temp();
    fn pca_cfg(relax: *mut c_double, niter: *mut c_int, ntol: *mut c_double);
    fn pca_firsttime(steady: c_int);
    fn pca_trip_temp() -> c_double;
    fn pca_passive(visit: c_int) -> c_int;
    fn pca_material();
    fn pca_trans_begin() -> c_int;
    fn pca_vel_zero();
    fn pca_vel_loads();
    fn pca_nactive() -> c_int;
    fn pca_vel_elem(t: c_int);
    fn pca_vel_asm();
    fn pca_vel_solve(norm: *mut c_double, relchange: *mut c_double);
    fn pca_vel_average();
    fn pca_pull_update();
    fn pca_speedup() -> c_double;
    fn pca_smooth() -> c_double;
    fn pca_simple_update(dtv: c_double);
    fn pca_sur_begin();
    fn pca_sur_elem(t: c_int, dtv: c_double);
    fn pca_sur_solve() -> c_double;
    fn pca_pullpos(dtv: c_double);
    fn pca_taverage(dtv: c_double);
    fn pca_newton_begin(visit: c_int, newton_in: c_int) -> c_int;
    fn pca_triple_melt(trip_temp: c_double);
    fn pca_createisotherm();
    fn pca_steady_map(newton: c_int);
    fn pca_steady_update(relax: c_double, visit: c_int, newton_in: c_int) -> c_int;
    fn pca_finish(trip_temp: c_double);
    fn pca_init_begin(h_model: *const c_void, h_solver: *const c_void);
    fn pca_init_body();
    fn pca_mh_eval(h_model: *const c_void, node: c_int, t: c_double) -> c_double;
}

/// Latching SAVE scalars (`FirstTime`, `SubroutineVisited`, `Newton`);
/// everything else SAVE lives wrapper-side with identical lifetimes.
struct PcsState {
    first_time: bool,
    subvisited: c_int,
    newton: bool,
}

impl PcsState {
    const fn new() -> Self {
        Self {
            first_time: true,
            subvisited: 0,
            newton: false,
        }
    }
}

static STATE: Mutex<PcsState> = Mutex::new(PcsState::new());

/// `phasechangesolve_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads (the body reads both).
#[export_name = "phasechangesolve_"]
pub unsafe extern "C" fn phasechange_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    let dtv: c_double = unsafe { *dt };
    let tr: c_int = unsafe { *transient };
    unsafe {
        pca_begin(model, solver);
        pca_info_start();
    }
    let steady = unsafe { pca_algo(dtv, tr) };

    state.subvisited += 1;
    let visit = state.subvisited;

    // IF(DIM /= 2) Fatal — the worker aborts inside on mismatch.
    if unsafe { pca_coord_dim() } != 2 {
        unsafe { pca_fatal_2d() };
        std::process::abort();
    }

    unsafe {
        pca_bind_primary();
        pca_bind_temp();
    }
    let (mut relax, mut niter, mut ntol) = (0.0, 0, 0.0);
    unsafe { pca_cfg(&mut relax, &mut niter, &mut ntol) };

    if state.first_time {
        unsafe { pca_firsttime(steady) };
    }

    let trip_temp = unsafe { pca_trip_temp() };

    // Passive-steps short-circuit straight to the label-200 tail.
    if unsafe { pca_passive(visit) } != 0 {
        unsafe { pca_finish(trip_temp) };
        state.first_time = false;
        return;
    }

    unsafe { pca_material() };

    if steady == 0 {
        // Transient branch: heat-flux velocity solve, then displacement.
        let useloads = unsafe { pca_trans_begin() };
        for _ in 1..=niter {
            unsafe { pca_vel_zero() };
            if useloads != 0 {
                unsafe { pca_vel_loads() };
            }
            let na = unsafe { pca_nactive() };
            for t in 1..=na {
                unsafe { pca_vel_elem(t) };
            }
            unsafe { pca_vel_asm() };
            let (mut norm, mut relchange) = (0.0, 0.0);
            unsafe { pca_vel_solve(&mut norm, &mut relchange) };
            let _ = norm;
            if relchange < ntol {
                break;
            }
        }
        unsafe {
            pca_vel_average();
            pca_pull_update();
        }
        let _speedup = unsafe { pca_speedup() };
        if unsafe { pca_smooth() } < 1.0e-20 {
            unsafe { pca_simple_update(dtv) };
        } else {
            unsafe { pca_sur_begin() };
            let na = unsafe { pca_nactive() };
            for t in 1..=na {
                unsafe { pca_sur_elem(t, dtv) };
            }
            let _norm: c_double = unsafe { pca_sur_solve() };
        }
        unsafe { pca_pullpos(dtv) };
    } else {
        // Steady branch: isotherm mapping formulation (or Newton).
        unsafe { pca_taverage(dtv) };
        let newton = unsafe { pca_newton_begin(visit, c_int::from(state.newton)) };
        unsafe { pca_triple_melt(trip_temp) };
        if newton == 0 {
            unsafe { pca_createisotherm() };
        }
        unsafe { pca_steady_map(newton) };
        state.newton = unsafe { pca_steady_update(relax, visit, newton) } != 0;
    }

    unsafe { pca_finish(trip_temp) };
    state.first_time = false;
}

/// `phasechangesolve_init_` — the solver init entry.
///
/// # Safety
///
/// Same handle contract as the main entry.
#[export_name = "phasechangesolve_init_"]
pub unsafe extern "C" fn phasechange_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        pca_init_begin(model, solver);
        pca_init_body();
    }
}

/// `meltingheat_` — the `Heat Flux` real-procedure UDF. `Model` crosses as
/// an opaque handle; `Node` and `t` cross by reference exactly as the F77
/// call passes them; the `REAL(KIND=dp)` result returns in the floating
/// register, matching gfortran's convention.
///
/// # Safety
///
/// `model` must be the Elmer address the caller passes; `node`/`t` must
/// be valid reads (they always are — the caller passes references).
#[export_name = "meltingheat_"]
pub unsafe extern "C" fn meltingheat_entry(
    model: *mut c_void,
    node: *const c_int,
    t: *const c_double,
) -> c_double {
    unsafe {
        // By-reference scalars, exactly as the F77 dummies arrive.
        pca_mh_eval(model, *node, *t)
    }
}
