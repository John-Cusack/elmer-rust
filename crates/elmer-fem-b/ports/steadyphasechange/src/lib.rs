//! Port of `fem/src/modules/SteadyPhaseChange.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/SteadyPhaseChange.F90` — two F77-global
//! entries, `steadyphasechange_` (the Lagrangian steady-state phase-change
//! solver) and `meltingheat_` (the `Heat Flux` real-procedure UDF called by
//! `HeatSolve`'s `LocalMatrixBC`). The emitted dynamic symbols are exactly
//! `steadyphasechange_` + `meltingheat_` (proven by `nm -D` on the
//! reference `.so`); there are no other globals.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/steadyphasechange_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MAXVAL`, `MATMUL`-free but order-sensitive isotherm mapping),
//! formatted I/O (`Info`/`Warn` with `Message`, `PRINT *`, `WRITE` to
//! `isotherm.dat`), and whole-array updates stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors SteadyPhaseChange.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the Rust side owns the
//!   branch-latching SAVE scalars (`FirstTime`, `SubroutineVisited`,
//!   `Newton`); every other SAVE item (first-mesh heap: `PrevTemp`,
//!   `NodeDone`, `Nodes`, `x/y/z/Basis/NodalTemp/TempDiff`,
//!   `BoundaryMarker`, `PhaseElements`; first-call scalars:
//!   `Trip_node`, `MeltPoint`, `ccum`, `prevtave`, ...; per-call
//!   accumulators; `IsoSurf` + flag; and the whole `MeltingHeat` SAVE
//!   state `Initialized/Nodes/NormalExist/NormalSol`) lives as
//!   wrapper-side SAVE with identical lifetimes (allocated once on the
//!   first call, `IsoSurf` allocated/freed around each `CreateIsotherm`
//!   pass pair, `MeltingHeat` initialized once — exactly as in the
//!   original).
//! - Dead SAVE weight is dropped, never faked: `ForceVector`,
//!   `FluxCorrect`, `Weights`, `SurfaceVelocitySet`, `Visited`, `pos0`,
//!   `prevpos0`, `MaxSurface`-family locals, `Eps`-host (recomputed per
//!   use), `NodalTemp` values (allocated but never filled — the
//!   allocation itself is kept for `Fatal` parity), and the unused `str`
//!   are never read on any output-affecting path.
//! - `MeltingHeat` crosses the C boundary as
//!   `double (void *, int *, double *)` (all dummies are interoperable:
//!   opaque `Model` handle + scalar `Node`/`t` by reference, matching the
//!   F77 calling convention; the function result travels in the floating
//!   register exactly as gfortran returns `REAL(KIND=dp)`). The serial
//!   proof tree has no OpenMP (`WITH_OpenMP=FALSE`), so the `!$` lines
//!   are comments there: `tid` is always 1 and a single `Nodes` slot is
//!   allocated — the wrapper mirrors that behavior.
//! - `Passive Steps` short-circuits to the label-200 tail exactly as in
//!   the original (via `sp_finish` + clearing `FirstTime`). Note the
//!   reference leaves `dpos` undefined on that path (uninitialized
//!   host local); the wrapper holds a deterministic module slot instead.
//!   The proof configs never set `Passive Steps`, so that path is
//!   excluded from the differential (documented, transcribed path kept).
//! - The single-element-crossing `n == 1` / empty-isotherm edges of
//!   `CreateIsotherm` (including the `imin`-undefined corner when no
//!   segment is found) are transcribed verbatim; the proof mesh produces
//!   a multi-segment isotherm, staying on defined inputs.
//! - `Model`, `dt`, `TransientSimulation` are unused by the solver body,
//!   exactly as in the original (received and ignored); `t` is likewise
//!   unused by `MeltingHeat`.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads. The
//!   `MeltingHeat` entry takes no lock (it owns no Rust state), so a
//!   `HeatSolve` re-entry into it can never deadlock against the solver
//!   entry's guard.
//! - The Rust side allocates no heap: every `STAT=istat` allocation check
//!   lives in the wrappers and calls the same `Fatal(...)` spellings as
//!   the original (`PhaseChangeSolver` for errors 1–3, `SteadyPhaseChange`
//!   for 4–5), instead of aborting.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sp_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sp_info_start();
    fn sp_dim() -> c_int;
    fn sp_fatal_2d();
    fn sp_surface_begin();
    fn sp_surfacemove_begin();
    fn sp_temperature_begin();
    fn sp_get_relax() -> c_double;
    fn sp_firsttime();
    fn sp_trip_temp() -> c_double;
    fn sp_passive_steps() -> c_int;
    fn sp_newton_update(visit: c_int, newton_in: c_int) -> c_int;
    fn sp_meltpoint(trip_temp: c_double);
    fn sp_createisotherm();
    fn sp_loopinit();
    fn sp_nbound() -> c_int;
    fn sp_elem_begin(t: c_int, n: *mut c_int, family: *mut c_int, doelem: *mut c_int);
    fn sp_node_update(nn: c_int, te: c_int, newton: c_int);
    fn sp_elem_integrate();
    fn sp_postloop(relax: c_double, visit: c_int, newton_in: c_int) -> c_int;
    fn sp_internal_mesh_movement() -> c_int;
    fn sp_boxmovemesh();
    fn sp_finish(trip_temp: c_double);
    fn sp_mh_init(h_model: *const c_void);
    fn sp_mh_eval(h_model: *const c_void, node: c_int, t: c_double) -> c_double;
}

/// First-call/branch-latching SAVE scalars (`FirstTime`,
/// `SubroutineVisited`, `Newton`); everything else SAVE lives
/// wrapper-side (see the module docs).
struct SteadyState {
    first_time: bool,
    subroutine_visited: c_int,
    newton: bool,
}

impl SteadyState {
    const fn new() -> Self {
        Self {
            first_time: true,
            subroutine_visited: 0,
            newton: false,
        }
    }
}

static STATE: Mutex<SteadyState> = Mutex::new(SteadyState::new());

/// `steadyphasechange_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "steadyphasechange_"]
pub unsafe extern "C" fn steadyphasechange_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sp_begin(model, solver) };
    unsafe { sp_info_start() };

    state.subroutine_visited += 1;
    let visit = state.subroutine_visited;

    // IF(DIM /= 2) Fatal — the worker query aborts inside on mismatch.
    if unsafe { sp_dim() } != 2 {
        unsafe { sp_fatal_2d() };
        return;
    }

    unsafe {
        sp_surface_begin();
        sp_surfacemove_begin();
        sp_temperature_begin();
    }
    let relax = unsafe { sp_get_relax() };

    if state.first_time {
        unsafe { sp_firsttime() };
    }

    let trip_temp = unsafe { sp_trip_temp() };

    // Passive-steps short-circuit straight to the label-200 tail.
    if unsafe { sp_passive_steps() } >= visit {
        unsafe { sp_finish(trip_temp) };
        state.first_time = false;
        return;
    }

    let newton = unsafe { sp_newton_update(visit, c_int::from(state.newton)) };
    unsafe { sp_meltpoint(trip_temp) };

    if newton == 0 {
        unsafe { sp_createisotherm() };
    }

    unsafe { sp_loopinit() };
    let nbnd = unsafe { sp_nbound() };
    for t in 1..=nbnd {
        let (mut n, mut family, mut doelem) = (0, 0, 0);
        unsafe { sp_elem_begin(t, &mut n, &mut family, &mut doelem) };
        if doelem == 0 {
            continue;
        }
        let _ = family;
        for nn in 1..=n {
            // Note the original's `TTemp = Temperature(TempPerm(t))`:
            // the *element* counter travels as the node index — kept.
            unsafe { sp_node_update(nn, t, newton) };
        }
        unsafe { sp_elem_integrate() };
    }

    state.newton = unsafe { sp_postloop(relax, visit, newton) } != 0;

    if unsafe { sp_internal_mesh_movement() } != 0 {
        unsafe { sp_boxmovemesh() };
    }

    unsafe { sp_finish(trip_temp) };
    state.first_time = false;
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
        sp_mh_init(model);
        // By-reference scalars, exactly as the F77 dummies arrive.
        sp_mh_eval(model, *node, *t)
    }
}
