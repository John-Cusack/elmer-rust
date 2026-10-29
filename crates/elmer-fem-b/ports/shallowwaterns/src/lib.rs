//! Port of `fem/src/modules/ShallowWaterNS.F90` (Elmer @ 9f6af2f85, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ShallowWaterNS.F90` — one F77-global
//! entry `shallowwaternssolver_` with three `CONTAINS` internals. The
//! emitted dynamic symbol is exactly `shallowwaternssolver_` (proven by
//! `nm -D` on the reference `.so`); the internals emit only locals.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/shallowwaterns_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`), whole-array fills/updates, and the two
//! `WRITE(Message, *)` + `Info` pairs stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors ShallowWaterNS.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the entry's `NewtonLinearization`
//!   flag carries the original's initialization-implies-`SAVE` lifetime
//!   — it latches process-wide once the
//!   `RelativeChange < NewtonTol .OR. iter >= NewtonIter` condition
//!   fires (set in a later nonlinear iteration or a later timestep,
//!   never reset), mirrored here by a process-wide latch that starts
//!   `.FALSE.`; every other local is a per-call automatic. The
//!   internal procedures' `Nodes`/`IP` locals live as wrapper-side
//!   `SAVE` state with identical lifetimes (consumed before the next
//!   fill, as in the original). Per-element work arrays are
//!   wrapper-side `SAVE` allocatables that only grow past the
//!   first-mesh sizing (allocation failure calls the same `Fatal`).
//! - `Model` and `dt` are unused by the body, exactly as in the
//!   original (received and ignored); `Transient` selects the
//!   `Default1stOrderTime` path in both local-matrix tails.
//! - The nonlinear loop (`NonlinIter`, defaults and early `EXIT`),
//!   the `MassFluxDiscretization` early return in
//!   `BoundaryLocalMatrix`, and the `IF (dim == 1)` / `MAX(j, nd-2)`
//!   bubble fixup are kept with the original's evaluation order.
//! - `BulkAssembly` / `BoundaryAssembly` run serially in index order,
//!   which is what the original's `!$omp` loops lower to in the
//!   serial (`-O0`, no OpenMP) proof build.
//! - No Rust threads: the solver calls entries serially, so the latch
//!   is a `static`-sharing guard, never contended across threads
//!   (`unsafe` only at the documented handle boundary).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::atomic::{AtomicBool, Ordering};

extern "C" {
    fn swa_begin(h_solver: *const c_void);
    fn swa_newton_iter() -> c_int;
    fn swa_newton_tol() -> c_double;
    fn swa_nonlin_iter() -> c_int;
    fn swa_nonlin_tol() -> c_double;
    fn swa_massflux() -> c_int;
    fn swa_nof_active() -> c_int;
    fn swa_nof_boundary() -> c_int;
    fn swa_default_initialize();
    fn swa_default_finish_bulk();
    fn swa_default_finish_bnd();
    fn swa_default_finish();
    fn swa_default_dirichlet();
    fn swa_default_solve() -> c_double;
    fn swa_nonlin_change() -> c_double;
    fn swa_info(norm: c_double, rel: c_double);
    fn swa_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn swa_local_ensure(n: c_int, nd: c_int);
    fn swa_local_fetch(
        n: c_int,
        nd: c_int,
        g: *mut c_double,
        icz: *mut c_double,
        wcoeff: *mut c_double,
        fcor: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn swa_local_assemble(
        n: c_int,
        nd: c_int,
        newton: c_int,
        massflux: c_int,
        g: c_double,
        icz: c_double,
        wcoeff: c_double,
        fcor: c_double,
    );
    fn swa_local_update(nd: c_int, transient: c_int);
    fn swa_begin_boundary(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn swa_bc() -> c_int;
    fn swa_boundary_assemble(n: c_int, nd: c_int, transient: c_int);
}

/// Process-wide `NewtonLinearization` latch. The original declares it
/// with initialization (`NewtonLinearization=.FALSE.`), which implies
/// the `SAVE` attribute: once the
/// `RelativeChange < NewtonTol .OR. iter >= NewtonIter` condition fires
/// — in a later nonlinear iteration or a later timestep — it stays
/// `.TRUE.` for the rest of the process. Elmer invokes solver entries
/// serially, so a relaxed atomic is a faithful sharing guard here.
static NEWTON_LATCH: AtomicBool = AtomicBool::new(false);
/// `shallowwaternssolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read if non-null (null reads as false).
#[export_name = "shallowwaternssolver_"]
pub unsafe extern "C" fn shallowwaterns_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let is_transient = !transient.is_null() && unsafe { transient.read() } != 0;
    let tflag: c_int = if is_transient { 1 } else { 0 };

    unsafe { swa_begin(solver) };

    // Params => GetSolverParams(Solver) lives in the wrapper; the
    // Found-defaults below mirror the original statements exactly.
    let newton_iter = unsafe { swa_newton_iter() };
    let newton_tol = unsafe { swa_newton_tol() };
    let nonlin_iter = unsafe { swa_nonlin_iter() };
    let nonlin_tol = unsafe { swa_nonlin_tol() };
    let massflux = unsafe { swa_massflux() };

    // NewtonLinearization: initialization-implies-SAVE latch — starts
    // from the process-wide latched value, never reset within the run.
    let mut newton_lin = NEWTON_LATCH.load(Ordering::Relaxed);

    // DO iter=1,NonlinIter
    for iter in 1..=nonlin_iter {
        unsafe { swa_default_initialize() };
        bulk_assembly(tflag, newton_lin, massflux);
        unsafe { swa_default_finish_bulk() };

        boundary_assembly(tflag, massflux);
        unsafe { swa_default_finish_bnd() };

        unsafe {
            swa_default_finish();
            swa_default_dirichlet();
        }

        // Norm = DefaultSolve(); RelativeChange = Solver % Variable % NonlinChange
        let norm: c_double = unsafe { swa_default_solve() };
        let rel: c_double = unsafe { swa_nonlin_change() };
        unsafe { swa_info(norm, rel) };

        if rel < newton_tol || iter >= newton_iter {
            newton_lin = true;
            NEWTON_LATCH.store(true, Ordering::Relaxed);
        }
        if rel < nonlin_tol {
            break;
        }
    }
}

/// The `BulkAssembly` internal procedure: serial active-element loop
/// (what the original's `!$omp` loop lowers to in the serial build).
fn bulk_assembly(tflag: c_int, newton_lin: bool, massflux: c_int) {
    let newton: c_int = if newton_lin { 1 } else { 0 };
    let active = unsafe { swa_nof_active() };
    for t in 1..=active {
        let (mut n, mut nd) = (0, 0);
        unsafe { swa_begin_active(t, &mut n, &mut nd) };
        unsafe { swa_local_ensure(n, nd) };
        let (mut g, mut icz, mut wcoeff, mut fcor) = (0.0, 0.0, 0.0, 0.0);
        unsafe { swa_local_fetch(n, nd, &mut g, &mut icz, &mut wcoeff, &mut fcor) };
        unsafe { swa_local_assemble(n, nd, newton, massflux, g, icz, wcoeff, fcor) };
        unsafe { swa_local_update(nd, tflag) };
    }
}

/// The `BoundaryAssembly` internal procedure: serial boundary-element
/// loop with the `ActiveBoundaryElement` skip, the
/// `MassFluxDiscretization` early return, and the `GetBC` guard, in
/// the original's order.
fn boundary_assembly(tflag: c_int, massflux: c_int) {
    let nbnd = unsafe { swa_nof_boundary() };
    for t in 1..=nbnd {
        let (mut n, mut nd, mut active_el) = (0, 0, 0);
        unsafe { swa_begin_boundary(t, &mut n, &mut nd, &mut active_el) };
        // IF(.NOT.ActiveBoundaryElement(Element)) CYCLE
        if active_el == 0 {
            continue;
        }
        // IF(.NOT.MassFluxDiscretization) RETURN (per element)
        if massflux == 0 {
            continue;
        }
        // BC => GetBC(Element); IF(.NOT.ASSOCIATED(BC)) RETURN
        if unsafe { swa_bc() } == 0 {
            continue;
        }
        unsafe { swa_local_ensure(n, nd) };
        unsafe { swa_boundary_assemble(n, nd, tflag) };
    }
}
