//! Port of `fem/src/modules/BeamSolver3D.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/BeamSolver3D.F90` — two F77-global
//! entries (`timoshenkosolver_`, `timoshenkosolver_init0_`) with no
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those
//! two (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/beamsolver3d_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the Timoshenko beam assembly in
//! `SolidMechanicsUtils`, `MATMUL` reductions, the second-order time
//! discretization) and the optional-arg `BeamStiffnessMatrix` calls stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors BeamSolver3D.F90 statement for statement):
//! - The nonlinear iteration loop maps to the same `DefaultInitialize` /
//!   per-element `BeamStiffnessMatrix` / finish / `DefaultSolve` /
//!   `DefaultConverged` sequence with identical exit semantics (`DO
//!   iter=1,maxiter` with `EXIT` on convergence; an empty range runs
//!   zero times, as in the original).
//! - `maxiter` comes from the same `ListGetInteger(..., minv=1)` call
//!   with the same `IF (.NOT. Found) maxiter = 1` fallback.
//! - `MassAssembly = TransientSimulation .OR. HarmonicAssembly` uses the
//!   entry's logical dummy, exactly as in the original.
//! - `Norm = DefaultSolve()` is kept and discarded, exactly as in the
//!   original. `Model` and `dt` are unused by the entry body
//!   (received and ignored), exactly as in the original.
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn b3d_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn b3d_default_start();
    fn b3d_check_axes();
    fn b3d_maxiter() -> c_int;
    fn b3d_modes(transient: c_int, harmonic: *mut c_int, mass: *mut c_int);
    fn b3d_displace_active() -> c_int;
    fn b3d_mesh_return();
    fn b3d_displace_forward();
    fn b3d_default_initialize();
    fn b3d_nof_active() -> c_int;
    fn b3d_elem(t: c_int, transient: c_int, mass: c_int, harmonic: c_int);
    fn b3d_finish_bulk();
    fn b3d_finish_boundary();
    fn b3d_finish_assembly();
    fn b3d_dirichlet();
    fn b3d_solve() -> c_double;
    fn b3d_converged() -> c_int;
    fn b3d_finish();
    fn b3d_init0(h_solver: *mut c_void, transient: c_int);
}

/// `timoshenkosolver_init0_` — the `_Init0` entry Elmer calls once at
/// solver setup. Mirrors `TimoshenkoSolver_Init0` (default `Variable`,
/// element definition, mass-matrix and timestepping defaults).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (only `transient`
/// is used, exactly as in the original).
#[export_name = "timoshenkosolver_init0_"]
pub unsafe extern "C" fn timoshenko_solver_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let t: c_int = if transient.is_null() {
        0
    } else {
        unsafe { *transient }
    };
    unsafe {
        b3d_init0(solver, t);
    }
}

/// `timoshenkosolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`timoshenko_solver_init0_entry`]; `model` is the
/// Elmer `Model_t` address the loader passes.
#[export_name = "timoshenkosolver_"]
pub unsafe extern "C" fn timoshenko_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let t: c_int = if transient.is_null() {
        0
    } else {
        unsafe { *transient }
    };
    unsafe {
        b3d_begin(model, solver);
        b3d_default_start();
        b3d_check_axes();
        let maxiter = b3d_maxiter();
        let (mut harmonic, mut mass) = (0, 0);
        b3d_modes(t, &mut harmonic, &mut mass);
        let displace = b3d_displace_active();
        if displace != 0 {
            b3d_mesh_return();
        }
        for _iter in 1..=maxiter {
            b3d_default_initialize();
            let active = b3d_nof_active();
            for k in 1..=active {
                b3d_elem(k, t, mass, harmonic);
            }
            b3d_finish_bulk();
            b3d_finish_boundary();
            b3d_finish_assembly();
            b3d_dirichlet();
            let _norm: c_double = b3d_solve();
            if b3d_converged() != 0 {
                break;
            }
        }
        b3d_finish();
        if displace != 0 {
            b3d_displace_forward();
        }
    }
}
