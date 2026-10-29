//! Port of `fem/src/modules/MagnetoDynamics2D.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/MagnetoDynamics2D.F90` — seven
//! F77-global entries (`magnetodynamics2d_` + `_init_`,
//! `magnetodynamics2dharmonic_` + `_init_` + `_init0_`, `bsolver_` +
//! `_init_`). The emitted dynamic symbols are exactly those seven
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/magnetodynamics2d_acc.F90` (three modules `mgdyn2d_t_acc`,
//! `mgdyn2d_h_acc`, `mgdyn2d_b_acc`, one per solver entry, plus init
//! wrappers), compiled against the built `.mod` files of the target
//! configure — the port itself contains no layout assumption and no
//! pinned offset. Order-sensitive reductions (`SUM`, `MATMUL`),
//! formatted I/O (`WRITE`/`Message`) and whole-array/strided updates stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors MagnetoDynamics2D.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: each solver entry owns its
//!   module-side SAVE state under the original names (`MassValues`,
//!   `BasisFunctionsAtIp`, `BasisFunctionsInUse`, `HandlesState`,
//!   harmonic `NewtonRaphson`/`RestartDone`/..., BSolver `ForceVector`/
//!   `SaveRHS`); non-`SAVE` entry locals are re-established by the
//!   per-invocation setup wrappers, exactly matching original lifetimes.
//! - The Gauss-point kernels (`LocalMatrix*`, `BulkAssembly`, lumping,
//!   hysteresis helpers) are transcribed statement-for-statement into
//!   the accessor modules; the Rust side only mirrors the entry-level
//!   control flow (nonlinear loops, element loops, convergence and
//!   torque-tolerance branches), exactly like the original's `IF`
//!   cascade and `CYCLE`/`EXIT` actions via integer action codes.
//! - Allocation failure calls the same `Fatal(...)` inside the wrappers
//!   (the `STAT=istat` checks), never abort directly.
//! - No Rust threads: the solver calls entries serially (the original's
//!   OpenMP element loops run with one thread on the serial tree), so
//!   the `Mutex` is only a `static`-sharing guard, never contended
//!   across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    // -- transient MagnetoDynamics2D phases (mgdyn2d_t_acc) --
    fn t_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    fn t_info_start();
    fn t_setup();
    fn t_nonliniter() -> c_int;
    fn t_use_handle() -> c_int;
    fn t_default_start();
    fn t_iter_begin(iter: c_int);
    fn t_nof_active() -> c_int;
    fn t_handle_elem(t: c_int) -> c_int;
    fn t_legacy_elem(t: c_int) -> c_int;
    fn t_finish_bulk();
    fn t_nof_bnd() -> c_int;
    fn t_bc_elem(t: c_int) -> c_int;
    fn t_finish_all();
    fn t_const_mass_capture();
    fn t_solve() -> c_double;
    fn t_torque_tol(iter: c_int) -> c_int;
    fn t_converged(iter: c_int) -> c_int;
    fn t_post();
    fn t_init_run(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    // -- harmonic phases (mgdyn2d_h_acc) --
    fn h_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    fn h_setup_ns() -> c_int;
    fn h_pre_iter() -> c_int;
    fn h_iter_begin(iter: c_int);
    fn h_nof_active() -> c_int;
    fn h_iter_elem(t: c_int);
    fn h_finish_bulk();
    fn h_nof_bnd() -> c_int;
    fn h_bc_elem(t: c_int) -> c_int;
    fn h_finish_all();
    fn h_solve() -> c_double;
    fn h_converged() -> c_int;
    fn h_post();
    fn h_restart();
    fn h_pop_ns();
    fn h_dorestart() -> c_int;
    fn h_init0_run(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    fn h_init_run(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    // -- BSolver phases (mgdyn2d_b_acc) --
    fn b_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
    fn b_warn_info();
    fn b_have_matrix() -> c_int;
    fn b_count_perm() -> c_int;
    fn b_slices_check();
    fn b_prep() -> c_int;
    fn b_bulk();
    fn b_finish();
    fn b_solve_stage() -> c_double;
    fn b_norm_info(totnorm: c_double);
    fn b_init_run(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, trans: c_int);
}

static STATE: Mutex<()> = Mutex::new(());
/// Read the Fortran `(dt, TransientSimulation)` by-reference dummies the
/// loader passes: gfortran `LOGICAL` is a 4-byte 0/1, `REAL(dp)` an
/// 8-byte float. The accessor `*_begin`/`*_init_run` wrappers take both
/// by value, so dereference here (null-tolerant, as the entries never
/// observe the pointers themselves).
unsafe fn load_dt_trans(dt: *const c_double, transient: *const c_int) -> (c_double, c_int) {
    let dtv: c_double = if dt.is_null() { 0.0 } else { *dt };
    let trv: c_int = if transient.is_null() {
        0
    } else if *transient != 0 {
        1
    } else {
        0
    };
    (dtv, trv)
}

/// `magnetodynamics2d_init_` — the transient solver init entry. The whole
/// body runs through the wrapper, which owns the block verbatim
/// (including the `GetSolverParams` / list-compat traffic that cannot
/// cross the boundary).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "magnetodynamics2d_init_"]
pub unsafe extern "C" fn magnetodynamics2d_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    let (dtv, trv) = unsafe { load_dt_trans(dt, transient) };
    unsafe { t_init_run(model, solver, dtv, trv) };
}

/// `magnetodynamics2d_` — the transient solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "magnetodynamics2d_"]
pub unsafe extern "C" fn magnetodynamics2d_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe {
        let (dtv, trv) = load_dt_trans(dt, transient);
        t_begin(model, solver, dtv, trv);
        t_info_start();
        t_setup();
        let nonliniter = t_nonliniter();
        t_default_start();
        'nonlin: for iter in 1..=nonliniter {
            t_iter_begin(iter);
            let active = t_nof_active();
            // HandleAsm selects the per-element path, exactly like the
            // original's IF/ELSE around the two OMP loops.
            // NOTE: the HandleAsm flag lives wrapper-side; the two query
            // wrappers below report which path the setup chose.
            if t_use_handle() != 0 {
                for t in 1..=active {
                    t_handle_elem(t);
                }
            } else {
                for t in 1..=active {
                    t_legacy_elem(t);
                }
            }
            t_finish_bulk();
            let bnd = t_nof_bnd();
            for t in 1..=bnd {
                t_bc_elem(t);
            }
            t_finish_all();
            t_const_mass_capture();
            let _norm: c_double = t_solve();
            if t_torque_tol(iter) == 0 {
                continue 'nonlin;
            }
            let conv = t_converged(iter);
            if conv == 1 {
                break 'nonlin;
            }
            if conv == 0 {
                continue 'nonlin;
            }
        }
        t_post();
    }
}

/// `magnetodynamics2dharmonic_init0_` — harmonic pre-init entry.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "magnetodynamics2dharmonic_init0_"]
pub unsafe extern "C" fn magnetodynamics2dharmonic_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    let (dtv, trv) = unsafe { load_dt_trans(dt, transient) };
    unsafe { h_init0_run(model, solver, dtv, trv) };
}

/// `magnetodynamics2dharmonic_init_` — harmonic init entry.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "magnetodynamics2dharmonic_init_"]
pub unsafe extern "C" fn magnetodynamics2dharmonic_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    let (dtv, trv) = unsafe { load_dt_trans(dt, transient) };
    unsafe { h_init_run(model, solver, dtv, trv) };
}

/// `magnetodynamics2dharmonic_` — the harmonic solver entry.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "magnetodynamics2dharmonic_"]
pub unsafe extern "C" fn magnetodynamics2dharmonic_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe {
        let (dtv, trv) = load_dt_trans(dt, transient);
        h_begin(model, solver, dtv, trv);
        if h_setup_ns() == 1 {
            return;
        }
        let nonliniter = h_pre_iter();
        for iter in 1..=nonliniter {
            h_iter_begin(iter);
            let active = h_nof_active();
            for t in 1..=active {
                h_iter_elem(t);
            }
            h_finish_bulk();
            let bnd = h_nof_bnd();
            for t in 1..=bnd {
                h_bc_elem(t);
            }
            h_finish_all();
            let _norm: c_double = h_solve();
            if h_converged() == 1 {
                break;
            }
        }
        h_post();
        if h_dorestart() != 0 {
            h_restart();
        }
        h_pop_ns();
    }
}

/// `bsolver_init_` — the BSolver init entry (keyword/default traffic
/// only; whole body in the wrapper).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "bsolver_init_"]
pub unsafe extern "C" fn bsolver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    let (dtv, trv) = unsafe { load_dt_trans(dt, transient) };
    unsafe { b_init_run(model, solver, dtv, trv) };
}

/// `bsolver_` — the B-field postprocessing entry.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "bsolver_"]
pub unsafe extern "C" fn bsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe {
        let (dtv, trv) = load_dt_trans(dt, transient);
        b_begin(model, solver, dtv, trv);
        b_warn_info();
        // IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN
        if b_have_matrix() == 0 {
            return;
        }
        // IF (COUNT(...) <= 0) RETURN
        if b_count_perm() <= 0 {
            return;
        }
        b_slices_check();
        let _totdofs = b_prep();
        b_bulk();
        b_finish();
        let totnorm: c_double = b_solve_stage();
        b_norm_info(totnorm);
    }
}
