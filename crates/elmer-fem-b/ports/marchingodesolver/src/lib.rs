//! Port of `fem/src/modules/MarchingODESolver.F90` (the campaign fork) for
//! lane T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib
//! with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/MarchingODESolver.F90` — two F77-global
//! entries `marchingodesolver_` and `marchingodesolver_init_` (proven by
//! `nm -D` on the reference `.so`); the four `CONTAINS` internals
//! (`FirstExtrudedMaterial`, `SetLayerValues`, `GetLayerValues`,
//! `GetCoefficients`) emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/marchingodesolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive whole-array
//! updates (`xvec = x0vec + ...`, `SUM` reductions, `MAX`/`SQRT` layer
//! scaling), the list-directed `PRINT`s, `WRITE(Message, ...)` + `Info`
//! lines, and the bare `STOP` stay in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors MarchingODESolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the accessor module holds
//!   every name from the original's `SAVE` list (plus the `=init`
//!   implied-`SAVE` locals `MaskExist`, `ParabolicModel`, `RequireBC`,
//!   `DoTransient`, `AnyDG`, `VectorSource`) with identical lifetimes —
//!   set once by `moa_initialize`, never reset per call. The
//!   `GetCoefficients` SAVE state (`Visited`, `Source_h`, `xloc`,
//!   `cdofs`) likewise persists. Per-call locals are re-set on every
//!   entry call, exactly as the original re-executes them.
//! - The original's `1 CONTINUE` / `GOTO 1` marching-step loop becomes an
//!   explicit `1..=npass` passage loop with `npass = max(dtn, 1)`; the
//!   `Taking marching step` info fires between passages exactly when the
//!   original's `IF (dti < dtn)` guard holds.
//! - `Model`, `dt`, and `TransientSimulation` travel into the wrappers
//!   (the body reads the draw velocity, transient flag, and timestep
//!   from them), exactly as in the original.
//! - Degenerate `Nonlinear System Max Iterations < 1` (zero-trip
//!   iteration loop, uninitialized `Norm`/`Change` in the original) is
//!   out of scope: no shipped `.sif` sets it, and the original prints
//!   uninitialized stack there.
//! - No Rust threads: the solver calls entries serially, and all shared
//!   state lives in the Fortran wrappers, so there is nothing to guard.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn moa_init_params();
    fn moa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn moa_initialized() -> c_int;
    fn moa_initialize();
    fn moa_set_dims();
    fn moa_beta() -> c_double;
    fn moa_parabolic();
    fn moa_velo();
    fn moa_maxiter() -> c_int;
    fn moa_nonlintol() -> c_double;
    fn moa_material_setup();
    fn moa_transient(dt: c_double, transient: c_int) -> c_int;
    fn moa_passage_begin() -> c_int;
    fn moa_nlayers() -> c_int;
    fn moa_advance_layer(layer: c_int, dtime: *mut c_double);
    fn moa_sub_head(subt: c_int);
    fn moa_iter_head(iter: c_int, maxiter: c_int);
    fn moa_solve_iter(beta: c_double, change: *mut c_double);
    fn moa_sub_tail();
    fn moa_info_step(dti: c_int);
    fn moa_info_done();
}

/// `marchingodesolver_init_` — the initializer entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; only the solver context
/// matters (via `GetSolverParams`), the rest are unused, exactly as in
/// the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "marchingodesolver_init_"]
pub unsafe extern "C" fn marchingodesolver_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { moa_init_params() };
}

/// `marchingodesolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads (they carry the timestep size
/// and the transient flag the marching logic needs).
#[export_name = "marchingodesolver_"]
pub unsafe extern "C" fn marchingodesolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { moa_begin(model, solver) };

    // IF (.NOT. Initialized) THEN <structured-mesh + ODE setup> END IF
    if unsafe { moa_initialized() } == 0 {
        unsafe { moa_initialize() };
    }

    // The variable on the layer: n = BotNodes; m = dofs * n
    unsafe { moa_set_dims() };

    // Beta / Timestepping Method branch.
    let beta = unsafe { moa_beta() };
    unsafe { moa_parabolic() };
    unsafe { moa_velo() };
    let maxiter = unsafe { moa_maxiter() };
    let nonlintol = unsafe { moa_nonlintol() };
    unsafe { moa_material_setup() };

    // dtn = 0; IF (Transient) THEN <match mesh step> END IF
    // DoTransient = (dtn > 0)
    let dtn = unsafe { moa_transient(*dt, *transient) };

    // dti = 1; 1 CONTINUE ... IF (dti < dtn) dti += 1, GOTO 1.
    // dtn <= 1 runs a single passage; the step info fires between
    // passages exactly when the original guard holds.
    let npass = dtn.max(1);
    let nlayers = unsafe { moa_nlayers() };
    for p in 1..=npass {
        let maxsubt = unsafe { moa_passage_begin() };
        for layer in 1..=nlayers {
            let mut dtime: c_double = 0.0;
            unsafe { moa_advance_layer(layer, &mut dtime) };
            for subt in 1..=maxsubt {
                unsafe { moa_sub_head(subt) };
                // DO iter = 1, MaxIter ... IF (Change < NonLinTol) EXIT
                for iter in 1..=maxiter {
                    unsafe { moa_iter_head(iter, maxiter) };
                    let mut change: c_double = 0.0;
                    unsafe { moa_solve_iter(beta, &mut change) };
                    if change < nonlintol {
                        break;
                    }
                }
                unsafe { moa_sub_tail() };
            }
        }
        if p < npass {
            unsafe { moa_info_step(p + 1) };
        }
    }

    unsafe { moa_info_done() };
}
