//! Port of `fem/src/modules/FacetShellSolve.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FacetShellSolve.F90` — two F77-global
//! entries (`shellsolver_`, `shellsolver_init0_`) with fourteen
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those
//! two (proven by `nm -D` on the reference `.so`); the internals emit
//! only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/facetshellsolve_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. All floating-point work (the
//! MITC shell assembly, `MATMUL`/`SUM` reductions, `SQRT`/`ACOS`
//! intrinsics, the `CPUTime`/`RealTime` reports) and all formatted/file
//! I/O stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FacetShellSolve.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the entry SAVE list
//!   (first-mesh `STIFF`/`MASS`/`DAMP`/`FORCE`/material arrays sized by
//!   `MIN(MaxElementNodes, 4)`, the `AllocationsDone` flag, the nodal
//!   stress pointers, the `LoadVector` staging) and the shared host
//!   variables of the internal procedures live as wrapper-side SAVE
//!   state with identical lifetimes. Non-SAVE entry locals retain their
//!   values across invocations here (module association) where the
//!   original left them undefined; every proof path assigns them before
//!   reading (`RelChange`/`RelUChange` are always set for loaded
//!   configs, `LoadScale` and `NumberOfElementNodes` are set whenever
//!   bulk elements exist).
//! - The `_Init0` entry transcribes the default `Deflection` variable
//!   (6 DOFs), the time derivative order, and the Shell/Drilling flags.
//! - The nonlinear loop bound and the `RelChange` exit test are driven
//!   from Rust through the `fss_maxit`/`fss_postsolve` queries; the loop
//!   body (bulk/boundary assembly, concentrated loads, `LoadScale`
//!   scaling, Dirichlet pass, `DefaultSolve`, timing reports) runs in
//!   the wrappers in original order.
//! - The `StabilityAnalysis` and transient (`Default2ndOrderTime`)
//!   branches are transcribed verbatim but are not exercised by the
//!   proof configs; the `HARMONIC`/`EigenOrHarmonicAnalysis` mass
//!   branches likewise.
//! - Allocation failure calls the same `Fatal('ShellSolver', ...)` via
//!   the original `STAT=istat` check. `Model`, `dt`, and the transient
//!   flag travel exactly as in the original: the model address feeds
//!   the wrappers that read `Model % Bodies`/`BodyForces`; `dt` is
//!   unused by the body, exactly as in the original, and is not passed;
//!   the transient flag feeds the density guard and the second-order
//!   time branch.
//! - The Rust side holds no state at all, so no `Mutex` is needed: the
//!   solver entries run serially, exactly as in the original.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fss_init_begin(h_solver: *mut c_void);
    fn fss_init_body();
    fn fss_begin(h_model: *mut c_void, h_solver: *mut c_void, transient: c_int);
    fn fss_head();
    fn fss_maxit() -> c_int;
    fn fss_iter_head(iter: c_int);
    fn fss_assemble();
    fn fss_solver_disabled() -> c_int;
    fn fss_presolve();
    fn fss_solve();
    fn fss_postsolve() -> c_int;
    fn fss_stress();
    fn fss_tail();
}

/// `shellsolver_init0_` — the `_Init0` entry Elmer calls once at solver
/// setup. Mirrors `ShellSolver_Init0` (default 6-DOF `Deflection`
/// variable, time derivative order, Shell/Drilling flags).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "shellsolver_init0_"]
pub unsafe extern "C" fn shellsolver_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fss_init_begin(solver);
        fss_init_body();
    }
}

/// `shellsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `transient` must be a valid read (it feeds the density guard and the
/// second-order time branch). `dt` may be null (it is ignored, exactly as
/// in the original).
#[export_name = "shellsolver_"]
pub unsafe extern "C" fn shellsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let tr: c_int = if transient.is_null() { 0 } else { *transient };
        fss_begin(model, solver, tr);
        fss_head();
        let maxit = fss_maxit();
        let mut iter: c_int = 0;
        while iter < maxit {
            iter += 1;
            fss_iter_head(iter);
            fss_assemble();
            if fss_solver_disabled() == 0 {
                fss_presolve();
                fss_solve();
                if fss_postsolve() != 0 {
                    break;
                }
            }
        }
        fss_stress();
        fss_tail();
    }
}
