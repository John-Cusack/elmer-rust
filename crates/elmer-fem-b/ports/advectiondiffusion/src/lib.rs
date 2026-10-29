//! Port of `fem/src/modules/AdvectionDiffusion.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/AdvectionDiffusion.F90` — two F77-global
//! entries (`advectiondiffusionsolver_`, `advectiondiffusionsolver_init_`)
//! with six `CONTAINS` internals. The emitted dynamic symbols are exactly
//! those two (proven by `nm -D` on the reference `.so`); the internals
//! emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/advectiondiffusion_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. All floating-point work (the
//! Galerkin assembly, `SUM` reductions, the `CPUTime`-based reports),
//! whole-array updates, and every formatted `WRITE`/`PRINT` stay in
//! Fortran wrappers (same compiler), so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors AdvectionDiffusion.F90 statement for
//! statement; the Rust side only owns the two entry loops and the
//! steady-iteration loop nest, every body statement lives in Fortran):
//! - `SAVE` replicated exactly: the whole entry local state (assembly
//!   arrays `LocalMassMatrix`/`LocalStiffMatrix`/`Load`/`C0`/`C1`/`CT`/
//!   `C2`/`Diffusivity`/`LocalForce`/`TimeForce`, nodal work arrays,
//!   `ElementNodes`/`ParentNodes`, `Mass`, `Hwrk`, `AllocationsDone`, and
//!   every per-invocation local that must survive across the per-element
//!   wrapper calls) lives as wrapper-module SAVE with identical
//!   lifetimes. The six `CONTAINS` internals are verbatim module
//!   procedures (their host-associated `Diffusivity` resolves to the same
//!   module SAVE array).
//! - `Relax`, `EquationName`, and the first `Norm` read are assigned but
//!   never read in the original; the reads are replicated and discarded.
//! - `ElementList` (automatic array sized by bulk elements) and
//!   `BackupForceVector` are declared but never referenced in the
//!   original; omitted (no executable effect).
//! - `Model`, `dt` are received and re-stored (`dt` feeds the
//!   `Add1stOrderTime` calls, `Model` serves every `Model %` access),
//!   exactly as in the original; `TransientSimulation` gates the same
//!   branches through the stored flag.
//! - `CPUTIme()` is the same external as `CPUTime()` (Fortran is
//!   case-insensitive); the wrapper keeps the original spelling.
//! - Allocation failure calls the same `Fatal('AdvectionDiffusion', ...)`
//!   via the `STAT=istat` check.
//! - No Rust threads: the solver entries run serially; the Rust side holds
//!   no shared state at all (all SAVE state is wrapper-side).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ada_init_begin(h_solver: *const c_void);
    fn ada_init_body();
    fn ada_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, transient: c_int);
    fn ada_preamble() -> c_int;
    fn ada_nonlinear_max() -> c_int;
    fn ada_iter_start(iter: c_int);
    fn ada_tmax() -> c_int;
    fn ada_bulk(t: c_int);
    fn ada_bulk_finish();
    fn ada_mixed();
    fn ada_boundary_start();
    fn ada_nbulk() -> c_int;
    fn ada_nboundary() -> c_int;
    fn ada_boundary(t: c_int);
    fn ada_assembly_finish();
    fn ada_solve_tail(iter: c_int, maxiter: c_int) -> c_int;
    fn ada_mass_check();
}

/// `advectiondiffusionsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model` and `solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null.
#[export_name = "advectiondiffusionsolver_"]
pub unsafe extern "C" fn advectiondiffusion_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let dtv = if dt.is_null() { 0.0 } else { *dt };
        let tr = if transient.is_null() { 0 } else { *transient };
        ada_begin(model, solver, dtv, tr);
        if ada_preamble() == 1 {
            return;
        }
        let niter = ada_nonlinear_max();
        for iter in 1..=niter {
            ada_iter_start(iter);
            let tmax = ada_tmax();
            for t in 1..=tmax {
                ada_bulk(t);
            }
            ada_bulk_finish();
            ada_mixed();
            ada_boundary_start();
            let nbulk = ada_nbulk();
            let nbound = ada_nboundary();
            for b in 1..=nbound {
                ada_boundary(nbulk + b);
            }
            ada_assembly_finish();
            if ada_solve_tail(iter, niter) == 1 {
                break;
            }
        }
        ada_mass_check();
    }
}

/// `advectiondiffusionsolver_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `AdvectionDiffusionSolver_init` (adds `Time
/// Derivative Order = 1`).
///
/// # Safety
///
/// Same contract as [`advectiondiffusion_solver_entry`].
#[export_name = "advectiondiffusionsolver_init_"]
pub unsafe extern "C" fn advectiondiffusion_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ada_init_begin(solver);
        ada_init_body();
    }
}
