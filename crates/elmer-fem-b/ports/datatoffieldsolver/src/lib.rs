//! Port of `fem/src/modules/DataToFieldSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/DataToFieldSolver.F90` — two F77-global
//! entries (`datatofieldsolver_`, `datatofieldsolver_init_`) with four
//! `CONTAINS` internals. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/datatoffieldsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. All floating-point work (basis
//! evaluation, the Galerkin assembly, `MATMUL`/`SUM` reductions,
//! `LOG`/`EXP`, the `ParallelReduction` weight correction, the
//! diagonal/rhs accumulation) and all file/formatted I/O (the
//! list-directed point-data read, `WRITE(Message, ...)`) stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors DataToFieldSolver.F90 statement for
//! statement):
//! - The Rust entry drives the structural branches in the original order:
//!   point-data file present → `AsciiPointsToMesh`; `Continue BC With
//!   Boundary Assembly` → `BoundaryAssembly`; `Revert Logarithmic
//!   Fitting` → revert-solve path (`EXP` + `ComputeChange`) vs plain
//!   `DefaultSolve`. String-keyword fallback chains (`Data Field Name` /
//!   `Target Variable`, `Weight Field Name` / `Weight Variable`, `Mask
//!   Field Name` / `Mask Variable`, `Normalize Data by Weight` alias)
//!   and every `Fatal` check stay inside the resolve wrapper in the
//!   original order.
//! - `SAVE` semantics replicated: `AsciiPointsToMesh` state
//!   (`ElementIndex`, `Basis`, `dBasisdx`, `ElementNodes`), `BulkAssembly`
//!   state (`Visited`, `Nodes`, `STIFF`, `FORCE`, `Basis`, `dBasisdx`,
//!   `NodalWeight`, `BoundaryPerm`, `BoundaryNodes`, `DataDiffusivity`),
//!   and `BoundaryLocalMatrix` state (`Nodes`, `ParentNodes`) live as
//!   wrapper-side SAVE state with identical lifetimes. One documented
//!   repair: the original declares `SAVE :: AllocationsDone, ...` but
//!   never guards the allocation, so a second in-process file read aborts;
//!   the wrapper allocates only when not yet associated (the author's
//!   evident intent) — the defined first-call path is statement-identical.
//! - Allocation failure calls the same `Fatal('DataToFieldSolver', ...)`
//!   via the `STAT=istat` checks, instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies, exactly
//!   as in the original (received and ignored).
//! - No Rust threads: the solver entries run serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn dts_begin(h_solver: *mut c_void);
    fn dts_info_start();
    fn dts_resolve();
    fn dts_has_pointfile() -> c_int;
    fn dts_ascii_points();
    fn dts_use_logscale() -> c_int;
    fn dts_continue_bc() -> c_int;
    fn dts_default_initialize();
    fn dts_bulk_assembly();
    fn dts_default_finish_bulk();
    fn dts_boundary_assembly();
    fn dts_default_finish();
    fn dts_default_dirichlet();
    fn dts_revert_log() -> c_int;
    fn dts_solve_revert() -> c_double;
    fn dts_solve() -> c_double;
    fn dts_info_done();
    fn dts_init_begin(h_solver: *mut c_void);
    fn dts_init_body();
}

/// `datatofieldsolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `DataToFieldSolver_init` (point-data exported variables,
/// default `Variable`, time derivative order, linear-system defaults).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "datatofieldsolver_init_"]
pub unsafe extern "C" fn datatofieldsolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        dts_init_begin(solver);
        dts_init_body();
    }
}

/// `datatofieldsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`datatofieldsolver_init_entry`].
#[export_name = "datatofieldsolver_"]
pub unsafe extern "C" fn datatofieldsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        dts_begin(solver);
        dts_info_start();
        dts_resolve();
        if dts_has_pointfile() != 0 {
            dts_ascii_points();
        }
        dts_use_logscale();
        let continue_bc = dts_continue_bc();
        dts_default_initialize();
        dts_bulk_assembly();
        dts_default_finish_bulk();
        if continue_bc != 0 {
            dts_boundary_assembly();
        }
        dts_default_finish();
        dts_default_dirichlet();
        if dts_revert_log() != 0 {
            let _norm: c_double = dts_solve_revert();
        } else {
            let _norm: c_double = dts_solve();
        }
        dts_info_done();
    }
}
