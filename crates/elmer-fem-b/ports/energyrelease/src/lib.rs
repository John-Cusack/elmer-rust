//! Port of `fem/src/modules/EnergyRelease.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/EnergyRelease.F90` — one F77-global entry
//! (`releaseratesolver_`) with one `CONTAINS` internal. The emitted dynamic
//! symbol is exactly `releaseratesolver_` (proven by `nm -D` on the
//! reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/energyrelease_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the per-integration-point tensor assembly,
//! `MATMUL`s, `SUM` reductions) and every report (list-directed `PRINT *`)
//! stays in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors EnergyRelease.F90 statement for statement):
//! - The SAVE work arrays (`LocalDisplacement`, `LocalStress`,
//!   `LocalPropagationShape`, `AllocationsDone`) live as wrapper-side
//!   SAVE state: allocated once at first-mesh `MaxElementNodes` sizing,
//!   kept across invocations, exactly as in the original.
//! - The internal procedure's `SAVE Nodes` + `SAVE IP`-equivalent state
//!   lives as wrapper-side SAVE state (consumed before the next fill, as
//!   in the original).
//! - Allocation failure calls the same `Fatal('PoissonSolve', ...)` the
//!   original calls (bug-compatible program name) via a `STAT=` check.
//! - `Norm` is declared-but-unused in the original and absent here;
//!   `Model`, `dt`, and the transient flag are unused by the body,
//!   exactly as in the original (received and ignored).
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ere_begin(h_solver: *mut c_void);
    fn ere_run();
}

/// `releaseratesolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "releaseratesolver_"]
pub unsafe extern "C" fn releaseratesolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ere_begin(solver);
        ere_run();
    }
}
