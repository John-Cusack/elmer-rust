//! Port of `fem/src/modules/ElementStats.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ElementStats.F90` — one F77-global entry
//! (`elementstats_`) with three `CONTAINS` internals. The emitted dynamic
//! symbol is exactly `elementstats_` (proven by `nm -D` on the reference
//! `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/elementstats_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the three element studies, `CEILING` binning,
//! `MIN`/`MAX` reductions, `ACOS`/`SQRT`) and every report (list-directed
//! `PRINT *`, `ES`-format `Info`, `Warn`/`Fatal`) stays in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors ElementStats.F90 statement for statement):
//! - The cohort second pass transcribes the original's `GOTO 100`
//!   re-sweep (including the repeated NaN reports) as an equivalent loop.
//! - The parallel-reduction branch is transcribed as written; in serial
//!   builds `IsParallel` is false by the same guard as the original.
//! - The studies' `SAVE Nodes` live as wrapper-side SAVE state (one set
//!   per study, consumed before the next fill, as in the original).
//! - The cohort table is allocated once per call when histograms are
//!   requested and released at the end, matching the original's
//!   per-call allocate/deallocate pair; allocation failure calls the same
//!   `Fatal('ElementStats', ...)`.
//! - The `Warn`/`Fatal` paths (`Element Warn Size`, `Stop At Nans`) run
//!   inside the wrappers with identical messages; a `Fatal` aborts the
//!   process exactly as in the original.
//! - `Model`, `dt`, and the transient flag are unused by the entry body
//!   itself (received and ignored), exactly as in the original.
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn est_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn est_bc_max() -> c_int;
    fn est_pass(is_bc: c_int);
    fn est_finish();
}

/// `elementstats_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer object addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "elementstats_"]
pub unsafe extern "C" fn elementstats_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        est_begin(model, solver);
        let bc_max = est_bc_max();
        for is_bc in 0..=bc_max {
            est_pass(is_bc);
        }
        est_finish();
    }
}
