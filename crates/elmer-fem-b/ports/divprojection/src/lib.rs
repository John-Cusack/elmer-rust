//! Port of `fem/src/modules/DivProjection.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/DivProjection.F90` — one F77-global entry
//! (`divprojection_`), no `CONTAINS` internals. The emitted dynamic symbol
//! is exactly `divprojection_` (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses.
//!
//! Faithfulness notes (mirrors DivProjection.F90 statement for statement):
//! the original body is empty apart from `USE` clauses (a dummy solver
//! stub "for generating pressure mass matrix" whose implementation is a
//! commented-out print). The port body is therefore empty too: it
//! receives the four entry arguments and returns without touching them,
//! exactly as the original touches nothing. No Elmer calls, no `SAVE`
//! state (the original has none), no Rust threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

/// `divprojection_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; all four are unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer object addresses the loader passes
/// (never dereferenced); `dt`/`transient` must be valid reads if non-null
/// (they are ignored).
#[export_name = "divprojection_"]
pub unsafe extern "C" fn divprojection_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}
