//! Port of `fem/src/modules/FluidicForce.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/FluidicForce.F90` — one F77-global entry
//! (`forcecompute_`) with one `CONTAINS` internal. The emitted dynamic
//! symbol is exactly `forcecompute_` (proven by `nm -D` on the reference
//! `.so`); the internal emits only local symbols. Note the campaign
//! library-vs-entry name difference: the proof invokes the port via the
//! campaign `Procedure` form `"FluidicForce" "ForceCompute"`.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/fluidicforce_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the per-integration-point stress assembly,
//! `MATMUL`/`TRANSPOSE`, `SUM` reductions, `SQRT`) and every report (the
//! `Level=4` `Info` lines, the list-directed shear-stress data files)
//! stays in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FluidicForce.F90 statement for statement):
//! - The whole entry body plus its `ForceIntegrate` internal is
//!   transcribed statement-for-statement into `ffc_run` (the T021-A
//!   EnergyRelease shape): same allocations (`MaxElementNodes`,
//!   `NumberOfBCs`), same loop nest, same `WRITE` formats, same file
//!   units for the shear-stress output. There is no `SAVE` state in the
//!   original (all work arrays are per-call allocatables), so re-entry
//!   is trivially clean.
//! - `Model`, `dt`, and the transient flag are unused by the entry
//!   itself beyond the `Model` handle traffic, exactly as in the
//!   original (received and ignored); `Mesh => GetMesh()` is kept as a
//!   touch call.
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ffc_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn ffc_run();
}

/// `forcecompute_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model` and `solver` must be the Elmer `Model_t` / `Solver_t` addresses
/// the loader passes. `dt`/`transient` must be valid reads if non-null
/// (they are ignored).
#[export_name = "forcecompute_"]
pub unsafe extern "C" fn forcecompute_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ffc_begin(model, solver);
        ffc_run();
    }
}
