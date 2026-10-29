//! Port of `fem/src/modules/SaveData/` (Elmer @ c0b4bc6) for lane T023
//! (rs-elmer-next Slice 1): whole-library Rust cdylib replacing
//! `SaveData.so`, with the same dynamic surface (`nm -D` equality over
//! the ten F77 entries).
//!
//! Boundary (solver-plugin, entry shape): each entry's derived dummies
//! (`Model`, `Solver`) cross as opaque addresses. Every component access
//! and every Elmer library call goes through the generated-style Fortran
//! `BIND(C)` wrappers in `accessors/sd_*_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset.
//!
//! Faithfulness notes per entry live atop each submodule. Shared rules:
//! - `SAVE` permanent storage is replicated as Rust `static Mutex<State>`
//!   with first-call sizing from the same queries; allocation failure
//!   calls the same `Fatal` via `try_reserve`.
//! - Dead declarations never allocated or touched by the source are
//!   documented, not replicated.
//! - Literally dead code (`IF (.FALSE.)` blocks, commented-out calls) is
//!   documented at the call site and omitted.
//! - `Model`, `dt`, `TransientSimulation` are received and ignored exactly
//!   where the original ignores them.
//! - No Rust threads: the solver entries run serially; each `Mutex` is only
//!   a `static`-sharing guard.
//! - The campaign builds serially here; parallel arms are transcribed and
//!   selected by the same runtime queries.

mod bvals;
mod depend;
mod line;
mod mats;
mod proj;
mod scalars;

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

/// `saveboundaryvalues_` — entry for SaveBoundaryValues.F90.
#[export_name = "saveboundaryvalues_"]
pub unsafe extern "C" fn saveboundaryvalues_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { bvals::entry(model, solver, dt, transient) }
}

/// `savedependence_` — entry for SaveDependence.F90.
#[export_name = "savedependence_"]
pub unsafe extern "C" fn savedependence_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { depend::entry(model, solver, dt, transient) }
}

/// `savedependence_init_` — init for SaveDependence.F90.
#[export_name = "savedependence_init_"]
pub unsafe extern "C" fn savedependence_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { depend::savedependence_init_entry(model, solver, dt, transient) }
}

/// `saveline_` — entry for SaveLine.F90.
#[export_name = "saveline_"]
pub unsafe extern "C" fn saveline_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { line::entry(model, solver, dt, transient) }
}

/// `saveline_init_` — init for SaveLine.F90.
#[export_name = "saveline_init_"]
pub unsafe extern "C" fn saveline_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { line::saveline_init_entry(model, solver, dt, transient) }
}

/// `savematerials_` — entry for SaveMaterials.F90.
#[export_name = "savematerials_"]
pub unsafe extern "C" fn savematerials_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { mats::entry(model, solver, dt, transient) }
}

/// `saveprojection_` — entry for SaveProjection.F90.
#[export_name = "saveprojection_"]
pub unsafe extern "C" fn saveprojection_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { proj::entry(model, solver, dt, transient) }
}

/// `saveprojection_init_` — init for SaveProjection.F90.
#[export_name = "saveprojection_init_"]
pub unsafe extern "C" fn saveprojection_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { proj::saveprojection_init_entry(model, solver, dt, transient) }
}

/// `savescalars_` — entry for SaveScalars.F90.
#[export_name = "savescalars_"]
pub unsafe extern "C" fn savescalars_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { scalars::entry(model, solver, dt, transient) }
}

/// `savescalars_init_` — init for SaveScalars.F90.
#[export_name = "savescalars_init_"]
pub unsafe extern "C" fn savescalars_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { scalars::savescalars_init_entry(model, solver, dt, transient) }
}
