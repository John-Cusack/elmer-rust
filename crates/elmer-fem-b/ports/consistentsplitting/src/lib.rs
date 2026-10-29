//! Port of `fem/src/modules/ConsistentSplitting/` (Elmer @ c0b4bc6) for lane
//! T023 (rs-elmer-next Slice 1): whole-library Rust cdylib replacing
//! `ConsistentSplitting.so`, with the same dynamic surface (`nm -D`
//! equality over the four F77 entries).
//!
//! Boundary (solver-plugin, entry shape): each entry's derived dummies
//! (`Model`, `Solver`) cross as opaque addresses. Every component access
//! and every Elmer library call goes through the generated-style Fortran
//! `BIND(C)` wrappers in `accessors/cs_*_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MATMUL`, `DOT_PRODUCT`) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes per entry live atop each submodule. Shared rules:
//! - `SAVE` permanent storage is replicated as Rust `static Mutex<State>`
//!   with first-call sizing from the same queries (`MaxElementDOFs`,
//!   `NumberOfRows`); allocation failure calls the same `Fatal` via
//!   `try_reserve`, never aborting differently.
//! - Dead declarations the original never allocates (`PVelo`/`ppres` in
//!   CompressibilitySolver, `LOAD` there) are documented, not replicated.
//! - Literally dead code (`IF (.FALSE.)` blocks, commented-out calls) is
//!   documented at the call site and omitted.
//! - `Model`, `dt`, `TransientSimulation` are received and ignored exactly
//!   where the original ignores them.
//! - No Rust threads: solver entries run serially; each `Mutex` is only a
//!   `static`-sharing guard.
//! - The campaign builds serially here (`ParEnv % PEs == 1`); parallel arms
//!   are transcribed and selected by the same runtime queries, so a
//!   parallel-tree differential would drive them — staged proofs run the
//!   serial tree the campaign pins.

mod comp;
mod opt;
mod pres;
mod velo;

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

/// `compressibilitysolver_` — entry for CompressibilitySolver.F90.
#[export_name = "compressibilitysolver_"]
pub unsafe extern "C" fn compressibility_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { comp::entry(model, solver, dt, transient) }
}

/// `pressuresolver_` — entry for PressureSolver.F90.
#[export_name = "pressuresolver_"]
pub unsafe extern "C" fn pressure_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { pres::entry(model, solver, dt, transient) }
}

/// `velocitysolver_` — entry for VelocitySolver.F90.
#[export_name = "velocitysolver_"]
pub unsafe extern "C" fn velocity_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { velo::entry(model, solver, dt, transient) }
}

/// `optimalsolutionupdate_` — entry for OptimalSolutionUpdate.F90.
#[export_name = "optimalsolutionupdate_"]
pub unsafe extern "C" fn optimal_solution_update_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe { opt::entry(model, solver, dt, transient) }
}
