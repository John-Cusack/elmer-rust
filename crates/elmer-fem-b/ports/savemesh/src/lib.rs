//! Port of `fem/src/modules/SaveMesh.F90` (Elmer @ 9f6af2f85, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/SaveMesh.F90` — one F77-global entry
//! `savemesh_` with no `CONTAINS` internals. The emitted dynamic symbol
//! is exactly `savemesh_` (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/savemesh_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. The
//! `WRITE` formatting (`i4.4`) and the `ParEnv % PEs` serial branch stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors SaveMesh.F90 statement for statement):
//! - `Mesh => Solver % Mesh` is a wrapper-side pointer rebind on every
//!   call, as in the original.
//! - The `"Save Mesh"` flag check (including the `Warn` on missing flag
//!   and the `.FALSE.` default) runs first; `.FALSE.` returns with
//!   nothing to do, exactly as in the original.
//! - `TimestepVar => VariableGet(..., "Timestep", .TRUE.)`,
//!   `"Mesh Name"` / `"Save Mesh Directory"` lookups (including the
//!   `Fatal` exits), both `MakeDirectory` calls, the `"Save All
//!   Timesteps"` lookup with `.FALSE.` default, the `WRITE(MeshDir,
//!   ...)` filename formatting, and the serial `WriteMeshToDisk2` call
//!   all run inside the `svm_do_save` wrapper in original order.
//! - `Model`, `dt`, and `TransientSimulation` are unused by the body
//!   except as the `WriteMeshToDisk2` model argument, exactly as in the
//!   original (`Model` crosses as an opaque handle and is only handed
//!   back to Fortran).
//! - Parallel arm out of scope: the serial proof tree (`ParEnv % PEs <=
//!   1`, no MPI) always takes the serial branch, which the wrapper
//!   replicates. On `PEs > 1` this port is a documented no-op where the
//!   original would `MPI_BARRIER` + write partitioned meshes.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn svm_begin(h_model: *const c_void, h_solver: *const c_void);
    fn svm_save_flag() -> c_int;
    fn svm_do_save();
}

/// Entry-state guard: the wrappers hold the `Model`/`Solver` addresses
/// as module-level pointers (mirroring the original's dummy arguments),
/// so rebinds serialize through this lock. Never contended: Elmer calls
/// solver entries serially.
static ENTRY: Mutex<()> = Mutex::new(());

/// `savemesh_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "savemesh_"]
pub unsafe extern "C" fn savemesh_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = ENTRY.lock().unwrap();
    unsafe { svm_begin(model, solver) };

    // IF(.NOT. SaveMeshLogical) RETURN — nothing to do.
    if unsafe { svm_save_flag() } == 0 {
        return;
    }

    // Remainder of the body: Timestep lookup, name/directory strings,
    // MakeDirectory calls, EveryTime branch + WRITE formatting, and the
    // serial WriteMeshToDisk2 — in original order inside the wrapper.
    unsafe { svm_do_save() };
}
