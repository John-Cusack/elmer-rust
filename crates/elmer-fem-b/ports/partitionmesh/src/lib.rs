//! Port of `fem/src/modules/PartitionMesh.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/PartitionMesh.F90` — one F77-global
//! entry `partitionmeshsolver_` with two `CONTAINS` internals. The
//! emitted dynamic symbol is exactly `partitionmeshsolver_` (proven by
//! `nm -D` on the reference `.so`); the internals emit only local
//! symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/partitionmesh_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive loops (`CreateNeighbourList` accumulation),
//! formatted I/O (`WRITE(Message,...)`), and whole-array updates
//! (`Var % Values = 0`, `ElementPart` fills) stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors PartitionMesh.F90 statement for statement):
//! - The entry-statement order is preserved: start `Info`, the
//!   `ParEnv % PEs > 1` early return (with the original's `PartitionMes`
//!   caller-name typo), `PartitionMeshSerial`, the neighbour-list build,
//!   the `Output Directory` branch, the `Partitioning Exit` STOP branch,
//!   the create-variable `Info`, `SetPartitionVariable`, and the done
//!   `Info`.
//! - The `CHARACTER(LEN=MAX_NAME_LEN)` directory scratch and the
//!   `ListGetString` hidden-length call never cross the boundary: the
//!   wrapper reads the value and the Rust side only sees the `Found`
//!   flag, exactly like the original's `IF( Found )` branch.
//! - `SAVE` semantics replicated: the original's host-associated
//!   `NeighbourList` starts unassociated on every invocation, so the
//!   wrapper releases any previous SAVE allocation before the fresh
//!   `ALLOCATE` (identical lifetimes, no leak across repeated calls);
//!   `ElementPart` re-points at `Mesh % RePartition` after every
//!   `PartitionMeshSerial` call, as in the original.
//! - `Partitioning Exit` calls the same `STOP` via a wrapper (the
//!   `pma_stop` call never returns, as in the original); the process
//!   proofs never set that keyword.
//! - The dead `n = ...` bulk/boundary count in the entry body is a pure
//!   local fill with no observable effect and is not materialised.
//! - No Rust threads: the solver calls entries serially, and the Rust
//!   side holds no `static` state at all (all lifetimes live in the
//!   wrapper module), so no sharing guard is needed.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn pma_begin(h_solver: *const c_void, h_model: *const c_void);
    fn pma_info_start();
    fn pma_pes() -> c_int;
    fn pma_info_serial_only();
    fn pma_partition();
    fn pma_create_neighbour_list();
    fn pma_output_found() -> c_int;
    fn pma_write_partitioned();
    fn pma_exit_flag() -> c_int;
    fn pma_stop();
    fn pma_info_create_var();
    fn pma_set_partition_variable();
    fn pma_info_done();
}

/// `partitionmeshsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "partitionmeshsolver_"]
pub unsafe extern "C" fn partition_mesh_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { pma_begin(solver, model) };

    unsafe { pma_info_start() };

    // Currently the mesh can only be partitioned in serial!
    if unsafe { pma_pes() } > 1 {
        unsafe { pma_info_serial_only() };
        return;
    }

    unsafe { pma_partition() };
    unsafe { pma_create_neighbour_list() };

    if unsafe { pma_output_found() } != 0 {
        unsafe { pma_write_partitioned() };
    }

    if unsafe { pma_exit_flag() } != 0 {
        unsafe { pma_stop() };
        std::process::abort();
    }

    unsafe { pma_info_create_var() };
    unsafe { pma_set_partition_variable() };

    unsafe { pma_info_done() };
}
