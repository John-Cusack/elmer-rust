//! Mesh utility solvers of Elmer's `fem/src/modules`: `MeshChecksum`
//! (MeshChecksum.F90), `RigidMeshMapper` (RigidMeshMapper.F90),
//! `StructuredMeshMapper` (StructuredMeshMapper.F90) and
//! `StructuredProjectToPlane` (StructuredProjectToPlane.F90).
//!
//! [`register`] binds all four to their `Procedure = "Lib" "Entry"` names.
#![forbid(unsafe_code)]

mod checksum;

pub use checksum::{checksum, pseudo_norm, MeshChecksum};

use elmer_solver::SolverRegistry;

/// Registers every solver of this crate.
pub fn register(registry: &mut SolverRegistry) {
    registry.register("MeshChecksum", "MeshChecksum", checksum::factory);
}
