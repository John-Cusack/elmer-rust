//! Mesh utility solvers of Elmer's `fem/src/modules`: `MeshChecksum`
//! (MeshChecksum.F90), `RigidMeshMapper` (RigidMeshMapper.F90),
//! `StructuredMeshMapper` (StructuredMeshMapper.F90) and
//! `StructuredProjectToPlane` (StructuredProjectToPlane.F90).
//!
//! [`register`] binds all four to their `Procedure = "Lib" "Entry"` names.
#![forbid(unsafe_code)]
// The ports index several parallel node arrays with the Fortran loop variable.
#![allow(clippy::needless_range_loop)]

mod checksum;
pub mod motion;
mod project;
mod rigid;
mod structure;
mod structured;
mod util;

pub use checksum::{checksum, pseudo_norm, MeshChecksum};
pub use project::StructuredProjectToPlane;
pub use rigid::RigidMeshMapper;
pub use structure::{detect_extruded_structure, Structure, Want};
pub use structured::StructuredMeshMapper;

use elmer_solver::SolverRegistry;

/// Registers every solver of this crate.
pub fn register(registry: &mut SolverRegistry) {
    registry.register("MeshChecksum", "MeshChecksum", checksum::factory);
    registry.register("RigidMeshMapper", "RigidMeshMapper", rigid::factory);
    registry.register(
        "StructuredProjectToPlane",
        "StructuredProjectToPlane",
        project::factory,
    );
    registry.register(
        "StructuredMeshMapper",
        "StructuredMeshMapper",
        structured::factory,
    );
}
