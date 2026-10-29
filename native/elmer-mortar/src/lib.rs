//! Periodic and mortar boundary coupling for the Rust-native Elmer rewrite.
//!
//! Ports the serial, nodal part of Elmer's boundary projectors:
//!
//! - [`create_interface_meshes`] (`CreateInterfaceMeshes`, MortarUtils.F90)
//!   cuts the two boundaries of a `Periodic BC` / `Mortar BC` pair out of the
//!   parent mesh as [`InterfaceMesh`]es, and [`map_interface_coordinate`]
//!   applies `Projector Coordinate Mapping`.
//! - [`overlay_interface_meshes`] (`OverlayIntefaceMeshes`) moves the target
//!   boundary onto this boundary with `Periodic BC Matrix`, or
//!   `Periodic BC Rotate`, `Periodic BC Scale` and `Periodic BC Translate`
//!   and their bounding-box defaults.
//! - [`nodal_projector`] (`NodalProjector`, via `MeshProjector` and
//!   `elmer-interp`) and [`periodic_projector`] (`PeriodicProjector`) build
//!   the projector matrix; [`generate_periodic_projectors`]
//!   (`GeneratePeriodicProjectors`) does it for every `Periodic BC`.
//! - [`save_projector`] writes the `Save Projector` files (`SaveProjector`).
//! - [`add_projector_topology`] (`ProjectorTopology`, ElementUtils.F90) adds
//!   the projector couplings to the matrix graph, and
//!   [`set_periodic_boundaries`] (`SetPeriodicBoundariesPass1/2`,
//!   BoundaryConditionUtils.F90) applies the constraint to an assembled
//!   system before the Dirichlet conditions.
//! - [`PeriodicAssembly`] is `elmer-solver`'s default assembly backend with
//!   both of those wired in: a `Simulation` with it in
//!   `backends.assembly` runs `Periodic BC` cases unchanged.
//!
//! For `Mortar BC` (Lagrange multipliers):
//!
//! - [`weighted_projector`] (`WeightedProjector`, InterpolateMeshToMesh.F90)
//!   is the Galerkin projector, and [`generate_mortar_projectors`]
//!   (`GenerateProjectors`, ProjectorUtils.F90) builds one per `Mortar BC`.
//! - [`generate_constraint_matrix`] (`GenerateConstraintMatrix`) stacks them
//!   into the solver's constraint matrix, one block per projector or, with
//!   `Mortar BCs Additive`, summed per node (`Projector Priority` choosing
//!   the owner of shared nodes); [`save_constraint_matrix`] writes it as
//!   `Save Constraint Matrix` does.
//! - [`collection_system`] (`SolveWithLinearRestriction`) builds the
//!   saddle-point system `[A C^T; C 0]` the linear solver then solves for
//!   the field and the multipliers.
//!
//! Refused with [`MortarError::Unsupported`]: level projectors, the
//! rotational, cylindrical, radial, axial, sliding, flat, plane and normal
//! projector families, discontinuous boundaries, Galerkin projectors with
//! edge DOFs or axisymmetric weighting, `Conforming BC`, `Periodic BC Rotate
//! Automatic`, explicit and Lagrange coefficient periodic conditions,
//! `Periodic BC Coefficient` jumps; for the constraints, vector fields,
//! `Mortar Diag` and mortar jumps, contact and
//! integral BCs, and eliminated, scaled or residual-mode restricted
//! systems. Parallel halo handling is out of scope.

#![forbid(unsafe_code)]

pub mod backend;
pub mod constrain;
pub mod constraint;
pub mod interface;
pub mod overlay;
pub mod projector;
pub mod restriction;
pub mod save;
pub mod weighted;

pub use backend::PeriodicAssembly;
pub use constrain::{
    add_projector_topology, periodic_active, set_periodic_boundaries, PeriodicTarget,
};
pub use constraint::{generate_constraint_matrix, ConstraintMatrix, ConstraintTarget};
pub use interface::{create_interface_meshes, map_interface_coordinate, InterfaceMesh};
pub use overlay::overlay_interface_meshes;
pub use projector::{
    generate_mortar_projectors, generate_periodic_projectors, nodal_projector, periodic_projector,
    save_request, GeneratedProjectors, Projector, ProjectorEnv, SaveRequest,
};
pub use restriction::{collection_system, CollectionSystem, Restriction};
pub use save::{save_constraint_matrix, save_projector};
pub use weighted::weighted_projector;

#[derive(Debug, thiserror::Error)]
pub enum MortarError {
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("{0}")]
    Fatal(String),
    #[error(transparent)]
    Interp(#[from] elmer_interp::InterpError),
    #[error("matrix has no entry ({row}, {col}) for the periodic constraint")]
    MissingEntry { row: usize, col: u32 },
}
