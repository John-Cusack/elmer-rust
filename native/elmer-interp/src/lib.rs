//! Mesh-to-mesh interpolation for the Rust-native Elmer rewrite.
//!
//! Ports the serial behaviour of `InterpolateMeshToMesh.F90`
//! (`InterpolateMeshToMeshQ`) and the point-location routines it calls in
//! `Interpolation.F90` (`BuildQuadrantTree`, `FindLeafElements`,
//! `PointInElement`) and `ElementGeometry.F90` (`GlobalToLocal`):
//!
//! - [`Locator`] finds the old-mesh bulk element and local coordinates of a
//!   point through a quadrant tree, with the Fortran's progressively sloppier
//!   tolerances and the exhaustive-scan fallback.
//! - [`interpolate`] moves same-named scalar variables (including previous
//!   time-step columns, DG variables, permutations and a body mask) from the
//!   old mesh to the new one.
//! - [`build_projector`] records the same interpolation as a sparse operator,
//!   with weights from the caller's [`elmer_core::ElementBasis`]; applying it
//!   later is [`Projector::apply`].
//!
//! Not ported (refused or out of scope): the partition exchange of
//! `InterpolateMeshToMesh` for `ParEnv % PEs > 1` (the MPI lane owns that),
//! edge-element (`Edge Basis`) projector rows, p-elements, and
//! `WeightedProjector`. Pyramids 613/614 and connectivity elements have no
//! nodal table and are rejected by [`Locator::new`].

#![forbid(unsafe_code)]

mod element_defs;
pub mod field;
pub mod geom;
pub mod locate;
pub mod octree;
pub mod options;
pub mod shape;
pub mod transfer;

pub use field::{DgIndexes, FieldData, Location};
pub use locate::{BodyMask, Located, Locator};
pub use options::InterpOptions;
pub use transfer::{build_projector, interpolate, Extras, InterpReport, Projector};

#[derive(Debug, thiserror::Error)]
pub enum InterpError {
    #[error("element type {0} has no nodal basis table")]
    UnsupportedElement(u16),
    #[error("{0} target nodes not found (Interpolation Abort Not Found)")]
    NotFound(usize),
    #[error("DG variable needs DG indexes for the {0} mesh")]
    MissingDgIndexes(&'static str),
    #[error("element basis failed: {0}")]
    Basis(String),
}
