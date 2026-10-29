//! Mesh operations of the Rust-native Elmer rewrite, numbered exactly as the
//! Fortran numbers them so results compare node by node.
//!
//! - [`increase_element_order`] ports `IncreaseElementOrder`
//!   (MeshBasics.F90:3564) for the linear elements 202, 303, 404, 504, 706
//!   and 808.
//! - [`create_discont_mesh`] ports `CreateDiscontMesh` (MeshBasics.F90:377),
//!   the mesh splitting of `Discontinuous Boundary`.
//! - [`read_target_names`] ports `ReadTargetNames` (MeshIO.F90:959), the
//!   `Use Mesh Names` targets that `DisContBoundaryDoubleNames` relies on.
//! - [`mesh2d`] ports the `Mesh2D` generator (meshgen2d) for the layer types
//!   of the adaptive-refinement quick tests.
//!
//! Operations take any [`elmer_core::MeshView`] and return an owned
//! [`Mesh`] (the one of `elmer-mesh-levels`, itself a `MeshView`).

#![forbid(unsafe_code)]

pub mod discont;
mod error;
pub mod mesh2d;
pub mod names;
pub mod order;

pub use discont::{create_discont_mesh, DiscontMesh};
pub use elmer_mesh_levels::Mesh;
pub use error::MeshOpsError;
pub use mesh2d::{mesh2d, Mesh2dFiles};
pub use names::read_target_names;
pub use order::increase_element_order;
