//! `elmer-mesh`: Rust-native Elmer mesh DB reader and load transforms.
//!
//! Reads the Elmer mesh database (`mesh.header`, `mesh.nodes`,
//! `mesh.elements`, `mesh.boundary`) written by ElmerGrid/ElmerSolver,
//! with the element-type table generated from the `elements.def`
//! snapshot. Implements `elmer_core::MeshView` and the S1 load-time
//! transforms (`MapCoordinates`, target mapping, DefaultIntExtBCs,
//! 2D IntersectionBCs, IncreaseElementOrder).

pub mod adapters;
pub mod edges;
pub mod elements;
pub mod mesh;
pub mod transform;

pub use mesh::{load, write_db, BoundElem, EdgeTable, Elem, Mesh, MeshError};
pub use transform::transform;
