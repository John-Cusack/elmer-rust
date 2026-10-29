//! Mesh multiplication and in-solver extrusion for the Rust-native Elmer
//! rewrite, numbered exactly as the Fortran numbers them so results compare
//! node by node.
//!
//! - [`split_mesh_equal`] ports `SplitMeshEqual` (MeshSplit.F90:52) for the
//!   linear elements 101, 202, 303, 404, 504, 706 and 808, with the edge and
//!   face numbering of `FindMeshEdges2D`, `FindMeshFaces3D` and
//!   `FindMeshEdges3D` (MeshBasics.F90:3875-4934) in [`topology`].
//! - [`mesh_levels`] is the Simulation `Mesh Levels` loop of `LoadModel`
//!   (ModelDescription.F90:3049-3116).
//! - [`extrude`] ports `MeshExtrude` (MeshExtrusion.F90:425) driven by the
//!   Simulation `Extruded Mesh *` keywords ([`ExtrudeOptions`]), including
//!   `ExtrudedDivision`, `GeometricUnitDivision` and `FunctionUnitDivision`
//!   in [`division`].
//!
//! Tags follow the Fortran's in-memory values after `MapBodiesAndBCs`
//! (MeshLoad.F90:330-455): a bulk tag is the body index and a boundary tag
//! the BC index (0 when no BC targets it). Splitting copies tags, so it works
//! on file tags too; extrusion drops boundary elements tagged 0 and numbers
//! the new boundaries as BC indices, so it needs the mapped tags.
//!
//! Not ported (refused): `Mesh Keep Grading`, `Follow * Boundary` curved
//! boundaries, pyramids, higher-order elements, parallel meshes,
//! `Extruded Mesh Slices`, `Extruded Mesh Rotational`, `Extruded Mesh End
//! Map`, `Extruded Elements`, rotor skew, and the `PrepareMesh` options that
//! rework the new mesh (`Split Mesh Quads`, `Increase Element Order`,
//! `Constant Stencil`). Solver-specific `Mesh Levels` and extrusion are the
//! caller's choice of keyword list.

#![forbid(unsafe_code)]

pub mod division;
mod error;
pub mod extrude;
mod mesh;
pub mod split;
pub mod topology;

pub use error::MeshLevelsError;
pub use extrude::{extrude, ExtrudeOptions};
pub use mesh::{Element, Mesh};
pub use split::{mesh_levels, split_mesh_equal};
