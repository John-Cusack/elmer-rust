//! Mesh operations of the Rust-native Elmer rewrite, numbered exactly as the
//! Fortran numbers them so results compare node by node.
//!
//! - [`increase_element_order`] ports `IncreaseElementOrder`
//!   (MeshBasics.F90:3564) for the linear elements 202, 303, 404, 504, 706
//!   and 808.
//!
//! Operations take any [`elmer_core::MeshView`] and return an owned
//! [`Mesh`] (the one of `elmer-mesh-levels`, itself a `MeshView`).

#![forbid(unsafe_code)]

mod error;
pub mod order;

pub use elmer_mesh_levels::Mesh;
pub use error::MeshOpsError;
pub use order::increase_element_order;
