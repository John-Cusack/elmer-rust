//! Linear elasticity kernels for the Rust-native Elmer rewrite: the small-strain
//! path of `ElasticSolve.F90` that `StressSolve` delegates to
//! (`StressSolve.F90:87`).
//!
//! Covers isotropic and anisotropic (orthotropic through a full or diagonal
//! `C`) Hooke laws in 2D plane strain / plane stress, axisymmetric and 3D, body
//! force, `Stress Pressure`, stress/strain loads, thermal strain, mass and
//! damping, and the traction, normal-traction and spring boundary conditions.
//! Large deflection, incompressible, Maxwell, UMAT, GPA and stability terms are
//! not implemented; callers must refuse a case that asks for them rather than
//! approximate it.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
pub mod material;
pub mod solver;

pub use boundary::{stress_boundary_local, Spring, StressBoundaryInput};
pub use bulk::{stress_bulk_local, CoordSystem, StressBulkInput};
pub use error::StressError;
pub use material::{diagonal_c, ElasticLaw};
pub use solver::{install_fe, register, FeServices, StressSolver};

/// Registry entry: `Procedure = "StressSolve" "StressSolver"`.
pub const REGISTRY_ENTRY: (&str, &str) = ("StressSolve", "StressSolver");
