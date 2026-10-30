//! Linear elasticity kernels for the Rust-native Elmer rewrite: the small-strain
//! path of `ElasticSolve.F90` that `StressSolve` delegates to
//! (`StressSolve.F90:87`).
//!
//! Covers isotropic and anisotropic (orthotropic through a full or diagonal
//! `C`) Hooke laws in 2D plane strain / plane stress, axisymmetric and 3D, body
//! force, `Stress Pressure`, stress/strain loads, thermal strain, mass and
//! damping, and the traction, normal-traction and spring boundary conditions.
//!
//! The solver evaluates every property through `elmer-expr` (tables, MATC,
//! procedures), at the element nodes or, with `Youngs Modulus at IP`,
//! `Poisson Ratio at IP`, `Heat Expansion Coefficient IP` and
//! `Stress Bodyforce at IP`, at the integration points. A transient run is
//! second order in time (`Default2ndOrderTime`, Newmark/Bossak) unless the sif
//! sets `Time Derivative Order = 1`, with `Damping` and `Rayleigh Damping
//! alpha/beta` in either. `Calculate Stresses/Strains/Principal` project
//! the elastic strain and stress onto the nodes ([`post`]).
//!
//! Large deflection, incompressible, Maxwell, UMAT, GPA, `Stress Load` /
//! `Strain Load` and stability terms are not implemented; callers must refuse a
//! case that asks for them rather than approximate it.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
pub mod material;
pub mod post;
mod props;
pub mod solver;

pub use boundary::{stress_boundary_local, stress_boundary_local_p, Spring, StressBoundaryInput};
pub use bulk::{shape_at_points, stress_bulk_local, stress_bulk_local_p, CoordSystem, IpProps, StressBulkInput};
pub use error::StressError;
pub use material::{diagonal_c, ElasticLaw};
pub use solver::{install_fe, register, FeServices, StressSolver};

/// Registry entry: `Procedure = "StressSolve" "StressSolver"`.
pub const REGISTRY_ENTRY: (&str, &str) = ("StressSolve", "StressSolver");
