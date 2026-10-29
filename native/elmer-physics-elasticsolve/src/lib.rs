//! `ElasticSolve` for the Rust-native Elmer rewrite: `Procedure = "ElasticSolve"
//! "ElasticSolver"` (`ElasticSolve.F90`), as an `elmer_solver::Solver`.
//!
//! Covers the small-strain and geometrically nonlinear (`Large Deflection`, the
//! default; St. Venant-Kirchhoff) linear elastic material in 2D plane strain /
//! plane stress, axisymmetric (`Axi Symmetric`, `Cylindric Symmetric`) and 3D
//! Cartesian, isotropic or anisotropic; the compressible neo-Hookean material
//! (`Neo-Hookean Material`); body, inertial and pressure loads, `Stress Load` /
//! `Strain Load`, thermal strain; the traction, normal-traction and spring
//! boundary conditions with `Pseudo-Traction`; and the Newton loop over
//! `Nonlinear System Max Iterations` that all of that needs. Transient runs are
//! first order (`Time Derivative Order = 1`, the damping matrix carries `d/dt`)
//! or second order (the default).
//!
//! `StressSolve` is the small-strain solver in `elmer-physics-elastic`; the
//! kernels there are not reused because `ElasticSolve` runs one Newton-form
//! assembly for every kinematics ([`bulk`]).
//!
//! Stress and strain postprocessing is done for the small-strain material only.
//! Incompressible and mixed formulations, Maxwell and UMAT materials, geometric
//! stiffness, eigen analysis, model lumping and staged construction are not
//! ported; callers get [`elmer_solver::SolverError::Unsupported`] for a case that
//! asks for them rather than an approximation.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
mod law;
pub mod neohooke;
mod props;
pub mod solver;
pub mod tensor;

pub use boundary::{elastic_boundary_local, node_uvw, BoundaryInput, ParentInput};
pub use bulk::{elastic_bulk_local, BulkInput};
pub use error::ElasticError;
pub use neohooke::{neo_hookean_local, NeoHookeInput};
pub use solver::{install_fe, register, ElasticSolver, FeServices};

/// Registry entry: `Procedure = "ElasticSolve" "ElasticSolver"`.
pub const REGISTRY_ENTRY: (&str, &str) = ("ElasticSolve", "ElasticSolver");
