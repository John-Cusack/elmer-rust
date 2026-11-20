//! `ShellSolver` (ShellSolver.F90): Reissner-Mindlin shells of the
//! curvature-adapted formulation. Each background quadrilateral is
//! replaced by a blending surface through its edges' Hermite curves (built
//! from the element directors), reparametrised along the lines of
//! curvature with a Taylor polynomial; the element matrix uses the reduced
//! membrane (curl kernel with edge DOFs) and MITC shear strains and the
//! large-deflection terms, solved incrementally as in the Fortran.
//!
//! Ported: four-node quadrilaterals with the directors of `mesh.director`
//! (`ReadSurfaceDirector`) or else the face normals, `Large Deflection` (default on, its nonlinear
//! iterations and residuals), `Nonlinear Bending Strains`, `Normal
//! Pressure`, the dead edge loads of `ShellBoundaryMatrix` (`Resultant
//! Force`, `Resultant Couple`, `Spring`), `Displace Mesh`, `Shear /
//! Membrane / Stretch Relaxation Alpha`, `Benchmark Problem` and `Spherical
//! Body`.
//!
//! Refused: triangles and other element types, p elements (`Cartesian
//! Formulation`), a `Director` field, beams, live edge loads, the solid
//! coupling, `Drilling
//! DOFs`, `Rotate DOFs`, bubbles, explicit strain reduction operators,
//! mass (transient, harmonic, eigen) and more than one partition.
#![forbid(unsafe_code)]
// The ports index several parallel arrays with the Fortran loop variables.
// The expressions keep the Fortran's shape (and so its rounding order).
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::assign_op_pattern,
    clippy::type_complexity,
    clippy::manual_memcpy
)]

mod geom;
mod kernel;
pub mod motion;
mod solver;

pub use solver::{factory, register, solver_defaults, ShellSolver, PROCEDURE};

#[derive(Debug, thiserror::Error)]
pub enum ShellError {
    #[error("{0}")]
    Unsupported(String),
    #[error("{0}")]
    Fatal(String),
}
