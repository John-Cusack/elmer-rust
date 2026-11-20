//! `ShellSolver` (ShellSolver.F90): Reissner-Mindlin shells of the
//! curvature-adapted formulation. Each background quadrilateral is
//! replaced by a blending surface through its edges' Hermite curves (built
//! from the element directors), reparametrised along the lines of
//! curvature with a Taylor polynomial; the element matrix uses the reduced
//! membrane (curl kernel with edge DOFs) and MITC shear strains and the
//! large-deflection terms, solved incrementally as in the Fortran.
//!
//! Ported: four-node quadrilaterals whose directors are the face normals
//! (no `mesh.director`), `Large Deflection` (default on, its nonlinear
//! iterations and residuals), `Nonlinear Bending Strains`, `Normal
//! Pressure`, `Shear / Membrane / Stretch Relaxation Alpha`, `Benchmark
//! Problem` and `Spherical Body`.
//!
//! Refused: triangles and other element types, p elements (`Cartesian
//! Formulation`), `mesh.director` and `Director` fields, beams, edge loads
//! (`ShellBoundaryMatrix`), the solid coupling, `Displace Mesh`, `Drilling
//! DOFs`, `Rotate DOFs`, bubbles, explicit strain reduction operators,
//! mass (transient, harmonic, eigen) and more than one partition.
#![forbid(unsafe_code)]
// The ports index several parallel arrays with the Fortran loop variables.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

mod geom;
mod kernel;
mod solver;

pub use solver::{factory, register, solver_defaults, ShellSolver, PROCEDURE};

#[derive(Debug, thiserror::Error)]
pub enum ShellError {
    #[error("{0}")]
    Unsupported(String),
    #[error("{0}")]
    Fatal(String),
}
