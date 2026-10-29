//! Magnetodynamics on Whitney edge elements: the steady-state branch of
//! `WhitneyAVSolver` (`fem/src/modules/MagnetoDynamics/WhitneyAVSolver.F90`,
//! [`steady`]) and `WhitneyAVHarmonicSolver` (`WhitneyAVHarmonicSolver.F90`,
//! [`harmonic`]) for the lowest-order (non-Piola) Whitney basis.
//!
//! The solver couples the vector potential `A` (one DOF per mesh edge) with
//! the electric scalar potential `V` (one DOF per node, only with `Static
//! Conductivity`, which a boundary `Electric Current Density` also turns on):
//!
//! ```text
//!   curl(nu curl A) + sigma grad V = J^s + curl M^s - sigma grad V^s
//!   -div(sigma grad V) = 0
//! ```
//!
//! The system layout follows Elmer's `Perm` without bandwidth optimisation:
//! the nodal DOFs of the active elements come first in node order, then one
//! DOF per edge in first-use order over the active elements. The steady
//! system is not gauged, so the `A` part of the solution depends on the
//! Krylov path; methods without a preconditioner see the same path under any
//! DOF permutation, an ILU preconditioner does not (see [`steady`]).
//!
//! Everything this port does not cover (Piola and second-order bases,
//! transient runs, imaginary material and source values, gauges, coils,
//! H-B curves, the `JFix`
//! source fixing, ...) is refused with [`MagnetoError::Unsupported`], never
//! approximated.
#![forbid(unsafe_code)]

pub mod case;
pub mod error;
pub mod harmonic;
pub mod layout;
pub mod local;
pub mod material;
pub mod solver;
pub mod steady;
pub mod topology;

pub use case::Case;
pub use error::MagnetoError;
pub use harmonic::{solve_harmonic, HarmonicSolution};
pub use solver::register;
pub use steady::{solve_steady, AvSolution};
