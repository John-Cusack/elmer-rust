//! FluxSolver: nodal flux and gradient recovery by Galerkin projection, the
//! port of Elmer's `FluxSolver` (`fem/src/modules/FluxSolver.F90`).
//!
//! For a scalar target field `u` (Temperature by default) the solver projects
//! `-C grad u` (`Calculate Flux`), `grad u` (`Calculate Grad`) and the pointwise
//! magnitudes `|C grad u|` and `|grad u|` (`Calculate Flux Magnitude`,
//! `Calculate Grad Magnitude`) onto the nodal basis. Every projected component
//! solves the same mass-matrix system `M x = b_i` with its own right-hand side,
//! by default with CG and ILU0 to `1e-10`.
//!
//! [`flux_bulk_local`] is the element kernel: a pure function of geometry and
//! nodal values, usable without a mesh or a solver context. [`FluxSolver`]
//! wires it into `elmer-solver`'s `Solver` trait and [`register`] binds it to
//! `Procedure = "FluxSolver" "FluxSolver"`.
//!
//! Supported: all four field requests, `Flux Abs` and `Grad Abs`, `Target
//! Variable`, `Flux Variable` and `Target Expression` (a material keyword
//! evaluated at the nodes), `Flux Coefficient` as scalar, diagonal or full
//! tensor material values, `Enforce Positive Magnitude`, axisymmetric weighting
//! and the solver's default linear-system settings. `Discontinuous Galerkin`
//! (elemental DG variables and face terms) is refused rather than approximated.
#![forbid(unsafe_code)]

mod error;
mod kernel;
mod solver;

pub use error::FluxError;
pub use kernel::{flux_bulk_local, FluxCoefficient, FluxConfig, FluxLocal};
pub use solver::{factory, install_fe, register, FeProvider, FluxSolver};

/// `Procedure = "FluxSolver" "FluxSolver"`.
pub const PROCEDURE: (&str, &str) = ("FluxSolver", "FluxSolver");
