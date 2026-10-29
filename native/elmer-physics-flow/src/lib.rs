//! FlowSolve for the Rust-native Elmer rewrite: the steady incompressible
//! Navier-Stokes equations with equal-order (P1/P1, Q1/Q1) velocity and
//! pressure and the residual-based stabilization of FlowSolve's default
//! `Stabilization Method = Stabilized` (`fem/src/modules/FlowSolve.F90`,
//! element matrices from `NavierStokesCompose` and `NavierStokesBoundary`
//! in `fem/src/NavierStokes.F90`).
//!
//! [`flow_bulk_local`] and [`flow_boundary_local`] are pure functions of the
//! element geometry and nodal values. [`FlowSolver`] drives them through
//! `elmer-solver`'s `Solver` trait with FlowSolve's nonlinear loop: Picard
//! iteration, switching to Newton after `Nonlinear System Newton After
//! Iterations` or `... Newton After Tolerance`. [`register`] binds it to
//! `Procedure = "FlowSolve" "FlowSolver"`.
//!
//! Supported: Stokes, no-convection and full Navier-Stokes flow models in
//! cartesian 2D and 3D, linear elements, constant or field-dependent density
//! and viscosity, body forces and hydrostatic gravity, traction, normal
//! pressure, mass flux and slip-coefficient boundaries, Dirichlet and
//! no-slip-wall conditions, and transient runs through `Default1stOrderTime`
//! (`elmer-timeint`). Bubbles, P2/P1, VMS, compressible models, porous
//! media, rotating frames, non-newtonian viscosity, normal-tangential and
//! free-surface boundaries are refused instead of approximated.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
pub mod solver;
pub mod stab;

pub use boundary::{flow_boundary_local, outward_normal, FlowBoundaryInput};
pub use bulk::{flow_bulk_local, Discretization, FlowBulkInput};
pub use error::FlowError;
pub use solver::{factory, install_fe, register, FeProvider, FlowSolver};
pub use stab::{element_diameter, stab_params, StabParams};

/// Registry entry: `Procedure = "FlowSolve" "FlowSolver"`.
pub const REGISTRY_ENTRY: (&str, &str) = ("FlowSolve", "FlowSolver");
