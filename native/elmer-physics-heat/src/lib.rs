//! HeatSolve element kernels: the local matrices of Elmer's `HeatSolver`
//! (`fem/src/modules/HeatSolve.F90`, vectorised path `LocalMatrixVec` and
//! `LocalMatrixBC`) for one scalar temperature degree of freedom per node.
//!
//! The kernels are pure functions of element geometry and *nodal* property
//! values, so they run without a mesh, a solver context or a property
//! evaluator: whoever drives them (the `Solver` impl, or a test-local glue
//! loop) evaluates keywords at the element nodes, as Elmer's
//! `ListGetElementReal` does, and the kernels interpolate to integration
//! points with the element basis.
//!
//! Time integration is not applied here. Transient elements return the mass
//! matrix in [`LocalSystem::mass`] next to the stiffness and force, and the
//! time-integration lane forms `M/dt` itself (Elmer's `Default1stOrderTime`).
//!
//! [`HeatSolver`] wires the kernels into `elmer-solver`'s `Solver` trait and
//! [`register`] binds it to `Procedure = "HeatSolve" "HeatSolver"`; it
//! evaluates properties with `elmer-expr` each nonlinear iteration, so
//! temperature-dependent properties converge by Picard iteration.
//!
//! Supported: scalar, diagonal and full-tensor conductivity, constant or
//! computed convection velocity (given at nodes) with plain Galerkin or SUPG
//! weighting, volumetric and per-mass heat source, heat transfer multiplier,
//! axisymmetric weighting, heat flux and Robin (heat transfer coefficient and
//! external temperature) boundary terms. With `Legacy Assembly = True` (the
//! HeatSolveLegacy path, see `legacy`) also `Smart Heater Control` and
//! `Radiation = Idealized`. `Radiation = Diffuse Gray` on axisymmetric
//! two-node boundaries runs through `elmer-viewfactors` (view factors,
//! normalisation, Gebhart factors, computed on the first solve) and widens
//! the matrix pattern by the factors' couplings. Everything else HeatSolve can do is refused by
//! [`check_supported`] instead of being approximated; `Convection = Computed`
//! and `Compressibility Model` wait for the solve, so a solver with `Exec
//! Solver = Never` passes setup as it does in Elmer.
#![forbid(unsafe_code)]

mod boundary;
mod bulk;
mod error;
mod legacy;
mod material;
mod preflight;
mod radiation;
mod solver;

pub use boundary::heat_boundary_local;
pub use bulk::heat_bulk_local;
pub use error::HeatError;
pub use material::{
    Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, HeatMaterialAtNodes,
    Stabilization,
};
pub use preflight::{check_supported, HeatSections, PROCEDURE};
pub use solver::{factory, install_fe, register, FeProvider, HeatSolver};
