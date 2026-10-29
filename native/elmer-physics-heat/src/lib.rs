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
//! [`register`] binds it to `Procedure = "HeatSolve" "HeatSolver"`.
//!
//! Supported: scalar, diagonal and full-tensor conductivity, constant or
//! computed convection velocity (given at nodes) with plain Galerkin or SUPG
//! weighting, volumetric and per-mass heat source, heat transfer multiplier,
//! axisymmetric weighting, heat flux and Robin (heat transfer coefficient and
//! external temperature) boundary terms. Everything else HeatSolve can do is
//! refused by [`check_supported`] instead of being approximated.
#![forbid(unsafe_code)]

mod boundary;
mod bulk;
mod error;
mod material;
mod preflight;
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
