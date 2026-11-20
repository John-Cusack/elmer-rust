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
//! Supported: scalar, diagonal and full-tensor conductivity, constant
//! convection velocity or `Convection = Computed` (the velocity of the flow
//! solution, less any `Mesh Velocity`) with plain Galerkin or SUPG weighting,
//! `Compressibility Model = Thermal` and `Perfect Gas` (the density law and,
//! with computed convection, the pressure work of the flow pressure),
//! volumetric and per-mass heat source, heat transfer multiplier,
//! axisymmetric weighting, heat flux and Robin (heat transfer coefficient and
//! external temperature) boundary terms. With `Legacy Assembly = True` (the
//! HeatSolveLegacy path, see `legacy`) also `Smart Heater Control` and
//! `Radiation = Idealized` and the `Heat Gap` boundary condition (nodes across
//! the gap are found by position and coupled by `-STIFF`). `Radiation = Diffuse Gray` on axisymmetric
//! two-node boundaries runs through `elmer-viewfactors` (view factors,
//! normalisation, Gebhart factors, computed on the first solve) and widens
//! the matrix pattern by the factors' couplings. Everything else HeatSolve can do is refused by
//! [`check_supported`] instead of being approximated; the legacy path's
//! `Convection = Computed` and `Compressibility Model` wait for the solve, so
//! a solver with `Exec Solver = Never` passes setup as it does in Elmer.
#![forbid(unsafe_code)]

mod adaptive;
mod boundary;
mod bulk;
mod bulk_vec;
mod condense;
mod error;
mod info_vec;
mod legacy;
mod legacy_compose;
mod linear_forms;
mod material;
mod preflight;
mod radiation;
mod radiators;
mod solver;

pub use boundary::heat_boundary_local;
pub use bulk::heat_bulk_local;
pub use bulk_vec::heat_bulk_local_vec;
pub use error::HeatError;
pub use info_vec::ElementInfoVec;
pub use material::{
    CompressModel, Compressibility, Conductivity, HeatBcAtNodes, HeatBoundaryConfig,
    HeatBulkConfig, HeatMaterialAtNodes, JouleHeat, Stabilization,
};
pub use preflight::{check_supported, HeatSections, PROCEDURE};
pub use solver::{factory, install_fe, register, solver_defaults, FeProvider, HeatSolver};
