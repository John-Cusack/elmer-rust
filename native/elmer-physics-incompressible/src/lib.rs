//! IncompressibleNSSolver for the Rust-native Elmer rewrite: the incompressible
//! Navier-Stokes equations with equal-order velocity and pressure held stable
//! by `Pressure Stabilization` (`fem/src/modules/IncompressibleNS.F90`, also
//! reached as `IncompressibleNSVec`), and on p elements with the dofs from
//! the solver's `SolverDofs`: the Taylor-Hood pair (`Element = "p:2"`) and,
//! the default MINI element of `_Init0` (`p:1 b:N`, the velocity bubbles
//! condensed by `LCondensate`, with their history in a transient run).
//!
//! [`incompressible_bulk_local`] and [`incompressible_boundary_local`] are
//! pure functions of the element geometry and nodal values; element sizes
//! come from `elmer-fe` (`StabParam`), as do the default basis and Gauss
//! rules, and the boundary integrals are `elmer-physics-flow`'s.
//! [`IncompressibleNsSolver`] drives them through `elmer-solver`'s `Solver`
//! trait; [`register`] binds it to `Procedure = "IncompressibleNS"
//! "IncompressibleNSSolver"` and `"IncompressibleNSVec" ...`.
//!
//! Compressibility models, non-newtonian and turbulent viscosity, axisymmetric coordinates,
//! friction and wall-law boundaries and the block preconditioner are
//! refused, not approximated.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
pub mod solver;

pub use boundary::{incompressible_boundary_local, IncompressibleBoundaryInput};
pub use bulk::{
    incompressible_bulk_local, incompressible_bulk_local_with, lcondensate, lcondensate_bdf1, pin_pressure, HMode,
    IncompressibleBulkInput, PressureStab,
};
pub use error::IncompressibleError;
pub use solver::{
    factory, install_fe, register, solver_defaults, solver_defaults_with, FeProvider, IncompressibleNsSolver,
};

/// Registry entries: `Procedure = "IncompressibleNS" "IncompressibleNSSolver"`
/// and the front `"IncompressibleNSVec" "IncompressibleNSSolver"`.
pub const REGISTRY_ENTRIES: [(&str, &str); 2] = [
    ("IncompressibleNS", "IncompressibleNSSolver"),
    ("IncompressibleNSVec", "IncompressibleNSSolver"),
];
