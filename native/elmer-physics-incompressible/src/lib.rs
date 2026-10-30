//! IncompressibleNSSolver for the Rust-native Elmer rewrite: the incompressible
//! Navier-Stokes equations with equal-order velocity and pressure held stable
//! by `Pressure Stabilization` (`fem/src/modules/IncompressibleNS.F90`, also
//! reached as `IncompressibleNSVec`), and the Taylor-Hood pair of a steady
//! run on p elements without bubbles (`Element = "p:2"` serendipity
//! elements, the dofs from the solver's `SolverDofs`).
//!
//! [`incompressible_bulk_local`] and [`incompressible_boundary_local`] are
//! pure functions of the element geometry and nodal values; element sizes
//! come from `elmer-fe` (`StabParam`), as do the default basis and Gauss
//! rules, and the boundary integrals are `elmer-physics-flow`'s.
//! [`IncompressibleNsSolver`] drives them through `elmer-solver`'s `Solver`
//! trait; [`register`] binds it to `Procedure = "IncompressibleNS"
//! "IncompressibleNSSolver"` and `"IncompressibleNSVec" ...`.
//!
//! The default MINI element (`Element = "p:1 b:1"`, condensed velocity
//! bubbles), other bubble elements, compressibility models, non-newtonian and
//! turbulent viscosity, axisymmetric coordinates, friction and wall-law
//! boundaries and the block preconditioner are refused, not approximated.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bulk;
pub mod error;
pub mod solver;

pub use boundary::{incompressible_boundary_local, IncompressibleBoundaryInput};
pub use bulk::{
    incompressible_bulk_local, incompressible_bulk_local_with, pin_pressure, HMode, IncompressibleBulkInput,
    PressureStab,
};
pub use error::IncompressibleError;
pub use solver::{factory, install_fe, register, solver_defaults, FeProvider, IncompressibleNsSolver};

/// Registry entries: `Procedure = "IncompressibleNS" "IncompressibleNSSolver"`
/// and the front `"IncompressibleNSVec" "IncompressibleNSSolver"`.
pub const REGISTRY_ENTRIES: [(&str, &str); 2] = [
    ("IncompressibleNS", "IncompressibleNSSolver"),
    ("IncompressibleNSVec", "IncompressibleNSSolver"),
];
