//! Two scalar physics solvers ElmerGUI offers: `StatElecSolve` (electrostatic
//! potential) and `HelmholtzSolve` (complex-valued acoustic pressure), ports
//! of `fem/src/modules/StatElecSolve.F90` and `HelmholtzSolve.F90` onto
//! `elmer-solver`'s `Solver` trait.
//!
//! The element kernels are pure functions of geometry and nodal values, usable
//! without a mesh or a solver context: [`stat_elec_bulk_local`],
//! [`stat_elec_boundary_local`], [`stat_elec_post_local`],
//! [`helmholtz_bulk_local`] and [`helmholtz_boundary_local`].
//! [`StatElecSolver`] and [`HelmholtzSolver`] wire them into the solver
//! framework and [`register`] binds both to their `Procedure` names.
//!
//! StatElecSolve: real potential with `Relative Permittivity`, `Charge
//! Density`, `Electric Flux` / `Surface Charge Density` and the electrostatic
//! layer condition; the postprocessing fields `Electric Energy Density`,
//! `Electric Flux`, `Electric Field`, `Nodal Energy Density` and `Nodal
//! Electric Flux`, and `res: Electric Energy`, `res: Potential Difference`,
//! `res: Effective Capacitance`, `Energy Control` / `Charge Control` and
//! `Calculate Average Potential`.
//!
//! HelmholtzSolve: `Sound Speed`, `Density`, `Sound Damping`, `Pressure
//! Source`, `Wave Flux`, `Wave Impedance`, `Plane Wave BC`, axisymmetric
//! coordinates and `Pressure 1` / `Pressure 2` Dirichlet values. The complex
//! system is assembled as the equivalent real one with interleaved real and
//! imaginary parts, which the framework solves; see the crate README for what
//! `Linear System Complex` still needs from `elmer-solver`.
//!
//! Refused rather than approximated (see the README): complex permittivity,
//! `Electric Infinity BC` (needs the boundary normal), capacitance matrices,
//! piezoelectric coupling, elemental (DG) fields, mesh adaptivity, bubbles,
//! mean-flow convection and the flow and structure interface conditions.
#![forbid(unsafe_code)]

mod common;
mod error;
mod helmholtz;
mod helmholtz_solver;
mod statelec;
mod statelec_solver;

pub use common::{install_fe, FeProvider};
pub use error::ScalarError;
pub use helmholtz::{
    helmholtz_boundary_local, helmholtz_bulk_local, HelmholtzBcAtNodes, HelmholtzBulkAtNodes, HelmholtzConfig,
    HelmholtzLocal,
};
pub use helmholtz_solver::{helmholtz_factory, HelmholtzSolver};
pub use statelec::{
    stat_elec_boundary_local, stat_elec_bulk_local, stat_elec_post_local, DenseLu, StatElecAtNodes, StatElecBcAtNodes,
    StatElecConfig, StatElecPostLocal, PHYSICAL_EPS0,
};
pub use statelec_solver::{stat_elec_factory, StatElecSolver};

use elmer_solver::SolverRegistry;

/// `Procedure = "StatElecSolve" "StatElecSolver"`.
pub const STAT_ELEC_PROCEDURE: (&str, &str) = ("StatElecSolve", "StatElecSolver");
/// `Procedure = "HelmholtzSolve" "HelmholtzSolver"`.
pub const HELMHOLTZ_PROCEDURE: (&str, &str) = ("HelmholtzSolve", "HelmholtzSolver");

/// Registers both solvers.
pub fn register(registry: &mut SolverRegistry) {
    statelec_solver::register(registry);
    helmholtz_solver::register(registry);
}
