//! ThermoElectricSolver for the Rust-native Elmer rewrite: temperature and
//! electric potential solved as one strongly coupled system with the Seebeck,
//! Peltier and Joule couplings (`fem/src/modules/ThermoElectricSolver.F90`).
//!
//! [`thermoelectric_bulk_local`] and [`thermoelectric_boundary_local`] are
//! pure functions of the element geometry and nodal values, on any
//! `ElementBasis` (nodal Lagrange, or a p-hierarchic basis with interior
//! bubbles whose condensation [`condense_bubbles`] performs).
//! [`ThermoElectricSolver`] drives them through `elmer-solver`'s `Solver`
//! trait; [`register`] binds it to `Procedure = "ThermoElectricSolver"
//! "ThermoElectricSolver"`.
#![forbid(unsafe_code)]

pub mod boundary;
pub mod bubble;
pub mod bulk;
pub mod condense;
pub mod error;
pub mod solver;

pub use bubble::{parse_element, BubbleBasis, BubbleConfig};
pub use boundary::{thermoelectric_boundary_local, ThermoElectricBoundaryInput};
pub use bulk::{thermoelectric_bulk_local, Stabilization, Tensor, ThermoElectricBulkInput};
pub use condense::{condense_bubbles, BubbleRecovery};
pub use error::ThermoElectricError;
pub use solver::{factory, install_fe, register, ThermoElectricSolver};

/// Registry entry: `Procedure = "ThermoElectricSolver" "ThermoElectricSolver"`.
pub const REGISTRY_ENTRIES: [(&str, &str); 1] = [("ThermoElectricSolver", "ThermoElectricSolver")];
