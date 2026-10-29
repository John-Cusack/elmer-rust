//! Errors of the thermoelectric kernels.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum ThermoElectricError {
    /// Something the Fortran supports and this crate refuses, not approximates.
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("{what}: expected {expected} nodal values, got {got}")]
    NodalLength {
        what: &'static str,
        expected: usize,
        got: usize,
    },
    /// The element basis could not be evaluated.
    #[error("element basis: {0}")]
    Basis(String),
}
