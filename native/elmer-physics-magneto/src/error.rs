//! Errors of the magnetodynamics solvers.

use elmer_core::LinalgError;
use elmer_solver::SolverError;

#[derive(Debug, thiserror::Error)]
pub enum MagnetoError {
    /// Bad or missing input (sif, mesh, keyword values).
    #[error("{0}")]
    Input(String),
    /// An option the Fortran supports that this port does not.
    #[error("not ported: {0}")]
    Unsupported(String),
    #[error(transparent)]
    Solver(#[from] SolverError),
    #[error(transparent)]
    Linalg(#[from] LinalgError),
}

impl MagnetoError {
    pub(crate) fn unsupported(msg: impl Into<String>) -> Self {
        Self::Unsupported(msg.into())
    }
}
