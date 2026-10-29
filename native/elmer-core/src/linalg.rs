//! Operator, preconditioner and solver traits shared by the linear-algebra lanes.

use crate::comm::Communicator;
use crate::csr::CsrMatrix;
use crate::scalar::Scalar;

#[derive(Debug, thiserror::Error)]
pub enum LinalgError {
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("breakdown in {method}: {detail}")]
    Breakdown {
        method: &'static str,
        detail: String,
    },
    #[error("singular matrix at row {0}")]
    Singular(usize),
    #[error("dimension mismatch: {0}")]
    Dim(String),
}

pub trait LinearOperator<T: Scalar = f64>: Sync {
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    /// `y = A x`.
    fn apply(&self, x: &[T], y: &mut [T]);
    fn apply_transpose(&self, _x: &[T], _y: &mut [T]) -> Result<(), LinalgError> {
        Err(LinalgError::Unsupported("transpose".into()))
    }
}

pub trait Preconditioner<T: Scalar = f64>: Sync {
    /// `z = M^-1 r`.
    fn apply(&self, r: &[T], z: &mut [T]);
    fn apply_transpose(&self, r: &[T], z: &mut [T]) -> Result<(), LinalgError> {
        self.apply(r, z);
        Ok(())
    }
}

/// `M = I`: copies `r` into `z`.
#[derive(Copy, Clone, Debug, Default)]
pub struct IdentityPrecond;

impl<T: Scalar> Preconditioner<T> for IdentityPrecond {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
    }
}

/// Elmer/HUTI stopping semantics are carried by the solver's own config type;
/// this is only the common result.
#[derive(Clone, Debug, PartialEq)]
pub enum SolveStatus {
    Converged,
    MaxIterations,
    Diverged,
    Breakdown(String),
}

#[derive(Clone, Debug)]
pub struct SolveStats {
    pub iterations: usize,
    pub residual: f64,
    pub status: SolveStatus,
}

pub trait IterativeSolver<T: Scalar = f64> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<T>,
        m: &dyn Preconditioner<T>,
        b: &[T],
        x: &mut [T],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError>;
}

pub trait DirectSolver<T: Scalar = f64> {
    fn factor(&mut self, a: &CsrMatrix<T>) -> Result<(), LinalgError>;
    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError>;
}
