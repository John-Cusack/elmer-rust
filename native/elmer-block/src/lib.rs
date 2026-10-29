//! Elmer's block linear system (`Linear System Solver = Block`,
//! `BlockSolve.F90`) for the Rust-native rewrite.
//!
//! [`block_solve`] is `BlockSolveExt`/`BlockSolveInt` (BlockSolve.F90:5007,
//! :5397) for the serial default route: the monolithic system is split by
//! component or by `Block Structure` ([`BlockSystem`], `BlockPickMatrix`);
//! the outer Krylov iteration (`BlockKrylovIter`) runs under the `outer:`
//! namespace over the block matrix-vector product ([`BlockOperator`],
//! `BlockMatrixVectorProd`), preconditioned by block Jacobi or Gauss-Seidel
//! sweeps ([`BlockPreconditioner`], `BlockMatrixPrec`) whose diagonal blocks
//! are solved by an [`InnerSolver`] under the `block:` and `block NN:`
//! namespaces. [`ElmerInner`] is that inner solve: `SolveLinearSystem`'s
//! scaling and zero right-hand side shortcut around the `elmer-krylov`,
//! `elmer-itermethods` and `elmer-direct` methods with the `elmer-sparse`
//! preconditioners, as IterSolve and DirectSolve select them.
//!
//! The caller (`SolveLinearSystem` in `elmer-solver`) dispatches here after
//! its own `Linear System Scaling`, passing that scaling
//! ([`BlockOptions::outer_scaling`]) and, for `Block Matrix Schur Variable`,
//! the Schur solver's matrix ([`BlockOptions::schur`]).
//!
//! Options not ported are refused with [`BlockError::Unsupported`], never
//! approximated: the other block splittings (A-V, Hcurl, Hdiv, Re-Im, nodal,
//! horizontal-vertical, Cartesian, domain, nested, complex), `Block Solvers`,
//! `Block Monolithic`, `Block Scaling`, `Block Preconditioner = False`
//! (`BlockStandardIter`), the `Prec Matrix ...` preconditioning matrices,
//! Schur operators and approximations, complex systems, constraint and add
//! matrices, and parallel runs.
#![forbid(unsafe_code)]

mod inner;
pub mod namespace;
mod operator;
mod precond;
mod solve;
pub mod split;

pub use inner::{ElmerInner, InnerSolver, InnerStats};
pub use operator::BlockOperator;
pub use precond::{BlockPreconditioner, InnerRecord, PrecBlock, PrecSettings};
pub use solve::{block_solve, block_solve_with, comp_norm, BlockOptions, BlockReport, OuterScaling};
pub use split::{BlockLayout, BlockSystem};

use elmer_core::LinalgError;

/// Errors of a block solve.
#[derive(Debug, thiserror::Error)]
pub enum BlockError {
    /// Elmer `CALL Fatal(caller, msg)`.
    #[error("{caller}: {msg}")]
    Fatal { caller: String, msg: String },
    /// Elmer `CALL NumericalError(caller, msg)` with `Global Abort Not
    /// Converged` in force (an iteration hit its maximum count).
    #[error("{caller}: Numerical Error: {msg}")]
    Numerical { caller: String, msg: String },
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error(transparent)]
    Linalg(#[from] LinalgError),
}

impl BlockError {
    pub fn fatal(caller: &str, msg: impl Into<String>) -> Self {
        BlockError::Fatal {
            caller: caller.into(),
            msg: msg.into(),
        }
    }
}
