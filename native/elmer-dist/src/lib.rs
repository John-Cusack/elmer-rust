//! Distributed linear algebra for the Rust-native Elmer: what
//! SParIterSolver.F90, SParIterComm.F90 and ParallelUtils.F90 do around one
//! parallel HUTI solve.
//!
//! * [`exchange`]: the [`Exchange`] transport (reductions from
//!   [`elmer_core::Communicator`] plus point-to-point traffic), for one rank
//!   ([`elmer_core::SerialComm`]), in-process ranks ([`ThreadComm`]) and MPI
//!   ([`MpiExchange`] over [`elmer_mpi::MpiComm`]).
//! * [`layout`]: [`DofLayout`], which rank owns each local dof and who holds
//!   copies of it (`ParallelInfo % NeighbourList`, `GlobalDOFs`).
//! * [`matrix`]: [`DistCsr`], Elmer's split matrix: owned rows by owned
//!   columns (`InsideMatrix`) plus interface blocks kept by the column owner
//!   (`IfMatrix`), and the [`DistOperator`] matvec (`SParMatrixVector`).
//! * [`solve`]: [`dist_solve`], diagonal scaling, the block-local
//!   preconditioner and the HUTI solve on owned dofs with global dot products.
//! * [`linsys`]: reader for the per-rank `Linear System Save` files a parallel
//!   ElmerSolver run writes (`linsys_a.dat.R`, `linsys_b.dat.R`).
//!
//! Every MPI process must run with `HWLOC_COMPONENTS=-gl`; the swarm limit is
//! `np <= 2`.
#![forbid(unsafe_code)]

pub mod exchange;
pub mod layout;
pub mod linsys;
pub mod matrix;
pub mod solve;

pub use exchange::{Exchange, MpiExchange, ThreadComm};
pub use layout::{owners_by_lowest_rank, DofLayout};
pub use linsys::{read_rank_system, RankSystem};
pub use matrix::{DistCsr, DistOperator};
pub use solve::{dist_solve, DistPrecond, DistSolveOptions};

/// Errors from distributed setup and solves. Setup is collective: when one
/// rank fails, every rank returns an error rather than waiting forever.
#[derive(Debug, thiserror::Error)]
pub enum DistError {
    #[error("rank {rank}: {msg}")]
    Local { rank: usize, msg: String },
    #[error("rank {rank}: another rank failed")]
    Remote { rank: usize },
    #[error(transparent)]
    Linalg(#[from] elmer_core::LinalgError),
    #[error("{path}: {msg}")]
    File {
        path: std::path::PathBuf,
        msg: String,
    },
}

/// Collective agreement on a fallible local step: every rank gets an error
/// when any rank's `local` is an error.
pub(crate) fn agree<T>(
    comm: &dyn elmer_core::Communicator,
    local: Result<T, String>,
) -> Result<T, DistError> {
    let mut flag = [if local.is_err() { 1.0 } else { 0.0 }];
    comm.max_f64(&mut flag);
    let rank = comm.rank();
    match local {
        Err(msg) => Err(DistError::Local { rank, msg }),
        Ok(_) if flag[0] != 0.0 => Err(DistError::Remote { rank }),
        Ok(v) => Ok(v),
    }
}
