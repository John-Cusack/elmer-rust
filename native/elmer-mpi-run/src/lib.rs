//! Partitioned runs of the Rust-native Elmer: what `ElmerSolver_mpi` does
//! with a mesh that `ElmerGrid ... -metis N` (or `-partdual -metiskway N`)
//! split into `partitioning.N`.
//!
//! Every rank reads its own partition ([`ParMesh`]), assembles its own
//! elements with `elmer-solver`'s solvers, and the ranks meet in three
//! places, all provided here as `elmer-solver` backends ([`backend`]):
//!
//! * Dirichlet conditions are communicated to every holder of a shared dof
//!   and scaled by the summed diagonal (`EnforceDirichletConditions`);
//! * linear systems are solved by `elmer-dist`'s distributed HUTI Krylov
//!   (`ParallelIter`: split matrix, block-local preconditioner, parallel
//!   diagonal scaling);
//! * norms and reductions go through the [`elmer_core::Communicator`]
//!   (`elmer-mpi`), with `Nonlinear System Consistent Norm` counting each
//!   shared node once, on its owner.
//!
//! `Mesh Levels` splits every rank's piece and numbers the new nodes across
//! ranks ([`levels`]).
//!
//! [`run_partitioned`] composes a case; [`MpiEngine`] puts it behind
//! `elmer-cli`'s [`Engine`](elmer_cli::driver::Engine); the
//! `elmer-mpi-run` binary is the `mpirun` entry point. Every MPI process
//! runs with `HWLOC_COMPONENTS=-gl`.
#![forbid(unsafe_code)]

pub mod backend;
pub mod context;
pub mod levels;
pub mod mesh;
pub mod run;
pub mod transport;

pub use backend::{linear_registry, prepare_slot, ParAssembly};
pub use context::RankContext;
pub use mesh::ParMesh;
pub use run::{run_partitioned, MpiEngine, ParOutcome};
pub use transport::WorldExchange;
