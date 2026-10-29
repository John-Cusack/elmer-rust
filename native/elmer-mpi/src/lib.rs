//! MPI support for the Rust-native Elmer rewrite.
//!
//! * [`comm`]: [`MpiRuntime`] / [`MpiComm`], an rsmpi (C MPI, no Fortran)
//!   implementation of [`elmer_core::Communicator`], plus the point-to-point
//!   halo exchange on shared nodes.
//! * [`partition`]: reader for the `partitioning.N/part.K.*` directories that
//!   `ElmerGrid ... -metis N` writes (owned/ghost nodes and elements, global
//!   ids, shared-node lists), producing a [`PartitionMesh`] that implements
//!   [`elmer_core::MeshView`].
//! * [`halo`]: the [`HaloPlan`] derived from a [`PartitionMesh`].
//!
//! Every MPI process of Elmer or of this crate's tests must run with
//! `HWLOC_COMPONENTS=-gl`, and the swarm limit is `np <= 2`.
#![deny(unsafe_op_in_unsafe_fn)]

pub mod comm;
pub mod halo;
pub mod partition;

pub use comm::{dot_owned, MpiComm, MpiError, MpiRuntime};
pub use halo::{HaloPlan, Neighbour};
pub use partition::{
    partition_subdir, PartElement, PartError, PartHeader, PartitionMesh, SharedNode,
};
