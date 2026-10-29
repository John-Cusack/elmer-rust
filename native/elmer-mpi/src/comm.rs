//! rsmpi-backed [`Communicator`] and halo exchange.

use std::sync::Mutex;

use elmer_core::Communicator;
use mpi::collective::SystemOperation;
use mpi::environment::{Threading, Universe};
use mpi::request::WaitGuard;
use mpi::topology::Communicator as _;
use mpi::topology::SimpleCommunicator;
use mpi::traits::*;

use crate::halo::HaloPlan;

#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    #[error("MPI was already initialized in this process")]
    AlreadyInitialized,
    #[error("MPI provides thread level {0:?}, need at least Serialized")]
    ThreadLevel(Threading),
    #[error("halo buffer too short: {len} values for {need} (stride {stride})")]
    ShortBuffer { len: usize, need: usize, stride: usize },
}

/// Owns the MPI environment; `MPI_Finalize` runs when it is dropped, so keep it
/// alive for as long as any [`MpiComm`] made from it is used.
pub struct MpiRuntime {
    universe: Universe,
}

impl MpiRuntime {
    /// `MPI_Init_thread(Serialized)`. At most one runtime per process.
    pub fn init() -> Result<Self, MpiError> {
        let (universe, provided) = mpi::initialize_with_threading(Threading::Serialized)
            .ok_or(MpiError::AlreadyInitialized)?;
        if provided < Threading::Serialized {
            return Err(MpiError::ThreadLevel(provided));
        }
        Ok(Self { universe })
    }

    /// `MPI_COMM_WORLD`.
    pub fn world(&self) -> MpiComm {
        MpiComm::new(self.universe.world())
    }
}

struct Handle(SimpleCommunicator);

// SAFETY: the handle is an opaque communicator id. MPI is initialized with
// `Threading::Serialized` and every call goes through the `Mutex` in `MpiComm`,
// so calls are never concurrent, which is what that level requires.
unsafe impl Send for Handle {}

/// [`Communicator`] over an MPI communicator.
///
/// Like every MPI collective, calls must be made by all ranks in the same
/// order. Methods panic on an MPI error only through rsmpi/MPI's own abort
/// (the default `MPI_ERRORS_ARE_FATAL`).
pub struct MpiComm {
    inner: Mutex<Handle>,
    rank: usize,
    size: usize,
}

impl MpiComm {
    fn new(comm: SimpleCommunicator) -> Self {
        let rank = comm.rank() as usize;
        let size = comm.size() as usize;
        Self { inner: Mutex::new(Handle(comm)), rank, size }
    }

    fn with<R>(&self, f: impl FnOnce(&SimpleCommunicator) -> R) -> R {
        let guard = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        f(&guard.0)
    }

    /// Sum the contributions of every rank that shares a node. On return each
    /// sharer holds the total (assembly of interface rows). `values` holds
    /// `stride` dofs per local node.
    pub fn halo_sum(
        &self,
        plan: &HaloPlan,
        values: &mut [f64],
        stride: usize,
    ) -> Result<(), MpiError> {
        plan.check_len(values.len(), stride)?;
        let send: Vec<Vec<f64>> =
            plan.neighbours().iter().map(|n| pack(values, &n.shared, stride)).collect();
        let recv = self.exchange(plan, send, |n| n.shared.len() * stride);
        for (n, buf) in plan.neighbours().iter().zip(&recv) {
            for (k, &node) in n.shared.iter().enumerate() {
                for d in 0..stride {
                    values[node as usize * stride + d] += buf[k * stride + d];
                }
            }
        }
        Ok(())
    }

    /// Copy owner values over non-owner copies of shared nodes.
    pub fn halo_update(
        &self,
        plan: &HaloPlan,
        values: &mut [f64],
        stride: usize,
    ) -> Result<(), MpiError> {
        plan.check_len(values.len(), stride)?;
        let send: Vec<Vec<f64>> =
            plan.neighbours().iter().map(|n| pack(values, &n.send_owned, stride)).collect();
        let recv = self.exchange(plan, send, |n| n.recv_ghost.len() * stride);
        for (n, buf) in plan.neighbours().iter().zip(&recv) {
            for (k, &node) in n.recv_ghost.iter().enumerate() {
                for d in 0..stride {
                    values[node as usize * stride + d] = buf[k * stride + d];
                }
            }
        }
        Ok(())
    }

    /// One nonblocking send and receive per neighbour, completed together, so
    /// the pairwise ordering of neighbours cannot deadlock.
    fn exchange(
        &self,
        plan: &HaloPlan,
        send: Vec<Vec<f64>>,
        recv_len: impl Fn(&crate::halo::Neighbour) -> usize,
    ) -> Vec<Vec<f64>> {
        let mut recv: Vec<Vec<f64>> =
            plan.neighbours().iter().map(|n| vec![0.0; recv_len(n)]).collect();
        self.with(|comm| {
            mpi::request::scope(|scope| {
                let mut guards = Vec::with_capacity(2 * recv.len());
                for (n, buf) in plan.neighbours().iter().zip(recv.iter_mut()) {
                    guards.push(WaitGuard::from(
                        comm.process_at_rank(n.rank as i32)
                            .immediate_receive_into(scope, &mut buf[..]),
                    ));
                }
                for (n, buf) in plan.neighbours().iter().zip(&send) {
                    guards.push(WaitGuard::from(
                        comm.process_at_rank(n.rank as i32).immediate_send(scope, &buf[..]),
                    ));
                }
                // Guards wait on drop, at the end of the scope closure.
            });
        });
        recv
    }
}

fn pack(values: &[f64], nodes: &[u32], stride: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(nodes.len() * stride);
    for &n in nodes {
        out.extend_from_slice(&values[n as usize * stride..(n as usize + 1) * stride]);
    }
    out
}

impl Communicator for MpiComm {
    fn rank(&self) -> usize {
        self.rank
    }

    fn size(&self) -> usize {
        self.size
    }

    fn sum_f64(&self, buf: &mut [f64]) {
        let local = buf.to_vec();
        self.with(|c| c.all_reduce_into(&local[..], &mut buf[..], SystemOperation::sum()));
    }

    fn max_f64(&self, buf: &mut [f64]) {
        let local = buf.to_vec();
        self.with(|c| c.all_reduce_into(&local[..], &mut buf[..], SystemOperation::max()));
    }

    fn barrier(&self) {
        self.with(|c| c.barrier());
    }
}

/// Global dot product of two distributed vectors: each shared node counts once,
/// on the rank that owns it (`owned[i]`), then the partial sums are added up.
pub fn dot_owned(comm: &dyn Communicator, x: &[f64], y: &[f64], owned: &[bool]) -> f64 {
    let mut s =
        [x.iter().zip(y).zip(owned).filter(|(_, &o)| o).map(|((a, b), _)| a * b).sum::<f64>()];
    comm.sum_f64(&mut s);
    s[0]
}
