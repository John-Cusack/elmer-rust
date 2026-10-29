//! [`WorldExchange`]: the MPI transport a partitioned run owns.

use elmer_core::Communicator;
use elmer_dist::{Exchange, MpiExchange};
use elmer_mpi::MpiComm;

/// [`Exchange`] over `MPI_COMM_WORLD` that owns its [`MpiComm`], so the
/// run's [`crate::RankContext`] can hold it. Drop it before the
/// [`elmer_mpi::MpiRuntime`].
pub struct WorldExchange(pub MpiComm);

impl Communicator for WorldExchange {
    fn rank(&self) -> usize {
        self.0.rank()
    }

    fn size(&self) -> usize {
        self.0.size()
    }

    fn sum_f64(&self, buf: &mut [f64]) {
        self.0.sum_f64(buf);
    }

    fn max_f64(&self, buf: &mut [f64]) {
        self.0.max_f64(buf);
    }

    fn barrier(&self) {
        self.0.barrier();
    }
}

impl Exchange for WorldExchange {
    fn all_to_all_u64(&self, send: &[Vec<u64>]) -> Vec<Vec<u64>> {
        MpiExchange::new(&self.0).all_to_all_u64(send)
    }

    fn exchange_f64(&self, send: &[(usize, Vec<f64>)], recv: &mut [(usize, Vec<f64>)]) {
        MpiExchange::new(&self.0).exchange_f64(send, recv);
    }
}
