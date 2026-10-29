//! Process communicator: distributed-ready (strategy report 3.5), serial now.

pub trait Communicator: Send + Sync {
    fn rank(&self) -> usize;
    fn size(&self) -> usize;
    /// In-place element-wise sum across ranks.
    fn sum_f64(&self, buf: &mut [f64]);
    /// In-place element-wise max across ranks.
    fn max_f64(&self, buf: &mut [f64]);
    fn barrier(&self);
}

/// One rank: reductions are the identity.
#[derive(Copy, Clone, Debug, Default)]
pub struct SerialComm;

impl Communicator for SerialComm {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        1
    }
    fn sum_f64(&self, _buf: &mut [f64]) {}
    fn max_f64(&self, _buf: &mut [f64]) {}
    fn barrier(&self) {}
}
