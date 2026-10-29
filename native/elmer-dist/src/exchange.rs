//! Point-to-point transport on top of [`Communicator`]'s reductions.

use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Barrier, Mutex};

use elmer_core::{Communicator, SerialComm};
use elmer_mpi::MpiComm;
use mpi::datatype::{Partition, PartitionMut};
use mpi::topology::Communicator as _;
use mpi::topology::SimpleCommunicator;
use mpi::traits::*;

/// MPI tag of [`Exchange::exchange_f64`] traffic (SParMatrixVector's 6001).
const TAG: i32 = 6001;

/// The traffic a distributed operator needs beyond sums and maxima. Every
/// method is collective: all ranks call it in the same order.
pub trait Exchange: Communicator {
    /// Personalised all-to-all: `send[r]` goes to rank `r` (`send.len()` is the
    /// communicator size). Returns what each rank sent here, by source rank.
    fn all_to_all_u64(&self, send: &[Vec<u64>]) -> Vec<Vec<u64>>;

    /// Exchange over a pattern both sides agreed on beforehand: `send` lists
    /// `(destination, values)`, `recv` lists `(source, buffer)` with each buffer
    /// already sized to what its source sends. A rank never names itself.
    fn exchange_f64(&self, send: &[(usize, Vec<f64>)], recv: &mut [(usize, Vec<f64>)]);
}

impl Exchange for SerialComm {
    fn all_to_all_u64(&self, send: &[Vec<u64>]) -> Vec<Vec<u64>> {
        send.to_vec()
    }

    fn exchange_f64(&self, send: &[(usize, Vec<f64>)], recv: &mut [(usize, Vec<f64>)]) {
        assert!(
            send.is_empty() && recv.is_empty(),
            "a single rank has no neighbours"
        );
    }
}

enum Msg {
    U64(Vec<u64>),
    F64(Vec<f64>),
}

/// In-process ranks on threads, for running a distributed algorithm at np > 1
/// without MPI. Messages travel on one FIFO channel per ordered rank pair,
/// which is enough because every operation is collective. Sums and maxima
/// combine contributions in rank order, so every rank gets the same bits.
pub struct ThreadComm {
    rank: usize,
    size: usize,
    to: Vec<Sender<Msg>>,
    from: Vec<Mutex<Receiver<Msg>>>,
    barrier: std::sync::Arc<Barrier>,
}

impl ThreadComm {
    /// Run `f` on `size` ranks, one thread each, and return each rank's result
    /// in rank order.
    pub fn run<R: Send>(size: usize, f: impl Fn(&ThreadComm) -> R + Sync) -> Vec<R> {
        assert!(size > 0);
        // senders[src][dst], receivers[dst][src].
        let mut senders: Vec<Vec<Sender<Msg>>> = (0..size).map(|_| Vec::new()).collect();
        let mut receivers: Vec<Vec<Mutex<Receiver<Msg>>>> = (0..size).map(|_| Vec::new()).collect();
        for s in senders.iter_mut() {
            for r in receivers.iter_mut() {
                let (tx, rx) = channel();
                s.push(tx);
                r.push(Mutex::new(rx));
            }
        }
        let barrier = std::sync::Arc::new(Barrier::new(size));
        let comms: Vec<ThreadComm> = senders
            .into_iter()
            .zip(receivers)
            .enumerate()
            .map(|(rank, (to, from))| ThreadComm {
                rank,
                size,
                to,
                from,
                barrier: barrier.clone(),
            })
            .collect();
        let f = &f;
        std::thread::scope(|s| {
            let handles: Vec<_> = comms.iter().map(|c| s.spawn(move || f(c))).collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("rank thread panicked"))
                .collect()
        })
    }

    fn send(&self, dst: usize, m: Msg) {
        self.to[dst].send(m).expect("peer rank exited");
    }

    fn recv(&self, src: usize) -> Msg {
        let rx = self.from[src].lock().unwrap_or_else(|e| e.into_inner());
        rx.recv().expect("peer rank exited")
    }

    fn recv_f64(&self, src: usize) -> Vec<f64> {
        match self.recv(src) {
            Msg::F64(v) => v,
            Msg::U64(_) => panic!("rank {}: expected f64 message from {src}", self.rank),
        }
    }

    /// Every rank's copy of `buf`, by rank.
    fn all_gather(&self, buf: &[f64]) -> Vec<Vec<f64>> {
        for r in (0..self.size).filter(|&r| r != self.rank) {
            self.send(r, Msg::F64(buf.to_vec()));
        }
        (0..self.size)
            .map(|r| {
                if r == self.rank {
                    buf.to_vec()
                } else {
                    self.recv_f64(r)
                }
            })
            .collect()
    }
}

impl Communicator for ThreadComm {
    fn rank(&self) -> usize {
        self.rank
    }

    fn size(&self) -> usize {
        self.size
    }

    fn sum_f64(&self, buf: &mut [f64]) {
        let all = self.all_gather(buf);
        for (i, b) in buf.iter_mut().enumerate() {
            *b = all.iter().fold(0.0, |s, v| s + v[i]);
        }
    }

    fn max_f64(&self, buf: &mut [f64]) {
        let all = self.all_gather(buf);
        for (i, b) in buf.iter_mut().enumerate() {
            *b = all.iter().map(|v| v[i]).fold(f64::NEG_INFINITY, f64::max);
        }
    }

    fn barrier(&self) {
        self.barrier.wait();
    }
}

impl Exchange for ThreadComm {
    fn all_to_all_u64(&self, send: &[Vec<u64>]) -> Vec<Vec<u64>> {
        assert_eq!(send.len(), self.size);
        for (r, v) in send.iter().enumerate().filter(|&(r, _)| r != self.rank) {
            self.send(r, Msg::U64(v.clone()));
        }
        (0..self.size)
            .map(|r| {
                if r == self.rank {
                    return send[r].clone();
                }
                match self.recv(r) {
                    Msg::U64(v) => v,
                    Msg::F64(_) => panic!("rank {}: expected u64 message from {r}", self.rank),
                }
            })
            .collect()
    }

    fn exchange_f64(&self, send: &[(usize, Vec<f64>)], recv: &mut [(usize, Vec<f64>)]) {
        for (dst, v) in send {
            self.send(*dst, Msg::F64(v.clone()));
        }
        for (src, buf) in recv.iter_mut() {
            let v = self.recv_f64(*src);
            assert_eq!(
                v.len(),
                buf.len(),
                "rank {}: message size from {src}",
                self.rank
            );
            buf.copy_from_slice(&v);
        }
    }
}

/// [`Exchange`] over an [`MpiComm`] on `MPI_COMM_WORLD` (what
/// [`elmer_mpi::MpiRuntime::world`] returns). Reductions go through the
/// `MpiComm`; point-to-point traffic uses the world communicator directly.
/// Like `MpiComm`, it assumes calls are never made from two threads at once.
pub struct MpiExchange<'a> {
    comm: &'a MpiComm,
}

impl<'a> MpiExchange<'a> {
    pub fn new(comm: &'a MpiComm) -> Self {
        MpiExchange { comm }
    }
}

impl Communicator for MpiExchange<'_> {
    fn rank(&self) -> usize {
        self.comm.rank()
    }

    fn size(&self) -> usize {
        self.comm.size()
    }

    fn sum_f64(&self, buf: &mut [f64]) {
        self.comm.sum_f64(buf);
    }

    fn max_f64(&self, buf: &mut [f64]) {
        self.comm.max_f64(buf);
    }

    fn barrier(&self) {
        self.comm.barrier();
    }
}

fn displs(counts: &[i32]) -> Vec<i32> {
    counts
        .iter()
        .scan(0, |acc, &c| {
            let d = *acc;
            *acc += c;
            Some(d)
        })
        .collect()
}

impl Exchange for MpiExchange<'_> {
    fn all_to_all_u64(&self, send: &[Vec<u64>]) -> Vec<Vec<u64>> {
        let world = SimpleCommunicator::world();
        let size = self.size();
        assert_eq!(send.len(), size);
        let send_counts: Vec<i32> = send.iter().map(|v| v.len() as i32).collect();
        let mut recv_counts = vec![0i32; size];
        world.all_to_all_into(&send_counts[..], &mut recv_counts[..]);

        let flat: Vec<u64> = send.concat();
        let sd = displs(&send_counts);
        let rd = displs(&recv_counts);
        let mut out = vec![0u64; recv_counts.iter().sum::<i32>() as usize];
        {
            let sp = Partition::new(&flat[..], &send_counts[..], &sd[..]);
            let mut rp = PartitionMut::new(&mut out[..], &recv_counts[..], &rd[..]);
            world.all_to_all_varcount_into(&sp, &mut rp);
        }
        recv_counts
            .iter()
            .zip(&rd)
            .map(|(&c, &d)| out[d as usize..(d + c) as usize].to_vec())
            .collect()
    }

    fn exchange_f64(&self, send: &[(usize, Vec<f64>)], recv: &mut [(usize, Vec<f64>)]) {
        let world = SimpleCommunicator::world();
        mpi::request::scope(|scope| {
            let mut guards = Vec::with_capacity(send.len() + recv.len());
            for (src, buf) in recv.iter_mut() {
                guards.push(mpi::request::WaitGuard::from(
                    world
                        .process_at_rank(*src as i32)
                        .immediate_receive_into_with_tag(scope, &mut buf[..], TAG),
                ));
            }
            for (dst, buf) in send {
                guards.push(mpi::request::WaitGuard::from(
                    world.process_at_rank(*dst as i32).immediate_send_with_tag(
                        scope,
                        &buf[..],
                        TAG,
                    ),
                ));
            }
            // Guards wait on drop, at the end of the scope closure.
        });
    }
}
