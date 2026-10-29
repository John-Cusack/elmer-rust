//! Halo-exchange plan: which local nodes each neighbouring rank shares.

use crate::comm::MpiError;
use elmer_core::MeshView;

use crate::partition::PartitionMesh;

/// Nodes shared with one neighbouring rank. All three lists are ordered by
/// global node id, so both sides of a pair agree on the message layout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Neighbour {
    pub rank: usize,
    /// Every local node shared with `rank`.
    pub shared: Vec<u32>,
    /// Shared nodes this rank owns (owner values are sent to `rank`).
    pub send_owned: Vec<u32>,
    /// Shared nodes `rank` owns (owner values arrive from `rank`).
    pub recv_ghost: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HaloPlan {
    n_nodes: usize,
    neighbours: Vec<Neighbour>,
}

impl HaloPlan {
    pub fn from_mesh(mesh: &PartitionMesh) -> Self {
        let me = mesh.rank();
        let mut by_rank: std::collections::BTreeMap<usize, Vec<(u64, u32, usize)>> =
            Default::default();
        for s in mesh.shared() {
            for &r in s.parts.iter().filter(|&&r| r != me) {
                by_rank.entry(r).or_default().push((s.gid, s.local, s.parts[0]));
            }
        }
        let neighbours = by_rank
            .into_iter()
            .map(|(rank, mut v)| {
                v.sort_unstable();
                Neighbour {
                    rank,
                    shared: v.iter().map(|t| t.1).collect(),
                    send_owned: v.iter().filter(|t| t.2 == me).map(|t| t.1).collect(),
                    recv_ghost: v.iter().filter(|t| t.2 == rank).map(|t| t.1).collect(),
                }
            })
            .collect();
        Self { n_nodes: mesh.n_nodes(), neighbours }
    }

    pub fn neighbours(&self) -> &[Neighbour] {
        &self.neighbours
    }

    pub(crate) fn check_len(&self, len: usize, stride: usize) -> Result<(), MpiError> {
        let need = self.n_nodes * stride;
        if len < need {
            return Err(MpiError::ShortBuffer { len, need, stride });
        }
        Ok(())
    }
}
