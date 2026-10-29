//! Which rank owns each local dof and which ranks hold copies of it.

use std::collections::HashMap;

use elmer_core::NO_DOF;
use elmer_mpi::PartitionMesh;

use crate::exchange::Exchange;
use crate::{agree, DistError};

/// One rank's dofs: global id and owning rank of every local dof, the owned
/// dofs in local order (the order `SolveHutiter` packs `TmpXVec` in), and the
/// copy lists that move values between owners and the ranks holding ghosts.
///
/// Built collectively by [`DofLayout::new`]. Copy lists are ordered by global
/// id, so both sides of a pair agree on message layout.
#[derive(Clone, Debug)]
pub struct DofLayout {
    rank: usize,
    size: usize,
    gid: Vec<u64>,
    owner: Vec<usize>,
    /// Local index of each owned dof.
    owned: Vec<u32>,
    /// Local dof -> position in `owned`, [`NO_DOF`] for ghosts.
    owned_pos: Vec<u32>,
    /// `(rank, owned positions)`: owned dofs that `rank` holds a copy of.
    copies_out: Vec<(usize, Vec<u32>)>,
    /// `(owner, local dofs)`: ghost dofs held here, by owner.
    copies_in: Vec<(usize, Vec<u32>)>,
}

impl DofLayout {
    /// Collective. `gid[i]` and `owner[i]` describe local dof `i`; every rank
    /// holding a dof must name the same owner, and the owner must hold it.
    pub fn new(comm: &dyn Exchange, gid: Vec<u64>, owner: Vec<usize>) -> Result<Self, DistError> {
        let (rank, size) = (comm.rank(), comm.size());
        let by_gid = agree(comm, check_local(rank, size, &gid, &owner))?;
        let mut owned_pos = vec![NO_DOF; gid.len()];
        let mut owned = Vec::new();
        for (i, &o) in owner.iter().enumerate() {
            if o == rank {
                owned_pos[i] = owned.len() as u32;
                owned.push(i as u32);
            }
        }

        // Ghost holders tell each owner which of its dofs they copy.
        let mut want: Vec<Vec<(u64, u32)>> = vec![Vec::new(); size];
        for (i, &o) in owner.iter().enumerate() {
            if o != rank {
                want[o].push((gid[i], i as u32));
            }
        }
        for w in want.iter_mut() {
            w.sort_unstable();
        }
        let requests: Vec<Vec<u64>> = want
            .iter()
            .map(|w| w.iter().map(|&(g, _)| g).collect())
            .collect();
        let asked = comm.all_to_all_u64(&requests);

        let mut copies_out = Vec::new();
        let mut bad = None;
        for (r, gids) in asked.iter().enumerate() {
            if gids.is_empty() {
                continue;
            }
            let mut pos = Vec::with_capacity(gids.len());
            for g in gids {
                match by_gid
                    .get(g)
                    .copied()
                    .filter(|&l| owner[l as usize] == rank)
                {
                    Some(l) => pos.push(owned_pos[l as usize]),
                    None => {
                        bad.get_or_insert(format!(
                            "rank {r} names this rank owner of dof {g}, which it does not own"
                        ));
                    }
                }
            }
            copies_out.push((r, pos));
        }
        agree(comm, bad.map_or(Ok(()), Err))?;

        let copies_in = want
            .into_iter()
            .enumerate()
            .filter(|(_, w)| !w.is_empty())
            .map(|(o, w)| (o, w.into_iter().map(|(_, l)| l).collect()))
            .collect();
        Ok(DofLayout {
            rank,
            size,
            gid,
            owner,
            owned,
            owned_pos,
            copies_out,
            copies_in,
        })
    }

    /// Nodal dofs of a partitioned mesh, `dofs` per node, node-major. Global
    /// dof id `dofs * (node_gid - 1) + k` (Elmer's `GlobalDOFs`, 0-based), the
    /// node's owner owns all its dofs. Local dof `dofs * node + k`.
    pub fn nodal(
        comm: &dyn Exchange,
        mesh: &PartitionMesh,
        dofs: usize,
    ) -> Result<Self, DistError> {
        let d = dofs as u64;
        let gid = mesh
            .node_gids()
            .iter()
            .flat_map(|&g| (0..d).map(move |k| d * (g - 1) + k))
            .collect();
        let owner = mesh
            .node_owners()
            .iter()
            .flat_map(|&o| std::iter::repeat_n(o, dofs))
            .collect();
        Self::new(comm, gid, owner)
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn n_local(&self) -> usize {
        self.gid.len()
    }

    pub fn n_owned(&self) -> usize {
        self.owned.len()
    }

    pub fn gids(&self) -> &[u64] {
        &self.gid
    }

    pub fn owners(&self) -> &[usize] {
        &self.owner
    }

    /// Local index of each owned dof, in local order.
    pub fn owned(&self) -> &[u32] {
        &self.owned
    }

    /// Local dof -> position among owned dofs, [`NO_DOF`] for ghosts.
    pub fn owned_pos(&self) -> &[u32] {
        &self.owned_pos
    }

    /// The owned entries of a local vector.
    pub fn gather_owned(&self, local: &[f64]) -> Vec<f64> {
        self.owned.iter().map(|&i| local[i as usize]).collect()
    }

    /// Owned totals of a local vector that holds partial sums at shared dofs:
    /// each owner adds the ghost copies' values to its own
    /// (`ExchangeRHSIf`). Collective.
    pub fn sum_to_owned(&self, comm: &dyn Exchange, local: &[f64]) -> Vec<f64> {
        let mut owned = self.gather_owned(local);
        let send: Vec<(usize, Vec<f64>)> = self
            .copies_in
            .iter()
            .map(|(o, l)| (*o, l.iter().map(|&i| local[i as usize]).collect()))
            .collect();
        let mut recv: Vec<(usize, Vec<f64>)> = self
            .copies_out
            .iter()
            .map(|(r, p)| (*r, vec![0.0; p.len()]))
            .collect();
        comm.exchange_f64(&send, &mut recv);
        for ((_, pos), (_, vals)) in self.copies_out.iter().zip(&recv) {
            for (&p, &v) in pos.iter().zip(vals) {
                owned[p as usize] += v;
            }
        }
        owned
    }

    /// Write owned values into a local vector and overwrite every ghost copy
    /// with its owner's value (`ExchangeResult`). Collective.
    pub fn owned_to_local(&self, comm: &dyn Exchange, owned: &[f64], local: &mut [f64]) {
        for (&i, &v) in self.owned.iter().zip(owned) {
            local[i as usize] = v;
        }
        let send: Vec<(usize, Vec<f64>)> = self
            .copies_out
            .iter()
            .map(|(r, p)| (*r, p.iter().map(|&k| owned[k as usize]).collect()))
            .collect();
        let mut recv: Vec<(usize, Vec<f64>)> = self
            .copies_in
            .iter()
            .map(|(o, l)| (*o, vec![0.0; l.len()]))
            .collect();
        comm.exchange_f64(&send, &mut recv);
        for ((_, dofs), (_, vals)) in self.copies_in.iter().zip(&recv) {
            for (&i, &v) in dofs.iter().zip(vals) {
                local[i as usize] = v;
            }
        }
    }

    /// Every copy of a shared dof gets the sum over its holders
    /// (`ParallelSumVector`). Collective.
    pub fn sum_shared(&self, comm: &dyn Exchange, local: &mut [f64]) {
        let owned = self.sum_to_owned(comm, local);
        self.owned_to_local(comm, &owned, local);
    }
}

/// Validate one rank's input; returns the gid -> local index map.
fn check_local(
    rank: usize,
    size: usize,
    gid: &[u64],
    owner: &[usize],
) -> Result<HashMap<u64, u32>, String> {
    if gid.len() != owner.len() {
        return Err(format!(
            "{} global ids but {} owners",
            gid.len(),
            owner.len()
        ));
    }
    if gid.len() >= NO_DOF as usize {
        return Err(format!("{} local dofs do not fit u32 indices", gid.len()));
    }
    if let Some(o) = owner.iter().find(|&&o| o >= size) {
        return Err(format!("owner rank {o} out of range for {size} ranks"));
    }
    let mut by_gid = HashMap::with_capacity(gid.len());
    for (i, &g) in gid.iter().enumerate() {
        if by_gid.insert(g, i as u32).is_some() {
            return Err(format!("global dof {g} appears twice on rank {rank}"));
        }
    }
    Ok(by_gid)
}

/// Owner of each local dof when the holders are known only by global id:
/// the lowest rank holding the dof (Elmer's own rule picks some holder; any
/// consistent choice gives the same assembled operator). Collective.
pub fn owners_by_lowest_rank(comm: &dyn Exchange, gid: &[u64]) -> Vec<usize> {
    let size = comm.size();
    // Each dof's "home" rank collects its holders and answers.
    let home = |g: u64| (g % size as u64) as usize;
    let mut ask: Vec<Vec<u64>> = vec![Vec::new(); size];
    for &g in gid {
        ask[home(g)].push(g);
    }
    let asked = comm.all_to_all_u64(&ask);
    let mut lowest: HashMap<u64, usize> = HashMap::new();
    for (r, gids) in asked.iter().enumerate() {
        for &g in gids {
            lowest
                .entry(g)
                .and_modify(|o| *o = (*o).min(r))
                .or_insert(r);
        }
    }
    let answers: Vec<Vec<u64>> = asked
        .iter()
        .map(|gids| gids.iter().map(|g| lowest[g] as u64).collect())
        .collect();
    let got = comm.all_to_all_u64(&answers);
    let mut next = vec![0usize; size];
    gid.iter()
        .map(|&g| {
            let h = home(g);
            let o = got[h][next[h]] as usize;
            next[h] += 1;
            o
        })
        .collect()
}
