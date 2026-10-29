//! Elmer's split parallel matrix (`SplitMatrix`, SParIterSolver.F90:143) and
//! its matvec (`SParMatrixVectorVals`, SParIterSolver.F90:2763).
//!
//! Each rank starts from its local assembled matrix: rows and columns are its
//! local dofs, and rows of shared dofs hold only this rank's contributions.
//! The global operator is the sum of the local ones. Every entry `(i, j)` moves
//! to the rank that owns column `j`, because that rank holds `x_j`:
//!
//! * row and column owned here: the `InsideMatrix` (owned positions on both
//!   axes), with contributions from all ranks added up;
//! * row owned by rank `r`, column owned here: the interface block for `r`.
//!
//! The matvec is then `y = Inside x` plus, from each neighbour, the partial row
//! sums its interface block computed against its own `x`. No ghost values of
//! `x` are needed and the Krylov vectors hold owned dofs only.

use std::collections::HashMap;

use elmer_core::{CsrMatrix, LinearOperator, NO_DOF};
use elmer_sparse::from_triplets;

use crate::exchange::Exchange;
use crate::layout::DofLayout;
use crate::{agree, DistError};

/// Interface rows this rank evaluates for the row owner `rank`: CSR over
/// owned columns, rows in the order the owner expects them.
#[derive(Clone, Debug)]
struct IfBlock {
    rank: usize,
    row_ptr: Vec<usize>,
    cols: Vec<u32>,
    vals: Vec<f64>,
}

/// One rank's part of a distributed matrix. See the module docs.
#[derive(Clone, Debug)]
pub struct DistCsr {
    inside: CsrMatrix<f64>,
    send_if: Vec<IfBlock>,
    /// `(rank, owned positions)`: where the partial sums from `rank` land.
    recv_if: Vec<(usize, Vec<u32>)>,
}

impl DistCsr {
    /// Collective. `local` is this rank's assembled matrix over the dofs of
    /// `layout` (square, `layout.n_local()` rows).
    pub fn from_local(
        comm: &dyn Exchange,
        layout: &DofLayout,
        local: &CsrMatrix<f64>,
    ) -> Result<Self, DistError> {
        let (me, size) = (comm.rank(), comm.size());
        let n = layout.n_local();
        agree(
            comm,
            if local.nrows == n && local.ncols == n && local.row_ptr.len() == n + 1 {
                Ok(())
            } else {
                Err(format!(
                    "local matrix {}x{} for {n} local dofs",
                    local.nrows, local.ncols
                ))
            },
        )?;
        let (gid, owner, pos) = (layout.gids(), layout.owners(), layout.owned_pos());

        let mut inside: Vec<(usize, u32, f64)> = Vec::new();
        // Interface entries kept here: (row owner, row gid, owned column, value).
        let mut iface: Vec<(usize, u64, u32, f64)> = Vec::new();
        // To column owners: row gid, row owner, column gid, value bits.
        let mut out: Vec<Vec<u64>> = vec![Vec::new(); size];
        for i in 0..n {
            for k in local.row_ptr[i]..local.row_ptr[i + 1] {
                let j = local.col_idx[k] as usize;
                let v = local.values[k];
                match (owner[i] == me, owner[j] == me) {
                    (true, true) => inside.push((pos[i] as usize, pos[j], v)),
                    (false, true) => iface.push((owner[i], gid[i], pos[j], v)),
                    (_, false) => {
                        out[owner[j]].extend([gid[i], owner[i] as u64, gid[j], v.to_bits()])
                    }
                }
            }
        }
        let got = comm.all_to_all_u64(&out);

        let owned_by_gid: HashMap<u64, u32> = layout
            .owned()
            .iter()
            .map(|&l| (gid[l as usize], pos[l as usize]))
            .collect();
        let mut bad = None;
        for (src, buf) in got.iter().enumerate() {
            for e in buf.chunks_exact(4) {
                let (gi, ro, gj, v) = (e[0], e[1] as usize, e[2], f64::from_bits(e[3]));
                let Some(&pj) = owned_by_gid.get(&gj) else {
                    bad.get_or_insert(format!("rank {src} sent column {gj}, not owned here"));
                    continue;
                };
                if ro == me {
                    match owned_by_gid.get(&gi) {
                        Some(&pi) => inside.push((pi as usize, pj, v)),
                        None => {
                            bad.get_or_insert(format!("rank {src} sent row {gi}, not owned here"));
                        }
                    }
                } else if ro < size {
                    iface.push((ro, gi, pj, v));
                } else {
                    bad.get_or_insert(format!("rank {src} sent row owner {ro}"));
                }
            }
        }
        agree(comm, bad.map_or(Ok(()), Err))?;

        let n_own = layout.n_owned();
        // Duplicates add up in arrival order: local entries, then by source
        // rank. Explicit zero diagonal slots for the preconditioners.
        inside.extend((0..n_own).map(|i| (i, i as u32, 0.0)));
        let inside = from_triplets(n_own, n_own, inside);

        // Interface blocks: rows by global id per row owner.
        let mut by_rank: Vec<Vec<(u64, u32, f64)>> = vec![Vec::new(); size];
        for (ro, gi, pj, v) in iface {
            by_rank[ro].push((gi, pj, v));
        }
        let mut send_if = Vec::new();
        let mut row_gids: Vec<Vec<u64>> = vec![Vec::new(); size];
        for (r, mut ents) in by_rank.into_iter().enumerate() {
            if ents.is_empty() {
                continue;
            }
            ents.sort_by_key(|&(g, c, _)| (g, c));
            let mut rows: Vec<u64> = ents.iter().map(|e| e.0).collect();
            rows.dedup();
            let row_of: HashMap<u64, u32> = rows
                .iter()
                .enumerate()
                .map(|(k, &g)| (g, k as u32))
                .collect();
            let m = from_triplets(
                rows.len(),
                n_own,
                ents.iter().map(|&(g, c, v)| (row_of[&g] as usize, c, v)),
            );
            send_if.push(IfBlock {
                rank: r,
                row_ptr: m.row_ptr,
                cols: m.col_idx,
                vals: m.values,
            });
            row_gids[r] = rows;
        }

        // Row owners learn which of their rows each neighbour will send.
        let rows_in = comm.all_to_all_u64(&row_gids);
        let mut recv_if = Vec::new();
        let mut bad = None;
        for (src, gids) in rows_in.iter().enumerate() {
            if gids.is_empty() {
                continue;
            }
            let p: Vec<u32> = gids
                .iter()
                .map(|g| {
                    owned_by_gid.get(g).copied().unwrap_or_else(|| {
                        bad.get_or_insert(format!("rank {src} computes row {g}, not owned here"));
                        NO_DOF
                    })
                })
                .collect();
            recv_if.push((src, p));
        }
        agree(comm, bad.map_or(Ok(()), Err))?;
        Ok(DistCsr {
            inside,
            send_if,
            recv_if,
        })
    }

    /// Owned rows by owned columns, with every contribution added up and a
    /// stored diagonal in each row. The block the preconditioner factors.
    pub fn inside(&self) -> &CsrMatrix<f64> {
        &self.inside
    }

    pub fn n_owned(&self) -> usize {
        self.inside.nrows
    }

    /// Ranks this one sends interface partial sums to and receives them from.
    pub fn neighbours(&self) -> (Vec<usize>, Vec<usize>) {
        (
            self.send_if.iter().map(|b| b.rank).collect(),
            self.recv_if.iter().map(|(r, _)| *r).collect(),
        )
    }

    /// `y = A x` on owned vectors. Collective.
    pub fn matvec(&self, comm: &dyn Exchange, x: &[f64], y: &mut [f64]) {
        let send: Vec<(usize, Vec<f64>)> = self
            .send_if
            .iter()
            .map(|b| {
                let part = (0..b.row_ptr.len() - 1)
                    .map(|j| {
                        (b.row_ptr[j]..b.row_ptr[j + 1])
                            .fold(0.0, |s, k| s + b.vals[k] * x[b.cols[k] as usize])
                    })
                    .collect();
                (b.rank, part)
            })
            .collect();
        let mut recv: Vec<(usize, Vec<f64>)> = self
            .recv_if
            .iter()
            .map(|(r, p)| (*r, vec![0.0; p.len()]))
            .collect();
        // The local product defines y; partial sums from neighbours add to it.
        elmer_sparse::spmv(&self.inside, x, y);
        comm.exchange_f64(&send, &mut recv);
        for ((_, pos), (_, vals)) in self.recv_if.iter().zip(&recv) {
            for (&p, &v) in pos.iter().zip(vals) {
                y[p as usize] += v;
            }
        }
    }

    /// The matrix as a [`LinearOperator`] on owned vectors, for HUTI.
    pub fn operator<'a>(&'a self, comm: &'a dyn Exchange) -> DistOperator<'a> {
        DistOperator { a: self, comm }
    }
}

/// [`DistCsr`] bound to its communicator. `apply` is collective, so every
/// rank must call it the same number of times, which the HUTI methods do.
pub struct DistOperator<'a> {
    a: &'a DistCsr,
    comm: &'a dyn Exchange,
}

impl LinearOperator<f64> for DistOperator<'_> {
    fn nrows(&self) -> usize {
        self.a.n_owned()
    }

    fn ncols(&self) -> usize {
        self.a.n_owned()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.a.matvec(self.comm, x, y);
    }
}
