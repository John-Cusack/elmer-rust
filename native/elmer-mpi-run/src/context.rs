//! The rank's parallel state for the duration of one run.
//!
//! `elmer-solver` reaches the linear solvers through plain function
//! factories and hands assembly backends only a matrix, so the transport,
//! the node ownership and the dof layout of every matrix it created live
//! here, one per thread (one per rank). [`RankContext::install`] sets it,
//! the returned guard clears it.

use std::cell::RefCell;
use std::rc::Rc;

use elmer_core::{CsrMatrix, NO_DOF};
use elmer_dist::{DistError, DofLayout, Exchange};

/// A matrix `elmer-solver` created, and the dofs behind its rows.
struct Registered {
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
    /// Global dof id and owner of every row.
    gid: Vec<u64>,
    owner: Vec<usize>,
    layout: Option<Rc<DofLayout>>,
}

/// Transport and ownership of one rank.
pub struct RankContext {
    pub ex: Box<dyn Exchange>,
    node_gid: Vec<u64>,
    node_owner: Vec<usize>,
    matrices: RefCell<Vec<Registered>>,
}

thread_local! {
    static CONTEXT: RefCell<Option<Rc<RankContext>>> = const { RefCell::new(None) };
}

/// Clears the thread's [`RankContext`] when dropped.
pub struct ContextGuard(());

impl Drop for ContextGuard {
    fn drop(&mut self) {
        CONTEXT.with(|c| c.borrow_mut().take());
    }
}

impl RankContext {
    /// Makes this the thread's context until the guard drops. `node_gid`
    /// and `node_owner` describe the local mesh nodes.
    pub fn install(ex: Box<dyn Exchange>, node_gid: Vec<u64>, node_owner: Vec<usize>) -> (Rc<Self>, ContextGuard) {
        let ctx = Rc::new(RankContext { ex, node_gid, node_owner, matrices: RefCell::new(Vec::new()) });
        CONTEXT.with(|c| *c.borrow_mut() = Some(Rc::clone(&ctx)));
        (ctx, ContextGuard(()))
    }

    /// The thread's context, if a parallel run is in progress.
    pub fn current() -> Option<Rc<Self>> {
        CONTEXT.with(|c| c.borrow().clone())
    }

    /// Records a matrix created for `perm` (node -> slot, [`NO_DOF`] when
    /// absent) with `dofs` dofs per slot. Global dof `dofs * (gid - 1) + k`
    /// is Elmer's nodal `GlobalDOFs`; the node's owner owns its dofs.
    pub fn register(&self, perm: &[u32], dofs: usize, a: &CsrMatrix<f64>) {
        let mut gid = vec![0u64; a.nrows];
        let mut owner = vec![0usize; a.nrows];
        let d = dofs as u64;
        for (node, &slot) in perm.iter().enumerate().take(self.node_gid.len()) {
            if slot == NO_DOF {
                continue;
            }
            for k in 0..dofs {
                let row = slot as usize * dofs + k;
                if row < a.nrows {
                    gid[row] = d * (self.node_gid[node] - 1) + k as u64;
                    owner[row] = self.node_owner[node];
                }
            }
        }
        let mut m = self.matrices.borrow_mut();
        m.retain(|r| r.row_ptr != a.row_ptr || r.col_idx != a.col_idx);
        m.push(Registered {
            row_ptr: a.row_ptr.clone(),
            col_idx: a.col_idx.clone(),
            gid,
            owner,
            layout: None,
        });
    }

    /// The dof layout of a matrix registered with [`register`](Self::register),
    /// found by its sparsity pattern. Collective on every call (one max
    /// reduction of a found flag), and again on first use of a matrix.
    pub fn layout_of(&self, a: &CsrMatrix<f64>) -> Result<Rc<DofLayout>, DistError> {
        let found = {
            let m = self.matrices.borrow();
            m.iter().position(|r| r.row_ptr == a.row_ptr && r.col_idx == a.col_idx)
        };
        // Collective, as elmer-dist's `agree`: a rank that does not know the
        // matrix must not return alone while the others go on into
        // `DofLayout::new` (or the solve's reductions) and wait for it
        // forever. Every rank learns that one failed, and all return an error.
        let mut missing = [if found.is_none() { 1.0 } else { 0.0 }];
        self.ex.max_f64(&mut missing);
        let rank = self.ex.rank();
        let Some(i) = found else {
            return Err(DistError::Local { rank, msg: "a linear system whose matrix the parallel assembly did not create".into() });
        };
        if missing[0] != 0.0 {
            return Err(DistError::Remote { rank });
        }
        if let Some(l) = &self.matrices.borrow()[i].layout {
            return Ok(Rc::clone(l));
        }
        let (gid, owner) = {
            let m = self.matrices.borrow();
            (m[i].gid.clone(), m[i].owner.clone())
        };
        let layout = Rc::new(DofLayout::new(&*self.ex, gid, owner)?);
        self.matrices.borrow_mut()[i].layout = Some(Rc::clone(&layout));
        Ok(layout)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Communicator;

    /// Rank 0 of two; the other rank's contribution to a max reduction is
    /// `other`. Point-to-point exchange is never reached by these tests.
    struct TwoRanks {
        other: f64,
    }

    impl Communicator for TwoRanks {
        fn rank(&self) -> usize {
            0
        }
        fn size(&self) -> usize {
            2
        }
        fn sum_f64(&self, buf: &mut [f64]) {
            buf.iter_mut().for_each(|x| *x += self.other);
        }
        fn max_f64(&self, buf: &mut [f64]) {
            buf.iter_mut().for_each(|x| *x = x.max(self.other));
        }
        fn barrier(&self) {}
    }

    impl Exchange for TwoRanks {
        fn all_to_all_u64(&self, _send: &[Vec<u64>]) -> Vec<Vec<u64>> {
            unreachable!("the agreement fails before any exchange")
        }
        fn exchange_f64(&self, _send: &[(usize, Vec<f64>)], _recv: &mut [(usize, Vec<f64>)]) {
            unreachable!("the agreement fails before any exchange")
        }
    }

    fn matrix() -> CsrMatrix<f64> {
        CsrMatrix { nrows: 1, ncols: 1, row_ptr: vec![0, 1], col_idx: vec![0], values: vec![1.0], diag: None }
    }

    /// The radiation_np2 hang: one rank's matrix was changed outside the
    /// parallel assembly. Every rank must get an error, not only that one.
    #[test]
    fn a_matrix_unknown_on_any_rank_fails_every_rank() {
        // Unknown here: this rank's own error.
        let (ctx, _g) = RankContext::install(Box::new(TwoRanks { other: 0.0 }), vec![1], vec![0]);
        assert!(matches!(ctx.layout_of(&matrix()), Err(DistError::Local { rank: 0, .. })));
        drop(_g);
        // Known here but unknown on the other rank: a remote failure, before
        // the collective DofLayout::new that the other rank would never join.
        let (ctx, _g) = RankContext::install(Box::new(TwoRanks { other: 1.0 }), vec![1], vec![0]);
        ctx.register(&[0], 1, &matrix());
        assert!(matches!(ctx.layout_of(&matrix()), Err(DistError::Remote { rank: 0 })));
    }
}
