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
    /// found by its sparsity pattern. Collective on first use.
    pub fn layout_of(&self, a: &CsrMatrix<f64>) -> Result<Rc<DofLayout>, DistError> {
        let found = {
            let m = self.matrices.borrow();
            m.iter().position(|r| r.row_ptr == a.row_ptr && r.col_idx == a.col_idx)
        };
        let local = |msg: &str| DistError::Local { rank: self.ex.rank(), msg: msg.into() };
        let Some(i) = found else {
            return Err(local("a linear system whose matrix the parallel assembly did not create"));
        };
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
