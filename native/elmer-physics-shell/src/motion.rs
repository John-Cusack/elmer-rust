//! `DisplaceMesh` for ShellSolver's `Displace Mesh`: the runner owns the
//! node coordinates and installs a [`NodeMover`] for the solvers of its
//! thread, as it does for StressSolve and the free surfaces.

use std::cell::RefCell;
use std::sync::Arc;

/// Replaces every node coordinate of the mesh the solvers read.
pub type NodeMover = Arc<dyn Fn(&[[f64; 3]])>;

thread_local! {
    static MOVER: RefCell<Option<NodeMover>> = const { RefCell::new(None) };
}

/// Installs the runner's node mover for the solvers run on this thread,
/// returning the previous one.
pub fn install_node_mover(m: Option<NodeMover>) -> Option<NodeMover> {
    MOVER.with(|c| std::mem::replace(&mut *c.borrow_mut(), m))
}

/// `DisplaceMesh( Mesh, Update, sgn, Perm, DOFs, .FALSE., 3 )`: the first
/// three components of each node's `DOFs` added `sgn` times. `None`
/// without an installed mover.
pub(crate) fn displace_mesh(
    coords: &[[f64; 3]],
    update: &[f64],
    sgn: f64,
    perm: &[u32],
    dofs: usize,
) -> Option<()> {
    let mover = MOVER.with(|c| c.borrow().clone())?;
    let mut x = coords.to_vec();
    for (i, c) in x.iter_mut().enumerate().take(perm.len()) {
        let k = perm[i];
        if k == elmer_core::NO_DOF {
            continue;
        }
        let k = dofs * k as usize;
        for d in 0..3 {
            c[d] += sgn * update[k + d];
        }
    }
    mover(&x);
    Some(())
}
