//! `DisplaceMesh` (MeshBasics.F90): StressSolve moves the mesh nodes by its
//! displacement after each solve and back at the start of the next one.
//! The runner owns the node coordinates and installs a [`NodeMover`] for
//! the solvers of its thread, as `elmer-physics-flow`'s free surfaces do.

use std::cell::RefCell;
use std::sync::Arc;

/// Replaces every node coordinate of the mesh the solvers read.
pub type NodeMover = Arc<dyn Fn(&[[f64; 3]])>;

thread_local! {
    static MOVER: RefCell<Option<NodeMover>> = const { RefCell::new(None) };
    static DONE: RefCell<Vec<&'static str>> = const { RefCell::new(Vec::new()) };
}

/// Installs the runner's node mover for the solvers run on this thread,
/// returning the previous one. It starts a new run: no routine has
/// displaced the mesh yet.
pub fn install_node_mover(m: Option<NodeMover>) -> Option<NodeMover> {
    DONE.with(|d| d.borrow_mut().clear());
    MOVER.with(|c| std::mem::replace(&mut *c.borrow_mut(), m))
}

/// The routine's `SAVE`d `AllocationsDone`, shared by every instance and
/// kept over `RUN` blocks as the Fortran's is, then set: whether an earlier
/// call of `routine` left its displacement on the mesh to be taken off.
pub fn allocations_done(routine: &'static str) -> bool {
    DONE.with(|d| {
        let mut d = d.borrow_mut();
        if d.contains(&routine) {
            true
        } else {
            d.push(routine);
            false
        }
    })
}

pub(crate) fn mover() -> Option<NodeMover> {
    MOVER.with(|c| c.borrow().clone())
}

/// `DisplaceMesh(Mesh, Update, sgn, Perm, DOFs)` without `UpdateDirs`:
/// `x(i) += sgn * Update(DOFs*(Perm(i)-1)+1)` and so on for the first
/// `DOFs` (at most three) coordinates of every node with a dof. `coords`
/// are the current coordinates; `None` without an installed mover.
pub fn displace_mesh(coords: &[[f64; 3]], update: &[f64], sgn: f64, perm: &[u32], dofs: usize, dirs: usize) -> Option<()> {
    let mover = mover()?;
    let mut x = coords.to_vec();
    for (i, c) in x.iter_mut().enumerate().take(perm.len()) {
        let k = perm[i];
        if k == elmer_core::NO_DOF {
            continue;
        }
        let k = dofs * k as usize;
        for (d, cd) in c.iter_mut().enumerate().take(dirs.min(3)) {
            *cd += sgn * update[k + d];
        }
    }
    mover(&x);
    Some(())
}
