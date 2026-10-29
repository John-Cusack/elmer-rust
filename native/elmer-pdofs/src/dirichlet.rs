//! The p-element part of `DefaultDirichletBCs` (DefUtils.F90:6355).
//!
//! Before the nodal values are set, Elmer walks the boundary elements of
//! every BC that gives the field: for an active p parent, the edge DOFs of
//! the boundary element (from `mGetBoundaryIndexesFromParent`) are
//! constrained to zero when the BC value is not a procedure. The nodal
//! rows are then set as for Lagrange elements
//! (`elmer_assembly::DirichletSet::set_boundaries`), the node slots
//! being the node indices when there is one DOF per node.
//!
//! Not covered: BC values given by a procedure (`ZeroRow` here, then the
//! `LocalBcBDOFs` projection of the edge DOFs), and 3D face DOFs.

use elmer_assembly::{DirichletSet, Perm};
use elmer_core::NO_DOF;

use crate::dofs::DofLayout;
use crate::Error;

/// Constrains the edge DOFs of the boundary elements `on_bc(e)` selects
/// (attached to a BC that gives this field, and active for the solver) to
/// zero, for component `dof` of `dofs`. `constant` is false when the BC
/// value is a procedure, which is refused.
pub fn constrain_edge_dofs(
    lay: &DofLayout<'_>,
    perm: &Perm,
    dof: usize,
    dofs: usize,
    on_bc: &dyn Fn(usize) -> bool,
    constant: bool,
    set: &mut DirichletSet,
) -> Result<(), Error> {
    let pm = lay.mesh;
    let mut ind = Vec::new();
    for (e, el) in pm.elements.iter().enumerate().skip(pm.n_bulk) {
        let Some(parent) = el.parents[0].or(el.parents[1]) else {
            continue;
        };
        if !on_bc(e) || !lay.is_active_p(parent as usize) {
            continue;
        }
        if !constant {
            return Err(Error::Unsupported(
                "p-element Dirichlet values given by a procedure (LocalBcBDOFs)".into(),
            ));
        }
        lay.boundary_indexes_from_parent(e, &mut ind);
        for &g in ind.iter().skip(el.nodes.len()) {
            let nb = perm.perm.get(g as usize).copied().unwrap_or(NO_DOF);
            if nb != NO_DOF {
                set.set(dofs * nb as usize + dof, 0.0);
            }
        }
    }
    Ok(())
}
