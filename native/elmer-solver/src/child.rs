//! Child solvers (`CreateChildSolver`, MainUtils.F90:2845): a container
//! solver that shares its parent's mesh, active elements and permutation
//! and owns nothing but a matrix and a right-hand side. It is never run
//! (`SOLVER_EXEC_NEVER`); the parent assembles into it with
//! `DefaultInitialize(USolver=Child)`, `DefaultUpdateEquations(...,
//! USolver=Child)` and `DefaultDirichletBCs(USolver=Child)`, and the block
//! solver finds its matrix by name through `Block Matrix Schur Variable`
//! (`BlockPrecMatrix`, BlockSolve.F90:2523). IncompressibleNS assembles the
//! pressure block of its block preconditioner this way.
//!
//! The child variable is kept with the child, not in the field registry:
//! Elmer creates it with `Variable Output = False` (the parent's `schur:
//! Variable Output`), so it reaches no result file.

use elmer_core::{CsrMatrix, ElemId, LocalSystem, Params, NO_DOF};

use crate::sim::LinearSystem;
use crate::{SolverCtx, SolverError};

/// A child solver's matrix and right-hand side, `dofs` rows per slot of the
/// parent's permutation.
#[derive(Clone, Debug)]
pub(crate) struct ChildSystem {
    /// The child variable's name, as `Block Matrix Schur Variable` names it.
    pub name: String,
    pub dofs: usize,
    pub system: LinearSystem,
}

/// `CRS_CreateChildMatrix` (CRSMatrix.F90) with `ColDofs = Dofs` and no
/// `Diagonal`: the topology of the parent's first row of each node, every
/// `parent_dofs`-th column, expanded to `dofs` rows and columns per node.
/// The parent must be node-major with a full `parent_dofs` x `parent_dofs`
/// block per coupled node pair, as `InitializeMatrix` makes it.
pub(crate) fn child_matrix(
    parent: &CsrMatrix<f64>,
    parent_dofs: usize,
    dofs: usize,
) -> CsrMatrix<f64> {
    if parent_dofs == dofs {
        // Reusing (or copying) the initial matrix topology.
        let mut m = parent.clone();
        m.values.iter_mut().for_each(|v| *v = 0.0);
        return m;
    }
    let nodes = parent.nrows / parent_dofs;
    let mut row_ptr = Vec::with_capacity(nodes * dofs + 1);
    let mut col_idx = Vec::new();
    row_ptr.push(0);
    for i in 0..nodes {
        let r = i * parent_dofs;
        let cols: Vec<u32> = (parent.row_ptr[r]..parent.row_ptr[r + 1])
            .step_by(parent_dofs)
            .map(|j| parent.col_idx[j] / parent_dofs as u32)
            .collect();
        for _ in 0..dofs {
            for &nn in &cols {
                col_idx.extend((0..dofs as u32).map(|l| nn * dofs as u32 + l));
            }
            row_ptr.push(col_idx.len());
        }
    }
    let n = nodes * dofs;
    let diag = (0..n)
        .map(|i| {
            (row_ptr[i]..row_ptr[i + 1])
                .find(|&j| col_idx[j] as usize == i)
                .unwrap_or(row_ptr[i])
        })
        .collect();
    CsrMatrix {
        nrows: n,
        ncols: n,
        values: vec![0.0; col_idx.len()],
        row_ptr,
        col_idx,
        diag: Some(diag),
    }
}

impl SolverCtx<'_> {
    fn child_index(&self, name: &str) -> Option<usize> {
        self.sim.slots[self.slot]
            .children
            .iter()
            .position(|c| c.name.eq_ignore_ascii_case(name))
    }

    fn child_mut(&mut self, name: &str, caller: &str) -> Result<&mut ChildSystem, SolverError> {
        let ix = self.child_index(name).ok_or_else(|| {
            SolverError::fatal(caller, format!("no child solver for variable: {name}"))
        })?;
        Ok(&mut self.sim.slots[self.slot].children[ix])
    }

    /// `CreateChildSolver(Solver, name, dofs, prefix)` followed by
    /// `DefaultInitialize(USolver=Child)`: on the first call the child's
    /// matrix is made from the parent's topology (`CreateChildMatrix`),
    /// so the parent's own [`default_initialize`](Self::default_initialize)
    /// must have run; every call zeroes the child's matrix and right-hand
    /// side. A solver whose dofs come from a [`crate::DofsBackend`] works the
    /// same way: the child shares the parent's slots.
    pub fn child_initialize(&mut self, name: &str, dofs: usize) -> Result<(), SolverError> {
        if dofs == 0 {
            return Err(SolverError::fatal(
                "CreateChildSolver",
                "a child solver of zero dofs",
            ));
        }
        if self.child_index(name).is_none() {
            let slot = &self.sim.slots[self.slot];
            let id = slot.variable.ok_or_else(|| {
                SolverError::fatal("CreateChildSolver", "Parent solver is missing variable!")
            })?;
            let parent_dofs = self.sim.vars.by_id(id).dofs;
            let parent = slot.system.as_ref().ok_or_else(|| {
                SolverError::fatal(
                    "CreateChildSolver",
                    "Parent matrix needed for child matrix!",
                )
            })?;
            if parent.matrix.nrows % parent_dofs != 0 {
                return Err(SolverError::fatal(
                    "CreateChildSolver",
                    "parent matrix is not node-major",
                ));
            }
            self.info(
                "CreateChildSolver",
                &format!(
                    "Creating solver of size {} for variable: {name}",
                    crate::fortran_fmt::i2s(dofs as i64)
                ),
                6,
            );
            let matrix = child_matrix(&parent.matrix, parent_dofs, dofs);
            let n = matrix.nrows;
            self.sim.slots[self.slot].children.push(ChildSystem {
                name: name.to_string(),
                dofs,
                system: LinearSystem {
                    matrix,
                    rhs: vec![0.0; n],
                    mass: None,
                },
            });
        }
        let child = self.child_mut(name, "DefaultInitialize")?;
        if child.dofs != dofs {
            return Err(SolverError::fatal(
                "CreateChildSolver",
                format!("child solver {name} exists with {} dofs", child.dofs),
            ));
        }
        child.system.matrix.values.iter_mut().for_each(|v| *v = 0.0);
        child.system.rhs.iter_mut().for_each(|v| *v = 0.0);
        Ok(())
    }

    /// `DefaultUpdateEquations(local, UElement=e, USolver=Child)`: adds a
    /// local system of the child's dofs per element slot at element `e`.
    pub fn child_update(
        &mut self,
        name: &str,
        local: &LocalSystem,
        e: ElemId,
    ) -> Result<(), SolverError> {
        if self.is_passive(e)? {
            return Ok(());
        }
        let mut idx = Vec::new();
        self.element_slots(e, &mut idx);
        let dofs = self.child_mut(name, "DefaultUpdateEquations")?.dofs;
        if local.n != idx.len() * dofs {
            return Err(SolverError::fatal(
                "DefaultUpdateEquations",
                format!(
                    "local system of {} rows for {} element dofs of child {name}",
                    local.n,
                    idx.len()
                ),
            ));
        }
        if idx.contains(&NO_DOF) {
            return Err(SolverError::fatal(
                "DefaultUpdateEquations",
                format!("element {} has nodes outside the child system", e.0),
            ));
        }
        let sim = &mut *self.sim;
        let backend = sim
            .backends
            .assembly
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
        let slot = &mut sim.slots[self.slot];
        let ix = slot
            .children
            .iter()
            .position(|c| c.name.eq_ignore_ascii_case(name))
            .expect("looked up above");
        let sys = &mut slot.children[ix].system;
        backend.add_local(&mut sys.matrix, &mut sys.rhs, &idx, local);
        Ok(())
    }

    /// `DefaultDirichletBCs(USolver=Child)`: the child variable's own
    /// Dirichlet conditions. A child is a container whose variable no
    /// section sets, so this is a no-op; a boundary condition that does name
    /// it is refused rather than ignored.
    pub fn child_dirichlet(&mut self, name: &str) -> Result<(), SolverError> {
        let dofs = self.child_mut(name, "DefaultDirichletBCs")?.dofs;
        let key = elmer_core::normalize_key(name);
        let mut keys = vec![key.clone()];
        keys.extend((1..=dofs).map(|k| format!("{key} {k}")));
        let named = self
            .sim
            .model
            .boundary_conditions
            .iter()
            .chain(self.sim.model.body_forces.iter())
            .any(|s| keys.iter().any(|k| s.get(k).is_some()));
        if named {
            return Err(SolverError::Unsupported(format!(
                "Dirichlet conditions of the child solver variable {name}"
            )));
        }
        Ok(())
    }

    /// The matrix of the child solver for variable `name`, if this solver has one.
    pub fn child_system(&self, name: &str) -> Option<&LinearSystem> {
        self.child_index(name)
            .map(|ix| &self.sim.slots[self.slot].children[ix].system)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn child_matrix_collapses_the_node_blocks() {
        // Two nodes coupled, 2 dofs each: a full 4x4 pattern.
        let parent = CsrMatrix {
            nrows: 4,
            ncols: 4,
            row_ptr: vec![0, 4, 8, 12, 16],
            col_idx: vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
            values: vec![1.0; 16],
            diag: None,
        };
        let c = child_matrix(&parent, 2, 1);
        assert_eq!(c.nrows, 2);
        assert_eq!(c.row_ptr, vec![0, 2, 4]);
        assert_eq!(c.col_idx, vec![0, 1, 0, 1]);
        assert_eq!(c.values, vec![0.0; 4]);
        assert_eq!(c.diag, Some(vec![0, 3]));
    }
}
