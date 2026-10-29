//! Local-to-global glue: `UpdateGlobalEquations` (SolverBasics.F90:667)
//! and `CRS_GlueLocalMatrix` (CRSMatrix.F90:717).

use elmer_core::{CsrMatrix, GlobalAssembler, LocalSystem, NO_DOF};

/// Adds element systems into a CSR matrix and right-hand side in Elmer's
/// order, so sums come out bitwise equal for the same element order.
///
/// `add_local(dof_index, local)` takes Elmer's `Indexes`: one entry per
/// element node, the node's [`Perm`](crate::Perm) slot (or [`NO_DOF`] to
/// skip it). The number of components per node is
/// `local.n / dof_index.len()`, and local row `dofs * i + c` is component
/// `c` of node `i` (Elmer's interleaved layout).
///
/// Entries missing from the pattern, and rows without a diagonal, are
/// dropped silently, as in Elmer. Normal-tangential rotation and passive
/// elements are the caller's business.
pub struct CsrAssembler<'a> {
    pub a: &'a mut CsrMatrix<f64>,
    pub rhs: &'a mut [f64],
    /// Optional mass and damping values on the pattern of `a` (Elmer
    /// `MassValues`, `DampValues`), fed from `local.mass` / `local.damp`.
    pub mass: Option<&'a mut [f64]>,
    pub damp: Option<&'a mut [f64]>,
}

impl<'a> CsrAssembler<'a> {
    pub fn new(a: &'a mut CsrMatrix<f64>, rhs: &'a mut [f64]) -> Self {
        assert_eq!(rhs.len(), a.nrows, "rhs length must match matrix rows");
        assert!(a.diag.is_some(), "glue needs the diagonal index");
        Self {
            a,
            rhs,
            mass: None,
            damp: None,
        }
    }

    pub fn with_mass(mut self, mass: &'a mut [f64]) -> Self {
        assert_eq!(
            mass.len(),
            self.a.nnz(),
            "mass values must match the pattern"
        );
        self.mass = Some(mass);
        self
    }

    pub fn with_damp(mut self, damp: &'a mut [f64]) -> Self {
        assert_eq!(
            damp.len(),
            self.a.nnz(),
            "damp values must match the pattern"
        );
        self.damp = Some(damp);
        self
    }
}

impl GlobalAssembler for CsrAssembler<'_> {
    fn add_local(&mut self, dof_index: &[u32], local: &LocalSystem) {
        let nd = dof_index.len();
        if nd == 0 {
            return;
        }
        assert!(
            local.n.is_multiple_of(nd),
            "local size {} is not a multiple of {nd} nodes",
            local.n
        );
        let dofs = local.n / nd;
        let CsrMatrix {
            row_ptr,
            col_idx,
            values,
            diag,
            ..
        } = &mut *self.a;
        let pat = Pattern {
            row_ptr,
            col_idx,
            diag: diag.as_deref().expect("glue needs the diagonal index"),
        };
        pat.glue(values, dof_index, dofs, local.n, &local.stiff);
        if let (Some(m), Some(lm)) = (self.mass.as_deref_mut(), local.mass.as_deref()) {
            pat.glue(m, dof_index, dofs, local.n, lm);
        }
        if let (Some(d), Some(ld)) = (self.damp.as_deref_mut(), local.damp.as_deref()) {
            pat.glue(d, dof_index, dofs, local.n, ld);
        }
        for (i, &idx) in dof_index.iter().enumerate() {
            if idx == NO_DOF {
                continue;
            }
            for j in 0..dofs {
                self.rhs[dofs * idx as usize + j] += local.force[dofs * i + j];
            }
        }
    }
}

struct Pattern<'p> {
    row_ptr: &'p [usize],
    col_idx: &'p [u32],
    diag: &'p [usize],
}

impl Pattern<'_> {
    /// Position of `(row, col)`, searching the half of the row on `col`'s
    /// side of the diagonal as Elmer's `BinarySearch` does.
    fn find(&self, row: usize, col: usize) -> Option<usize> {
        let d = self.diag[row];
        let (lo, hi) = if col >= row {
            (d, self.row_ptr[row + 1])
        } else {
            (self.row_ptr[row], d)
        };
        self.col_idx[lo..hi]
            .binary_search(&(col as u32))
            .ok()
            .map(|p| lo + p)
    }

    /// `CRS_GlueLocalMatrix` loop order: rows by node, then component from
    /// last to first; columns likewise.
    fn glue(&self, values: &mut [f64], idx: &[u32], dofs: usize, n: usize, local: &[f64]) {
        for (i, &ii) in idx.iter().enumerate() {
            if ii == NO_DOF {
                continue;
            }
            for k in 0..dofs {
                let row = dofs * ii as usize + (dofs - 1 - k);
                if self.diag[row] == usize::MAX {
                    continue;
                }
                let lrow = dofs * i + (dofs - 1 - k);
                for (j, &jj) in idx.iter().enumerate() {
                    if jj == NO_DOF {
                        continue;
                    }
                    for l in 0..dofs {
                        let col = dofs * jj as usize + (dofs - 1 - l);
                        if let Some(c) = self.find(row, col) {
                            values[c] += local[lrow * n + dofs * j + (dofs - 1 - l)];
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::list::ListMatrix;
    use crate::perm::initialize_matrix;

    #[test]
    fn glue_two_dofs_skips_absent_nodes() {
        let mut l = ListMatrix::new(2);
        l.add_clique(&mut [0, 1]);
        let mut a = initialize_matrix(&l, None, 2);
        let mut b = vec![0.0; 4];
        let mut local = LocalSystem::zeros(6);
        for r in 0..6 {
            local.force[r] = r as f64;
            for c in 0..6 {
                local.stiff[r * 6 + c] = (10 * r + c) as f64;
            }
        }
        // Element nodes: slot 1, absent, slot 0.
        CsrAssembler::new(&mut a, &mut b).add_local(&[1, NO_DOF, 0], &local);
        // Global row 0 (slot 0, comp 0) = local row 4; cols slot0 c0,c1 =
        // local 4,5; slot1 c0,c1 = local 0,1.
        let row0: Vec<f64> = a.values[a.row_ptr[0]..a.row_ptr[1]].to_vec();
        assert_eq!(row0, vec![44.0, 45.0, 40.0, 41.0]);
        assert_eq!(b, vec![4.0, 5.0, 0.0, 1.0]);
    }
}
