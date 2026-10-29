//! Dirichlet conditions with Elmer semantics.
//!
//! Elmer splits this in two. `DefaultDirichletBCs` (DefUtils.F90:6200)
//! and `SetDirichletBoundaries` (BoundaryConditionUtils.F90:297) only
//! *record* constrained rows and their values (`ConstrainedDOF`,
//! `DValues`); [`DirichletSet`] does that for nodal BCs on boundary
//! elements, given values already resolved per node. Then
//! `EnforceDirichletConditions` (BoundaryConditionUtils.F90:3269) edits
//! the matrix and right-hand side; that is [`apply_dirichlet`].
//!
//! Not covered (callers must refuse such cases): normal-tangential and
//! periodic BCs, body-force and point (`Target Nodes`/`Target
//! Coordinates`) Dirichlet values, `Passive` BCs, edge/face DOFs, mass and
//! damping zeroing for eigen/harmonic analysis, block off-diagonal
//! matrices, and the `row equilibration`/`constant` scaling methods.

use elmer_core::{CsrMatrix, ElemId, LinalgError, MeshView};

use crate::perm::Perm;

/// How `EnforceDirichletConditions` scales the diagonal it writes
/// (`Linear System Dirichlet Scaling`, falling back to `Linear System
/// Scaling`, default true; `Linear System Scaling Method`, default
/// `diagonal`).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DirichletScaling {
    /// Scaling off, or method `none`: diagonal 1, rhs = value.
    None,
    /// Method `diagonal` (the default): the diagonal stays (to rounding)
    /// `|a_kk|`, computed as Elmer does through `1/sqrt|a_kk|`.
    Diagonal,
}

/// Matrix-side options of `EnforceDirichletConditions`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct DirichletMode {
    /// `Linear System Symmetric`: also eliminate the constrained columns
    /// from unconstrained rows, moving them to the rhs
    /// (`CRS_ElimSymmDirichlet`). Off: only the constrained rows change.
    pub symmetric: bool,
    pub scaling: DirichletScaling,
}

impl Default for DirichletMode {
    /// Elmer defaults: not symmetric, diagonal scaling.
    fn default() -> Self {
        Self {
            symmetric: false,
            scaling: DirichletScaling::Diagonal,
        }
    }
}

/// One BC's contribution to one component, with values resolved per node.
pub struct DirichletBc<'a> {
    /// Boundary elements this BC is attached to (Elmer: `Constraint ==`
    /// the BC's tag after `Target Boundaries` mapping).
    pub elements: &'a dyn Fn(ElemId) -> bool,
    /// Value at a node of an element (`ListGetReal(BC, Name, n, Indexes)`).
    pub value: &'a dyn Fn(ElemId, u32) -> f64,
    /// `<Name> Condition`: the node is skipped where this is negative.
    pub condition: Option<&'a dyn Fn(ElemId, u32) -> f64>,
}

/// Recorded constraints, Elmer's `ConstrainedDOF` and `DValues`.
#[derive(Clone, Debug, PartialEq)]
pub struct DirichletSet {
    pub constrained: Vec<bool>,
    pub values: Vec<f64>,
}

impl DirichletSet {
    pub fn new(nrows: usize) -> Self {
        Self {
            constrained: vec![false; nrows],
            values: vec![0.0; nrows],
        }
    }

    /// `SetSinglePoint`: constrain `row` to `value`; a later call wins.
    pub fn set(&mut self, row: usize, value: f64) {
        self.constrained[row] = true;
        self.values[row] = value;
    }

    /// `SetDirichletBoundaries` + `SetElementValues` for component `dof`
    /// (0-based) of a field with `dofs` components, BCs given in BC-number
    /// order. Elmer's default order is element-major (for each boundary
    /// element, each BC), so on shared nodes the later boundary element
    /// wins; with `by_bc_numbering` (`Set Dirichlet BCs by BC Numbering`)
    /// it is BC-major and the later BC wins.
    pub fn set_boundaries(
        &mut self,
        mesh: &dyn MeshView,
        perm: &Perm,
        dof: usize,
        dofs: usize,
        bcs: &[DirichletBc<'_>],
        by_bc_numbering: bool,
    ) {
        let nb = mesh.n_bulk() as u32;
        let boundary = (nb..nb + mesh.n_boundary() as u32).map(ElemId);
        if by_bc_numbering {
            for bc in bcs {
                for e in boundary.clone() {
                    self.set_element(mesh, perm, dof, dofs, bc, e);
                }
            }
        } else {
            for e in boundary {
                for bc in bcs {
                    self.set_element(mesh, perm, dof, dofs, bc, e);
                }
            }
        }
    }

    fn set_element(
        &mut self,
        mesh: &dyn MeshView,
        perm: &Perm,
        dof: usize,
        dofs: usize,
        bc: &DirichletBc<'_>,
        e: ElemId,
    ) {
        if !(bc.elements)(e) {
            return;
        }
        for &node in mesh.element(e).nodes {
            if let Some(cond) = bc.condition {
                if cond(e, node) < 0.0 {
                    continue;
                }
            }
            if let Some(row) = perm.row(node, dof, dofs) {
                self.set(row, (bc.value)(e, node));
            }
        }
    }

    /// Constrained rows and values, ascending by row.
    pub fn to_list(&self) -> Vec<(u32, f64)> {
        self.constrained
            .iter()
            .zip(&self.values)
            .enumerate()
            .filter(|(_, (&c, _))| c)
            .map(|(i, (_, &v))| (i as u32, v))
            .collect()
    }
}

/// `EnforceDirichletConditions` on a CSR system: constrain each
/// `(row, value)` (a repeated row takes its last value).
///
/// In order: compute the diagonal scaling of the unconstrained matrix;
/// with `mode.symmetric`, move constrained columns of the other rows to
/// the rhs (`b_n -= a_nk * v_k`, column order); then zero each
/// constrained row and set its diagonal to `s` and its rhs to `s * v`.
///
/// Fails, before changing anything, if a constrained row has no diagonal
/// entry (Elmer would fall back to a list matrix there) or dimensions
/// disagree.
pub fn apply_dirichlet(
    a: &mut CsrMatrix<f64>,
    b: &mut [f64],
    dofs: &[(u32, f64)],
    mode: DirichletMode,
) -> Result<(), LinalgError> {
    let n = a.nrows;
    if b.len() != n || a.row_ptr.len() != n + 1 {
        return Err(LinalgError::Dim(format!(
            "matrix has {n} rows, rhs {} entries, row_ptr {}",
            b.len(),
            a.row_ptr.len()
        )));
    }
    let mut set = DirichletSet::new(n);
    for &(row, v) in dofs {
        let row = row as usize;
        if row >= n {
            return Err(LinalgError::Dim(format!("Dirichlet row {row} >= {n}")));
        }
        set.set(row, v);
    }
    if !set.constrained.iter().any(|&c| c) {
        return Ok(());
    }
    let diag_pos = diagonal_positions(a);
    if let Some(k) = (0..n).find(|&k| set.constrained[k] && diag_pos[k].is_none()) {
        return Err(LinalgError::Unsupported(format!(
            "Dirichlet row {k} has no diagonal entry"
        )));
    }

    let scaling = match mode.scaling {
        DirichletScaling::None => None,
        DirichletScaling::Diagonal => Some(diagonal_scaling(a, &diag_pos)),
    };

    if mode.symmetric {
        for r in 0..n {
            if set.constrained[r] {
                continue;
            }
            for k in a.row_ptr[r]..a.row_ptr[r + 1] {
                let c = a.col_idx[k] as usize;
                if set.constrained[c] {
                    b[r] -= a.values[k] * set.values[c];
                    a.values[k] = 0.0;
                }
            }
        }
    }

    for k in 0..n {
        if !set.constrained[k] {
            continue;
        }
        let s = match &scaling {
            None => 1.0,
            Some(d) => {
                let mut s = d[k];
                if s.abs() <= f64::MIN_POSITIVE {
                    s = 1.0;
                }
                1.0 / (s * s)
            }
        };
        a.values[a.row_ptr[k]..a.row_ptr[k + 1]].fill(0.0);
        a.values[diag_pos[k].expect("checked above")] = s;
        b[k] = s * set.values[k];
    }
    Ok(())
}

fn diagonal_positions(a: &CsrMatrix<f64>) -> Vec<Option<usize>> {
    match &a.diag {
        Some(d) => d.iter().map(|&p| (p != usize::MAX).then_some(p)).collect(),
        None => (0..a.nrows)
            .map(|i| {
                let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
                a.col_idx[lo..hi]
                    .iter()
                    .position(|&c| c as usize == i)
                    .map(|p| lo + p)
            })
            .collect(),
    }
}

/// `ScaleLinearSystemDiagonal` with `ApplyScaling=.FALSE.`
/// (MatrixScaling.F90:122): `1/sqrt|a_ii|`, with the absolute row sum
/// standing in for a (near) zero diagonal, and 1 for an all-zero row.
fn diagonal_scaling(a: &CsrMatrix<f64>, diag_pos: &[Option<usize>]) -> Vec<f64> {
    let tiny = f64::MIN_POSITIVE;
    let mut d: Vec<f64> = diag_pos
        .iter()
        .map(|p| p.map_or(0.0, |p| a.values[p]))
        .collect();
    if d.iter().any(|x| x.abs() <= tiny) {
        for (i, di) in d.iter_mut().enumerate() {
            if di.abs() <= tiny {
                let row = &a.values[a.row_ptr[i]..a.row_ptr[i + 1]];
                *di = row.iter().fold(0.0, |acc, v| acc + v.abs());
            }
        }
    }
    for di in &mut d {
        *di = if di.abs() > tiny {
            1.0 / di.abs().sqrt()
        } else {
            1.0
        };
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 3x3 full pattern with the given row-major values.
    fn full3(v: [f64; 9]) -> CsrMatrix<f64> {
        CsrMatrix {
            nrows: 3,
            ncols: 3,
            row_ptr: vec![0, 3, 6, 9],
            col_idx: vec![0, 1, 2, 0, 1, 2, 0, 1, 2],
            values: v.to_vec(),
            diag: Some(vec![0, 4, 8]),
        }
    }

    #[test]
    fn unsymmetric_unscaled_only_touches_constrained_rows() {
        let mut a = full3([4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]);
        let mut b = vec![1.0, 2.0, 3.0];
        let mode = DirichletMode {
            symmetric: false,
            scaling: DirichletScaling::None,
        };
        apply_dirichlet(&mut a, &mut b, &[(0, 5.0)], mode).unwrap();
        assert_eq!(
            a.values,
            vec![1.0, 0.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]
        );
        assert_eq!(b, vec![5.0, 2.0, 3.0]);
    }

    #[test]
    fn symmetric_scaled_moves_columns_to_rhs() {
        let mut a = full3([4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]);
        let mut b = vec![1.0, 2.0, 3.0];
        let mode = DirichletMode {
            symmetric: true,
            scaling: DirichletScaling::Diagonal,
        };
        apply_dirichlet(&mut a, &mut b, &[(2, 1.0), (0, 9.0), (0, 2.0)], mode).unwrap();
        let s = 1.0 / ((1.0 / 4f64.sqrt()) * (1.0 / 4f64.sqrt()));
        assert_eq!(a.values, vec![s, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, s]);
        assert_eq!(b, vec![2.0 * s, 2.0 + 2.0 + 1.0, s]);
    }

    #[test]
    fn zero_diagonal_falls_back_to_row_sum() {
        let mut a = full3([0.0, -1.0, 3.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]);
        let mut b = vec![0.0; 3];
        apply_dirichlet(&mut a, &mut b, &[(0, 1.0)], DirichletMode::default()).unwrap();
        let d = 1.0 / 4f64.sqrt();
        assert_eq!(a.values[0], 1.0 / (d * d));
        assert_eq!(b[0], 1.0 / (d * d));
    }

    #[test]
    fn missing_diagonal_is_an_error() {
        let mut a = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 1, 2],
            col_idx: vec![1, 1],
            values: vec![1.0, 1.0],
            diag: Some(vec![usize::MAX, 1]),
        };
        let mut b = vec![0.0; 2];
        let mode = DirichletMode {
            symmetric: false,
            scaling: DirichletScaling::None,
        };
        assert!(apply_dirichlet(&mut a, &mut b, &[(0, 1.0)], mode).is_err());
    }
}
