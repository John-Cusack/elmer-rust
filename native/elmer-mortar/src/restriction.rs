//! The saddle-point system of a constrained solve
//! (`SolveWithLinearRestriction`, SolveCore.F90): the stiffness matrix with
//! the constraint matrix and its transpose around it,
//!
//! ```text
//! [ A  C^T ] [ x ]   [ b ]
//! [ C  0   ] [ l ] = [ g ]
//! ```
//!
//! solved as one system by the solver's linear solver, `l` being the
//! Lagrange multipliers.

use std::collections::BTreeMap;

use elmer_core::{CsrMatrix, Params, NO_DOF};

use crate::constraint::ConstraintMatrix;
use crate::MortarError;

/// The collection system (`StiffMatrix % CollectionMatrix`) and its start
/// vector.
#[derive(Clone, Debug)]
pub struct CollectionSystem {
    pub matrix: CsrMatrix<f64>,
    pub rhs: Vec<f64>,
    /// The current solution, then the multipliers (zero, or the exported
    /// ones from the previous solve).
    pub x: Vec<f64>,
    /// The size of the stiffness matrix.
    pub n: usize,
}

impl CollectionSystem {
    /// The solution of the stiffness matrix's unknowns.
    pub fn solution(&self) -> &[f64] {
        &self.x[..self.n]
    }

    /// The Lagrange multipliers (`Export Lagrange Multiplier`).
    pub fn multipliers(&self) -> &[f64] {
        &self.x[self.n..]
    }
}

/// The inputs of [`collection_system`].
pub struct Restriction<'a> {
    /// The solver section.
    pub params: &'a dyn Params,
    /// The stiffness matrix after the Dirichlet conditions, and its
    /// right-hand side and current solution.
    pub a: &'a CsrMatrix<f64>,
    pub b: &'a [f64],
    pub x: &'a [f64],
    pub c: &'a ConstraintMatrix,
    /// The Dirichlet rows `(row, value)` (`ConstrainedDOF`, `DValues`).
    pub constrained: &'a [(u32, f64)],
    /// The exported multipliers of the previous solve, if any.
    pub multipliers: Option<&'a [f64]>,
    /// Parent node -> node slot (`Solver % Variable % Perm`).
    pub perm: &'a [u32],
    pub dofs: usize,
    /// `Mesh % NumberOfNodes`.
    pub n_nodes: usize,
}

const REFUSED: [&str; 7] = [
    "restriction system residual mode",
    "no explicit constrained matrix",
    "eliminate linear constraints",
    "constraint scaling",
    "linear system complex",
    "eigen analysis",
    "harmonic analysis",
];

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `AddToMatrixElement` on a list matrix.
fn add(list: &mut [BTreeMap<u32, f64>], i: usize, j: usize, v: f64) {
    list[i].entry(j as u32).and_modify(|x| *x += v).or_insert(v);
}

/// The collection system of `SolveWithLinearRestriction` without an
/// `AddMatrix`, elimination, residual mode, constraint scaling or complex
/// arithmetic. With `Enforce Exact Dirichlet BCs` (default on) a constraint
/// does not touch Dirichlet rows: their columns move to the constraint's
/// right-hand side, and a constraint row of a Dirichlet node becomes
/// `l = 0`. A constraint row with no entries on the unknowns gets a zero
/// right-hand side (and a unit diagonal with `Use Tree Gauge`).
pub fn collection_system(r: &Restriction<'_>) -> Result<CollectionSystem, MortarError> {
    let p = r.params;
    if let Some(k) = REFUSED.iter().find(|k| flag(p, k)) {
        return Err(MortarError::Unsupported(format!(
            "SolveWithLinearRestriction: {k}"
        )));
    }
    if p.get("lagrange multiplier relaxation factor").is_some() {
        return Err(MortarError::Unsupported(
            "Lagrange Multiplier Relaxation Factor".into(),
        ));
    }
    let n = r.a.nrows;
    let m = r.c.n_rows();
    let total = n + m;
    let mut list: Vec<BTreeMap<u32, f64>> = vec![BTreeMap::new(); total];
    let mut rhs = vec![0.0; total];

    let enforce =
        p.logical("enforce exact dirichlet bcs").unwrap_or(true) && !r.constrained.is_empty();
    let mut dvalue: Vec<Option<f64>> = vec![None; n];
    for &(row, v) in r.constrained {
        dvalue[row as usize] = Some(v);
    }
    let tree_gauge = flag(p, "use tree gauge");
    let mut iperm = vec![NO_DOF; n.div_ceil(r.dofs.max(1))];
    for (node, &s) in r.perm.iter().enumerate() {
        if s != NO_DOF {
            iperm[s as usize] = node as u32;
        }
    }

    for i in (0..m).rev() {
        let k = n + i;
        add(&mut list, k, k, 0.0);
        let range = r.c.rows[i]..r.c.rows[i + 1];
        let mut rowsum = 0.0;
        for j in range.clone().rev() {
            rowsum += r.c.values[j].abs();
        }
        let mut non_empty = false;
        if rowsum > f64::EPSILON {
            for j in range.rev() {
                let colj = r.c.cols[j] as usize;
                let val = r.c.values[j];
                if val.abs() < f64::EPSILON * rowsum {
                    continue;
                }
                let fixed = if enforce && colj < n {
                    dvalue[colj]
                } else {
                    None
                };
                match fixed {
                    None => {
                        add(&mut list, colj, k, val);
                        if colj < n {
                            add(&mut list, k, colj, val);
                            non_empty |= val != 0.0;
                        }
                    }
                    Some(d) => rhs[k] -= val * d,
                }
            }
        }
        let mut free = true;
        if enforce {
            let l = r.c.inv_perm[i] % n;
            if dvalue[l].is_some() {
                let node = iperm[l / r.dofs];
                if node == NO_DOF || (node as usize) < r.n_nodes {
                    free = false;
                    for v in list[k].values_mut() {
                        *v = 0.0;
                    }
                    rhs[k] = 0.0;
                    list[k].insert(k as u32, 1.0);
                }
            }
        }
        if free {
            if !non_empty {
                rhs[k] = 0.0;
                if tree_gauge {
                    list[k].insert(k as u32, 1.0);
                }
            } else {
                rhs[k] += r.c.rhs[i];
            }
        }
    }

    for i in (0..n).rev() {
        for j in (r.a.row_ptr[i]..r.a.row_ptr[i + 1]).rev() {
            add(&mut list, i, r.a.col_idx[j] as usize, r.a.values[j]);
        }
        rhs[i] += r.b[i];
    }

    let mut matrix = CsrMatrix {
        nrows: total,
        ncols: total,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    for row in list {
        for (c, v) in row {
            matrix.col_idx.push(c);
            matrix.values.push(v);
        }
        matrix.row_ptr.push(matrix.col_idx.len());
    }
    matrix.diag = Some(elmer_sparse::ops::find_diag(&matrix));

    let mut x = vec![0.0; total];
    x[..n].copy_from_slice(&r.x[..n]);
    if flag(p, "export lagrange multiplier") {
        if let Some(l) = r.multipliers {
            let k = l.len().min(m);
            x[n..n + k].copy_from_slice(&l[..k]);
        }
    }
    Ok(CollectionSystem { matrix, rhs, x, n })
}
