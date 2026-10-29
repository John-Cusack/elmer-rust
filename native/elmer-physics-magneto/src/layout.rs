//! The DOF layout of the AV solvers: Elmer's `Perm` without bandwidth
//! optimisation for `Element = "n:1 e:1"` (or `"n:0 e:1"`): the nodal
//! slots of the active elements in node order, then one slot per edge in
//! first-use order over the active elements.

use elmer_core::{CsrMatrix, Scalar};

use crate::case::Case;
use crate::error::MagnetoError;
use crate::topology::Edges;

pub struct AvLayout {
    /// Active bulk elements (indices into `mesh.bulk`).
    pub active: Vec<usize>,
    /// Slot of each mesh node, `u32::MAX` when not in the system.
    pub node_slot: Vec<u32>,
    pub n_nodal: usize,
    pub edges: Edges,
    /// DOFs of each active element: nodal slots (with `V`), then edges.
    pub elem_dofs: Vec<Vec<u32>>,
    /// Whether the nodal scalar potential is part of the system.
    pub with_v: bool,
}

impl AvLayout {
    pub fn build(case: &Case, solver_id: usize, with_v: bool) -> Result<Self, MagnetoError> {
        let active: Vec<usize> = (0..case.mesh.bulk.len())
            .filter(|&e| case.solver_active_in(solver_id, case.mesh.bulk[e].body))
            .collect();
        if active.is_empty() {
            return Err(MagnetoError::Input("no active elements".into()));
        }
        let mut node_slot = vec![u32::MAX; case.mesh.nodes.len()];
        let mut edges = Edges::default();
        let mut elem_edges = Vec::with_capacity(active.len());
        for &e in &active {
            let el = &case.mesh.bulk[e];
            if !matches!(el.code, 504 | 706 | 808) {
                return Err(MagnetoError::unsupported(format!("bulk element type {}", el.code)));
            }
            if with_v {
                for &i in &el.nodes {
                    node_slot[i as usize] = 0;
                }
            }
            elem_edges.push(edges.add_element(el.code, &el.nodes));
        }
        let mut n_nodal = 0u32;
        for s in node_slot.iter_mut().filter(|s| **s == 0) {
            *s = n_nodal;
            n_nodal += 1;
        }
        let mut layout = Self {
            active,
            node_slot,
            n_nodal: n_nodal as usize,
            edges,
            elem_dofs: Vec::new(),
            with_v,
        };
        layout.elem_dofs = layout
            .active
            .iter()
            .zip(&elem_edges)
            .map(|(&e, eids)| layout.dofs_of(&case.mesh.bulk[e].nodes, eids))
            .collect();
        Ok(layout)
    }

    /// Number of slots.
    pub fn n(&self) -> usize {
        self.n_nodal + self.edges.n
    }

    /// The slots of an element with nodes `nodes` and edges `eids`.
    pub fn dofs_of(&self, nodes: &[u32], eids: &[u32]) -> Vec<u32> {
        let mut d = Vec::with_capacity(nodes.len() + eids.len());
        if self.with_v {
            d.extend(nodes.iter().map(|&i| self.node_slot[i as usize]));
        }
        d.extend(eids.iter().map(|&e| self.n_nodal as u32 + e));
        d
    }

    /// The zero matrix on the element pattern.
    pub fn matrix<T: Scalar>(&self) -> CsrMatrix<T> {
        pattern(self.n(), &self.elem_dofs)
    }
}

/// The zero `n x n` matrix on the pattern of element DOF lists: sorted
/// columns per row, with the diagonal.
pub fn pattern<T: Scalar>(n: usize, elem_dofs: &[Vec<u32>]) -> CsrMatrix<T> {
    let mut rows: Vec<Vec<u32>> = (0..n).map(|i| vec![i as u32]).collect();
    for dofs in elem_dofs {
        for &r in dofs {
            rows[r as usize].extend_from_slice(dofs);
        }
    }
    let mut row_ptr = vec![0];
    let mut col_idx = Vec::new();
    let mut diag = Vec::with_capacity(n);
    for (i, r) in rows.iter_mut().enumerate() {
        r.sort_unstable();
        r.dedup();
        diag.push(col_idx.len() + r.iter().position(|&c| c as usize == i).expect("diagonal"));
        col_idx.extend_from_slice(r);
        row_ptr.push(col_idx.len());
    }
    CsrMatrix {
        nrows: n,
        ncols: n,
        values: vec![T::zero(); col_idx.len()],
        row_ptr,
        col_idx,
        diag: Some(diag),
    }
}

/// A real zero matrix on the element pattern and a zero right-hand side.
pub fn nodal_matrix(n: usize, elem_dofs: &[Vec<u32>]) -> (CsrMatrix<f64>, Vec<f64>) {
    (pattern(n, elem_dofs), vec![0.0; n])
}

/// Adds a local `nd x nd` row-major matrix and vector at `dofs`.
pub fn glue<T: Scalar>(a: &mut CsrMatrix<T>, rhs: &mut [T], dofs: &[u32], stiff: &[T], force: &[T]) {
    let nd = dofs.len();
    for (p, &r) in dofs.iter().enumerate() {
        let r = r as usize;
        rhs[r] = rhs[r] + force[p];
        let (lo, hi) = (a.row_ptr[r], a.row_ptr[r + 1]);
        for (q, &c) in dofs.iter().enumerate() {
            let k = lo + a.col_idx[lo..hi].binary_search(&c).expect("pattern");
            a.values[k] = a.values[k] + stiff[p * nd + q];
        }
    }
}

/// The nodal coordinates and global (1-based) node numbers of an element.
pub fn geometry(case: &Case, nodes: &[u32]) -> (Vec<[f64; 3]>, Vec<i64>) {
    let xyz = nodes.iter().map(|&i| case.mesh.nodes[i as usize]).collect();
    let gindex = nodes.iter().map(|&i| i as i64 + 1).collect();
    (xyz, gindex)
}
