//! DOF permutation and matrix topology: Elmer's `CreateMatrix`
//! (ElementUtils.F90:1963) for nodal solvers.
//!
//! The pipeline is the same as Elmer's:
//! 1. `InitialPermutation` (Lists.F90:234) numbers the nodes of the active
//!    elements in element order, or uses the identity when every node is
//!    active;
//! 2. `MakeListMatrix` (ElementUtils.F90:381) builds the node graph in that
//!    numbering;
//! 3. `OptimizeBandwidth` (BandwidthOptimize.F90:182) runs Elmer's RCM and
//!    keeps it unless it made the bandwidth worse;
//! 4. `InitializeMatrix` (ElementUtils.F90:1857) expands nodes to `dofs`
//!    rows each and sorts the columns (`CRS_SortMatrix`).
//!
//! Supported: nodal Lagrange DOFs, one slot per node, on bulk and boundary
//! elements the caller marks active, plus Elmer's `102` connection
//! elements. Not supported (callers must refuse such cases): edge, face
//! and bubble DOFs, DG numbering, radiation and projector (periodic,
//! mortar) couplings, `Apply Conforming BCs`, and parallel halo pruning.

use elmer_core::{CsrMatrix, ElemId, MeshView, NO_DOF};

use crate::bandwidth::{bandwidth, bandwidth_reordered, rcm_order, BandwidthReport};
use crate::list::ListMatrix;

/// Node to DOF-slot permutation (Elmer `Perm`, made 0-based).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Perm {
    /// `perm[node]` is the node's slot, or [`NO_DOF`] when it is not in the
    /// system (Elmer `Perm(node) == 0`). Length `mesh.n_nodes()`.
    pub perm: Vec<u32>,
    /// `inv[slot]` is the node in that slot. Length `n_dofs`.
    pub inv: Vec<u32>,
    /// Number of slots; the matrix has `n_dofs * dofs` rows.
    pub n_dofs: usize,
}

impl Perm {
    /// Builds `inv` and `n_dofs` from a forward table.
    ///
    /// # Panics
    /// If the used slots are not exactly `0..n` for some `n`.
    pub fn from_perm(perm: Vec<u32>) -> Self {
        let n_dofs = perm.iter().filter(|&&p| p != NO_DOF).count();
        let mut inv = vec![NO_DOF; n_dofs];
        for (node, &p) in perm.iter().enumerate() {
            if p != NO_DOF {
                assert!((p as usize) < n_dofs, "perm slot {p} out of range");
                assert!(inv[p as usize] == NO_DOF, "perm slot {p} used twice");
                inv[p as usize] = node as u32;
            }
        }
        Self { perm, inv, n_dofs }
    }

    /// Matrix row of component `dof` (0-based) of `node`, with `dofs`
    /// components per node, or `None` when the node is not in the system.
    pub fn row(&self, node: u32, dof: usize, dofs: usize) -> Option<usize> {
        match self.perm[node as usize] {
            NO_DOF => None,
            p => Some(dofs * p as usize + dof),
        }
    }
}

/// Bandwidth keywords of a solver section.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PermOptions {
    /// `Optimize Bandwidth`, default true (MainUtils.F90:1724). Forced off
    /// by Elmer for `DG Reduced Basis` and `Apply Conforming BCs`; the
    /// caller applies those overrides.
    pub optimize_bandwidth: bool,
    /// `Optimize Bandwidth Use Always`, default false: keep the RCM order
    /// even when it is worse than the initial one.
    pub use_optimized_always: bool,
}

impl Default for PermOptions {
    fn default() -> Self {
        Self {
            optimize_bandwidth: true,
            use_optimized_always: false,
        }
    }
}

/// Permutation, matrix pattern and bandwidth decision of one solver.
#[derive(Clone, Debug)]
pub struct System {
    pub perm: Perm,
    /// Pattern with zero values, `perm.n_dofs * dofs` square.
    pub matrix: CsrMatrix<f64>,
    pub bandwidth: BandwidthReport,
}

fn all_elements(mesh: &dyn MeshView) -> impl Iterator<Item = ElemId> {
    (0..(mesh.n_bulk() + mesh.n_boundary()) as u32).map(ElemId)
}

fn boundary_elements(mesh: &dyn MeshView) -> impl Iterator<Item = ElemId> {
    let nb = mesh.n_bulk() as u32;
    (nb..nb + mesh.n_boundary() as u32).map(ElemId)
}

/// `InitialPermutation` for nodal DOFs. `active(e)` must be Elmer's
/// `CheckElementEquation` for the solver's equation, over bulk *and*
/// boundary elements (a boundary element counts only through a BC's
/// `Body Id`; otherwise its body id is 0 and it is inactive).
pub fn initial_permutation(mesh: &dyn MeshView, active: &dyn Fn(ElemId) -> bool) -> Perm {
    let n_nodes = mesh.n_nodes();
    let mut perm = vec![NO_DOF; n_nodes];
    let mut k = 0u32;
    let mut number = |perm: &mut [u32], nodes: &[u32]| {
        for &node in nodes {
            let p = &mut perm[node as usize];
            if *p == NO_DOF {
                *p = k;
                k += 1;
            }
        }
    };
    for e in all_elements(mesh) {
        if active(e) {
            number(&mut perm, mesh.element(e).nodes);
        }
    }
    // Connection elements join the system whatever their equation.
    for e in boundary_elements(mesh) {
        let el = mesh.element(e);
        if el.etype.0 == 102 {
            number(&mut perm, el.nodes);
        }
    }
    // CreateMatrix: when every slot is used Elmer discards the element
    // order and takes the identity.
    if k as usize == n_nodes {
        for (i, p) in perm.iter_mut().enumerate() {
            *p = i as u32;
        }
    }
    Perm::from_perm(perm)
}

/// `MakeListMatrix` for nodal DOFs: the node graph of the active elements
/// in the numbering of `perm`.
pub fn list_matrix(
    mesh: &dyn MeshView,
    active: &dyn Fn(ElemId) -> bool,
    perm: &Perm,
) -> ListMatrix {
    let mut list = ListMatrix::new(perm.n_dofs);
    let mut ind = Vec::new();
    for e in all_elements(mesh) {
        if !active(e) {
            continue;
        }
        ind.clear();
        ind.extend(
            mesh.element(e)
                .nodes
                .iter()
                .map(|&n| perm.perm[n as usize])
                .filter(|&p| p != NO_DOF),
        );
        list.add_clique(&mut ind);
    }
    // Point-family boundary elements couple their first node to the rest.
    for e in boundary_elements(mesh) {
        let el = mesh.element(e);
        if !(102..200).contains(&el.etype.0) {
            continue;
        }
        let k1 = perm.perm[el.nodes[0] as usize];
        if k1 != NO_DOF {
            for &node in el.nodes {
                let k2 = perm.perm[node as usize];
                if k2 != NO_DOF {
                    list.add(k1 as usize, k2);
                    list.add(k2 as usize, k1);
                }
            }
        }
        if el.etype.0 == 102 {
            let k2 = perm.perm[el.nodes[1] as usize];
            if k2 != NO_DOF {
                list.add(k2 as usize, k2);
            }
        }
    }
    list
}

/// The permutation a solver ends up with: initial numbering, then RCM when
/// `opt_bw` (Elmer's `Optimize Bandwidth`) is set and it does not make the
/// bandwidth worse.
pub fn build_perm(mesh: &dyn MeshView, active: &dyn Fn(ElemId) -> bool, opt_bw: bool) -> Perm {
    let opts = PermOptions {
        optimize_bandwidth: opt_bw,
        ..PermOptions::default()
    };
    build_perm_with(mesh, active, opts).0
}

/// [`build_perm`] with all options, also returning the bandwidth decision.
pub fn build_perm_with(
    mesh: &dyn MeshView,
    active: &dyn Fn(ElemId) -> bool,
    opts: PermOptions,
) -> (Perm, BandwidthReport) {
    let init = initial_permutation(mesh, active);
    let list = list_matrix(mesh, active, &init);
    let (perm, report, _) = optimize(init, &list, opts);
    (perm, report)
}

/// [`build_perm_with`] with extra list-matrix rows: each `(node, nodes)`
/// couples `node`'s row to every one of `nodes`, as `MakeListMatrix` adds
/// the diffuse gray radiation couplings (ElementUtils.F90:744-830) before
/// `OptimizeBandwidth`.
pub fn build_perm_coupled(
    mesh: &dyn MeshView,
    active: &dyn Fn(ElemId) -> bool,
    opts: PermOptions,
    extra: &[(u32, Vec<u32>)],
) -> (Perm, BandwidthReport) {
    let init = initial_permutation(mesh, active);
    let mut list = list_matrix(mesh, active, &init);
    for (node, cols) in extra {
        let k1 = init.perm[*node as usize];
        if k1 == NO_DOF {
            continue;
        }
        for &c in cols {
            let k2 = init.perm[c as usize];
            if k2 != NO_DOF {
                list.add(k1 as usize, k2);
            }
        }
    }
    let (perm, report, _) = optimize(init, &list, opts);
    (perm, report)
}

/// Returns the final permutation, the report, and the initial-to-final
/// slot map when the order changed.
fn optimize(
    init: Perm,
    list: &ListMatrix,
    opts: PermOptions,
) -> (Perm, BandwidthReport, Option<Vec<u32>>) {
    let initial = bandwidth(list) + 1;
    if !opts.optimize_bandwidth || init.n_dofs == 0 {
        let report = BandwidthReport {
            initial,
            optimized: None,
            accepted: false,
        };
        return (init, report, None);
    }
    let new_index = rcm_order(list);
    let optimized = bandwidth_reordered(list, &new_index) + 1;
    let accepted = !(initial < optimized && !opts.use_optimized_always);
    let report = BandwidthReport {
        initial,
        optimized: Some(optimized),
        accepted,
    };
    if !accepted {
        return (init, report, None);
    }
    let perm = init
        .perm
        .iter()
        .map(|&p| {
            if p == NO_DOF {
                NO_DOF
            } else {
                new_index[p as usize]
            }
        })
        .collect();
    (Perm::from_perm(perm), report, Some(new_index))
}

/// Matrix pattern for `dofs` components per node, rows and columns in the
/// numbering of `perm` (which must come from the same `active` set).
pub fn create_matrix(
    mesh: &dyn MeshView,
    active: &dyn Fn(ElemId) -> bool,
    perm: &Perm,
    dofs: usize,
) -> CsrMatrix<f64> {
    initialize_matrix(&list_matrix(mesh, active, perm), None, dofs)
}

/// `CreateMatrix` for a nodal solver: permutation, bandwidth decision and
/// matrix pattern in one pass over the mesh.
pub fn create_system(
    mesh: &dyn MeshView,
    active: &dyn Fn(ElemId) -> bool,
    opts: PermOptions,
    dofs: usize,
) -> System {
    let init = initial_permutation(mesh, active);
    let list = list_matrix(mesh, active, &init);
    let (perm, bandwidth, map) = optimize(init, &list, opts);
    let matrix = initialize_matrix(&list, map.as_deref(), dofs);
    System {
        perm,
        matrix,
        bandwidth,
    }
}

/// `InitializeMatrix` + `CRS_SortMatrix`: expand list rows (renumbered by
/// `map` when given) to `dofs` rows each, columns ascending.
///
/// `diag[i]` is the position of `(i, i)`, or `usize::MAX` when the row has
/// no diagonal entry (Elmer's `Diag(i) == 0`).
pub fn initialize_matrix(list: &ListMatrix, map: Option<&[u32]>, dofs: usize) -> CsrMatrix<f64> {
    assert!(dofs > 0, "dofs must be positive");
    let k = list.nrows();
    let n = k * dofs;
    let slot = |i: usize| map.map_or(i, |m| m[i] as usize);

    let mut row_len = vec![0usize; n];
    for i in 0..k {
        let j = slot(i);
        for l in 0..dofs {
            row_len[dofs * j + l] = dofs * list.degree(i);
        }
    }
    let mut row_ptr = Vec::with_capacity(n + 1);
    row_ptr.push(0);
    for len in &row_len {
        row_ptr.push(row_ptr.last().unwrap() + len);
    }
    let nnz = row_ptr[n];

    let mut col_idx = vec![0u32; nnz];
    for i in 0..k {
        let j = slot(i);
        for l in 0..dofs {
            let start = row_ptr[dofs * j + l];
            let cols = &mut col_idx[start..start + dofs * list.degree(i)];
            let mut c = 0;
            for &kk in list.row(i) {
                let base = dofs * slot(kk as usize);
                for m in 0..dofs {
                    cols[c] = (base + m) as u32;
                    c += 1;
                }
            }
            if map.is_some() {
                cols.sort_unstable();
            }
        }
    }

    let diag = (0..n)
        .map(|i| {
            let row = &col_idx[row_ptr[i]..row_ptr[i + 1]];
            row.binary_search(&(i as u32))
                .map_or(usize::MAX, |p| row_ptr[i] + p)
        })
        .collect();

    CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr,
        col_idx,
        values: vec![0.0; nnz],
        diag: Some(diag),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{ElementRef, ElementType};

    /// Bulk elements given as node lists of one type, plus boundary ones.
    struct Mesh {
        n_nodes: usize,
        bulk: Vec<(u16, Vec<u32>)>,
        boundary: Vec<(u16, Vec<u32>)>,
    }

    impl MeshView for Mesh {
        fn dim(&self) -> usize {
            1
        }
        fn n_nodes(&self) -> usize {
            self.n_nodes
        }
        fn coord(&self, node: u32) -> [f64; 3] {
            [node as f64, 0.0, 0.0]
        }
        fn n_bulk(&self) -> usize {
            self.bulk.len()
        }
        fn n_boundary(&self) -> usize {
            self.boundary.len()
        }
        fn element(&self, e: ElemId) -> ElementRef<'_> {
            let i = e.0 as usize;
            let (code, nodes) = if i < self.bulk.len() {
                &self.bulk[i]
            } else {
                &self.boundary[i - self.bulk.len()]
            };
            ElementRef {
                etype: ElementType(*code),
                nodes,
                tag: 1,
                parents: [None, None],
            }
        }
    }

    fn lines(edges: &[(u32, u32)], n_nodes: usize) -> Mesh {
        Mesh {
            n_nodes,
            bulk: edges.iter().map(|&(a, b)| (202, vec![a, b])).collect(),
            boundary: Vec::new(),
        }
    }

    /// Elmer's RCM makes this graph worse (half bandwidth 3 -> 4 with the
    /// `+1`); traced by hand through `Levelize`/`Renumber`.
    fn worse_graph() -> Mesh {
        lines(&[(2, 4), (2, 0), (3, 2), (1, 0), (2, 1)], 5)
    }

    fn bulk_only(m: &Mesh) -> impl Fn(ElemId) -> bool + '_ {
        move |e: ElemId| (e.0 as usize) < m.bulk.len()
    }

    #[test]
    fn worse_bandwidth_is_rejected() {
        let m = worse_graph();
        let (perm, report) = build_perm_with(&m, &bulk_only(&m), PermOptions::default());
        assert_eq!(
            report,
            BandwidthReport {
                initial: 3,
                optimized: Some(4),
                accepted: false
            }
        );
        assert_eq!(report.half_bandwidth(), 3);
        assert_eq!(perm.perm, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn use_always_keeps_the_worse_order() {
        let m = worse_graph();
        let opts = PermOptions {
            optimize_bandwidth: true,
            use_optimized_always: true,
        };
        let sys = create_system(&m, &bulk_only(&m), opts, 1);
        assert!(sys.bandwidth.accepted);
        assert_eq!(sys.perm.perm, vec![2, 1, 3, 4, 0]);
        assert_eq!(sys.perm.inv, vec![4, 1, 0, 2, 3]);
        // Pattern rows come out sorted in the new numbering.
        let row = |r: usize| &sys.matrix.col_idx[sys.matrix.row_ptr[r]..sys.matrix.row_ptr[r + 1]];
        assert_eq!(row(3), &[0, 1, 2, 3, 4]);
        assert_eq!(row(0), &[0, 3]);
    }

    #[test]
    fn no_optimisation_reports_initial_only() {
        let m = worse_graph();
        let opts = PermOptions {
            optimize_bandwidth: false,
            use_optimized_always: false,
        };
        let (perm, report) = build_perm_with(&m, &bulk_only(&m), opts);
        assert_eq!(report.optimized, None);
        assert!(!report.accepted);
        assert_eq!(perm.perm, vec![0, 1, 2, 3, 4]);
        assert_eq!(build_perm(&m, &bulk_only(&m), false), perm);
    }

    #[test]
    fn inactive_nodes_are_left_out_in_element_order() {
        // Elements 0 and 2 active; node 0 and 1 belong only to element 1.
        let m = lines(&[(4, 3), (0, 1), (3, 2)], 5);
        let active = |e: ElemId| e.0 != 1;
        let init = initial_permutation(&m, &active);
        assert_eq!(init.perm, vec![NO_DOF, NO_DOF, 2, 1, 0]);
        assert_eq!(init.inv, vec![4, 3, 2]);
        let a = create_matrix(&m, &active, &init, 2);
        assert_eq!(a.nrows, 6);
        // Slot 1 (node 3) touches slots 0, 1, 2: 3 nodes x 2 components.
        assert_eq!(a.row_ptr[3] - a.row_ptr[2], 6);
        assert_eq!(a.diag.as_ref().unwrap()[3], a.row_ptr[3] + 3);
    }

    #[test]
    fn connection_elements_join_and_couple() {
        // Node 3 sits only on a 102 boundary element linking it to node 0.
        let m = Mesh {
            n_nodes: 4,
            bulk: vec![(202, vec![0, 1]), (202, vec![1, 2])],
            boundary: vec![(102, vec![0, 3]), (101, vec![2])],
        };
        let active = |e: ElemId| e.0 < 2;
        let init = initial_permutation(&m, &active);
        // Every node numbered: identity.
        assert_eq!(init.perm, vec![0, 1, 2, 3]);
        let list = list_matrix(&m, &active, &init);
        assert_eq!(list.row(0), &[0, 1, 3]);
        assert_eq!(list.row(3), &[0, 3]);
        assert_eq!(list.row(2), &[1, 2]);
    }

    #[test]
    fn empty_system() {
        let m = lines(&[(0, 1)], 2);
        let sys = create_system(&m, &|_| false, PermOptions::default(), 1);
        assert_eq!(sys.perm.n_dofs, 0);
        assert_eq!(sys.matrix.nrows, 0);
        assert_eq!(sys.matrix.row_ptr, vec![0]);
    }
}
