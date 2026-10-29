//! Permutation and matrix pattern of a p-element solver: `CreateMatrix`
//! (ElementUtils.F90:1963) over the slot space of [`crate::dofs`].
//!
//! 1. The size of `Perm` (`Ndeg`, MainUtils.F90:1657);
//! 2. `CreateMatrix` recomputes the mesh edge/bubble maxima;
//! 3. `InitialPermutation` (Lists.F90:234) numbers nodal, edge and bubble
//!    slots of the active elements in element order, or the identity when
//!    every slot is used;
//! 4. `MakeListMatrix` (ElementUtils.F90:381) couples all slots of each
//!    active element (`PackSortAdd`);
//! 5. `OptimizeBandwidth` and `InitializeMatrix` as in `elmer-assembly`.
//!
//! Not covered (callers must refuse them): DG, radiation, projectors,
//! `Apply Conforming BCs`, `Optimize Bandwidth` on boundary-only solvers
//! (`BCMode`) and parallel halo pruning.

use elmer_assembly::{
    bandwidth, bandwidth_reordered, initialize_matrix, rcm_order, BandwidthReport, ListMatrix,
    Perm, PermOptions,
};
use elmer_core::{CsrMatrix, NO_DOF};

use crate::defs::{ElementDefs, Kind};
use crate::dofs::DofLayout;
use crate::pmesh::PMesh;

/// Permutation, pattern and bandwidth decision of one p solver.
#[derive(Clone, Debug)]
pub struct PSystem {
    /// Slot to Perm index; `perm.perm.len()` is `Ndeg`.
    pub perm: Perm,
    /// Zero-valued pattern, `perm.n_dofs * dofs` square.
    pub matrix: CsrMatrix<f64>,
    pub bandwidth: BandwidthReport,
}

/// `Ndeg`: the size of `Perm` for a non-DG solver.
pub fn perm_size(pm: &PMesh, global_bubbles: bool) -> usize {
    let bulk = &pm.elements[..pm.n_bulk];
    let max_bdofs = bulk.iter().map(|e| e.bdofs).max().unwrap_or(0);
    let dofs_per_node = bulk
        .iter()
        .map(|e| e.ndofs / e.nodes.len())
        .max()
        .unwrap_or(0);
    let mut max_edofs = pm.edges.iter().map(|e| e.bdofs).max().unwrap_or(0);
    let mut max_fdofs = pm.faces.iter().map(|f| f.bdofs).max().unwrap_or(0);
    let mut ndeg = dofs_per_node * pm.n_nodes();
    if global_bubbles {
        ndeg += max_bdofs * pm.n_bulk;
        for el in &pm.elements[pm.n_bulk..] {
            if el.code >= 300 {
                max_fdofs = max_fdofs.max(el.bdofs);
            } else {
                max_edofs = max_edofs.max(el.bdofs);
            }
        }
    }
    ndeg + max_edofs * pm.edges.len() + max_fdofs * pm.faces.len()
}

/// The mesh maxima as `CreateMatrix` resets them before numbering.
pub fn update_mesh_maxima(pm: &mut PMesh, global_bubbles: bool) {
    let mut min_e = usize::MAX;
    let mut max_e = 0;
    let mut min_f = usize::MAX;
    let mut max_f = 0;
    for e in &pm.edges {
        min_e = min_e.min(e.bdofs);
        max_e = max_e.max(e.bdofs);
    }
    for f in &pm.faces {
        min_f = min_f.min(f.bdofs);
        max_f = max_f.max(f.bdofs);
    }
    if global_bubbles {
        for el in &pm.elements[pm.n_bulk..] {
            if el.code >= 300 {
                min_f = min_f.min(el.bdofs);
                max_f = max_f.max(el.bdofs);
            } else {
                min_e = min_e.min(el.bdofs);
                max_e = max_e.max(el.bdofs);
            }
        }
    }
    pm.max_edge_dofs = max_e;
    pm.min_edge_dofs = if min_e <= max_e { min_e } else { max_e };
    pm.max_face_dofs = max_f;
    pm.min_face_dofs = if min_f <= max_f { min_f } else { max_f };
    pm.max_bdofs = pm.elements[..pm.n_bulk]
        .iter()
        .map(|e| e.bdofs)
        .max()
        .unwrap_or(0);
}

/// `InitialPermutation` (non-DG, equation mode). `active(e)` is
/// `CheckElementEquation` over bulk and boundary elements (the elements
/// must have a body, as Elmer indexes `Def_Dofs` with it). Returns the
/// forward table of length [`perm_size`] (0-based, [`NO_DOF`] unused) and
/// the number of slots used.
pub fn initial_permutation(
    lay: &DofLayout<'_>,
    active: &dyn Fn(usize) -> bool,
) -> (Vec<u32>, usize) {
    let pm = lay.mesh;
    let defs = lay.defs;
    let gb = lay.global_bubbles;
    let mut perm = vec![NO_DOF; perm_size(pm, gb)];
    let mut k = 0u32;
    let mut number = |perm: &mut [u32], j: usize| {
        if perm[j] == NO_DOF {
            perm[j] = k;
            k += 1;
        }
    };
    let n_all = pm.elements.len();
    let body = |e: usize| pm.elements[e].body.max(1) as usize;

    let mut edge_dofs = vec![0usize; pm.edges.len()];
    let mut face_dofs = vec![0usize; pm.faces.len()];
    if defs.any_p() {
        for e in (0..n_all).filter(|&e| active(e)) {
            let el = &pm.elements[e];
            let p = defs.get(el.family(), body(e), Kind::P);
            if el.code >= 300 {
                for &g in &el.edges {
                    let d = &mut edge_dofs[g as usize];
                    *d = (*d).max(pm.edge_dofs_of(e, p));
                }
            }
            if el.code >= 500 {
                for (i, &g) in el.faces.iter().enumerate() {
                    let d = &mut face_dofs[g as usize];
                    *d = (*d).max(pm.face_dofs_of(e, p, i + 1));
                }
            }
        }
    }
    let face_base = lay.nodal_offset() + lay.edge_offset();

    let bubble_base = lay.nodal_offset() + lay.edge_offset() + lay.face_offset();
    for e in (0..n_all).filter(|&e| active(e)) {
        let el = &pm.elements[e];
        let row = defs.row(el.family(), body(e));
        let ndofs = row[Kind::Nodal as usize - 1];
        if ndofs > 0 {
            for &node in &el.nodes {
                for j in 0..ndofs as usize {
                    number(&mut perm, pm.max_ndofs * node as usize + j);
                }
            }
        }
        for &g in &el.edges {
            if el.code == 202 && !gb {
                continue;
            }
            let e_def = row[Kind::Edge as usize - 1];
            let n = if e_def >= 0 {
                e_def as usize
            } else if row[Kind::P as usize - 1] > 1 {
                edge_dofs[g as usize]
            } else {
                0
            };
            for i in 1..=n {
                number(&mut perm, lay.edge_slot(g, i) as usize);
            }
        }
        for &g in &el.faces {
            let face = &pm.faces[g as usize];
            if el.code == face.code && !gb {
                continue;
            }
            let mut l = row[Kind::Face as usize - 1].max(0);
            let fam = usize::from(face.code / 100) + 6;
            if l == 0 {
                if let Some(left) = face.left {
                    l = defs.get(fam, body(left as usize), Kind::Bubble).max(0);
                }
                if let Some(right) = face.right {
                    l = l.max(defs.get(fam, body(right as usize), Kind::Bubble));
                }
            }
            let n = if l > 0 {
                l as usize
            } else if row[Kind::P as usize - 1] > 1 {
                face_dofs[g as usize]
            } else {
                0
            };
            for i in 0..n {
                number(&mut perm, face_base + pm.max_face_dofs * g as usize + i);
            }
        }
        if gb && el.bubble_base > 0 {
            let nb = row[Kind::Bubble as usize - 1];
            let p = row[Kind::P as usize - 1];
            let n = if nb >= 0 || p >= 1 {
                let b = if p > 1 { pm.bubble_dofs_of(e, p) } else { 0 };
                b.max(nb.max(0) as usize)
            } else if lay.bubbles {
                el.bdofs
            } else {
                0
            };
            for i in 1..=n {
                number(&mut perm, bubble_base + el.bubble_index(i) - 1);
            }
        }
    }

    // Connection elements: their node numbers index Perm directly.
    for el in &pm.elements[pm.n_bulk..] {
        if el.code == 102 {
            for &node in &el.nodes {
                number(&mut perm, node as usize);
            }
        }
    }

    let k = k as usize;
    if k == perm.len() {
        for (i, p) in perm.iter_mut().enumerate() {
            *p = i as u32;
        }
    }
    (perm, k)
}

/// `MakeListMatrix` (non-DG): the slot graph of the active elements in
/// the numbering of `perm`.
pub fn list_matrix(lay: &DofLayout<'_>, active: &dyn Fn(usize) -> bool, perm: &Perm) -> ListMatrix {
    let pm = lay.mesh;
    let gb = lay.global_bubbles;
    let edofs = if pm.edges.is_empty() {
        0
    } else {
        pm.max_edge_dofs
    };
    let nodal_offset = lay.nodal_offset();
    let bubble_base = nodal_offset + lay.edge_offset() + lay.face_offset();
    let fdofs = if pm.faces.is_empty() {
        0
    } else {
        pm.max_face_dofs
    };
    let face_base = nodal_offset + edofs * pm.edges.len();
    let mut list = ListMatrix::new(perm.n_dofs);
    let mut ind: Vec<u32> = Vec::new();
    let at = |s: usize| perm.perm[s];
    for e in (0..pm.elements.len()).filter(|&e| active(e)) {
        let el = &pm.elements[e];
        ind.clear();
        let per_node = el.ndofs / el.nodes.len();
        for &node in &el.nodes {
            for j in 0..per_node {
                ind.push(at(pm.max_ndofs * node as usize + j));
            }
        }
        if edofs > 0 && !el.edges.is_empty() {
            for &g in &el.edges {
                for i in 0..pm.edges[g as usize].bdofs {
                    ind.push(at(nodal_offset + edofs * g as usize + i));
                }
            }
            if gb && el.code == 202 {
                let g = el.edges[0] as usize;
                for i in 0..el.bdofs {
                    ind.push(at(nodal_offset + edofs * g + i));
                }
            }
        }
        if fdofs > 0 && !el.faces.is_empty() {
            for &g in &el.faces {
                for i in 0..pm.faces[g as usize].bdofs {
                    ind.push(at(face_base + fdofs * g as usize + i));
                }
            }
            if gb && el.code == pm.faces[el.faces[0] as usize].code {
                let g = el.faces[0] as usize;
                for i in 0..el.bdofs {
                    ind.push(at(face_base + fdofs * g + i));
                }
            }
        }
        if gb && el.bubble_base > 0 {
            for i in 1..=el.bdofs {
                ind.push(at(bubble_base + el.bubble_index(i) - 1));
            }
        }
        ind.retain(|&p| p != NO_DOF);
        list.add_clique(&mut ind);
    }
    for el in &pm.elements[pm.n_bulk..] {
        if !(102..200).contains(&el.code) {
            continue;
        }
        let k1 = perm.perm[el.nodes[0] as usize];
        if k1 != NO_DOF {
            for &node in &el.nodes {
                let k2 = perm.perm[node as usize];
                if k2 != NO_DOF {
                    list.add(k1 as usize, k2);
                    list.add(k2 as usize, k1);
                }
            }
        }
        if el.code == 102 {
            let k2 = perm.perm[el.nodes[1] as usize];
            if k2 != NO_DOF {
                list.add(k2 as usize, k2);
            }
        }
    }
    list
}

/// `CreateMatrix` for a p solver with `dofs` components: resets the mesh
/// maxima, numbers the slots, builds the pattern and applies Elmer's
/// bandwidth decision.
pub fn create_p_system(
    pm: &mut PMesh,
    defs: &ElementDefs,
    global_bubbles: bool,
    bubbles: bool,
    active: &dyn Fn(usize) -> bool,
    opts: PermOptions,
    dofs: usize,
) -> PSystem {
    update_mesh_maxima(pm, global_bubbles);
    let lay = DofLayout {
        mesh: pm,
        defs,
        global_bubbles,
        bubbles,
    };
    let (perm, _) = initial_permutation(&lay, active);
    let init = Perm::from_perm(perm);
    let list = list_matrix(&lay, active, &init);

    // OptimizeBandwidth, as elmer-assembly's `create_system`.
    let initial = bandwidth(&list) + 1;
    let (perm, report, map) = if !opts.optimize_bandwidth || init.n_dofs == 0 {
        let r = BandwidthReport {
            initial,
            optimized: None,
            accepted: false,
        };
        (init, r, None)
    } else {
        let new_index = rcm_order(&list);
        let optimized = bandwidth_reordered(&list, &new_index) + 1;
        let accepted = !(initial < optimized && !opts.use_optimized_always);
        let r = BandwidthReport {
            initial,
            optimized: Some(optimized),
            accepted,
        };
        if accepted {
            let p = init
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
            (Perm::from_perm(p), r, Some(new_index))
        } else {
            (init, r, None)
        }
    };
    let matrix = initialize_matrix(&list, map.as_deref(), dofs);
    PSystem {
        perm,
        matrix,
        bandwidth: report,
    }
}
