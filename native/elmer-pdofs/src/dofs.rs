//! Global DOF slots of one element: `mGetElementDOFs` (ElemInfo.F90:9375)
//! and `mGetBoundaryIndexesFromParent` (ElementUtils.F90:4229).
//!
//! A slot is an index into the solver's `Perm` table (Elmer's `Indexes`,
//! made 0-based). The slot space is laid out as
//!
//! * nodes: `MaxNDOFs * node + j`,
//! * edges: `NodalOffset + MaxEdgeDOFs * edge + i`,
//! * faces: `+ EdgeOffset` (none in 1D/2D),
//! * bubbles: `+ FaceOffset + BubbleIndexes(i) - 1`,
//!
//! with `NodalOffset = MaxNDOFs * NumberOfNodes` and so on, all from the
//! [`PMesh`] maxima.

use elmer_core::NO_DOF;

use crate::defs::{ElementDefs, Kind};
use crate::pmesh::PMesh;

/// One solver's view of the p mesh: its `Def_Dofs` and DOF keywords.
#[derive(Clone, Copy, Debug)]
pub struct DofLayout<'a> {
    pub mesh: &'a PMesh,
    pub defs: &'a ElementDefs,
    /// `Solver % GlobalBubbles` (`Bubbles in Global System`, which
    /// defaults to true when the `Element` keyword has `p:` or `b:`).
    pub global_bubbles: bool,
    /// `Bubbles` (solver keyword, default false): bubble count when the
    /// definition gives neither `b:` nor `p:`.
    pub bubbles: bool,
}

impl<'a> DofLayout<'a> {
    pub fn new(mesh: &'a PMesh, defs: &'a ElementDefs, global_bubbles: bool) -> Self {
        Self {
            mesh,
            defs,
            global_bubbles,
            bubbles: false,
        }
    }

    /// `MaxNDOFs * NumberOfNodes`.
    pub fn nodal_offset(&self) -> usize {
        self.mesh.max_ndofs * self.mesh.n_nodes()
    }

    /// `MaxEdgeDOFs * NumberOfEdges`.
    pub fn edge_offset(&self) -> usize {
        self.mesh.max_edge_dofs * self.mesh.edges.len()
    }

    /// `MaxFaceDOFs * NumberOfFaces`.
    pub fn face_offset(&self) -> usize {
        self.mesh.max_face_dofs * self.mesh.n_faces
    }

    /// Slot of edge DOF `i` (1-based) of global edge `edge`.
    pub fn edge_slot(&self, edge: u32, i: usize) -> u32 {
        (self.nodal_offset() + self.mesh.max_edge_dofs * edge as usize + i - 1) as u32
    }

    /// `isActivePElement`: the element has p definitions and some body
    /// gives its family a degree.
    pub fn is_active_p(&self, e: usize) -> bool {
        let el = &self.mesh.elements[e];
        el.pdefs.is_some() && self.defs.p_active(el.family())
    }

    /// The body whose `Def_Dofs` an element uses in `mGetElementDOFs`:
    /// its own, else the left parent's, else (when that is 0 or out of
    /// range) the right parent's, else 1.
    pub fn dof_body(&self, e: usize) -> usize {
        let el = &self.mesh.elements[e];
        let nb = self.defs.n_bodies() as i32;
        let body_of = |p: Option<u32>| p.map(|p| self.mesh.elements[p as usize].body);
        let mut id = el.body;
        if id == 0 && e >= self.mesh.n_bulk {
            if let Some(b) = body_of(el.parents[0]) {
                id = b;
            }
            if id == 0 || id > nb {
                if let Some(b) = body_of(el.parents[1]) {
                    id = b;
                }
            }
        }
        if id == 0 {
            1
        } else {
            id as usize
        }
    }

    /// Bubble count of a bulk element: `MAX(b:, getBubbleDOFs(p))`, or all
    /// of `BubbleIndexes` with `Bubbles` when neither is defined.
    fn bubble_count(&self, e: usize, fam: usize, id: usize) -> usize {
        let nb = self.defs.get(fam, id, Kind::Bubble);
        let p = self.defs.get(fam, id, Kind::P);
        if nb >= 0 || p >= 1 {
            let b = if p > 1 {
                self.mesh.bubble_dofs_of(e, p)
            } else {
                0
            };
            b.max(nb.max(0) as usize)
        } else if self.bubbles {
            self.mesh.elements[e].bdofs
        } else {
            0
        }
    }

    /// `mGetElementDOFs` for a non-DG solver: the element's slots, nodal
    /// first, then edges, then the parent's edge (boundary line elements)
    /// or the bubbles (bulk elements with global bubbles).
    pub fn element_dofs(&self, e: usize, out: &mut Vec<u32>) {
        out.clear();
        let pm = self.mesh;
        let el = &pm.elements[e];
        let fam = el.family();
        let id = self.dof_body(e);
        let defs = self.defs;

        let ndofs = defs.get(fam, id, Kind::Nodal);
        if ndofs > 0 {
            for &node in &el.nodes {
                for j in 0..ndofs as usize {
                    out.push((pm.max_ndofs * node as usize + j) as u32);
                }
            }
        }

        let row = defs.row(fam, id);
        let mut need_edges = row[1..].iter().any(|&v| v >= 0);
        if !need_edges && (fam == 3 || fam == 4) {
            need_edges = defs.get(6 + fam, id, Kind::Bubble) >= 0;
        }
        if !need_edges {
            return;
        }

        let edge_count = |elem: usize, p_fam: usize| -> usize {
            let e_def = defs.get(fam, id, Kind::Edge);
            if e_def >= 0 {
                e_def as usize
            } else if defs.get(fam, id, Kind::P) > 1 {
                pm.edge_dofs_of(elem, defs.get(p_fam, id, Kind::P))
            } else {
                0
            }
        };

        let boundary = e >= pm.n_bulk;
        let mut edges_done = false;
        if !el.edges.is_empty() {
            edges_done = true;
            for &edge in &el.edges {
                // Every mesh edge is a 202 element.
                if el.code == 202 && !(self.global_bubbles && el.body > 0 && boundary) {
                    edges_done = false;
                    continue;
                }
                for i in 1..=edge_count(e, fam) {
                    out.push(self.edge_slot(edge, i));
                }
            }
        }

        if boundary {
            let parent = if self.is_active_p(e) {
                el.pdefs.as_ref().and_then(|pd| pd.local_parent)
            } else {
                el.parents[0].or(el.parents[1])
            };
            let Some(parent) = parent else { return };
            let par = &pm.elements[parent as usize];
            if fam == 2 && !edges_done && !par.edges.is_empty() {
                let ind = if self.is_active_p(e) {
                    el.pdefs.as_ref().map_or(0, |pd| pd.local_number)
                } else {
                    // First parent edge whose nodes all belong to the element.
                    par.edges
                        .iter()
                        .position(|&g| {
                            let n = pm.edges[g as usize].nodes;
                            n.iter().filter(|x| el.nodes.contains(x)).count() == el.nodes.len()
                        })
                        .map_or(0, |k| k + 1)
                };
                if ind == 0 {
                    return;
                }
                let edge = par.edges[ind - 1];
                for i in 1..=edge_count(parent as usize, par.family()) {
                    out.push(self.edge_slot(edge, i));
                }
            }
        } else if el.bdofs > 0 && self.global_bubbles {
            let base = self.nodal_offset() + self.edge_offset() + self.face_offset();
            for i in 1..=self.bubble_count(e, fam, id) {
                out.push((base + el.bubble_index(i) - 1) as u32);
            }
        }
    }

    /// [`DofLayout::element_dofs`] mapped through `perm` (slot to Perm
    /// index, [`NO_DOF`] when absent): the `dof_index` that
    /// `elmer_assembly::CsrAssembler::add_local` expects.
    pub fn element_rows(&self, e: usize, perm: &[u32], out: &mut Vec<u32>) {
        self.element_dofs(e, out);
        for s in out.iter_mut() {
            *s = perm.get(*s as usize).copied().unwrap_or(NO_DOF);
        }
    }

    /// `mGetBoundaryIndexesFromParent` for a boundary element of a 1D or
    /// 2D mesh: its nodes, then the edge DOF slots of the parent's local
    /// edge (`Element % BDOFs` of them). Empty when the element has no p
    /// definitions or no local parent.
    ///
    /// Note the Fortran offsets the edge slots by `NumberOfNodes`, not by
    /// `MaxNDOFs * NumberOfNodes`; kept.
    pub fn boundary_indexes_from_parent(&self, e: usize, out: &mut Vec<u32>) {
        out.clear();
        let pm = self.mesh;
        let el = &pm.elements[e];
        let Some(pd) = &el.pdefs else { return };
        let Some(parent) = pd.local_parent else {
            return;
        };
        out.extend_from_slice(&el.nodes);
        let par = &pm.elements[parent as usize];
        if par.dim() == 2 && !pm.edges.is_empty() && pd.local_number > 0 {
            let edge = par.edges[pd.local_number - 1] as usize;
            for i in 1..=el.bdofs {
                out.push((pm.n_nodes() + edge * pm.max_edge_dofs + i - 1) as u32);
            }
        }
    }
}
