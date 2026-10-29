//! Mesh-level p-element setup: what `PrepareMesh` (MeshLoad.F90:643) adds to
//! a loaded mesh before any solver builds its matrix.
//!
//! * `NonNodalElements` (MeshLoad.F90:713): per-element nodal and bubble DOF
//!   counts and the p definitions (`Element % PDefs`), including the
//!   MATC-dependent degree `p:%f` evaluated element by element;
//! * `FindMeshEdges2D` (MeshBasics.F90:3875): the global edge table, edges
//!   numbered in first-visit order over bulk then boundary elements;
//! * `SetMeshEdgeFaceDOFs` (MeshBasics.F90:2123) and `AssignLocalNumber`
//!   (MeshBasics.F90:5826): edge degrees and DOF counts, and the local edge
//!   (`PDefs % LocalNumber`) and parent of every p boundary element, whose
//!   nodes are **reordered** to the parent's edge orientation;
//! * `SetMeshMaxDOFs` (MeshBasics.F90:2271): integration point counts
//!   (`getNumberOfGaussPoints`, PElementMaps.F90:1135), the mesh maxima and
//!   the bubble index blocks.
//!
//! Supported: meshes of first-order 1D and 2D elements (`202`, `303`,
//! `404`; boundary `101`, `202`). 3D meshes (face tables, `ConvertToACTetra`)
//! are refused with [`Error::Unsupported`]. Element definitions given in
//! equation sections are not handled (the caller must refuse them); solver
//! sections and body-wise `Solver i: Element` keywords are.

use elmer_core::{ElemId, MeshView};
use elmer_fe_hp::maps::{bubble_dofs, bubble_max_degree_1d, edge_map, face_dofs};

use crate::defs::{is_mesh_dependent, ElementDefs, MaxDefs};
use crate::Error;

/// `Element % PDefs`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PDefs {
    /// `P`: for bulk elements the mesh-wide maximum degree seen so far (see
    /// [`PMesh::new`]); for edges and boundary elements the maximum degree of
    /// the adjacent p elements.
    pub p: i32,
    pub tetra_type: i32,
    /// `isEdge`: the element is an edge or a boundary element.
    pub is_edge: bool,
    /// `GaussPoints`: integration points for `GaussPoints(Element)`.
    pub gauss_points: usize,
    /// `localNumber`: 1-based local edge of [`PDefs::local_parent`], 0 if
    /// none was found.
    pub local_number: usize,
    /// `LocalParent`: bulk element index.
    pub local_parent: Option<u32>,
    pub serendipity: bool,
}

impl PDefs {
    /// `AllocatePDefinitions` (MeshIO.F90:88).
    fn new(serendipity: bool) -> Self {
        Self {
            p: 0,
            tetra_type: 0,
            is_edge: false,
            gauss_points: 0,
            local_number: 0,
            local_parent: None,
            serendipity,
        }
    }
}

/// One mesh element (bulk or boundary) with its p-element data.
#[derive(Clone, Debug, PartialEq)]
pub struct PElement {
    pub code: u16,
    /// Node indices, 0-based. For p boundary elements these are in the
    /// orientation of the parent's local edge (`AssignLocalNumber`), which
    /// may differ from the mesh file.
    pub nodes: Vec<u32>,
    /// `BodyId`: the mesh tag for bulk elements, 0 for boundary elements
    /// unless set with [`PMeshOptions::boundary_body`].
    pub body: i32,
    /// Boundary tag (boundary elements), 0 for bulk elements.
    pub boundary_tag: i32,
    /// `BoundaryInfo % Left / Right`, as bulk element indices.
    pub parents: [Option<u32>; 2],
    /// `EdgeIndexes`, 0-based global edges; empty when not associated.
    pub edges: Vec<u32>,
    pub pdefs: Option<PDefs>,
    /// `NDOFs`: nodal DOFs of the element (nodes times DOFs per node).
    pub ndofs: usize,
    /// `BDOFs`: bubble DOFs (bulk), or the edge DOF count copied from the
    /// parent's edge (boundary).
    pub bdofs: usize,
    /// `DGDOFs`.
    pub dg_dofs: usize,
    /// First of the `bdofs` consecutive 1-based `BubbleIndexes`; 0 when the
    /// element has none.
    pub bubble_base: usize,
}

impl PElement {
    pub fn family(&self) -> usize {
        usize::from(self.code / 100)
    }

    /// `Element % TYPE % DIMENSION`.
    pub fn dim(&self) -> usize {
        match self.code / 100 {
            1 => 0,
            2 => 1,
            3 | 4 => 2,
            _ => 3,
        }
    }

    /// `Element % TYPE % NumberOfEdges` of the first-order types.
    pub fn n_edges(&self) -> usize {
        n_edges(self.code)
    }

    /// `BubbleIndexes(i)`, `i` 1-based.
    pub fn bubble_index(&self, i: usize) -> usize {
        self.bubble_base + i - 1
    }
}

/// A global edge (`Mesh % Edges(i)`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MeshEdge {
    /// End nodes, in the orientation of the left parent's local edge.
    pub nodes: [u32; 2],
    /// `BoundaryInfo % Left / Right`: the first two bulk elements visiting
    /// the edge.
    pub left: Option<u32>,
    pub right: Option<u32>,
    pub pdefs: Option<PDefs>,
    /// `BDOFs`: edge DOF count (the maximum over the adjacent elements).
    pub bdofs: usize,
}

fn n_edges(code: u16) -> usize {
    match code / 100 {
        2 => 1,
        3 => 3,
        4 => 4,
        5 => 6,
        6 => 8,
        7 => 9,
        8 => 12,
        _ => 0,
    }
}

fn first_order(code: u16) -> bool {
    matches!(code, 101 | 202 | 303 | 404 | 504 | 605 | 706 | 808)
}

/// Local (1-based) end nodes of edge `k` (1-based) of a 1D or 2D element,
/// as `FindMeshEdges2D` walks them.
fn edge_nodes_2d(code: u16, k: usize) -> [usize; 2] {
    let n = n_edges(code);
    if n == 1 {
        [1, 2]
    } else if k < n {
        [k, k + 1]
    } else {
        [k, 1]
    }
}

/// What a solver contributes to the mesh setup.
pub struct SolverElements<'a> {
    /// `Solver % Def_Dofs`, from [`ElementDefs::from_solver_element`];
    /// updated in place by body-wise and MATC definitions.
    pub defs: &'a mut ElementDefs,
    /// The solver-section `Element` keyword.
    pub element: Option<&'a str>,
    /// `Solver % DG` (`Discontinuous Galerkin`).
    pub dg: bool,
    /// `Solver i: Element` in a body section (1-based body), when present.
    pub body_element: &'a dyn Fn(usize) -> Option<String>,
}

/// A MATC degree function `p_of(function, [body, x, y, z])`.
pub type MatcDegree<'a> = dyn FnMut(&str, [f64; 4]) -> Result<f64, Error> + 'a;

/// Mesh-level keywords and model facts that [`PMesh::new`] needs.
pub struct PMeshOptions<'a> {
    /// `Serendipity P Elements` (Simulation), default true.
    pub serendipity: bool,
    /// `Need Edges 2D` / `Need Edges 3D` (Simulation) for the mesh
    /// dimension, when given.
    pub need_edges: Option<bool>,
    /// Number of bodies in the model.
    pub n_bodies: usize,
    /// Whether body `b` (1-based) has an `Equation` keyword.
    pub body_has_equation: &'a dyn Fn(usize) -> bool,
    /// MATC degree function: `p_of(function, [body, x, y, z])` with the
    /// element midpoint, as `GetMatcReal(function//'(cx)', 4, ...)`.
    pub matc_degree: &'a mut MatcDegree<'a>,
    /// `BodyId` of boundary element `i` (0-based among boundary elements);
    /// Elmer leaves it 0 unless the BC gives a `Body Id`.
    pub boundary_body: &'a dyn Fn(usize) -> i32,
}

/// A mesh prepared for p elements.
#[derive(Clone, Debug, PartialEq)]
pub struct PMesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub n_bulk: usize,
    /// Bulk then boundary elements, as `Mesh % Elements`.
    pub elements: Vec<PElement>,
    pub edges: Vec<MeshEdge>,
    /// `Mesh % NumberOfFaces` (always 0: 3D is not supported).
    pub n_faces: usize,
    pub max_ndofs: usize,
    pub max_edge_dofs: usize,
    pub min_edge_dofs: usize,
    pub max_face_dofs: usize,
    pub max_bdofs: usize,
    pub max_element_dofs: usize,
    pub max_element_nodes: usize,
}

impl PMesh {
    pub fn n_nodes(&self) -> usize {
        self.coords.len()
    }

    pub fn n_boundary(&self) -> usize {
        self.elements.len() - self.n_bulk
    }

    pub fn element(&self, e: ElemId) -> &PElement {
        &self.elements[e.0 as usize]
    }

    pub fn is_boundary(&self, e: ElemId) -> bool {
        e.0 as usize >= self.n_bulk
    }

    /// `PrepareMesh` p-element steps for the solvers on this mesh.
    ///
    /// `max_defs` is the `Def_Dofs(10, 6)` that `LoadModel` accumulated
    /// with [`ElementDefs::from_solver_element`] (start from
    /// [`crate::defs::initial_max_defs`]); it is raised further here.
    ///
    /// Elmer's quirk is kept: when the degree comes from a MATC function,
    /// each bulk element's `PDefs % P` is the maximum degree over the
    /// elements processed **so far**, so it depends on the element order;
    /// the solver's per-body degree (used by the basis) is the value of the
    /// last element of that body.
    pub fn new(
        mesh: &dyn MeshView,
        solvers: &mut [SolverElements<'_>],
        max_defs: &mut MaxDefs,
        opts: PMeshOptions<'_>,
    ) -> Result<Self, Error> {
        let n_bulk = mesh.n_bulk();
        let n_all = n_bulk + mesh.n_boundary();
        let mut elements = Vec::with_capacity(n_all);
        for i in 0..n_all {
            let el = mesh.element(ElemId(i as u32));
            let boundary = i >= n_bulk;
            elements.push(PElement {
                code: el.etype.0,
                nodes: el.nodes.to_vec(),
                body: if boundary {
                    (opts.boundary_body)(i - n_bulk)
                } else {
                    el.tag
                },
                boundary_tag: if boundary { el.tag } else { 0 },
                parents: el.parents.map(|p| p.map(|e| e.0)),
                edges: Vec::new(),
                pdefs: None,
                ndofs: 0,
                bdofs: 0,
                dg_dofs: 0,
                bubble_base: 0,
            });
        }
        let coords = (0..mesh.n_nodes() as u32).map(|n| mesh.coord(n)).collect();
        let dim = elements[..n_bulk]
            .iter()
            .map(PElement::dim)
            .max()
            .unwrap_or(0);
        let mut pm = PMesh {
            dim,
            coords,
            n_bulk,
            elements,
            edges: Vec::new(),
            n_faces: 0,
            max_ndofs: 0,
            max_edge_dofs: 0,
            min_edge_dofs: usize::MAX,
            max_face_dofs: 0,
            max_bdofs: 0,
            max_element_dofs: 0,
            max_element_nodes: 0,
        };
        pm.non_nodal_elements(solvers, max_defs, opts)?;
        Ok(pm)
    }

    fn midpoint(&self, e: usize) -> [f64; 3] {
        let nodes = &self.elements[e].nodes;
        let n = nodes.len() as f64;
        let mut c = [0.0; 3];
        for (k, ck) in c.iter_mut().enumerate() {
            // SUM(Mesh % Nodes % x(Element % NodeIndexes)) / n
            let mut s = 0.0;
            for &node in nodes {
                s += self.coords[node as usize][k];
            }
            *ck = s / n;
        }
        c
    }

    fn non_nodal_elements(
        &mut self,
        solvers: &mut [SolverElements<'_>],
        in_dofs: &mut MaxDefs,
        opts: PMeshOptions<'_>,
    ) -> Result<(), Error> {
        let n_bulk = self.n_bulk;
        let mut edge_dofs = vec![0usize; n_bulk];
        let mut face_dofs = vec![0usize; n_bulk];

        let first_order_elements = self.elements[..n_bulk].iter().all(|e| first_order(e.code));
        let mesh_deps = solvers
            .iter()
            .any(|s| s.element.is_some_and(is_mesh_dependent));
        // UpdateDefDofs: only MATC-dependent solver definitions are re-read
        // (no equation-section definitions are supported).
        let update: Vec<bool> = solvers
            .iter()
            .map(|s| s.element.is_some_and(is_mesh_dependent))
            .collect();

        let mut matc = |f: &str, cx: [f64; 4]| -> Result<i32, Error> {
            // NINT
            Ok((opts.matc_degree)(f, cx)?.round() as i32)
        };

        if !mesh_deps {
            for body in 1..=opts.n_bodies {
                if !(opts.body_has_equation)(body) {
                    continue;
                }
                for s in solvers.iter_mut() {
                    if let Some(def) = (s.body_element)(body) {
                        let dg = s.dg;
                        s.defs.get_max_defs(&def, body, dg, in_dofs, &mut |_| {
                            Err(Error::Unsupported(
                                "MATC degree in a body-wise element definition".into(),
                            ))
                        })?;
                    }
                }
            }
        }

        let mut need_edges = in_dofs[8][4] > 0 || in_dofs[9][4] > 0;
        let mut body0 = i32::MIN;
        let mut found_eq = false;
        for i in 0..n_bulk {
            let body = self.elements[i].body;
            if mesh_deps {
                if body != body0 {
                    found_eq = body >= 1 && (opts.body_has_equation)(body as usize);
                }
                if found_eq {
                    let mid = self.midpoint(i);
                    for (s, &upd) in solvers.iter_mut().zip(&update) {
                        if !upd {
                            continue;
                        }
                        let def = s.element.unwrap_or("");
                        let dg = s.dg;
                        let cx = [f64::from(body), mid[0], mid[1], mid[2]];
                        s.defs
                            .get_max_defs(def, body as usize, dg, in_dofs, &mut |f| matc(f, cx))?;
                    }
                }
                body0 = body;
            }

            let el = &mut self.elements[i];
            let n = el.nodes.len();
            let fam = el.family() - 1;
            el.ndofs = n * in_dofs[fam][0].max(0) as usize;
            edge_dofs[i] = in_dofs[fam][1].max(0) as usize;
            face_dofs[i] = in_dofs[fam][2].max(0) as usize;

            let dg = in_dofs[fam][3] == 0;
            if dg {
                in_dofs[fam][3] = n as i32;
            }
            if in_dofs[fam][3] > 0 {
                el.dg_dofs = in_dofs[fam][3] as usize;
            }
            if dg {
                need_edges = true;
            } else {
                need_edges = need_edges || in_dofs[fam][1] > 0 || in_dofs[fam][2] > 0;
            }

            if first_order_elements && in_dofs[fam][5] > 0 {
                let mut pd = PDefs::new(opts.serendipity);
                if !dg {
                    need_edges = in_dofs[fam][5] > 1;
                }
                pd.p = in_dofs[fam][5];
                el.bdofs = if in_dofs[fam][4] > 0 {
                    in_dofs[fam][4] as usize
                } else {
                    bubble_dofs(el.code, pd.p, pd.serendipity)
                };
                if el.code == 504 {
                    return Err(Error::Unsupported(
                        "tetrahedral p elements (ConvertToACTetra)".into(),
                    ));
                }
                el.pdefs = Some(pd);
            } else {
                el.pdefs = None;
                el.bdofs = in_dofs[fam][4].max(0) as usize;
            }
            self.max_element_nodes = self.max_element_nodes.max(n);
        }

        for i in n_bulk..self.elements.len() {
            let n = self.elements[i].nodes.len();
            let fam = self.elements[i].family();
            let mut ndofs = n * in_dofs[fam - 1][0].max(0) as usize;
            let mut bdofs = 0usize;
            let mut pdefs = None;
            let parents = self.elements[i].parents;
            let edim = self.elements[i].dim();
            for parent in parents.into_iter().flatten() {
                let j = parent as usize;
                if j >= n_bulk {
                    // "Element indexes beyond face or edge table": only warned.
                    continue;
                }
                if self.elements[j].ndofs == 0 {
                    ndofs = 0;
                }
                if edim == 1 {
                    bdofs = edge_dofs[j];
                } else {
                    bdofs = face_dofs[j];
                    if !need_edges && in_dofs[fam - 1][4] > 0 {
                        let mut pd = PDefs::new(opts.serendipity);
                        pd.is_edge = true;
                        pdefs = Some(pd);
                        bdofs = bdofs.max(in_dofs[fam - 1][4].max(0) as usize);
                    } else if fam + 6 <= 10 {
                        bdofs = bdofs.max(in_dofs[fam + 5][4].max(0) as usize);
                    }
                }
            }
            let el = &mut self.elements[i];
            el.ndofs = ndofs;
            el.bdofs = bdofs;
            if pdefs.is_some() {
                el.pdefs = pdefs;
            }
        }

        if let Some(ne) = opts.need_edges {
            need_edges = ne;
        }
        if need_edges {
            self.set_mesh_edge_face_dofs(&edge_dofs, &face_dofs)?;
        }
        if self.max_element_dofs == 0 {
            self.max_element_dofs = self.max_element_nodes;
        }
        self.set_mesh_max_dofs();

        // Edges are released again when nothing needs them (PrepareMesh,
        // "Why the heck did we allocate the edges and faces?!").
        let del = self.max_face_dofs + self.max_edge_dofs == 0
            && self.max_element_dofs == self.max_element_nodes
            && !solvers.iter().any(|s| s.dg);
        if del && !self.edges.is_empty() {
            self.edges.clear();
            for el in &mut self.elements {
                el.edges.clear();
            }
        }
        Ok(())
    }

    /// `FindMeshEdges2D` without a mask.
    fn find_mesh_edges_2d(&mut self) {
        let n_nodes = self.coords.len();
        let mut table: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n_nodes];
        let mut edges: Vec<MeshEdge> = Vec::new();
        for i in 0..self.elements.len() {
            let code = self.elements[i].code;
            if code < 200 {
                continue;
            }
            let ne = n_edges(code);
            let mut idx = vec![0u32; ne];
            for (k, slot) in idx.iter_mut().enumerate() {
                let [a, b] = edge_nodes_2d(code, k + 1);
                let (na, nb) = (self.elements[i].nodes[a - 1], self.elements[i].nodes[b - 1]);
                let (n1, n2) = if nb < na { (nb, na) } else { (na, nb) };
                // The hash list is searched head first; at most one entry
                // per (n1, n2), so the order does not matter.
                let found = table[n1 as usize]
                    .iter()
                    .find(|&&(node, _)| node == n2)
                    .map(|&(_, e)| e);
                let edge = match found {
                    Some(e) => e,
                    None => {
                        let e = edges.len() as u32;
                        table[n1 as usize].push((n2, e));
                        let pdefs = self.elements[i].pdefs.as_ref().map(|pd| {
                            let mut d = PDefs::new(pd.serendipity);
                            d.p = 0;
                            d
                        });
                        edges.push(MeshEdge {
                            nodes: [na, nb],
                            left: None,
                            right: None,
                            pdefs,
                            bdofs: 0,
                        });
                        e
                    }
                };
                *slot = edge;
                if i < self.n_bulk {
                    let ed = &mut edges[edge as usize];
                    if ed.left.is_some() {
                        ed.right = Some(i as u32);
                    } else {
                        ed.left = Some(i as u32);
                    }
                }
            }
            self.elements[i].edges = idx;
        }
        self.edges = edges;
    }

    /// `SetMeshEdgeFaceDOFs` with `EdgeDOFs`/`FaceDOFs`/`inDOFs` present and
    /// `NeedEdges` absent.
    fn set_mesh_edge_face_dofs(
        &mut self,
        edge_dofs: &[usize],
        _face_dofs: &[usize],
    ) -> Result<(), Error> {
        match self.dim {
            1 | 2 => self.find_mesh_edges_2d(),
            3 => {
                return Err(Error::Unsupported(
                    "3D meshes for p elements (FindMeshFaces3D / FindMeshEdges3D)".into(),
                ))
            }
            _ => {
                return Err(Error::Unsupported(
                    "mesh without 1D/2D bulk elements".into(),
                ))
            }
        }
        // FindMeshEdges2D also gives boundary line elements of a 2D mesh
        // their (single) edge index; AssignConstraints only copies BC tags.

        for i in 0..self.n_bulk {
            let el_p = self.elements[i].pdefs.as_ref().map(|pd| pd.p);
            let el_edges = self.elements[i].edges.clone();
            for &edge in &el_edges {
                let e = edge as usize;
                if let Some(p) = el_p {
                    let ser = self.elements[i].pdefs.as_ref().unwrap().serendipity;
                    let ed = &mut self.edges[e];
                    let pd = ed.pdefs.get_or_insert_with(|| PDefs::new(ser));
                    pd.p = p.max(pd.p);
                    ed.bdofs = ed.bdofs.max((pd.p - 1).max(0) as usize);
                    pd.is_edge = true;
                    // (BDOFs+2)**DIMENSION, an edge being 1D.
                    pd.gauss_points = ed.bdofs + 2;
                    let parent = ed.left.or(ed.right);
                    if let Some(parent) = parent {
                        self.assign_local_number_edge(e, parent as usize);
                    }
                }
                let ed = &mut self.edges[e];
                ed.bdofs = ed.bdofs.max(edge_dofs[i]);
                self.min_edge_dofs = self.min_edge_dofs.min(ed.bdofs);
                self.max_edge_dofs = self.max_edge_dofs.max(ed.bdofs);
            }
            if self.min_edge_dofs > self.max_edge_dofs {
                self.min_edge_dofs = self.max_edge_dofs;
            }
            // Faces: none in 1D/2D.
        }

        for i in self.n_bulk..self.elements.len() {
            let [left, right] = self.elements[i].parents;
            let mut p_alloc = false;
            if let Some(l) = left {
                if let Some(pd) = &self.elements[l as usize].pdefs {
                    p_alloc = true;
                    let mut d = PDefs::new(pd.serendipity);
                    d.is_edge = true;
                    self.elements[i].pdefs = Some(d);
                    self.assign_local_number_boundary(i, l as usize);
                }
            }
            if let Some(r) = right {
                if let Some(pd) = &self.elements[r as usize].pdefs {
                    if !p_alloc {
                        let mut d = PDefs::new(pd.serendipity);
                        d.is_edge = true;
                        self.elements[i].pdefs = Some(d);
                        self.assign_local_number_boundary(i, r as usize);
                    }
                }
            }
        }
        Ok(())
    }

    /// Local edge of bulk element `parent` (1-based) whose nodes all belong
    /// to `nodes`, with the parent's nodes along that edge.
    fn find_local_edge(&self, nodes: &[u32], parent: usize) -> Option<(usize, [u32; 2])> {
        let pel = &self.elements[parent];
        if pel.dim() != 2 {
            // AssignLocalNumber: nothing to do for 0D/1D parents; 3D refused.
            return None;
        }
        let map = edge_map(pel.code)?;
        for edge_number in 1..=pel.n_edges() {
            let edge = &self.edges[pel.edges[edge_number - 1] as usize];
            let n = edge.nodes.iter().filter(|g| nodes.contains(g)).count();
            if n == nodes.len() {
                let [a, b] = map[edge_number - 1];
                return Some((edge_number, [pel.nodes[a - 1], pel.nodes[b - 1]]));
            }
        }
        None
    }

    /// `AssignLocalNumber`'s guard before reordering: how many of the found
    /// entity's nodes appear among the first `n` parent nodes of the edge.
    fn reorder_hits(entity: &[u32; 2], bindex: &[u32; 2], n: usize) -> usize {
        entity
            .iter()
            .filter(|g| bindex[..n.min(2)].contains(g))
            .count()
    }

    /// `AssignLocalNumber(Edge, Parent)` for a mesh edge.
    fn assign_local_number_edge(&mut self, e: usize, parent: usize) {
        let nodes = self.edges[e].nodes;
        if let Some((num, bindex)) = self.find_local_edge(&nodes, parent) {
            let ed = &mut self.edges[e];
            if Self::reorder_hits(&nodes, &bindex, 2) == 2 {
                ed.nodes = bindex;
            }
            let pd = ed.pdefs.as_mut().expect("p edge");
            pd.local_number = num;
            pd.local_parent = Some(parent as u32);
            // Copying the edge's own attributes onto itself is a no-op.
        }
    }

    /// `AssignLocalNumber(Element, Parent)` for a boundary element.
    fn assign_local_number_boundary(&mut self, i: usize, parent: usize) {
        if self.elements[parent].dim() < 2 {
            return;
        }
        let nodes = self.elements[i].nodes.clone();
        let Some((num, bindex)) = self.find_local_edge(&nodes, parent) else {
            return;
        };
        let edge = &self.edges[self.elements[parent].edges[num - 1] as usize];
        let (is_edge, gp, p) = match &edge.pdefs {
            Some(pd) => (pd.is_edge, pd.gauss_points, pd.p),
            None => (false, 0, 0),
        };
        let edge_bdofs = edge.bdofs;
        let k = Self::reorder_hits(&edge.nodes, &bindex, nodes.len());
        let n = nodes.len().min(2);
        let el = &mut self.elements[i];
        if k == nodes.len() {
            el.nodes[..n].copy_from_slice(&bindex[..n]);
        }
        let pd = el.pdefs.as_mut().expect("p boundary element");
        pd.local_number = num;
        pd.local_parent = Some(parent as u32);
        pd.is_edge = is_edge;
        pd.gauss_points = gp;
        pd.p = p;
        el.bdofs = el.bdofs.max(edge_bdofs);
    }

    /// `getNumberOfGaussPoints(Element, Mesh)` for a bulk p element.
    fn number_of_gauss_points(&self, i: usize) -> usize {
        let el = &self.elements[i];
        let pd = el.pdefs.as_ref().expect("p element");
        let mut edge_p = 0;
        if el.dim() >= 2 {
            for &e in &el.edges {
                let ed = &self.edges[e as usize];
                if ed.bdofs == 0 {
                    continue;
                }
                edge_p = edge_p.max(ed.pdefs.as_ref().map_or(0, |d| d.p));
            }
        }
        let face_p = 0;
        let mut bubble_p = effective_bubble_p(el.code, pd.p, el.bdofs, pd.serendipity);
        if pd.serendipity && el.family() == 4 {
            let maxp = 1.max(edge_p).max(face_p).max(bubble_p);
            if !(pd.p < 4 && el.bdofs > 0) && maxp > 1 && maxp <= 7 {
                return match maxp {
                    2 => 8,
                    3 => 12,
                    4 => 20,
                    5 => 36,
                    6 => 45,
                    _ => 60,
                };
            }
        }
        bubble_p = bubble_max_degree_1d(el.code, bubble_p, pd.serendipity);
        let maxp = (1.max(edge_p).max(face_p).max(bubble_p) + 1) as usize;
        maxp.pow(el.dim() as u32)
    }

    /// `SetMeshMaxDOFs`.
    fn set_mesh_max_dofs(&mut self) {
        for i in 0..self.n_bulk {
            if self.elements[i].pdefs.is_some() {
                let gp = self.number_of_gauss_points(i);
                self.elements[i].pdefs.as_mut().unwrap().gauss_points = gp;
            }
            let el = &self.elements[i];
            self.max_bdofs = self.max_bdofs.max(el.bdofs);
            self.max_ndofs = self.max_ndofs.max(el.ndofs / el.nodes.len());
        }
        for i in 0..self.n_bulk {
            let el = &self.elements[i];
            self.max_element_dofs = self.max_element_dofs.max(
                (el.nodes.len() * self.max_ndofs + el.n_edges() * self.max_edge_dofs + el.bdofs)
                    .max(el.dg_dofs),
            );
            if el.bdofs > 0 {
                self.elements[i].bubble_base = self.max_bdofs * i + 1;
            }
        }
        if self.min_edge_dofs == usize::MAX {
            self.min_edge_dofs = self.max_edge_dofs;
        }
    }

    /// `getEdgeDOFs(Element, p)`: 0 for elements without p definitions.
    pub fn edge_dofs_of(&self, e: usize, p: i32) -> usize {
        if self.elements[e].pdefs.is_none() {
            0
        } else {
            (p - 1).max(0) as usize
        }
    }

    /// `getFaceDOFs(Element, p, face)`.
    pub fn face_dofs_of(&self, e: usize, p: i32, face: usize) -> usize {
        match &self.elements[e].pdefs {
            None => 0,
            Some(pd) => face_dofs(self.elements[e].code, p, face, pd.serendipity),
        }
    }

    /// `getBubbleDOFs(Element, p)`.
    pub fn bubble_dofs_of(&self, e: usize, p: i32) -> usize {
        match &self.elements[e].pdefs {
            None => 0,
            Some(pd) => bubble_dofs(self.elements[e].code, p, pd.serendipity),
        }
    }
}

/// `getEffectiveBubbleP(Element, p, bdofs)`: the smallest degree from `p`
/// up whose bubble count reaches `bdofs`.
pub fn effective_bubble_p(code: u16, p: i32, bdofs: usize, serendipity: bool) -> i32 {
    let mut p = p;
    let mut curr = bubble_dofs(code, p, serendipity);
    let nb = curr.max(bdofs);
    while curr < nb {
        p += 1;
        curr = bubble_dofs(code, p, serendipity);
    }
    p
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::{initial_max_defs, Kind};
    use elmer_core::{ElementRef, ElementType};

    /// Two quads side by side (nodes 0..6) with their outer boundary lines.
    struct TwoQuads {
        bulk: Vec<Vec<u32>>,
        boundary: Vec<(Vec<u32>, u32)>,
    }

    impl TwoQuads {
        fn new() -> Self {
            // 3 4 5
            // 0 1 2
            Self {
                bulk: vec![vec![0, 1, 4, 3], vec![1, 2, 5, 4]],
                boundary: vec![(vec![1, 0], 0), (vec![2, 5], 1), (vec![5, 4], 1)],
            }
        }
    }

    impl MeshView for TwoQuads {
        fn dim(&self) -> usize {
            2
        }
        fn n_nodes(&self) -> usize {
            6
        }
        fn coord(&self, node: u32) -> [f64; 3] {
            [f64::from(node % 3), f64::from(node / 3), 0.0]
        }
        fn n_bulk(&self) -> usize {
            2
        }
        fn n_boundary(&self) -> usize {
            self.boundary.len()
        }
        fn element(&self, e: ElemId) -> ElementRef<'_> {
            let i = e.0 as usize;
            if i < 2 {
                ElementRef {
                    etype: ElementType(404),
                    nodes: &self.bulk[i],
                    tag: i as i32 + 1,
                    parents: [None, None],
                }
            } else {
                let (nodes, p) = &self.boundary[i - 2];
                ElementRef {
                    etype: ElementType(202),
                    nodes,
                    tag: 1,
                    parents: [Some(ElemId(*p)), None],
                }
            }
        }
    }

    fn prepare(def: &str, f: &mut MatcDegree<'_>) -> (PMesh, ElementDefs) {
        let mesh = TwoQuads::new();
        let mut max = initial_max_defs();
        let mut defs = ElementDefs::from_solver_element(Some(def), 2, false, &mut max).unwrap();
        let none = |_: usize| None;
        let pm = {
            let mut solvers = [SolverElements {
                defs: &mut defs,
                element: Some(def),
                dg: false,
                body_element: &none,
            }];
            PMesh::new(
                &mesh,
                &mut solvers,
                &mut max,
                PMeshOptions {
                    serendipity: false,
                    need_edges: None,
                    n_bodies: 2,
                    body_has_equation: &|_| true,
                    matc_degree: f,
                    boundary_body: &|_| 0,
                },
            )
            .unwrap()
        };
        (pm, defs)
    }

    #[test]
    fn edges_in_first_visit_order() {
        let (pm, _) = prepare("p:3", &mut |_, _| unreachable!());
        let nodes: Vec<[u32; 2]> = pm.edges.iter().map(|e| e.nodes).collect();
        // Each edge in its left parent's edge-map order (quad edge 3 is local 4 -> 3).
        assert_eq!(
            nodes,
            vec![[0, 1], [1, 4], [3, 4], [0, 3], [1, 2], [2, 5], [4, 5]]
        );
        assert_eq!(pm.elements[1].edges, vec![4, 5, 6, 1]);
        assert_eq!(pm.edges[1].left, Some(0));
        assert_eq!(pm.edges[1].right, Some(1));
        assert!(pm.edges.iter().all(|e| e.bdofs == 2));
        // Non-serendipity quad p=3: 4 bubbles; p+1 = 4 points per direction.
        assert_eq!(pm.elements[0].bdofs, 4);
        assert_eq!(pm.elements[0].pdefs.as_ref().unwrap().gauss_points, 16);
        assert_eq!(pm.max_bdofs, 4);
        assert_eq!(pm.elements[1].bubble_base, 5);
        assert_eq!(pm.max_element_dofs, 4 + 4 * 2 + 4);
    }

    #[test]
    fn boundary_nodes_follow_the_parent_edge() {
        let (pm, _) = prepare("p:2", &mut |_, _| unreachable!());
        let b = &pm.elements[2];
        // Mesh file order (1, 0) becomes the parent's edge 1 order (0, 1).
        assert_eq!(b.nodes, vec![0, 1]);
        let pd = b.pdefs.as_ref().unwrap();
        assert_eq!(
            (pd.local_number, pd.local_parent, pd.p, pd.is_edge),
            (1, Some(0), 2, true)
        );
        assert_eq!(b.bdofs, 1);
        assert_eq!(b.edges, vec![0]);
        // Edge 2 of the right quad runs 2 -> 5 already.
        assert_eq!(pm.elements[3].nodes, vec![2, 5]);
        // Mesh file (5, 4) becomes edge 3 of the right quad, local 4 -> 3 = (4, 5).
        assert_eq!(pm.elements[4].nodes, vec![4, 5]);
        assert_eq!(pm.elements[4].pdefs.as_ref().unwrap().local_number, 3);
    }

    #[test]
    fn matc_degree_is_a_running_maximum() {
        let (pm, defs) = prepare("p:%f", &mut |f, cx| {
            assert_eq!(f, "f");
            Ok(if cx[0] > 1.0 { 4.0 } else { 2.0 })
        });
        assert_eq!(defs.get(4, 1, Kind::P), 2);
        assert_eq!(defs.get(4, 2, Kind::P), 4);
        assert_eq!(pm.elements[0].pdefs.as_ref().unwrap().p, 2);
        assert_eq!(pm.elements[1].pdefs.as_ref().unwrap().p, 4);
        // The shared edge takes the larger degree.
        assert_eq!(pm.edges[1].bdofs, 3);
        assert_eq!(pm.edges[0].bdofs, 1);
        assert_eq!(pm.max_edge_dofs, 3);
        assert_eq!(pm.elements[0].pdefs.as_ref().unwrap().gauss_points, 25);
    }

    #[test]
    fn effective_bubble_degree() {
        assert_eq!(effective_bubble_p(404, 2, 0, false), 2);
        assert_eq!(effective_bubble_p(404, 2, 4, false), 3);
        assert_eq!(effective_bubble_p(303, 2, 1, false), 3);
    }
}
