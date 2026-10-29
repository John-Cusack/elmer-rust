//! Mesh-level p-element setup: what `PrepareMesh` (MeshLoad.F90:643) adds to
//! a loaded mesh before any solver builds its matrix.
//!
//! * `NonNodalElements` (MeshLoad.F90:713): per-element nodal and bubble DOF
//!   counts and the p definitions (`Element % PDefs`), including the
//!   MATC-dependent degree `p:%f` evaluated element by element;
//! * `FindMeshEdges2D` (MeshBasics.F90:3875), or `FindMeshFaces3D`
//!   (MeshBasics.F90:4135) and `FindMeshEdges3D` (MeshBasics.F90:4558):
//!   the global edge and face tables, numbered in first-visit order over
//!   bulk then boundary elements;
//! * `SetMeshEdgeFaceDOFs` (MeshBasics.F90:2123) and `AssignLocalNumber`
//!   (MeshBasics.F90:5826): edge and face degrees and DOF counts, and the
//!   local edge or face (`PDefs % LocalNumber`) and parent of every p
//!   boundary element, whose nodes are **reordered** to the parent's edge
//!   or face orientation (in 3D its edges become the parent's face edges);
//! * `SetMeshMaxDOFs` (MeshBasics.F90:2271): integration point counts
//!   (`getNumberOfGaussPoints`, PElementMaps.F90:1135), the mesh maxima and
//!   the bubble index blocks.
//!
//! Supported: meshes of first-order elements; DOF tables for every type,
//! the p basis (see [`crate::basis`]) for `202`, `303`, `404`, `504`, `605`,
//! `706` and `808`.
//! p tetrahedra are reordered by `ConvertToACTetra` (MeshBasics.F90:5747),
//! which sets their `TetraType`. Element definitions come from the solver
//! sections, the equations (`Element`, `Element{id}`) and the bodies
//! (`Solver id: Element`), with `NonNodalElements`' precedence and its
//! per-element re-reading for MATC degrees.

use elmer_core::{ElemId, MeshView};
use elmer_fe_hp::maps::{
    bubble_dofs, bubble_max_degree_1d, edge_map, face_dofs, face_edge_map, face_map,
    pelement_edge_map, tetra_face_direction,
};

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
    /// `FaceIndexes`, 0-based global faces (3D meshes); empty when not
    /// associated.
    pub faces: Vec<u32>,
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

/// A global face of a 3D mesh (`Mesh % Faces(i)`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MeshFace {
    /// `303` or `404`.
    pub code: u16,
    /// Nodes in the face map order of the first element visiting it.
    pub nodes: Vec<u32>,
    /// `BoundaryInfo % Left / Right`: the first two bulk elements visiting
    /// the face.
    pub left: Option<u32>,
    pub right: Option<u32>,
    /// `EdgeIndexes`, 0-based, in face-edge-map order (0 entries of the
    /// Fortran are `u32::MAX` until set).
    pub edges: Vec<u32>,
    pub pdefs: Option<PDefs>,
    /// `BDOFs`: face DOF count.
    pub bdofs: usize,
}

fn n_faces(code: u16) -> usize {
    match code / 100 {
        5 => 4,
        6 | 7 => 5,
        8 => 6,
        _ => 0,
    }
}

/// Nodes of face `k` (0-based) of a 3D element type.
fn face_node_count(code: u16, k: usize) -> usize {
    match code / 100 {
        3 | 5 => 3,
        4 | 8 => 4,
        6 => {
            if k == 0 {
                4
            } else {
                3
            }
        }
        7 => {
            if k <= 1 {
                3
            } else {
                4
            }
        }
        _ => 0,
    }
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

/// An element's edge map (1-based end nodes) and face edge map (1-based
/// local edges per face).
type EdgeMaps = (Vec<[usize; 2]>, Vec<Vec<usize>>);

/// `GetElementFaceMap` / `getElementBoundaryMap` of a p element: the
/// tetra-type tables for tetrahedra.
fn p_face_map(el: &PElement) -> Result<elmer_fe_hp::FaceDir, Error> {
    if el.code == 504 {
        return Ok(tetra_face_direction(pd_tetra_type(el)));
    }
    face_map(el.code).ok_or_else(|| Error::Unsupported(format!("face map of type {}", el.code)))
}

/// `ConvertToACTetra` (MeshBasics.F90:5747): reorders a p tetrahedron's
/// nodes so that face 1 (face 2 when node 4 is the lowest) starts at the
/// lowest node and face 3 ends at the highest, and returns the resulting
/// `TetraType` (1 when node 2 < node 3, else 2).
pub fn convert_to_ac_tetra(nodes: &mut [u32]) -> Result<i32, Error> {
    let type1 = tetra_face_direction(1);
    let face = |k: usize| [type1[k][0] - 1, type1[k][1] - 1, type1[k][2] - 1];
    let mut gmin = nodes[0];
    let mut gmin_i = 0;
    let mut gmax = nodes[0];
    for (i, &n) in nodes.iter().enumerate().take(4).skip(1) {
        if gmin > n {
            gmin = n;
            gmin_i = i;
        } else if gmax < n {
            gmax = n;
        }
    }
    let f = if gmin_i == 3 { face(1) } else { face(0) };
    let mut g = f.map(|l| nodes[l]);
    while g[0] != gmin {
        g.rotate_left(1);
    }
    for (k, &l) in f.iter().enumerate() {
        nodes[l] = g[k];
    }
    let f = face(2);
    let mut g = f.map(|l| nodes[l]);
    while g[2] != gmax {
        g.rotate_left(1);
    }
    for (k, &l) in f.iter().enumerate() {
        nodes[l] = g[k];
    }
    match nodes[1].cmp(&nodes[2]) {
        std::cmp::Ordering::Less => Ok(1),
        std::cmp::Ordering::Greater => Ok(2),
        std::cmp::Ordering::Equal => Err(Error::BadDefinition(
            "ConvertToACTetra: Corrupt element type".into(),
        )),
    }
}

fn pd_tetra_type(el: &PElement) -> i32 {
    el.pdefs.as_ref().map_or(1, |pd| pd.tetra_type.max(1))
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
    /// `Element{i}` in an equation section (1-based equation), when present.
    pub equation_element: &'a dyn Fn(usize) -> Option<String>,
}

/// MATC for `p:%function` degrees: `p_of(command, [body, x, y, z])`
/// evaluates `command` (`function(cx)`) with the MATC variable `cx` set to
/// the body id and element midpoint, as `GetMatcReal(command, 4, cx, 'cx')`.
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
    /// Number of equation sections.
    pub n_equations: usize,
    /// The `Equation` of body `b` (1-based), when it has one.
    pub body_equation: &'a dyn Fn(usize) -> Option<usize>,
    /// `Element` of equation `q` (1-based), when present.
    pub equation_element: &'a dyn Fn(usize) -> Option<String>,
    /// MATC degree evaluation (see [`MatcDegree`]); called with the
    /// command `function(cx)`.
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
    pub faces: Vec<MeshFace>,
    /// `Mesh % NumberOfFaces`.
    pub n_faces: usize,
    pub max_ndofs: usize,
    pub max_edge_dofs: usize,
    pub min_edge_dofs: usize,
    pub max_face_dofs: usize,
    pub min_face_dofs: usize,
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
                faces: Vec::new(),
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
            faces: Vec::new(),
            n_faces: 0,
            max_ndofs: 0,
            max_edge_dofs: 0,
            min_edge_dofs: usize::MAX,
            max_face_dofs: 0,
            min_face_dofs: usize::MAX,
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
        // NonNodalElements (MeshLoad.F90:787): which definitions exist and
        // whether any degree comes from MATC.
        let ns = solvers.len();
        let mut mesh_deps = false;
        let mut found_eq_defs = false;
        let mut found_solver_defs = vec![false; ns];
        for q in 1..=opts.n_equations {
            match (opts.equation_element)(q) {
                Some(def) => {
                    found_eq_defs = true;
                    mesh_deps |= is_mesh_dependent(&def);
                }
                None => {
                    for (si, s) in solvers.iter().enumerate() {
                        if let Some(def) = (s.equation_element)(q) {
                            found_solver_defs[si] = true;
                            mesh_deps |= is_mesh_dependent(&def);
                        }
                    }
                }
            }
        }
        // UpdateDefDofs: solvers re-read per body (or per element).
        let mut update = vec![true; ns];
        for (si, s) in solvers.iter().enumerate() {
            match s.element {
                Some(def) => {
                    found_solver_defs[si] = true;
                    if is_mesh_dependent(def) {
                        mesh_deps = true;
                    } else {
                        update[si] = false;
                    }
                }
                None => {
                    if !found_eq_defs {
                        update[si] = false;
                    }
                }
            }
        }

        // str = TRIM(ElementDef(j+3:))//'(cx)'; x = GetMatcReal(str, 4,
        // [BodyId, x, y, z], 'cx'); NINT(x).
        let mut matc = |f: &str, cx: [f64; 4]| -> Result<i32, Error> {
            Ok((opts.matc_degree)(&format!("{f}(cx)"), cx)?.round() as i32)
        };

        if !mesh_deps {
            let no_matc = |_: &str| -> Result<i32, Error> {
                Err(Error::Unsupported(
                    "MATC degree outside the per-element loop".into(),
                ))
            };
            for body in 1..=opts.n_bodies {
                let Some(eq) = (opts.body_equation)(body) else {
                    continue;
                };
                let (mut def0, found_def0) =
                    match found_eq_defs.then(|| (opts.equation_element)(eq)).flatten() {
                        Some(d) => (d, true),
                        None => (String::new(), false),
                    };
                for (si, s) in solvers.iter_mut().enumerate() {
                    let dg = s.dg;
                    if let Some(def) = (s.body_element)(body) {
                        s.defs
                            .get_max_defs(&def, body, dg, in_dofs, &mut |f| no_matc(f))?;
                        continue;
                    }
                    let def = if found_solver_defs[si] {
                        (s.equation_element)(eq)
                    } else {
                        None
                    };
                    if let Some(def) = def {
                        s.defs
                            .get_max_defs(&def, body, dg, in_dofs, &mut |f| no_matc(f))?;
                    } else if update[si] {
                        let own = !found_def0 && found_solver_defs[si];
                        if own {
                            def0 = s.element.unwrap_or("").to_string();
                        }
                        s.defs
                            .get_max_defs(&def0, body, dg, in_dofs, &mut |f| no_matc(f))?;
                        if own {
                            def0.clear();
                        }
                    }
                }
            }
        }

        let mut need_edges = in_dofs[8][4] > 0 || in_dofs[9][4] > 0;
        let mut body0 = i32::MIN;
        let mut found_eq: Option<usize> = None;
        let mut def0 = String::new();
        for i in 0..n_bulk {
            let body = self.elements[i].body;
            if mesh_deps {
                if body != body0 {
                    found_eq = if body >= 1 {
                        (opts.body_equation)(body as usize)
                    } else {
                        None
                    };
                    def0.clear();
                }
                if let Some(eq) = found_eq {
                    // FoundDef0 is reset for every element; the equation's
                    // Element is read only on the first element of a body.
                    let mut found_def0 = false;
                    if found_eq_defs && body != body0 {
                        if let Some(d) = (opts.equation_element)(eq) {
                            def0 = d;
                            found_def0 = true;
                        }
                    }
                    let mid = self.midpoint(i);
                    let cx = [f64::from(body), mid[0], mid[1], mid[2]];
                    for (si, s) in solvers.iter_mut().enumerate() {
                        let dg = s.dg;
                        let def = if found_solver_defs[si] {
                            (s.equation_element)(eq)
                        } else {
                            None
                        };
                        if let Some(def) = def {
                            s.defs
                                .get_max_defs(&def, body as usize, dg, in_dofs, &mut |f| {
                                    matc(f, cx)
                                })?;
                        } else if update[si] {
                            let own = !found_def0 && found_solver_defs[si];
                            if own {
                                def0 = s.element.unwrap_or("").to_string();
                            }
                            s.defs
                                .get_max_defs(&def0, body as usize, dg, in_dofs, &mut |f| {
                                    matc(f, cx)
                                })?;
                            if own {
                                def0.clear();
                            }
                        }
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
                    // The Ainsworth & Coyle trick.
                    pd.tetra_type = convert_to_ac_tetra(&mut el.nodes)?;
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
            self.set_mesh_edge_face_dofs(&edge_dofs, &face_dofs, in_dofs)?;
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
            self.faces.clear();
            self.n_faces = 0;
            for el in &mut self.elements {
                el.edges.clear();
                el.faces.clear();
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

    /// The face map `FindMeshFaces3D` uses for element `i` (1-based local
    /// nodes per face) and its face count: `GetElementFaceMap` for p
    /// elements, the routine's own tables otherwise.
    fn element_face_map(&self, i: usize) -> Result<Vec<Vec<usize>>, Error> {
        let el = &self.elements[i];
        let code = el.code;
        if code / 100 == 3 || code / 100 == 4 {
            let n = face_node_count(code, 0);
            return Ok(vec![(1..=n).collect()]);
        }
        let rows: Vec<Vec<usize>> = if el.pdefs.is_some() {
            let map = p_face_map(el)?;
            (0..n_faces(code))
                .map(|k| map[k][..face_node_count(code, k)].to_vec())
                .collect()
        } else {
            let t: &[&[usize]] = match code / 100 {
                5 => &[&[1, 2, 3], &[1, 2, 4], &[2, 3, 4], &[3, 1, 4]],
                6 => &[
                    &[1, 2, 3, 4],
                    &[1, 2, 5],
                    &[2, 3, 5],
                    &[3, 4, 5],
                    &[4, 1, 5],
                ],
                7 => &[
                    &[1, 2, 3],
                    &[4, 5, 6],
                    &[1, 2, 5, 4],
                    &[3, 2, 5, 6],
                    &[3, 1, 4, 6],
                ],
                8 => &[
                    &[1, 2, 3, 4],
                    &[5, 6, 7, 8],
                    &[1, 2, 6, 5],
                    &[2, 3, 7, 6],
                    &[3, 4, 8, 7],
                    &[4, 1, 5, 8],
                ],
                _ => {
                    return Err(Error::Unsupported(format!(
                        "FindMeshFaces: element type {code}"
                    )))
                }
            };
            t.iter().map(|r| r.to_vec()).collect()
        };
        Ok(rows)
    }

    /// `FindMeshFaces3D` without a mask: faces numbered in first-visit
    /// order over bulk then boundary elements (types 300 and up), keyed by
    /// their three lowest node numbers.
    fn find_mesh_faces_3d(&mut self) -> Result<(), Error> {
        let mut table: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); self.coords.len()];
        let mut faces: Vec<MeshFace> = Vec::new();
        for i in 0..self.elements.len() {
            let code = self.elements[i].code;
            if code < 300 {
                continue;
            }
            let map = self.element_face_map(i)?;
            let mut idx = vec![0u32; map.len()];
            for (k, local) in map.iter().enumerate() {
                let nodes: Vec<u32> = local
                    .iter()
                    .map(|&l| self.elements[i].nodes[l - 1])
                    .collect();
                let mut nf = nodes.clone();
                nf.sort_unstable();
                let found = table[nf[0] as usize]
                    .iter()
                    .find(|&&(a, b, _)| a == nf[1] && b == nf[2])
                    .map(|&(_, _, f)| f);
                let face = match found {
                    Some(f) => f,
                    None => {
                        let f = faces.len() as u32;
                        table[nf[0] as usize].push((nf[1], nf[2], f));
                        let pdefs = self.elements[i]
                            .pdefs
                            .as_ref()
                            .map(|pd| PDefs::new(pd.serendipity));
                        faces.push(MeshFace {
                            code: if nodes.len() == 3 { 303 } else { 404 },
                            nodes,
                            left: None,
                            right: None,
                            edges: Vec::new(),
                            pdefs,
                            bdofs: 0,
                        });
                        f
                    }
                };
                idx[k] = face;
                if i < self.n_bulk {
                    let fc = &mut faces[face as usize];
                    if fc.left.is_some() {
                        fc.right = Some(i as u32);
                    } else {
                        fc.left = Some(i as u32);
                    }
                }
            }
            self.elements[i].faces = idx;
        }
        self.faces = faces;
        self.n_faces = self.faces.len();
        Ok(())
    }

    /// The edge map (1-based local end nodes) and, for 3D types, the face
    /// edge map (1-based local edges per face) `FindMeshEdges3D` uses for
    /// element `i`.
    fn element_edge_maps(&self, i: usize) -> Result<EdgeMaps, Error> {
        let el = &self.elements[i];
        let code = el.code;
        let ne = n_edges(code);
        if el.pdefs.is_some() {
            let em = pelement_edge_map(code, pd_tetra_type(el))
                .ok_or_else(|| Error::Unsupported(format!("edge map of type {code}")))?;
            let fem = if code >= 500 {
                let m = face_edge_map(code, pd_tetra_type(el))
                    .ok_or_else(|| Error::Unsupported(format!("face edge map of type {code}")))?;
                (0..n_faces(code))
                    .map(|k| m[k].iter().copied().filter(|&x| x > 0).collect())
                    .collect()
            } else {
                Vec::new()
            };
            return Ok((em[..ne].to_vec(), fem));
        }
        let (em, fem): (&[[usize; 2]], &[&[usize]]) = match code / 100 {
            2 => (&[[1, 2]], &[]),
            3 => (&[[1, 2], [2, 3], [3, 1]], &[]),
            4 => (&[[1, 2], [2, 3], [3, 4], [4, 1]], &[]),
            5 => (
                &[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]],
                &[&[1, 2, 3], &[1, 5, 4], &[2, 6, 5], &[3, 4, 6]],
            ),
            6 => (
                &[
                    [1, 2],
                    [2, 3],
                    [3, 4],
                    [4, 1],
                    [1, 5],
                    [2, 5],
                    [3, 5],
                    [4, 5],
                ],
                &[
                    &[1, 2, 3, 4],
                    &[1, 6, 5],
                    &[2, 7, 6],
                    &[3, 8, 7],
                    &[4, 5, 8],
                ],
            ),
            7 => (
                &[
                    [1, 2],
                    [2, 3],
                    [1, 3],
                    [4, 5],
                    [5, 6],
                    [6, 4],
                    [1, 4],
                    [2, 5],
                    [3, 6],
                ],
                &[
                    &[1, 2, 3],
                    &[4, 5, 6],
                    &[1, 8, 4, 7],
                    &[2, 9, 5, 8],
                    &[3, 7, 6, 9],
                ],
            ),
            8 => (
                &[
                    [1, 2],
                    [2, 3],
                    [4, 3],
                    [1, 4],
                    [5, 6],
                    [6, 7],
                    [8, 7],
                    [5, 8],
                    [1, 5],
                    [2, 6],
                    [3, 7],
                    [4, 8],
                ],
                &[
                    &[1, 2, 3, 4],
                    &[5, 6, 7, 8],
                    &[1, 10, 5, 9],
                    &[2, 11, 6, 10],
                    &[3, 12, 7, 11],
                    &[4, 9, 8, 12],
                ],
            ),
            _ => {
                return Err(Error::Unsupported(format!(
                    "FindMeshEdges3D: element type {code}"
                )))
            }
        };
        Ok((em.to_vec(), fem.iter().map(|r| r.to_vec()).collect()))
    }

    /// `FindMeshEdges3D`: edges numbered in first-visit order over bulk then
    /// boundary elements; the faces of 3D elements get their edges through
    /// the face edge maps. (Edge `BoundaryInfo` is not kept: in 3D Elmer
    /// fills it with elements and faces alike and only uses it to orient
    /// the edge's own node order, which no DOF or basis reads.)
    fn find_mesh_edges_3d(&mut self) -> Result<(), Error> {
        let mut table: Vec<Vec<(u32, u32)>> = vec![Vec::new(); self.coords.len()];
        let mut edges: Vec<MeshEdge> = Vec::new();
        for f in &mut self.faces {
            f.edges = vec![u32::MAX; if f.code == 303 { 3 } else { 4 }];
        }
        for i in 0..self.elements.len() {
            if self.elements[i].code < 200 {
                continue;
            }
            let (em, fem) = self.element_edge_maps(i)?;
            let mut idx = vec![0u32; em.len()];
            for (k, &[a, b]) in em.iter().enumerate() {
                let (na, nb) = (self.elements[i].nodes[a - 1], self.elements[i].nodes[b - 1]);
                let (n1, n2) = if na < nb { (na, nb) } else { (nb, na) };
                let found = table[n1 as usize]
                    .iter()
                    .find(|&&(n, _)| n == n2)
                    .map(|&(_, e)| e);
                let edge = match found {
                    Some(e) => e,
                    None => {
                        let e = edges.len() as u32;
                        table[n1 as usize].push((n2, e));
                        let pdefs = self.elements[i]
                            .pdefs
                            .as_ref()
                            .map(|pd| PDefs::new(pd.serendipity));
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
                idx[k] = edge;
                for (ii, row) in fem.iter().enumerate() {
                    let face = self.elements[i].faces[ii] as usize;
                    if let Some(jj) = row.iter().position(|&x| x == k + 1) {
                        self.faces[face].edges[jj] = edge;
                    }
                }
            }
            self.elements[i].edges = idx;
        }
        self.edges = edges;
        Ok(())
    }

    /// `SetMeshEdgeFaceDOFs` with `EdgeDOFs`/`FaceDOFs`/`inDOFs` present and
    /// `NeedEdges` absent.
    fn set_mesh_edge_face_dofs(
        &mut self,
        edge_dofs: &[usize],
        face_dofs: &[usize],
        in_dofs: &MaxDefs,
    ) -> Result<(), Error> {
        match self.dim {
            1 | 2 => self.find_mesh_edges_2d(),
            3 => {
                self.find_mesh_faces_3d()?;
                self.find_mesh_edges_3d()?;
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
                    // In 3D the edge's parent is a face or an element, and
                    // AssignLocalNumber only reorders the edge's own nodes.
                    let parent = ed.left.or(ed.right);
                    if let (Some(parent), true) = (parent, self.dim < 3) {
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
            let el_faces = self.elements[i].faces.clone();
            for (j, &face) in el_faces.iter().enumerate() {
                let f = face as usize;
                if self.faces[f].edges.contains(&u32::MAX) {
                    continue;
                }
                if let Some(pd) = self.elements[i].pdefs.clone() {
                    let code = self.elements[i].code;
                    let fc = &mut self.faces[f];
                    let fp = fc.pdefs.get_or_insert_with(|| PDefs::new(pd.serendipity));
                    fp.p = pd.p.max(fp.p);
                    let p = fp.p;
                    fp.is_edge = true;
                    fc.bdofs =
                        fc.bdofs
                            .max(elmer_fe_hp::maps::face_dofs(code, p, j + 1, pd.serendipity));
                    let gp = self.face_gauss_points(f);
                    self.faces[f]
                        .pdefs
                        .as_mut()
                        .expect("set above")
                        .gauss_points = gp;
                    // AssignLocalNumber(Face, Left) only reorders the
                    // face's own nodes, which no DOF or basis reads.
                }
                let fc = &mut self.faces[f];
                let fam = usize::from(fc.code / 100);
                fc.bdofs = fc
                    .bdofs
                    .max(face_dofs[i])
                    .max(in_dofs[fam + 5][4].max(0) as usize);
                self.min_face_dofs = self.min_face_dofs.min(fc.bdofs);
                self.max_face_dofs = self.max_face_dofs.max(fc.bdofs);
            }
        }
        if self.min_face_dofs > self.max_face_dofs {
            self.min_face_dofs = self.max_face_dofs;
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
                    self.assign_local_number_boundary(i, l as usize)?;
                }
            }
            if let Some(r) = right {
                if let Some(pd) = &self.elements[r as usize].pdefs {
                    if !p_alloc {
                        let mut d = PDefs::new(pd.serendipity);
                        d.is_edge = true;
                        self.elements[i].pdefs = Some(d);
                        self.assign_local_number_boundary(i, r as usize)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// `getNumberOfGaussPointsFace(Face, Mesh)`.
    fn face_gauss_points(&self, f: usize) -> usize {
        let fc = &self.faces[f];
        let edge_p = fc
            .edges
            .iter()
            .map(|&e| self.edges[e as usize].pdefs.as_ref().map_or(0, |d| d.p))
            .max()
            .unwrap_or(0);
        let pd = fc.pdefs.as_ref().expect("p face");
        let maxp = if fc.bdofs == 0 {
            edge_p
        } else {
            edge_p.max(pd.p)
        };
        if pd.serendipity && fc.code / 100 == 4 && maxp > 1 && maxp <= 8 {
            return [8, 12, 20, 25, 36, 45, 60][maxp as usize - 2];
        }
        ((maxp + 1) as usize).pow(2)
    }

    /// `AssignLocalNumber(Element, Parent)` for a boundary face element of
    /// a 3D parent: local face and parent, nodes reordered to the parent's
    /// face map, the face's p definitions and DOF count, and the parent's
    /// edges of that face (`getFaceEdgeMap`) as the element's edges.
    fn assign_local_number_face(&mut self, i: usize, parent: usize) -> Result<(), Error> {
        let pel = &self.elements[parent];
        let code = pel.code;
        let map = p_face_map(pel)?;
        let fem = face_edge_map(code, pd_tetra_type(pel))
            .ok_or_else(|| Error::Unsupported(format!("face edge map of type {code}")))?;
        let nodes = self.elements[i].nodes.clone();
        for k in 0..n_faces(code) {
            let face = &self.faces[pel.faces[k] as usize];
            let n = face.nodes.iter().filter(|g| nodes.contains(g)).count();
            if n != nodes.len() {
                continue;
            }
            let bindex: Vec<u32> = map[k][..n].iter().map(|&l| pel.nodes[l - 1]).collect();
            let hits = face.nodes.iter().filter(|g| bindex.contains(g)).count();
            let (is_edge, gp, p) = match &face.pdefs {
                Some(pd) => (pd.is_edge, pd.gauss_points, pd.p),
                None => (false, 0, 0),
            };
            let face_bdofs = face.bdofs;
            let edges: Vec<u32> = if face.edges.is_empty() {
                Vec::new()
            } else {
                (0..face.edges.len())
                    .map(|j| pel.edges[fem[k][j] - 1])
                    .collect()
            };
            let el = &mut self.elements[i];
            if hits == n {
                el.nodes[..n].copy_from_slice(&bindex);
            }
            let pd = el.pdefs.as_mut().expect("p boundary element");
            pd.local_number = k + 1;
            pd.local_parent = Some(parent as u32);
            pd.is_edge = is_edge;
            pd.gauss_points = gp;
            pd.p = p;
            el.bdofs = el.bdofs.max(face_bdofs);
            if !edges.is_empty() {
                el.edges = edges;
            }
            return Ok(());
        }
        Ok(())
    }

    /// Local edge of bulk element `parent` (1-based) whose nodes all belong
    /// to `nodes`, with the parent's nodes along that edge.
    fn find_local_edge(&self, nodes: &[u32], parent: usize) -> Option<(usize, [u32; 2])> {
        let pel = &self.elements[parent];
        if pel.dim() != 2 {
            // AssignLocalNumber: nothing to do for 0D/1D parents; 3D parents
            // go through `assign_local_number_face`.
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
    fn assign_local_number_boundary(&mut self, i: usize, parent: usize) -> Result<(), Error> {
        match self.elements[parent].dim() {
            0 | 1 => return Ok(()),
            3 => return self.assign_local_number_face(i, parent),
            _ => {}
        }
        let nodes = self.elements[i].nodes.clone();
        let Some((num, bindex)) = self.find_local_edge(&nodes, parent) else {
            return Ok(());
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
        Ok(())
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
        let face_p = if el.dim() == 3 {
            el.faces
                .iter()
                .map(|&f| &self.faces[f as usize])
                .filter(|f| f.bdofs > 0)
                .map(|f| f.pdefs.as_ref().map_or(0, |d| d.p))
                .max()
                .unwrap_or(0)
        } else {
            0
        };
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
                (el.nodes.len() * self.max_ndofs
                    + el.n_edges() * self.max_edge_dofs
                    + n_faces(el.code) * self.max_face_dofs
                    + el.bdofs)
                    .max(el.dg_dofs),
            );
            if el.bdofs > 0 {
                self.elements[i].bubble_base = self.max_bdofs * i + 1;
            }
        }
        if self.min_face_dofs == usize::MAX {
            self.min_face_dofs = self.max_face_dofs;
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
                equation_element: &none,
            }];
            PMesh::new(
                &mesh,
                &mut solvers,
                &mut max,
                PMeshOptions {
                    serendipity: false,
                    need_edges: None,
                    n_bodies: 2,
                    n_equations: 1,
                    body_equation: &|_| Some(1),
                    equation_element: &|_| None,
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
            assert_eq!(f, "f(cx)");
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
    fn ac_tetra_orders_faces_by_global_index() {
        // Face 1 (nodes 1 2 3) starts at the lowest node; face 3 (2 3 4)
        // ends at the highest.
        // [5, 9, 2, 7]: face 1 rotates to [2, 5, 9], then face 3 (5, 9, 7)
        // to [7, 5, 9]; node 2 > node 3 gives type 2.
        let mut n = [5, 9, 2, 7];
        assert_eq!(convert_to_ac_tetra(&mut n).unwrap(), 2);
        assert_eq!(n, [2, 7, 5, 9]);
        let mut n = [5, 9, 2, 7];
        n.swap(1, 3);
        // [5, 7, 2, 9]: face 1 rotates to [2, 5, 7], then face 3 to [5, 7, 9].
        assert_eq!(convert_to_ac_tetra(&mut n).unwrap(), 1);
        assert_eq!(n, [2, 5, 7, 9]);
        // Lowest node in position 4: face 2 (1 2 4) is rotated instead.
        let mut n = [5, 9, 7, 1];
        convert_to_ac_tetra(&mut n).unwrap();
        assert_eq!(n[0], 1);
    }

    #[test]
    fn effective_bubble_degree() {
        assert_eq!(effective_bubble_p(404, 2, 0, false), 2);
        assert_eq!(effective_bubble_p(404, 2, 4, false), 3);
        assert_eq!(effective_bubble_p(303, 2, 1, false), 3);
    }
}
