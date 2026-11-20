//! The meshes of the refinement: an owned triangle mesh (`Mesh_t` as
//! SplitOneLevel builds it) and the 2-D edge table of `FindMeshEdges2D`
//! (MeshBasics.F90).

use std::collections::HashMap;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

/// One element of an owned mesh: type, nodes (0-based), body or boundary
/// tag, and the bulk parents of a boundary element.
#[derive(Clone, Debug)]
pub struct Elem {
    pub etype: ElementType,
    pub nodes: Vec<u32>,
    pub tag: i32,
    pub parents: [Option<ElemId>; 2],
}

/// A mesh the refinement made. Bulk elements come first.
#[derive(Clone, Debug)]
pub struct OwnedMesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub elems: Vec<Elem>,
    pub n_bulk: usize,
}

impl OwnedMesh {
    /// A copy of any mesh.
    pub fn from_view(m: &dyn MeshView) -> OwnedMesh {
        let elems = (0..(m.n_bulk() + m.n_boundary()) as u32)
            .map(|e| {
                let r = m.element(ElemId(e));
                Elem {
                    etype: r.etype,
                    nodes: r.nodes.to_vec(),
                    tag: r.tag,
                    parents: r.parents,
                }
            })
            .collect();
        OwnedMesh {
            dim: m.dim(),
            coords: (0..m.n_nodes() as u32).map(|k| m.coord(k)).collect(),
            elems,
            n_bulk: m.n_bulk(),
        }
    }
}

impl MeshView for OwnedMesh {
    fn dim(&self) -> usize {
        self.dim
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let el = &self.elems[e.0 as usize];
        ElementRef {
            etype: el.etype,
            nodes: &el.nodes,
            tag: el.tag,
            parents: el.parents,
        }
    }
}

/// One edge of `FindMeshEdges2D`: nodes in the order of the first element
/// with it, and its bulk parents (`BoundaryInfo % Left`, `Right`).
#[derive(Clone, Debug)]
pub struct Edge {
    pub nodes: [u32; 2],
    pub left: Option<usize>,
    pub right: Option<usize>,
}

/// `Mesh % Edges` and every bulk element's `EdgeIndexes`.
pub struct Edges {
    pub edges: Vec<Edge>,
    pub of_element: Vec<Vec<usize>>,
}

impl Edges {
    /// `FindMeshEdges2D` over the bulk elements (triangles and
    /// quadrilaterals) and then the boundary elements, which only meet
    /// edges the bulk already made.
    pub fn find(mesh: &dyn MeshView) -> Edges {
        let nb = mesh.n_bulk();
        let mut edges: Vec<Edge> = Vec::new();
        let mut hash: HashMap<(u32, u32), usize> = HashMap::new();
        let mut of_element = Vec::with_capacity(nb);
        for t in 0..nb + mesh.n_boundary() {
            let el = mesh.element(ElemId(t as u32));
            let n = match el.etype.0 / 100 {
                2 => 1,
                3 => 3,
                4 => 4,
                _ => 0,
            };
            let mut idx = Vec::with_capacity(n);
            for k in 0..n {
                let a = el.nodes[k];
                let b = if n == 1 { el.nodes[1] } else { el.nodes[(k + 1) % n] };
                let key = (a.min(b), a.max(b));
                let e = *hash.entry(key).or_insert_with(|| {
                    edges.push(Edge {
                        nodes: [a, b],
                        left: None,
                        right: None,
                    });
                    edges.len() - 1
                });
                if t < nb {
                    if edges[e].left.is_some() {
                        edges[e].right = Some(t);
                    } else {
                        edges[e].left = Some(t);
                    }
                }
                idx.push(e);
            }
            if t < nb {
                of_element.push(idx);
            }
        }
        Edges { edges, of_element }
    }
}
