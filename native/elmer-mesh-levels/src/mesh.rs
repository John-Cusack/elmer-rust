use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

#[derive(Clone, Debug, PartialEq)]
pub struct Element {
    pub etype: ElementType,
    /// 0-based node indices, Elmer local order.
    pub nodes: Vec<u32>,
    pub tag: i32,
    /// `BoundaryInfo % Left` and `Right`; boundary elements only.
    pub parents: [Option<ElemId>; 2],
}

/// An owned mesh: bulk elements `0..n_bulk`, then boundary elements.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub elements: Vec<Element>,
    pub n_bulk: usize,
}

impl Mesh {
    pub fn from_view(mesh: &dyn MeshView) -> Self {
        let n = mesh.n_bulk() + mesh.n_boundary();
        let elements = (0..n as u32)
            .map(|e| {
                let el = mesh.element(ElemId(e));
                Element {
                    etype: el.etype,
                    nodes: el.nodes.to_vec(),
                    tag: el.tag,
                    parents: el.parents,
                }
            })
            .collect();
        Mesh {
            dim: mesh.dim(),
            coords: (0..mesh.n_nodes() as u32).map(|i| mesh.coord(i)).collect(),
            elements,
            n_bulk: mesh.n_bulk(),
        }
    }

    pub fn bulk(&self) -> &[Element] {
        &self.elements[..self.n_bulk]
    }

    pub fn boundary(&self) -> &[Element] {
        &self.elements[self.n_bulk..]
    }
}

impl MeshView for Mesh {
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
        self.elements.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let el = &self.elements[e.0 as usize];
        ElementRef {
            etype: el.etype,
            nodes: &el.nodes,
            tag: el.tag,
            parents: el.parents,
        }
    }
}
