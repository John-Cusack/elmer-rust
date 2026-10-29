//! A `MeshView` over the test mesh (until `elmer-mesh` merges).

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

use super::meshio::Mesh;

pub struct SimMesh {
    pub mesh: Mesh,
    nodes: Vec<Vec<u32>>,
}

impl SimMesh {
    pub fn new(mesh: Mesh) -> Self {
        let nodes = mesh
            .bulk
            .iter()
            .chain(&mesh.boundary)
            .map(|e| e.nodes.iter().map(|&n| n as u32).collect())
            .collect();
        Self { mesh, nodes }
    }
}

impl MeshView for SimMesh {
    fn dim(&self) -> usize {
        self.mesh.dim
    }
    fn n_nodes(&self) -> usize {
        self.mesh.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.mesh.bulk.len()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.boundary.len()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        let nb = self.mesh.bulk.len();
        let (el, parents) = if i < nb {
            (&self.mesh.bulk[i], [None, None])
        } else {
            let b = &self.mesh.boundary[i - nb];
            (b, b.parents.map(|p| p.map(|p| ElemId(p as u32))))
        };
        ElementRef {
            etype: ElementType(el.etype.0),
            nodes: &self.nodes[i],
            tag: el.tag,
            parents,
        }
    }
}
