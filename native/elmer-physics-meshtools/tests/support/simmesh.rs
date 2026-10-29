//! A `MeshView` over the test mesh whose node coordinates a mapper can move
//! (through `motion::SharedCoords`).

use std::sync::Arc;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};
use elmer_physics_meshtools::motion::SharedCoords;

use super::meshio::Mesh;

pub struct SimMesh {
    pub mesh: Mesh,
    pub coords: Arc<SharedCoords>,
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
        let coords = Arc::new(SharedCoords::new(mesh.coords.clone()));
        Self { mesh, coords, nodes }
    }
}

impl MeshView for SimMesh {
    fn dim(&self) -> usize {
        self.mesh.dim
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords.get(node)
    }
    fn n_bulk(&self) -> usize {
        self.mesh.bulk.len()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.boundary.len()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        let el = if i < self.mesh.bulk.len() {
            &self.mesh.bulk[i]
        } else {
            &self.mesh.boundary[i - self.mesh.bulk.len()]
        };
        ElementRef {
            etype: ElementType(el.etype.0),
            nodes: &self.nodes[i],
            tag: el.tag,
            parents: [None, None],
        }
    }
}
