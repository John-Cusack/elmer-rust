//! Adapter from the native `elmer-mesh` model to `elmer_core::MeshView`.

use crate::mesh::Mesh;
use elmer_core::mesh::{ElemId, ElementRef, ElementType, MeshView};

/// Thin view over [`Mesh`]; no copying.
pub struct MeshViewAdapter<'a> {
    mesh: &'a Mesh,
}

impl<'a> MeshViewAdapter<'a> {
    pub fn new(mesh: &'a Mesh) -> Self {
        MeshViewAdapter { mesh }
    }
}

impl MeshView for MeshViewAdapter<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim
    }

    fn n_nodes(&self) -> usize {
        self.mesh.nodes.len()
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.nodes[node as usize]
    }

    fn n_bulk(&self) -> usize {
        self.mesh.bulk.len()
    }

    fn n_boundary(&self) -> usize {
        self.mesh.boundary.len()
    }

    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        if i < self.mesh.bulk.len() {
            let el = &self.mesh.bulk[i];
            ElementRef {
                etype: ElementType(el.code),
                nodes: &el.nodes,
                tag: el.body,
                parents: [None, None],
            }
        } else {
            let b = &self.mesh.boundary[i - self.mesh.bulk.len()];
            ElementRef {
                etype: ElementType(b.code),
                nodes: &b.nodes,
                tag: b.constraint,
                parents: [
                    b.parents[0].map(|p| ElemId(p as u32)),
                    b.parents[1].map(|p| ElemId(p as u32)),
                ],
            }
        }
    }

    fn element_edges(&self, e: ElemId) -> Option<&[u32]> {
        let i = e.0 as usize;
        if i < self.mesh.bulk.len() {
            Some(&self.mesh.bulk[i].edges)
        } else {
            Some(&self.mesh.boundary[i - self.mesh.bulk.len()].edges)
        }
    }
}

impl Mesh {
    /// Borrow this mesh as an `elmer_core::MeshView`.
    pub fn view(&self) -> MeshViewAdapter<'_> {
        MeshViewAdapter::new(self)
    }
}
