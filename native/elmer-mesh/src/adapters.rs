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
        self.mesh.dim()
    }

    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }

    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }

    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }

    fn element(&self, e: ElemId) -> ElementRef<'_> {
        self.mesh.element(e)
    }

    fn element_edges(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_edges(e)
    }
}

/// `Mesh` itself is a [`MeshView`]: bulk elements come first, then boundary
/// elements, with Elmer's 1-based file order mapped to 0-based [`ElemId`].
/// This is the form the swarm's `s1_adapters` check asserts.
impl MeshView for Mesh {
    fn dim(&self) -> usize {
        self.dim
    }

    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.nodes[node as usize]
    }

    fn n_bulk(&self) -> usize {
        self.bulk.len()
    }

    fn n_boundary(&self) -> usize {
        self.boundary.len()
    }

    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        if i < self.bulk.len() {
            let el = &self.bulk[i];
            ElementRef {
                etype: ElementType(el.code),
                nodes: &el.nodes,
                tag: el.body,
                parents: [None, None],
            }
        } else {
            let b = &self.boundary[i - self.bulk.len()];
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
        if i < self.bulk.len() {
            Some(&self.bulk[i].edges)
        } else {
            Some(&self.boundary[i - self.bulk.len()].edges)
        }
    }
}

impl Mesh {
    /// Borrow this mesh as an `elmer_core::MeshView`.
    pub fn view(&self) -> MeshViewAdapter<'_> {
        MeshViewAdapter::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{BoundElem, Elem};

    fn two_tri() -> Mesh {
        Mesh {
            dim: 2,
            max_dim: 2,
            nodes: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            node_tags: vec![1, 2, 3, 4],
            bulk: vec![
                Elem { code: 303, body: 1, nodes: vec![0, 1, 2], tag: 1, edges: vec![0, 1, 2] },
                Elem { code: 303, body: 1, nodes: vec![0, 2, 3], tag: 2, edges: vec![3, 4, 5] },
            ],
            boundary: vec![BoundElem {
                code: 202,
                constraint: 1,
                parents: [Some(0), None],
                nodes: vec![0, 1],
                tag: 1,
                edges: vec![0],
            }],
            edges: None,
        }
    }

    fn check_view<V: MeshView>(v: &V) {
        assert_eq!(v.n_bulk(), 2);
        assert_eq!(v.n_boundary(), 1);
        assert_eq!(v.coord(2), [1.0, 1.0, 0.0]);
        // Bulk elements come first, then boundary elements.
        let b = v.element(ElemId(2));
        assert_eq!(b.tag, 1);
        assert_eq!(b.nodes, &[0, 1]);
        assert_eq!(b.parents, [Some(ElemId(0)), None]);
        assert_eq!(v.element_edges(ElemId(0)), Some(&[0, 1, 2][..]));
    }

    /// The swarm's `s1_adapters` check asserts these bounds verbatim.
    #[test]
    fn mesh_and_adapter_are_views() {
        let m = two_tri();
        check_view(&m);
        check_view(&m.view());
    }
}
