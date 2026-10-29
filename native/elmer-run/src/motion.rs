//! Node coordinates the mesh mappers can move.
//!
//! RigidMeshMapper and StructuredMeshMapper assign `Mesh % Nodes` in place;
//! `elmer-physics-meshtools` reaches the coordinates through its
//! thread-local [`motion::install`]. [`MovableMesh`] keeps them in a
//! [`SharedCoords`] that every solver and the VTU writer read, and
//! [`MotionGuard`] installs it for the run.

use std::sync::Arc;

use elmer_core::{ElemId, ElementRef, MeshView};
use elmer_physics_meshtools::motion::{self, MeshMotion, SharedCoords};

/// `mesh` with its node coordinates in a [`SharedCoords`].
pub(crate) struct MovableMesh<'a> {
    mesh: &'a dyn MeshView,
    coords: Arc<SharedCoords>,
}

impl<'a> MovableMesh<'a> {
    pub(crate) fn new(mesh: &'a dyn MeshView) -> Self {
        let coords = (0..mesh.n_nodes() as u32).map(|n| mesh.coord(n)).collect();
        MovableMesh { mesh, coords: Arc::new(SharedCoords::new(coords)) }
    }

    /// Installs the coordinates as this thread's mesh motion until the
    /// guard drops.
    pub(crate) fn install(&self) -> MotionGuard {
        let coords: Arc<dyn MeshMotion> = self.coords.clone();
        MotionGuard(motion::install(Some(coords)))
    }
}

/// Restores the previously installed mesh motion on drop.
pub(crate) struct MotionGuard(Option<Arc<dyn MeshMotion>>);

impl Drop for MotionGuard {
    fn drop(&mut self) {
        motion::install(self.0.take());
    }
}

impl MeshView for MovableMesh<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords.get(node)
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
    fn element_faces(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_faces(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::ElementType;

    struct OneSegment;

    impl MeshView for OneSegment {
        fn dim(&self) -> usize {
            1
        }
        fn n_nodes(&self) -> usize {
            2
        }
        fn coord(&self, node: u32) -> [f64; 3] {
            [node as f64, 0.0, 0.0]
        }
        fn n_bulk(&self) -> usize {
            1
        }
        fn n_boundary(&self) -> usize {
            0
        }
        fn element(&self, _e: ElemId) -> ElementRef<'_> {
            ElementRef { etype: ElementType(202), nodes: &[0, 1], tag: 1, parents: [None, None] }
        }
    }

    /// A mapper's `set_coords` is what the mesh reads while the guard
    /// lives; dropping it uninstalls the motion.
    #[test]
    fn installed_motion_moves_the_nodes() {
        let mesh = MovableMesh::new(&OneSegment);
        assert_eq!(mesh.coord(1), [1.0, 0.0, 0.0]);
        {
            let _guard = mesh.install();
            let current = motion::install(None).unwrap();
            current.set_coords(&[[0.0; 3], [2.0, 0.5, 0.0]]);
            motion::install(Some(current));
            assert_eq!(mesh.coord(1), [2.0, 0.5, 0.0]);
        }
        assert!(motion::install(None).is_none());
    }
}
