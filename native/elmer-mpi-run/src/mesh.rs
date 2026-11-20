//! [`ParMesh`]: one rank's piece of an ElmerGrid partitioning as the
//! [`MeshView`] the solver assembles on, with the global id and owner of
//! every local node.
//!
//! Local nodes keep their partition-file order, as Elmer reads them, so the
//! bandwidth optimisation gives Elmer's per-partition `Perm`. The first rank
//! in a node's `part.K.shared` list owns it ([`ParMesh::owned_nodes`]).

use std::path::Path;

use elmer_core::{ElemId, ElementRef, MeshView};
use elmer_mpi::{PartError, PartitionMesh};

pub(crate) struct Elem {
    pub(crate) etype: elmer_core::ElementType,
    pub(crate) nodes: Vec<u32>,
    pub(crate) tag: i32,
    pub(crate) parents: [Option<ElemId>; 2],
}

/// One rank's mesh: owned bulk elements, then the boundary elements.
pub struct ParMesh {
    rank: usize,
    n_parts: usize,
    dim: usize,
    coords: Vec<[f64; 3]>,
    node_gid: Vec<u64>,
    /// Ranks holding each node; the first owns it.
    node_parts: Vec<Vec<usize>>,
    node_owner: Vec<usize>,
    elems: Vec<Elem>,
    n_bulk: usize,
    /// Global id (`GElementIndex`) and owning partition (`PartIndex`) of
    /// every element, bulk then boundary; empty for a serial mesh.
    elem_gid: Vec<u64>,
    elem_owner: Vec<usize>,
    /// [`MeshView::edge_order`] ([`crate::edges`]).
    edge_keys: Option<Vec<u64>>,
}

/// Node data in local order.
pub(crate) struct Nodes {
    pub(crate) coords: Vec<[f64; 3]>,
    pub(crate) gid: Vec<u64>,
    pub(crate) parts: Vec<Vec<usize>>,
}

impl ParMesh {
    /// One rank's mesh over `nodes` in their given order.
    pub(crate) fn new(rank: usize, n_parts: usize, dim: usize, nodes: Nodes, elems: Vec<Elem>, n_bulk: usize) -> Self {
        let node_owner = nodes.parts.iter().map(|p| p[0]).collect();
        ParMesh {
            rank,
            n_parts,
            dim,
            coords: nodes.coords,
            node_gid: nodes.gid,
            node_parts: nodes.parts,
            node_owner,
            elems,
            n_bulk,
            elem_gid: Vec::new(),
            edge_keys: None,
            elem_owner: Vec::new(),
        }
    }

    /// Rank `rank` of `<mesh_dir>/partitioning.<n_parts>`, coordinates
    /// multiplied by `scale` (`Coordinate Scaling`). Halo elements
    /// (`ElmerGrid -halo`) are kept, marked by [`ParMesh::halo_elements`]:
    /// no solver assembles them, a DG solver numbers their DG nodes.
    pub fn read(mesh_dir: &Path, rank: usize, n_parts: usize, scale: f64) -> Result<Self, String> {
        let part = PartitionMesh::read(mesh_dir, rank, n_parts).map_err(|e: PartError| e.to_string())?;
        Ok(Self::from_partition(&part, scale))
    }

    /// The elements (bulk then boundary) of another partition, `None` when
    /// every element is this rank's.
    pub fn halo_elements(&self) -> Option<Vec<bool>> {
        let h: Vec<bool> = self.elem_owner.iter().map(|&o| o != self.rank).collect();
        h.iter().any(|&x| x).then_some(h)
    }

    /// Sets [`MeshView::edge_order`].
    pub fn set_edge_order(&mut self, keys: Vec<u64>) {
        self.edge_keys = Some(keys);
    }

    /// Global element ids (`GElementIndex`), bulk then boundary.
    pub fn element_gids(&self) -> &[u64] {
        &self.elem_gid
    }

    /// Owning partition of every element, bulk then boundary.
    pub fn element_owners(&self) -> &[usize] {
        &self.elem_owner
    }

    /// The whole mesh DB in `mesh_dir` as the only rank's piece, as
    /// `ElmerSolver_mpi` reads it at np=1: every node owned, global id
    /// `index + 1`.
    pub fn read_serial(mesh_dir: &Path, scale: f64) -> Result<Self, String> {
        let dir = mesh_dir.to_str().ok_or("non-UTF-8 mesh path")?;
        let mesh = elmer_mesh::load(dir).map_err(|e| format!("{dir}: {e:?}"))?;
        Ok(Self::from_view(&elmer_mesh::adapters::MeshViewAdapter::new(&mesh), scale))
    }

    /// One rank holding all of `view`.
    pub fn from_view(view: &dyn MeshView, scale: f64) -> Self {
        let n = view.n_nodes();
        let nodes = Nodes {
            coords: (0..n).map(|i| view.coord(i as u32).map(|c| scale * c)).collect(),
            gid: (1..=n as u64).collect(),
            parts: vec![vec![0]; n],
        };
        Self::new(0, 1, view.dim(), nodes, elems_of(view), view.n_bulk())
    }

    /// [`ParMesh::read`] from an already read partition without halo
    /// elements.
    pub fn from_partition(part: &PartitionMesh, scale: f64) -> Self {
        let rank = part.rank();
        let n = part.n_nodes();
        let mut parts = vec![vec![rank]; n];
        for s in part.shared() {
            parts[s.local as usize] = s.parts.clone();
        }
        let nodes = Nodes {
            coords: (0..n).map(|i| part.coord(i as u32).map(|c| scale * c)).collect(),
            gid: part.node_gids().to_vec(),
            parts,
        };
        let mut m = Self::new(rank, part.n_parts(), part.dim(), nodes, elems_of(part), part.n_bulk());
        let all = part.bulk_elements().iter().chain(part.boundary_elements());
        (m.elem_gid, m.elem_owner) = all.map(|e| (e.gid, e.owner)).unzip();
        m
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn n_parts(&self) -> usize {
        self.n_parts
    }

    /// Global (1-based, as in the mesh files) id of every local node.
    pub fn node_gids(&self) -> &[u64] {
        &self.node_gid
    }

    /// Owning rank of every local node.
    pub fn node_owners(&self) -> &[usize] {
        &self.node_owner
    }

    /// Ranks holding each local node (this one included); the first owns it.
    pub fn node_parts(&self) -> &[Vec<usize>] {
        &self.node_parts
    }

    /// Whether this rank owns each local node.
    pub fn owned_nodes(&self) -> Vec<bool> {
        self.node_owner.iter().map(|&o| o == self.rank).collect()
    }

    /// Bulk and boundary element types present.
    pub fn element_types(&self) -> std::collections::BTreeSet<u16> {
        self.elems.iter().map(|e| e.etype.0).collect()
    }
}

pub(crate) fn elems_of(view: &dyn MeshView) -> Vec<Elem> {
    (0..view.n_bulk() + view.n_boundary())
        .map(|e| {
            let el = view.element(ElemId(e as u32));
            Elem { etype: el.etype, nodes: el.nodes.to_vec(), tag: el.tag, parents: el.parents }
        })
        .collect()
}

impl MeshView for ParMesh {
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
        ElementRef { etype: el.etype, nodes: &el.nodes, tag: el.tag, parents: el.parents }
    }

    fn edge_order(&self) -> Option<&[u64]> {
        self.edge_keys.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/angle/angle")
    }

    #[test]
    fn nodes_keep_the_partition_file_order() {
        for rank in 0..3 {
            let part = PartitionMesh::read(&fixture(), rank, 3).unwrap();
            let mesh = ParMesh::from_partition(&part, 1.0);
            let owned = mesh.owned_nodes();
            assert_eq!(owned.iter().filter(|&&o| o).count(), part.n_owned_nodes());
            assert_eq!(mesh.node_gids(), part.node_gids());
            for (i, &o) in owned.iter().enumerate() {
                assert_eq!(o, mesh.node_owners()[i] == rank);
                assert_eq!(part.coord(i as u32), mesh.coord(i as u32));
            }
            // The same elements over the same nodes.
            for e in 0..part.n_bulk() + part.n_boundary() {
                let (a, b) = (part.element(ElemId(e as u32)), mesh.element(ElemId(e as u32)));
                assert_eq!((a.etype, a.tag, a.parents, a.nodes), (b.etype, b.tag, b.parents, b.nodes));
            }
        }
    }
}
