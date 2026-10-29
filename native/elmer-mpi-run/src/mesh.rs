//! [`ParMesh`]: one rank's piece of an ElmerGrid partitioning as the
//! [`MeshView`] the solver assembles on, with the global id and owner of
//! every local node.
//!
//! Local nodes are renumbered so that the nodes another rank owns come
//! first and the owned ones last, each group in partition-file order (the
//! first rank in a node's `part.K.shared` list owns it). The
//! owned nodes are then the tail `n_ghost..` of the node range, which is
//! what lets [`crate::backend`] express Elmer's consistent parallel norm
//! (`Nonlinear System Consistent Norm`) through `elmer-solver`'s `Norm
//! Permutation`.

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
    n_ghost: usize,
    elems: Vec<Elem>,
    n_bulk: usize,
}

/// Node data in any order, before [`ParMesh::new`] puts the ghosts first.
pub(crate) struct Nodes {
    pub(crate) coords: Vec<[f64; 3]>,
    pub(crate) gid: Vec<u64>,
    pub(crate) parts: Vec<Vec<usize>>,
}

impl ParMesh {
    /// Renumbers `nodes` ghosts first, owned last (each group in the given
    /// order) and the element node lists with them.
    pub(crate) fn new(rank: usize, n_parts: usize, dim: usize, nodes: Nodes, mut elems: Vec<Elem>, n_bulk: usize) -> Self {
        let n = nodes.gid.len();
        let owner = |i: usize| nodes.parts[i][0];
        let mut order: Vec<usize> = (0..n).filter(|&i| owner(i) != rank).collect();
        let n_ghost = order.len();
        order.extend((0..n).filter(|&i| owner(i) == rank));
        let mut new_of = vec![0u32; n];
        for (k, &i) in order.iter().enumerate() {
            new_of[i] = k as u32;
        }
        for e in &mut elems {
            for v in &mut e.nodes {
                *v = new_of[*v as usize];
            }
        }
        ParMesh {
            rank,
            n_parts,
            dim,
            coords: order.iter().map(|&i| nodes.coords[i]).collect(),
            node_gid: order.iter().map(|&i| nodes.gid[i]).collect(),
            node_parts: order.iter().map(|&i| nodes.parts[i].clone()).collect(),
            node_owner: order.iter().map(|&i| owner(i)).collect(),
            n_ghost,
            elems,
            n_bulk,
        }
    }

    /// Rank `rank` of `<mesh_dir>/partitioning.<n_parts>`, coordinates
    /// multiplied by `scale` (`Coordinate Scaling`). Partitions written with
    /// halo elements (`ElmerGrid -halo`) are an error: their ghost elements
    /// would be assembled twice.
    pub fn read(mesh_dir: &Path, rank: usize, n_parts: usize, scale: f64) -> Result<Self, String> {
        let part = PartitionMesh::read(mesh_dir, rank, n_parts).map_err(|e: PartError| e.to_string())?;
        if part.n_owned_bulk() != part.bulk_elements().len() {
            return Err("partitions with halo elements (ElmerGrid -halo)".into());
        }
        Ok(Self::from_partition(&part, scale))
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
        Self::new(rank, part.n_parts(), part.dim(), nodes, elems_of(part), part.n_bulk())
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

    /// Nodes `0..n_ghost` are owned by other ranks, the rest by this one.
    pub fn n_ghost(&self) -> usize {
        self.n_ghost
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
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/angle/angle")
    }

    #[test]
    fn owned_nodes_are_numbered_last() {
        for rank in 0..3 {
            let part = PartitionMesh::read(&fixture(), rank, 3).unwrap();
            let mesh = ParMesh::from_partition(&part, 1.0);
            let k = mesh.n_ghost();
            assert_eq!(mesh.n_nodes() - k, part.n_owned_nodes());
            assert!(mesh.node_owners()[..k].iter().all(|&o| o != rank));
            assert!(mesh.node_owners()[k..].iter().all(|&o| o == rank));
            // Same elements over the same global nodes and coordinates.
            for e in 0..part.n_bulk() + part.n_boundary() {
                let (a, b) = (part.element(ElemId(e as u32)), mesh.element(ElemId(e as u32)));
                assert_eq!((a.etype, a.tag, a.parents), (b.etype, b.tag, b.parents));
                for (&i, &j) in a.nodes.iter().zip(b.nodes) {
                    assert_eq!(part.node_gids()[i as usize], mesh.node_gids()[j as usize]);
                    assert_eq!(part.coord(i), mesh.coord(j));
                }
            }
        }
    }
}
