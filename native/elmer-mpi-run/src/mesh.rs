//! [`ParMesh`]: one rank's piece of an ElmerGrid partitioning as the
//! [`MeshView`] the solver assembles on, with the global id and owner of
//! every local node.
//!
//! Local nodes keep their partition-file order, as Elmer reads them, so the
//! bandwidth optimisation gives Elmer's per-partition `Perm`. The first rank
//! in a node's `part.K.shared` list owns it ([`ParMesh::owned_nodes`]).

use std::path::Path;

use elmer_core::{ElemId, ElementRef, MeshView};
use elmer_dist::Exchange;
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
    /// The last `n_imported` elements came from other partitions
    /// ([`ParMesh::radiation_distribute`]); they are not halo elements.
    n_imported: usize,
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
            n_imported: 0,
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
        let local = self.elem_owner.len() - self.n_imported;
        let h: Vec<bool> = self.elem_owner.iter().enumerate().map(|(e, &o)| e < local && o != self.rank).collect();
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

    /// `RadiationParallelMeshDistribute` (MeshBasics.F90:1411-1771): every
    /// radiating boundary element of the other partitions, with its nodes,
    /// is appended to this one, so that each rank sees the whole radiation
    /// surface. `radiating(tag)` is `RadiationCheck` of a boundary element's
    /// condition. Collective.
    ///
    /// * The nodes of this rank's radiating elements become shared with
    ///   every rank (`NeighbourList` = this rank, then the others).
    /// * A received node already here gains the sender and this rank among
    ///   its holders; a new one is appended, held by the sender and this
    ///   rank. A received element keeps the sender as owner (`PartIndex`)
    ///   and its global id, without parents.
    /// * Each radiating node is then owned by the owner of the radiating
    ///   element with the largest global id touching it.
    ///
    /// Returns the imported elements' flags (bulk then boundary).
    pub fn radiation_distribute(&mut self, ex: &dyn Exchange, radiating: &dyn Fn(i32) -> bool) -> Vec<bool> {
        let (np, me) = (self.n_parts, self.rank);
        let mut imported = vec![false; self.elems.len()];
        if np <= 1 || ex.size() <= 1 {
            return imported;
        }
        let nb = self.n_bulk;
        let n_nodes = self.coords.len();
        let own: Vec<usize> = (nb..self.elems.len()).filter(|&e| radiating(self.elems[e].tag)).collect();
        let mut flag = vec![false; n_nodes];
        for &e in &own {
            for &n in &self.elems[e].nodes {
                flag[n as usize] = true;
            }
        }
        // cPerm and the coordinates sent; the own radiating nodes are
        // shared with every rank.
        let mut cperm = vec![0u64; n_nodes];
        let mut coords: Vec<u64> = Vec::new();
        let mut n_coord = 0u64;
        for i in 0..n_nodes {
            if !flag[i] {
                continue;
            }
            n_coord += 1;
            cperm[i] = n_coord;
            let c = self.coords[i];
            coords.extend([self.node_gid[i], c[0].to_bits(), c[1].to_bits(), c[2].to_bits()]);
            self.node_parts[i] = std::iter::once(me).chain((0..ex.size()).filter(|&j| j != me)).collect();
        }
        let mut msg: Vec<u64> = vec![own.len() as u64, n_coord];
        for &e in &own {
            let el = &self.elems[e];
            msg.extend([u64::from(el.etype.0), el.tag as i64 as u64, self.elem_gid[e], el.nodes.len() as u64]);
            msg.extend(el.nodes.iter().map(|&n| cperm[n as usize]));
            self.elem_owner[e] = me;
        }
        msg.extend(coords);
        let send: Vec<Vec<u64>> = (0..ex.size()).map(|j| if j == me { Vec::new() } else { msg.clone() }).collect();
        let recv = ex.all_to_all_u64(&send);
        let total: u64 = recv.iter().enumerate().filter(|&(i, _)| i != me).map(|(_, r)| r.first().copied().unwrap_or(0)).sum();
        if own.is_empty() || total == 0 {
            self.node_owner = self.node_parts.iter().map(|p| p[0]).collect();
            return imported;
        }
        for (i, r) in recv.iter().enumerate() {
            if i == me || r.first().copied().unwrap_or(0) == 0 {
                continue;
            }
            let (n_surf, n_c) = (r[0] as usize, r[1] as usize);
            let mut pos = 2;
            let mut surfs = Vec::with_capacity(n_surf);
            for _ in 0..n_surf {
                let (code, tag, gid, nn) = (r[pos] as u16, r[pos + 1] as i64 as i32, r[pos + 2], r[pos + 3] as usize);
                surfs.push((code, tag, gid, r[pos + 4..pos + 4 + nn].to_vec()));
                pos += 4 + nn;
            }
            // SearchNode over the nodes present before this sender's.
            let mut by_gid: Vec<(u64, usize)> = self.node_gid.iter().copied().zip(0..).collect();
            by_gid.sort_unstable();
            let mut local = vec![0usize; n_c + 1];
            for j in 0..n_c {
                let q = &r[pos + 4 * j..pos + 4 * j + 4];
                let k = by_gid.binary_search_by_key(&q[0], |&(g, _)| g).ok().map(|x| by_gid[x].1);
                local[j + 1] = match k {
                    Some(k) => {
                        let parts = &mut self.node_parts[k];
                        if !parts.contains(&i) {
                            parts.push(i);
                        }
                        if !parts.contains(&me) {
                            parts.push(me);
                        }
                        k
                    }
                    None => {
                        self.coords.push([f64::from_bits(q[1]), f64::from_bits(q[2]), f64::from_bits(q[3])]);
                        self.node_gid.push(q[0]);
                        self.node_parts.push(vec![i, me]);
                        self.coords.len() - 1
                    }
                };
            }
            for (code, tag, gid, nodes) in surfs {
                self.elems.push(Elem {
                    etype: elmer_core::ElementType(code),
                    nodes: nodes.iter().map(|&c| local[c as usize] as u32).collect(),
                    tag,
                    parents: [None, None],
                });
                self.elem_gid.push(gid);
                self.elem_owner.push(i);
                imported.push(true);
                self.n_imported += 1;
            }
        }
        // The owner of a radiating node: the owner of the radiating element
        // of largest global id on it.
        let mut gmax = vec![0u64; self.coords.len()];
        let rad: Vec<usize> = (nb..self.elems.len()).filter(|&e| radiating(self.elems[e].tag)).collect();
        for &e in &rad {
            for &n in &self.elems[e].nodes {
                gmax[n as usize] = gmax[n as usize].max(self.elem_gid[e]);
            }
        }
        for &e in &rad {
            let owner = self.elem_owner[e];
            for &n in &self.elems[e].nodes {
                if self.elem_gid[e] != gmax[n as usize] {
                    continue;
                }
                let parts = &mut self.node_parts[n as usize];
                if let Some(k) = parts.iter().position(|&p| p == owner) {
                    parts.swap(0, k);
                }
            }
        }
        self.node_owner = self.node_parts.iter().map(|p| p[0]).collect();
        imported
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
