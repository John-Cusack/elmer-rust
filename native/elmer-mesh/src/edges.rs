//! Global edge numbering, replicating `FindMeshEdges2D` and the linear
//! part of `FindMeshEdges3D`.
//!
//! Edges are numbered in first-use order over bulk then boundary
//! elements; each element records its edges in Elmer's canonical
//! per-type order (consecutive node pairs in 2D, the `TetraEdgeMap` /
//! `BrickEdgeMap` order in 3D). `IncreaseElementOrder` consumes this
//! table, so the order must match Elmer's or mid-edge nodes land in
//! the wrong slots.

use crate::elements::linear_edges;
use crate::mesh::{EdgeTable, Mesh};
use std::collections::HashMap;

impl Mesh {
    /// Build the global edge table and per-element edge indexes.
    /// 1D/2D use consecutive-pair enumeration; 3D uses the canonical
    /// family maps. Elements without a linear map keep empty indexes.
    pub fn find_edges(&mut self) {
        // (min_node, max_node) -> edge id, in first-use order.
        let mut map: HashMap<(u32, u32), u32> = HashMap::new();
        let mut edges: Vec<(u32, u32)> = Vec::new();
        let n_elem = self.bulk.len() + self.boundary.len();

        // Per-element edge endpoint pairs first (borrows end before insert).
        let mut elem_pairs: Vec<Vec<(u32, u32)>> = Vec::with_capacity(n_elem);
        let push_elem = |code: u16, nodes: &[u32], pairs: &mut Vec<Vec<(u32, u32)>>| {
            let mut p = Vec::new();
            if let Some(def) = linear_edges(code) {
                for &(a, b) in def {
                    if a < nodes.len() && b < nodes.len() {
                        p.push((nodes[a], nodes[b]));
                    }
                }
            }
            pairs.push(p);
        };
        for e in &self.bulk {
            push_elem(e.code, &e.nodes, &mut elem_pairs);
        }
        for e in &self.boundary {
            push_elem(e.code, &e.nodes, &mut elem_pairs);
        }

        let mut elem_edge_ids: Vec<Vec<u32>> = Vec::with_capacity(n_elem);
        for pairs in &elem_pairs {
            let mut ids = Vec::with_capacity(pairs.len());
            for &(n1, n2) in pairs {
                let key = if n1 < n2 { (n1, n2) } else { (n2, n1) };
                let id = *map.entry(key).or_insert_with(|| {
                    let id = edges.len() as u32;
                    edges.push(key);
                    id
                });
                ids.push(id);
            }
            elem_edge_ids.push(ids);
        }

        let n_bulk = self.bulk.len();
        for (i, ids) in elem_edge_ids.into_iter().enumerate() {
            if i < n_bulk {
                self.bulk[i].edges = ids;
            } else {
                self.boundary[i - n_bulk].edges = ids;
            }
        }
        self.edges = Some(EdgeTable { edges });
    }

    pub fn n_edges(&self) -> usize {
        self.edges.as_ref().map(|e| e.edges.len()).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use crate::mesh::{BoundElem, Elem, Mesh};

    fn two_tri_mesh() -> Mesh {
        // Two triangles sharing edge (1,2): nodes 0..3.
        Mesh {
            dim: 2,
            max_dim: 2,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            node_tags: vec![1, 2, 3, 4],
            bulk: vec![
                Elem { code: 303, body: 1, nodes: vec![0, 1, 2], tag: 1, edges: vec![] },
                Elem { code: 303, body: 1, nodes: vec![0, 2, 3], tag: 2, edges: vec![] },
            ],
            boundary: vec![BoundElem {
                code: 202,
                constraint: 1,
                parents: [Some(0), None],
                nodes: vec![0, 1],
                tag: 1,
                edges: vec![],
            }],
            edges: None,
        }
    }

    #[test]
    fn shared_edge_numbered_once_in_canonical_order() {
        let mut m = two_tri_mesh();
        m.find_edges();
        // Triangle edges in (n,n+1) order; shared edge (2,0)~(0,2)
        // gets one id regardless of orientation.
        assert_eq!(m.bulk[0].edges.len(), 3);
        assert_eq!(m.bulk[1].edges.len(), 3);
        assert_eq!(m.bulk[0].edges[2], m.bulk[1].edges[0]);
        // Boundary 202 on nodes (0,1) shares the first triangle edge.
        assert_eq!(m.boundary[0].edges, vec![m.bulk[0].edges[0]]);
        assert_eq!(m.n_edges(), 5);
    }

    #[test]
    fn increase_order_midpoints_and_types() {
        let mut m = two_tri_mesh();
        m.find_edges();
        let eid = m.bulk[0].edges[0];
        crate::transform::increase_element_order(&mut m).unwrap();
        assert_eq!(m.nodes.len(), 4 + 5);
        // Mid-edge node of (0,1) is at (0.5, 0, 0) and stored after corners.
        let mid = m.nodes[4 + eid as usize];
        assert!((mid[0] - 0.5).abs() < 1e-15 && mid[1] == 0.0);
        assert!(m.bulk.iter().all(|e| e.code == 306 && e.nodes.len() == 6));
        assert_eq!(m.boundary[0].code, 203);
    }
}
