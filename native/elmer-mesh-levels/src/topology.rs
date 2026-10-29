//! Global edge and face numbering of `FindMeshEdges` (MeshBasics.F90:3720).
//!
//! The Fortran keys a hash table on the smallest node(s), but the numbers it
//! hands out are just the order of first appearance while walking every
//! element (bulk, then boundary) and its local edges or faces. An edge keeps
//! the node order of the element that first produced it, and a face the
//! node order of its first element's face map; `SplitMeshEqual` sums the
//! face nodes in that order.

use std::collections::HashMap;

use crate::{Mesh, MeshLevelsError};

/// The element types the edge and face maps below cover.
pub const LINEAR: [u16; 7] = [101, 202, 303, 404, 504, 706, 808];

/// Local edges of the linear elements as 0-based node pairs
/// (`FindMeshEdges3D` edge maps; `FindMeshEdges2D` walks the same cycle).
fn edge_map(code: u16) -> Option<&'static [[usize; 2]]> {
    Some(match code / 100 {
        2 => &[[0, 1]],
        3 => &[[0, 1], [1, 2], [2, 0]],
        4 => &[[0, 1], [1, 2], [2, 3], [3, 0]],
        5 => &[[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]],
        7 => &[
            [0, 1],
            [1, 2],
            [0, 2],
            [3, 4],
            [4, 5],
            [5, 3],
            [0, 3],
            [1, 4],
            [2, 5],
        ],
        8 => &[
            [0, 1],
            [1, 2],
            [3, 2],
            [0, 3],
            [4, 5],
            [5, 6],
            [7, 6],
            [4, 7],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
        ],
        _ => return None,
    })
}

/// Local faces of the linear elements as 0-based node lists
/// (`FindMeshFaces3D` face maps).
fn face_map(code: u16) -> Option<&'static [&'static [usize]]> {
    Some(match code / 100 {
        3 => &[&[0, 1, 2]],
        4 => &[&[0, 1, 2, 3]],
        5 => &[&[0, 1, 2], &[0, 1, 3], &[1, 2, 3], &[2, 0, 3]],
        7 => &[
            &[0, 1, 2],
            &[3, 4, 5],
            &[0, 1, 4, 3],
            &[2, 1, 4, 5],
            &[2, 0, 3, 5],
        ],
        8 => &[
            &[0, 1, 2, 3],
            &[4, 5, 6, 7],
            &[0, 1, 5, 4],
            &[1, 2, 6, 5],
            &[2, 3, 7, 6],
            &[3, 0, 4, 7],
        ],
        _ => return None,
    })
}

/// Edges and faces of a mesh, with each element's local-to-global maps
/// (`Element % EdgeIndexes`, `Element % FaceIndexes`). Elements without
/// edges (101) or faces (2D meshes) get empty lists.
#[derive(Clone, Debug, Default)]
pub struct Topology {
    pub edges: Vec<[u32; 2]>,
    pub faces: Vec<Vec<u32>>,
    pub element_edges: Vec<Vec<u32>>,
    pub element_faces: Vec<Vec<u32>>,
}

/// `FindMeshEdges` without a mask: edges only when the bulk elements are at
/// most 2D (`FindMeshEdges2D`), faces then edges in 3D (`FindMeshFaces3D`,
/// `FindMeshEdges3D`). Only linear elements are accepted.
pub fn find_mesh_edges(mesh: &Mesh) -> Result<Topology, MeshLevelsError> {
    if let Some(el) = mesh.elements.iter().find(|e| !LINEAR.contains(&e.etype.0)) {
        return Err(MeshLevelsError::Unsupported(format!(
            "element type {} in mesh edge numbering",
            el.etype.0
        )));
    }
    let elem_dim = mesh
        .bulk()
        .iter()
        .map(|e| elem_dim(e.etype.0))
        .max()
        .unwrap_or(0);
    let n = mesh.elements.len();
    let mut topo = Topology {
        element_edges: vec![Vec::new(); n],
        element_faces: vec![Vec::new(); n],
        ..Default::default()
    };
    if elem_dim == 3 {
        let mut table: HashMap<[u32; 3], u32> = HashMap::new();
        for (e, el) in mesh.elements.iter().enumerate() {
            let Some(map) = face_map(el.etype.0) else {
                continue;
            };
            for f in map {
                let nodes: Vec<u32> = f.iter().map(|&k| el.nodes[k]).collect();
                let mut key = nodes.clone();
                key.sort_unstable();
                let id = *table.entry([key[0], key[1], key[2]]).or_insert_with(|| {
                    topo.faces.push(nodes);
                    topo.faces.len() as u32 - 1
                });
                topo.element_faces[e].push(id);
            }
        }
    }
    let mut table: HashMap<[u32; 2], u32> = HashMap::new();
    for (e, el) in mesh.elements.iter().enumerate() {
        let Some(map) = edge_map(el.etype.0) else {
            continue;
        };
        for &[a, b] in map {
            let (n1, n2) = (el.nodes[a], el.nodes[b]);
            let id = *table.entry([n1.min(n2), n1.max(n2)]).or_insert_with(|| {
                topo.edges.push([n1, n2]);
                topo.edges.len() as u32 - 1
            });
            topo.element_edges[e].push(id);
        }
    }
    Ok(topo)
}

/// `EnsureElemDim` (MeshBasics.F90:3800) of one element.
fn elem_dim(code: u16) -> usize {
    match code {
        c if c > 500 => 3,
        c if c > 300 => 2,
        c if c > 200 => 1,
        _ => 0,
    }
}
