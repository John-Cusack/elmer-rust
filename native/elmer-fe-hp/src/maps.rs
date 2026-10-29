//! Local edge/face numbering, global-parity directions and p-DOF counts.
//!
//! Transliterates `H1Basis_GetEdgeMap`, `H1Basis_GetFaceMap`,
//! `H1Basis_GetEdgeDirection`, `H1Basis_GetFaceDirection`,
//! `H1Basis_GetTetra{Edge,Face}Direction` (H1Basis.F90), the mapping tables of
//! `PElementMaps.F90` (`InitializeMappings`) and its DOF-count functions
//! `getEdgeDOFs`, `getFaceDOFs`, `getBubbleDOFs`, `getBubbleMaxDegree1D`.
//!
//! Local node numbers are 1-based as in the Fortran.

#![allow(clippy::chunks_exact_to_as_chunks)]

use crate::types::{EdgeDir, FaceDir};

/// `H1Basis_GetEdgeMap`: local node pairs of each edge for element code
/// `ecode`, or `None` when the Fortran reports "Not fully implemented yet".
pub fn edge_map(ecode: u16) -> Option<EdgeDir> {
    let flat: &[usize] = match ecode {
        202 => &[1, 2],
        303 => &[1, 2, 2, 3, 3, 1],
        404 => &[1, 2, 2, 3, 4, 3, 1, 4],
        504 => return Some(tetra_edge_direction(1)),
        605 => &[1, 2, 2, 3, 4, 3, 1, 4, 1, 5, 2, 5, 3, 5, 4, 5],
        706 => &[1, 2, 2, 3, 3, 1, 4, 5, 5, 6, 6, 4, 1, 4, 2, 5, 3, 6],
        808 => &[
            1, 2, 2, 3, 4, 3, 1, 4, 5, 6, 6, 7, 8, 7, 5, 8, 1, 5, 2, 6, 3, 7, 4, 8,
        ],
        _ => return None,
    };
    let mut map = [[0usize; 2]; 12];
    for (i, e) in flat.chunks_exact(2).enumerate() {
        map[i] = [e[0], e[1]];
    }
    Some(map)
}

/// `H1Basis_GetFaceMap`: local nodes of each face (`0` in slot 4 for a
/// triangle) for element code `ecode`.
pub fn face_map(ecode: u16) -> Option<FaceDir> {
    let flat: &[usize] = match ecode {
        303 => &[1, 2, 3, 0],
        404 => &[1, 2, 3, 4],
        504 => return Some(tetra_face_direction(1)),
        605 => &[
            1, 2, 3, 4, //
            1, 2, 5, 0, //
            2, 3, 5, 0, //
            3, 4, 5, 0, //
            4, 1, 5, 0,
        ],
        706 => &[
            1, 2, 3, 0, //
            4, 5, 6, 0, //
            1, 2, 5, 4, //
            2, 3, 6, 5, //
            3, 1, 4, 6,
        ],
        808 => &[
            1, 2, 3, 4, //
            5, 6, 7, 8, //
            1, 2, 6, 5, //
            2, 3, 7, 6, //
            4, 3, 7, 8, //
            1, 4, 8, 5,
        ],
        _ => return None,
    };
    let mut map = [[0usize; 4]; 6];
    for (i, f) in flat.chunks_exact(4).enumerate() {
        map[i] = [f[0], f[1], f[2], f[3]];
    }
    Some(map)
}

/// `H1Basis_GetTetraEdgeDirection`. Panics on an unknown tetra type (the
/// Fortran calls `Fatal`).
pub fn tetra_edge_direction(ttype: i32) -> EdgeDir {
    let flat: [usize; 12] = match ttype {
        1 => [1, 2, 2, 3, 1, 3, 1, 4, 2, 4, 3, 4],
        2 => [1, 2, 3, 2, 1, 3, 1, 4, 2, 4, 3, 4],
        _ => panic!("tetra_edge_direction: unknown tetra type {ttype}"),
    };
    let mut map = [[0usize; 2]; 12];
    for (i, e) in flat.chunks_exact(2).enumerate() {
        map[i] = [e[0], e[1]];
    }
    map
}

/// `H1Basis_GetTetraFaceDirection`. Panics on an unknown tetra type.
pub fn tetra_face_direction(ttype: i32) -> FaceDir {
    let flat: [usize; 16] = match ttype {
        1 => [1, 2, 3, 0, 1, 2, 4, 0, 2, 3, 4, 0, 1, 3, 4, 0],
        2 => [1, 3, 2, 0, 1, 2, 4, 0, 3, 2, 4, 0, 1, 3, 4, 0],
        _ => panic!("tetra_face_direction: unknown tetra type {ttype}"),
    };
    let mut map = [[0usize; 4]; 6];
    for (i, f) in flat.chunks_exact(4).enumerate() {
        map[i] = [f[0], f[1], f[2], f[3]];
    }
    map
}

/// `H1Basis_GetEdgeDirection`: the edge map with each pair swapped where
/// needed so the global node index increases along the edge. `globalind` is
/// indexed by 0-based local node (`globalind[k-1]` is the index of local node
/// `k`). Only the first `nedges` entries are touched, like the Fortran.
pub fn edge_direction(ecode: u16, nedges: usize, globalind: &[i64]) -> Option<EdgeDir> {
    let mut direction = edge_map(ecode)?;
    for d in direction.iter_mut().take(nedges) {
        if globalind[d[1] - 1] < globalind[d[0] - 1] {
            d.swap(0, 1);
        }
    }
    Some(direction)
}

/// `H1Basis_GetFaceDirection`: the face map re-ordered so the local nodes of
/// each face follow ascending global index (triangles) or start at the
/// smallest global index and continue towards its smaller neighbour (quads).
pub fn face_direction(ecode: u16, nfaces: usize, globalind: &[i64]) -> Option<FaceDir> {
    let mut direction = face_map(ecode)?;
    let g = |n: usize| globalind[n - 1];
    for d in direction.iter_mut().take(nfaces) {
        if d[3] == 0 {
            // Triangle face: sort the three local nodes by global index.
            if g(d[0]) > g(d[1]) {
                d.swap(0, 1);
            }
            if g(d[0]) > g(d[2]) {
                d.swap(0, 2);
            }
            if g(d[1]) > g(d[2]) {
                d.swap(1, 2);
            }
        } else {
            // Square face: rotate the smallest global index to the front.
            let mut min_i = 0;
            for i in 1..4 {
                if g(d[i]) < g(d[min_i]) {
                    min_i = i;
                }
            }
            match min_i {
                0 => {}
                1 => {
                    let tmp = d[0];
                    for i in 0..3 {
                        d[i] = d[i + 1];
                    }
                    d[3] = tmp;
                }
                2 => {
                    d.swap(0, 2);
                    d.swap(1, 3);
                }
                _ => {
                    let tmp = d[3];
                    for i in (0..3).rev() {
                        d[i + 1] = d[i];
                    }
                    d[0] = tmp;
                }
            }
            // The smaller of the two neighbours of the minimum goes second;
            // [A,B,D,C] is then a globally consistent ordering.
            if g(d[1]) > g(d[3]) {
                d.swap(1, 3);
            }
        }
    }
    Some(direction)
}

/// Edge map tables of `PElementMaps::InitializeMappings` (`getXxxEdgeMap`),
/// by 1-based edge number. Tetra tables are type 1.
pub fn pelement_edge_map(ecode: u16, tetra_type: i32) -> Option<EdgeDir> {
    if ecode == 504 {
        return Some(tetra_edge_direction(tetra_type));
    }
    edge_map(ecode)
}

/// `TetraFaceEdgeMap{1,2}` / `BrickFaceEdgeMap` / `WedgeFaceEdgeMap` /
/// `PyramidFaceEdgeMap` of `PElementMaps`: the (1-based) edges bounding each
/// face, `0` padding the unused slots.
pub fn face_edge_map(ecode: u16, tetra_type: i32) -> Option<[[usize; 4]; 6]> {
    let rows: &[[usize; 4]] = match ecode {
        504 => match tetra_type {
            1 => &[[1, 2, 3, 0], [1, 5, 4, 0], [2, 6, 5, 0], [3, 6, 4, 0]],
            2 => &[[3, 2, 1, 0], [1, 5, 4, 0], [2, 5, 6, 0], [3, 6, 4, 0]],
            _ => return None,
        },
        808 => &[
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [1, 10, 5, 9],
            [2, 11, 6, 10],
            [3, 11, 7, 12],
            [4, 12, 8, 9],
        ],
        706 => &[
            [1, 2, 3, 0],
            [4, 5, 6, 0],
            [1, 8, 4, 7],
            [2, 9, 5, 8],
            [3, 7, 6, 9],
        ],
        605 => &[
            [1, 2, 3, 4],
            [1, 6, 5, 0],
            [2, 7, 6, 0],
            [3, 8, 7, 0],
            [4, 5, 8, 0],
        ],
        _ => return None,
    };
    let mut map = [[0usize; 4]; 6];
    map[..rows.len()].copy_from_slice(rows);
    Some(map)
}

/// `getEdgeDOFs(Element, p)` for an element that has p definitions.
pub fn edge_dofs(p: i32) -> usize {
    (p - 1).max(0) as usize
}

/// `getFaceDOFs(Element, p, faceNumber)`; `face_number` is 1-based and only
/// matters for pyramids and wedges. `serendipity` is `PDefs % Serendipity`.
/// Unsupported types return `p` like the Fortran's default branch.
pub fn face_dofs(ecode: u16, p: i32, face_number: usize, serendipity: bool) -> usize {
    if p <= 1 {
        return 0;
    }
    let quad = |sd: bool| {
        if sd {
            (p - 2) * (p - 3) / 2
        } else {
            (p - 1) * (p - 1)
        }
    };
    let tri = (p - 1) * (p - 2) / 2;
    let n = match ecode / 100 {
        3 | 5 => tri,
        4 | 8 => quad(serendipity),
        6 => match face_number {
            1 => (p - 1) * (p - 1),
            2..=5 => tri,
            _ => 0,
        },
        7 => match face_number {
            1 | 2 => tri,
            3..=5 => quad(serendipity),
            _ => 0,
        },
        _ => p,
    };
    n.max(0) as usize
}

/// `getBubbleDOFs(Element, p)`.
pub fn bubble_dofs(ecode: u16, p: i32, serendipity: bool) -> usize {
    if p <= 1 {
        return 0;
    }
    let n = match ecode / 100 {
        2 => p - 1,
        3 => (p - 2) * (p - 1) / 2,
        4 => {
            if serendipity {
                (p - 2) * (p - 3) / 2
            } else {
                (p - 1) * (p - 1)
            }
        }
        5 => (p - 3) * (p - 2) * (p - 1) / 6,
        6 => (p - 2) * (p - 1) * p / 6,
        7 => {
            if serendipity {
                (p - 2) * (p - 3) * (p - 4) / 6
            } else {
                (p - 1) * (p - 1) * (p - 2) / 2
            }
        }
        8 => {
            if serendipity {
                (p - 3) * (p - 4) * (p - 5) / 6
            } else {
                (p - 1) * (p - 1) * (p - 1)
            }
        }
        _ => p,
    };
    n.max(0) as usize
}

/// `getBubbleMaxDegree1D(Element, p)`.
pub fn bubble_max_degree_1d(ecode: u16, p: i32, serendipity: bool) -> i32 {
    let mut deg = p;
    if serendipity {
        match ecode / 100 {
            4 | 7 => deg = p - 2,
            8 => deg = p - 4,
            _ => {}
        }
    }
    deg.max(0)
}

/// Number of nodal functions, edges and faces per element code.
fn topology(ecode: u16) -> Option<(usize, usize, usize)> {
    Some(match ecode {
        // A line's p functions are its "bubbles"; it has no separate edge set.
        202 => (2, 0, 0),
        303 => (3, 3, 1),
        404 => (4, 4, 1),
        504 => (4, 6, 4),
        605 => (5, 8, 5),
        706 => (6, 9, 5),
        808 => (8, 12, 6),
        _ => return None,
    })
}

/// `(nodal, edge, face, bubble)` basis-function counts for polynomial degree
/// `p` with the full (non-serendipity) space: `edge` is per edge and `face`
/// per (first) face, the others as in `getEdgeDOFs`, `getFaceDOFs` and
/// `getBubbleDOFs`. For wedges and pyramids use [`face_dofs`] to get the
/// count of each particular face.
pub fn p_dof_counts(ecode: u16, p: u8) -> Option<(usize, usize, usize, usize)> {
    let (nodes, _, _) = topology(ecode)?;
    let p = i32::from(p);
    Some((
        nodes,
        edge_dofs(p),
        face_dofs(ecode, p, 1, false),
        bubble_dofs(ecode, p, false),
    ))
}

/// Total number of basis functions of the element at uniform degree `p`
/// (`ElementInfoVec` counts: nodes + all edge, face and bubble functions).
pub fn p_total_dofs(ecode: u16, p: u8, serendipity: bool) -> Option<usize> {
    let (nodes, nedges, nfaces) = topology(ecode)?;
    let p = i32::from(p);
    let mut n = nodes + nedges * edge_dofs(p);
    for f in 1..=nfaces {
        n += face_dofs(ecode, p, f, serendipity);
    }
    // A 2D element's single "face" is its interior, already counted as bubble.
    if ecode / 100 == 3 || ecode / 100 == 4 {
        n -= face_dofs(ecode, p, 1, serendipity);
    }
    Some(n + bubble_dofs(ecode, p, serendipity))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brick_edge_direction_follows_global_index() {
        // Global indices decreasing with local number: every edge must come
        // out with ascending global index from first to second node.
        let g: Vec<i64> = (1..=8).rev().collect();
        let d = edge_direction(808, 12, &g).unwrap();
        for e in d.iter() {
            assert!(g[e[0] - 1] < g[e[1] - 1]);
        }
        // Same numbering as the map: nothing swaps.
        let id: Vec<i64> = (1..=8).collect();
        let d = edge_direction(808, 12, &id).unwrap();
        assert_eq!(d[2], [3, 4]);
    }

    #[test]
    fn quad_face_direction_rotates_smallest_first() {
        // Local node 3 holds the smallest global index.
        let g = [40, 30, 10, 20];
        let d = face_direction(404, 1, &g).unwrap();
        // 3 first, then the smaller neighbour (4: 20 < 2: 30), then 1, 2.
        assert_eq!(d[0], [3, 4, 1, 2]);
    }

    #[test]
    fn triangle_face_direction_sorts() {
        let g = [30, 10, 20];
        let d = face_direction(303, 1, &g).unwrap();
        assert_eq!(d[0], [2, 3, 1, 0]);
    }

    #[test]
    fn dof_counts_match_closed_forms() {
        // brick p=4: 8 + 12*3 + 6*9 + 27
        assert_eq!(p_total_dofs(808, 4, false), Some(8 + 36 + 54 + 27));
        // tetra p=4: 4 + 6*3 + 4*3 + 1
        assert_eq!(p_total_dofs(504, 4, false), Some(4 + 18 + 12 + 1));
        // quad p=3 (interior counted as bubble): 4 + 4*2 + 4
        assert_eq!(p_total_dofs(404, 3, false), Some(4 + 8 + 4));
        assert_eq!(p_dof_counts(706, 3), Some((6, 2, 1, 2)));
        assert_eq!(face_dofs(706, 3, 3, false), 4);
        assert_eq!(bubble_max_degree_1d(808, 6, true), 2);
    }
}
