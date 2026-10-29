//! Element-type table generated from the `elements.def` snapshot.
//!
//! Lookup mirrors `GetElementType`: unknown codes are an error at load
//! (Elmer warns in `GetElementType` and its callers fatal).

include!(concat!(env!("OUT_DIR"), "/element_table.rs"));

/// Type info for an Elmer element code, or `None` when unknown.
pub fn element_info(code: u16) -> Option<ElementTypeInfo> {
    ELEMENT_TYPES.iter().find(|e| e.code == code).copied()
}

/// Reference node coordinates (U, V, W) for a type, or `None`.
pub fn reference_nodes(code: u16) -> Option<&'static [[f64; 3]]> {
    ELEMENT_NODES.iter().find(|(c, _)| *c == code).map(|(_, n)| *n)
}

/// Canonical linear edge endpoints per element family, 0-based node
/// positions. Matches `FindMeshEdges2D` (consecutive pairs) and the
/// `TetraEdgeMap`/`BrickEdgeMap` tables in `FindMeshEdges3D`.
pub fn linear_edges(code: u16) -> Option<&'static [(usize, usize)]> {
    Some(match code / 100 {
        2 => &[(0, 1)],
        3 => &[(0, 1), (1, 2), (2, 0)],
        4 => &[(0, 1), (1, 2), (2, 3), (3, 0)],
        5 => &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
        // NOTE: Elmer numbers 6xx = pyramids (8 edges) and 7xx = wedges
        // (9 edges); the elements.def names are the other way round.
        6 => &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (0, 4),
            (1, 4),
            (2, 4),
            (3, 4),
        ],
        7 => &[
            (0, 1),
            (1, 2),
            (0, 2),
            (3, 4),
            (4, 5),
            (5, 3),
            (0, 3),
            (1, 4),
            (2, 5),
        ],
        8 => &[
            (0, 1),
            (1, 2),
            (3, 2),
            (0, 3),
            (4, 5),
            (5, 6),
            (7, 6),
            (4, 7),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ],
        _ => return None,
    })
}
