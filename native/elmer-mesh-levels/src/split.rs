//! `SplitMeshEqual` (MeshSplit.F90:52) and the `Mesh Levels` loop of
//! `LoadModel` (ModelDescription.F90:3049).
//!
//! New nodes are numbered: old nodes, then one per old edge (edge order),
//! then one per quadrilateral face (face order), then one per 404/808 bulk
//! element (element order). Each old bulk element becomes 2, 4 or 8 children
//! in the tables below; each old boundary element with a parent becomes
//! 1, 2 or 4 children whose only parent (`Left`) is the child of the old
//! parent that contains it. Boundary elements without a parent are dropped,
//! and the old second parent (`Right`) is not carried over.

use elmer_core::{ElemId, MeshView, ParamMap, Params};

use crate::topology::find_mesh_edges;
use crate::{Element, Mesh, MeshLevelsError};

/// A child node in Fortran terms (1-based): old node `NodeIndexes(k)`,
/// midpoint of `EdgeIndexes(k)`, centre of `FaceIndexes(k)`, element centre.
#[derive(Copy, Clone)]
enum P {
    V(usize),
    E(usize),
    F(usize),
    C,
}
use P::{C, E, F, V};

const LINE: [[P; 2]; 2] = [[V(1), E(1)], [E(1), V(2)]];
const TRI: [[P; 3]; 4] = [
    [V(1), E(1), E(3)],
    [E(1), V(2), E(2)],
    [E(1), E(2), E(3)],
    [E(2), V(3), E(3)],
];
const QUAD: [[P; 4]; 4] = [
    [V(1), E(1), C, E(4)],
    [E(1), V(2), E(2), C],
    [C, E(2), V(3), E(3)],
    [E(4), C, E(3), V(4)],
];
const TET_CORNERS: [[P; 4]; 4] = [
    [V(1), E(1), E(3), E(4)],
    [V(2), E(1), E(5), E(2)],
    [V(3), E(2), E(6), E(3)],
    [V(4), E(4), E(6), E(5)],
];
/// The inner octahedron split along the diagonal between edges 2-4, 3-5 or
/// 1-6 (Delaunay choice: the shortest).
const TET_INNER: [[[P; 4]; 4]; 3] = [
    [
        [E(6), E(4), E(5), E(2)],
        [E(6), E(4), E(2), E(3)],
        [E(4), E(5), E(2), E(1)],
        [E(4), E(3), E(1), E(2)],
    ],
    [
        [E(5), E(4), E(6), E(3)],
        [E(5), E(1), E(2), E(3)],
        [E(3), E(2), E(6), E(5)],
        [E(3), E(1), E(4), E(5)],
    ],
    [
        [E(6), E(4), E(5), E(1)],
        [E(6), E(2), E(3), E(1)],
        [E(1), E(2), E(5), E(6)],
        [E(1), E(3), E(4), E(6)],
    ],
];
const WEDGE: [[P; 6]; 8] = [
    [V(1), E(1), E(3), E(7), F(3), F(5)],
    [E(1), V(2), E(2), F(3), E(8), F(4)],
    [E(3), E(2), V(3), F(5), F(4), E(9)],
    [E(1), E(2), E(3), F(3), F(4), F(5)],
    [E(7), F(3), F(5), V(4), E(4), E(6)],
    [F(3), E(8), F(4), E(4), V(5), E(5)],
    [F(5), F(4), E(9), E(6), E(5), V(6)],
    [F(3), F(4), F(5), E(4), E(5), E(6)],
];
const BRICK: [[P; 8]; 8] = [
    [V(1), E(1), F(1), E(4), E(9), F(3), C, F(6)],
    [E(1), V(2), E(2), F(1), F(3), E(10), F(4), C],
    [E(4), F(1), E(3), V(4), F(6), C, F(5), E(12)],
    [F(1), E(2), V(3), E(3), C, F(4), E(11), F(5)],
    [E(9), F(3), C, F(6), V(5), E(5), F(2), E(8)],
    [F(3), E(10), F(4), C, E(5), V(6), E(6), F(2)],
    [F(6), C, F(5), E(12), E(8), F(2), E(7), V(8)],
    [C, F(4), E(11), F(5), F(2), E(6), V(7), E(7)],
];

/// One level of uniform refinement, `SplitMeshEqual(Mesh)` without the
/// nodal element size `h` (`Mesh Keep Grading`), and without
/// `FollowCurvedBoundary` (callers refuse `Follow * Boundary` BCs; see
/// [`mesh_levels`]).
pub fn split_mesh_equal(mesh: &dyn MeshView) -> Result<Mesh, MeshLevelsError> {
    let old = Mesh::from_view(mesh);
    if old.bulk().iter().any(|e| e.etype.family() == 6) {
        return Err(MeshLevelsError::Unsupported(
            "pyramids in SplitMeshEqual".into(),
        ));
    }
    let topo = find_mesh_edges(&old)?;
    let n_old = old.coords.len();
    let n_edges = topo.edges.len();

    // Quad face centres get node numbers after the edge midpoints.
    let mut face_node = vec![None; topo.faces.len()];
    let mut node_cnt = n_old + n_edges;
    for (f, nodes) in topo.faces.iter().enumerate() {
        if nodes.len() == 4 {
            face_node[f] = Some(node_cnt as u32);
            node_cnt += 1;
        }
    }
    let first_centre = node_cnt;
    node_cnt += old
        .bulk()
        .iter()
        .filter(|e| matches!(e.etype.family(), 4 | 8))
        .count();

    // Fortran SUM starts from zero, so a pair of -0.0 averages to +0.0.
    let mean = |nodes: &[u32], c: usize| {
        nodes
            .iter()
            .fold(0.0, |s, &n| s + old.coords[n as usize][c])
            / nodes.len() as f64
    };
    let mean3 = |nodes: &[u32]| [mean(nodes, 0), mean(nodes, 1), mean(nodes, 2)];
    let mut coords = Vec::with_capacity(node_cnt);
    coords.extend_from_slice(&old.coords);
    coords.extend(topo.edges.iter().map(|e| mean3(e)));
    coords.extend(topo.faces.iter().filter(|f| f.len() == 4).map(|f| mean3(f)));
    coords.extend(
        old.bulk()
            .iter()
            .filter(|e| matches!(e.etype.family(), 4 | 8))
            .map(|e| mean3(&e.nodes)),
    );

    let mut elements = Vec::new();
    let mut children = vec![Vec::new(); old.n_bulk];
    let mut centre = first_centre;
    for (i, el) in old.bulk().iter().enumerate() {
        let edges = &topo.element_edges[i];
        let faces = &topo.element_faces[i];
        let mut c = None;
        if matches!(el.etype.0, 404 | 808) {
            c = Some(centre as u32);
            centre += 1;
        }
        let node = |p: P| -> Result<u32, MeshLevelsError> {
            Ok(match p {
                V(k) => el.nodes[k - 1],
                E(k) => (n_old + edges[k - 1] as usize) as u32,
                F(k) => face_node[faces[k - 1] as usize].ok_or_else(|| {
                    MeshLevelsError::Invalid("face centre of a triangular face".into())
                })?,
                C => c.expect("centre node of a 404/808"),
            })
        };
        let mut push = |rows: &[&[P]]| -> Result<(), MeshLevelsError> {
            for row in rows {
                let nodes = row.iter().map(|&p| node(p)).collect::<Result<_, _>>()?;
                children[i].push(elements.len() as u32);
                elements.push(Element {
                    etype: el.etype,
                    nodes,
                    tag: el.tag,
                    parents: [None, None],
                });
            }
            Ok(())
        };
        match el.etype.0 {
            101 => push(&[&[V(1)]])?,
            202 => push(&LINE.each_ref().map(|r| &r[..]))?,
            303 => push(&TRI.each_ref().map(|r| &r[..]))?,
            404 => push(&QUAD.each_ref().map(|r| &r[..]))?,
            504 => {
                push(&TET_CORNERS.each_ref().map(|r| &r[..]))?;
                let mid = |k: usize| coords[n_old + edges[k - 1] as usize];
                let dist = |a: usize, b: usize| {
                    let (p, q) = (mid(a), mid(b));
                    let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
                };
                let d = [dist(4, 2), dist(5, 3), dist(6, 1)];
                let diag = if d[1] < d[0] && d[1] < d[2] {
                    1
                } else if d[2] < d[0] && d[2] < d[1] {
                    2
                } else {
                    0
                };
                push(&TET_INNER[diag].each_ref().map(|r| &r[..]))?;
            }
            706 => push(&WEDGE.each_ref().map(|r| &r[..]))?,
            808 => push(&BRICK.each_ref().map(|r| &r[..]))?,
            c => {
                return Err(MeshLevelsError::Unsupported(format!(
                    "element type {c} in SplitMeshEqual"
                )))
            }
        }
    }
    let n_bulk = elements.len();

    for (j, el) in old.boundary().iter().enumerate() {
        let Some(parent) = el.parents[0].or(el.parents[1]) else {
            continue;
        };
        let v = &el.nodes;
        if el.etype.0 == 101 {
            // The parents are the latest two lines through the point.
            let mut parents = [None; 2];
            let mut k = 0;
            for (e, cand) in elements.iter().enumerate().rev() {
                if cand.etype.family() == 2 && cand.nodes.contains(&v[0]) {
                    parents[k] = Some(ElemId(e as u32));
                    k += 1;
                    if k == 2 {
                        break;
                    }
                }
            }
            elements.push(Element {
                etype: el.etype,
                nodes: v.clone(),
                tag: el.tag,
                parents,
            });
            continue;
        }
        let pid = parent.0 as usize;
        if pid >= old.n_bulk {
            return Err(MeshLevelsError::Invalid(format!(
                "boundary element {} has a non-bulk parent",
                j + 1
            )));
        }
        let pe = &topo.element_edges[pid];
        let edge_node = |a: u32, b: u32| -> Result<u32, MeshLevelsError> {
            pe.iter()
                .find(|&&e| {
                    let [p, q] = topo.edges[e as usize];
                    p.min(q) == a.min(b) && p.max(q) == a.max(b)
                })
                .map(|&e| (n_old + e as usize) as u32)
                .ok_or_else(|| {
                    MeshLevelsError::Invalid(format!(
                        "no parent edge with nodes {} {}",
                        a + 1,
                        b + 1
                    ))
                })
        };
        let rows: Vec<Vec<u32>> = match el.etype.0 {
            202 => {
                let m = edge_node(v[0], v[1])?;
                vec![vec![v[0], m], vec![m, v[1]]]
            }
            303 => {
                let (m12, m23, m31) = (
                    edge_node(v[0], v[1])?,
                    edge_node(v[1], v[2])?,
                    edge_node(v[2], v[0])?,
                );
                vec![
                    vec![v[0], m12, m31],
                    vec![m12, v[1], m23],
                    vec![m12, m23, m31],
                    vec![m31, m23, v[2]],
                ]
            }
            404 => {
                let mut key = v.clone();
                key.sort_unstable();
                let face = topo.element_faces[pid]
                    .iter()
                    .find(|&&f| {
                        let mut nodes = topo.faces[f as usize].clone();
                        nodes.sort_unstable();
                        nodes == key
                    })
                    .and_then(|&f| face_node[f as usize])
                    .ok_or_else(|| {
                        MeshLevelsError::Invalid(format!(
                            "no parent face for boundary element {}",
                            j + 1
                        ))
                    })?;
                let (m12, m23, m34, m41) = (
                    edge_node(v[0], v[1])?,
                    edge_node(v[1], v[2])?,
                    edge_node(v[2], v[3])?,
                    edge_node(v[3], v[0])?,
                );
                vec![
                    vec![v[0], m12, face, m41],
                    vec![m12, v[1], m23, face],
                    vec![m41, face, m34, v[3]],
                    vec![face, m23, v[2], m34],
                ]
            }
            c => {
                return Err(MeshLevelsError::Unsupported(format!(
                    "boundary element type {c} in SplitMeshEqual"
                )))
            }
        };
        for nodes in rows {
            // Lines need both nodes in the child, faces three of theirs.
            let need = if nodes.len() == 2 { 2 } else { 3 };
            let child = children[pid]
                .iter()
                .copied()
                .find(|&c| {
                    nodes
                        .iter()
                        .filter(|n| elements[c as usize].nodes.contains(n))
                        .count()
                        >= need
                })
                .ok_or_else(|| {
                    MeshLevelsError::Invalid(format!(
                        "parent element not found for boundary {}",
                        j + 1
                    ))
                })?;
            elements.push(Element {
                etype: el.etype,
                nodes,
                tag: el.tag,
                parents: [Some(ElemId(child)), None],
            });
        }
    }
    Ok(Mesh {
        dim: old.dim,
        coords,
        elements,
        n_bulk,
    })
}

/// The Simulation `Mesh Levels` loop: `Mesh Levels - 1` splits of the loaded
/// mesh, returning the finest (the active mesh), or `None` when there is
/// nothing to split. `bcs` are the Boundary Condition sections, checked for
/// the `Follow` keywords that `FollowCurvedBoundary` would apply.
pub fn mesh_levels(
    mesh: &dyn MeshView,
    simulation: &ParamMap,
    bcs: &[ParamMap],
) -> Result<Option<Mesh>, MeshLevelsError> {
    let levels = simulation.integer("mesh levels").unwrap_or(1);
    if levels <= 1 {
        return Ok(None);
    }
    check_prepare_mesh(simulation)?;
    if simulation.logical("mesh keep grading") == Some(true) {
        return Err(MeshLevelsError::Unsupported("Mesh Keep Grading".into()));
    }
    if bcs
        .iter()
        .any(|bc| bc.0.keys().any(|k| k.starts_with("follow")))
    {
        return Err(MeshLevelsError::Unsupported(
            "Follow * Boundary with Mesh Levels".into(),
        ));
    }
    let mut cur = split_mesh_equal(mesh)?;
    for _ in 2..levels {
        cur = split_mesh_equal(&cur)?;
    }
    Ok(Some(cur))
}

/// `PrepareMesh` options that rework a freshly built mesh again.
pub(crate) fn check_prepare_mesh(simulation: &ParamMap) -> Result<(), MeshLevelsError> {
    for key in [
        "split mesh quads",
        "increase element order",
        "constant stencil",
    ] {
        if simulation.logical(key) == Some(true) {
            return Err(MeshLevelsError::Unsupported(format!(
                "{key} on a refined or extruded mesh"
            )));
        }
    }
    Ok(())
}
