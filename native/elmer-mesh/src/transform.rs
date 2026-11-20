//! Mesh-load transforms, replicating `MeshLoad.F90` / `MeshTagging.F90`.
//!
//! Implemented (S1 scope, report §5.3):
//! - `MapCoordinates` (coordinate mapping + scaling, `MeshLoad.F90:558`),
//! - `MapBodiesAndBCs` (Target Bodies/Boundaries, `MeshLoad.F90:230`),
//! - `MapInternalExternalBCs` (`MeshLoad.F90:466`),
//! - `CreateIntersectionBCs` in 2D, both modes (`MeshTagging.F90:52`),
//! - `IncreaseElementOrder` (`MeshBasics.F90:3564`),
//! - `SplitMeshQuads` (`MeshSplit.F90:2058`): `Split Mesh Quads` (2D) and
//!   `Split Mesh Prisms` (3D).
//!
//! Explicitly out of scope for this PR: 3D intersection BCs (need face
//! tables), `TagBCsUsingRule`, curved-boundary
//! following, and all parallel handling. Those are capability refusals
//! in a later lane, not silent skips.

use crate::mesh::{BoundElem, Elem, Mesh, MeshError};
use elmer_sif::model::{Entry, Model, Section, Value};

fn fatal(msg: String) -> MeshError {
    MeshError { msg }
}

fn int_entry(s: &Section, name: &str) -> Option<i32> {
    match s.get(name)?.value {
        Value::IntArray { ref values, .. } => values.first().copied(),
        _ => None,
    }
}

fn int_array(s: &Section, name: &str) -> Option<Vec<i32>> {
    match s.get(name)?.value {
        Value::IntArray { ref values, .. } => Some(values.clone()),
        _ => None,
    }
}

fn logical(s: &Section, name: &str) -> bool {
    matches!(s.get(name), Some(Entry { value: Value::Logical(true), .. }))
}

fn present(s: &Section, name: &str) -> bool {
    s.get(name).is_some()
}

/// Full S1 load pipeline: read is done by [`crate::mesh::load`]; this
/// applies the model-dependent transforms in Elmer's order.
pub fn transform(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    map_coordinates(mesh, model)?;
    map_bodies_and_bcs(mesh, model)?;
    map_internal_external_bcs(mesh, model)?;
    split_mesh_quads(mesh, model)?;
    create_intersection_bcs(mesh, model)?;
    if sim_logical(model, "increase element order") {
        increase_element_order(mesh)?;
    }
    Ok(())
}

fn sim_logical(model: &Model, name: &str) -> bool {
    model.simulation.as_ref().map(|s| logical(s, name)).unwrap_or(false)
}

fn sim_section(model: &Model) -> Option<&Section> {
    model.simulation.as_ref()
}

/// `WedgeFaceMap(3:5, 1:4)` (MeshBasics.F90:4182-4184): the quadrilateral
/// faces of a prism, 0-based.
const WEDGE_QUAD_FACES: [[usize; 4]; 3] = [[0, 1, 4, 3], [2, 1, 4, 5], [2, 0, 3, 5]];

/// The split of a prism by the cut corners of its three quadrilateral faces
/// (MeshSplit.F90:2329-2405): `CutComb` and the three tetrahedra.
fn cut_comb(nodes: &[u32], cut: [u32; 3]) -> u32 {
    let is24 = cut[0] == nodes[1] || cut[0] == nodes[3];
    let is26 = cut[1] == nodes[1] || cut[1] == nodes[5];
    let is16 = cut[2] == nodes[0] || cut[2] == nodes[3];
    111 + if is24 { 0 } else { 100 } + if is26 { 0 } else { 10 } + if is16 { 0 } else { 1 }
}

fn prism_cut(nodes: &[u32], cut: [u32; 3]) -> Result<(u32, [[usize; 4]; 3]), MeshError> {
    let comb = cut_comb(nodes, cut);
    let maps: [[usize; 4]; 3] = match comb {
        111 => [[1, 2, 4, 6], [2, 5, 4, 6], [1, 2, 3, 6]],
        112 => [[1, 2, 4, 3], [2, 5, 4, 6], [2, 3, 6, 4]],
        122 => [[1, 2, 4, 3], [2, 5, 4, 3], [4, 5, 3, 6]],
        221 => [[1, 2, 5, 3], [1, 5, 4, 6], [1, 3, 5, 6]],
        222 => [[1, 2, 5, 3], [1, 5, 4, 3], [4, 5, 3, 6]],
        211 => [[1, 2, 5, 6], [1, 5, 4, 6], [1, 2, 3, 6]],
        c => return Err(fatal(format!("SplitMeshQuads: Unknown case for split: {c}"))),
    };
    Ok((comb, maps))
}

/// `SplitMeshQuads(Mesh, Model % Simulation)` (MeshSplit.F90:2058): with
/// `Split Mesh Quads` (2D) the quadrilaterals become two triangles each;
/// with `Split Mesh Prisms` (3D) each prism becomes three tetrahedra and
/// each boundary quadrilateral two triangles. The diagonal of a
/// quadrilateral face goes through its cut corner: the smallest node
/// (`Split Mesh Prisms Min/Max`: the shorter or longer diagonal), changed
/// prism by prism until every prism can be cut. A boundary quadrilateral is
/// split in the node order of its face as `FindMeshFaces3D` stores it (that
/// of the first element with the face), and its parents become the new
/// elements that contain it. Pyramids and hexahedra are refused, as in
/// Elmer.
fn split_mesh_quads(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    let dim = mesh.dim;
    if dim < 2 {
        return Ok(());
    }
    if !sim_logical(model, if dim == 2 { "split mesh quads" } else { "split mesh prisms" }) {
        return Ok(());
    }
    let mut cnt = [0usize; 9];
    for code in mesh.bulk.iter().map(|e| e.code).chain(mesh.boundary.iter().map(|e| e.code)) {
        cnt[(code / 100) as usize] += 1;
    }
    if cnt[6] + cnt[8] > 0 {
        return Err(fatal("SplitMeshQuads: Not implemented yet for pyramids and hexahedrons!".into()));
    }
    if dim == 3 && cnt[7] == 0 {
        return Ok(());
    }
    if mesh.bulk.iter().chain(std::iter::empty()).any(|e| ![101, 202, 303, 404, 504, 706].contains(&e.code))
        || mesh.boundary.iter().any(|e| ![101, 202, 303, 404].contains(&e.code))
    {
        return Err(fatal("SplitMeshQuads: higher-order elements are not split".into()));
    }
    // FindMeshFaces3D: the quadrilateral faces, in the node order of the
    // first element that has them; MinCorner and CutCorner.
    let mut face_of: std::collections::HashMap<[u32; 4], usize> = std::collections::HashMap::new();
    let mut face_nodes: Vec<[u32; 4]> = Vec::new();
    let mut elem_faces: Vec<[usize; 3]> = vec![[0; 3]; mesh.bulk.len()];
    for (i, e) in mesh.bulk.iter().enumerate() {
        if e.code != 706 {
            continue;
        }
        for (f, map) in WEDGE_QUAD_FACES.iter().enumerate() {
            let fnodes = map.map(|k| e.nodes[k]);
            let mut key = fnodes;
            key.sort_unstable();
            let id = *face_of.entry(key).or_insert_with(|| {
                face_nodes.push(fnodes);
                face_nodes.len() - 1
            });
            elem_faces[i][f] = id;
        }
    }
    let min_corner: Vec<u32> = face_nodes.iter().map(|f| *f.iter().min().unwrap()).collect();
    let mut cut = min_corner.clone();
    let do_max = if sim_logical(model, "split mesh prisms max") {
        -1.0
    } else if sim_logical(model, "split mesh prisms min") {
        1.0
    } else {
        0.0
    };
    if do_max != 0.0 {
        let d2 = |a: u32, b: u32| {
            let (p, q) = (mesh.nodes[a as usize], mesh.nodes[b as usize]);
            (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)
        };
        for (f, n) in face_nodes.iter().enumerate() {
            let (d13, d24) = (d2(n[0], n[2]), d2(n[1], n[3]));
            cut[f] = if do_max * d13 < do_max * d24 { n[0].min(n[2]) } else { n[1].min(n[3]) };
        }
    }
    if dim == 3 {
        loop {
            let mut changes = 0;
            for (i, e) in mesh.bulk.iter().enumerate() {
                if e.code != 706 {
                    continue;
                }
                let lc = elem_faces[i].map(|f| cut[f]);
                let comb = cut_comb(&e.nodes, lc);
                if comb == 121 || comb == 212 {
                    let mx = *lc.iter().max().unwrap();
                    let k = lc.iter().position(|&c| c == mx).unwrap() + 1;
                    let icut = match (k, comb) {
                        (1, 121) => 1,
                        (1, _) => 2,
                        (2, 121) => 2,
                        (2, _) => 3,
                        (3, 121) => 3,
                        (_, _) => 1,
                    };
                    cut[elem_faces[i][k - 1]] = e.nodes[icut - 1];
                    changes += 1;
                }
            }
            if changes == 0 {
                break;
            }
        }
    }
    // The new elements: bulk first, then boundary, in the old order.
    let mut bulk: Vec<Elem> = Vec::with_capacity(mesh.bulk.len() + 2 * cnt[7] + cnt[4]);
    let mut offset = vec![0usize; mesh.bulk.len()];
    for (i, e) in mesh.bulk.iter().enumerate() {
        offset[i] = bulk.len();
        let mut child = |code: u16, nodes: Vec<u32>| {
            bulk.push(Elem { code, body: e.body, nodes, tag: e.tag, edges: Vec::new() });
        };
        match e.code {
            404 => {
                for map in [[0, 1, 2], [2, 3, 0]] {
                    child(303, map.iter().map(|&k| e.nodes[k]).collect());
                }
            }
            706 => {
                let (_, maps) = prism_cut(&e.nodes, elem_faces[i].map(|f| cut[f]))?;
                for map in maps {
                    child(504, map.iter().map(|&k| e.nodes[k - 1]).collect());
                }
            }
            c => child(c, e.nodes.clone()),
        }
    }
    let n_cands = |code: u16| match code {
        404 => 2,
        706 => 3,
        _ => 1,
    };
    let mut boundary: Vec<BoundElem> = Vec::with_capacity(mesh.boundary.len() + cnt[4]);
    for b in &mesh.boundary {
        let pieces: Vec<(u16, Vec<u32>)> = match b.code {
            404 => {
                // The face of the parent with these nodes, and its cut.
                let mut key = [b.nodes[0], b.nodes[1], b.nodes[2], b.nodes[3]];
                key.sort_unstable();
                let f = *face_of
                    .get(&key)
                    .ok_or_else(|| fatal("SplitMeshQuads: boundary quadrilateral without a prism face".into()))?;
                let fnodes = face_nodes[f];
                let odd = cut[f] == fnodes[0] || cut[f] == fnodes[2];
                let maps = if odd { [[0, 1, 2], [2, 3, 0]] } else { [[3, 0, 1], [1, 2, 3]] };
                maps.iter().map(|m| (303, m.iter().map(|&k| fnodes[k]).collect())).collect()
            }
            c => vec![(c, b.nodes.clone())],
        };
        for (code, nodes) in pieces {
            let mut parents = [None, None];
            for (j, p) in b.parents.iter().enumerate() {
                let Some(p) = *p else { continue };
                let off = offset[p];
                let hit = (off..off + n_cands(mesh.bulk[p].code))
                    .find(|&c| nodes.iter().all(|n| bulk[c].nodes.contains(n)))
                    .ok_or_else(|| {
                        fatal(format!("UpdateParentElements: Could not find parent for type {}", mesh.bulk[p].code))
                    })?;
                parents[j] = Some(hit);
            }
            boundary.push(BoundElem { code, constraint: b.constraint, parents, nodes, tag: b.tag, edges: Vec::new() });
        }
    }
    mesh.bulk = bulk;
    mesh.boundary = boundary;
    if mesh.edges.is_some() {
        mesh.edges = None;
        mesh.find_edges();
    }
    Ok(())
}

fn map_coordinates(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    let sim = match sim_section(model) {
        Some(s) => s,
        None => return Ok(()),
    };
    if let Some(map) = int_array(sim, "coordinate mapping") {
        if map.len() != 3 || ![1, 2, 3].iter().all(|c| map.contains(c)) {
            return Err(fatal(format!("inconsistent coordinate mapping: {map:?}")));
        }
        // Elmer swaps whole coordinate vectors through pointers.
        let xs: Vec<f64> = mesh.nodes.iter().map(|n| n[0]).collect();
        let ys: Vec<f64> = mesh.nodes.iter().map(|n| n[1]).collect();
        let zs: Vec<f64> = mesh.nodes.iter().map(|n| n[2]).collect();
        let cols = [&xs, &ys, &zs];
        for (i, n) in mesh.nodes.iter_mut().enumerate() {
            for (k, &c) in map.iter().enumerate() {
                n[k] = cols[(c - 1) as usize][i];
            }
        }
    }
    if let Some(scale) = sim_real_array(sim, "coordinate scaling") {
        // Wrk(j,1) with j = MIN(i, SIZE(Wrk,1)): a short vector repeats
        // its last entry for the remaining dimensions.
        let dim = mesh.dim.max(1);
        let mut factors = [1.0; 3];
        if !scale.is_empty() {
            for (f, j) in factors.iter_mut().zip(0..dim) {
                let k = j.min(scale.len() - 1);
                *f = scale[k];
            }
        }
        for n in &mut mesh.nodes {
            for (k, f) in factors.iter().enumerate() {
                // Elmer scales x always, y when dim>1, z when dim>2.
                if k == 0 || (k == 1 && dim > 1) || (k == 2 && dim > 2) {
                    n[k] *= *f;
                }
            }
        }
    }
    // Mesh dimension is set at read; remap here as Elmer re-derives it
    // after mapping (SetMeshDimension is idempotent on the values).
    crate::mesh::set_mesh_dimension(mesh);
    Ok(())
}

fn sim_real_array(sim: &Section, name: &str) -> Option<Vec<f64>> {
    match sim.get(name)?.value {
        Value::RealArray { ref values, .. } => Some(values.clone()),
        _ => None,
    }
}

// ---- MapBodiesAndBCs ----

fn map_bodies_and_bcs(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    // Bodies: Target Bodies remap + Default Target body.
    let found = model.bodies.iter().any(|b| {
        present(b, "target bodies")
            || logical(b, "default target")
    });
    if found {
        if mesh.bulk.is_empty() {
            return Err(fatal("body indexes are screwed!".to_string()));
        }
        let minid = mesh.bulk.iter().map(|e| e.body).min().unwrap_or(1).min(1);
        let maxid = mesh
            .bulk
            .iter()
            .map(|e| e.body)
            .max()
            .unwrap_or(0)
            .max(model.bodies.len() as i32);
        let off = minid;
        let mut index_map = vec![0i32; (maxid - minid + 1) as usize];
        let at = |id: i32| (id - off) as usize;
        let mut default_target = 0;
        for (i, b) in model.bodies.iter().enumerate() {
            let id = i as i32 + 1;
            if logical(b, "default target") {
                default_target = id;
            }
        }
        for (i, b) in model.bodies.iter().enumerate() {
            let id = i as i32 + 1;
            if let Some(blist) = int_array(b, "target bodies") {
                for body in blist {
                    if body > maxid || body < minid {
                        continue;
                    } else if index_map[at(body)] != 0 {
                        // Multiple bodies share a target: warn, keep first.
                    } else {
                        index_map[at(body)] = id;
                    }
                }
            } else if default_target == 0 {
                if index_map[at(id)] != 0 {
                    // already set by a target list; keep
                } else {
                    index_map[at(id)] = id;
                }
            }
        }
        for e in &mut mesh.bulk {
            let id = e.body;
            if id >= minid && id <= maxid {
                if index_map[at(id)] == 0 && default_target != 0 {
                    index_map[at(id)] = default_target;
                }
                e.body = index_map[at(id)];
            }
        }
    }
    let _ = found;

    if mesh.boundary.is_empty() {
        return Ok(());
    }
    // Boundaries: Boundary N -> Boundary Condition id, then Target
    // Boundaries / Default Target BC remap.
    let minid = mesh
        .boundary
        .iter()
        .map(|e| e.constraint)
        .min()
        .unwrap_or(1)
        .min(1);
    let maxid = mesh
        .boundary
        .iter()
        .map(|e| e.constraint)
        .max()
        .unwrap_or(0)
        .max(model.bcs.len() as i32);
    if minid > maxid {
        return Err(fatal("boundary indexes are screwed".to_string()));
    }
    let off = minid;
    let mut index_map = vec![0i32; (maxid - minid + 1) as usize];
    let at = |id: i32| (id - off) as usize;
    for nb in &model.boundaries {
        let id = int_entry(&nb.section, "boundary condition").unwrap_or(0);
        if id == 0 {
            continue;
        }
        let bndry = nb.id as i32;
        if bndry > maxid {
            // Warn: BoundaryId exceeds range (kept as warning).
        } else if bndry == 0 {
            // Warn: BoundaryId is zero.
        } else {
            index_map[at(bndry)] = id;
        }
    }
    let mut default_target_bc = 0;
    for (i, bc) in model.bcs.iter().enumerate() {
        let id = i as i32 + 1;
        if logical(bc, "default target") || logical(bc, "default bc") {
            default_target_bc = id;
        }
        if let Some(blist) = int_array(bc, "target boundaries") {
            for bndry in blist {
                if bndry > maxid {
                    continue;
                } else if index_map[at(bndry)] != 0 {
                    // Multiple BCs share a target: warn, keep first.
                } else {
                    index_map[at(bndry)] = id;
                }
            }
        } else if !(present(bc, "target nodes") || present(bc, "target coordinates"))
            && index_map[at(id)] != 0
            && id == default_target_bc
        {
            // Default target already listed: warn, keep.
        }
    }
    for e in &mut mesh.boundary {
        let bndry = e.constraint;
        if bndry < minid || bndry > maxid {
            // Warn: not in range; Elmer then treats via IndexMap default.
        }
        let mapped = if bndry >= minid && bndry <= maxid {
            index_map[at(bndry)]
        } else {
            0
        };
        if mapped < 0 {
            e.constraint = 0;
            continue;
        } else if mapped == 0 {
            if default_target_bc != 0 {
                if bndry >= minid && bndry <= maxid {
                    index_map[at(bndry)] = default_target_bc;
                }
                e.constraint = default_target_bc;
            } else {
                if bndry >= minid && bndry <= maxid {
                    index_map[at(bndry)] = -1;
                }
                e.constraint = 0;
                continue;
            }
        } else {
            e.constraint = mapped;
        }
        // Body Id / Normal Target Body of the owning BC.
        if e.constraint >= 1 && (e.constraint as usize) <= model.bcs.len() {
            let bc = &model.bcs[e.constraint as usize - 1];
            // Body Id (default 0 when absent).
            // NOTE: Elmer also sets OutBody (Normal Target Body); the S1
            // mesh keeps parents only — recorded for the solver lane.
            let _ = bc;
        }
    }
    Ok(())
}

// ---- MapInternalExternalBCs ----

fn map_internal_external_bcs(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    if mesh.boundary.is_empty() {
        return Ok(());
    }
    let mut default_int = 0;
    let mut default_ext = 0;
    for (i, bc) in model.bcs.iter().enumerate() {
        let id = i as i32 + 1;
        if logical(bc, "default internal bc") {
            default_int = id;
        }
        if logical(bc, "default external bc") {
            default_ext = id;
        }
    }
    if default_int == 0 && default_ext == 0 {
        return Ok(());
    }
    for e in &mut mesh.boundary {
        if e.constraint != 0 {
            continue;
        }
        // Only dim-1 elements: 3D needs code>=300, 2D needs code>=200.
        if mesh.dim == 3 && e.code < 300 {
            continue;
        }
        if mesh.dim == 2 && e.code < 200 {
            continue;
        }
        let mut cnt = 0;
        if e.parents[0].is_some() {
            cnt += 1;
        }
        if e.parents[1].is_some() {
            cnt += 1;
        }
        // Serial: no parallel shared-node correction.
        if cnt == 2 && default_int > 0 {
            e.constraint = default_int;
        } else if cnt == 1 && default_ext > 0 {
            e.constraint = default_ext;
        }
    }
    Ok(())
}

// ---- CreateIntersectionBCs (2D) ----

fn create_intersection_bcs(mesh: &mut Mesh, model: &Model) -> Result<(), MeshError> {
    // Collect (new_tag, parent_a, parent_b, is_body_mode) rows.
    struct Row {
        new_tag: i32,
        pa: i32,
        pb: i32,
        body_mode: bool,
    }
    let mut rows: Vec<Row> = Vec::new();
    for (i, bc) in model.bcs.iter().enumerate() {
        let id = i as i32 + 1;
        let tag = bc.id as i32;
        if let Some(v) = int_array(bc, "intersection bc") {
            if v.len() != 2 {
                return Err(fatal("intersection BCs only available for two parents!".to_string()));
            }
            rows.push(Row { new_tag: tag, pa: v[0], pb: v[1], body_mode: false });
        } else if let Some(v) = int_array(bc, "intersection body") {
            if v.len() != 2 {
                return Err(fatal("intersection BCs only available for two parents!".to_string()));
            }
            rows.push(Row { new_tag: tag, pa: v[0], pb: v[1], body_mode: true });
        }
        let _ = id;
    }
    if rows.is_empty() {
        return Ok(());
    }
    if mesh.dim == 3 {
        if rows.iter().any(|r| r.body_mode) {
            return Err(fatal("3D Intersection Body needs face tables (later lane)".to_string()));
        }
        // `CreateIntersection3D`: an edge shared by two boundary faces of
        // the two parent BCs becomes a line element (once per edge over all
        // rows), with the edge's nodes as `Mesh % Edges` has them and the
        // two faces as its parents.
        if mesh.edges.is_none() {
            mesh.find_edges();
        }
        let n_edges = mesh.n_edges();
        let mut first: Vec<Option<(u32, u32)>> = vec![None; n_edges];
        for b in &mesh.bulk {
            let def = crate::elements::linear_edges(b.code).unwrap_or(&[]);
            for (k, &eid) in b.edges.iter().enumerate() {
                let e = eid as usize;
                if first[e].is_none() {
                    if let Some(&(a, c)) = def.get(k) {
                        first[e] = Some((b.nodes[a], b.nodes[c]));
                    }
                }
            }
        }
        let nbc = mesh.boundary.len();
        let mut owner_of = vec![0i32; nbc];
        for (t, e) in mesh.boundary.iter().enumerate() {
            if e.code < 300 {
                continue;
            }
            let Some(i) = model.bcs.iter().position(|bc| e.constraint == bc.id as i32) else { continue };
            let owner = i as i32 + 1;
            if rows.iter().any(|r| r.pa == owner || r.pb == owner) {
                owner_of[t] = owner;
            }
        }
        let mut done = vec![false; n_edges];
        let mut new3d: Vec<BoundElem> = Vec::new();
        for row in &rows {
            for t in 0..nbc {
                let j = owner_of[t];
                if j == 0 || (j != row.pa && j != row.pb) {
                    continue;
                }
                for t2 in t + 1..nbc {
                    let j2 = owner_of[t2];
                    if j2 == 0 || j2 == j || (j2 != row.pa && j2 != row.pb) {
                        continue;
                    }
                    let (a, b) = (&mesh.boundary[t], &mesh.boundary[t2]);
                    if a.nodes.iter().filter(|n| b.nodes.contains(n)).count() < 2 {
                        continue;
                    }
                    let Some(&e) = a.edges.iter().find(|e| !done[**e as usize] && b.edges.contains(e)) else { continue };
                    done[e as usize] = true;
                    let (n1, n2) = first[e as usize].unwrap_or(mesh.edges.as_ref().expect("found").edges[e as usize]);
                    new3d.push(BoundElem {
                        code: 202,
                        constraint: row.new_tag,
                        parents: [Some(mesh.bulk.len() + t), Some(mesh.bulk.len() + t2)],
                        nodes: vec![n1, n2],
                        tag: 0,
                        edges: Vec::new(),
                    });
                }
            }
        }
        mesh.boundary.extend(new3d);
        return Ok(());
    }
    if mesh.dim != 2 {
        return Ok(());
    }
    // BulkMode (Intersection Body) in 2D needs the global edge table.
    let need_edges = rows.iter().any(|r| r.body_mode);
    if need_edges && mesh.edges.is_none() {
        mesh.find_edges();
    }
    // Candidate parent boundary elements (2D boundary elements only).
    let nbc = mesh.boundary.len();
    let mut parent_of = vec![0i32; nbc];
    for (t, e) in mesh.boundary.iter().enumerate() {
        if e.code < 200 {
            continue;
        }
        // Resolve the BC id owning this constraint tag.
        let mut owner = 0;
        for (i, bc) in model.bcs.iter().enumerate() {
            if e.constraint == bc.id as i32 {
                owner = i as i32 + 1;
                break;
            }
        }
        if owner == 0 {
            continue;
        }
        if rows.iter().any(|r| r.pa == owner || r.pb == owner) {
            parent_of[t] = owner;
        }
    }
    // Two passes: count, then append (Elmer allocates once).
    struct New {
        node: u32,
        tag: i32,
        /// A BulkMode line element (`Edge % TYPE`, its nodes): the edge
        /// between the two bodies, with its bulk parents.
        line: Option<(u16, Vec<u32>)>,
        /// Boundary-element positions (into `mesh.boundary`).
        left_bpos: usize,
        right_bpos: usize,
        /// Bulk-element indices (BulkMode).
        left_bulk: Option<usize>,
        right_bulk: Option<usize>,
    }
    let mut new_elems: Vec<New> = Vec::new();
    // NodeDone is per-row in Elmer (reset for each l).
    for row in &rows {
        if row.body_mode {
            // `Mesh % Edges` in FindMeshEdges2D's order: each edge with the
            // nodes of the first bulk element that has it, `Left` that
            // element and `Right` the second. An edge between the two
            // bodies (0: no element on that side) becomes a boundary
            // element of the edge's type. (2D BulkMode; 3D is refused above.)
            let n_edges = mesh.n_edges();
            let mut first: Vec<Option<(u32, u32)>> = vec![None; n_edges];
            let mut left: Vec<Option<usize>> = vec![None; n_edges];
            let mut right: Vec<Option<usize>> = vec![None; n_edges];
            for (bi, b) in mesh.bulk.iter().enumerate() {
                let def = crate::elements::linear_edges(b.code).unwrap_or(&[]);
                for (k, &eid) in b.edges.iter().enumerate() {
                    let e = eid as usize;
                    if first[e].is_none() {
                        if let Some(&(a, c)) = def.get(k) {
                            first[e] = Some((b.nodes[a], b.nodes[c]));
                        }
                    }
                    if left[e].is_some() {
                        right[e] = Some(bi);
                    } else {
                        left[e] = Some(bi);
                    }
                }
            }
            for e in 0..n_edges {
                let j = left[e].map_or(0, |b| mesh.bulk[b].body);
                let j2 = right[e].map_or(0, |b| mesh.bulk[b].body);
                if j == j2 || (j != row.pa && j != row.pb) || (j2 != row.pa && j2 != row.pb) {
                    continue;
                }
                let Some((a, c)) = first[e] else { continue };
                new_elems.push(New {
                    node: a,
                    tag: row.new_tag,
                    line: Some((202, vec![a, c])),
                    left_bpos: 0,
                    right_bpos: 0,
                    left_bulk: left[e],
                    right_bulk: right[e],
                });
            }
        } else {
            let mut node_done = vec![false; mesh.nodes.len()];
            for t in 0..nbc {
                let j = parent_of[t];
                if j == 0 || (j != row.pa && j != row.pb) {
                    continue;
                }
                for t2 in t + 1..nbc {
                    let j2 = parent_of[t2];
                    if j2 == 0 || j2 == j || (j2 != row.pa && j2 != row.pb) {
                        continue;
                    }
                    // Shared nodes; require exactly one.
                    let (a, b) = {
                        let ea = &mesh.boundary[t];
                        let eb = &mesh.boundary[t2];
                        (ea.nodes.clone(), eb.nodes.clone())
                    };
                    let mut e = 0u32;
                    let mut k = 0;
                    for &na in &a {
                        if b.contains(&na) {
                            e = na;
                            k += 1;
                        }
                    }
                    if k != 1 || node_done[e as usize] {
                        continue;
                    }
                    node_done[e as usize] = true;
                    new_elems.push(New {
                        node: e,
                        tag: row.new_tag,
                        line: None,
                        left_bpos: t,
                        right_bpos: t2,
                        left_bulk: None,
                        right_bulk: None,
                    });
                }
            }
        }
    }
    // Append 101 point elements. Boundary-mode parents are the two
    // parent boundary elements (as bulk positions); bulk-mode parents
    // are the two bulk elements.
    let n_bulk = mesh.bulk.len();
    for ne in new_elems {
        if let Some((code, nodes)) = ne.line {
            mesh.boundary.push(BoundElem {
                code,
                constraint: ne.tag,
                parents: [ne.left_bulk, ne.right_bulk],
                nodes,
                tag: 0,
                edges: Vec::new(),
            });
            continue;
        }
        let parents = [Some(n_bulk + ne.left_bpos), Some(n_bulk + ne.right_bpos)];
        mesh.boundary.push(BoundElem {
            code: 101,
            constraint: ne.tag,
            parents,
            nodes: vec![ne.node],
            tag: 0,
            edges: Vec::new(),
        });
    }
    Ok(())
}

// ---- IncreaseElementOrder ----

pub(crate) fn increase_element_order(mesh: &mut Mesh) -> Result<(), MeshError> {
    if mesh.edges.is_none() {
        mesh.find_edges();
    }
    let n0 = mesh.nodes.len();
    let _n1 = mesh.n_edges();
    // New nodes at edge midpoints.
    let edges = mesh.edges.as_ref().cloned().unwrap_or_default();
    for &(n1n, n2n) in &edges.edges {
        let (a, b) = (mesh.nodes[n1n as usize], mesh.nodes[n2n as usize]);
        mesh.nodes.push([
            0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2]),
        ]);
    }
    // Rewrite elements: linear nodes + mid-edge nodes in canonical order.
    let rewrite = |nodes: &[u32], edge_ids: &[u32], code: u16| -> Result<(u16, Vec<u32>), MeshError> {
        let m1 = nodes.len();
        let m2 = edge_ids.len();
        let new_code = match code {
            202 | 303 | 404 | 504 | 808 => code + m2 as u16,
            101 => return Ok((101, nodes.to_vec())),
            _ => {
                return Err(fatal(format!(
                    "cannot increase element order for: {code}"
                )));
            }
        };
        let expect_m1 = match code {
            202 => 2,
            303 => 3,
            404 => 4,
            504 => 4,
            808 => 8,
            _ => m1,
        };
        if m1 != expect_m1 {
            return Err(fatal(format!(
                "cannot increase element order for: {code} (got {m1} nodes)"
            )));
        }
        let mut new_nodes = nodes.to_vec();
        for &eid in edge_ids {
            new_nodes.push(n0 as u32 + eid);
        }
        if code == 808 {
            // Elmer's 820 orders the top/bottom face mids swapped vs the
            // edge order; replicate the swap verbatim.
            if new_nodes.len() != 20 {
                return Err(fatal("cannot increase element order for: 808".to_string()));
            }
            let t12 = new_nodes[12..16].to_vec();
            let t16 = new_nodes[16..20].to_vec();
            new_nodes[16..20].copy_from_slice(&t12);
            new_nodes[12..16].copy_from_slice(&t16);
        }
        Ok((new_code, new_nodes))
    };
    for e in mesh.bulk.iter_mut() {
        let (code, nodes) = rewrite(&e.nodes.clone(), &e.edges.clone(), e.code)?;
        e.code = code;
        e.nodes = nodes;
    }
    for e in mesh.boundary.iter_mut() {
        let (code, nodes) = rewrite(&e.nodes.clone(), &e.edges.clone(), e.code)?;
        e.code = code;
        e.nodes = nodes;
    }
    // The global edge table is consumed by the upgrade.
    mesh.edges = None;
    for e in mesh.bulk.iter_mut() {
        e.edges.clear();
    }
    for e in mesh.boundary.iter_mut() {
        e.edges.clear();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::Elem;
    use elmer_sif::model::Value as SifValue;

    fn int_section(id: u32, name: &str, values: Vec<i32>) -> Section {
        Section {
            id,
            entries: vec![Entry {
                name: name.to_string(),
                value: SifValue::IntArray { values },
                partag: None,
                distribute: false,
            }],
        }
    }

    fn two_tri(dim: usize) -> Mesh {
        Mesh {
            dim,
            max_dim: dim,
            nodes: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            node_tags: vec![1, 2, 3, 4],
            bulk: vec![
                Elem { code: 303, body: 1, nodes: vec![0, 1, 2], tag: 1, edges: vec![] },
                Elem { code: 303, body: 1, nodes: vec![0, 2, 3], tag: 2, edges: vec![] },
            ],
            boundary: vec![],
            edges: None,
        }
    }

    fn sim_model(sections: Vec<Section>, sim: Option<Section>) -> Model {
        Model { bcs: sections, simulation: sim, ..Model::default() }
    }

    #[test]
    fn intersection_needs_exactly_two_parents() {
        let mut mesh = two_tri(2);
        let model = sim_model(vec![int_section(1, "intersection bc", vec![1, 2, 3])], None);
        let err = transform(&mut mesh, &model).unwrap_err();
        assert!(err.msg.contains("two parents"), "unexpected: {err:?}");
    }

    #[test]
    fn intersection_3d_is_later_lane() {
        let mut mesh = two_tri(3);
        let model = sim_model(vec![int_section(5, "intersection body", vec![1, 2])], None);
        let err = transform(&mut mesh, &model).unwrap_err();
        assert!(err.msg.contains("later lane"), "unexpected: {err:?}");
    }

    fn logical_section(name: &str) -> Section {
        Section {
            id: 0,
            entries: vec![Entry {
                name: name.to_string(),
                value: SifValue::Logical(true),
                partag: None,
                distribute: false,
            }],
        }
    }

    #[test]
    fn increase_rejects_nonlinear_type() {
        let mut mesh = two_tri(2);
        // Quadratic tet (already second order) can never be increased.
        mesh.bulk.push(Elem {
            code: 510,
            body: 1,
            nodes: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            tag: 3,
            edges: vec![],
        });
        let model = sim_model(vec![], Some(logical_section("increase element order")));
        let err = transform(&mut mesh, &model).unwrap_err();
        assert!(err.msg.contains("cannot increase element order for: 510"), "unexpected: {err:?}");
    }

    #[test]
    fn increase_rejects_wrong_node_count() {
        let mut mesh = two_tri(2);
        // A 303 must carry exactly 3 corner nodes.
        mesh.bulk[0].nodes.push(3);
        let model = sim_model(vec![], Some(logical_section("increase element order")));
        let err = transform(&mut mesh, &model).unwrap_err();
        assert!(err.msg.contains("got 4 nodes"), "unexpected: {err:?}");
    }

    fn tet_volume(m: &Mesh, n: &[u32]) -> f64 {
        let p = |i: usize| m.nodes[n[i] as usize];
        let d = |a: [f64; 3], b: [f64; 3]| [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let (u, v, w) = (d(p(0), p(1)), d(p(0), p(2)), d(p(0), p(3)));
        (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]))
            .abs()
            / 6.0
    }

    /// `Split Mesh Prisms`: a prism becomes three tetrahedra filling it,
    /// and its boundary quadrilateral two triangles, each with the new
    /// tetrahedron that contains it as its parent.
    #[test]
    fn split_mesh_prisms_fills_the_prism() {
        let mut mesh = Mesh {
            dim: 3,
            max_dim: 3,
            nodes: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
                [1.0, 0.0, 2.0],
                [0.0, 1.0, 2.0],
            ],
            node_tags: (1..=6).collect(),
            bulk: vec![Elem { code: 706, body: 1, nodes: (0..6).collect(), tag: 1, edges: vec![] }],
            boundary: vec![BoundElem {
                code: 404,
                constraint: 1,
                parents: [Some(0), None],
                nodes: vec![0, 1, 4, 3],
                tag: 1,
                edges: vec![],
            }],
            edges: None,
        };
        let model = Model { simulation: Some(logical_section("split mesh prisms")), ..Model::default() };
        split_mesh_quads(&mut mesh, &model).unwrap();
        assert_eq!(mesh.bulk.len(), 3);
        assert!(mesh.bulk.iter().all(|e| e.code == 504));
        let vol: f64 = mesh.bulk.iter().map(|e| tet_volume(&mesh, &e.nodes)).sum();
        assert!((vol - 1.0).abs() < 1e-12, "{vol}");
        assert_eq!(mesh.boundary.len(), 2);
        for b in &mesh.boundary {
            assert_eq!(b.code, 303);
            let p = b.parents[0].expect("parent");
            assert!(b.nodes.iter().all(|n| mesh.bulk[p].nodes.contains(n)));
        }
    }
}
