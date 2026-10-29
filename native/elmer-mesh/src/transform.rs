//! Mesh-load transforms, replicating `MeshLoad.F90` / `MeshTagging.F90`.
//!
//! Implemented (S1 scope, report §5.3):
//! - `MapCoordinates` (coordinate mapping + scaling, `MeshLoad.F90:558`),
//! - `MapBodiesAndBCs` (Target Bodies/Boundaries, `MeshLoad.F90:230`),
//! - `MapInternalExternalBCs` (`MeshLoad.F90:466`),
//! - `CreateIntersectionBCs` in 2D, both modes (`MeshTagging.F90:52`),
//! - `IncreaseElementOrder` (`MeshBasics.F90:3564`).
//!
//! Explicitly out of scope for this PR: 3D intersection BCs (need face
//! tables), `SplitMeshQuads`, `TagBCsUsingRule`, curved-boundary
//! following, and all parallel handling. Those are capability refusals
//! in a later lane, not silent skips.

use crate::mesh::{BoundElem, Mesh, MeshError};
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
        return Err(fatal(
            "3D intersection BCs need face tables (later lane)".to_string(),
        ));
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
            // Edges shared by bulk elements whose bodies are exactly the
            // two parents. (2D BulkMode; 3D is refused above.)
            let mut edge_users: std::collections::HashMap<u32, Vec<usize>> =
                std::collections::HashMap::new();
            for (bi, b) in mesh.bulk.iter().enumerate() {
                for &eid in &b.edges {
                    edge_users.entry(eid).or_default().push(bi);
                }
            }
            let table = mesh.edges.as_ref().cloned().unwrap_or_default();
            for (eid, users) in &edge_users {
                if users.len() != 2 {
                    continue;
                }
                let (b0, b1) = (mesh.bulk[users[0]].body, mesh.bulk[users[1]].body);
                if b0 == b1 {
                    continue;
                }
                if !((b0 == row.pa || b0 == row.pb) && (b1 == row.pa || b1 == row.pb)) {
                    continue;
                }
                let (n1, _) = table.edges[*eid as usize];
                new_elems.push(New {
                    node: n1,
                    tag: row.new_tag,
                    left_bpos: 0,
                    right_bpos: 0,
                    left_bulk: Some(users[0]),
                    right_bulk: Some(users[1]),
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
        let parents = if let (Some(a), Some(b)) = (ne.left_bulk, ne.right_bulk) {
            [Some(a), Some(b)]
        } else {
            [Some(n_bulk + ne.left_bpos), Some(n_bulk + ne.right_bpos)]
        };
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
