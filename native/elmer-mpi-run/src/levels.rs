//! `Mesh Levels` on a partitioned mesh: `SplitMeshEqual` on every rank's
//! piece (`elmer-mesh-levels`), then the parallel information of the new
//! nodes, which Elmer derives in `SParEdgeNumbering` /
//! `SParFaceNumbering` (SParIterComm.F90:750, :1414).
//!
//! A new node sits at the midpoint of an old edge, the centre of an old
//! quadrilateral face, or the centre of an old element. It can be shared
//! only when all its parent nodes are shared, and only with ranks that hold
//! all of them; such a rank shares it exactly when it has the same edge or
//! face, which both sides recognise by the sorted global ids of the parent
//! nodes. The lowest holding rank owns the node. New global ids continue
//! after the largest old one: owners number their new nodes in local order,
//! rank by rank, and send the ids to the other holders.

use std::collections::HashMap;

use elmer_core::{ParamMap, Params};
use elmer_dist::Exchange;
use elmer_mesh_levels::{split_mesh_equal, topology::find_mesh_edges, Mesh, MeshLevelsError};
use elmer_solver::SolverError;

use crate::mesh::{Elem, Nodes, ParMesh};

fn levels_err(e: MeshLevelsError) -> SolverError {
    match e {
        MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshLevelsError::Invalid(m) => SolverError::fatal("LoadModel", m),
    }
}

/// The `Mesh Levels` keywords `elmer_mesh_levels::mesh_levels` refuses,
/// checked before any rank splits. Returns the number of splits.
pub fn splits(simulation: &ParamMap, bcs: &[ParamMap]) -> Result<usize, SolverError> {
    let levels = simulation.integer("mesh levels").unwrap_or(1);
    if levels <= 1 {
        return Ok(0);
    }
    for key in ["split mesh quads", "increase element order", "constant stencil"] {
        if simulation.logical(key) == Some(true) {
            return Err(SolverError::Unsupported(format!("{key} on a refined mesh")));
        }
    }
    if simulation.logical("mesh keep grading") == Some(true) {
        return Err(SolverError::Unsupported("Mesh Keep Grading".into()));
    }
    if bcs.iter().any(|bc| bc.0.keys().any(|k| k.starts_with("follow"))) {
        return Err(SolverError::Unsupported("Follow * Boundary with Mesh Levels".into()));
    }
    Ok(levels as usize - 1)
}

/// The old nodes each new node of `split_mesh_equal(old)` averages, in the
/// split's numbering: edges, quadrilateral faces, 404/808 centres.
fn new_node_parents(old: &Mesh) -> Result<Vec<Vec<u32>>, SolverError> {
    let topo = find_mesh_edges(old).map_err(levels_err)?;
    let mut parents: Vec<Vec<u32>> = topo.edges.iter().map(|e| e.to_vec()).collect();
    parents.extend(topo.faces.iter().filter(|f| f.len() == 4).cloned());
    parents.extend(old.bulk().iter().filter(|e| matches!(e.etype.family(), 4 | 8)).map(|e| e.nodes.clone()));
    Ok(parents)
}

/// Message layout of a key list: `len, gid_1 .. gid_len` per key.
fn push_key(buf: &mut Vec<u64>, key: &[u64]) {
    buf.push(key.len() as u64);
    buf.extend_from_slice(key);
}

fn keys(buf: &[u64]) -> Vec<&[u64]> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < buf.len() {
        let n = buf[i] as usize;
        out.push(&buf[i + 1..i + 1 + n]);
        i += 1 + n;
    }
    out
}

/// The rank-local split and the parents of its new nodes, checked against
/// the coordinates the split gave them.
fn local_split(old: &Mesh, mesh: &ParMesh) -> Result<(Vec<Vec<u32>>, Mesh), SolverError> {
    let parents = new_node_parents(old)?;
    let new = split_mesh_equal(mesh).map_err(levels_err)?;
    let n_old = old.coords.len();
    if new.coords.len() != n_old + parents.len() {
        return Err(SolverError::fatal(
            "SplitMeshEqual",
            format!("{} new nodes, {} parent sets", new.coords.len() - n_old, parents.len()),
        ));
    }
    // The split averages the parents the same way (sum from zero, divide).
    for (k, p) in parents.iter().enumerate() {
        let c = new.coords[n_old + k];
        let mean = |d: usize| p.iter().fold(0.0, |s, &n| s + old.coords[n as usize][d]) / p.len() as f64;
        if c != [mean(0), mean(1), mean(2)] {
            return Err(SolverError::fatal("SplitMeshEqual", format!("new node {} is not its parents' centre", n_old + k + 1)));
        }
    }
    Ok((parents, new))
}

/// One `SplitMeshEqual` of a partitioned mesh, with global ids, holders and
/// owners for the new nodes. Collective.
pub fn split(mesh: &ParMesh, ex: &dyn Exchange) -> Result<ParMesh, SolverError> {
    let (me, np) = (ex.rank(), ex.size());
    let old = Mesh::from_view(mesh);
    let local = local_split(&old, mesh);
    // Every rank learns whether any rank failed before the collectives.
    let mut failed = [if local.is_err() { 1.0 } else { 0.0 }];
    ex.max_f64(&mut failed);
    let (parents, new) = match local {
        Err(e) => return Err(e),
        Ok(_) if failed[0] != 0.0 => return Err(SolverError::fatal("SplitMeshEqual", "another rank could not split its partition")),
        Ok(v) => v,
    };

    let gids = mesh.node_gids();
    let holders = mesh.node_parts();
    // Ranks other than this one that hold every parent, and the sorted
    // parent gids, of each new node.
    let candidates: Vec<Vec<usize>> = parents
        .iter()
        .map(|p| {
            let mut c: Vec<usize> = holders[p[0] as usize].iter().copied().filter(|&r| r != me).collect();
            for &n in &p[1..] {
                c.retain(|r| holders[n as usize].contains(r));
            }
            c.sort_unstable();
            c
        })
        .collect();
    let key_of = |k: usize| {
        let mut key: Vec<u64> = parents[k].iter().map(|&n| gids[n as usize]).collect();
        key.sort_unstable();
        key
    };
    let mut by_key: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut ask: Vec<Vec<u64>> = vec![Vec::new(); np];
    let mut asked_for: Vec<Vec<usize>> = vec![Vec::new(); np];
    for (k, c) in candidates.iter().enumerate() {
        if c.is_empty() {
            continue;
        }
        let key = key_of(k);
        for &r in c {
            push_key(&mut ask[r], &key);
            asked_for[r].push(k);
        }
        by_key.insert(key, k);
    }

    // Does the asking rank's edge or face exist here too?
    let asked = ex.all_to_all_u64(&ask);
    let answers: Vec<Vec<u64>> = asked
        .iter()
        .enumerate()
        .map(|(r, buf)| {
            keys(buf)
                .into_iter()
                .map(|key| u64::from(by_key.get(key).is_some_and(|&k| candidates[k].contains(&r))))
                .collect()
        })
        .collect();
    let got = ex.all_to_all_u64(&answers);
    let mut node_parts: Vec<Vec<usize>> = vec![vec![me]; parents.len()];
    for (r, yes) in got.iter().enumerate() {
        for (&k, &y) in asked_for[r].iter().zip(yes) {
            if y == 1 {
                node_parts[k].push(r);
            }
        }
    }
    for p in &mut node_parts {
        p.sort_unstable();
    }

    // Owners number their new nodes after the largest old id.
    let mut top = [gids.iter().copied().max().unwrap_or(0) as f64];
    ex.max_f64(&mut top);
    let mut counts = vec![0.0; np];
    counts[me] = node_parts.iter().filter(|p| p[0] == me).count() as f64;
    ex.sum_f64(&mut counts);
    let mut next = top[0] as u64 + counts[..me].iter().sum::<f64>() as u64 + 1;
    let mut new_gid = vec![0u64; parents.len()];
    for (k, p) in node_parts.iter().enumerate() {
        if p[0] == me {
            new_gid[k] = next;
            next += 1;
        }
    }
    // The other holders ask the owner for the id.
    let mut req: Vec<Vec<u64>> = vec![Vec::new(); np];
    let mut req_for: Vec<Vec<usize>> = vec![Vec::new(); np];
    for (k, p) in node_parts.iter().enumerate() {
        if p[0] != me {
            push_key(&mut req[p[0]], &key_of(k));
            req_for[p[0]].push(k);
        }
    }
    let reqs = ex.all_to_all_u64(&req);
    let mut missing = None;
    let replies: Vec<Vec<u64>> = reqs
        .iter()
        .map(|buf| {
            keys(buf)
                .into_iter()
                .map(|key| match by_key.get(key) {
                    Some(&k) if node_parts[k][0] == me => new_gid[k],
                    _ => {
                        missing.get_or_insert_with(|| format!("asked for the id of a new node with parents {key:?}"));
                        0
                    }
                })
                .collect()
        })
        .collect();
    let ids = ex.all_to_all_u64(&replies);
    for (r, v) in ids.iter().enumerate() {
        for (&k, &g) in req_for[r].iter().zip(v) {
            new_gid[k] = g;
        }
    }
    let mut bad = [if missing.is_some() { 1.0 } else { 0.0 }];
    ex.max_f64(&mut bad);
    if bad[0] != 0.0 {
        return Err(SolverError::fatal("SParEdgeNumbering", missing.unwrap_or_else(|| "another rank lost a new node".into())));
    }

    let nodes = Nodes {
        coords: new.coords.clone(),
        gid: gids.iter().copied().chain(new_gid).collect(),
        parts: holders.iter().cloned().chain(node_parts).collect(),
    };
    let elems = new
        .elements
        .iter()
        .map(|e| Elem { etype: e.etype, nodes: e.nodes.clone(), tag: e.tag, parents: e.parents })
        .collect();
    Ok(ParMesh::new(me, mesh.n_parts(), new.dim, nodes, elems, new.n_bulk))
}
