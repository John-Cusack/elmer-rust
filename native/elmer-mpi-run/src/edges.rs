//! The parallel numbering of a 2D mesh's edges (`SParEdgeNumbering`,
//! SParIterComm.F90:750-1409), for the order of `Mesh % Edges`.
//!
//! A partitioned `ElmerSolver_mpi` numbers the edges globally while loading
//! the mesh and then sorts its edge table by those numbers, so a DG solver,
//! which numbers its dofs over the edges (`InitialPermutation`) and
//! assembles the jumps edge by edge, sees the edges in global order rather
//! than in the order `FindMeshEdges2D` created them. An edge is an
//! interface edge when both its nodes are shared and it lacks one of its
//! two parents here; it belongs to the greatest partition both nodes list.
//! Each rank numbers the edges it owns after the lower ranks' counts; the
//! owners send the numbers of their interface edges to the other holders,
//! who find them by the global ids of the edge nodes. An edge the owner
//! sends but the receiver does not take as an interface edge is dropped
//! from the owner's holders; an edge left without a number passes its
//! ownership to the next partition of its list and the round repeats.

use elmer_core::ElemId;
use elmer_dist::Exchange;

use crate::mesh::ParMesh;

/// One edge of the local first-encounter list.
struct Edgen {
    n: u64,
    interface: bool,
    /// Holder ranks, the owner first; a holder that did not take an edge
    /// is marked `-rank - 1` until it is dropped.
    neighbours: Vec<i64>,
    /// The global ids of the two edge nodes.
    parent: [u64; 2],
}

/// The global number of each edge of `edges` (`dg_edges`' list, in
/// `FindMeshEdges2D` order). Collective.
pub fn edge_numbers(ex: &dyn Exchange, mesh: &ParMesh, edges: &[(ElemId, Option<ElemId>, Vec<u32>)]) -> Result<Vec<u64>, String> {
    edge_numbers_with_holders(ex, mesh, edges).map(|(n, _)| n)
}

/// [`edge_numbers`] with each edge's holder ranks, the owner first
/// (`Mesh % ParallelInfo % EdgeNeighbourList`); a private edge lists this
/// rank alone. Collective.
#[allow(clippy::type_complexity)]
pub fn edge_numbers_with_holders(
    ex: &dyn Exchange,
    mesh: &ParMesh,
    edges: &[(ElemId, Option<ElemId>, Vec<u32>)],
) -> Result<(Vec<u64>, Vec<Vec<usize>>), String> {
    let me = ex.rank() as i64;
    let np = ex.size();
    let gid = mesh.node_gids();
    let parts = mesh.node_parts();
    let shared = |n: u32| parts[n as usize].len() > 1;
    let mut e: Vec<Edgen> = edges
        .iter()
        .map(|(_, right, nodes)| {
            let (p, q) = (nodes[0], nodes[1]);
            let mut ed = Edgen { n: 0, interface: false, neighbours: vec![me], parent: [gid[p as usize], gid[q as usize]] };
            if shared(p) && shared(q) && right.is_none() {
                // The node of the greater global id lists first.
                let (a, b) = if gid[p as usize] > gid[q as usize] { (p, q) } else { (q, p) };
                let l2 = &parts[b as usize];
                let mut common: Vec<i64> =
                    parts[a as usize].iter().filter(|r| l2.contains(r)).map(|&r| r as i64).collect();
                if common.len() >= 2 {
                    let q = *common.iter().max().expect("two entries");
                    let at = common.iter().position(|&r| r == q).expect("present");
                    common.swap(0, at);
                    ed.interface = true;
                    ed.neighbours = common;
                }
            }
            ed
        })
        .collect();
    let matches = |ed: &Edgen, a: u64, b: u64| (ed.parent[0] == a && ed.parent[1] == b) || (ed.parent[0] == b && ed.parent[1] == a);
    let enc = |v: i64| v as u64;
    let dec = |v: u64| v as i64;

    for _round in 0..64 {
        // Owned edges after the lower ranks' counts.
        let mut counts = vec![0.0f64; np];
        counts[me as usize] = e.iter().filter(|x| x.neighbours[0] == me).count() as f64;
        ex.sum_f64(&mut counts);
        let mut k = counts[..me as usize].iter().sum::<f64>() as u64;
        for x in e.iter_mut() {
            x.n = 0;
            if x.neighbours[0] == me {
                k += 1;
                x.n = k;
            }
        }
        let iedges: Vec<usize> = (0..e.len()).filter(|&i| e[i].interface).collect();

        // Owners send (number, node ids) to every other holder.
        let mut send = vec![Vec::new(); np];
        for (r, out) in send.iter_mut().enumerate() {
            for &l in &iedges {
                if e[l].neighbours[0] != me {
                    continue;
                }
                for &m in &e[l].neighbours[1..] {
                    if m == r as i64 {
                        out.extend([e[l].n, e[l].parent[0], e[l].parent[1]]);
                    }
                }
            }
        }
        let got = ex.all_to_all_u64(&send);
        let mut missing = vec![Vec::new(); np];
        for (proc, msg) in got.iter().enumerate() {
            for t in msg.as_chunks::<3>().0 {
                match iedges.iter().find(|&&k| matches(&e[k], t[1], t[2])) {
                    Some(&k) => e[k].n = t[0],
                    None => missing[proc].push(t[0]),
                }
            }
        }
        // The holders that did not take an edge.
        let back = ex.all_to_all_u64(&missing);
        for (proc, msg) in back.iter().enumerate() {
            for &j in &iedges {
                if msg.contains(&e[j].n) {
                    if let Some(m) = e[j].neighbours[1..].iter().position(|&r| r == proc as i64) {
                        e[j].neighbours[m + 1] = -(proc as i64) - 1;
                    }
                }
            }
        }
        for &j in &iedges {
            if e[j].n == 0 || e[j].neighbours[0] != me {
                continue;
            }
            e[j].neighbours.retain(|&r| r >= 0);
            if e[j].neighbours.len() == 1 {
                e[j].interface = false;
            }
        }
        let iedges: Vec<usize> = (0..e.len()).filter(|&i| e[i].interface).collect();

        // Owners send each holder the final holder lists.
        let mut send = vec![Vec::new(); np];
        for (r, out) in send.iter_mut().enumerate() {
            if r as i64 == me {
                continue;
            }
            for &j in &iedges {
                if e[j].n == 0 || e[j].neighbours[0] != me || !e[j].neighbours.contains(&(r as i64)) {
                    continue;
                }
                out.push(enc(-(e[j].n as i64)));
                out.extend(e[j].neighbours[1..].iter().map(|&m| enc(m)));
            }
        }
        let got = ex.all_to_all_u64(&send);
        for (proc, msg) in got.iter().enumerate() {
            let g: Vec<i64> = msg.iter().map(|&v| dec(v)).collect();
            let mut j = 0;
            while j < g.len() {
                let k = -g[j];
                if let Some(&l) = iedges.iter().find(|&&l| e[l].n != 0 && e[l].n as i64 == k) {
                    let m = (j + 1..g.len()).find(|&m| g[m] < 0).unwrap_or(g.len());
                    let mut nb = vec![proc as i64];
                    nb.extend_from_slice(&g[j + 1..m]);
                    e[l].neighbours = nb;
                    j = m;
                } else {
                    j += 1;
                }
            }
        }

        // An edge left without a number goes to the next partition listed.
        let mut unnumbered = [0.0f64];
        for x in e.iter_mut().filter(|x| x.n == 0) {
            unnumbered[0] += 1.0;
            let k = x.neighbours.len();
            x.neighbours.remove(0);
            if x.neighbours.is_empty() {
                x.neighbours.push(me);
            }
            if k <= 2 {
                x.interface = false;
            }
        }
        ex.sum_f64(&mut unnumbered);
        if unnumbered[0] == 0.0 {
            let holders = e
                .iter()
                .map(|x| if x.interface { x.neighbours.iter().filter(|&&r| r >= 0).map(|&r| r as usize).collect() } else { vec![me as usize] })
                .collect();
            return Ok((e.iter().map(|x| x.n).collect(), holders));
        }
    }
    Err("SParEdgeNumbering: edges left unnumbered after 64 rounds".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{Elem, Nodes};
    use elmer_core::{Communicator, ElementType};
    use elmer_dist::ThreadComm;

    /// A 2x3 grid of nodes, global ids 1 2 3 / 4 5 6, cut in two triangles
    /// per square: rank 0 owns 1-2-5 and 1-5-4, rank 1 owns 2-3-6 and 2-6-5,
    /// and each holds the other's triangle on the 2-5 seam as halo.
    fn rank_mesh(rank: usize) -> ParMesh {
        let (gids, elems): (Vec<u64>, Vec<[u64; 3]>) = if rank == 0 {
            (vec![1, 2, 4, 5, 6], vec![[1, 2, 5], [1, 5, 4], [2, 6, 5]])
        } else {
            (vec![1, 2, 3, 5, 6], vec![[2, 3, 6], [2, 6, 5], [1, 2, 5]])
        };
        let parts = gids
            .iter()
            .map(|&g| match g {
                4 => vec![0],
                3 => vec![1],
                6 => vec![1, 0],
                _ => vec![0, 1],
            })
            .collect();
        let local = |g: u64| gids.iter().position(|&x| x == g).expect("held") as u32;
        let elems = elems
            .iter()
            .map(|e| Elem { etype: ElementType(303), nodes: e.iter().map(|&g| local(g)).collect(), tag: 1, parents: [None, None] })
            .collect();
        let nodes = Nodes { coords: vec![[0.0; 3]; gids.len()], gid: gids.clone(), parts };
        ParMesh::new(rank, 2, 2, nodes, elems, 3)
    }

    #[test]
    fn every_edge_gets_a_distinct_global_number() {
        let out = ThreadComm::run(2, |ex| {
            let mesh = rank_mesh(ex.rank());
            let edges = elmer_solver::dg_edges(&mesh).expect("triangles");
            let n = edge_numbers(ex, &mesh, &edges).expect("numbered");
            let gid = mesh.node_gids();
            let named: Vec<([u64; 2], u64)> = edges
                .iter()
                .zip(&n)
                .map(|((_, _, nodes), &k)| {
                    let mut g = [gid[nodes[0] as usize], gid[nodes[1] as usize]];
                    g.sort_unstable();
                    (g, k)
                })
                .collect();
            named
        });
        for named in &out {
            let mut k: Vec<u64> = named.iter().map(|x| x.1).collect();
            assert!(k.iter().all(|&x| x > 0), "{named:?}");
            k.sort_unstable();
            k.dedup();
            assert_eq!(k.len(), named.len(), "{named:?}");
        }
        // 2-6 has one parent on rank 0, but both on its first owner, rank 1,
        // which therefore never sends it: rank 0 numbers it itself.
        let on0 = out[0].iter().find(|x| x.0 == [2, 6]).expect("edge 2-6").1;
        let owned0 = out[0].len() as u64;
        assert!(on0 <= owned0, "{:?}", out[0]);
    }
}
