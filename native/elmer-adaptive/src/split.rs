//! `SplitMesh`, `SplitOneLevel`, `RGBRefinement` and `SetParents`
//! (Adaptive.F90): longest-edge bisection of marked linear triangles with
//! the red-green-blue closure, repeated while the children stay marked.

use elmer_core::{ElemId, ElementType, MeshView};

use crate::mesh::{Edges, Elem, OwnedMesh};
use crate::AdaptiveError;

/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// `SplitMesh`: marks each triangle `INT( MIN( MaxChange, 1 / Lambda ) )`
/// times, `Lambda = ( ErrorLimit / ErrorIndicator ) ** ( 1 / EhConvergence )`,
/// and splits until no child is marked. `None` when nothing is marked.
#[allow(clippy::too_many_arguments)]
pub fn split_mesh(
    ref_mesh: &dyn MeshView,
    adaptive_depth: usize,
    error_indicator: &[f64],
    error_limit: f64,
    h_convergence: &[f64],
    max_change: f64,
) -> Result<Option<OwnedMesh>, AdaptiveError> {
    let nb = ref_mesh.n_bulk();
    let mut splitted = vec![0i32; nb];
    let mut marked = 0;
    for i in 0..nb {
        let el = ref_mesh.element(ElemId(i as u32));
        if el.etype.0 != 303 {
            return Err(AdaptiveError::Fatal(
                "SplitMesh: Internal splitting implemented only for linear triangles.".into(),
            ));
        }
        let n = el.nodes.len();
        let eh = if adaptive_depth < 1 {
            1.0
        } else {
            let mut s = 0.0;
            for &k in el.nodes {
                s += h_convergence[k as usize];
            }
            s / n as f64
        };
        if error_indicator[i] > 100.0 * AEPS {
            let lambda = (error_limit / error_indicator[i]).powf(1.0 / eh);
            // `Splitted` is an INTEGER: the assignment truncates.
            splitted[i] = max_change.min(1.0 / lambda) as i32;
        }
        if splitted[i] > 0 {
            marked += 1;
        }
    }
    if marked == 0 {
        return Ok(None);
    }
    let (mut mesh, mut sp) = split_one_level(ref_mesh, &splitted);
    loop {
        if !sp.iter().any(|&s| s > 0) {
            break;
        }
        let (m2, s2) = split_one_level(&mesh, &sp);
        mesh = m2;
        sp = s2;
    }
    Ok(Some(mesh))
}

fn length(m: &dyn MeshView, e: [u32; 2]) -> f64 {
    let (a, b) = (m.coord(e[0]), m.coord(e[1]));
    ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt()
}

/// `RGBRefinement`: the longest edge of each marked triangle, then the
/// longest edge of every triangle with a marked edge until none is added.
fn rgb_refinement(m: &dyn MeshView, edges: &Edges, splitted: &[i32]) -> (Vec<bool>, usize) {
    let mut es = vec![false; edges.edges.len()];
    let longest = |i: usize| -> usize {
        let mut max_len = 0.0f64;
        let mut longest = 0usize;
        for &e in &edges.of_element[i] {
            let l = length(m, edges.edges[e].nodes);
            if l >= max_len {
                max_len = l;
                longest = e;
            }
        }
        longest
    };
    for i in 0..m.n_bulk() {
        if splitted[i] > 0 {
            es[longest(i)] = true;
        }
    }
    let mut marked = es.iter().filter(|&&x| x).count();
    loop {
        let mut hanging = 0;
        for i in 0..m.n_bulk() {
            let found = edges.of_element[i].iter().any(|&e| es[e]);
            let l = longest(i);
            if found && !es[l] {
                hanging += 1;
                es[l] = true;
            }
        }
        if hanging > 0 {
            marked += hanging;
        } else {
            break;
        }
    }
    (es, marked)
}

/// `SplitOneLevel`: splits every triangle along its marked edges (one,
/// two or three of them) and the boundary elements on split edges.
/// Returns the new mesh and its triangles' `Splitted` (the parent's
/// minus one).
fn split_one_level(m: &dyn MeshView, splitted: &[i32]) -> (OwnedMesh, Vec<i32>) {
    let edges = Edges::find(m);
    let (es, marked_edges) = rgb_refinement(m, &edges, splitted);
    let nn = m.n_nodes();
    let nb = m.n_bulk();
    let mut coords: Vec<[f64; 3]> = (0..nn as u32).map(|k| m.coord(k)).collect();
    coords.reserve(marked_edges);
    for (i, e) in edges.edges.iter().enumerate() {
        if es[i] {
            let (a, b) = (m.coord(e.nodes[0]), m.coord(e.nodes[1]));
            coords.push([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, 0.0]);
        }
    }
    // MarkedOrder: the new node of each split edge, numbered in element order.
    let mut marked_order = vec![0u32; edges.edges.len()];
    let mut k = 0u32;
    for i in 0..nb {
        for &e in &edges.of_element[i] {
            if es[e] && marked_order[e] == 0 {
                k += 1;
                marked_order[e] = k + nn as u32;
            }
        }
    }
    // MarkedOrder is 1-based in the Fortran; our node ids are 0-based.
    let mo = |e: usize| marked_order[e] - 1;
    let mut elems: Vec<Elem> = Vec::new();
    let mut children = vec![[usize::MAX; 4]; nb + m.n_boundary()];
    let mut new_split: Vec<i32> = Vec::new();
    let push = |elems: &mut Vec<Elem>, new_split: &mut Vec<i32>, nodes: [u32; 3], tag: i32, s: i32| -> usize {
        elems.push(Elem {
            etype: ElementType(303),
            nodes: nodes.to_vec(),
            tag,
            parents: [None, None],
        });
        new_split.push(s);
        elems.len() - 1
    };
    for i in 0..nb {
        let el = m.element(ElemId(i as u32));
        let rn: Vec<u32> = el.nodes.to_vec();
        let tag = el.tag;
        let ei = &edges.of_element[i];
        let marked = ei.iter().filter(|&&e| es[e]).count();
        let s = splitted[i] - 1;
        let en = |e: usize| edges.edges[e].nodes;
        let opposite = |e: usize| -> u32 {
            let n = en(e);
            *rn.iter().find(|&&x| x != n[0] && x != n[1]).expect("a triangle")
        };
        match marked {
            0 => {
                children[i][0] = push(&mut elems, &mut new_split, [rn[0], rn[1], rn[2]], tag, s);
            }
            1 => {
                let j = (0..3).find(|&j| es[ei[j]]).expect("one");
                let e = ei[j];
                let kn = opposite(e);
                children[i][0] = push(&mut elems, &mut new_split, [kn, en(e)[0], mo(e)], tag, s);
                children[i][1] = push(&mut elems, &mut new_split, [kn, mo(e), en(e)[1]], tag, s);
            }
            2 => {
                let j = (0..3).find(|&j| !es[ei[j]]).expect("one unsplit");
                let e = ei[j];
                let kn = opposite(e);
                let others: Vec<usize> = (0..3).filter(|&k| k != j).collect();
                children[i][0] =
                    push(&mut elems, &mut new_split, [kn, mo(ei[others[0]]), mo(ei[others[1]])], tag, s);
                let mut max_len = 0.0f64;
                let mut longest = 0usize;
                for &k in &others {
                    let l = length(m, en(ei[k]));
                    if l >= max_len {
                        max_len = l;
                        longest = k;
                    }
                }
                let ej = en(e);
                let ek = en(ei[longest]);
                let third = if ej[0] == ek[0] || ej[0] == ek[1] { ej[1] } else { ej[0] };
                children[i][1] =
                    push(&mut elems, &mut new_split, [mo(ei[others[0]]), mo(ei[others[1]]), third], tag, s);
                children[i][2] = push(&mut elems, &mut new_split, [ej[0], ej[1], mo(ei[longest])], tag, s);
            }
            _ => {
                children[i][0] = push(&mut elems, &mut new_split, [rn[0], mo(ei[0]), mo(ei[2])], tag, s);
                children[i][1] = push(&mut elems, &mut new_split, [rn[1], mo(ei[1]), mo(ei[0])], tag, s);
                children[i][2] = push(&mut elems, &mut new_split, [rn[2], mo(ei[2]), mo(ei[1])], tag, s);
                children[i][3] = push(&mut elems, &mut new_split, [mo(ei[0]), mo(ei[1]), mo(ei[2])], tag, s);
            }
        }
    }
    let n_bulk = elems.len();
    // SetParents: the child of the old parent that holds all the nodes.
    let parent_child = |elems: &[Elem], parent: usize, nodes: &[u32]| -> Option<ElemId> {
        for &c in &children[parent] {
            if c == usize::MAX {
                continue;
            }
            let cn = &elems[c].nodes;
            if nodes.iter().all(|n| cn.contains(n)) {
                return Some(ElemId(c as u32));
            }
        }
        None
    };
    for j in nb..nb + m.n_boundary() {
        let el = m.element(ElemId(j as u32));
        let rn: Vec<u32> = el.nodes.to_vec();
        let tag = el.tag;
        let parent = el.parents[0].or(el.parents[1]).map(|p| p.0 as usize);
        let Some(p) = parent else {
            elems.push(Elem {
                etype: el.etype,
                nodes: rn,
                tag,
                parents: [None, None],
            });
            continue;
        };
        let e = *edges.of_element[p]
            .iter()
            .find(|&&e| {
                let n = edges.edges[e].nodes;
                (n[0] == rn[0] && n[1] == rn[1]) || (n[1] == rn[0] && n[0] == rn[1])
            })
            .expect("the boundary element is an edge of its parent");
        let edge = &edges.edges[e];
        let add = |elems: &mut Vec<Elem>, nodes: Vec<u32>| {
            let left = edge.left.and_then(|l| parent_child(elems, l, &nodes));
            let right = edge.right.and_then(|r| parent_child(elems, r, &nodes));
            elems.push(Elem {
                etype: el.etype,
                nodes,
                tag,
                parents: [left, right],
            });
        };
        if es[e] {
            add(&mut elems, vec![rn[0], mo(e)]);
            add(&mut elems, vec![mo(e), rn[1]]);
        } else {
            add(&mut elems, rn.clone());
        }
    }
    (
        OwnedMesh {
            dim: m.dim(),
            coords,
            elems,
            n_bulk,
        },
        new_split,
    )
}
