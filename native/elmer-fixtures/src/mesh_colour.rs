//! `MeshColour/MeshColour.F90`: a unit test of the mesh graph utilities as
//! a solver. `MeshColour_init` makes the variable a hidden global
//! (`Variable = -nooutput -global meshcolour_var`, given a `Norm Variable
//! Index`); `MeshColour` builds the dual graph of the bulk elements
//! (`ElmerMeshToDualGraph`: elements are adjacent when they share a node),
//! colours it (`ElmerGraphColour`), gathers the colour index lists
//! (`ElmerColouringToGraph`) and checks both (`GraphColourVerify`,
//! `GraphColourListVerify`). The variable and its norm become `1 + nerror`,
//! the number of errors the checks found.
//!
//! The METIS cross-check (`HAVE_METIS`) is not part of the port: it compares
//! the dual graph with `METIS_MeshToDual`'s and only adds errors.

use elmer_core::{ElemId, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "MeshColour";

/// `MeshColour_init`: with a positive `Norm Variable Index` and no
/// `Variable`, the variable is `-nooutput -global meshcolour_var`.
pub fn init_section(p: &mut ParamMap) {
    if p.integer("norm variable index").unwrap_or(0) > 0 && p.get("variable").is_none() {
        p.insert("Variable", Value::Str("-nooutput -global meshcolour_var".into()));
    }
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MeshColour))
}

struct MeshColour;

/// A graph in compressed rows: the neighbours of `v` are
/// `ind[ptr[v]..ptr[v + 1]]`.
struct Graph {
    ptr: Vec<usize>,
    ind: Vec<usize>,
}

/// `ElmerMeshToDualGraph`: the bulk elements, adjacent when they share a
/// node, each element's neighbours in increasing order and without itself.
fn dual_graph(elements: &[Vec<u32>], n_nodes: usize) -> Graph {
    let mut node_elems: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (e, nodes) in elements.iter().enumerate() {
        for &k in nodes {
            node_elems[k as usize].push(e);
        }
    }
    let mut ptr = vec![0];
    let mut ind = Vec::new();
    for (e, nodes) in elements.iter().enumerate() {
        let mut adj: Vec<usize> = nodes.iter().flat_map(|&k| node_elems[k as usize].iter().copied()).filter(|&f| f != e).collect();
        adj.sort_unstable();
        adj.dedup();
        ind.extend(adj);
        ptr.push(ind.len());
    }
    Graph { ptr, ind }
}

/// `ElmerGraphColour`: greedy colouring in vertex order, each vertex taking
/// the smallest colour (1-based) none of its coloured neighbours has.
/// Returns the colours and their count.
fn colour(g: &Graph) -> (Vec<usize>, usize) {
    let n = g.ptr.len() - 1;
    let mut c = vec![0usize; n];
    let mut nc = 0;
    let mut used = Vec::new();
    for v in 0..n {
        used.clear();
        used.resize(nc + 2, false);
        for &w in &g.ind[g.ptr[v]..g.ptr[v + 1]] {
            used[c[w]] = true;
        }
        let k = (1..).find(|&k| !used[k]).unwrap_or(1);
        c[v] = k;
        nc = nc.max(k);
    }
    (c, nc)
}

/// `ElmerColouringToGraph`: the vertices of each colour, in vertex order.
fn colour_lists(c: &[usize], nc: usize) -> Graph {
    let mut ptr = vec![0usize; nc + 1];
    for &k in c {
        ptr[k] += 1;
    }
    for k in 0..nc {
        ptr[k + 1] += ptr[k];
    }
    let mut ind = vec![0usize; c.len()];
    let mut at = ptr.clone();
    for (v, &k) in c.iter().enumerate() {
        ind[at[k - 1]] = v;
        at[k - 1] += 1;
    }
    Graph { ptr, ind }
}

/// `GraphColourVerify`: each colour in `1..=nc`, no edge between vertices
/// of one colour.
fn verify_colouring(g: &Graph, c: &[usize], nc: usize) -> usize {
    let mut err = 0;
    for v in 0..c.len() {
        if c[v] < 1 || c[v] > nc {
            err += 1;
            continue;
        }
        err += g.ind[g.ptr[v]..g.ptr[v + 1]].iter().filter(|&&w| c[w] == c[v]).count();
    }
    err
}

/// `GraphColourListVerify`: every vertex in the list of its colour.
fn verify_lists(c: &[usize], nc: usize, lists: &Graph) -> usize {
    let mut err = 0;
    for k in 1..=nc {
        for &v in &lists.ind[lists.ptr[k - 1]..lists.ptr[k]] {
            if c[v] != k {
                err += 1;
            }
        }
    }
    if lists.ptr[nc] != c.len() {
        err += 1;
    }
    err
}

impl Solver for MeshColour {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let elements: Vec<Vec<u32>> = (0..ctx.mesh.n_bulk()).map(|t| ctx.mesh.element(ElemId(t as u32)).nodes.to_vec()).collect();
        let g = dual_graph(&elements, ctx.mesh.n_nodes());
        let (c, nc) = colour(&g);
        let mut nerror = verify_colouring(&g, &c, nc);
        let lists = colour_lists(&c, nc);
        nerror += verify_lists(&c, nc, &lists);
        let var = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
        var.values.iter_mut().for_each(|x| *x = (1 + nerror) as f64);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_strip_of_quads_takes_two_colours() {
        // Three quads in a row: 0-1 and 1-2 share nodes, 0-2 do not.
        let els = vec![vec![0, 1, 5, 4], vec![1, 2, 6, 5], vec![2, 3, 7, 6]];
        let g = dual_graph(&els, 8);
        assert_eq!(g.ind, vec![1, 0, 2, 1]);
        let (c, nc) = colour(&g);
        assert_eq!((c.clone(), nc), (vec![1, 2, 1], 2));
        let lists = colour_lists(&c, nc);
        assert_eq!(lists.ind, vec![0, 2, 1]);
        assert_eq!(verify_colouring(&g, &c, nc) + verify_lists(&c, nc, &lists), 0);
        assert_eq!(verify_colouring(&g, &[1, 1, 1], 1), 4);
    }
}
