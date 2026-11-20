//! ResultOutputSolver's STL format (`STL Format = True` or `Output Format =
//! stl`): `SaveSTLSurface` (SaveUtils.F90:1439), serially.
//!
//! The surface goes to `STL Filename` (default `mesh.stl`, relative to the
//! case directory) as ASCII STL, one `solid` per surface:
//!
//! * a 3D mesh with `STL Separate Bodies`: per body, the boundary faces with
//!   a parent in it, named by the body's `Name` (else `body<i>`);
//! * a 3D mesh with `STL Separate BCs`: per boundary condition, its faces,
//!   named by the BC's `Name` (else `boundary<i>`);
//! * otherwise the bulk triangles and quadrilaterals, as one `body`.
//!
//! Quadrilaterals become the triangles (1,2,3) and (1,3,4). Each facet gets
//! the face's normal and a vote: its agreement with the parent-oriented
//! normal (bodies), one (BCs), or with the direction from the origin to the
//! face's centre (bulk; the mesh centre the Fortran computes is reset to
//! zero before use). `OrientFacets` then joins facets that share an edge
//! (vertices merged on a grid of 1e-7 of the extent), flips facets so each
//! connected patch is consistently wound and turns a patch over when its
//! votes say so; a flip negates the normal and swaps the first and third
//! vertex.
//!
//! With `Reference Values`, the solver's norm is `AscBinCompareNorm` of
//! them against `[surfaces, facets, sum of normal components, sum of first
//! vertices' x, y, z]`. Refused: parallel runs and faces with more than
//! their corner nodes (the Fortran reads past its four-node buffer).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use elmer_core::{ElemId, MeshView, ParamMap, Params};
use elmer_fe::basis::{nodal_basis, nodal_first_derivatives};
use elmer_solver::{Model, Solver, SolverCtx, SolverError};

use crate::output::{asc_bin_compare_norm, is_result_output, reference_values};

const CALLER: &str = "SaveSTLSurface";
const FACET_VALS: usize = 13;

/// Whether a ResultOutputSolver section writes STL.
pub(crate) fn is_stl(params: &ParamMap) -> bool {
    params.logical("stl format") == Some(true) || params.string("output format").is_some_and(|f| f.eq_ignore_ascii_case("stl"))
}

/// One STL output section.
pub(crate) struct StlOutput {
    pub index: usize,
    file: PathBuf,
    reference: Option<[f64; 6]>,
}

/// The model's STL-only ResultOutputSolver sections.
pub(crate) fn prepare_stl_output(model: &Model, case_dir: &Path, parallel: bool) -> Result<Vec<StlOutput>, SolverError> {
    let mut out = Vec::new();
    for (index, params) in model.solvers.iter().enumerate() {
        if !is_result_output(params) || !is_stl(params) {
            continue;
        }
        if params.string("exec solver").is_some_and(|e| e.eq_ignore_ascii_case("never")) {
            continue;
        }
        let id = index + 1;
        if parallel {
            return Err(SolverError::Unsupported(format!("Solver {id}: STL output in parallel")));
        }
        let reference = match params.get("reference values") {
            None => None,
            Some(v) => Some(
                reference_values(v)
                    .ok_or_else(|| SolverError::Unsupported(format!("Solver {id}: Reference Values with fewer than 6 values")))?,
            ),
        };
        let name = params.string("stl filename").unwrap_or("mesh.stl");
        out.push(StlOutput { index, file: case_dir.join(name), reference });
    }
    Ok(out)
}

impl StlOutput {
    pub(crate) fn into_solver(self) -> StlSolver {
        StlSolver { out: self }
    }
}

pub(crate) struct StlSolver {
    out: StlOutput,
}

fn def(code: u16) -> Result<&'static elmer_fe::ElementTypeDef, SolverError> {
    elmer_fe::element_type(code).map_err(|e| SolverError::fatal(CALLER, e.to_string()))
}

/// `NormalVector(Element, Nodes)` of a triangle or quadrilateral at its
/// centre, `(1/3, 1/3)` or `(0, 0)`.
fn face_normal(code: u16, x: &[[f64; 3]]) -> Result<[f64; 3], SolverError> {
    let d = def(code)?;
    let (u, v) = if code / 100 == 3 { (1.0 / 3.0, 1.0 / 3.0) } else { (0.0, 0.0) };
    let mut dl = vec![[0.0f64; 3]; d.n_nodes];
    nodal_first_derivatives(d, u, v, 0.0, &mut dl);
    let dot = |c: usize, a: usize| (0..x.len()).fold(0.0, |s, k| s + dl[k][c] * x[k][a]);
    let (dxdu, dydu, dzdu) = (dot(0, 0), dot(0, 1), dot(0, 2));
    let (dxdv, dydv, dzdv) = (dot(1, 0), dot(1, 1), dot(1, 2));
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let det = 1.0 / (auu * avv - auv * auv).sqrt();
    Ok([(dydu * dzdv - dydv * dzdu) * det, (dxdv * dzdu - dxdu * dzdv) * det, (dxdu * dydv - dxdv * dydu) * det])
}

/// `NodalBasisFunctions` interpolation of `x` at `uvw`.
fn interpolate(code: u16, x: &[[f64; 3]], uvw: [f64; 3]) -> Result<[f64; 3], SolverError> {
    let d = def(code)?;
    let mut b = vec![0.0f64; d.n_nodes];
    nodal_basis(d, uvw[0], uvw[1], uvw[2], &mut b);
    Ok(std::array::from_fn(|a| (0..x.len()).fold(0.0, |s, k| s + b[k] * x[k][a])))
}

/// `NormalVector(Element, Nodes, Parent = Parent)`: the face normal turned
/// to point out of the parent (`CheckNormalDirectionParent`).
fn parent_normal(mesh: &dyn MeshView, code: u16, x: &[[f64; 3]], parent: ElemId) -> Result<[f64; 3], SolverError> {
    let n = face_normal(code, x)?;
    let centre = if code / 100 == 3 { [1.0 / 3.0, 1.0 / 3.0, 0.0] } else { [0.0; 3] };
    let c = interpolate(code, x, centre)?;
    let p = mesh.element(parent);
    let pc = match p.etype.0 / 100 {
        2 | 4 | 8 => [0.0; 3],
        3 | 7 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        5 => [0.25; 3],
        6 => [0.0, 0.0, 1.0 / 3.0],
        _ => return Err(SolverError::fatal("CheckNormalDirection", "Invalid elementcode for parent element!")),
    };
    let px: Vec<[f64; 3]> = p.nodes.iter().map(|&k| mesh.coord(k)).collect();
    let q = interpolate(p.etype.0, &px, pc)?;
    let d = [q[0] - c[0], q[1] - c[1], q[2] - c[2]];
    Ok(if d[0] * n[0] + d[1] * n[1] + d[2] * n[2] > 0.0 { n.map(|v| -v) } else { n })
}

/// `OrientFacets`: returns the number of flipped facets.
fn orient_facets(f: &mut [[f64; FACET_VALS]], keep_dir: bool) -> usize {
    let nf = f.len();
    if nf == 0 {
        return 0;
    }
    let ie = [[0usize, 1], [1, 2], [2, 0]];
    let mut cmin = [f[0][3], f[0][4], f[0][5]];
    let mut cmax = cmin;
    for fa in f.iter() {
        for j in 0..3 {
            for a in 0..3 {
                cmin[a] = cmin[a].min(fa[3 * j + 3 + a]);
                cmax[a] = cmax[a].max(fa[3 * j + 3 + a]);
            }
        }
    }
    let mut q = (0..3).map(|a| cmax[a] - cmin[a]).fold(f64::MIN, f64::max);
    if q <= 0.0 {
        q = 1.0;
    }
    q *= 1.0e-7;
    // Vertices merged on the grid, numbered in first-seen order.
    let mut verts: HashMap<[i64; 3], usize> = HashMap::new();
    let mut vid = vec![[0usize; 3]; nf];
    for (i, fa) in f.iter().enumerate() {
        for j in 0..3 {
            let iv: [i64; 3] = std::array::from_fn(|a| ((fa[3 * j + 3 + a] - cmin[a]) / q).round() as i64);
            let next = verts.len() + 1;
            vid[i][j] = *verts.entry(iv).or_insert(next);
        }
    }
    // Edge neighbours: the first two facets on an edge.
    let mut nbr = vec![[0usize; 3]; nf]; // 1-based, 0 = none
    let mut same = vec![[false; 3]; nf];
    let mut edges: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
    for i in 0..nf {
        for j in 0..3 {
            let (va, vb) = (vid[i][ie[j][0]], vid[i][ie[j][1]]);
            let key = (va.min(vb), va.max(vb));
            match edges.get(&key) {
                None => {
                    edges.insert(key, (i, j));
                }
                Some(&(k, l)) => {
                    if nbr[k][l] == 0 {
                        nbr[i][j] = k + 1;
                        nbr[k][l] = i + 1;
                        let da = vid[i][ie[j][0]] < vid[i][ie[j][1]];
                        let db = vid[k][ie[l][0]] < vid[k][ie[l][1]];
                        same[i][j] = da == db;
                        same[k][l] = same[i][j];
                    }
                }
            }
        }
    }
    let mut flip = vec![false; nf];
    let mut visited = vec![false; nf];
    for i in 0..nf {
        if visited[i] {
            continue;
        }
        let mut stack = vec![i];
        visited[i] = true;
        let mut comp = Vec::new();
        while let Some(k) = stack.pop() {
            comp.push(k);
            for j in 0..3 {
                let l = nbr[k][j];
                if l == 0 || visited[l - 1] {
                    continue;
                }
                flip[l - 1] = same[k][j] != flip[k];
                visited[l - 1] = true;
                stack.push(l - 1);
            }
        }
        let vote = comp.iter().fold(0.0, |s, &k| if flip[k] { s - f[k][FACET_VALS - 1] } else { s + f[k][FACET_VALS - 1] });
        if vote < 0.0 && !keep_dir {
            for &k in &comp {
                flip[k] = !flip[k];
            }
        }
    }
    let mut nrev = 0;
    for (fa, &fl) in f.iter_mut().zip(&flip) {
        if !fl {
            continue;
        }
        nrev += 1;
        for v in &mut fa[0..3] {
            *v = -*v;
        }
        for a in 0..3 {
            fa.swap(3 + a, 9 + a);
        }
    }
    nrev
}

impl Solver for StlSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let params = ctx.sim.model.solvers[self.out.index].clone();
        let mesh = ctx.mesh;
        let model = &ctx.sim.model;
        let nb = mesh.n_bulk();
        let three_d = nb > 0 && mesh.element(ElemId(0)).etype.0 > 500;
        let do_bodies = three_d && params.logical("stl separate bodies") == Some(true);
        let do_bcs = three_d && params.logical("stl separate bcs") == Some(true);
        let (range, min_body, max_body) = if do_bodies || do_bcs {
            let min = params.integer("stl min entity").unwrap_or(0).max(1) as usize;
            // `STL Max Entity` found: the Fortran takes the number of BCs.
            let max = match params.integer("stl max entity") {
                Some(_) => model.boundary_conditions.len(),
                None if do_bodies => model.bodies.len(),
                None => 0,
            };
            (nb..nb + mesh.n_boundary(), min, max)
        } else {
            (0..nb, 1, 1)
        };
        let mut text = String::new();
        let mut this = [0.0f64; 6];
        let mut body = min_body;
        loop {
            let name = if do_bodies {
                model.bodies.get(body - 1).and_then(|b| b.string("name")).map(str::to_string).unwrap_or(format!("body{body}"))
            } else if do_bcs {
                model.boundary_conditions.get(body - 1).and_then(|b| b.string("name")).map(str::to_string).unwrap_or(format!("boundary{body}"))
            } else {
                "body".to_string()
            };
            let mut facets: Vec<[f64; FACET_VALS]> = Vec::new();
            for t in range.clone() {
                let e = ElemId(t as u32);
                let el = mesh.element(e);
                let kmax = match el.etype.0 / 100 {
                    3 => 1,
                    4 => 2,
                    _ => continue,
                };
                let n = (el.etype.0 / 100) as usize;
                if el.nodes.len() != n {
                    return Err(SolverError::Unsupported(format!("{CALLER}: faces of element type {} (more than their corners)", el.etype.0)));
                }
                let x: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                let (normal, vote) = if do_bodies {
                    let is_body = |p: Option<ElemId>| p.filter(|&p| model.body_of_tag(mesh.element(p).tag) == Some(body));
                    let Some(parent) = is_body(el.parents[0]).or_else(|| is_body(el.parents[1])) else { continue };
                    let np = parent_normal(mesh, el.etype.0, &x, parent)?;
                    let n0 = face_normal(el.etype.0, &x)?;
                    (n0, n0[0] * np[0] + n0[1] * np[1] + n0[2] * np[2])
                } else if do_bcs {
                    if model.bc_of_tag(el.tag) != Some(body) {
                        continue;
                    }
                    let Some(parent) = el.parents[0].or(el.parents[1]) else {
                        return Err(SolverError::fatal(CALLER, "a boundary face without a parent"));
                    };
                    (parent_normal(mesh, el.etype.0, &x, parent)?, 1.0)
                } else {
                    let c: [f64; 3] = std::array::from_fn(|a| x.iter().map(|p| p[a]).sum::<f64>() / n as f64);
                    let n0 = face_normal(el.etype.0, &x)?;
                    (n0, n0[0] * c[0] + n0[1] * c[1] + n0[2] * c[2])
                };
                for k in 0..kmax {
                    let inds = if k == 0 { [0, 1, 2] } else { [0, 2, 3] };
                    let mut fa = [0.0f64; FACET_VALS];
                    fa[..3].copy_from_slice(&normal);
                    for (j, &i) in inds.iter().enumerate() {
                        fa[3 * j + 3..3 * j + 6].copy_from_slice(&x[i]);
                    }
                    fa[FACET_VALS - 1] = vote;
                    facets.push(fa);
                }
            }
            if !facets.is_empty() {
                orient_facets(&mut facets, false);
                let es = |v: f64| elmer_solver::fortran_fmt::es(v, 15, 6);
                text.push_str(&format!("solid {name}\n"));
                for fa in &facets {
                    text.push_str(&format!("  facet normal{}{}{}\n", es(fa[0]), es(fa[1]), es(fa[2])));
                    text.push_str("    outer loop\n");
                    for j in 0..3 {
                        let v = &fa[3 * j + 3..3 * j + 6];
                        text.push_str(&format!("      vertex{}{}{}\n", es(v[0]), es(v[1]), es(v[2])));
                    }
                    text.push_str("    endloop\n  endfacet\n");
                }
                text.push_str(&format!("endsolid {name}\n"));
                this[0] += 1.0;
                for fa in &facets {
                    this[1] += 1.0;
                    this[2] += fa[0] + fa[1] + fa[2];
                    this[3] += fa[3];
                    this[4] += fa[4];
                    this[5] += fa[5];
                }
            }
            body += 1;
            if body > max_body {
                break;
            }
        }
        std::fs::write(&self.out.file, text)
            .map_err(|e| SolverError::fatal(CALLER, format!("Could not open file: {}: {e}", self.out.file.display())))?;
        if let Some(reference) = self.out.reference {
            let nrm = asc_bin_compare_norm(&reference, &this);
            let id = ctx.solver_id();
            if let Some(slot) = ctx.sim.slots.iter_mut().find(|s| s.id == id) {
                match slot.variable {
                    Some(v) => ctx.sim.vars.by_id_mut(v).conv.norm = nrm,
                    None => slot.null_norm = nrm,
                }
            }
            ctx.info("ResultOutputSolver", &format!("Calculate Pseudonorm:{}", elmer_solver::fortran_fmt::es(nrm, 15, 6)), 4);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn facet(v: [[f64; 3]; 3], vote: f64) -> [f64; FACET_VALS] {
        let mut f = [0.0; FACET_VALS];
        f[2] = 1.0;
        for j in 0..3 {
            f[3 * j + 3..3 * j + 6].copy_from_slice(&v[j]);
        }
        f[12] = vote;
        f
    }

    #[test]
    fn inconsistent_neighbours_are_flipped_and_votes_turn_patches() {
        // Two triangles of the unit square sharing the diagonal a-c:
        // consistently wound (f), and the second reversed (g).
        let (a, b, c, d) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        let mut f = vec![facet([a, b, c], 1.0), facet([a, c, d], 1.0)];
        let mut g = vec![facet([a, b, c], 1.0), facet([d, c, a], 1.0)];
        assert_eq!(orient_facets(&mut f, false), 0);
        assert_eq!(orient_facets(&mut g, false), 1);
        // The flip negates the normal and swaps the first and third vertex.
        assert_eq!(g[1][2], -1.0);
        assert_eq!(&g[1][3..6], &a);
        assert_eq!(&g[1][9..12], &d);
        // A patch whose votes are negative is turned over as a whole.
        let mut h = vec![facet([a, b, c], -1.0), facet([a, c, d], -1.0)];
        assert_eq!(orient_facets(&mut h, false), 2);
        assert_eq!(orient_facets(&mut [facet([a, b, c], -1.0)], true), 0);
    }
}
