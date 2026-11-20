//! `SaveSTLSurface` (SaveUtils.F90:1439): ResultOutputSolver's `STL Format`.
//! The triangles and quadrilaterals (split in two) of the surface are
//! written as ASCII STL, one `solid` per body with `STL Separate Bodies`
//! (a 3D mesh's boundary elements by the body of the parent on that side).
//! `OrientFacets` winds every connected patch consistently and turns it so
//! that it agrees with its elements' votes.
//!
//! Ported: serial runs, `STL Separate Bodies`, the whole-mesh form, `STL
//! Filename`, and `Reference Values` (the six checksums of the written
//! facets, compared by `AscBinCompareNorm` into the solver's norm).
//! Refused: `STL Separate BCs` and `STL Min/Max Entity`.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::PathBuf;

use elmer_core::{ElemId, MeshView, ParamMap, Params};
use elmer_solver::{Solver, SolverCtx, SolverError};

/// The values kept per facet: normal, three vertices, vote (`FacetVals`).
type Facet = [f64; 13];

/// gfortran `ESw.d`: `d.ddddddE+xx`, right-aligned, the sign of a
/// negative zero kept.
fn es(x: f64, w: usize, d: usize) -> String {
    let s = format!("{:.*e}", d, x);
    let (mant, exp) = s.split_once('e').expect("exponent");
    let e: i32 = exp.parse().expect("exponent");
    let body = if e.abs() >= 100 {
        format!("{mant}{}{:03}", if e < 0 { '-' } else { '+' }, e.abs())
    } else {
        format!("{mant}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    };
    format!("{body:>w$}")
}

/// `NormalVector( Boundary, Nodes [, Parent] )` of a surface element (no
/// orientation check without `parent`).
fn normal(ctx: &SolverCtx<'_>, e: ElemId, parent: Option<ElemId>) -> Result<[f64; 3], SolverError> {
    let el = ctx.mesh.element(e);
    let def = elmer_fe::element_type(el.etype.0).map_err(|m| SolverError::fatal("NormalVector", m.to_string()))?;
    let n = def.n_nodes;
    let xyz: Vec<[f64; 3]> = el.nodes[..n].iter().map(|&k| ctx.mesh.coord(k)).collect();
    let (u, v) = if el.etype.0 / 100 == 3 { (1.0 / 3.0, 1.0 / 3.0) } else { (0.0, 0.0) };
    let mut dl = vec![[0.0f64; 3]; n];
    elmer_fe::basis::nodal_first_derivatives(def, u, v, 0.0, &mut dl);
    let dot = |k: usize, c: usize| {
        let mut s = 0.0;
        for i in 0..n {
            s += dl[i][k] * xyz[i][c];
        }
        s
    };
    let (dxdu, dydu, dzdu) = (dot(0, 0), dot(0, 1), dot(0, 2));
    let (dxdv, dydv, dzdv) = (dot(1, 0), dot(1, 1), dot(1, 2));
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let det = 1.0 / (auu * avv - auv * auv).sqrt();
    let mut nrm = [
        (dydu * dzdv - dydv * dzdu) * det,
        (dxdv * dzdu - dxdu * dzdv) * det,
        (dxdu * dydv - dxdv * dydu) * det,
    ];
    let Some(p) = parent else { return Ok(nrm) };
    // `CheckNormalDirectionParent`: the point at the element's centre,
    // the parent's centre.
    let mut b = vec![0.0f64; n];
    elmer_fe::basis::nodal_basis(def, u, v, 0.0, &mut b);
    let mut c = [0.0f64; 3];
    for (k, ck) in c.iter_mut().enumerate() {
        for i in 0..n {
            *ck += b[i] * xyz[i][k];
        }
    }
    let pel = ctx.mesh.element(p);
    let pdef = elmer_fe::element_type(pel.etype.0).map_err(|m| SolverError::fatal("NormalVector", m.to_string()))?;
    let (pu, pv, pw) = match pel.etype.0 / 100 {
        2 | 4 | 8 => (0.0, 0.0, 0.0),
        3 | 7 => (1.0 / 3.0, 1.0 / 3.0, 0.0),
        5 => (1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0),
        6 => (0.0, 0.0, 1.0 / 3.0),
        _ => return Err(SolverError::fatal("CheckNormalDirection", "Invalid elementcode for parent element!")),
    };
    let pn = pdef.n_nodes;
    let mut pb = vec![0.0f64; pn];
    elmer_fe::basis::nodal_basis(pdef, pu, pv, pw, &mut pb);
    let mut d = [0.0f64; 3];
    for (k, dk) in d.iter_mut().enumerate() {
        let mut s = 0.0;
        for i in 0..pn {
            s += pb[i] * ctx.mesh.coord(pel.nodes[i])[k];
        }
        *dk = s - c[k];
    }
    if d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2] > 0.0 {
        nrm = [-nrm[0], -nrm[1], -nrm[2]];
    }
    Ok(nrm)
}

/// `OrientFacets`: consistent winding per connected patch, each patch
/// turned to agree with the sum of its votes. Returns the facets reversed.
fn orient_facets(f: &mut [Facet]) -> usize {
    let nf = f.len();
    if nf == 0 {
        return 0;
    }
    let ie = [[0usize, 1], [1, 2], [2, 0]];
    let mut cmin = [f[0][3], f[0][4], f[0][5]];
    let mut cmax = cmin;
    for fa in f.iter() {
        for j in 0..3 {
            for k in 0..3 {
                let c = fa[3 * j + 3 + k];
                cmin[k] = cmin[k].min(c);
                cmax[k] = cmax[k].max(c);
            }
        }
    }
    let mut q = (0..3).map(|k| cmax[k] - cmin[k]).fold(f64::MIN, f64::max);
    if q <= 0.0 {
        q = 1.0;
    }
    q *= 1.0e-7;
    // Vertices by their rounded position, numbered in first-seen order.
    let mut verts: HashMap<[i64; 3], usize> = HashMap::new();
    let mut vid = vec![[0usize; 3]; nf];
    for (i, fa) in f.iter().enumerate() {
        for j in 0..3 {
            let key: [i64; 3] = std::array::from_fn(|k| ((fa[3 * j + 3 + k] - cmin[k]) / q).round() as i64);
            let n = verts.len() + 1;
            vid[i][j] = *verts.entry(key).or_insert(n);
        }
    }
    let mut nbr = vec![[usize::MAX; 3]; nf];
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
                    if nbr[k][l] == usize::MAX {
                        nbr[i][j] = k;
                        nbr[k][l] = i;
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
                if l == usize::MAX || visited[l] {
                    continue;
                }
                flip[l] = same[k][j] != flip[k];
                visited[l] = true;
                stack.push(l);
            }
        }
        let mut vote = 0.0;
        for &k in &comp {
            if flip[k] {
                vote -= f[k][12];
            } else {
                vote += f[k][12];
            }
        }
        if vote < 0.0 {
            for &k in &comp {
                flip[k] = !flip[k];
            }
        }
    }
    let mut rev = 0;
    for (fa, &fl) in f.iter_mut().zip(&flip) {
        if !fl {
            continue;
        }
        rev += 1;
        for k in 0..3 {
            fa[k] = -fa[k];
            fa.swap(3 + k, 9 + k);
        }
    }
    rev
}

/// The ResultOutputSolver section with `STL Format`.
pub struct StlSolver {
    pub params: ParamMap,
    pub case_dir: PathBuf,
}

impl StlSolver {
    pub fn check(params: &ParamMap) -> Result<(), SolverError> {
        for k in ["stl separate bcs", "stl min entity", "stl max entity"] {
            if params.get(k).is_some() {
                return Err(SolverError::Unsupported(format!("STL output: '{k}'")));
            }
        }
        Ok(())
    }
}

impl Solver for StlSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("STL output in parallel".into()));
        }
        let file = self.params.string("stl filename").unwrap_or("mesh.stl").to_string();
        let nb = ctx.mesh.n_bulk();
        let three_d = nb > 0 && ctx.mesh.element(ElemId(0)).etype.0 > 500;
        let bodies = three_d && self.params.logical("stl separate bodies") == Some(true);
        let n_bodies = ctx.sim.model.bodies.len();
        let (range, groups): (std::ops::Range<usize>, Vec<usize>) = if bodies {
            (nb..nb + ctx.mesh.n_boundary(), (1..=n_bodies).collect())
        } else {
            (0..nb, vec![1])
        };
        let body_of = |e: ElemId| ctx.sim.model.body_of_tag(ctx.mesh.element(e).tag);
        let mut out = String::new();
        // `This Values`: solids, facets, sum of the normals, sums of the
        // first vertex's coordinates.
        let mut this = [0.0f64; 6];
        for &body in &groups {
            let name = if bodies {
                match ctx.sim.model.bodies[body - 1].string("name") {
                    Some(s) => s.to_lowercase(),
                    None => format!("body{body}"),
                }
            } else {
                "body".to_string()
            };
            let mut facets: Vec<Facet> = Vec::new();
            for i in range.clone() {
                let e = ElemId(i as u32);
                let el = ctx.mesh.element(e);
                let n = (el.etype.0 / 100) as usize;
                let kmax = match n {
                    3 => 1,
                    4 => 2,
                    _ => continue,
                };
                let xyz: Vec<[f64; 3]> = el.nodes[..n].iter().map(|&k| ctx.mesh.coord(k)).collect();
                let (nrm, vote) = if bodies {
                    let [l, r] = el.parents;
                    let parent = match (l, r) {
                        (Some(l), _) if body_of(l) == Some(body) => l,
                        (_, Some(r)) if body_of(r) == Some(body) => r,
                        _ => continue,
                    };
                    let np = normal(ctx, e, Some(parent))?;
                    let nrm = normal(ctx, e, None)?;
                    (nrm, nrm[0] * np[0] + nrm[1] * np[1] + nrm[2] * np[2])
                } else {
                    // `MeshCenter` is zero by then (SaveUtils.F90).
                    let mut c = [0.0f64; 3];
                    for (k, ck) in c.iter_mut().enumerate() {
                        let mut s = 0.0;
                        for p in &xyz {
                            s += p[k];
                        }
                        *ck = s / n as f64;
                    }
                    let nrm = normal(ctx, e, None)?;
                    (nrm, nrm[0] * c[0] + nrm[1] * c[1] + nrm[2] * c[2])
                };
                for k in 0..kmax {
                    let ix = if k == 0 { [0, 1, 2] } else { [0, 2, 3] };
                    let mut fa: Facet = [0.0; 13];
                    fa[..3].copy_from_slice(&nrm);
                    for (j, &v) in ix.iter().enumerate() {
                        fa[3 * j + 3..3 * j + 6].copy_from_slice(&xyz[v]);
                    }
                    fa[12] = vote;
                    facets.push(fa);
                }
            }
            if facets.is_empty() {
                continue;
            }
            orient_facets(&mut facets);
            this[0] += 1.0;
            for fa in &facets {
                this[1] += 1.0;
                let mut s = 0.0;
                for k in 0..3 {
                    s += fa[k];
                }
                this[2] += s;
                this[3] += fa[3];
                this[4] += fa[4];
                this[5] += fa[5];
            }
            let _ = writeln!(out, "solid {name}");
            for fa in &facets {
                let _ = writeln!(out, "  facet normal{}{}{}", es(fa[0], 15, 6), es(fa[1], 15, 6), es(fa[2], 15, 6));
                let _ = writeln!(out, "    outer loop");
                for j in 0..3 {
                    let v = &fa[3 * j + 3..3 * j + 6];
                    let _ = writeln!(out, "      vertex{}{}{}", es(v[0], 15, 6), es(v[1], 15, 6), es(v[2], 15, 6));
                }
                let _ = writeln!(out, "    endloop");
                let _ = writeln!(out, "  endfacet");
            }
            let _ = writeln!(out, "endsolid {name}");
        }
        std::fs::write(self.case_dir.join(&file), out).map_err(|e| SolverError::fatal("SaveSTLSurface", format!("{file}: {e}")))?;
        // ResultOutputSolver.F90:351: `Solver % Variable % Norm =
        // AscBinCompareNorm(RefResults, This Values)`.
        if let Some(v) = self.params.get("reference values") {
            let reference = crate::output::reference_values(v)
                .ok_or_else(|| SolverError::Unsupported("STL output: Reference Values with fewer than 6 values".into()))?;
            let nrm = crate::output::asc_bin_compare_norm(&reference, &this);
            let id = ctx.solver_id();
            if let Some(slot) = ctx.sim.slots.iter_mut().find(|s| s.id == id) {
                match slot.variable {
                    Some(v) => ctx.sim.vars.by_id_mut(v).conv.norm = nrm,
                    None => slot.null_norm = nrm,
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es_matches_gfortran() {
        assert_eq!(es(0.0, 15, 6), "   0.000000E+00");
        assert_eq!(es(-0.0, 15, 6), "  -0.000000E+00");
        assert_eq!(es(-1.0, 15, 6), "  -1.000000E+00");
        assert_eq!(es(0.25, 15, 6), "   2.500000E-01");
    }
}

#[allow(dead_code)]
fn _mesh(_: &dyn MeshView) {}
