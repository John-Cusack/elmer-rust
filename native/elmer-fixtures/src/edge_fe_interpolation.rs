//! `EdgeFEInterpolationTest.F90` (fem/src/modules) as the
//! `SecondOrderEdgeElement2D_BCs` test runs it (`Procedure =
//! "EdgeFEInterpolationTest" "BestApproximationSolver"`): the H(curl)
//! best approximation `(u, v) + MatPar (curl u, curl v)` of a manufactured
//! field in the Piola edge element space of the first kind, with the
//! boundary edge dofs fixed by `<Variable> {e} 1..3` Dirichlet conditions
//! (`VectorElementEdgeDOFs`, the tangential moments), and with `Error
//! Computation` (the default) the relative L2 (or curl) error printed as
//! `Error Norm = ...`.
//!
//! [`init_section`] is `BestApproximationSolver_Init0`: the `Element` of
//! the basis `EdgeElementStyle` picks. The Gauss rule is `GaussPoints(
//! Element, PReferenceElement=Piola, EdgeBasisDegree=k)`: of degree 2 the p
//! triangle's 6 points. Ported: second-order (`Quadratic Approximation`)
//! Piola elements on triangles. Refused: the second kind, cubic and
//! first-order elements, 3D, `Simplicial Mesh`.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe_hp::edge::{edge_element_info, p_reference_nodal_basis};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "BestApproximationSolver";

/// `EdgeElementStyle(SolverParams, PiolaVersion, SecondFamily, BasisDegree)`.
fn edge_style(p: &ParamMap) -> (bool, bool, usize) {
    let quadratic = p.logical("quadratic approximation") == Some(true);
    let second = p.logical("second kind basis") == Some(true);
    let piola = p.logical("use piola transform") == Some(true) || quadratic || second;
    (piola, second, if quadratic { 2 } else { 1 })
}

/// `BestApproximationSolver_Init0`.
pub fn init_section(p: &mut ParamMap) {
    if p.get("element").is_some() {
        return;
    }
    let (piola, second, k) = edge_style(p);
    let e = match (k, second) {
        (2, true) => "n:0 e:3 -tri b:3 -tri_face b:3",
        (2, false) => "n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2",
        (_, true) => "n:0 e:2",
        (_, false) if piola => "n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2",
        _ => "n:0 e:1",
    };
    p.insert("Element", Value::Str(e.into()));
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(BestApproximation))
}

struct BestApproximation;

/// The test field and its curl.
fn field(mode: i32, xq: f64, yq: f64) -> Result<([f64; 3], [f64; 3]), SolverError> {
    match mode {
        5 => Ok(([-xq * (yq * yq), yq * (xq * xq), 0.0], [0.0, 0.0, 4.0 * xq * yq])),
        m => Err(SolverError::Unsupported(format!("{CALLER}: Test Mode {m}"))),
    }
}

struct Point {
    det: f64,
    s: f64,
    basis: Vec<f64>,
    e: Vec<[f64; 3]>,
    c: Vec<[f64; 3]>,
}

fn points(ctx: &SolverCtx<'_>, e: ElemId, nd: usize) -> Result<(Vec<Point>, Vec<[f64; 3]>), SolverError> {
    let el = ctx.mesh.element(e);
    let etype = el.etype;
    if etype.0 != 303 {
        return Err(SolverError::Unsupported(format!("{CALLER} on element type {}", etype.0)));
    }
    let rule = elmer_pdofs::quadrature::gauss_points_triangle_p(6).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
    let mut out = Vec::with_capacity(rule.len());
    for t in 0..rule.len() {
        let uvw = [rule.u[t], rule.v[t], rule.w[t]];
        let ev = edge_element_info(etype, &xyz, &gindex, uvw, 2, ctx.mesh.dim()).map_err(|m| SolverError::fatal(CALLER, m))?;
        if ev.dofs != nd {
            return Err(SolverError::fatal(CALLER, format!("{} edge basis functions for {nd} element dofs", ev.dofs)));
        }
        out.push(Point { det: ev.det, s: rule.s[t], basis: p_reference_nodal_basis(etype, uvw), e: ev.edge, c: ev.curl });
    }
    Ok((out, xyz))
}

fn sum(a: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for x in a {
        s += x;
    }
    s
}

/// `DefaultDirichletBCs` for `<name> {e}` on the boundary edges of
/// triangles (`PickActiveFace` gives the mesh edge, which has its only
/// triangle's edge-map node order, `FindMeshEdges2D`), and
/// `VectorElementEdgeDOFs(BC, Edge, n, Parent, np, Name//' {e}', Work,
/// EDOFs=2, SecondFamily=.FALSE., BasisDegree=2)`: `Integral(1) = int
/// (L + VL.e)`, `Integral(2) = int (L + VL.e)(Basis(2) - Basis(1))`, the
/// first negated when the edge runs against the global node order. The
/// values go to the parent's dofs of that edge.
fn edge_dirichlet(ctx: &SolverCtx<'_>, name: &str) -> Result<Vec<(u32, f64)>, SolverError> {
    const EDGES: [[usize; 2]; 3] = [[0, 1], [1, 2], [2, 0]];
    let prefix = format!("{name} {{e}}");
    let condition = format!("{name} condition");
    let nb = ctx.mesh.n_bulk();
    let active: std::collections::HashSet<u32> = ctx.active_elements().map(|e| e.0).collect();
    let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
    let rule = elmer_fe::quadrature::gauss_points(ElementType(202), None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let mut fixed = Vec::new();
    let mut slots = Vec::new();
    let mut ev = BasisEval::default();
    for &e in &boundary {
        let Some(bc) = ctx.bc_params(e) else { continue };
        if !bc.0.keys().any(|k| k.starts_with(&prefix)) {
            continue;
        }
        if let Some(c) = ctx.real_at_nodes(bc, &condition, e)? {
            if c.iter().sum::<f64>() > 0.0 {
                continue;
            }
        }
        let el = ctx.mesh.element(e);
        if el.etype.0 != 202 {
            return Err(SolverError::Unsupported(format!("{CALLER}: {{e}} conditions on element type {}", el.etype.0)));
        }
        let Some(parent) = el.parents.iter().flatten().copied().find(|p| (p.0 as usize) < nb && active.contains(&p.0)) else {
            continue;
        };
        let pel = ctx.mesh.element(parent);
        if pel.etype.0 != 303 {
            return Err(SolverError::Unsupported(format!("{CALLER}: {{e}} conditions on a parent of type {}", pel.etype.0)));
        }
        let bnodes = el.nodes.to_vec();
        let Some(l) = EDGES.iter().position(|&[j, k]| {
            (pel.nodes[j] == bnodes[0] && pel.nodes[k] == bnodes[1]) || (pel.nodes[j] == bnodes[1] && pel.nodes[k] == bnodes[0])
        }) else {
            continue;
        };
        let [j, k] = EDGES[l];
        let (nj, nk) = (pel.nodes[j], pel.nodes[k]);
        // The condition at the mesh edge's nodes (nj, nk).
        let order: [usize; 2] = if bnodes[0] == nj { [0, 1] } else { [1, 0] };
        let at = |key: &str| -> Result<[f64; 2], SolverError> {
            Ok(match ctx.real_at_nodes(bc, key, e)? {
                Some(v) => [v[order[0]], v[order[1]]],
                None => [0.0; 2],
            })
        };
        let load = at(&prefix)?;
        let vload = [at(&format!("{prefix} 1"))?, at(&format!("{prefix} 2"))?, at(&format!("{prefix} 3"))?];
        let (pj, pk) = (ctx.mesh.coord(nj), ctx.mesh.coord(nk));
        let d = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let dir = d.map(|x| x / len);
        let xyz = [pj, pk];
        let mut integral = [0.0f64; 2];
        for t in 0..rule.len() {
            elmer_fe::LagrangeBasis::new(3)
                .eval(ElementType(202), &xyz, [rule.u[t], 0.0, 0.0], &mut ev)
                .map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let s = rule.s[t] * ev.det_j;
            let (b0, b1) = (ev.basis[0], ev.basis[1]);
            let lv = load[0] * b0 + load[1] * b1;
            let vl: [f64; 3] = std::array::from_fn(|c| vload[c][0] * b0 + vload[c][1] * b1);
            let ve = (vl[0] * dir[0] + vl[1] * dir[1]) + vl[2] * dir[2];
            integral[0] += s * (lv + ve);
            let v = b1 - b0;
            integral[1] += s * (lv + ve) * v;
        }
        if nk < nj {
            integral[0] = -integral[0];
        }
        ctx.element_slots(parent, &mut slots);
        for (jj, val) in integral.iter().enumerate() {
            let s = slots[2 * l + jj];
            if s != NO_DOF {
                fixed.push((s, *val));
            }
        }
    }
    Ok(fixed)
}

impl Solver for BestApproximation {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the edge element dofs (no DofsBackend)"));
        }
        let p = ctx.params().clone();
        let (piola, second, k) = edge_style(&p);
        if !piola || second || k != 2 {
            return Err(SolverError::Unsupported(format!("{CALLER}: only second-order Piola elements of the first kind")));
        }
        if p.logical("simplicial mesh") == Some(true) || ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported(format!("{CALLER}: 3D or Simplicial Mesh")));
        }
        let error_estimation = p.logical("error computation").unwrap_or(true);
        let curl_norm = p.logical("compute curl error").unwrap_or(false);
        let mode = p.integer("test mode").unwrap_or(2);
        let name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        let dim = ctx.mesh.dim();
        ctx.default_initialize()?;
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
        let mut slots = Vec::new();
        for &e in &active {
            ctx.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, format!("element {} lacks dofs", e.0)));
            }
            let nd = slots.len();
            let (pts, xyz) = points(ctx, e, nd)?;
            let n = xyz.len();
            let acoef = match ctx.body_params(e).and_then(|b| b.material) {
                Some(m) => ctx.real_at_nodes(m, "material param", e)?.unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            let mut stiff = vec![0.0; nd * nd];
            let mut force = vec![0.0; nd];
            for pt in &pts {
                let (dj, s) = (pt.det, pt.s);
                let xq = sum((0..n).map(|i| xyz[i][0] * pt.basis[i]));
                let yq = sum((0..n).map(|i| xyz[i][1] * pt.basis[i]));
                let mat_par = sum((0..n).map(|i| acoef[i] * pt.basis[i]));
                for pp in 0..nd {
                    for q in 0..nd {
                        let t1 = 1.0 * sum((0..dim).map(|c| pt.e[q][c] * pt.e[pp][c])) * dj * s;
                        let t2 = mat_par * sum((0..3).map(|c| pt.c[q][c] * pt.c[pp][c])) * dj * s;
                        let a = &mut stiff[pp * nd + q];
                        *a = (*a + t1) + t2;
                    }
                    let (e, c) = (&pt.e[pp], &pt.c[pp]);
                    force[pp] = match mode {
                        5 => {
                            force[pp]
                                + (-xq * (yq * yq)) * e[0] * dj * s
                                + (yq * (xq * xq)) * e[1] * dj * s
                                + mat_par * (0.0) * c[0] * dj * s
                                + mat_par * (0.0) * c[1] * dj * s
                                + mat_par * (4.0 * xq * yq) * c[2] * dj * s
                        }
                        m => return Err(SolverError::Unsupported(format!("{CALLER}: Test Mode {m}"))),
                    };
                }
            }
            let ls = LocalSystem { n: nd, stiff, mass: None, damp: None, force };
            ctx.default_update(&ls, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        let fixed = edge_dirichlet(ctx, &name)?;
        ctx.default_dirichlet_after(&fixed)?;
        ctx.default_solve()?;

        if error_estimation {
            let vals = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.values.to_vec();
            let (mut err, mut sol_norm) = (0.0f64, 0.0f64);
            for &e in &active {
                ctx.element_slots(e, &mut slots);
                let nd = slots.len();
                let load: Vec<f64> = slots.iter().map(|&s| vals[s as usize]).collect();
                let (pts, xyz) = points(ctx, e, nd)?;
                let n = xyz.len();
                for pt in &pts {
                    let xq = sum((0..n).map(|i| xyz[i][0] * pt.basis[i]));
                    let yq = sum((0..n).map(|i| xyz[i][1] * pt.basis[i]));
                    let mut u = [0.0; 3];
                    for i in 0..dim {
                        u[i] = sum((0..nd).map(|j| load[j] * pt.e[j][i]));
                    }
                    let rotu: [f64; 3] = std::array::from_fn(|c| sum((0..nd).map(|j| load[j] * pt.c[j][c])));
                    let (sol, rotsol) = field(mode, xq, yq)?;
                    let e3: [f64; 3] = std::array::from_fn(|c| sol[c] - u[c]);
                    let rote: [f64; 3] = std::array::from_fn(|c| rotsol[c] - rotu[c]);
                    let ss = |v: &[f64; 3]| sum(v.iter().map(|x| x * x));
                    if curl_norm {
                        sol_norm += ss(&rotsol) * pt.det * pt.s;
                        err += ss(&rote) * pt.det * pt.s;
                    } else {
                        sol_norm += ss(&sol) * pt.det * pt.s;
                        err += ss(&e3) * pt.det * pt.s;
                    }
                }
            }
            let msg = format!("Error Norm = {}", crate::poisson_ml::fortran_e(err.sqrt() / sol_norm.sqrt(), 16, 8));
            ctx.info(CALLER, &msg, 1);
        }
        Ok(())
    }
}
