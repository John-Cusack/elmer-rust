//! The test-local `Hrot.f90` of `HelmholtzEdge` (`Procedure = "Hrot"
//! "HrotSolver"`): the 2D vector Helmholtz equation `rot rot u + A u = f`
//! with complex `A = (A 1, A 2)`, on lowest-order Whitney (edge) elements of
//! triangles (`Element = "n:0 e:1"`) and a complex variable of two real
//! components (`Variable DOFs = 2`; the solver sets `Solver % Matrix %
//! Complex`, here `Linear System Complex` from [`init_section`]).
//!
//! The basis is the plugin's own `EdgeBasis`: `W_k = B_i grad B_j - B_j
//! grad B_i` on the triangle's edges `(1,2), (2,3), (3,1)`, negated when the
//! edge runs against the global node order. Boundary elements add `g` and
//! `B` (per unit length) on their parent's edge; `MyDirichletBCs` sets the
//! edge DOF of `Tangent Component 1|2` to `+-SUM(data)/2 * h` (the sign by
//! the mesh edge's node order) by zeroing the row. The solve is
//! `DefaultSolve`. `WriteResults` (the `ratkaisu.ep` ElmerPost file the
//! plugin writes beside the results) is not written: the field is in the
//! result file.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::metric::MetricStatus;
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "HrotSolver";

/// `Solver % Matrix % Complex = .TRUE.`
pub fn init_section(p: &mut ParamMap) {
    p.insert("Linear System Complex", Value::Logical(true));
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Hrot))
}

struct Hrot;

const EDGES: [[usize; 2]; 3] = [[0, 1], [1, 2], [2, 0]];

/// `SUM(Basis(1:n) * c(1:n))`.
fn interp(basis: &[f64], c: &[f64]) -> f64 {
    let mut s = 0.0;
    for (b, x) in basis.iter().zip(c) {
        s += b * x;
    }
    s
}

/// The plugin's `EdgeBasis(WhitneyBasis, RotWhitneyBasis, u, v, w, Element)`.
fn edge_basis(be: &BasisEval, g: &[u32]) -> ([[f64; 2]; 3], [f64; 3]) {
    let (b, d) = (&be.basis, &be.dbasis_dx);
    let mut w = [[0.0; 2]; 3];
    let mut r = [0.0; 3];
    for (k, &[i, j]) in EDGES.iter().enumerate() {
        w[k][0] = b[i] * d[j][0] - b[j] * d[i][0];
        w[k][1] = b[i] * d[j][1] - b[j] * d[i][1];
        r[k] = 2.0 * (d[i][1] * d[j][0] - d[j][1] * d[i][0]);
        if g[j] < g[i] {
            w[k] = [-w[k][0], -w[k][1]];
            r[k] = -r[k];
        }
    }
    (w, r)
}

/// The complex local system as the interleaved real one.
fn to_real(n: usize, re: &[f64], im: &[f64], fre: &[f64], fim: &[f64]) -> LocalSystem {
    let m = 2 * n;
    let mut ls = LocalSystem { n: m, stiff: vec![0.0; m * m], mass: None, damp: None, force: vec![0.0; m] };
    for p in 0..n {
        for q in 0..n {
            let (a, b) = (re[p * n + q], im[p * n + q]);
            ls.stiff[2 * p * m + 2 * q] = a;
            ls.stiff[2 * p * m + 2 * q + 1] = -b;
            ls.stiff[(2 * p + 1) * m + 2 * q] = b;
            ls.stiff[(2 * p + 1) * m + 2 * q + 1] = a;
        }
        ls.force[2 * p] = fre[p];
        ls.force[2 * p + 1] = fim[p];
    }
    ls
}

fn at_nodes(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    Ok(match p {
        Some(p) => ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]),
        None => vec![0.0; n],
    })
}

/// The parent's local edge (`EDGES`) with both nodes of the boundary
/// element `b`.
fn parent_edge(pnodes: &[u32], b: &[u32]) -> Option<usize> {
    EDGES.iter().position(|&[i, j]| b.contains(&pnodes[i]) && b.contains(&pnodes[j]))
}

impl Solver for Hrot {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the edge element DOFs (no DofsBackend)"));
        }
        if ctx.params().logical("eigen analysis") == Some(true) {
            return Err(SolverError::Unsupported(format!("{CALLER}: Eigen Analysis")));
        }
        if ctx.primary().map(|v| v.dofs) != Some(2) {
            return Err(SolverError::fatal(CALLER, "the variable needs two components"));
        }
        let name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        let cdim = ctx.mesh.dim();
        let nbulk = ctx.mesh.n_bulk();
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
        let def = elmer_fe::element_type(303).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
        let mut status = MetricStatus::default();
        let mut slots = Vec::new();
        for &e in &active {
            let el = ctx.mesh.element(e);
            if el.etype.0 != 303 {
                return Err(SolverError::Unsupported(format!("{CALLER} on element type {}", el.etype.0)));
            }
            ctx.element_slots(e, &mut slots);
            if slots.len() != 3 || slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, "a triangle without its three edge DOFs"));
            }
            let g: Vec<u32> = el.nodes.to_vec();
            let xyz: Vec<[f64; 3]> = g.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let bp = ctx.body_params(e);
            let load1 = at_nodes(ctx, bp.as_ref().and_then(|b| b.body_force), "f 1", e, 3)?;
            let load2 = at_nodes(ctx, bp.as_ref().and_then(|b| b.body_force), "f 2", e, 3)?;
            let a1 = at_nodes(ctx, bp.as_ref().and_then(|b| b.material), "a 1", e, 3)?;
            let a2 = at_nodes(ctx, bp.as_ref().and_then(|b| b.material), "a 2", e, 3)?;
            let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let (mut sre, mut sim) = (vec![0.0; 9], vec![0.0; 9]);
            let mut fre = vec![0.0; 3];
            let fim = vec![0.0; 3];
            for t in 0..rule.len() {
                let mut be = BasisEval::default();
                elmer_fe::info::element_info(def, &xyz, [rule.u[t], rule.v[t], rule.w[t]], cdim, &mut status, &mut be)
                    .map_err(|m| SolverError::fatal("ElementInfo", m.to_string()))?;
                let (dj, s) = (be.det_j, rule.s[t]);
                let basis = &be.basis[..3];
                let (l1, l2) = (interp(basis, &load1), interp(basis, &load2));
                let (ca1, ca2) = (interp(basis, &a1), interp(basis, &a2));
                let (w, rot) = edge_basis(&be, &g);
                for i in 0..3 {
                    fre[i] += (l1 * w[i][0] + l2 * w[i][1]) * dj * s;
                    for j in 0..3 {
                        sre[i * 3 + j] += rot[i] * rot[j] * dj * s;
                        let ww = w[i][0] * w[j][0] + w[i][1] * w[j][1];
                        sre[i * 3 + j] += ca1 * ww * dj * s;
                        sim[i * 3 + j] += ca2 * ww * dj * s;
                    }
                }
            }
            ctx.default_update(&to_real(3, &sre, &sim, &fre, &fim), e)?;
        }
        ctx.default_finish_bulk_assembly()?;

        // Robin terms `B` and loads `g` on the parent's boundary edge.
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let coord = |k: u32| ctx.mesh.coord(k);
        let edge_h = |a: u32, b: u32| {
            let (pa, pb) = (coord(a), coord(b));
            let h = (pb[0] - pa[0]) * (pb[0] - pa[0]) + (pb[1] - pa[1]) * (pb[1] - pa[1]);
            h.sqrt()
        };
        for &b in &boundary {
            let el = ctx.mesh.element(b);
            if el.etype.0 / 100 == 1 {
                continue;
            }
            let Some(bc) = ctx.bc_params(b) else { continue };
            let keys = ["g 1", "g 2", "b 1", "b 2"];
            if !keys.iter().any(|k| bc.get(k).is_some()) {
                continue;
            }
            let Some(parent) = el.parents[0].or(el.parents[1]) else { continue };
            let n = el.nodes.len();
            let v: Vec<Vec<f64>> = keys.iter().map(|k| at_nodes(ctx, Some(bc), k, b, n)).collect::<Result<_, _>>()?;
            let pnodes = ctx.mesh.element(parent).nodes.to_vec();
            let Some(j) = parent_edge(&pnodes, el.nodes) else { continue };
            let h = edge_h(pnodes[EDGES[j][0]], pnodes[EDGES[j][1]]);
            let half = |x: &[f64]| (x[0] + x[1]) / 2.0;
            let (mut sre, mut sim) = (vec![0.0; 9], vec![0.0; 9]);
            let (mut fre, mut fim) = (vec![0.0; 3], vec![0.0; 3]);
            fre[j] = half(&v[0]);
            fim[j] = half(&v[1]);
            sre[j * 3 + j] = half(&v[2]) / h;
            sim[j * 3 + j] = half(&v[3]) / h;
            ctx.default_update(&to_real(3, &sre, &sim, &fre, &fim), parent)?;
        }
        ctx.default_finish_boundary_assembly()?;

        // MyDirichletBCs.
        let mut rows: Vec<(usize, Option<f64>)> = Vec::new();
        for dof in 0..2 {
            let comp = format!("{name} {}", dof + 1);
            for &b in &boundary {
                let el = ctx.mesh.element(b);
                if el.etype.0 / 100 == 1 {
                    continue;
                }
                let Some(bc) = ctx.bc_params(b) else { continue };
                if bc.get(&comp).is_none() {
                    continue;
                }
                let data = at_nodes(ctx, Some(bc), &comp, b, el.nodes.len())?;
                let Some(parent) = el.parents[0].or(el.parents[1]) else { continue };
                let pnodes = ctx.mesh.element(parent).nodes.to_vec();
                let Some(j) = parent_edge(&pnodes, el.nodes) else { continue };
                ctx.element_slots(parent, &mut slots);
                let s = slots[j];
                if s == NO_DOF {
                    continue;
                }
                // The mesh edge has its (only) triangle's node order.
                let (n1, n2) = (pnodes[EDGES[j][0]], pnodes[EDGES[j][1]]);
                let h = edge_h(n1, n2);
                let val = (data[0] + data[1]) / 2.0 * h;
                let val = if n2 > n1 { val } else { -val };
                for l in 0..2 {
                    rows.push((2 * s as usize + l, (l == dof).then_some(val)));
                }
            }
        }
        {
            let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no linear system"))?;
            let a = &mut sys.matrix;
            let diag = a.diag.clone().ok_or_else(|| SolverError::fatal(CALLER, "matrix without a diagonal"))?;
            for (r, v) in rows {
                a.values[a.row_ptr[r]..a.row_ptr[r + 1]].fill(0.0);
                if let Some(v) = v {
                    sys.rhs[r] = v;
                }
                a.values[diag[r]] = 1.0;
            }
        }
        ctx.default_solve()?;
        Ok(())
    }
}
