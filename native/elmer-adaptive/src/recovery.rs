//! `FluxRecovery` (Adaptive.F90): RefineMesh's error indicator of `Flux
//! Recovery`. Each triangle's flux is equilibrated locally in the H(div)
//! space of the `RTFlux` variable (the solver whose `Variable = RTFlux`,
//! e.g. `VectorHelmholtz_Dummy` with `Element = "n:0 e:2"`). The fluxes are
//! averaged over the elements sharing a face and replaced by the nodal
//! reactions on Dirichlet boundaries. The indicator of an element is the
//! norm of `EA grad u_h - N_h` there, relative to the norm of `N_h`.
//!
//! Ported: the BDM(1) space (`Second Kind Basis = True`) on linear
//! triangles in Cartesian 2D, serial. Refused: RT (the first-kind space
//! with its first-order equilibration), other elements, and partitions.

use elmer_core::{ElemId, Params, NO_DOF};
use elmer_fe_hp::face::{d_triangle_nodal_p_basis, face_element_info_triangle, triangle_face_orientation, FaceInfo};
use elmer_solver::{SolverCtx, SolverError};

/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// `GetEdgeMap(3)`, 0-based.
const EDGE_MAP: [[usize; 2]; 3] = [[0, 1], [1, 2], [2, 0]];

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("FluxRecovery: {what}"))
}

fn fatal(what: impl Into<String>) -> SolverError {
    SolverError::fatal("FluxRecovery", what)
}

/// The `RTFlux` solver: its slot and keywords.
struct RtSolver {
    slot: usize,
    matpar: String,
    force: String,
}

/// A triangle of the heat solver: nodes, coordinates and solution.
struct Tri {
    g: [u32; 3],
    x: [f64; 3],
    y: [f64; 3],
    z: [f64; 3],
}

impl Tri {
    fn of(ctx: &SolverCtx<'_>, e: ElemId) -> Tri {
        let el = ctx.mesh.element(e);
        let mut t = Tri { g: [0; 3], x: [0.0; 3], y: [0.0; 3], z: [0.0; 3] };
        for i in 0..3 {
            let c = ctx.mesh.coord(el.nodes[i]);
            t.g[i] = el.nodes[i];
            (t.x[i], t.y[i], t.z[i]) = (c[0], c[1], c[2]);
        }
        t
    }

    fn info(&self, u: f64, v: f64) -> FaceInfo {
        face_element_info_triangle(&self.x, &self.y, &self.z, &self.g, u, v, true, true)
    }
}

/// `detF` before its absolute value: `LeftHanded = detF < 0`.
fn signed_det(fi: &FaceInfo) -> f64 {
    fi.f[0][0] * fi.f[1][1] - fi.f[0][1] * fi.f[1][0]
}

/// `dBasisdx` of FaceElementInfo for a triangle: the nodal p basis
/// gradients through `LG`, the transposed inverse of the map gradient.
fn dbasis(fi: &FaceInfo) -> [[f64; 2]; 3] {
    let lf = &fi.f;
    let det = signed_det(fi);
    let lg0 = [[1.0 / det * lf[1][1], -1.0 / det * lf[0][1]], [-1.0 / det * lf[1][0], 1.0 / det * lf[0][0]]];
    let lg = [[lg0[0][0], lg0[1][0]], [lg0[0][1], lg0[1][1]]];
    let mut out = [[0.0; 2]; 3];
    for (i, o) in out.iter_mut().enumerate() {
        let dl = d_triangle_nodal_p_basis(i + 1);
        for j in 0..2 {
            for k in 0..2 {
                o[j] += dl[k] * lg[j][k];
            }
        }
    }
    out
}

/// The sign of face `p` (`FaceElementOrientation` and `LeftHanded`).
fn face_sign(left: bool, reverse: bool) -> f64 {
    match (left, reverse) {
        (true, true) => 1.0,
        (true, false) => -1.0,
        (false, true) => -1.0,
        (false, false) => 1.0,
    }
}

fn weights() -> (f64, f64) {
    (0.5 * (1.0 + 1.0 / 3f64.sqrt()), 0.5 * (1.0 - 1.0 / 3f64.sqrt()))
}

fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// What `EstimateError` needs of an element.
struct Local {
    tri: Tri,
    uk: [f64; 3],
    diff: Vec<f64>,
    load: Vec<f64>,
}

fn rule() -> Result<elmer_fe::quadrature::GaussRule, SolverError> {
    elmer_pdofs::quadrature::gauss_points_triangle_p(6).map_err(|e| fatal(e.to_string()))
}

/// `EstimateError(..., LinFun=)`: the local BDM(1) flux with the Lagrange
/// multiplier of the element balance.
fn lin_fun(l: &Local) -> Result<[f64; 6], SolverError> {
    const F: usize = 6;
    const D: usize = 7;
    let ip = rule()?;
    let rev = triangle_face_orientation(&l.tri.g);
    let mut mass = [[0.0f64; D]; D];
    let mut rhs = [0.0f64; D];
    let (w1, w2) = weights();
    let mut left = false;
    for t in 0..ip.len() {
        let fi = l.tri.info(ip.u[t], ip.v[t]);
        left = signed_det(&fi) < 0.0;
        let w = ip.s[t] * fi.det_f;
        let uh = sum(&l.uk, &fi.basis[..3]);
        let ea = sum(&l.diff, &fi.basis[..3]);
        for p in 0..F {
            for q in 0..F {
                let mut s = 0.0;
                for k in 0..2 {
                    s += fi.fbasis[q][k] * fi.fbasis[p][k];
                }
                mass[p][q] += s * w / ea;
            }
            rhs[p] -= w * uh * fi.div[p];
        }
        let f = sum(&l.load, &fi.basis[..3]);
        rhs[F] += w * f;
    }
    for p in 0..3 {
        let s = face_sign(left, rev[p]);
        let [i, j] = EDGE_MAP[p];
        let (ni, nj) = (l.tri.g[i], l.tri.g[j]);
        if nj < ni {
            rhs[2 * p] += s * (w2 * l.uk[i] + w1 * l.uk[j]);
            rhs[2 * p + 1] += s * (w1 * l.uk[i] + w2 * l.uk[j]);
        } else {
            rhs[2 * p] += s * (w1 * l.uk[i] + w2 * l.uk[j]);
            rhs[2 * p + 1] += s * (w2 * l.uk[i] + w1 * l.uk[j]);
        }
        mass[F][2 * p] = -s;
        mass[F][2 * p + 1] = -s;
        mass[2 * p][F] = -s;
        mass[2 * p + 1][F] = -s;
    }
    // `InvertMatrix(Mass(1:DOFs,1:DOFs), DOFs)`, column-major.
    let mut a = vec![0.0f64; D * D];
    for i in 0..D {
        for j in 0..D {
            a[i + j * D] = mass[i][j];
        }
    }
    elmer_fe::invert_matrix(&mut a, D).ok_or_else(|| fatal("InvertMatrix: inversion needs successful LU-decomposition"))?;
    let mut c = [0.0f64; F];
    for (i, ci) in c.iter_mut().enumerate() {
        let mut s = 0.0;
        for j in 0..D {
            s += a[i + j * D] * rhs[j];
        }
        *ci = s;
    }
    Ok(c)
}

/// `EstimateError(..., PostLinFun=, APostEst_K=, SolNormEst=)` with the BDM
/// flux `c`: the element's estimate, and its `SolNormEst` share.
fn post(l: &Local, c: &[f64; 6]) -> Result<(f64, f64), SolverError> {
    let ip = rule()?;
    let (mut est, mut sol) = (0.0, 0.0);
    for t in 0..ip.len() {
        let fi = l.tri.info(ip.u[t], ip.v[t]);
        let w = ip.s[t] * fi.det_f;
        let db = dbasis(&fi);
        let mut grad = [0.0f64; 2];
        for (i, g) in grad.iter_mut().enumerate() {
            let d: Vec<f64> = db.iter().map(|r| r[i]).collect();
            *g = sum(&l.uk, &d);
        }
        let ea = sum(&l.diff, &fi.basis[..3]);
        let mut nh = [0.0f64; 2];
        for (i, n) in nh.iter_mut().enumerate() {
            let fb: Vec<f64> = fi.fbasis[..6].iter().map(|r| r[i]).collect();
            *n = sum(c, &fb);
        }
        let mut s2 = 0.0;
        for k in 0..2 {
            s2 += nh[k].powi(2);
        }
        sol += s2 * w;
        let mut e2 = 0.0;
        for k in 0..2 {
            e2 += (ea * grad[k] - nh[k]).powi(2);
        }
        est += e2 * w;
    }
    Ok((est, sol))
}

fn rt_solver(ctx: &SolverCtx<'_>) -> Result<RtSolver, SolverError> {
    let id = ctx.sim.vars.id("rtflux").ok_or_else(|| fatal("RTFlux variable cannot be found"))?;
    let slot = ctx
        .sim
        .slots
        .iter()
        .position(|s| s.variable == Some(id))
        .ok_or_else(|| fatal("RTFlux variable is not associated with any solver"))?;
    let p = &ctx.sim.slots[slot].params;
    if p.logical("second kind basis") != Some(true) {
        return Err(unsupported("the Raviart-Thomas space (Second Kind Basis = False)"));
    }
    if ctx.sim.slots[slot].solver_dofs().is_none() {
        return Err(fatal("the RTFlux solver has no face dofs"));
    }
    Ok(RtSolver {
        slot,
        matpar: p.string("material parameter name").unwrap_or("").to_string(),
        force: p.string("body force name").unwrap_or("").to_string(),
    })
}

fn local(ctx: &SolverCtx<'_>, rt: &RtSolver, e: ElemId) -> Result<Option<Local>, SolverError> {
    let el = ctx.mesh.element(e);
    let prim = ctx.primary().ok_or_else(|| fatal("solver has no variable"))?;
    let perm = prim.perm.ok_or_else(|| fatal("solver variable has no permutation"))?;
    let mut uk = [0.0; 3];
    for i in 0..3 {
        let k = perm[el.nodes[i] as usize];
        if k == NO_DOF {
            return Ok(None);
        }
        uk[i] = prim.values[k as usize];
    }
    let bp = ctx.body_params(e).ok_or_else(|| fatal("element without a body"))?;
    let mat = bp.material.ok_or_else(|| fatal("element without a material"))?;
    let diff = ctx.real_at_nodes(mat, &rt.matpar, e)?.unwrap_or_else(|| vec![0.0; 3]);
    let load = match bp.body_force {
        Some(bf) => ctx.real_at_nodes(bf, &rt.force, e)?.unwrap_or_else(|| vec![0.0; 3]),
        None => vec![0.0; 3],
    };
    Ok(Some(Local { tri: Tri::of(ctx, e), uk, diff, load }))
}

fn rt_rows(ctx: &SolverCtx<'_>, rt: &RtSolver, e: ElemId) -> Vec<u32> {
    let mut out = Vec::new();
    ctx.sim.slots[rt.slot].solver_dofs().expect("checked").element_slots(e, &mut out);
    out
}

/// `FluxRecovery(Model, Solver, Mesh, ErrorIndicator, MaxError)`: fills
/// `ErrorIndicator` (one per active element, the bulk in order) and
/// returns `MaxError`.
pub(crate) fn flux_recovery(ctx: &mut SolverCtx<'_>, ei: &mut [f64]) -> Result<f64, SolverError> {
    if ctx.comm.size() > 1 {
        return Err(unsupported("partitioned runs"));
    }
    let rt = rt_solver(ctx)?;
    let active: Vec<ElemId> = ctx.active_elements().collect();
    if active.iter().any(|&e| ctx.mesh.element(e).etype.0 != 303) {
        return Err(unsupported("elements other than linear triangles"));
    }
    let n_rt = ctx.sim.vars.get("rtflux").expect("checked").values.len();
    let mut rt_vals = vec![0.0f64; n_rt];
    let mut visits = vec![0u32; n_rt];
    for &e in &active {
        let Some(l) = local(ctx, &rt, e)? else { continue };
        let lin = lin_fun(&l)?;
        for (i, &j) in rt_rows(ctx, &rt, e).iter().enumerate().take(6) {
            if j == NO_DOF {
                continue;
            }
            rt_vals[j as usize] += lin[i];
            visits[j as usize] += 1;
        }
    }
    for (v, &c) in rt_vals.iter_mut().zip(&visits) {
        if c > 1 {
            *v /= f64::from(c);
        }
    }

    // The reactions on the Dirichlet boundary.
    let prim = ctx.primary().ok_or_else(|| fatal("solver has no variable"))?;
    let name = prim.name.to_string();
    let perm: Vec<u32> = prim.perm.map(|p| p.to_vec()).ok_or_else(|| fatal("no permutation"))?;
    let n_vals = prim.values.len();
    let constrained = ctx.constrained_rows();
    let loads = ctx.sim.vars.get(&format!("{name} loads")).map(|v| v.values.clone());
    if let (Some(lv), false) = (&loads, constrained.is_empty()) {
        let mut is_c = vec![false; n_vals];
        for &r in &constrained {
            if (r as usize) < is_c.len() {
                is_c[r as usize] = true;
            }
        }
        // `NodalLoads % Values(Solver % Variable % Perm(node))`.
        let load_at = |node: u32| -> f64 { lv[perm[node as usize] as usize] };
        let nb = ctx.mesh.n_bulk();
        let bnd: Vec<ElemId> = (nb..nb + ctx.mesh.n_boundary()).map(|e| ElemId(e as u32)).collect();
        let full = |e: ElemId| -> bool {
            let el = ctx.mesh.element(e);
            el.etype.0 / 100 == 2
                && el.nodes.iter().all(|&n| perm[n as usize] != NO_DOF)
                && el.nodes.iter().all(|&n| is_c[perm[n as usize] as usize])
        };
        let h = |e: ElemId| -> f64 {
            let el = ctx.mesh.element(e);
            let (a, b) = (ctx.mesh.coord(el.nodes[0]), ctx.mesh.coord(el.nodes[1]));
            ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt()
        };
        let mut rw = vec![0.0f64; n_vals.max(ctx.mesh.n_nodes())];
        for &e in &bnd {
            if full(e) {
                let hk = h(e);
                for &n in ctx.mesh.element(e).nodes {
                    rw[n as usize] += hk;
                }
            }
        }
        let (w1, w2) = weights();
        let detw = w1.powi(2) - w2.powi(2);
        for &e in &bnd {
            if !full(e) {
                continue;
            }
            let el = ctx.mesh.element(e);
            let parent = el.parents[0].or(el.parents[1]).ok_or_else(|| fatal("A parent element is not defined"))?;
            let pel = ctx.mesh.element(parent);
            if pel.etype.0 != 303 {
                return Err(unsupported("a parent other than a linear triangle"));
            }
            let face = (0..3)
                .find(|&f| {
                    let [a, b] = EDGE_MAP[f];
                    let pn = [pel.nodes[a], pel.nodes[b]];
                    el.nodes.iter().all(|n| pn.contains(n))
                })
                .ok_or_else(|| fatal("Cannot determine an element face"))?;
            let tri = Tri::of(ctx, parent);
            let rev = triangle_face_orientation(&tri.g);
            // IsLeftHanded(Parent): FaceElementInfo at (0, sqrt(3)/3).
            let left = signed_det(&tri.info(0.0, 3f64.sqrt() / 3.0)) < 0.0;
            let s = face_sign(left, rev[face]);
            let (ir, jr) = (el.nodes[0], el.nodes[1]);
            let i = if rev[face] { pel.nodes[EDGE_MAP[face][1]] } else { pel.nodes[EDGE_MAP[face][0]] };
            let hk = h(e);
            let (r1, r2) = if ir == i {
                (s * load_at(ir) * hk / rw[ir as usize], s * load_at(jr) * hk / rw[jr as usize])
            } else {
                (s * load_at(jr) * hk / rw[jr as usize], s * load_at(ir) * hk / rw[ir as usize])
            };
            let rows = rt_rows(ctx, &rt, parent);
            let (a, b) = (rows[2 * face], rows[2 * face + 1]);
            if a == NO_DOF || b == NO_DOF {
                continue;
            }
            rt_vals[a as usize] = 1.0 / detw * (w1 * r1 - w2 * r2);
            rt_vals[b as usize] = 1.0 / detw * (-w2 * r1 + w1 * r2);
        }
    }

    // The element-wise estimates with the equilibrated fluxes.
    let mut sol_norm = 0.0;
    for (k, &e) in active.iter().enumerate() {
        let Some(l) = local(ctx, &rt, e)? else { continue };
        let rows = rt_rows(ctx, &rt, e);
        let mut c = [0.0f64; 6];
        for (i, ci) in c.iter_mut().enumerate() {
            if rows[i] != NO_DOF {
                *ci = rt_vals[rows[i] as usize];
            }
        }
        let (est, sol) = post(&l, &c)?;
        sol_norm += sol;
        if k < ei.len() {
            ei[k] = est.sqrt();
        }
    }
    if sol_norm > AEPS {
        let f = 1.0 / sol_norm.sqrt();
        for x in ei.iter_mut() {
            *x *= f;
        }
    }
    {
        let v = ctx.sim.vars.get_mut("rtflux").expect("checked");
        v.values.copy_from_slice(&rt_vals);
    }
    Ok(ei.iter().copied().fold(f64::MIN, f64::max))
}
