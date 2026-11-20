//! The residual estimators FlowSolve gives RefineMesh
//! (`FlowSolver_Inside_Residual`, `FlowSolver_Edge_Residual`,
//! `FlowSolver_Boundary_Residual`, FlowSolve.F90) over the flow solution
//! (velocity and pressure), in Cartesian coordinates on linear triangles
//! (SplitMesh refines only those).
//!
//! Refused: other coordinate systems and element types, transient runs,
//! the compressible (perfect gas) density, Boussinesq forces, surface
//! tension and `Circle Boundary`.

use elmer_core::{ElemId, ElementType, ParamMap, Params};
use elmer_solver::{SolverCtx, SolverError};

use crate::{element_hk, Edge, Estimator};

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("FlowSolve adaptivity: {what}"))
}

/// `SUM( a(1:n) * b(1:n) )`.
fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// FlowSolve's estimator over `Flow Solution`.
pub struct FlowEstimator;

/// Velocity (3 components, zero beyond the dimension) and pressure at
/// `nodes`; `None` when a node has no dof (`ANY( Perm <= 0 )`).
fn flow_at(ctx: &SolverCtx<'_>, nodes: &[u32], dim: usize) -> Option<(Vec<[f64; 3]>, Vec<f64>)> {
    let v = ctx.primary()?;
    let perm = v.perm?;
    let dofs = dim + 1;
    let mut vel = Vec::with_capacity(nodes.len());
    let mut p = Vec::with_capacity(nodes.len());
    for &n in nodes {
        let k = perm[n as usize];
        if k == elmer_core::NO_DOF {
            return None;
        }
        let base = dofs * k as usize;
        let mut u = [0.0; 3];
        for (c, x) in u.iter_mut().enumerate().take(dim) {
            *x = v.values[base + c];
        }
        vel.push(u);
        p.push(v.values[base + dim]);
    }
    Some((vel, p))
}

fn material(ctx: &SolverCtx<'_>, e: ElemId) -> Result<ParamMap, SolverError> {
    ctx.body_params(e)
        .and_then(|b| b.material)
        .cloned()
        .ok_or_else(|| SolverError::fatal("RefineMesh", "the element has no material"))
}

/// `ListGetReal( list, key, n, NodeIndexes )` at the nodes of `e`, zeros
/// when absent.
fn at_nodes(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str, e: ElemId) -> Result<Vec<f64>, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    Ok(ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]))
}

fn basis_at(ctx: &SolverCtx<'_>, etype: ElementType, xyz: &[[f64; 3]], uvw: [f64; 3]) -> Result<elmer_core::BasisEval, SolverError> {
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let mut ev = elmer_core::BasisEval::default();
    fe.basis
        .eval(etype, xyz, uvw, &mut ev)
        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
    Ok(ev)
}

fn gauss(ctx: &SolverCtx<'_>, etype: ElementType) -> Result<Vec<([f64; 3], f64)>, SolverError> {
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    Ok(fe.quadrature.gauss_points(etype, None).iter().map(|q| (q.uvw, q.weight)).collect())
}

/// The parent-local point of the edge point with basis `eb` (linear
/// triangle parents).
fn parent_uvw(pnodes: &[u32], edge_nodes: &[u32], eb: &[f64]) -> [f64; 3] {
    const UVW: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let loc: Vec<[f64; 3]> = edge_nodes
        .iter()
        .map(|n| UVW[pnodes.iter().position(|q| q == n).expect("an edge of its parent")])
        .collect();
    std::array::from_fn(|d| {
        let f: Vec<f64> = loc.iter().map(|c| c[d]).collect();
        sum(eb, &f)
    })
}

fn coords(ctx: &SolverCtx<'_>, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&k| ctx.mesh.coord(k)).collect()
}

/// The values of a parent-node keyword at the edge nodes.
fn on_edge(values: &[f64], pnodes: &[u32], en: &[u32]) -> Vec<f64> {
    en.iter().map(|n| values[pnodes.iter().position(|q| q == n).expect("edge node")]).collect()
}

/// `Grad = MATMUL( Velocity(:,1:Pn), dBasisdx(1:Pn,:) )`.
fn grad(vel: &[[f64; 3]], db: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let mut g = [[0.0; 3]; 3];
    for (j, row) in g.iter_mut().enumerate() {
        for (k, x) in row.iter_mut().enumerate() {
            let mut s = 0.0;
            for (v, d) in vel.iter().zip(db) {
                s += v[j] * d[k];
            }
            *x = s;
        }
    }
    g
}

impl FlowEstimator {
    /// Refuses what this port does not cover.
    pub fn check(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
        let sim = ctx.simulation_params();
        if let Some(cs) = sim.string("coordinate system") {
            let cs = elmer_core::normalize_key(cs);
            if !cs.starts_with("cartesian") {
                return Err(unsupported(&format!("Coordinate System = {cs}")));
            }
        }
        if sim.string("simulation type").is_some_and(|s| elmer_core::normalize_key(s) == "transient") {
            return Err(unsupported("transient runs"));
        }
        let model = &ctx.sim.model;
        if model.materials.iter().any(|m| m.get("compressibility model").is_some()) {
            return Err(unsupported("Compressibility Model"));
        }
        if model.body_forces.iter().any(|b| b.logical("boussinesq") == Some(true)) {
            return Err(unsupported("Boussinesq"));
        }
        if model.boundary_conditions.iter().any(|b| {
            b.get("surface tension coefficient").is_some()
                || b.get("surface tension expansion coefficient").is_some()
                || b.logical("circle boundary") == Some(true)
        }) {
            return Err(unsupported("surface tension or Circle Boundary"));
        }
        for e in 0..ctx.mesh.n_bulk() {
            let t = ctx.mesh.element(ElemId(e as u32)).etype.0;
            if t != 303 {
                return Err(unsupported(&format!("element type {t}")));
            }
        }
        Ok(())
    }
}

impl Estimator for FlowEstimator {
    /// `FlowSolver_Inside_Residual`, steady: the time term vanishes and
    /// `ddBasisddx` is zero on linear triangles.
    fn inside(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let el = ctx.mesh.element(e);
        let nodes = el.nodes.to_vec();
        let n = nodes.len();
        let dim = ctx.mesh.dim();
        let Some((vel, pres)) = flow_at(ctx, &nodes, dim) else { return Ok((0.0, [0.0; 2])) };
        let xyz = coords(ctx, &nodes);
        let mat = material(ctx, e)?;
        let dens = at_nodes(ctx, &mat, "density", e)?;
        let visc = at_nodes(ctx, &mat, "viscosity", e)?;
        let convect = ctx
            .body_params(e)
            .and_then(|b| b.equation)
            .and_then(|q| q.logical("ns convect"))
            .unwrap_or(true);
        let mut force = vec![vec![0.0; n]; 3];
        if let Some(bf) = ctx.body_params(e).and_then(|b| b.body_force).cloned() {
            for (i, f) in force.iter_mut().enumerate() {
                *f = at_nodes(ctx, &bf, &format!("flow bodyforce {}", i + 1), e)?;
            }
        }
        let hk = element_hk(ctx, e)?;
        let (mut fnorm, mut rnorm) = (0.0, 0.0);
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &xyz, uvw)?;
            let b = &ev.basis[..n];
            let db = &ev.dbasis_dx[..n];
            let s = w * ev.det_j;
            let density = sum(&dens, b);
            let viscosity = sum(&visc, b);
            let col = |j: usize| -> Vec<f64> { db.iter().map(|d| d[j]).collect() };
            let comp = |i: usize| -> Vec<f64> { vel.iter().map(|v| v[i]).collect() };
            let mut residual = [0.0f64; 4];
            for i in 0..dim {
                residual[i] = -density * sum(&force[i], b);
                residual[i] += sum(&pres, &col(i));
                for j in 0..dim {
                    // `- Viscosity * SUM( Velocity(i) * ddBasisddx(:,j,j) )` and
                    // `- Viscosity * SUM( Velocity(j) * ddBasisddx(:,i,j) )`: zero.
                    residual[i] -= sum(&visc, &col(j)) * sum(&comp(i), &col(j));
                    residual[i] -= sum(&visc, &col(j)) * sum(&comp(j), &col(i));
                }
                if convect {
                    // `Density * SUM((Velocity-PrevVelo)*Basis) / dt`: zero.
                    for j in 0..dim {
                        residual[i] += density * sum(&comp(j), b) * sum(&comp(i), &col(j));
                    }
                }
            }
            for j in 0..dim {
                residual[dim] += density * sum(&comp(j), &col(j));
            }
            for f in force.iter().take(dim) {
                fnorm += s * (density * sum(f, b)).powi(2);
            }
            let _ = viscosity;
            let mut r2 = 0.0;
            for r in residual.iter().take(dim) {
                r2 += r * r;
            }
            rnorm += s * (hk * hk * r2 + residual[dim].powi(2));
        }
        Ok((fnorm, [rnorm, rnorm]))
    }

    /// `FlowSolver_Edge_Residual`: the jump of the stress across the edge.
    fn edge(&mut self, ctx: &SolverCtx<'_>, edge: &Edge) -> Result<[f64; 2], SolverError> {
        let dim = ctx.mesh.dim();
        let en = edge.nodes.to_vec();
        let exyz = coords(ctx, &en);
        // `NormalVector( Edge, EdgeNodes, u, v, .FALSE. )` of a line.
        let (dxdu, dydu) = (0.0 + -0.5 * exyz[0][0] + 0.5 * exyz[1][0], 0.0 + -0.5 * exyz[0][1] + 0.5 * exyz[1][1]);
        let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
        let normal = [-dydu * det, dxdu * det, 0.0];
        let (mut len, mut rnorm) = (0.0, 0.0);
        for (uvw, w) in gauss(ctx, ElementType(202))? {
            let ev = basis_at(ctx, ElementType(202), &exyz, uvw)?;
            let eb = ev.basis[..2].to_vec();
            let s = w * ev.det_j;
            let mut stress = [[[0.0f64; 3]; 3]; 2];
            for (i, side) in [edge.left, edge.right].into_iter().enumerate() {
                let Some(p) = side else { continue };
                let pe = ElemId(p as u32);
                let pel = ctx.mesh.element(pe);
                let pnodes = pel.nodes.to_vec();
                let Some((vel, pres)) = flow_at(ctx, &pnodes, dim) else { continue };
                let pxyz = coords(ctx, &pnodes);
                let pev = basis_at(ctx, pel.etype, &pxyz, parent_uvw(&pnodes, &en, &eb))?;
                let pn = pnodes.len();
                let mat = material(ctx, pe)?;
                let visc_e = on_edge(&at_nodes(ctx, &mat, "viscosity", pe)?, &pnodes, &en);
                let viscosity = sum(&visc_e, &eb);
                let g = grad(&vel, &pev.dbasis_dx[..pn]);
                let st = &mut stress[i];
                for j in 0..3 {
                    for k in 0..3 {
                        st[j][k] = viscosity * (g[j][k] + g[k][j]);
                    }
                }
                let p_ip = sum(&pres, &pev.basis[..pn]);
                for j in 0..dim {
                    st[j][j] -= p_ip;
                    for k in 0..dim {
                        st[j][j] -= (2.0 / 3.0) * viscosity * g[k][k];
                    }
                }
            }
            len += s;
            let mut jump = [0.0f64; 3];
            for (j, x) in jump.iter_mut().enumerate() {
                for k in 0..3 {
                    *x += (stress[0][j][k] - stress[1][j][k]) * normal[k];
                }
            }
            let mut j2 = 0.0;
            for x in jump.iter().take(dim) {
                j2 += x * x;
            }
            rnorm += s * j2;
        }
        let ind = len * rnorm;
        Ok([ind, ind])
    }

    /// `FlowSolver_Boundary_Residual`: the boundary traction against the
    /// given forces, in the directions the BC does not fix.
    fn boundary(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let dim = ctx.mesh.dim();
        let el = ctx.mesh.element(e);
        let en = el.nodes.to_vec();
        let tag = el.tag;
        let [left, right] = el.parents;
        let in_system = |p: ElemId| flow_at(ctx, ctx.mesh.element(p).nodes, dim).is_some();
        let parent = match left {
            None => right,
            Some(l) if !in_system(l) => right,
            Some(l) => Some(l),
        };
        let Some(pe) = parent else { return Ok((0.0, [0.0; 2])) };
        let pel = ctx.mesh.element(pe);
        let pnodes = pel.nodes.to_vec();
        let Some((vel, pres)) = flow_at(ctx, &pnodes, dim) else { return Ok((0.0, [0.0; 2])) };
        let Some(bc) = ctx.sim.model.bc_of_tag(tag).map(|i| ctx.sim.model.boundary_conditions[i - 1].clone()) else {
            return Ok((0.0, [0.0; 2]));
        };
        let pn = pnodes.len();
        let exyz = coords(ctx, &en);
        let pxyz = coords(ctx, &pnodes);
        let mat = material(ctx, pe)?;
        let visc_e = on_edge(&at_nodes(ctx, &mat, "viscosity", pe)?, &pnodes, &en);
        let mut force = vec![vec![0.0; en.len()]; 3];
        for (i, f) in force.iter_mut().enumerate() {
            *f = at_nodes(ctx, &bc, &format!("pressure {}", i + 1), e)?;
        }
        let ext = at_nodes(ctx, &bc, "external pressure", e)?;
        let mut slip = vec![vec![0.0; en.len()]; 3];
        for (i, f) in slip.iter_mut().enumerate() {
            *f = at_nodes(ctx, &bc, &format!("slip coefficient {}", i + 1), e)?;
        }
        let dir: [f64; 3] = std::array::from_fn(|i| if bc.get(&format!("velocity {}", i + 1)).is_some() { 0.0 } else { 1.0 });
        let (mut gnorm, mut len, mut rnorm) = (0.0, 0.0, 0.0);
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &exyz, uvw)?;
            let eb = ev.basis[..en.len()].to_vec();
            let s = w * ev.det_j;
            let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
            let pev = basis_at(ctx, pel.etype, &pxyz, parent_uvw(&pnodes, &en, &eb))?;
            let pb = &pev.basis[..pn];
            let viscosity = sum(&visc_e, &eb);
            let ext_ip = sum(&ext, &eb);
            let mut residual: [f64; 3] = std::array::from_fn(|i| sum(&force[i], &eb) - ext_ip * normal[i]);
            for i in 0..dim {
                let sl = sum(&slip[i], &eb);
                let vi: Vec<f64> = vel.iter().map(|v| v[i]).collect();
                residual[i] -= sl * sum(&vi, pb);
            }
            // `SUM( dEdgeBasisdx(:,i) * Tension )`: no surface tension.
            let g = grad(&vel, &pev.dbasis_dx[..pn]);
            let p_ip = sum(&pres, pb);
            let mut st = [[0.0f64; 3]; 3];
            for j in 0..3 {
                for k in 0..3 {
                    st[j][k] = viscosity * (g[j][k] + g[k][j]);
                }
                st[j][j] -= p_ip;
            }
            let mut fs = [0.0f64; 3];
            for (j, x) in fs.iter_mut().enumerate() {
                for k in 0..3 {
                    *x += st[j][k] * normal[k];
                }
            }
            for i in 0..3 {
                residual[i] -= fs[i] * dir[i];
            }
            len += s;
            gnorm += s * (fs[0] * fs[0] + fs[1] * fs[1] + fs[2] * fs[2]);
            let mut r2 = 0.0;
            for r in residual.iter().take(dim) {
                r2 += r * r;
            }
            rnorm += s * r2;
        }
        if dim == 3 {
            len = len.sqrt();
        }
        let ind = len * rnorm;
        Ok((gnorm, [ind, ind]))
    }
}
