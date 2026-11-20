//! The residual estimators HeatSolve gives RefineMesh
//! (`HeatSolver_Inside_Residual`, `HeatSolver_Edge_Residual`,
//! `HeatSolver_Boundary_Residual`, HeatSolve.F90) for adaptive mesh
//! refinement, in Cartesian coordinates on linear elements.
//!
//! Refused: other coordinate systems, convection, the compressible
//! (perfect gas) density, radiation boundaries, transient runs and p
//! elements.

use elmer_adaptive::{element_hk, Edge, Estimator};
use elmer_core::{ElemId, ElementType, MeshView, ParamMap, Params};
use elmer_solver::{SolverCtx, SolverError};

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("HeatSolve adaptivity: {what}"))
}

/// HeatSolve's estimator over its temperature field.
pub(crate) struct HeatEstimator;

/// `SUM( a(1:n) * b(1:n) )`.
fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// The temperature at the nodes of `e`, `None` where one has no dof.
fn temperature(ctx: &SolverCtx<'_>, nodes: &[u32]) -> Option<Vec<f64>> {
    let v = ctx.primary()?;
    let perm = v.perm?;
    nodes
        .iter()
        .map(|&n| {
            let k = perm[n as usize];
            (k != elmer_core::NO_DOF).then(|| v.values[k as usize])
        })
        .collect()
}

fn material(ctx: &SolverCtx<'_>, e: ElemId) -> Result<ParamMap, SolverError> {
    ctx.body_params(e)
        .and_then(|b| b.material)
        .cloned()
        .ok_or_else(|| SolverError::fatal("HeatSolve", "no material"))
}

/// A keyword at the nodes of element `e`, zero when absent.
fn at_nodes(
    ctx: &SolverCtx<'_>,
    p: &ParamMap,
    key: &str,
    e: ElemId,
) -> Result<(Vec<f64>, bool), SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    Ok(match ctx.real_at_nodes(p, key, e)? {
        Some(v) => (v, true),
        None => (vec![0.0; n], false),
    })
}

fn basis_at(
    ctx: &SolverCtx<'_>,
    etype: ElementType,
    xyz: &[[f64; 3]],
    uvw: [f64; 3],
) -> Result<elmer_core::BasisEval, SolverError> {
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
    Ok(fe
        .quadrature
        .gauss_points(etype, None)
        .iter()
        .map(|q| (q.uvw, q.weight))
        .collect())
}

fn coords(mesh: &dyn MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&k| mesh.coord(k)).collect()
}

/// `Element % TYPE % NodeU/V/W` of a linear triangle or quadrilateral.
fn node_uvw(etype: ElementType, k: usize) -> [f64; 3] {
    match etype.0 {
        303 => [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]][k],
        _ => [
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ][k],
    }
}

/// The parent-local coordinates of the edge point with `edge_basis`.
fn parent_uvw(
    parent_type: ElementType,
    parent_nodes: &[u32],
    edge_nodes: &[u32],
    eb: &[f64],
) -> [f64; 3] {
    let loc: Vec<[f64; 3]> = edge_nodes
        .iter()
        .map(|n| {
            let k = parent_nodes
                .iter()
                .position(|p| p == n)
                .expect("an edge of its parent");
            node_uvw(parent_type, k)
        })
        .collect();
    std::array::from_fn(|d| {
        let f: Vec<f64> = loc.iter().map(|c| c[d]).collect();
        sum(eb, &f)
    })
}

impl HeatEstimator {
    pub(crate) fn check(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
        let sim = ctx.simulation_params();
        if let Some(cs) = sim.string("coordinate system") {
            let cs = elmer_core::normalize_key(cs);
            if !cs.starts_with("cartesian") {
                return Err(unsupported(&format!("Coordinate System = {cs}")));
            }
        }
        if sim
            .string("simulation type")
            .is_some_and(|s| elmer_core::normalize_key(s) == "transient")
        {
            return Err(unsupported("transient runs"));
        }
        let model = &ctx.sim.model;
        if model.equations.iter().any(|e| {
            e.get("convection")
                .is_some_and(|v| *v != elmer_core::Value::Str("none".into()))
        }) {
            return Err(unsupported("convection"));
        }
        if model
            .materials
            .iter()
            .any(|m| m.get("compressibility model").is_some())
        {
            return Err(unsupported("Compressibility Model"));
        }
        if model
            .boundary_conditions
            .iter()
            .any(|b| b.get("radiation").is_some())
        {
            return Err(unsupported("radiation boundaries"));
        }
        for e in 0..ctx.mesh.n_bulk() {
            let t = ctx.mesh.element(ElemId(e as u32)).etype.0;
            if t != 303 && t != 404 {
                return Err(unsupported(&format!("element type {t}")));
            }
        }
        Ok(())
    }
}

impl Estimator for HeatEstimator {
    fn inside(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let el = ctx.mesh.element(e);
        let nodes = el.nodes.to_vec();
        let n = nodes.len();
        let Some(t) = temperature(ctx, &nodes) else {
            return Ok((0.0, [0.0; 2]));
        };
        let dim = ctx.mesh.dim();
        let xyz = coords(ctx.mesh, &nodes);
        let mat = material(ctx, e)?;
        let (cond, _) = at_nodes(ctx, &mat, "heat conductivity", e)?;
        let (dens, _) = at_nodes(ctx, &mat, "density", e)?;
        let bf = ctx.body_params(e).and_then(|b| b.body_force).cloned();
        let (src, vol) = match &bf {
            Some(bf) => match at_nodes(ctx, bf, "volumetric heat source", e)? {
                (v, true) => (v, true),
                _ => (at_nodes(ctx, bf, "heat source", e)?.0, false),
            },
            None => (vec![0.0; n], false),
        };
        let mut fnorm = 0.0;
        let mut rnorm = 0.0;
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &xyz, uvw)?;
            let b = &ev.basis[..n];
            let s = w * ev.det_j;
            let density = sum(&dens, b);
            let mut residual = if vol {
                -sum(&src, b)
            } else {
                -density * sum(&src, b)
            };
            for j in 0..dim {
                let dbj: Vec<f64> = ev.dbasis_dx[..n].iter().map(|d| d[j]).collect();
                residual -= sum(&t, &dbj) * sum(&cond, &dbj);
                // `Conductivity * SUM( T * ddBasisddx(:,j,j) )`: zero for
                // linear elements.
            }
            for _ in 0..dim {
                fnorm += s * (density * sum(&src, b)).powi(2);
            }
            rnorm += s * residual.powi(2);
        }
        let hk = element_hk(ctx, e)?;
        let ind = hk.powi(2) * rnorm;
        Ok((fnorm, [ind, ind]))
    }

    fn edge(&mut self, ctx: &SolverCtx<'_>, edge: &Edge) -> Result<[f64; 2], SolverError> {
        let dim = ctx.mesh.dim();
        let en = edge.nodes.to_vec();
        let exyz = coords(ctx.mesh, &en);
        let mut grad = [[0.0f64; 2]; 3];
        let mut len = 0.0;
        let mut rnorm = 0.0;
        // The unchecked normal of the line (NormalVector, Check = .FALSE.).
        let (dxdu, dydu) = (
            0.0 + -0.5 * exyz[0][0] + 0.5 * exyz[1][0],
            0.0 + -0.5 * exyz[0][1] + 0.5 * exyz[1][1],
        );
        let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
        let normal = [-dydu * det, dxdu * det, 0.0];
        for (uvw, w) in gauss(ctx, ElementType(202))? {
            let ev = basis_at(ctx, ElementType(202), &exyz, uvw)?;
            let eb = ev.basis[..2].to_vec();
            let s = w * ev.det_j;
            for (i, side) in [edge.left, edge.right].into_iter().enumerate() {
                let Some(p) = side else { continue };
                let pe = ElemId(p as u32);
                let pel = ctx.mesh.element(pe);
                let pnodes = pel.nodes.to_vec();
                let Some(t) = temperature(ctx, &pnodes) else {
                    continue;
                };
                let uvw = parent_uvw(pel.etype, &pnodes, &en, &eb);
                let pxyz = coords(ctx.mesh, &pnodes);
                let pev = basis_at(ctx, pel.etype, &pxyz, uvw)?;
                let mat = material(ctx, pe)?;
                let (cond_p, _) = at_nodes(ctx, &mat, "heat conductivity", pe)?;
                let cond_e: Vec<f64> = en
                    .iter()
                    .map(|n| cond_p[pnodes.iter().position(|q| q == n).expect("edge node")])
                    .collect();
                let conductivity = sum(&cond_e, &eb);
                for j in 0..dim {
                    let dbj: Vec<f64> =
                        pev.dbasis_dx[..pnodes.len()].iter().map(|d| d[j]).collect();
                    grad[j][i] = conductivity * sum(&dbj, &t);
                }
            }
            len += s;
            let mut jump = 0.0;
            for k in 0..dim {
                jump += (grad[k][0] - grad[k][1]) * normal[k];
            }
            rnorm += s * jump.powi(2);
        }
        let ind = len * rnorm;
        Ok([ind, ind])
    }

    fn boundary(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let el = ctx.mesh.element(e);
        let en = el.nodes.to_vec();
        let tag = el.tag;
        let [left, right] = el.parents;
        let in_system = |p: ElemId| temperature(ctx, ctx.mesh.element(p).nodes).is_some();
        let parent = match left {
            None => right,
            Some(l) if !in_system(l) => right,
            Some(l) => Some(l),
        };
        let Some(pe) = parent else {
            return Ok((0.0, [0.0; 2]));
        };
        let pel = ctx.mesh.element(pe);
        let pnodes = pel.nodes.to_vec();
        let Some(t) = temperature(ctx, &pnodes) else {
            return Ok((0.0, [0.0; 2]));
        };
        let Some(bc) = ctx
            .sim
            .model
            .bc_of_tag(tag)
            .map(|i| ctx.sim.model.boundary_conditions[i - 1].clone())
        else {
            return Ok((0.0, [0.0; 2]));
        };
        let dim = ctx.mesh.dim();
        let exyz = coords(ctx.mesh, &en);
        let pxyz = coords(ctx.mesh, &pnodes);
        let dirichlet = bc.get("temperature").is_some();
        let (flux, _) = at_nodes(ctx, &bc, "heat flux", e)?;
        let (tc, _) = at_nodes(ctx, &bc, "heat transfer coefficient", e)?;
        let (ext, _) = at_nodes(ctx, &bc, "external temperature", e)?;
        let mat = material(ctx, pe)?;
        let (cond_p, _) = at_nodes(ctx, &mat, "heat conductivity", pe)?;
        let cond_e: Vec<f64> = en
            .iter()
            .map(|n| cond_p[pnodes.iter().position(|q| q == n).expect("edge node")])
            .collect();
        let mut gnorm = 0.0;
        let mut len = 0.0;
        let mut rnorm = 0.0;
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &exyz, uvw)?;
            let eb = ev.basis[..en.len()].to_vec();
            let normal =
                elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
            let s = w * ev.det_j;
            let puvw = parent_uvw(pel.etype, &pnodes, &en, &eb);
            let pev = basis_at(ctx, pel.etype, &pxyz, puvw)?;
            let pb = &pev.basis[..pnodes.len()];
            let conductivity = sum(&cond_e, &eb);
            let mut residual = -sum(&flux, &eb);
            residual += sum(&tc, &eb) * (sum(&t, pb) - sum(&ext, &eb));
            for k in 0..dim {
                let dbk: Vec<f64> = pev.dbasis_dx[..pnodes.len()].iter().map(|d| d[k]).collect();
                let g = sum(&dbk, &t);
                residual += conductivity * g * normal[k];
                gnorm += s * (conductivity * g * normal[k]).powi(2);
            }
            len += s;
            if !dirichlet {
                rnorm += s * residual.powi(2);
            }
        }
        let ind = len * rnorm;
        Ok((gnorm, [ind, ind]))
    }
}
