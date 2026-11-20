//! The residual estimators of the scalar diffusion solvers that give
//! RefineMesh their own `_Inside_Residual`, `_Edge_Residual` and
//! `_Boundary_Residual`: StatElecSolve (permittivity, charge density) and
//! StatCurrentSolve (electric conductivity, current source). They are the
//! same code with other keywords, in Cartesian coordinates on linear
//! triangles and quadrilaterals.

use elmer_core::{ElemId, ElementType, ParamMap, Params};
use elmer_solver::{SolverCtx, SolverError};

use crate::{element_hk, Edge, Estimator};

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("diffusion adaptivity: {what}"))
}

fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// The keywords of one solver.
pub struct DiffusionEstimator {
    /// The material coefficient, the first one present (`Relative
    /// Permittivity`, else `Permittivity`).
    pub coefficient: &'static [&'static str],
    /// The body force source, the first one present.
    pub source: &'static [&'static str],
}

impl DiffusionEstimator {
    /// StatElecSolve.F90's estimators.
    pub const STAT_ELEC: DiffusionEstimator = DiffusionEstimator {
        coefficient: &["relative permittivity", "permittivity"],
        source: &["charge density", "source"],
    };

    /// StatCurrentSolve.F90's estimators.
    pub const STAT_CURRENT: DiffusionEstimator = DiffusionEstimator {
        coefficient: &["electric conductivity"],
        source: &["current source"],
    };

    /// What the port refuses, before any residual.
    pub fn check(&self, ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
        let sim = ctx.simulation_params();
        if let Some(cs) = sim.string("coordinate system") {
            if !elmer_core::normalize_key(cs).starts_with("cartesian") {
                return Err(unsupported("a non-Cartesian coordinate system"));
            }
        }
        if sim.string("simulation type").is_some_and(|s| elmer_core::normalize_key(s) == "transient") {
            return Err(unsupported("transient runs"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("3D meshes"));
        }
        for e in 0..ctx.mesh.n_bulk() {
            let t = ctx.mesh.element(ElemId(e as u32)).etype.0;
            if t != 303 && t != 404 {
                return Err(unsupported(&format!("element type {t}")));
            }
        }
        for m in &ctx.sim.model.materials {
            for k in self.coefficient {
                if matches!(m.get(k), Some(elmer_core::Value::Real { rows, cols, .. }) if *rows * *cols > 1) {
                    return Err(unsupported(&format!("a tensor {k}")));
                }
            }
        }
        Ok(())
    }
}

fn potential(ctx: &SolverCtx<'_>, nodes: &[u32]) -> Option<Vec<f64>> {
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

/// `ListGetReal` of the first of `keys` present in `p`, at the nodes of `e`
/// (zeros when none is).
fn first_at_nodes(ctx: &SolverCtx<'_>, p: &ParamMap, keys: &[&str], e: ElemId) -> Result<Vec<f64>, SolverError> {
    for k in keys {
        if let Some(v) = ctx.real_at_nodes(p, k, e)? {
            return Ok(v);
        }
    }
    Ok(vec![0.0; ctx.mesh.element(e).nodes.len()])
}

fn material(ctx: &SolverCtx<'_>, e: ElemId) -> Result<ParamMap, SolverError> {
    ctx.body_params(e)
        .and_then(|b| b.material)
        .cloned()
        .ok_or_else(|| SolverError::fatal("RefineMesh", "the element has no material"))
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

/// `Element % TYPE % NodeU/V/W` of a linear triangle or quadrilateral.
fn node_uvw(etype: ElementType, k: usize) -> [f64; 3] {
    match etype.0 {
        303 => [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]][k],
        _ => [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]][k],
    }
}

/// The parent-local point of the edge point with basis `eb`:
/// `u = SUM( EdgeBasis(1:en) * x(1:en) )`.
fn parent_uvw(etype: ElementType, pnodes: &[u32], edge_nodes: &[u32], eb: &[f64]) -> [f64; 3] {
    let loc: Vec<[f64; 3]> = edge_nodes
        .iter()
        .map(|n| node_uvw(etype, pnodes.iter().position(|q| q == n).expect("an edge of its parent")))
        .collect();
    std::array::from_fn(|d| {
        let f: Vec<f64> = loc.iter().map(|c| c[d]).collect();
        sum(eb, &f)
    })
}

fn coords(ctx: &SolverCtx<'_>, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&k| ctx.mesh.coord(k)).collect()
}

impl Estimator for DiffusionEstimator {
    fn inside(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let el = ctx.mesh.element(e);
        let nodes = el.nodes.to_vec();
        let n = nodes.len();
        let Some(pot) = potential(ctx, &nodes) else { return Ok((0.0, [0.0; 2])) };
        let dim = ctx.mesh.dim();
        let xyz = coords(ctx, &nodes);
        let mat = material(ctx, e)?;
        let coef = first_at_nodes(ctx, &mat, self.coefficient, e)?;
        let src = match ctx.body_params(e).and_then(|b| b.body_force).cloned() {
            Some(bf) => first_at_nodes(ctx, &bf, self.source, e)?,
            None => vec![0.0; n],
        };
        let (mut fnorm, mut rnorm) = (0.0, 0.0);
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &xyz, uvw)?;
            let b = &ev.basis[..n];
            let s = w * ev.det_j;
            let mut residual = -sum(&src, b);
            for j in 0..dim {
                let dbj: Vec<f64> = ev.dbasis_dx[..n].iter().map(|d| d[j]).collect();
                residual -= sum(&pot, &dbj) * sum(&coef, &dbj);
                // `Coefficient * SUM( Potential * ddBasisddx(:,j,j) )`:
                // zero on linear elements.
            }
            for _ in 0..dim {
                fnorm += s * sum(&src, b).powi(2);
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
        let exyz = coords(ctx, &en);
        // `Grad = 0` once, before the points: a side out of the system
        // keeps the previous point's gradient.
        let mut grad = [[0.0f64; 2]; 3];
        let (mut len, mut rnorm) = (0.0, 0.0);
        // `NormalVector( Edge, EdgeNodes, u, v, .FALSE. )` of a line.
        let (dxdu, dydu) = (0.0 + -0.5 * exyz[0][0] + 0.5 * exyz[1][0], 0.0 + -0.5 * exyz[0][1] + 0.5 * exyz[1][1]);
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
                let Some(pot) = potential(ctx, &pnodes) else { continue };
                let pxyz = coords(ctx, &pnodes);
                let pev = basis_at(ctx, pel.etype, &pxyz, parent_uvw(pel.etype, &pnodes, &en, &eb))?;
                let mat = material(ctx, pe)?;
                let coef_p = first_at_nodes(ctx, &mat, self.coefficient, pe)?;
                let coef_e: Vec<f64> = en
                    .iter()
                    .map(|n| coef_p[pnodes.iter().position(|q| q == n).expect("edge node")])
                    .collect();
                let coefficient = sum(&coef_e, &eb);
                for j in 0..dim {
                    let dbj: Vec<f64> = pev.dbasis_dx[..pnodes.len()].iter().map(|d| d[j]).collect();
                    grad[j][i] = coefficient * sum(&dbj, &pot);
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
        let [left, right] = el.parents;
        let in_system = |p: ElemId| potential(ctx, ctx.mesh.element(p).nodes).is_some();
        let pe = match left {
            None => right,
            Some(l) if !in_system(l) => right,
            Some(l) => Some(l),
        };
        let Some(pe) = pe else { return Ok((0.0, [0.0; 2])) };
        let pel = ctx.mesh.element(pe);
        let pnodes = pel.nodes.to_vec();
        let Some(pot) = potential(ctx, &pnodes) else { return Ok((0.0, [0.0; 2])) };
        let Some(bc) = ctx.sim.model.bc_of_tag(el.tag).map(|i| ctx.sim.model.boundary_conditions[i - 1].clone()) else {
            return Ok((0.0, [0.0; 2]));
        };
        let var = ctx.primary().map(|v| v.name.to_string()).unwrap_or_default();
        let dirichlet = bc.get(&var).is_some()
            || bc.0.keys().any(|k| elmer_core::normalize_key(k).starts_with("constraint mode"));
        let flux = ctx.real_at_nodes(&bc, "electric flux", e)?.unwrap_or_else(|| vec![0.0; en.len()]);
        // The parent's material at the edge's nodes.
        let mat = material(ctx, pe)?;
        let coef_p = first_at_nodes(ctx, &mat, self.coefficient, pe)?;
        let coef: Vec<f64> = en
            .iter()
            .map(|n| coef_p[pnodes.iter().position(|q| q == n).expect("boundary node of its parent")])
            .collect();
        let dim = ctx.mesh.dim();
        let exyz = coords(ctx, &en);
        let pxyz = coords(ctx, &pnodes);
        let (mut gnorm, mut len, mut rnorm) = (0.0, 0.0, 0.0);
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &exyz, uvw)?;
            let eb = ev.basis[..en.len()].to_vec();
            let s = w * ev.det_j;
            let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
            let pev = basis_at(ctx, pel.etype, &pxyz, parent_uvw(pel.etype, &pnodes, &en, &eb))?;
            let coefficient = sum(&coef, &eb);
            let mut residual = -sum(&flux, &eb);
            for k in 0..dim {
                let dbk: Vec<f64> = pev.dbasis_dx[..pnodes.len()].iter().map(|d| d[k]).collect();
                let g = coefficient * sum(&dbk, &pot) * normal[k];
                residual += g;
                gnorm += s * g.powi(2);
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
