//! The elasticity residual estimators of Stress.F90
//! (`ElasticityInsideResidual`, `ElasticityEdgeResidual`,
//! `ElasticityBoundaryResidual`) that ElasticSolve and StressSolve give
//! RefineMesh, in Cartesian coordinates on linear elements with isotropic
//! Youngs modulus.
//!
//! Ported: the large-deflection stress (`LameParameters`,
//! `LargeDeflectionResidualStress`), which ElasticSolve uses by default, and
//! the small-strain `LocalStress` of StressSolve (isotropic, compressible,
//! no thermal strain: the estimators pass a zero temperature).
//! Refused: anisotropic or tensor moduli, `Incompressible`, transient runs,
//! other coordinate systems and nodal `Stress 1..6` fields.

use elmer_core::{ElemId, ElementType, MeshView, ParamMap, Params};
use elmer_solver::{SolverCtx, SolverError};

use crate::{Edge, Estimator};

type M3 = [[f64; 3]; 3];

/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("elasticity adaptivity: {what}"))
}

fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `LameParameters`.
fn lame(young: f64, poisson: f64, plane_stress: bool) -> (f64, f64) {
    let l1 = if plane_stress {
        young * poisson / (1.0 - poisson.powi(2))
    } else {
        young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    };
    (l1, young / (2.0 * (1.0 + poisson)))
}

/// `LargeDeflectionResidualStress`: `(Stress, Strain)` from the gradient.
fn ld_stress(g: &M3, l1: f64, l2: f64, dim: usize) -> (M3, M3) {
    let mut strain = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut m = 0.0;
            for k in 0..3 {
                m += g[k][i] * g[k][j];
            }
            strain[i][j] = (g[j][i] + g[i][j] + m) / 2.0;
        }
    }
    let mut tr = 0.0;
    for i in 0..dim {
        tr += strain[i][i];
    }
    let mut s2 = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let id = if i == j { 1.0 } else { 0.0 };
            s2[i][j] = 2.0 * l2 * strain[i][j] + l1 * tr * id;
        }
    }
    let mut stress = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                let defg = if i == k { 1.0 } else { 0.0 } + g[i][k];
                s += defg * s2[k][j];
            }
            stress[i][j] = s;
        }
    }
    (stress, strain)
}

/// `BuildIsotropicC` (Stress.F90) for Cartesian coordinates, with the
/// out-of-plane row `C(4,1:2) = C(1,2)` LocalStress adds in plane strain.
fn isotropic_c(young: f64, poisson: f64, dim: usize, plane_stress: bool) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    let f = if dim == 2 && plane_stress {
        c[0][0] = 1.0;
        c[0][1] = poisson;
        c[1][0] = poisson;
        c[1][1] = 1.0;
        c[2][2] = 0.5 * (1.0 - poisson);
        1.0 - poisson * poisson
    } else if dim == 2 {
        c[0][0] = 1.0 - poisson;
        c[0][1] = poisson;
        c[1][0] = poisson;
        c[1][1] = 1.0 - poisson;
        c[2][2] = 0.5 - poisson;
        (1.0 + poisson) * (1.0 - 2.0 * poisson)
    } else {
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] = if i == j { 1.0 - poisson } else { poisson };
            }
            c[3 + i][3 + i] = 0.5 - poisson;
        }
        (1.0 + poisson) * (1.0 - 2.0 * poisson)
    };
    // `C = C * Young / ( ... )`.
    for row in c.iter_mut() {
        for x in row.iter_mut() {
            *x = *x * young / f;
        }
    }
    if dim == 2 && !plane_stress {
        c[3][0] = c[0][1];
        c[3][1] = c[0][1];
    }
    c
}

/// `LocalStress` (Stress.F90) of an isotropic compressible material at zero
/// temperature in Cartesian coordinates: `(Stress, Strain)` from the
/// gradient `LGrad(1:dim,1:dim)`.
fn local_stress(g: &M3, young: f64, poisson: f64, plane_stress: bool, dim: usize) -> (M3, M3) {
    let c = isotropic_c(young, poisson, dim, plane_stress);
    let mut lg = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            lg[i][j] = g[i][j];
        }
    }
    let mut strain = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            strain[i][j] = (lg[i][j] + lg[j][i]) / 2.0;
        }
    }
    // `Strain2Stress`.
    let (sv, map): (Vec<f64>, &[(usize, usize)]) = if dim == 2 {
        (vec![strain[0][0], strain[1][1], strain[0][1] * 2.0], &[(0, 0), (1, 1), (0, 1)])
    } else {
        (
            vec![
                strain[0][0],
                strain[1][1],
                strain[2][2],
                strain[0][1] * 2.0,
                strain[1][2] * 2.0,
                strain[0][2] * 2.0,
            ],
            &[(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)],
        )
    };
    let mut stress = [[0.0f64; 3]; 3];
    for (i, &(p, q)) in map.iter().enumerate() {
        let mut csum = 0.0;
        for (j, x) in sv.iter().enumerate() {
            csum += c[i][j] * x;
        }
        stress[p][q] = csum;
        stress[q][p] = csum;
    }
    if dim == 2 {
        if !plane_stress {
            // `PlaneStrainStressZZ`.
            let s = [strain[0][0], strain[1][1], 2.0 * strain[0][1]];
            let mut szz = 0.0;
            for j in 0..3 {
                szz += c[3][j] * s[j];
            }
            stress[2][2] += szz;
        } else {
            strain[2][2] = -(poisson / (1.0 - poisson)) * (strain[0][0] + strain[1][1]);
        }
    }
    (stress, strain)
}

/// `ResGrad = MATMUL( NodalDisplacement(:,1:nd), dBasisdx(1:nd,:) )`.
fn grad(u: &[[f64; 3]], db: &[[f64; 3]]) -> M3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut s = 0.0;
            for (uk, dk) in u.iter().zip(db) {
                s += uk[i] * dk[j];
            }
            s
        })
    })
}

/// The estimator of a displacement solver.
pub struct ElasticityEstimator {
    pub large_deflection: bool,
}

struct Parent {
    nodes: Vec<u32>,
    etype: ElementType,
    xyz: Vec<[f64; 3]>,
    young: Vec<f64>,
    poisson: Vec<f64>,
    plane_stress: bool,
    disp: Vec<[f64; 3]>,
}

fn displacement(ctx: &SolverCtx<'_>, nodes: &[u32]) -> Option<Vec<[f64; 3]>> {
    let v = ctx.primary()?;
    let perm = v.perm?;
    let dofs = v.dofs;
    nodes
        .iter()
        .map(|&n| {
            let k = perm[n as usize];
            (k != elmer_core::NO_DOF).then(|| {
                let mut u = [0.0f64; 3];
                for d in 0..dofs.min(3) {
                    u[d] = v.values[dofs * k as usize + d];
                }
                u
            })
        })
        .collect()
}

fn at_nodes(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str, e: ElemId) -> Result<Vec<f64>, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    Ok(ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]))
}

fn parent(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Option<Parent>, SolverError> {
    let el = ctx.mesh.element(e);
    let nodes = el.nodes.to_vec();
    let Some(disp) = displacement(ctx, &nodes) else { return Ok(None) };
    let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal("Elasticity", "no body"))?;
    let mat = bp.material.cloned().ok_or_else(|| SolverError::fatal("Elasticity", "no material"))?;
    let plane_stress = bp.equation.and_then(|q| q.logical("plane stress")).unwrap_or(false);
    if matches!(mat.get("youngs modulus"), Some(elmer_core::Value::Real { rows, cols, .. }) if *rows * *cols > 1) {
        return Err(unsupported("a tensor Youngs Modulus (InputTensor)"));
    }
    Ok(Some(Parent {
        young: at_nodes(ctx, &mat, "youngs modulus", e)?,
        poisson: at_nodes(ctx, &mat, "poisson ratio", e)?,
        xyz: nodes.iter().map(|&k| ctx.mesh.coord(k)).collect(),
        nodes,
        etype: el.etype,
        plane_stress,
        disp,
    }))
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

fn node_uvw(etype: ElementType, k: usize) -> [f64; 3] {
    match etype.0 {
        303 => [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]][k],
        _ => [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]][k],
    }
}

fn parent_uvw(p: &Parent, edge_nodes: &[u32], eb: &[f64]) -> [f64; 3] {
    let loc: Vec<[f64; 3]> = edge_nodes
        .iter()
        .map(|n| node_uvw(p.etype, p.nodes.iter().position(|q| q == n).expect("an edge of its parent")))
        .collect();
    std::array::from_fn(|d| {
        let f: Vec<f64> = loc.iter().map(|c| c[d]).collect();
        sum(eb, &f)
    })
}

impl ElasticityEstimator {
    /// The stress of `p` with `ev` (its basis at a point).
    fn stress(&self, p: &Parent, ev: &elmer_core::BasisEval, dim: usize) -> Result<(M3, M3), SolverError> {
        let n = p.nodes.len();
        let b = &ev.basis[..n];
        let g = grad(&p.disp, &ev.dbasis_dx[..n]);
        if !self.large_deflection {
            // `Poisson = SUM( Basis(1:n) * PoissonRatio(1:n) )`, likewise Young.
            return Ok(local_stress(&g, sum(b, &p.young), sum(b, &p.poisson), p.plane_stress, dim));
        }
        let (l1, l2) = lame(sum(&p.young, b), sum(&p.poisson, b), p.plane_stress);
        Ok(ld_stress(&g, l1, l2, dim))
    }

    pub fn check(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
        let sim = ctx.simulation_params();
        if let Some(cs) = sim.string("coordinate system") {
            if !elmer_core::normalize_key(cs).starts_with("cartesian") {
                return Err(unsupported("a non-Cartesian coordinate system"));
            }
        }
        if sim.string("simulation type").is_some_and(|s| elmer_core::normalize_key(s) == "transient") {
            return Err(unsupported("transient runs"));
        }
        if ctx.params().logical("incompressible") == Some(true) {
            return Err(unsupported("Incompressible"));
        }
        if ctx.sim.vars.get("stress 1").is_some() {
            return Err(unsupported("nodal Stress 1..6 fields"));
        }
        Ok(())
    }
}

impl Estimator for ElasticityEstimator {
    fn inside(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let Some(p) = parent(ctx, e)? else { return Ok((0.0, [0.0; 2])) };
        let dim = ctx.mesh.dim();
        let n = p.nodes.len();
        let bf = ctx.body_params(e).and_then(|b| b.body_force).cloned();
        let mut nf = [vec![0.0f64; n], vec![0.0f64; n], vec![0.0f64; n]];
        if let Some(bf) = &bf {
            let f1 = at_nodes(ctx, bf, "stress bodyforce 1", e)?;
            let f2 = at_nodes(ctx, bf, "stress bodyforce 2", e)?;
            let f3 = at_nodes(ctx, bf, "stress bodyforce 3", e)?;
            // `NodalForce(2) = NodalForce(1) + ...`, `(3) = NodalForce(1) + ...`.
            for i in 0..n {
                nf[0][i] += f1[i];
            }
            for i in 0..n {
                nf[1][i] = nf[0][i] + f2[i];
            }
            for i in 0..n {
                nf[2][i] = nf[0][i] + f3[i];
            }
        }
        // The stress at the element's nodes.
        let mut nodal = Vec::with_capacity(n);
        for i in 0..n {
            let ev = basis_at(ctx, p.etype, &p.xyz, node_uvw(p.etype, i))?;
            nodal.push(self.stress(&p, &ev, dim)?.0);
        }
        let (mut rnorm, mut area, mut energy, mut young) = (0.0, 0.0, 0.0, 0.0);
        for (uvw, w) in gauss(ctx, p.etype)? {
            let ev = basis_at(ctx, p.etype, &p.xyz, uvw)?;
            let b = &ev.basis[..n];
            let s = w * ev.det_j;
            let mut res = [0.0f64; 3];
            for i in 0..3 {
                res[i] = -sum(&nf[i], b);
                for j in 0..3 {
                    let sij: Vec<f64> = nodal.iter().map(|st| st[i][j]).collect();
                    let dbj: Vec<f64> = ev.dbasis_dx[..n].iter().map(|d| d[j]).collect();
                    res[i] -= sum(&sij, &dbj);
                }
            }
            young += s * sum(&p.young, b);
            let (stress, strain) = self.stress(&p, &ev, dim)?;
            let mut dd = 0.0;
            for i in 0..dim {
                for j in 0..dim {
                    dd += strain[i][j] * stress[i][j];
                }
            }
            energy += s * dd / 2.0;
            area += s;
            let mut r2 = 0.0;
            for x in &res[..dim] {
                r2 += x.powi(2);
            }
            rnorm += s * r2;
        }
        let young = young / area;
        let ind = area * rnorm / young;
        Ok((energy, [ind, ind]))
    }

    fn edge(&mut self, ctx: &SolverCtx<'_>, edge: &Edge) -> Result<[f64; 2], SolverError> {
        let dim = ctx.mesh.dim();
        let en = edge.nodes.to_vec();
        let exyz: Vec<[f64; 3]> = en.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let (dxdu, dydu) = (0.0 + -0.5 * exyz[0][0] + 0.5 * exyz[1][0], 0.0 + -0.5 * exyz[0][1] + 0.5 * exyz[1][1]);
        let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
        let normal = [-dydu * det, dxdu * det, 0.0];
        let sides: Vec<Option<Parent>> = [edge.left, edge.right]
            .into_iter()
            .map(|s| s.map_or(Ok(None), |p| parent(ctx, ElemId(p as u32))))
            .collect::<Result<_, _>>()?;
        let (mut len, mut rnorm, mut young) = (0.0, 0.0, 0.0);
        // The last side's modulus and basis carry over (the Fortran's locals).
        let mut last: Option<(Vec<f64>, Vec<f64>)> = None;
        for (uvw, w) in gauss(ctx, ElementType(202))? {
            let ev = basis_at(ctx, ElementType(202), &exyz, uvw)?;
            let eb = ev.basis[..2].to_vec();
            let s = w * ev.det_j;
            let mut st = [[[0.0f64; 3]; 3]; 2];
            for (i, side) in sides.iter().enumerate() {
                let Some(p) = side else { continue };
                let pev = basis_at(ctx, p.etype, &p.xyz, parent_uvw(p, &en, &eb))?;
                st[i] = self.stress(p, &pev, dim)?.0;
                last = Some((p.young.clone(), pev.basis[..p.nodes.len()].to_vec()));
            }
            len += s;
            let mut j2 = 0.0;
            for i in 0..dim {
                let mut jump = 0.0;
                for j in 0..3 {
                    jump += (st[0][i][j] - st[1][i][j]) * normal[j];
                }
                j2 += jump.powi(2);
            }
            rnorm += s * j2;
            if let Some((yv, b)) = &last {
                young += s * sum(yv, b);
            }
        }
        let young = young / len;
        let ind = len * rnorm / young;
        Ok([ind, ind])
    }

    fn boundary(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError> {
        let el = ctx.mesh.element(e);
        let en = el.nodes.to_vec();
        let [left, right] = el.parents;
        let in_system = |p: ElemId| displacement(ctx, ctx.mesh.element(p).nodes).is_some();
        let pe = match left {
            None => right,
            Some(l) if !in_system(l) => right,
            Some(l) => Some(l),
        };
        let Some(pe) = pe else { return Ok((0.0, [0.0; 2])) };
        let Some(p) = parent(ctx, pe)? else { return Ok((0.0, [0.0; 2])) };
        let Some(bc) = ctx.sim.model.bc_of_tag(el.tag).map(|i| ctx.sim.model.boundary_conditions[i - 1].clone()) else {
            return Ok((0.0, [0.0; 2]));
        };
        let dim = ctx.mesh.dim();
        let exyz: Vec<[f64; 3]> = en.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let force = [
            at_nodes(ctx, &bc, "force 1", e)?,
            at_nodes(ctx, &bc, "force 2", e)?,
            at_nodes(ctx, &bc, "force 3", e)?,
        ];
        let ext = at_nodes(ctx, &bc, "normal force", e)?;
        let mut dir = [1.0f64; 3];
        if bc.get("displacement").is_some() {
            dir = [0.0; 3];
        }
        for k in 0..3 {
            if bc.get(&format!("displacement {}", k + 1)).is_some() {
                dir[k] = 0.0;
            }
        }
        let (mut len, mut rnorm, mut young) = (0.0, 0.0, 0.0);
        for (uvw, w) in gauss(ctx, el.etype)? {
            let ev = basis_at(ctx, el.etype, &exyz, uvw)?;
            let eb = ev.basis[..en.len()].to_vec();
            let s = w * ev.det_j;
            let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
            let pev = basis_at(ctx, p.etype, &p.xyz, parent_uvw(&p, &en, &eb))?;
            let (stress, _) = self.stress(&p, &pev, dim)?;
            let pext = sum(&ext, &eb);
            let mut res = [0.0f64; 3];
            for i in 0..3 {
                res[i] = sum(&force[i], &eb) - pext * normal[i];
                let mut fs = 0.0;
                for j in 0..3 {
                    fs += stress[i][j] * normal[j];
                }
                res[i] -= fs * dir[i];
            }
            len += s;
            let mut r2 = 0.0;
            for x in &res[..dim] {
                r2 += x.powi(2);
            }
            rnorm += s * r2;
            young += s * sum(&p.young, &pev.basis[..p.nodes.len()]);
        }
        let mut ind = 0.0;
        if young > AEPS {
            let young = young / len;
            ind = len * rnorm / young;
        }
        Ok((0.0, [ind, ind]))
    }
}

#[allow(dead_code)]
fn _mesh(_: &dyn MeshView) {}
