//! `StokesProj/Velocity.f90` and `StokesProj/Pressure.f90`: a projection
//! scheme for steady Navier-Stokes split into two solver procedures that
//! the coupled-system loop alternates.
//!
//! * `PressureSolver` solves `-div grad dp = -div u` for its own variable
//!   (the pressure correction) from the velocity field `VelocityTot`, then
//!   adds `Pressure Relax` times it to the exported `PressureTot`.
//! * `VelocitySolver` solves the momentum equation for its variable with
//!   `PressureTot` on the right-hand side, Picard convection from its own
//!   previous iterate (`Convect`, default on), and switches to Newton
//!   linearisation for good once two successive norms are within 1 %
//!   (`SAVE Newton, Norm`).
//!
//! Both integrate with `GaussPoints(Element, ntot)` and keep the Fortran's
//! association order. Bubbles (`LCondensate`) are refused: only the `n:1`
//! path is ported.

use elmer_core::{BasisEval, ElemId, ElementBasis, FieldRef, LocalSystem, ParamMap, Params, Quadrature, NO_DOF};
use elmer_fe::ElmerQuadrature;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{element_nodes, lagrange_basis};

fn check_nodal(params: &ParamMap, which: &str) -> Result<(), SolverError> {
    if params.string("element").is_some_and(|e| elmer_core::normalize_key(e) != "n:1")
        || params.logical("bubbles") == Some(true)
    {
        return Err(SolverError::Unsupported(format!("StokesProj {which}: bubbles or p-elements (LCondensate)")));
    }
    Ok(())
}

/// `GetVectorLocalSolution` / `GetScalarLocalSolution`: the field's
/// components at the element's nodes, zero where it has no dof.
fn local_solution(field: Option<FieldRef<'_>>, nodes: &[u32], dim: usize) -> Vec<[f64; 3]> {
    let mut out = vec![[0.0; 3]; nodes.len()];
    let Some(f) = field else { return out };
    for (o, &n) in out.iter_mut().zip(nodes) {
        let slot = match f.perm {
            Some(perm) => match perm.get(n as usize) {
                Some(&s) if s != NO_DOF => s as usize,
                _ => continue,
            },
            None => n as usize,
        };
        for (k, ok) in o.iter_mut().enumerate().take(dim.min(f.dofs)) {
            *ok = f.values.get(slot * f.dofs + k).copied().unwrap_or(0.0);
        }
    }
    out
}

/// `GetReal(list, key)` at the nodes, zeros when absent.
fn nodal(ctx: &SolverCtx<'_>, list: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    match list {
        Some(p) => Ok(ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n])),
        None => Ok(vec![0.0; n]),
    }
}

/// `SUM(a * b)`, accumulated from zero in order.
fn dot(a: impl IntoIterator<Item = f64>, b: impl IntoIterator<Item = f64>) -> f64 {
    a.into_iter().zip(b).fold(0.0, |s, (x, y)| s + x * y)
}

/// The pressure correction solver.
pub struct PressureSolver {
    local: LocalSystem,
}

pub fn pressure_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    check_nodal(params, "PressureSolver")?;
    Ok(Box::new(PressureSolver { local: LocalSystem::zeros(0) }))
}

/// `PressureSolver`'s `LocalMatrix`: `STIFF += s grad.grad`,
/// `FORCE -= s divU Basis`.
pub fn pressure_local(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    velo: &[[f64; 3]],
    dim: usize,
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let nd = xyz.len();
    *local = LocalSystem::zeros(nd);
    let mut ev = BasisEval::default();
    for qp in ElmerQuadrature.gauss_points(etype, Some(nd)) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let (b, g) = (&ev.basis[..nd], &ev.dbasis_dx[..nd]);
        let s = qp.weight * ev.det_j;
        let mut div_u = 0.0;
        for i in 0..dim {
            div_u += dot(velo.iter().map(|v| v[i]), g.iter().map(|d| d[i]));
        }
        for p in 0..nd {
            for q in 0..nd {
                let m = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                local.stiff[p * nd + q] += s * m;
            }
        }
        for (f, bp) in local.force.iter_mut().zip(b) {
            *f -= s * div_u * bp;
        }
    }
    Ok(())
}

impl Solver for PressureSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::Unsupported("StokesProj PressureSolver: a scalar Variable".into()));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let basis = lagrange_basis(ctx)?;
        let dim = basis.cdim;
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        for e in active {
            let el = ctx.mesh.element(e);
            let etype = el.etype;
            let velo = local_solution(ctx.variable("velocitytot"), el.nodes, dim);
            let xyz = element_nodes(ctx.mesh, e);
            pressure_local(etype, &xyz, &velo, dim, &basis, &mut self.local)
                .map_err(|m| SolverError::fatal("PressureSolve", m))?;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;

        let relax = ctx.creal(ctx.params(), "pressure relax")?.unwrap_or(1.0);
        let delta: Vec<f64> = ctx.primary().map(|v| v.values.to_vec()).unwrap_or_default();
        // `PressureVariable % Values = PressureVariable % Values + Relax * Solver % Variable % Values`
        let total = ctx
            .sim
            .vars
            .get_mut("pressuretot")
            .ok_or_else(|| SolverError::fatal("PressureSolve", "PressureTot not found"))?;
        if total.values.len() != delta.len() {
            return Err(SolverError::fatal("PressureSolve", "PressureTot and the solver variable differ in size"));
        }
        for (t, d) in total.values.iter_mut().zip(&delta) {
            *t += relax * d;
        }
        total.values_changed = true;
        Ok(())
    }
}

/// The momentum solver, with the `SAVE`d `Newton` switch and `Norm`.
pub struct VelocitySolver {
    newton: bool,
    norm: f64,
    local: LocalSystem,
}

pub fn velocity_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    check_nodal(params, "VelocitySolver")?;
    Ok(Box::new(VelocitySolver { newton: false, norm: 0.0, local: LocalSystem::zeros(0) }))
}

/// Nodal inputs of `VelocitySolver`'s `LocalMatrix`.
pub struct VelocityNodal<'a> {
    pub load: &'a [[f64; 3]],
    pub rho: &'a [f64],
    pub mu: &'a [f64],
    pub velo: &'a [[f64; 3]],
    pub pressure: &'a [f64],
}

/// `VelocitySolver`'s `LocalMatrix`, interleaved `dim` dofs per node.
pub fn velocity_local(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    v: &VelocityNodal<'_>,
    newton: bool,
    dim: usize,
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let nd = xyz.len();
    let nn = nd * dim;
    *local = LocalSystem::zeros(nn);
    let mut ev = BasisEval::default();
    for qp in ElmerQuadrature.gauss_points(etype, Some(nd)) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let (b, g) = (&ev.basis[..nd], &ev.dbasis_dx[..nd]);
        let s = qp.weight * ev.det_j;
        let mu = dot(b.iter().copied(), v.mu.iter().copied());
        let rho = dot(b.iter().copied(), v.rho.iter().copied());
        let pressure = dot(v.pressure.iter().copied(), b.iter().copied());
        let mut velo = [0.0; 3];
        let mut grad = [[0.0; 3]; 3];
        let mut load_ip = [0.0; 3];
        for i in 0..dim {
            velo[i] = dot(v.velo.iter().map(|x| x[i]), b.iter().copied());
            for j in 0..dim {
                grad[i][j] = dot(v.velo.iter().map(|x| x[i]), g.iter().map(|d| d[j]));
            }
            load_ip[i] = dot(b.iter().copied(), v.load.iter().map(|x| x[i]));
        }
        if newton {
            for i in 0..dim {
                load_ip[i] += rho * (0..dim).fold(0.0, |acc, j| acc + grad[i][j] * velo[j]);
            }
        }
        for p in 0..nd {
            for q in 0..nd {
                let at = |i: usize, j: usize| (dim * p + i) * nn + dim * q + j;
                for i in 0..dim {
                    for j in 0..dim {
                        local.stiff[at(i, i)] += s * mu * g[q][j] * g[p][j];
                        local.stiff[at(i, j)] += s * mu * g[q][i] * g[p][j];
                        local.stiff[at(i, i)] += s * rho * velo[j] * g[q][j] * b[p];
                        if newton {
                            local.stiff[at(i, j)] += s * rho * grad[i][j] * b[q] * b[p];
                        }
                    }
                }
            }
            for i in 0..dim {
                local.force[dim * p + i] += s * (load_ip[i] * b[p] + pressure * g[p][i]);
            }
        }
    }
    Ok(())
}

impl Solver for VelocitySolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let basis = lagrange_basis(ctx)?;
        let dim = basis.cdim;
        if ctx.primary().is_none_or(|v| v.dofs != dim) {
            return Err(SolverError::Unsupported(format!("StokesProj VelocitySolver: a Variable of {dim} dofs")));
        }
        let convect = ctx.params().logical("convect").unwrap_or(true);
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        for e in active {
            let el = ctx.mesh.element(e);
            let (etype, n) = (el.etype, el.nodes.len());
            let xyz = element_nodes(ctx.mesh, e);
            let (m, bf) = ctx.body_params(e).map_or((None, None), |b| (b.material, b.body_force));
            let mut load = vec![[0.0; 3]; n];
            for (k, key) in ["source x", "source y", "source z"].into_iter().enumerate() {
                for (l, x) in load.iter_mut().zip(nodal(ctx, bf, key, e, n)?) {
                    l[k] = x;
                }
            }
            let rho = nodal(ctx, m, "density", e, n)?;
            let mu = nodal(ctx, m, "viscosity", e, n)?;
            let nodes = ctx.mesh.element(e).nodes;
            let velo = if convect { local_solution(ctx.primary(), nodes, dim) } else { vec![[0.0; 3]; n] };
            let pressure: Vec<f64> =
                local_solution(ctx.variable("pressuretot"), nodes, 1).iter().map(|p| p[0]).collect();
            let v = VelocityNodal { load: &load, rho: &rho, mu: &mu, velo: &velo, pressure: &pressure };
            velocity_local(etype, &xyz, &v, self.newton, dim, &basis, &mut self.local)
                .map_err(|m| SolverError::fatal("StokesSolve", m))?;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        let prev = self.norm;
        self.norm = ctx.default_solve()?;
        let relc = 2.0 * (prev - self.norm).abs() / (prev + self.norm);
        if relc < 1.0e-2 {
            self.newton = true;
        }
        Ok(())
    }
}
