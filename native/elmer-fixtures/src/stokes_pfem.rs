//! `StokesPFEM/Stokes.f90`: `SUBROUTINE StokesSolver`, steady
//! Navier-Stokes on p elements with the velocity components and the
//! pressure in one `dim + 1` dof variable (`Procedure = "Stokes1"
//! "StokesSolver"`).
//!
//! Velocity uses all the element's dofs; pressure only the nodal ones:
//! the pressure rows and columns of the edge and face dofs are replaced by
//! the identity (P(k+1)/P1 pairs). Convection is Picard from the variable's
//! current values (`Convect`, default on), Newton for good once successive
//! norms are within 1 % (`SAVE Newton, Norm`). The integration rule is
//! `GaussPoints(Element, nd + nb)`. Local bubbles (`LCondensate` of the
//! velocity bubbles) and transient runs are refused.

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, QuadPoint};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{check_element_definitions, nodal};

const CALLER: &str = "StokesSolve";

/// Nodal inputs of `LocalMatrix`: `n` values per property, velocity at the
/// element's `nd` dofs.
pub struct StokesNodal<'a> {
    pub load: &'a [[f64; 4]],
    pub rho: &'a [f64],
    pub mu: &'a [f64],
    pub velo: &'a [[f64; 3]],
}

/// `LocalMatrix(MASS, STIFF, FORCE, LOAD, rho, mu, Velocity, Element, n,
/// nd, ntot, dim)` into `local` (`(dim+1) * ntot` square, interleaved).
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub fn local_matrix(
    etype: ElementType,
    xyz: &[[f64; 3]],
    ntot: usize,
    v: &StokesNodal<'_>,
    convect: bool,
    newton: bool,
    dim: usize,
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    local: &mut LocalSystem,
) -> Result<(), String> {
    let n = xyz.len();
    let nd = v.velo.len();
    let c = dim + 1;
    let nn = c * ntot;
    *local = LocalSystem::zeros(nn);
    let mut mass = vec![0.0; nn * nn];
    let mut ev = BasisEval::default();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).fold(0.0, |s, (x, y)| s + x * y);
    for qp in rule {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < ntot {
            return Err(format!("ElementInfo gave {} basis functions for {ntot} dofs", ev.basis.len()));
        }
        let (b, g) = (&ev.basis[..ntot], &ev.dbasis_dx[..ntot]);
        let s = qp.weight * ev.det_j;
        let mu = dot(&b[..n], v.mu);
        let rho = dot(&b[..n], v.rho);
        let mut velo = [0.0; 3];
        let mut grad = [[0.0; 3]; 3];
        for i in 0..dim {
            velo[i] = (0..nd).fold(0.0, |acc, k| acc + v.velo[k][i] * b[k]);
            for j in 0..dim {
                grad[i][j] = (0..nd).fold(0.0, |acc, k| acc + v.velo[k][i] * g[k][j]);
            }
        }
        let mut load_ip = [0.0; 4];
        for (j, l) in load_ip.iter_mut().enumerate().take(c) {
            *l = (0..n).fold(0.0, |acc, k| acc + b[k] * v.load[k][j]);
        }
        if convect && newton {
            for i in 0..dim {
                load_ip[i] += rho * (0..dim).fold(0.0, |acc, j| acc + grad[i][j] * velo[j]);
            }
        }
        for p in 0..ntot {
            for q in 0..ntot {
                let at = |i: usize, j: usize| (c * p + i) * nn + c * q + j;
                for i in 0..dim {
                    mass[at(i, i)] += s * rho * b[q] * b[p];
                    for j in 0..dim {
                        local.stiff[at(i, i)] += s * mu * g[q][j] * g[p][j];
                        local.stiff[at(i, j)] += s * mu * g[q][i] * g[p][j];
                        if convect {
                            local.stiff[at(i, i)] += s * rho * velo[j] * g[q][j] * b[p];
                            if newton {
                                local.stiff[at(i, j)] += s * rho * grad[i][j] * b[q] * b[p];
                            }
                        }
                    }
                    local.stiff[at(i, dim)] -= s * b[q] * g[p][i];
                    local.stiff[at(dim, i)] += s * rho * g[q][i] * b[p];
                }
            }
            for k in 0..c {
                local.force[c * p + k] += s * load_ip[k] * b[p];
            }
        }
    }
    // Pressure only on the nodes: identity rows for the other dofs.
    for p in n..ntot {
        let i = c * (p + 1) - 1;
        local.force[i] = 0.0;
        for k in 0..nn {
            mass[k * nn + i] = 0.0;
            mass[i * nn + k] = 0.0;
            local.stiff[i * nn + k] = 0.0;
            local.stiff[k * nn + i] = 0.0;
        }
        local.stiff[i * nn + i] = 1.0;
    }
    local.mass = Some(mass);
    Ok(())
}

pub struct StokesSolver {
    newton: bool,
    norm: f64,
    local: LocalSystem,
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StokesSolver { newton: false, norm: 0.0, local: LocalSystem::zeros(0) }))
}

impl Solver for StokesSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_element_definitions(ctx)?;
        let dim = crate::hooks::coordinate_system(ctx)?.dim;
        if ctx.primary().is_none_or(|v| v.dofs != dim + 1) {
            return Err(SolverError::Unsupported(format!("Stokes1 StokesSolver: a Variable of {} dofs", dim + 1)));
        }
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::NotWired(
                "Stokes1 StokesSolver: its p elements need a dof layout (elmer-pdofs in Backends::dofs)".into(),
            ));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("Stokes1 StokesSolver: transient runs".into()));
        }
        let dim = crate::hooks::coordinate_system(ctx)?.dim;
        let convect = ctx.params().logical("convect").unwrap_or(true);
        let fatal = |m: String| SolverError::fatal(CALLER, m);
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for e in active {
            let d = ctx.solver_dofs().ok_or_else(|| fatal("no dof layout".into()))?;
            d.element_slots(e, &mut slots);
            let (nd, nb) = (slots.len(), d.local_bubbles(e));
            if nb > 0 {
                return Err(SolverError::Unsupported("Stokes1 StokesSolver: local bubbles (LCondensate)".into()));
            }
            let order = d.element_nodes(e);
            let basis = d.element_basis(e).map_err(fatal)?;
            let rule = d.element_rule_np(e, nd + nb).map_err(fatal)?;
            let etype = ctx.mesh.element(e).etype;
            let xyz: Vec<[f64; 3]> = order.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let n = order.len();
            let (m, bf) = ctx.body_params(e).map_or((None, None), |b| (b.material, b.body_force));
            let mut load = vec![[0.0; 4]; n];
            for (k, key) in ["source x", "source y", "source z"].into_iter().enumerate() {
                for (l, x) in load.iter_mut().zip(nodal(ctx, bf, key, e, &order)?) {
                    l[k] = x;
                }
            }
            let rho = nodal(ctx, m, "density", e, &order)?;
            let mu = nodal(ctx, m, "viscosity", e, &order)?;
            // GetVectorLocalSolution: the variable at the element's dofs.
            let mut velo = vec![[0.0; 3]; nd];
            if convect {
                let v = ctx.primary().ok_or_else(|| fatal("no solver variable".into()))?;
                for (o, &sl) in velo.iter_mut().zip(&slots) {
                    for (k, ok) in o.iter_mut().enumerate().take(dim) {
                        *ok = v.values[sl as usize * v.dofs + k];
                    }
                }
            }
            let nodal_in = StokesNodal { load: &load, rho: &rho, mu: &mu, velo: &velo };
            local_matrix(etype, &xyz, nd + nb, &nodal_in, convect, self.newton, dim, basis.as_ref(), &rule, &mut self.local)
                .map_err(fatal)?;
            self.local.mass = None;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        let prev = self.norm;
        self.norm = ctx.default_solve()?;
        if 2.0 * (prev - self.norm).abs() / (prev + self.norm) < 1.0e-2 {
            self.newton = true;
        }
        Ok(())
    }
}
