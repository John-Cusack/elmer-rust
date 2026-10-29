//! `ModelPDEhandle/ModelPDEhandle.F90`: `SUBROUTINE AdvDiffSolver`, the
//! ModelPDE advection-diffusion-reaction equation with its properties read
//! through keyword handles (`ListInitElementKeyword` /
//! `ListGetElementReal`) at the integration points.
//!
//! A handle gives a constant keyword as it is and interpolates any other
//! value from the element's nodes with the basis (`SUM(Basis * values)`),
//! and an absent keyword is zero. Bulk: `D grad.grad + C (a.grad u) v + R u v`
//! and the source; boundary: `Robin Coefficient` and `Field Flux +
//! Robin Coefficient * External Field`. The nonlinear loop stops on
//! `NonlinConverged`. Bubbles (`LCondensate`) are not on the `n:1` path and
//! are refused.

use elmer_core::{BasisEval, ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_fe::ElmerQuadrature;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{element_nodes, lagrange_basis};

const CALLER: &str = "AdvDiffSolver";

/// One handle's value on an element.
enum Prop {
    Absent,
    Const(f64),
    Nodal(Vec<f64>),
}

impl Prop {
    fn read(ctx: &SolverCtx<'_>, list: Option<&ParamMap>, key: &str, e: ElemId) -> Result<Self, SolverError> {
        let Some(list) = list else { return Ok(Prop::Absent) };
        match list.get(key) {
            None => Ok(Prop::Absent),
            Some(Value::Real { rows: 1, cols: 1, data }) => Ok(Prop::Const(data[0])),
            Some(_) => Ok(ctx.real_at_nodes(list, key, e)?.map_or(Prop::Absent, Prop::Nodal)),
        }
    }

    /// `ListGetElementReal(handle, Basis, Element, Found)`.
    fn at(&self, basis: &[f64]) -> f64 {
        match self {
            Prop::Absent => 0.0,
            Prop::Const(c) => *c,
            Prop::Nodal(v) => basis.iter().zip(v).fold(0.0, |s, (b, x)| s + b * x),
        }
    }
}

pub struct ModelPdeHandle {
    local: LocalSystem,
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`ModelPdeHandle`].
pub fn factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if params.string("element").is_some_and(|e| elmer_core::normalize_key(e) != "n:1")
        || params.logical("bubbles") == Some(true)
    {
        return Err(SolverError::Unsupported("ModelPDEhandle AdvDiffSolver: bubbles or p-elements (LCondensate)".into()));
    }
    Ok(Box::new(ModelPdeHandle { local: LocalSystem::zeros(0) }))
}

/// The bulk properties of one element.
struct Bulk {
    load: Prop,
    rho: Prop,
    velo: [Prop; 3],
    d: Prop,
    c: Prop,
    r: Prop,
}

/// `LocalMatrix`: stiffness, mass and force over the Gauss points, in the
/// Fortran's order of accumulation.
fn local_matrix(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    props: &Bulk,
    dim: usize,
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let nd = xyz.len();
    *local = LocalSystem::zeros(nd);
    let mut mass = vec![0.0; nd * nd];
    let mut ev = BasisEval::default();
    for qp in ElmerQuadrature.gauss_points(etype, None) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let (b, g) = (&ev.basis[..nd], &ev.dbasis_dx);
        let load = props.load.at(b);
        let rho = props.rho.at(b);
        let mut a = [0.0; 3];
        a[0] = props.velo[0].at(b);
        a[1] = props.velo[1].at(b);
        if dim == 3 {
            a[2] = props.velo[2].at(b);
        }
        let (d, c, r) = (props.d.at(b), props.c.at(b), props.r.at(b));
        let w = qp.weight * ev.det_j;
        for p in 0..nd {
            for q in 0..nd {
                let m = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                local.stiff[p * nd + q] += w * d * m;
            }
        }
        for p in 0..nd {
            for q in 0..nd {
                let adv = (0..dim).fold(0.0, |s, k| s + a[k] * g[q][k]);
                local.stiff[p * nd + q] += w * c * adv * b[p];
                local.stiff[p * nd + q] += w * r * b[q] * b[p];
                mass[p * nd + q] += w * rho * b[q] * b[p];
            }
        }
        for (f, bp) in local.force.iter_mut().zip(b) {
            *f += w * load * bp;
        }
    }
    local.mass = Some(mass);
    Ok(())
}

/// `LocalMatrixBC`.
fn local_matrix_bc(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    (flux, coeff, ext): (&Prop, &Prop, &Prop),
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let nd = xyz.len();
    *local = LocalSystem::zeros(nd);
    let mut ev = BasisEval::default();
    for qp in ElmerQuadrature.gauss_points(etype, None) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let b = &ev.basis[..nd];
        let w = qp.weight * ev.det_j;
        let (f, c, x) = (flux.at(b), coeff.at(b), ext.at(b));
        for p in 0..nd {
            for q in 0..nd {
                local.stiff[p * nd + q] += w * c * b[q] * b[p];
            }
        }
        for (fp, bp) in local.force.iter_mut().zip(b) {
            *fp += w * (f + c * x) * bp;
        }
    }
    Ok(())
}

impl Solver for ModelPdeHandle {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        crate::hooks::check_element_definitions(ctx)?;
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::Unsupported("ModelPDEhandle AdvDiffSolver: a scalar Variable".into()));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let basis = lagrange_basis(ctx)?;
        let dim = basis.cdim;
        let maxiter = ctx.params().integer("nonlinear system max iterations").map_or(1, |n| n.max(1));
        let bulk: Vec<ElemId> = ctx.active_elements().collect();
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();
        for _ in 0..maxiter {
            ctx.default_initialize()?;
            for &e in &bulk {
                let etype = ctx.mesh.element(e).etype;
                let xyz = element_nodes(ctx.mesh, e);
                let (m, bf) = ctx.body_params(e).map_or((None, None), |b| (b.material, b.body_force));
                let props = Bulk {
                    load: Prop::read(ctx, bf, "field source", e)?,
                    rho: Prop::read(ctx, m, "time derivative coefficient", e)?,
                    velo: [
                        Prop::read(ctx, m, "convection velocity 1", e)?,
                        Prop::read(ctx, m, "convection velocity 2", e)?,
                        Prop::read(ctx, m, "convection velocity 3", e)?,
                    ],
                    d: Prop::read(ctx, m, "diffusion coefficient", e)?,
                    c: Prop::read(ctx, m, "convection coefficient", e)?,
                    r: Prop::read(ctx, m, "reaction coefficient", e)?,
                };
                local_matrix(etype, &xyz, &props, dim, &basis, &mut self.local)
                    .map_err(|m| SolverError::fatal(CALLER, m))?;
                if transient {
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            for &e in &boundary {
                // ActiveBoundaryElement: every node in the solver's system.
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let etype = ctx.mesh.element(e).etype;
                let xyz = element_nodes(ctx.mesh, e);
                let bc = ctx.bc_params(e);
                let flux = Prop::read(ctx, bc, "field flux", e)?;
                let coeff = Prop::read(ctx, bc, "robin coefficient", e)?;
                let ext = Prop::read(ctx, bc, "external field", e)?;
                local_matrix_bc(etype, &xyz, (&flux, &coeff, &ext), &basis, &mut self.local)
                    .map_err(|m| SolverError::fatal(CALLER, m))?;
                ctx.default_update(&self.local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}
