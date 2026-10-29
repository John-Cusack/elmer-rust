//! `ModelPDEsplitted/ModelPDEsplitted.f90`: ModelPDE written as two step
//! routines, `AdvDiffSolver_bulk` and `AdvDiffSolver_boundary`, for
//! `Procedure = "ModelPDEsplitted" "AdvDiffSolver"`.
//!
//! The library has no `AdvDiffSolver` itself. Because `<proc>_bulk` exists,
//! Elmer runs the section in steps (`SOLVER_MODE_STEPS`,
//! MainUtils.F90:1402) through `ExecSolverInSteps`: per nonlinear
//! iteration `DefaultInitialize`, the bulk routine,
//! `DefaultFinishBulkAssembly`, the boundary routine,
//! `DefaultFinishBoundaryAssembly`, `DefaultFinishAssembly`,
//! `DefaultDirichletBCs` and `DefaultSolve`, until `NonlinConverged`.
//! [`ModelPdeSplitted`] is that driver around the two routines.
//!
//! The routines read nodal property values (`GetReal`) and interpolate them
//! with the nodal basis functions, assemble all `nd + nb` dofs of the
//! element and condense the local bubbles (`LCondensate`). On p elements
//! (the cases' `Element = "p:2"`) basis, Gauss rule, node order and dof
//! counts come from the solver's dof layout (`elmer-pdofs`, through
//! `SolverCtx::solver_dofs`); without one the element is nodal Lagrange.

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, QuadPoint, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{coordinate_system, element, nodal};
use crate::linalg::condensate;

const CALLER: &str = "AdvDiffSolver";

/// `SUM(Basis(1:n) * values(1:n))`.
fn interp(basis: &[f64], values: &[f64]) -> f64 {
    basis.iter().zip(values).fold(0.0, |s, (b, v)| s + b * v)
}

/// The bulk properties of one element at its nodes.
pub struct BulkNodal {
    pub load: Vec<f64>,
    pub diff: Vec<f64>,
    pub react: Vec<f64>,
    pub conv: Vec<f64>,
    pub time: Vec<f64>,
    pub velo: [Vec<f64>; 3],
}

/// `AdvDiffSolver_bulk`'s `LocalMatrix` over `m = nd + nb` dofs: stiffness,
/// mass and force.
#[allow(clippy::too_many_arguments)]
pub fn bulk_local(
    etype: ElementType,
    xyz: &[[f64; 3]],
    m: usize,
    props: &BulkNodal,
    dim: usize,
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    local: &mut LocalSystem,
) -> Result<(), String> {
    let n = xyz.len();
    *local = LocalSystem::zeros(m);
    let mut mass = vec![0.0; m * m];
    let mut ev = BasisEval::default();
    for qp in rule {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < m {
            return Err(format!("ElementInfo gave {} basis functions for {m} dofs", ev.basis.len()));
        }
        let (b, g) = (&ev.basis[..m], &ev.dbasis_dx[..m]);
        let bn = &b[..n];
        let load = interp(bn, &props.load);
        let rho = interp(bn, &props.time);
        // a = MATMUL(Velo(:,1:n), Basis(1:n))
        let a: [f64; 3] = std::array::from_fn(|i| props.velo[i].iter().zip(bn).fold(0.0, |s, (v, bp)| s + v * bp));
        let d = interp(bn, &props.diff);
        let c = interp(bn, &props.conv);
        let r = interp(bn, &props.react);
        let w = qp.weight * ev.det_j;
        for p in 0..m {
            for q in 0..m {
                let mm = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                local.stiff[p * m + q] += w * d * mm;
            }
        }
        for p in 0..m {
            for q in 0..m {
                let adv = (0..dim).fold(0.0, |s, k| s + a[k] * g[q][k]);
                local.stiff[p * m + q] += w * c * adv * b[p];
                local.stiff[p * m + q] += w * r * b[q] * b[p];
                mass[p * m + q] += w * rho * b[q] * b[p];
            }
        }
        for (f, bp) in local.force.iter_mut().zip(b) {
            *f += w * load * bp;
        }
    }
    local.mass = Some(mass);
    Ok(())
}

/// `AdvDiffSolver_boundary`'s `LocalMatrixBC` over `nd` dofs.
#[allow(clippy::too_many_arguments)]
pub fn boundary_local(
    etype: ElementType,
    xyz: &[[f64; 3]],
    nd: usize,
    (flux, coeff, ext): (&[f64], &[f64], &[f64]),
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    local: &mut LocalSystem,
) -> Result<(), String> {
    let n = xyz.len();
    *local = LocalSystem::zeros(nd);
    let mut ev = BasisEval::default();
    for qp in rule {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < nd {
            return Err(format!("ElementInfo gave {} basis functions for {nd} dofs", ev.basis.len()));
        }
        let b = &ev.basis[..nd];
        let w = qp.weight * ev.det_j;
        let (f, c, x) = (interp(&b[..n], flux), interp(&b[..n], coeff), interp(&b[..n], ext));
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

pub struct ModelPdeSplitted {
    local: LocalSystem,
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for the stepped
/// `AdvDiffSolver` (`_bulk` + `_boundary`).
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ModelPdeSplitted { local: LocalSystem::zeros(0) }))
}

impl ModelPdeSplitted {
    /// `AdvDiffSolver_bulk`.
    fn bulk(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let dim = coordinate_system(ctx)?.dim;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        for e in active {
            let el = element(ctx, e)?;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let (m, bf) = ctx.body_params(e).map_or((None, None), |b| (b.material, b.body_force));
            let o = &el.nodes;
            let props = BulkNodal {
                load: nodal(ctx, bf, "field source", e, o)?,
                diff: nodal(ctx, m, "diffusion coefficient", e, o)?,
                react: nodal(ctx, m, "reaction coefficient", e, o)?,
                conv: nodal(ctx, m, "convection coefficient", e, o)?,
                time: nodal(ctx, m, "time derivative coefficient", e, o)?,
                velo: [
                    if dim >= 1 { nodal(ctx, m, "convection velocity 1", e, o)? } else { vec![0.0; o.len()] },
                    if dim >= 2 { nodal(ctx, m, "convection velocity 2", e, o)? } else { vec![0.0; o.len()] },
                    if dim >= 3 { nodal(ctx, m, "convection velocity 3", e, o)? } else { vec![0.0; o.len()] },
                ],
            };
            let total = el.nd + el.nb;
            bulk_local(el.etype, &xyz, total, &props, dim, el.basis.as_ref(), &el.rule, &mut self.local)
                .map_err(|m| SolverError::fatal(CALLER, m))?;
            if transient {
                if el.nb > 0 {
                    return Err(SolverError::Unsupported(
                        "ModelPDEsplitted: transient run with local bubbles (Default1stOrderTime before LCondensate)".into(),
                    ));
                }
                ctx.default_first_order_time(&mut self.local, e, dt)?;
            }
            if el.nb > 0 {
                let l = &mut self.local;
                condensate(el.nd, el.nb, &mut l.stiff, &mut l.force).map_err(|m| SolverError::fatal(CALLER, m))?;
                l.n = el.nd;
                l.mass = None;
            }
            ctx.default_update(&self.local, e)?;
        }
        Ok(())
    }

    /// `AdvDiffSolver_boundary`.
    fn boundary(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();
        for e in boundary {
            // ActiveBoundaryElement()
            ctx.element_dofs(e, &mut dofs);
            if dofs.contains(&NO_DOF) {
                continue;
            }
            let el = element(ctx, e)?;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let bc = ctx.bc_params(e);
            let o = &el.nodes;
            let flux = nodal(ctx, bc, "field flux", e, o)?;
            let coeff = nodal(ctx, bc, "robin coefficient", e, o)?;
            let ext = nodal(ctx, bc, "external field", e, o)?;
            boundary_local(el.etype, &xyz, el.nd, (&flux, &coeff, &ext), el.basis.as_ref(), &el.rule, &mut self.local)
                .map_err(|m| SolverError::fatal(CALLER, m))?;
            ctx.default_update(&self.local, e)?;
        }
        Ok(())
    }
}

impl Solver for ModelPdeSplitted {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        crate::hooks::check_element_definitions(ctx)?;
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::Unsupported("ModelPDEsplitted AdvDiffSolver: a scalar Variable".into()));
        }
        let element = ctx.params().string("element");
        if ctx.solver_dofs().is_none() && element.is_some_and(|e| elmer_core::normalize_key(e) != "n:1") {
            return Err(SolverError::NotWired(format!(
                "ModelPDEsplitted AdvDiffSolver: Element = {} without a dof layout (elmer-pdofs in Backends::dofs)",
                element.unwrap_or_default()
            )));
        }
        Ok(())
    }

    /// `ExecSolverInSteps`.
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let maxiter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1);
        for _ in 0..maxiter {
            ctx.default_initialize()?;
            self.bulk(ctx, dt, transient)?;
            ctx.default_finish_bulk_assembly()?;
            self.boundary(ctx)?;
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
