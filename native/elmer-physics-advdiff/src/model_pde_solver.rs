//! `AdvDiffSolver` of ModelPDE.F90: the `Procedure = "ModelPDE"
//! "AdvDiffSolver"` solver.

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::fe::Fe;
use crate::model_pde::{
    model_pde_boundary_local, model_pde_bulk_local, ModelPdeBcNodal, ModelPdeNodal,
};
use crate::preflight::{check_model_pde_supported, ModelPdeSections, MODEL_PDE_PROCEDURE};

/// Registers `("ModelPDE", "AdvDiffSolver")`.
pub fn register_model_pde(registry: &mut SolverRegistry) {
    registry.register(
        MODEL_PDE_PROCEDURE.0,
        MODEL_PDE_PROCEDURE.1,
        model_pde_factory,
    );
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`ModelPdeSolver`].
pub fn model_pde_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ModelPdeSolver {
        local: LocalSystem::zeros(0),
    }))
}

pub struct ModelPdeSolver {
    local: LocalSystem,
}

fn nodal_or_zero(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    key: &str,
    e: ElemId,
    n: usize,
) -> Result<Vec<f64>, SolverError> {
    match p {
        Some(p) => Ok(ctx
            .real_at_nodes(p, key, e)?
            .unwrap_or_else(|| vec![0.0; n])),
        None => Ok(vec![0.0; n]),
    }
}

impl Solver for ModelPdeSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|f| f.dofs != 1) {
            return Err(SolverError::Unsupported(
                "ModelPDE solves a scalar variable (Variable DOFs = 1)".into(),
            ));
        }
        check_model_pde_supported(&ModelPdeSections {
            solver: ctx.params(),
        })?;
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        let dim = ctx.mesh.dim();
        let basis = LagrangeBasis::new(dim);
        let quad = ElmerQuadrature;
        let fe = Fe::new(&basis, &quad);
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();

        for _ in 0..maxiter {
            ctx.default_initialize()?;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let nodal = {
                    let bp = ctx.body_params(e).ok_or_else(|| {
                        SolverError::fatal("AdvDiffSolver", format!("element {} has no body", e.0))
                    })?;
                    let (m, bf) = (bp.material, bp.body_force);
                    let mut velocity = vec![[0.0; 3]; n];
                    for i in 0..dim {
                        let v =
                            nodal_or_zero(ctx, m, &format!("convection velocity {}", i + 1), e, n)?;
                        for (vel, x) in velocity.iter_mut().zip(v) {
                            vel[i] = x;
                        }
                    }
                    ModelPdeNodal {
                        diffusion: nodal_or_zero(ctx, m, "diffusion coefficient", e, n)?,
                        convection: nodal_or_zero(ctx, m, "convection coefficient", e, n)?,
                        reaction: nodal_or_zero(ctx, m, "reaction coefficient", e, n)?,
                        time: nodal_or_zero(ctx, m, "time derivative coefficient", e, n)?,
                        velocity,
                        source: nodal_or_zero(ctx, bf, "field source", e, n)?,
                    }
                };
                model_pde_bulk_local(&fe, etype, &xyz, &nodal, dim, None, &mut self.local)?;
                if transient {
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let p = ctx.bc_params(e);
                let bc = ModelPdeBcNodal {
                    flux: nodal_or_zero(ctx, p, "field flux", e, n)?,
                    robin: nodal_or_zero(ctx, p, "robin coefficient", e, n)?,
                    external: nodal_or_zero(ctx, p, "external field", e, n)?,
                };
                model_pde_boundary_local(&fe, etype, &xyz, &bc, &mut self.local)?;
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}
