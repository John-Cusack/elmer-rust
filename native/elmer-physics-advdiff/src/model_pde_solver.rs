//! `AdvDiffSolver` of ModelPDE.F90: the `Procedure = "ModelPDE"
//! "AdvDiffSolver"` solver.

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_solver::SolverDofs;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::fe::{condensate, Fe};
use crate::model_pde::{
    model_pde_boundary_local, model_pde_boundary_local_with, model_pde_bulk_local,
    model_pde_bulk_local_with, ModelPdeBcNodal, ModelPdeNodal,
};
use crate::preflight::{
    check_model_pde_supported, check_model_pde_supported_p, ModelPdeSections, MODEL_PDE_PROCEDURE,
};

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

/// `values` (at `mesh_nodes`, the mesh's node order) at `nodes` instead:
/// p boundary elements list their nodes in the parent edge or face order.
fn reorder(values: Vec<f64>, mesh_nodes: &[u32], nodes: &[u32]) -> Vec<f64> {
    if mesh_nodes == nodes {
        return values;
    }
    nodes
        .iter()
        .map(|n| {
            values[mesh_nodes
                .iter()
                .position(|m| m == n)
                .expect("same node set")]
        })
        .collect()
}

fn dofs_error(what: &str, msg: String) -> SolverError {
    SolverError::fatal("AdvDiffSolver", format!("{what}: {msg}"))
}

/// The element's basis, rule, node order, dof count and local bubble count
/// from the solver's p dofs.
struct PElement {
    basis: Box<dyn elmer_core::ElementBasis + Send + Sync>,
    rule: Vec<elmer_core::QuadPoint>,
    nodes: Vec<u32>,
    nd: usize,
    nb: usize,
}

fn p_element(d: &dyn SolverDofs, e: ElemId) -> Result<PElement, SolverError> {
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    Ok(PElement {
        basis: d
            .element_basis(e)
            .map_err(|m| dofs_error("ElementInfo", m))?,
        rule: d
            .element_rule(e)
            .map_err(|m| dofs_error("GaussPoints", m))?,
        nodes: d.element_nodes(e),
        nd: slots.len(),
        nb: d.local_bubbles(e),
    })
}

/// Leading `n x n` block and `n` entries of an `m`-square local system.
fn shrink(local: &mut LocalSystem, n: usize) {
    let m = local.n;
    if m == n {
        return;
    }
    let mut k = vec![0.0; n * n];
    for i in 0..n {
        k[i * n..(i + 1) * n].copy_from_slice(&local.stiff[i * m..i * m + n]);
    }
    local.stiff = k;
    local.force.truncate(n);
    local.mass = None;
    local.n = n;
}

impl Solver for ModelPdeSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|f| f.dofs != 1) {
            return Err(SolverError::Unsupported(
                "ModelPDE solves a scalar variable (Variable DOFs = 1)".into(),
            ));
        }
        let sections = ModelPdeSections {
            solver: ctx.params(),
        };
        if ctx.solver_dofs().is_some() {
            check_model_pde_supported_p(&sections)?;
        } else {
            check_model_pde_supported(&sections)?;
        }
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
                let pe = match ctx.solver_dofs() {
                    Some(d) => Some(p_element(d, e)?),
                    None => None,
                };
                let nodes: Vec<u32> = pe
                    .as_ref()
                    .map_or_else(|| el.nodes.to_vec(), |p| p.nodes.clone());
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
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
                match &pe {
                    None => {
                        model_pde_bulk_local(&fe, etype, &xyz, &nodal, dim, None, &mut self.local)?
                    }
                    Some(p) => {
                        // LocalMatrix(Element, n, nd+nb), then CondensateP(nd, nb).
                        let m = p.nd + p.nb;
                        model_pde_bulk_local_with(
                            &*p.basis,
                            &p.rule,
                            etype,
                            &xyz,
                            &nodal,
                            dim,
                            m,
                            &mut self.local,
                        )?;
                        if p.nb > 0 {
                            if transient {
                                return Err(SolverError::Unsupported(
                                    "ModelPDE: transient p elements with local bubbles".into(),
                                ));
                            }
                            let l = &mut self.local;
                            condensate(p.nd, p.nb, &mut l.stiff, &mut l.force, None)?;
                            shrink(l, p.nd);
                        }
                    }
                }
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
                let pe = match ctx.solver_dofs() {
                    Some(d) => Some(p_element(d, e)?),
                    None => None,
                };
                let nodes: Vec<u32> = pe
                    .as_ref()
                    .map_or_else(|| el.nodes.to_vec(), |p| p.nodes.clone());
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let p = ctx.bc_params(e);
                let at = |key: &str| -> Result<Vec<f64>, SolverError> {
                    Ok(reorder(nodal_or_zero(ctx, p, key, e, n)?, el.nodes, &nodes))
                };
                let bc = ModelPdeBcNodal {
                    flux: at("field flux")?,
                    robin: at("robin coefficient")?,
                    external: at("external field")?,
                };
                match &pe {
                    None => model_pde_boundary_local(&fe, etype, &xyz, &bc, &mut self.local)?,
                    Some(pe) => model_pde_boundary_local_with(
                        &*pe.basis,
                        &pe.rule,
                        etype,
                        &xyz,
                        &bc,
                        dofs.len(),
                        &mut self.local,
                    )?,
                }
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
