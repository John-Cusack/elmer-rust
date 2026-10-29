//! `StatCurrentSolver`: the `Procedure = "StatCurrentSolve"
//! "StatCurrentSolver"` solver, driving the kernels of [`crate::statcurrent`]
//! through [`SolverCtx`] the way `StatCurrentSolver` (StatCurrentSolve.F90:211)
//! and `StatCurrentSolver_post` (:830) do.

use elmer_core::{
    ElemId, ElementBasis, FieldKind, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF,
};
use elmer_solver::sim::node_order_perm;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};
use std::sync::Arc;

use crate::common::{
    axisymmetric, coefficient_at_nodes, fe_provider, init_dirichlet, real_at_nodes, slot_mut,
    unsupported,
};
use crate::post::{self, Derived, PostField, PostSpec};
use crate::potential::{Coefficient, PostConfig, PotentialPostLocal};
use crate::statcurrent::{
    stat_current_boundary_local, stat_current_bulk_local, StatCurrentAtNodes, StatCurrentBcAtNodes,
    StatCurrentConfig,
};
use crate::statelec::PHYSICAL_EPS0;
use crate::STAT_CURRENT_PROCEDURE;

const CALLER: &str = "StatCurrentSolver";

/// Registers `("StatCurrentSolve", "StatCurrentSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(
        STAT_CURRENT_PROCEDURE.0,
        STAT_CURRENT_PROCEDURE.1,
        stat_current_factory,
    );
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`StatCurrentSolver`].
/// Fails with `NotWired` until [`install_fe`](crate::install_fe) has been called.
pub fn stat_current_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let (basis, quad) = fe_provider()?();
    Ok(Box::new(StatCurrentSolver::new(basis, quad)))
}

pub struct StatCurrentSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
    local: LocalSystem,
    post: PotentialPostLocal,
}

impl StatCurrentSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: LocalSystem::zeros(0),
            post: PotentialPostLocal::default(),
        }
    }
}

const POST_FIELDS: [PostField; 5] = [
    PostField {
        name: "nodal joule heating",
        flag: "calculate nodal heating",
        field_type: 1,
        nodal: true,
    },
    PostField {
        name: "joule heating",
        flag: "calculate joule heating",
        field_type: 1,
        nodal: false,
    },
    PostField {
        name: "nodal current",
        flag: "calculate nodal current",
        field_type: 2,
        nodal: true,
    },
    PostField {
        name: "volume current",
        flag: "calculate volume current",
        field_type: 2,
        nodal: false,
    },
    PostField {
        name: "electric field",
        flag: "calculate electric field",
        field_type: 3,
        nodal: false,
    },
];

const POST: PostSpec = PostSpec {
    caller: "StatCurrentSolver_post",
    fields: &POST_FIELDS,
    total_factor: 1.0,
    total_key: "RES: Total Joule Heating",
    potential_difference_key: None,
    derived: Derived::Resistance("RES: Effective Resistance"),
    weight_eps: f64::EPSILON,
    control_total: "power control",
    control_flux: "current control",
    scaling_key: "RES: CurrentSolver Scaling",
    scale_on_control: true,
};

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// Refuses what this port does not implement, instead of dropping it.
fn check_supported(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for (key, reason) in [
        (
            "legacy assembly",
            "the legacy StatCurrentSolveLegacy assembly is not ported",
        ),
        (
            "harmonic mode",
            "the complex harmonic assembly is not ported",
        ),
        (
            "calculate conductivity matrix",
            "Constraint Modes analysis is not ported",
        ),
        (
            "constraint modes analysis",
            "Constraint Modes analysis is not ported",
        ),
        ("adaptive mesh refinement", "mesh refinement is not ported"),
        (
            "calculate elemental fields",
            "elemental (-dg) fields are not ported",
        ),
    ] {
        if flag(p, key) {
            return Err(unsupported(format!("Solver keyword \"{key}\""), reason));
        }
    }
    if p.integer("relative integration order")
        .is_some_and(|n| n != 0)
    {
        return Err(unsupported(
            "Solver keyword \"Relative Integration Order\"",
            "the quadrature takes a point count, not an order",
        ));
    }
    for bc in &ctx.sim.model.boundary_conditions {
        if bc.get("farfield potential").is_some() {
            return Err(unsupported(
                "Boundary Condition \"Farfield Potential\"",
                "needs NormalVector(Element, Nodes, u, v, check orientation), which the element basis lacks",
            ));
        }
        for key in ["conductivity body", "constraint mode potential"] {
            if bc.get(key).is_some() {
                return Err(unsupported(
                    format!("Boundary Condition keyword \"{key}\""),
                    "Constraint Modes analysis is not ported",
                ));
            }
        }
    }
    Ok(())
}

fn bulk_coefficient(ctx: &SolverCtx<'_>, e: ElemId, n: usize) -> Result<Coefficient, SolverError> {
    let bp = ctx
        .body_params(e)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
    coefficient_at_nodes(
        ctx,
        bp.material,
        "Electric Conductivity",
        e,
        n,
        ctx.mesh.dim(),
    )?
    .ok_or_else(|| {
        SolverError::fatal(
            "ListGetElementReal",
            "Failed to find \"Electric Conductivity\"",
        )
    })
}

impl Solver for StatCurrentSolver {
    /// `StatCurrentSolver_init` (StatCurrentSolve.F90:95): the `Variable`
    /// default, the exported postprocessing fields on the solver's
    /// permutation and `Time Derivative Order = 1`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_supported(ctx)?;
        axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();

        if ctx.params().get("variable").is_none() {
            let mesh = ctx.mesh;
            let id = ctx.solver_id();
            let sim = &mut *ctx.sim;
            let slot = sim
                .slots
                .iter()
                .find(|s| s.id == id)
                .expect("the running solver has a slot");
            let perm = match &sim.backends.assembly {
                Some(b) => b.build_perm_for(mesh, &slot.active, &slot.params)?,
                None => node_order_perm(mesh, &slot.active),
            };
            let perm_arc: Arc<[u32]> = perm.perm.into();
            let mut v = Variable::new(
                "potential",
                FieldKind::Nodal,
                1,
                Some(perm_arc.clone()),
                vec![0.0; perm.n_slots],
            );
            if sim.is_transient() {
                v.prev_values = vec![vec![0.0; v.values.len()]; 1];
            }
            let vid = sim.vars.add(v)?;
            let slot = slot_mut(ctx);
            slot.variable = Some(vid);
            slot.perm = Some(perm_arc);
        }
        {
            let v = ctx
                .primary()
                .ok_or_else(|| SolverError::fatal(CALLER, "solver has no potential variable"))?;
            if v.dofs != 1 {
                return Err(unsupported(
                    format!("variable \"{}\" with {} dofs", v.name, v.dofs),
                    "only the real single-potential path is ported",
                ));
            }
        }
        post::create_fields(ctx, &POST, dim)?;
        slot_mut(ctx)
            .params
            .insert("Time derivative order", Value::Integer(vec![1]));
        init_dirichlet(ctx)?;
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        check_supported(ctx)?;
        let axi = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let eps0 = ctx
            .constants()
            .real("permittivity of vacuum")
            .unwrap_or(PHYSICAL_EPS0);
        let cfg = StatCurrentConfig {
            dim,
            axisymmetric: axi,
            transient,
            eps0,
            quad_points: None,
        };
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();

        ctx.info(
            CALLER,
            "------------------------------------------------",
            4,
        );
        ctx.info(CALLER, "Solving static current conduction solver", 4);

        let mut dofs = Vec::new();
        for _ in 1..=maxiter {
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = {
                    let bp = ctx.body_params(e).ok_or_else(|| {
                        SolverError::fatal(CALLER, format!("element {} has no body", e.0))
                    })?;
                    StatCurrentAtNodes {
                        conductivity: bulk_coefficient(ctx, e, n)?,
                        relative_permittivity: if transient {
                            real_at_nodes(ctx, bp.material, "Relative Permittivity", e)?
                        } else {
                            None
                        },
                        current_source: real_at_nodes(ctx, bp.body_force, "Current Source", e)?,
                    }
                };
                stat_current_bulk_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &cfg,
                    &mut self.local,
                )?;
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
                let p = ctx.bc_params(e);
                let bc = StatCurrentBcAtNodes {
                    current_density: real_at_nodes(ctx, p, "Current Density", e)?,
                    external_conductivity: real_at_nodes(ctx, p, "External Conductivity", e)?,
                    external_potential: real_at_nodes(ctx, p, "External Potential", e)?,
                };
                if bc.current_density.is_none() && bc.external_conductivity.is_none() {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                stat_current_boundary_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &bc,
                    &cfg,
                    &mut self.local,
                )?;
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }

        if post::post_active(ctx.params(), &POST) {
            let pcfg = PostConfig {
                dim,
                axisymmetric: axi,
                quad_points: None,
            };
            post::postprocess(
                ctx,
                &*self.basis,
                &*self.quad,
                &mut self.post,
                &POST,
                &pcfg,
                1.0,
                &bulk_coefficient,
            )?;
        }
        Ok(())
    }
}
