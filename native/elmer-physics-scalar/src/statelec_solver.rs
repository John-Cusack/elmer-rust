//! `StatElecSolver`: the `Procedure = "StatElecSolve" "StatElecSolver"`
//! solver, driving the kernels of [`crate::statelec`] through [`SolverCtx`]
//! the way `StatElecSolver` (StatElecSolve.F90:252) drives its assembly and
//! `StatElecSolver_post` (:1010) drives the field recovery.

use std::sync::Arc;

use elmer_core::{
    ElemId, ElementBasis, FieldKind, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF,
};
use elmer_solver::sim::node_order_perm;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

use crate::common::{
    axisymmetric, fe_provider, init_dirichlet, real_at_nodes, slot_mut, unsupported,
};
use crate::post::{self, Derived, PostField, PostSpec};
use crate::potential::{Coefficient, PostConfig, PotentialPostLocal};
use crate::statelec::{
    stat_elec_boundary_local, stat_elec_boundary_local_with, stat_elec_bulk_local,
    stat_elec_bulk_local_with, StatElecAtNodes, StatElecBcAtNodes, StatElecConfig, PHYSICAL_EPS0,
};
use crate::{STAT_ELEC_PROCEDURE, STAT_ELEC_VEC_PROCEDURE};

const CALLER: &str = "StatElecSolver";

/// Registers `("StatElecSolve", "StatElecSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(
        STAT_ELEC_PROCEDURE.0,
        STAT_ELEC_PROCEDURE.1,
        stat_elec_factory,
    );
    // StatElecSolveVec is a thin front for StatElecSolve that injects the
    // physical vacuum permittivity, which this solver defaults to anyway.
    registry.register(
        STAT_ELEC_VEC_PROCEDURE.0,
        STAT_ELEC_VEC_PROCEDURE.1,
        stat_elec_factory,
    );
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`StatElecSolver`]. Fails
/// with `NotWired` until [`install_fe`](crate::install_fe) has been called.
pub fn stat_elec_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let (basis, quad) = fe_provider()?();
    Ok(Box::new(StatElecSolver::new(basis, quad)))
}

pub struct StatElecSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
    local: LocalSystem,
    post: PotentialPostLocal,
}

impl StatElecSolver {
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
        name: "nodal energy density",
        flag: "calculate nodal energy",
        field_type: 1,
        nodal: true,
    },
    PostField {
        name: "electric energy density",
        flag: "calculate electric energy",
        field_type: 1,
        nodal: false,
    },
    PostField {
        name: "nodal electric flux",
        flag: "calculate nodal flux",
        field_type: 2,
        nodal: true,
    },
    PostField {
        name: "electric flux",
        flag: "calculate electric flux",
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
    caller: "StatElecSolver_post",
    fields: &POST_FIELDS,
    total_factor: 0.5,
    total_key: "res: Electric Energy",
    potential_difference_key: Some("res: Potential Difference"),
    derived: Derived::Capacitance("RES: Effective Capacitance"),
    weight_eps: f64::MIN_POSITIVE,
    control_total: "energy control",
    control_flux: "charge control",
    scaling_key: "RES: StatElec Scaling",
    scale_on_control: false,
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
            "the legacy StatElecSolveLegacy assembly is not ported",
        ),
        (
            "calculate capacitance matrix",
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
        (
            "average within materials",
            "DG body averaging is not ported",
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
    if p.integer("variable dofs") == Some(2) {
        return Err(unsupported(
            "complex-valued potential (Variable DOFs = 2)",
            "the complex assembly needs Linear System Complex in DefaultSolve",
        ));
    }
    let model = &ctx.sim.model;
    for m in &model.materials {
        if m.get("relative permittivity im").is_some() {
            return Err(unsupported(
                "Material keyword \"Relative Permittivity Im\"",
                "the complex assembly needs Linear System Complex in DefaultSolve",
            ));
        }
    }
    for bf in &model.body_forces {
        if flag(bf, "piezo material") {
            return Err(unsupported(
                "Body Force keyword \"Piezo Material\"",
                "piezoelectric coupling is not ported",
            ));
        }
    }
    for bc in &model.boundary_conditions {
        if flag(bc, "electric infinity bc") || bc.get("infinity bc").is_some() {
            return Err(unsupported(
                "Boundary Condition \"Electric Infinity BC\"",
                "needs NormalVector(Element, Nodes, u, v, check orientation), which the element basis lacks",
            ));
        }
        for key in ["capacitance body", "constraint mode potential"] {
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

fn eps0(ctx: &SolverCtx<'_>) -> f64 {
    ctx.constants()
        .real("permittivity of vacuum")
        .unwrap_or(PHYSICAL_EPS0)
}

impl Solver for StatElecSolver {
    /// `StatElecSolver_init` (StatElecSolve.F90:129): the `Variable` default,
    /// the exported postprocessing fields on the solver's permutation, and
    /// `Time Derivative Order = 0`. The framework has no hook ahead of
    /// variable creation, so the variable and the time order are settled here.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_supported(ctx)?;
        axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();

        if ctx.params().get("variable").is_none() {
            if ctx
                .params()
                .string("element")
                .is_some_and(|e| elmer_core::normalize_key(e) != "n:1")
            {
                return Err(unsupported(
                    "Element without a Variable keyword",
                    "the default potential variable is created nodal",
                ));
            }
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
            let v = Variable::new(
                "potential",
                FieldKind::Nodal,
                1,
                Some(perm_arc.clone()),
                vec![0.0; perm.n_slots],
            );
            let id = ctx.sim.vars.add(v)?;
            let slot = slot_mut(ctx);
            slot.variable = Some(id);
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
        init_dirichlet(ctx)?;

        let slot = slot_mut(ctx);
        slot.params
            .insert("Time derivative order", Value::Integer(vec![0]));
        slot.time_order = 0;
        if let Some(id) = slot.variable {
            ctx.sim.vars.by_id_mut(id).prev_values.clear();
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        check_supported(ctx)?;
        let axi = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let cfg = StatElecConfig {
            dim,
            axisymmetric: axi,
            eps0: eps0(ctx),
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
        ctx.info(CALLER, "Solving static electric field for insulators", 4);

        if ctx.solver_dofs().is_some() && post::post_active(ctx.params(), &POST) {
            return Err(unsupported(
                "field and energy postprocessing on p elements",
                "the potential post-processing integrates with the nodal basis only",
            ));
        }
        let mut dofs = Vec::new();
        for _ in 1..=maxiter {
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let mat = bulk_inputs(ctx, e, xyz.len())?;
                match ctx.solver_dofs() {
                    Some(d) => {
                        // LocalMatrix(Element, n, nd+nb, nb), then CondensateP(nd, nb).
                        let p = p_element(d, e)?;
                        stat_elec_bulk_local_with(
                            &*p.basis,
                            &p.rule,
                            el.etype,
                            &xyz,
                            &mat,
                            &cfg,
                            p.nd + p.nb,
                            &mut self.local,
                        )?;
                        condense(&mut self.local, p.nd, p.nb)?;
                    }
                    None => stat_elec_bulk_local(
                        &*self.basis,
                        &*self.quad,
                        el.etype,
                        &xyz,
                        &mat,
                        &cfg,
                        &mut self.local,
                    )?,
                }
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let Some(bc) = bc_inputs(ctx, e)? else {
                    continue;
                };
                let el = ctx.mesh.element(e);
                match ctx.solver_dofs() {
                    Some(d) => {
                        // A p boundary element lists its nodes in the parent's
                        // edge order (AssignLocalNumber); its values follow it.
                        let p = p_element(d, e)?;
                        let xyz: Vec<[f64; 3]> =
                            p.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                        let bc = bc.reordered(el.nodes, &p.nodes);
                        stat_elec_boundary_local_with(
                            &*p.basis,
                            &p.rule,
                            el.etype,
                            &xyz,
                            &bc,
                            &cfg,
                            dofs.len(),
                            &mut self.local,
                        )?;
                    }
                    None => {
                        let xyz: Vec<[f64; 3]> =
                            el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                        stat_elec_boundary_local(
                            &*self.basis,
                            &*self.quad,
                            el.etype,
                            &xyz,
                            &bc,
                            &cfg,
                            &mut self.local,
                        )?;
                    }
                }
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
            let coef_of =
                |c: &SolverCtx<'_>, e: ElemId, n: usize| -> Result<Coefficient, SolverError> {
                    Ok(Coefficient::Scalar(
                        bulk_inputs(c, e, n)?.relative_permittivity,
                    ))
                };
            post::postprocess(
                ctx,
                &*self.basis,
                &*self.quad,
                &mut self.post,
                &POST,
                &pcfg,
                cfg.eps0,
                &coef_of,
            )?;
        }
        Ok(())
    }
}

/// An element's basis, rule, node order, dof count and local bubble count
/// from the solver's p dofs.
struct PElement {
    basis: Box<dyn elmer_core::ElementBasis + Send + Sync>,
    rule: Vec<elmer_core::QuadPoint>,
    nodes: Vec<u32>,
    nd: usize,
    nb: usize,
}

fn p_element(d: &dyn elmer_solver::SolverDofs, e: ElemId) -> Result<PElement, SolverError> {
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    Ok(PElement {
        basis: d
            .element_basis(e)
            .map_err(|m| SolverError::fatal("ElementInfo", m))?,
        rule: d
            .element_rule(e)
            .map_err(|m| SolverError::fatal("GaussPoints", m))?,
        nodes: d.element_nodes(e),
        nd: slots.len(),
        nb: d.local_bubbles(e),
    })
}

/// `CondensateP(nd, nb, STIFF, FORCE)`: eliminates the last `nb` of the
/// `nd + nb` local dofs, leaving an `nd`-square system.
fn condense(local: &mut LocalSystem, nd: usize, nb: usize) -> Result<(), SolverError> {
    if nb == 0 {
        return Ok(());
    }
    let m = nd + nb;
    let k = &local.stiff;
    let f = &local.force;
    // Kbb^-1 by Gauss-Jordan with partial pivoting.
    let mut a: Vec<f64> = (0..nb * nb)
        .map(|ij| k[(nd + ij / nb) * m + nd + ij % nb])
        .collect();
    let mut inv = vec![0.0; nb * nb];
    for i in 0..nb {
        inv[i * nb + i] = 1.0;
    }
    for c in 0..nb {
        let piv = (c..nb)
            .max_by(|&x, &y| a[x * nb + c].abs().total_cmp(&a[y * nb + c].abs()))
            .unwrap_or(c);
        if a[piv * nb + c] == 0.0 {
            return Err(SolverError::fatal(
                CALLER,
                "CondensateP: singular bubble block",
            ));
        }
        for j in 0..nb {
            a.swap(c * nb + j, piv * nb + j);
            inv.swap(c * nb + j, piv * nb + j);
        }
        let dg = a[c * nb + c];
        for j in 0..nb {
            a[c * nb + j] /= dg;
            inv[c * nb + j] /= dg;
        }
        for i in 0..nb {
            if i != c {
                let fac = a[i * nb + c];
                if fac != 0.0 {
                    for j in 0..nb {
                        a[i * nb + j] -= fac * a[c * nb + j];
                        inv[i * nb + j] -= fac * inv[c * nb + j];
                    }
                }
            }
        }
    }
    let mut ks = vec![0.0; nd * nd];
    let mut fs = vec![0.0; nd];
    for i in 0..nd {
        // Klb Kbb^-1, row i.
        let g: Vec<f64> = (0..nb)
            .map(|j| (0..nb).map(|l| k[i * m + nd + l] * inv[l * nb + j]).sum())
            .collect();
        fs[i] = f[i] - (0..nb).map(|j| g[j] * f[nd + j]).sum::<f64>();
        for c in 0..nd {
            ks[i * nd + c] =
                k[i * m + c] - (0..nb).map(|j| g[j] * k[(nd + j) * m + c]).sum::<f64>();
        }
    }
    local.n = nd;
    local.stiff = ks;
    local.force = fs;
    Ok(())
}

fn bulk_inputs(ctx: &SolverCtx<'_>, e: ElemId, n: usize) -> Result<StatElecAtNodes, SolverError> {
    let bp = ctx
        .body_params(e)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
    Ok(StatElecAtNodes {
        relative_permittivity: real_at_nodes(ctx, bp.material, "Relative Permittivity", e)?
            .unwrap_or_else(|| vec![0.0; n]),
        charge_density: real_at_nodes(ctx, bp.body_force, "Charge Density", e)?,
    })
}

/// The boundary inputs of element `e`, or `None` when its condition sets
/// none of the keywords `LocalMatrixBC` reads.
fn bc_inputs(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Option<StatElecBcAtNodes>, SolverError> {
    let p = ctx.bc_params(e);
    let mut electric_flux = real_at_nodes(ctx, p, "Electric Flux", e)?;
    if electric_flux.is_none() {
        electric_flux = real_at_nodes(ctx, p, "Surface Charge Density", e)?;
    }
    let bc = StatElecBcAtNodes {
        electric_flux,
        layer_relative_permittivity: real_at_nodes(ctx, p, "Layer Relative Permittivity", e)?,
        layer_thickness: real_at_nodes(ctx, p, "Layer Thickness", e)?,
        electrode_potential: real_at_nodes(ctx, p, "Electrode Potential", e)?,
        layer_charge_density: real_at_nodes(ctx, p, "Layer Charge Density", e)?,
    };
    Ok((bc.electric_flux.is_some() || bc.layer_relative_permittivity.is_some()).then_some(bc))
}
