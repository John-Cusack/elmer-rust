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

use crate::common::{axisymmetric, fe_provider, real_at_nodes, slot_mut, unsupported};
use crate::statelec::{
    stat_elec_bulk_local, stat_elec_boundary_local, stat_elec_post_local, DenseLu,
    StatElecAtNodes, StatElecBcAtNodes, StatElecConfig, StatElecPostLocal, PHYSICAL_EPS0,
};
use crate::STAT_ELEC_PROCEDURE;

const CALLER: &str = "StatElecSolver";

/// Registers `("StatElecSolve", "StatElecSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(STAT_ELEC_PROCEDURE.0, STAT_ELEC_PROCEDURE.1, stat_elec_factory);
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
    post: StatElecPostLocal,
}

impl StatElecSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: LocalSystem::zeros(0),
            post: StatElecPostLocal::default(),
        }
    }
}

/// A postprocessing field: `PostVars` in StatElecSolver_post. Nodal fields
/// (`Nodal Energy Density`, `Nodal Electric Flux`) come before the
/// distributed ones of the same type, as the Fortran requires.
struct PostField {
    name: &'static str,
    flag: &'static str,
    /// 1 energy density (one component), 2 flux, 3 field (`dim` components).
    field_type: u8,
    nodal: bool,
}

const POST_FIELDS: [PostField; 5] = [
    PostField { name: "nodal energy density", flag: "calculate nodal energy", field_type: 1, nodal: true },
    PostField { name: "electric energy density", flag: "calculate electric energy", field_type: 1, nodal: false },
    PostField { name: "nodal electric flux", flag: "calculate nodal flux", field_type: 2, nodal: true },
    PostField { name: "electric flux", flag: "calculate electric flux", field_type: 2, nodal: false },
    PostField { name: "electric field", flag: "calculate electric field", field_type: 3, nodal: false },
];

impl PostField {
    fn dofs(&self, dim: usize) -> usize {
        if self.field_type == 1 {
            1
        } else {
            dim
        }
    }

    fn components(&self, dim: usize) -> Vec<String> {
        if self.dofs(dim) == 1 {
            Vec::new()
        } else {
            (1..=dim).map(|k| format!("{} {k}", self.name)).collect()
        }
    }
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// Refuses what this port does not implement, instead of dropping it.
fn check_supported(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for (key, reason) in [
        ("legacy assembly", "the legacy StatElecSolveLegacy assembly is not ported"),
        ("calculate capacitance matrix", "Constraint Modes analysis is not ported"),
        ("constraint modes analysis", "Constraint Modes analysis is not ported"),
        ("adaptive mesh refinement", "mesh refinement is not ported"),
        ("calculate elemental fields", "elemental (-dg) fields are not ported"),
        ("average within materials", "DG body averaging is not ported"),
    ] {
        if flag(p, key) {
            return Err(unsupported(format!("Solver keyword \"{key}\""), reason));
        }
    }
    if p.integer("relative integration order").is_some_and(|n| n != 0) {
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
            return Err(unsupported("Body Force keyword \"Piezo Material\"", "piezoelectric coupling is not ported"));
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
    ctx.constants().real("permittivity of vacuum").unwrap_or(PHYSICAL_EPS0)
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
            let mesh = ctx.mesh;
            let id = ctx.solver_id();
            let sim = &mut *ctx.sim;
            let slot = sim.slots.iter().find(|s| s.id == id).expect("the running solver has a slot");
            let perm = match &sim.backends.assembly {
                Some(b) => b.build_perm_for(mesh, &slot.active, &slot.params)?,
                None => node_order_perm(mesh, &slot.active),
            };
            let perm_arc: Arc<[u32]> = perm.perm.into();
            let v = Variable::new("potential", FieldKind::Nodal, 1, Some(perm_arc.clone()), vec![0.0; perm.n_slots]);
            let id = ctx.sim.vars.add(v)?;
            let slot = slot_mut(ctx);
            slot.variable = Some(id);
            slot.perm = Some(perm_arc);
        }

        let perm = slot_mut(ctx).perm.clone();
        let n_slots = {
            let v = ctx
                .primary()
                .ok_or_else(|| SolverError::fatal(CALLER, "solver has no potential variable"))?;
            if v.dofs != 1 {
                return Err(unsupported(
                    format!("variable \"{}\" with {} dofs", v.name, v.dofs),
                    "only the real single-potential path is ported",
                ));
            }
            v.values.len()
        };
        let wanted: Vec<&PostField> = POST_FIELDS.iter().filter(|f| flag(ctx.params(), f.flag)).collect();
        for f in wanted {
            if ctx.sim.vars.id(f.name).is_some() {
                continue;
            }
            let dofs = f.dofs(dim);
            let mut v = Variable::new(f.name, FieldKind::Nodal, dofs, perm.clone(), vec![0.0; n_slots * dofs]);
            v.secondary = true;
            ctx.sim.vars.add_with_components(v, &f.components(dim))?;
        }

        let slot = slot_mut(ctx);
        slot.params.insert("Time derivative order", Value::Integer(vec![0]));
        slot.time_order = 0;
        if let Some(id) = slot.variable {
            ctx.sim.vars.by_id_mut(id).prev_values.clear();
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
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

        ctx.info(CALLER, "------------------------------------------------", 4);
        ctx.info(CALLER, "Solving static electric field for insulators", 4);

        let mut dofs = Vec::new();
        for _ in 1..=maxiter {
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let mat = self.bulk_inputs(ctx, e, xyz.len())?;
                stat_elec_bulk_local(&*self.basis, &*self.quad, el.etype, &xyz, &mat, &cfg, &mut self.local)?;
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let Some(bc) = bc_inputs(ctx, e)? else { continue };
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                stat_elec_boundary_local(&*self.basis, &*self.quad, el.etype, &xyz, &bc, &cfg, &mut self.local)?;
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }

        if POST_FIELDS.iter().any(|f| flag(ctx.params(), f.flag)) || flag(ctx.params(), "calculate average potential") {
            self.postprocess(ctx, &cfg)?;
        }
        Ok(())
    }
}

impl StatElecSolver {
    fn bulk_inputs(&self, ctx: &SolverCtx<'_>, e: ElemId, n: usize) -> Result<StatElecAtNodes, SolverError> {
        let bp = ctx
            .body_params(e)
            .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
        Ok(StatElecAtNodes {
            relative_permittivity: real_at_nodes(ctx, bp.material, "Relative Permittivity", e)?
                .unwrap_or_else(|| vec![0.0; n]),
            charge_density: real_at_nodes(ctx, bp.body_force, "Charge Density", e)?,
        })
    }

    /// `StatElecSolver_post`: the element integrals of the requested fields,
    /// their nodal projection, the weighted nodal averages and, with any
    /// distributed field, the total energy, the potential difference and the
    /// effective capacitance.
    fn postprocess(&mut self, ctx: &mut SolverCtx<'_>, cfg: &StatElecConfig) -> Result<(), SolverError> {
        const POST: &str = "StatElecSolver_post";
        let dim = cfg.dim;
        let constant_weights = flag(ctx.params(), "constant weights");
        let calc_ave = flag(ctx.params(), "calculate average potential");

        ctx.info(POST, "------------------------------------------------", 4);
        ctx.info(POST, "Calculating postprocessing fields", 4);

        // The fields that exist, with their accumulators (Values = 0 first).
        let mut active: Vec<(&PostField, usize, Vec<f64>)> = Vec::new();
        for f in &POST_FIELDS {
            if let Some(v) = ctx.sim.vars.get(f.name) {
                active.push((f, v.dofs, vec![0.0; v.values.len()]));
            }
        }
        let need_scaling = active.iter().any(|(f, _, _)| !f.nodal);
        let n_slots = ctx
            .primary()
            .ok_or_else(|| SolverError::fatal(POST, "solver has no variable"))?
            .values
            .len();
        let mut weight_vector = vec![0.0; n_slots];
        let n_bodies = ctx.sim.model.bodies.len();
        let (mut pot_int, mut pot_vol) = (vec![0.0; n_bodies], vec![0.0; n_bodies]);
        let mut energy_tot = 0.0;

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for e in elems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nd| ctx.mesh.coord(nd)).collect();
            let eps_r = self.bulk_inputs(ctx, e, n)?.relative_permittivity;
            ctx.element_slots(e, &mut slots);
            let pot: Vec<f64> = {
                let p = ctx.primary().expect("checked above");
                slots
                    .iter()
                    .map(|&s| if s == NO_DOF { 0.0 } else { p.values[s as usize] })
                    .collect()
            };
            stat_elec_post_local(&*self.basis, &*self.quad, el.etype, &xyz, &pot, &eps_r, cfg, &mut self.post)?;
            let post = &self.post;
            energy_tot += post.energy_integral;
            if calc_ave {
                if let Some(b) = ctx.body_params(e) {
                    pot_vol[b.body_id - 1] += post.volume;
                    pot_int[b.body_id - 1] += post.potential_integral;
                }
            }
            if need_scaling {
                for (i, &s) in slots.iter().enumerate().filter(|(_, &s)| s != NO_DOF) {
                    weight_vector[s as usize] += if constant_weights { 1.0 } else { post.weights[i] };
                }
            }

            let lu = DenseLu::new(post.mass.clone(), n)?;
            for field_type in 1..=3u8 {
                let comps = if field_type == 1 { 1 } else { dim };
                for m in 0..comps {
                    let raw: &[f64] = match field_type {
                        1 => &post.energy,
                        2 => &post.flux[m],
                        _ => &post.field[m],
                    };
                    let mut x = raw.to_vec();
                    let mut solved = false;
                    for (f, vdofs, acc) in active.iter_mut().filter(|(f, _, _)| f.field_type == field_type) {
                        if m >= *vdofs {
                            continue;
                        }
                        if !f.nodal && !solved {
                            lu.solve(&mut x);
                            solved = true;
                        }
                        for (i, &s) in slots.iter().enumerate().filter(|(_, &s)| s != NO_DOF) {
                            acc[s as usize * *vdofs + m] += if constant_weights { x[i] } else { post.weights[i] * x[i] };
                        }
                    }
                }
            }
        }

        if need_scaling {
            let half = energy_tot / 2.0;
            ctx.info(POST, &format!("Total Electric Energy: {half:12.3E}"), 6);
            let sim_params = &mut ctx.sim.model.simulation;
            sim_params.insert("res: Electric Energy", Value::real(half));
            let fixed = ctx.dirichlet_values()?;
            let pot_diff = if fixed.is_empty() {
                0.0
            } else {
                let max = fixed.iter().fold(f64::NEG_INFINITY, |m, &(_, v)| m.max(v));
                let min = fixed.iter().fold(f64::INFINITY, |m, &(_, v)| m.min(v));
                max - min
            };
            if pot_diff > f64::MIN_POSITIVE {
                ctx.sim.model.simulation.insert("res: Potential Difference", Value::real(pot_diff));
                let capacitance = 2.0 * half / (pot_diff * pot_diff);
                ctx.info(POST, &format!("Effective Capacitance: {capacitance:12.3E}"), 6);
                ctx.sim.model.simulation.insert("RES: Effective Capacitance", Value::real(capacitance));
            }
            for (_, vdofs, acc) in active.iter_mut().filter(|(f, _, _)| !f.nodal) {
                for slot in 0..acc.len() / *vdofs {
                    let w = weight_vector[slot];
                    if w.abs() > f64::MIN_POSITIVE {
                        for m in 0..*vdofs {
                            acc[slot * *vdofs + m] /= w;
                        }
                    }
                }
            }

            // Physical scaling of the solution, if requested.
            let p = ctx.params();
            let mut control = p.real("energy control").map(|t| (t / half).sqrt());
            if control.is_none() {
                if let Some(q) = p.real("charge control") {
                    if pot_diff < f64::MIN_POSITIVE {
                        return Err(SolverError::fatal(POST, "Charge cannot be controlled without pot. difference"));
                    }
                    control = Some(q / (half / pot_diff));
                }
            }
            if let Some(s) = control {
                ctx.info(POST, &format!("Control Scaling: {s:12.3E}"), 4);
                ctx.sim.model.simulation.insert("RES: StatElec Scaling", Value::real(s));
                if let Some(p) = ctx.primary_mut() {
                    p.values.iter_mut().for_each(|v| *v *= s);
                }
                for (f, _, acc) in active.iter_mut() {
                    let k = if f.field_type == 1 { s * s } else { s };
                    acc.iter_mut().for_each(|v| *v *= k);
                }
            }
        }

        for (f, _, acc) in active {
            let v = ctx.sim.vars.get_mut(f.name).expect("the field exists");
            v.values.copy_from_slice(&acc);
            v.values_changed = true;
        }

        if calc_ave {
            for i in 0..n_bodies {
                if pot_vol[i] < f64::EPSILON {
                    continue;
                }
                let ave = pot_int[i] / pot_vol[i];
                ctx.info(POST, &format!("Average body{} potential: {ave:12.5E}", i + 1), 7);
                ctx.sim
                    .model
                    .simulation
                    .insert(&format!("res: Average body{} potential", i + 1), Value::real(ave));
            }
        }
        Ok(())
    }
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
