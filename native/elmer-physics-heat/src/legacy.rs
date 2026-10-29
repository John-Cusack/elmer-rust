//! The `Legacy Assembly = True` driver: `HeatSolve` delegating to
//! HeatSolveLegacy.F90 (`DelegateToHeatLegacy`), restricted to what the
//! kernels here implement.
//!
//! Beyond the vectorised path this adds
//! * `Smart Heater Control`: the heater body force does not enter the right
//!   hand side. It is assembled into a second vector and the system is solved
//!   for both, `T = Y + P X` with `A X = heater load`, `A Y = rhs` and the
//!   power scaling `P` chosen so a control node reaches `Smart Heater
//!   Temperature` (HeatSolveLegacy.F90:1176-1270);
//! * `Radiation = Idealized`: grey radiation to `Radiation External
//!   Temperature` (or `External Temperature`), linearised by Newton's method
//!   once `Nonlinear System Newton After Iterations/Tolerance` allows it and
//!   by the Picard form before that (`AddHeatFluxBC`, :1477).
//!
//! The legacy assembly forms the heat source and the boundary coefficients
//! from nodal values and hands them to the same element integrals, so the
//! kernels in [`crate::bulk`] and [`crate::boundary`] serve both paths.
//! Transient simulations, the `Smart Heater Boundary` control node and
//! integral heat control are refused.

use elmer_core::{ElemId, Params, Value, NO_DOF};
use elmer_solver::{SolverCtx, SolverError};

use crate::material::{HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, Stabilization};
use crate::preflight::stabilized;
use crate::solver::{axisymmetric, bulk_inputs, real_at_nodes, HeatSolver};
use crate::{heat_boundary_local, heat_bulk_local};

const CALLER: &str = "HeatSolve";

/// What the `Smart Heater Control` body force asks for.
struct Smart {
    /// `Smart Heater Temperature`: the value the control node is driven to.
    melt_point: f64,
    /// Mesh node (0-based) whose temperature is controlled.
    node: u32,
    /// `Smart Heater Control After Tolerance`, 1.0 when absent.
    tol: f64,
}

fn smart_heater(ctx: &SolverCtx<'_>) -> Result<Option<Smart>, SolverError> {
    let model = &ctx.sim.model;
    if !model
        .body_forces
        .iter()
        .any(|b| b.get("Smart Heater Control").is_some())
    {
        return Ok(None);
    }
    let Some(bf) = model
        .body_forces
        .iter()
        .rev()
        .find(|b| b.logical("Smart Heater Control") == Some(true))
    else {
        return Err(SolverError::Unsupported(
            "Smart Heater Control present but never True".into(),
        ));
    };
    let melt = ctx.creal(bf, "Smart Heater Temperature")?;
    let mut node = bf.integer("Smart Heater Control Node").map(|n| n - 1);
    if node.is_none() {
        if let Some(Value::Real { data, .. }) = bf.get("Smart Heater Control Point") {
            if data.len() >= 3 {
                let cp = [data[0], data[1], data[2]];
                let perm = ctx.primary().and_then(|p| p.perm.map(<[u32]>::to_vec));
                let mut best = f64::MAX;
                for l in 0..ctx.mesh.n_nodes() as u32 {
                    if let Some(p) = &perm {
                        if p[l as usize] == NO_DOF {
                            continue;
                        }
                    }
                    let x = ctx.mesh.coord(l);
                    let d: f64 = (0..3).map(|i| (cp[i] - x[i]).powi(2)).sum();
                    if d < best {
                        best = d;
                        node = Some(l as i32);
                    }
                }
                ctx.info(
                    CALLER,
                    &format!("Found Control Point at distance: {}", best.sqrt()),
                    5,
                );
            }
        }
    }
    let (Some(melt_point), Some(node)) = (melt, node) else {
        return Err(SolverError::Unsupported(
            "Smart Heater Boundary: a control node found from the boundary condition".into(),
        ));
    };
    if ctx.params().logical("Smart Heater Average") == Some(true) {
        return Err(SolverError::Unsupported("Smart Heater Average".into()));
    }
    let tol = ctx
        .params()
        .real("Smart Heater Control After Tolerance")
        .unwrap_or(1.0);
    Ok(Some(Smart {
        melt_point,
        node: node as u32,
        tol,
    }))
}

/// `TabulateBoundaryAverages`: the element mean of the boundary condition's
/// `Emissivity`, else of the parent material's.
fn emissivity(ctx: &SolverCtx<'_>, e: ElemId, n: usize) -> Result<f64, SolverError> {
    if let Some(bc) = ctx.bc_params(e) {
        if let Some(v) = real_at_nodes(ctx, Some(bc), "Emissivity", e)? {
            return Ok(v.iter().sum::<f64>() / n as f64);
        }
    }
    for parent in ctx.mesh.element(e).parents.into_iter().flatten() {
        let mat = ctx.body_params(parent).and_then(|b| b.material);
        if let Some(v) = real_at_nodes(ctx, mat, "Emissivity", e)? {
            return Ok(v.iter().sum::<f64>() / n as f64);
        }
    }
    Ok(0.0)
}

impl HeatSolver {
    pub(crate) fn solve_legacy(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        transient: bool,
    ) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported(
                "Legacy Assembly with a transient simulation".into(),
            ));
        }
        let axisym = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let params = ctx.params();
        let maxiter = params
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let supg = stabilized(params);
        let long_edge = ctx
            .simulation_params()
            .logical("Stabilization Use Longest Element Edge")
            .unwrap_or(false);

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let is_radiation = ctx
            .sim
            .model
            .boundary_conditions
            .iter()
            .any(|b| b.get("Radiation").is_some());
        let stefan = if is_radiation {
            Some(
                ctx.sim
                    .model
                    .constants
                    .real("Stefan Boltzmann")
                    .ok_or_else(|| {
                        SolverError::fatal(CALLER, "Stefan Boltzmann is not defined".to_string())
                    })?,
            )
        } else {
            None
        };
        let (newton_tol, newton_iter) = if is_radiation {
            (
                ctx.params()
                    .real("Nonlinear System Newton After Tolerance")
                    .unwrap_or(0.0),
                ctx.params()
                    .integer("Nonlinear System Newton After Iterations")
                    .unwrap_or(0),
            )
        } else {
            (1.0, 0)
        };
        let mut newton = newton_iter == 0;

        let smart = smart_heater(ctx)?;
        // Without `Smart Heater Control After Tolerance` the control is active
        // from the start.
        let mut smart_tol_reached = smart.is_some()
            && ctx
                .params()
                .real("Smart Heater Control After Tolerance")
                .is_none();
        if ctx
            .params()
            .real("Smart Heater Time Scale")
            .is_some_and(|t| t > 0.0)
        {
            return Err(SolverError::Unsupported("Smart Heater Time Scale".into()));
        }

        let var_name = ctx.variable_name().unwrap_or("temperature").to_string();
        let slot_ix = ctx
            .sim
            .slots
            .iter()
            .position(|s| s.id == ctx.solver_id())
            .ok_or_else(|| SolverError::fatal(CALLER, "solver slot not found".to_string()))?;
        let relax_saved = ctx.sim.slots[slot_ix]
            .params
            .get("Nonlinear System Relaxation Factor")
            .cloned();

        let mut bulk_cfg = HeatBulkConfig::new(dim);
        bulk_cfg.axisymmetric = axisym;
        let bnd_cfg = HeatBoundaryConfig {
            axisymmetric: axisym,
            quad_points: None,
        };
        let mut dofs = Vec::new();
        let mut xx: Vec<f64> = Vec::new();
        let mut yy: Vec<f64> = Vec::new();
        // The last bulk element's flag outlives the element loop, as the
        // Fortran's `HeaterControlLocal` does when the boundary loop reads it.
        let mut heater_local = false;

        for iter in 1..=maxiter {
            ctx.info(CALLER, &format!("TEMPERATURE ITERATION {iter}"), 4);
            ctx.default_initialize()?;
            let n_dof = ctx.system().map_or(0, |s| s.rhs.len());
            let mut force_heater = vec![0.0; n_dof];
            let control = smart.is_some() && newton && smart_tol_reached;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = bulk_inputs(ctx, e, n, true)?;
                bulk_cfg.stabilization = if supg && mat.convection_velocity.is_some() {
                    let (m_k, h_k) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                        .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
                    Stabilization::Supg { h_k, m_k }
                } else {
                    Stabilization::None
                };
                heat_bulk_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &bulk_cfg,
                    &mut self.local,
                )?;
                heater_local = control
                    && ctx
                        .body_params(e)
                        .and_then(|b| b.body_force)
                        .is_some_and(|b| b.logical("Smart Heater Control") == Some(true));
                if heater_local {
                    ctx.element_dofs(e, &mut dofs);
                    for (p, &d) in dofs.iter().enumerate() {
                        if d != NO_DOF {
                            force_heater[d as usize] += self.local.force[p];
                        }
                    }
                    self.local.force.iter_mut().for_each(|f| *f = 0.0);
                }
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                let Some(bc) = ctx.bc_params(e) else {
                    continue;
                };
                // A Dirichlet node carries no heater load.
                if heater_local && bc.get(&var_name).is_some() {
                    ctx.element_dofs(e, &mut dofs);
                    for &d in dofs.iter().filter(|&&d| d != NO_DOF) {
                        force_heater[d as usize] = 0.0;
                    }
                }
                if bc.logical("Heat Flux BC") == Some(false) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let mut coeff = vec![0.0; n];
                let mut load = vec![0.0; n];
                let mut any = false;
                let radiation = bc.string("Radiation").map(str::to_ascii_lowercase);
                if radiation.as_deref() == Some("idealized") {
                    let sigma = stefan.expect("Stefan Boltzmann is read when a BC radiates");
                    let emis = emissivity(ctx, e, n)?;
                    let text =
                        match real_at_nodes(ctx, Some(bc), "Radiation External Temperature", e)? {
                            Some(t) => t,
                            None => real_at_nodes(ctx, Some(bc), "External Temperature", e)?
                                .unwrap_or_else(|| vec![0.0; n]),
                        };
                    let temps: Vec<f64> = {
                        let f = ctx.primary().expect("heat solver has a variable");
                        el.nodes
                            .iter()
                            .map(|&node| {
                                let s = f.perm.map_or(node, |p| p[node as usize]);
                                f.values[s as usize]
                            })
                            .collect()
                    };
                    for j in 0..n {
                        let (t, x) = (temps[j], text[j]);
                        if newton {
                            coeff[j] = emis * 4.0 * t.powi(3) * sigma;
                            load[j] = emis * (3.0 * t.powi(4) + x.powi(4)) * sigma;
                        } else {
                            coeff[j] =
                                emis * (t.powi(3) + t * t * x + t * x * x + x.powi(3)) * sigma;
                            load[j] = coeff[j] * x;
                        }
                    }
                    any = true;
                }
                if let Some(h) = real_at_nodes(ctx, Some(bc), "Heat Transfer Coefficient", e)? {
                    let text = real_at_nodes(ctx, Some(bc), "External Temperature", e)?
                        .unwrap_or_else(|| vec![0.0; n]);
                    for j in 0..n {
                        load[j] += h[j] * text[j];
                        coeff[j] += h[j];
                    }
                    any = true;
                }
                if let Some(q) = real_at_nodes(ctx, Some(bc), "Heat Flux", e)? {
                    for j in 0..n {
                        load[j] += q[j];
                    }
                    any = true;
                }
                if !any {
                    continue;
                }
                let at_nodes = HeatBcAtNodes {
                    heat_flux: Some(load),
                    heat_transfer_coefficient: Some(coeff),
                    external_temperature: None,
                };
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                heat_boundary_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &at_nodes,
                    &bnd_cfg,
                    &mut self.local,
                )?;
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;

            if let (Some(smart), true) = (&smart, control) {
                self.solve_controlled(
                    ctx,
                    slot_ix,
                    smart,
                    &force_heater,
                    &mut xx,
                    &mut yy,
                    relax_saved.as_ref(),
                )?;
            } else {
                ctx.default_solve()?;
            }

            let (change, converged) = {
                let c = ctx.convergence().expect("heat solver has a variable");
                (c.nonlin_change, c.nonlin_converged > 0)
            };
            ctx.info(CALLER, &format!("Relative Change : {change}"), 4);
            if change < newton_tol || iter >= newton_iter {
                newton = true;
            }
            if converged && (smart.is_none() || smart_tol_reached) {
                break;
            }
            if let Some(s) = &smart {
                if change < s.tol {
                    smart_tol_reached = true;
                    yy = ctx.primary().expect("variable").values.to_vec();
                }
            }
        }
        Ok(())
    }

    /// The control branch of the nonlinear iteration: solve for the heater
    /// load (`X`) and the remaining load (`Y`) without the nonlinear
    /// bookkeeping, then scale the heater so the control node hits its
    /// target and compute the change of the combined field.
    #[allow(clippy::too_many_arguments)]
    fn solve_controlled(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        slot_ix: usize,
        smart: &Smart,
        force_heater: &[f64],
        xx: &mut Vec<f64>,
        yy: &mut Vec<f64>,
        relax_saved: Option<&Value>,
    ) -> Result<(), SolverError> {
        let prev = ctx.primary().expect("variable").values.to_vec();
        let n = prev.len();
        xx.resize(n, 0.0);
        yy.resize(n, 0.0);

        let slot = &mut ctx.sim.slots[slot_ix];
        slot.params
            .insert("Skip Compute Nonlinear Change", Value::Logical(true));
        slot.params
            .insert("Skip Advance Nonlinear Iter", Value::Logical(true));
        let relaxed = slot
            .params
            .real("Nonlinear System Relaxation Factor")
            .is_some_and(|r| r != 1.0);
        if relaxed {
            slot.params
                .insert("Nonlinear System Relaxation Factor", Value::real(1.0));
        }

        let rhs0 = ctx.system().expect("assembled system").rhs.clone();
        ctx.system_mut().expect("assembled system").rhs[..force_heater.len()]
            .copy_from_slice(force_heater);
        ctx.primary_mut()
            .expect("variable")
            .values
            .copy_from_slice(xx);
        ctx.default_solve()?;
        xx.copy_from_slice(ctx.primary().expect("variable").values);

        ctx.system_mut().expect("assembled system").rhs = rhs0;
        ctx.primary_mut()
            .expect("variable")
            .values
            .copy_from_slice(yy);
        ctx.default_solve()?;
        yy.copy_from_slice(ctx.primary().expect("variable").values);

        let slot = &mut ctx.sim.slots[slot_ix];
        slot.params
            .insert("Skip Compute Nonlinear Change", Value::Logical(false));
        slot.params
            .insert("Skip Advance Nonlinear Iter", Value::Logical(false));

        let k = {
            let f = ctx.primary().expect("variable");
            f.perm.map_or(smart.node, |p| p[smart.node as usize]) as usize
        };
        let power = (smart.melt_point - yy[k]) / xx[k];
        ctx.info(CALLER, &format!("Heater scaling: {power}"), 4);
        let id = ctx.sim.slots[slot_ix].variable.expect("variable");
        let var = ctx.sim.vars.by_id_mut(id);
        var.nonlin_values = Some(prev);
        for (t, (x, y)) in var.values.iter_mut().zip(xx.iter().zip(yy.iter())) {
            *t = y + power * x;
        }
        var.values_changed = true;
        ctx.nonlinear_change()?;

        // The relaxation applies to the iteration again from here on.
        if relaxed {
            if let Some(r) = relax_saved {
                ctx.sim.slots[slot_ix]
                    .params
                    .insert("Nonlinear System Relaxation Factor", r.clone());
            }
        }
        Ok(())
    }
}
