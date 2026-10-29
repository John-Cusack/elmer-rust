//! The postprocessing solver `StatElecSolver_post` (StatElecSolve.F90:1010)
//! and `StatCurrentSolver_post` (StatCurrentSolve.F90:830) share one
//! structure: element integrals of the requested fields, their local nodal
//! projection, weighted nodal averaging, the totals and the control scaling.
//! [`PostSpec`] holds what differs between the two modules.

use std::sync::Arc;

use elmer_core::{ElemId, ElementBasis, FieldKind, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{SolverCtx, SolverError, Variable};

use crate::common::slot_mut;
use crate::potential::{
    potential_post_local, Coefficient, DenseLu, PostConfig, PotentialPostLocal,
};

/// A postprocessing field. Nodal fields (`Nodal Energy Density`, `Nodal
/// Electric Flux`) come before the distributed ones of the same type, as the
/// Fortran requires.
pub(crate) struct PostField {
    pub name: &'static str,
    pub flag: &'static str,
    /// 1 energy density or heating (one component), 2 flux or current, 3 field.
    pub field_type: u8,
    pub nodal: bool,
}

impl PostField {
    pub fn dofs(&self, dim: usize) -> usize {
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

/// What the total integral is divided by to give the module's second result.
pub(crate) enum Derived {
    /// `2 W / dV^2`, the effective capacitance.
    Capacitance(&'static str),
    /// `dV^2 / P`, the effective resistance.
    Resistance(&'static str),
}

pub(crate) struct PostSpec {
    pub caller: &'static str,
    pub fields: &'static [PostField],
    /// 0.5 turns the stored `int E.D` into the energy; 1 keeps the power.
    pub total_factor: f64,
    pub total_key: &'static str,
    pub potential_difference_key: Option<&'static str>,
    pub derived: Derived,
    /// Nodal weights at or below this are not divided by.
    pub weight_eps: f64,
    /// Keywords scaling the total (`Energy Control`) and the flux through the
    /// potential difference (`Charge Control`) to a target.
    pub control_total: &'static str,
    pub control_flux: &'static str,
    pub scaling_key: &'static str,
    /// Run the global scaling for a control keyword even with no
    /// distributed field (StatCurrentSolve).
    pub scale_on_control: bool,
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// Whether the post solver does anything for this solver section.
pub(crate) fn post_active(p: &ParamMap, spec: &PostSpec) -> bool {
    spec.fields.iter().any(|f| flag(p, f.flag))
        || flag(p, "calculate average potential")
        || (spec.scale_on_control
            && (p.get(spec.control_total).is_some() || p.get(spec.control_flux).is_some()))
}

/// The exported fields of `_init`, on the permutation of the primary variable.
pub(crate) fn create_fields(
    ctx: &mut SolverCtx<'_>,
    spec: &PostSpec,
    dim: usize,
) -> Result<(), SolverError> {
    let perm: Option<Arc<[u32]>> = slot_mut(ctx).perm.clone();
    let n_slots = ctx
        .primary()
        .ok_or_else(|| SolverError::fatal(spec.caller, "solver has no variable"))?
        .values
        .len();
    let wanted: Vec<&PostField> = spec
        .fields
        .iter()
        .filter(|f| flag(ctx.params(), f.flag))
        .collect();
    for f in wanted {
        if ctx.sim.vars.id(f.name).is_some() {
            continue;
        }
        let dofs = f.dofs(dim);
        let mut v = Variable::new(
            f.name,
            FieldKind::Nodal,
            dofs,
            perm.clone(),
            vec![0.0; n_slots * dofs],
        );
        v.secondary = true;
        ctx.sim.vars.add_with_components(v, &f.components(dim))?;
    }
    Ok(())
}

pub(crate) type CoefOf<'a> =
    &'a dyn Fn(&SolverCtx<'_>, ElemId, usize) -> Result<Coefficient, SolverError>;

/// The postprocessing pass after a solve.
#[allow(clippy::too_many_arguments)]
pub(crate) fn postprocess(
    ctx: &mut SolverCtx<'_>,
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    scratch: &mut PotentialPostLocal,
    spec: &PostSpec,
    cfg: &PostConfig,
    scale: f64,
    coef_of: CoefOf<'_>,
) -> Result<(), SolverError> {
    let caller = spec.caller;
    let dim = cfg.dim;
    let constant_weights = flag(ctx.params(), "constant weights");
    let calc_ave = flag(ctx.params(), "calculate average potential");
    let control = spec.scale_on_control
        && (ctx.params().get(spec.control_total).is_some()
            || ctx.params().get(spec.control_flux).is_some());

    ctx.info(
        caller,
        "------------------------------------------------",
        4,
    );
    ctx.info(caller, "Calculating postprocessing fields", 4);

    // The fields that exist, with their accumulators (Values = 0 first).
    let mut active: Vec<(&PostField, usize, Vec<f64>)> = Vec::new();
    for f in spec.fields {
        if let Some(v) = ctx.sim.vars.get(f.name) {
            active.push((f, v.dofs, vec![0.0; v.values.len()]));
        }
    }
    let need_scaling = active.iter().any(|(f, _, _)| !f.nodal);
    let n_slots = ctx
        .primary()
        .ok_or_else(|| SolverError::fatal(caller, "solver has no variable"))?
        .values
        .len();
    let mut weight_vector = vec![0.0; n_slots];
    let n_bodies = ctx.sim.model.bodies.len();
    let (mut pot_int, mut pot_vol) = (vec![0.0; n_bodies], vec![0.0; n_bodies]);
    let mut total = 0.0;

    let elems: Vec<ElemId> = ctx.active_elements().collect();
    let mut slots = Vec::new();
    for e in elems {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nd| ctx.mesh.coord(nd)).collect();
        let coef = coef_of(ctx, e, n)?;
        ctx.element_slots(e, &mut slots);
        let pot: Vec<f64> = {
            let p = ctx.primary().expect("checked above");
            slots
                .iter()
                .map(|&s| {
                    if s == NO_DOF {
                        0.0
                    } else {
                        p.values[s as usize]
                    }
                })
                .collect()
        };
        potential_post_local(
            basis, quad, el.etype, &xyz, &pot, &coef, scale, cfg, scratch,
        )?;
        let post = &*scratch;
        total += post.heat_integral;
        if calc_ave {
            if let Some(b) = ctx.body_params(e) {
                pot_vol[b.body_id - 1] += post.volume;
                pot_int[b.body_id - 1] += post.potential_integral;
            }
        }
        if need_scaling {
            for (i, &s) in slots.iter().enumerate().filter(|(_, &s)| s != NO_DOF) {
                weight_vector[s as usize] += if constant_weights {
                    1.0
                } else {
                    post.weights[i]
                };
            }
        }

        let lu = DenseLu::new(post.mass.clone(), n)?;
        for field_type in 1..=3u8 {
            let comps = if field_type == 1 { 1 } else { dim };
            for m in 0..comps {
                let raw: &[f64] = match field_type {
                    1 => &post.heat,
                    2 => &post.flux[m],
                    _ => &post.field[m],
                };
                let mut x = raw.to_vec();
                let mut solved = false;
                for (f, vdofs, acc) in active
                    .iter_mut()
                    .filter(|(f, _, _)| f.field_type == field_type)
                {
                    if m >= *vdofs {
                        continue;
                    }
                    if !f.nodal && !solved {
                        lu.solve(&mut x);
                        solved = true;
                    }
                    for (i, &s) in slots.iter().enumerate().filter(|(_, &s)| s != NO_DOF) {
                        acc[s as usize * *vdofs + m] += if constant_weights {
                            x[i]
                        } else {
                            post.weights[i] * x[i]
                        };
                    }
                }
            }
        }
    }

    if need_scaling || control {
        let total = total * spec.total_factor;
        ctx.info(
            caller,
            &format!("Total {}: {total:12.3E}", spec.total_key),
            6,
        );
        ctx.sim
            .model
            .simulation
            .insert(spec.total_key, Value::real(total));
        let fixed = ctx.dirichlet_values()?;
        let pot_diff = if fixed.is_empty() {
            0.0
        } else {
            let max = fixed.iter().fold(f64::NEG_INFINITY, |m, &(_, v)| m.max(v));
            let min = fixed.iter().fold(f64::INFINITY, |m, &(_, v)| m.min(v));
            max - min
        };
        if pot_diff > f64::MIN_POSITIVE {
            if let Some(k) = spec.potential_difference_key {
                ctx.sim.model.simulation.insert(k, Value::real(pot_diff));
            }
            let (key, v) = match spec.derived {
                Derived::Capacitance(k) => (k, 2.0 * total / (pot_diff * pot_diff)),
                Derived::Resistance(k) => (k, pot_diff * pot_diff / total),
            };
            ctx.info(caller, &format!("{key}: {v:12.3E}"), 6);
            ctx.sim.model.simulation.insert(key, Value::real(v));
        }
        for (f, vdofs, acc) in active.iter_mut() {
            if f.nodal {
                continue;
            }
            for slot in 0..acc.len() / *vdofs {
                let w = weight_vector[slot];
                if w.abs() > spec.weight_eps {
                    for m in 0..*vdofs {
                        acc[slot * *vdofs + m] /= w;
                    }
                }
            }
        }

        // Physical scaling of the solution, if requested.
        let p = ctx.params();
        let mut factor = p.real(spec.control_total).map(|t| (t / total).sqrt());
        if factor.is_none() {
            if let Some(q) = p.real(spec.control_flux) {
                if pot_diff < f64::MIN_POSITIVE {
                    return Err(SolverError::fatal(
                        caller,
                        "Flux cannot be controlled without pot. difference",
                    ));
                }
                factor = Some(q / (total / pot_diff));
            }
        }
        if let Some(s) = factor {
            ctx.info(caller, &format!("Control Scaling: {s:12.3E}"), 4);
            ctx.sim
                .model
                .simulation
                .insert(spec.scaling_key, Value::real(s));
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
            ctx.info(
                caller,
                &format!("Average body{} potential: {ave:12.5E}", i + 1),
                7,
            );
            ctx.sim.model.simulation.insert(
                &format!("res: Average body{} potential", i + 1),
                Value::real(ave),
            );
        }
    }
    Ok(())
}
