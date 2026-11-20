//! `Apply Limiter`: the soft limiter `DefaultDirichletBCs` sets before the
//! boundary conditions (`DetermineSoftLimiter`, SolveCore.F90:7912). The
//! dofs whose value passes the `<var> Upper Limit` / `<var> Lower Limit` of
//! their element's body force are fixed at the limit, and leave the set
//! again when their nodal load (`CalculateLoads` of the assembled system,
//! kept in `<var> Contact Load`) pushes back; in the first nonlinear
//! iteration of the first timestep the set starts from the values, and in
//! the first nonlinear iteration of a later timestep it is kept as it is.
//!
//! The load of dof `ind` is read as `LoadValues(ind)`, `ind` the dof's
//! index in the solver variable, from the `Contact Load` field that
//! `CalculateLoads` fills in the field's own numbering: the two are the
//! same unless the variable's permutation has changed since the field was
//! made (CutFEM), and then the Fortran reads another node's load. That is
//! kept.
//!
//! Refused: limits from boundary conditions, `Linear System Limiter`,
//! `Apply Limiter Loads After`, the conservative and downstream options,
//! relative tolerances, the global active/passive conditions, given
//! initial or active sets and `Save Limiter`.

use std::sync::Arc;

use elmer_core::{normalize_key, ElemId, Params, NO_DOF};

use crate::error::SolverError;
use crate::solver::SolverCtx;
use crate::variable::Variable;

const CALLER: &str = "DetermineSoftLimiter";

/// `Var % LowerLimitActive` / `Var % UpperLimitActive` of the solver's
/// variable, kept between calls.
#[derive(Clone, Debug, Default)]
pub struct LimiterState {
    active: [Vec<bool>; 2],
}

impl SolverCtx<'_> {
    /// `DetermineSoftLimiter(Solver)`: the dofs to fix at their limits.
    pub(crate) fn determine_soft_limiter(&mut self) -> Result<Vec<(u32, f64)>, SolverError> {
        let p = self.params().clone();
        for k in [
            "linear system limiter",
            "apply limiter loads after",
            "apply limiter conservative add after iterations",
            "apply limiter conservative remove after iterations",
            "apply limiter remove downstream",
            "limiter relative tolerance",
            "limiter global active condition",
            "limiter global passive condition",
            "save limiter",
        ] {
            if p.get(k).is_some_and(|v| !matches!(v, elmer_core::Value::Logical(false))) {
                return Err(SolverError::Unsupported(format!("Apply Limiter with {k}")));
            }
        }
        let id = self.sim.slots[self.slot].variable.ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
        let (name, dofs) = {
            let v = self.sim.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        if dofs != 1 {
            return Err(SolverError::Unsupported("Apply Limiter for a variable of more than one dof".into()));
        }
        let it = |k: &str| self.sim.vars.scalar(k).map_or(1, |x| x.round() as i64);
        let (nonlin, coupled, time) = (it("nonlin iter"), it("coupled iter"), it("timestep"));
        let first_time = nonlin <= 1 && coupled <= 1 && time == 1;
        let mut freeze = time > 1 && nonlin <= 1 && coupled <= 1;
        if time == 1 && !first_time {
            freeze = p.logical("limiter passive first timestep") == Some(true);
        }
        let load_eps = p.real("limiter load tolerance").unwrap_or(1.0e-8);
        let val_eps = p.real("limiter value tolerance").unwrap_or(1.0e-8);
        let eq_sign = if p.logical("limiter load sign negative") == Some(true) { -1.0 } else { 1.0 };

        // CalculateLoads(Solver, Solver % Matrix, Var % Values, ...): A x - b
        // of the system as it is now, into `<var> Contact Load`.
        let load_name = normalize_key(&format!("{name} contact load"));
        let load_id = match self.sim.vars.id(&load_name) {
            Some(l) => l,
            None => {
                // Made with the variable (AddEquationBasics): its own
                // permutation and size, before any CutFEM swap.
                let (perm, len) = match self.cutfem_original() {
                    Some((perm, len)) => (perm, len),
                    None => {
                        let v = self.sim.vars.by_id(id);
                        (v.perm.clone(), v.values.len())
                    }
                };
                let output = self.sim.vars.by_id(id).output;
                let mut lv = Variable::new(&load_name, elmer_core::FieldKind::Nodal, 1, perm, vec![0.0; len]);
                lv.output = output;
                self.sim.vars.add(lv)?
            }
        };
        let temp = {
            let sys = self.sim.slots[self.slot].system.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no system"))?;
            let x = &self.sim.vars.by_id(id).values;
            let mut t = vec![0.0; sys.matrix.nrows];
            crate::matvec::crs_matrix_vector_multiply(&sys.matrix, 0, x, &mut t);
            for (ti, bi) in t.iter_mut().zip(&sys.rhs) {
                *ti -= bi;
            }
            t
        };
        {
            let var_perm: Option<Arc<[u32]>> = self.sim.vars.by_id(id).perm.clone();
            let lv = self.sim.vars.by_id_mut(load_id);
            let lperm = lv.perm.clone();
            let n = lperm.as_ref().map_or(lv.values.len(), |p| p.len());
            for i in 0..n {
                let j = match &lperm {
                    Some(p) => p[i],
                    None => i as u32,
                };
                let k = match &var_perm {
                    Some(p) => p.get(i).copied().unwrap_or(NO_DOF),
                    None => i as u32,
                };
                if j != NO_DOF && k != NO_DOF {
                    lv.values[j as usize] = temp[k as usize];
                }
            }
        }
        let loads = self.sim.vars.by_id(load_id).values.clone();

        let totsize = self.sim.vars.by_id(id).values.len();
        let nb = self.mesh.n_bulk();
        let mut fixed: Vec<(u32, f64)> = Vec::new();
        for upper in [0usize, 1] {
            let limit_sign = if upper == 0 { -eq_sign } else { eq_sign };
            let (limit_name, init_name, active_name) = if upper == 0 {
                (format!("{name} lower limit"), format!("{name} lower initial"), format!("{name} lower active"))
            } else {
                (format!("{name} upper limit"), format!("{name} upper initial"), format!("{name} upper active"))
            };
            let model = &self.sim.model;
            if model.boundary_conditions.iter().any(|b| b.get(&limit_name).is_some()) {
                return Err(SolverError::Unsupported(format!("{limit_name} in a boundary condition")));
            }
            if !model.body_forces.iter().any(|b| b.get(&limit_name).is_some()) {
                continue;
            }
            for sec in model.body_forces.iter().chain(&model.initial_conditions) {
                for k in [&init_name, &active_name, &format!("{active_name} initial")] {
                    if sec.get(k).is_some() {
                        return Err(SolverError::Unsupported(format!("Apply Limiter with {k}")));
                    }
                }
            }
            if model.initial_conditions.iter().any(|ic| ic.get(&limit_name).is_some()) {
                return Err(SolverError::Unsupported(format!("{limit_name} in an initial condition")));
            }
            self.info(CALLER, &format!("Applying limit: {limit_name}"), 8);
            let mut state = self.sim.slots[self.slot].limiter.take().unwrap_or_default();
            if state.active[upper].len() != totsize {
                state.active[upper] = vec![false; totsize];
            }
            let mut done = vec![false; totsize];
            let (mut added, mut removed) = (0usize, 0usize);
            for t in 0..nb {
                let e = ElemId(t as u32);
                let body = match self.cutfem_body(e) {
                    Ok(b) => b,
                    Err(_) => self.sim.model.body_of_tag(self.mesh.element(e).tag).map_or(0, |b| b as i32),
                };
                if body <= 0 {
                    continue;
                }
                let Some(bf) = self.sim.model.body_params(body as usize).and_then(|b| b.body_force).cloned() else { continue };
                let Some(lim) = self.real_at_nodes(&bf, &limit_name, e)? else { continue };
                let nodes = self.mesh.element(e).nodes.to_vec();
                for (i, &node) in nodes.iter().enumerate() {
                    let Some(ind) = self.sim.vars.by_id(id).slot(node) else { continue };
                    if done[ind] {
                        continue;
                    }
                    let act = &mut state.active[upper];
                    let val = self.sim.vars.by_id(id).values[ind];
                    let over = |v: f64| if upper == 0 { v < lim[i] - val_eps } else { v > lim[i] + val_eps };
                    if first_time {
                        if over(val) {
                            added += 1;
                            act[ind] = true;
                        } else {
                            act[ind] = false;
                        }
                    } else if freeze {
                    } else if act[ind] {
                        let load = *loads.get(ind).ok_or_else(|| {
                            SolverError::fatal(CALLER, "a dof beyond the Contact Load field (out-of-bounds read in the Fortran)")
                        })?;
                        if limit_sign * load > limit_sign * load_eps {
                            removed += 1;
                            act[ind] = false;
                        }
                    } else if over(val) {
                        added += 1;
                        act[ind] = true;
                    }
                    if act[ind] {
                        fixed.push((ind as u32, lim[i]));
                        let v = &mut self.sim.vars.by_id_mut(id).values[ind];
                        *v = if upper == 0 { v.max(lim[i]) } else { v.min(lim[i]) };
                    }
                    done[ind] = true;
                }
            }
            let count = state.active[upper].iter().filter(|&&a| a).count();
            let which = if upper == 0 { "lower" } else { "upper" };
            self.info(CALLER, &format!("Number of {which} limited dofs for {name}: {count}"), 5);
            self.info(CALLER, &format!("Added {added} and removed {removed} dofs in contact set"), 6);
            self.sim.slots[self.slot].limiter = Some(state);
        }
        Ok(fixed)
    }
}
