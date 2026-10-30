//! Nodal loads: `<var> Load` boundary and body-force values added to the
//! right-hand side (`SetNodalLoads`, SolverBasics.F90:1822), the bulk
//! system snapshot (`DefaultFinishBulkAssembly`, DefUtils.F90:7756) and the
//! `<var> Loads` / `<var> Residual` fields computed from it after the solve
//! (`CalculateLoads`, SolveCore.F90:364, called at :1739 and :1769).

use elmer_core::{CsrMatrix, ElemId, LinearOperator, Logger, MeshView, ParamMap, Params};

use crate::error::SolverError;
use crate::sim::{Simulation, VarSpec};
use crate::variable::{VarId, Variable};

/// Per-assembly bookkeeping of one slot, reset by `DefaultInitialize`.
#[derive(Clone, Debug, Default)]
pub(crate) struct AssemblyState {
    /// `A % BulkValues`, `A % BulkRHS`.
    pub bulk: Option<(Vec<f64>, Vec<f64>)>,
    /// `DefaultFinishBulkAssembly` has run for this assembly.
    pub bulk_done: bool,
    /// `DefaultFinishBoundaryAssembly` has run for this assembly.
    pub boundary_done: bool,
    /// A boundary element added a nonzero contribution.
    pub boundary_touched: bool,
    /// The Dirichlet dofs of this assembly (`A % ConstrainedDOF`).
    pub constrained: Vec<(u32, f64)>,
}

/// Whether `DefaultFinishBulkAssembly` keeps the bulk system: `Calculate
/// Loads` (unless `Calculate Loads Slot` names another stage) or `Control
/// Use Loads` (DefUtils.F90:7797-7808). The other triggers (`Constant Bulk
/// System`, ...) are refused where they are read.
pub(crate) fn keeps_bulk(p: &ParamMap) -> Result<bool, SolverError> {
    if p.logical("control use loads") == Some(true) || p.logical("constraint modes analysis") == Some(true) {
        return Ok(true);
    }
    if p.logical("calculate loads") != Some(true) {
        return Ok(false);
    }
    match p.string("calculate loads slot") {
        None => Ok(true),
        Some(s) if elmer_core::normalize_key(s) == "bulk assembly" => Ok(true),
        Some(s) => Err(SolverError::Unsupported(format!("Calculate Loads Slot = {s}"))),
    }
}

/// `A x - b` over the first `n` rows.
pub(crate) fn residual(a: &CsrMatrix<f64>, values: &[f64], x: &[f64], b: &[f64]) -> Vec<f64> {
    let m = CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: values.to_vec(),
        diag: None,
    };
    let mut r = vec![0.0; a.nrows];
    m.apply(x, &mut r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri -= bi;
    }
    r
}

impl Simulation {
    /// The `<var> Loads` and `<var> Residual` fields of a solver that asks
    /// for them (`Calculate Loads`, `Calculate Residual`; `Calculate
    /// Boundary Fluxes` implies loads), with the solver variable's perm
    /// (AddEquationBasics, MainUtils.F90:2342-2390).
    pub(crate) fn add_load_variables(&mut self, params: &mut ParamMap, spec: &VarSpec, primary: VarId) -> Result<(), SolverError> {
        // MainUtils.F90:2345: fluxes are summed from the nodal loads.
        if params.logical("calculate boundary fluxes") == Some(true) {
            if params.logical("calculate boundary weights") == Some(true) {
                return Err(SolverError::Unsupported("Calculate Boundary Weights".into()));
            }
            params.insert("Calculate Loads", elmer_core::Value::Logical(true));
        }
        for what in ["loads", "residual"] {
            if params.logical(&format!("calculate {what}")) != Some(true) {
                continue;
            }
            let name = format!("{} {what}", spec.name);
            if self.vars.id(&name).is_some() {
                continue;
            }
            let p = self.vars.by_id(primary);
            let mut v = Variable::new(&name, p.kind, p.dofs, p.perm.clone(), vec![0.0; p.values.len()]);
            v.output = params.logical(&format!("save {what}")).unwrap_or(p.output);
            self.vars.add(v)?;
        }
        Ok(())
    }

    /// `SetNodalLoads` for slot `ix`'s variable into its right-hand side:
    /// each node gets the first boundary condition's `<name> Load` that
    /// reaches it, then (separately) the first body force's.
    pub(crate) fn set_nodal_loads(&mut self, ix: usize, mesh: &dyn MeshView) -> Result<(), SolverError> {
        let Some(id) = self.slots[ix].variable else { return Ok(()) };
        let (name, dofs) = {
            let v = self.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        let mut rhs_add: Vec<(usize, f64)> = Vec::new();
        for dof in 0..dofs {
            let comp = if dofs == 1 { name.clone() } else { format!("{name} {}", dof + 1) };
            let key = format!("{comp} load");
            let point = |p: &ParamMap| p.get("target nodes").is_some() || p.get("target coordinates").is_some();
            for p in self.model.boundary_conditions.iter().chain(&self.model.body_forces) {
                if p.get(&format!("{key} dofs")).is_some() || (p.get(&key).is_some() && point(p)) {
                    return Err(SolverError::Unsupported(format!("{key} on target nodes or as DOFs")));
                }
            }
            let nslots = self.vars.by_id(id).values.len() / dofs;
            // Boundary conditions (one DoneLoad over all of them), then body forces.
            for bulk in [false, true] {
                let mut done = vec![false; nslots];
                let elements = if bulk {
                    0..mesh.n_bulk() as u32
                } else {
                    mesh.n_bulk() as u32..(mesh.n_bulk() + mesh.n_boundary()) as u32
                };
                let n_sections = if bulk { self.model.body_forces.len() } else { self.model.boundary_conditions.len() };
                for sec in 0..n_sections {
                    let section = if bulk { &self.model.body_forces[sec] } else { &self.model.boundary_conditions[sec] };
                    let Some(val) = section.get(&key) else { continue };
                    if !bulk && section.get("target boundaries").is_none() && section.get("name").is_none() {
                        continue;
                    }
                    for e in elements.clone().map(ElemId) {
                        let el = mesh.element(e);
                        let owner = if bulk {
                            self.model
                                .body_of_tag(el.tag)
                                .and_then(|b| self.model.bodies[b - 1].integer("body force"))
                                .map(|b| b as usize)
                        } else {
                            self.model.bc_of_tag(el.tag)
                        };
                        if owner != Some(sec + 1) {
                            continue;
                        }
                        let w = crate::expr::eval_at_nodes_with(
                            val,
                            el.nodes,
                            &self.fields(mesh),
                            self.backends.matc.as_deref(),
                            &self.backends.udfs,
                        )?;
                        let v = self.vars.by_id(id);
                        for (&n, &x) in el.nodes.iter().zip(&w) {
                            let Some(k) = v.slot(n) else { continue };
                            if done[k] {
                                continue;
                            }
                            done[k] = true;
                            rhs_add.push((dofs * k + dof, x));
                        }
                    }
                }
            }
        }
        if let Some(sys) = self.slots[ix].system.as_mut() {
            for (k, x) in rhs_add {
                if let Some(b) = sys.rhs.get_mut(k) {
                    *b += x;
                }
            }
        }
        Ok(())
    }

    /// `CalculateLoads` into `<var> Loads` from the kept bulk system (after
    /// the solve, before the nonlinear change) or into `<var> Residual`
    /// from the solved system (after it), over the first `n` dofs.
    pub(crate) fn calculate_loads(
        &mut self,
        ix: usize,
        n: usize,
        from_bulk: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
    ) -> Result<(), SolverError> {
        let slot = &self.slots[ix];
        let Some(id) = slot.variable else { return Ok(()) };
        let what = if from_bulk { "loads" } else { "residual" };
        let name = format!("{} {what}", self.vars.by_id(id).name);
        let Some(lid) = self.vars.id(&name) else { return Ok(()) };
        if slot.params.logical(&format!("calculate {what}")) == Some(false) {
            return Ok(());
        }
        for k in ["calculate energy inner product", "calculate energy norm"] {
            if slot.params.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(k.into()));
            }
        }
        let Some(sys) = slot.system.as_ref() else { return Ok(()) };
        let x = &self.vars.by_id(id).values[..n];
        let r = if from_bulk {
            let (values, rhs) = slot
                .asm
                .bulk
                .as_ref()
                .ok_or_else(|| SolverError::fatal("CalculateLoads", "\"BulkValues\" are not associated!"))?;
            residual(&sys.matrix, values, x, &rhs[..n])
        } else {
            residual(&sys.matrix, &sys.matrix.values, x, &sys.rhs[..n])
        };
        let loads = self.vars.by_id_mut(lid);
        let m = loads.values.len().min(r.len());
        loads.values[..m].copy_from_slice(&r[..m]);
        loads.values_changed = true;
        if from_bulk && self.slots[ix].params.logical("calculate boundary fluxes") == Some(true) {
            self.boundary_fluxes(id, lid, mesh, log)?;
        }
        Ok(())
    }

    /// The `Calculate Boundary Fluxes` part of `CalculateLoads`
    /// (SolveCore.F90:536-640): each boundary condition's sum of the nodal
    /// loads of its boundary nodes, a node counted for the first boundary
    /// that reaches it (`Calculate Boundary Flux` on some boundary limits the
    /// sums to those), stored as the Simulation entry `res: <var> Flux over
    /// BC i` (`<var> j Flux ...` per component) for SaveScalars.
    fn boundary_fluxes(&mut self, id: VarId, lid: VarId, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        let (name, dofs) = {
            let v = self.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        log.info("CalculateLoads", "Computing boundary fluxes from nodal loads", 6);
        let nbc = self.model.boundary_conditions.len();
        let only_given = self.model.boundary_conditions.iter().any(|b| b.get("calculate boundary flux").is_some());
        let loads = self.vars.by_id(lid);
        let nslots = loads.values.len() / dofs.max(1);
        let mut sum = vec![0.0; nbc * dofs];
        let mut active = vec![false; nbc];
        let mut shared = vec![false; nbc];
        let mut summed = vec![0usize; nslots];
        let nb = mesh.n_bulk() as u32;
        for e in (nb..nb + mesh.n_boundary() as u32).map(ElemId) {
            let el = mesh.element(e);
            let Some(bc) = self.model.bc_of_tag(el.tag) else { continue };
            if only_given && self.model.boundary_conditions[bc - 1].logical("calculate boundary flux") != Some(true) {
                continue;
            }
            for &n in el.nodes {
                let Some(ind) = loads.slot(n) else { continue };
                if summed[ind] == 0 {
                    for j in 0..dofs {
                        sum[dofs * (bc - 1) + j] += loads.values[dofs * ind + j];
                    }
                    summed[ind] = bc;
                    active[bc - 1] = true;
                } else if summed[ind] != bc {
                    shared[bc - 1] = true;
                    shared[summed[ind] - 1] = true;
                }
            }
        }
        let mut n_active = 0;
        for i in 0..nbc {
            if !active[i] {
                continue;
            }
            if shared[i] {
                let msg = format!("Boundary {} includes inseparable dofs!", i + 1);
                if self.flux_shared_warned {
                    log.info("CalculateLoads", &msg, 8);
                } else {
                    log.warn("CalculateLoads", &msg);
                }
            }
            n_active += 1;
            for j in 0..dofs {
                let what = if dofs == 1 {
                    format!("{name} Flux over BC {}", i + 1)
                } else {
                    format!("{name} {} Flux over BC {}", j + 1, i + 1)
                };
                let v = sum[dofs * i + j];
                self.model.simulation.insert(&format!("res: {what}"), elmer_core::Value::real(v));
                log.info("CalculateLoads", &format!("{what}: {}", crate::fortran_fmt::es(v, 12, 5)), 6);
            }
        }
        self.flux_shared_warned = true;
        if n_active > 1 {
            for j in 0..dofs {
                let what = if dofs == 1 {
                    format!("{name} Flux over all BCs")
                } else {
                    format!("{name} {} Flux over all BCs", j + 1)
                };
                let total: f64 = sum.iter().skip(j).step_by(dofs).sum();
                log.info("CalculateLoads", &format!("{what}: {}", crate::fortran_fmt::es(total, 12, 5)), 6);
            }
        }
        Ok(())
    }
}

/// Whether a local system adds anything (a zero boundary contribution
/// leaves the bulk system as it was).
pub(crate) fn is_nonzero(local: &elmer_core::LocalSystem) -> bool {
    local.stiff.iter().chain(&local.force).any(|&v| v != 0.0)
}
