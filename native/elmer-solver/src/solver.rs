//! The `Solver` trait every physics crate implements, the registry that maps
//! a `Procedure = "Lib" "Entry"` to one, and [`SolverCtx`], the view a solver
//! gets of the running simulation (Elmer's `Solver_t` plus DefUtils.F90).

use std::collections::BTreeMap;

use elmer_core::{
    normalize_key, Communicator, ElemId, FieldRef, LocalSystem, Logger, MeshView, ParamMap, Params, Value, NO_DOF,
};

use crate::backend::{LinearSolverChoice, Perm, TimeTerms};
use crate::error::SolverError;
use crate::model::BodyParams;
use crate::sim::{LinearSystem, Simulation};
use crate::variable::{FieldMut, Strided, StridedMut};

/// One `Solver n` section's procedure.
///
/// `solve` is one visit of the solver in the coupled-system loop: Elmer's
/// `CALL Proc(Model, Solver, dt, Transient)` from `SingleSolver`
/// (MainUtils.F90:5235). A typical implementation loops
/// `Nonlinear System Max Iterations` times over
/// `default_initialize` / `default_update` / `default_dirichlet` /
/// `default_solve` and stops early when `default_converged` is true.
pub trait Solver: Send {
    /// Once, after the solver's variable exists (`<Proc>_init`).
    fn init(&mut self, _ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError>;
}

/// Builds a solver from its section.
pub type SolverFactory = fn(&ParamMap) -> Result<Box<dyn Solver>, SolverError>;

/// `(library, entry)` -> factory. Names are matched case-insensitively and
/// a library suffix such as `.so` is ignored.
#[derive(Default)]
pub struct SolverRegistry {
    map: BTreeMap<(String, String), SolverFactory>,
}

fn proc_key(library: &str, entry: &str) -> (String, String) {
    let lib = normalize_key(library);
    let lib = lib.strip_suffix(".so").or_else(|| lib.strip_suffix(".dll")).unwrap_or(&lib);
    let lib = lib.rsplit('/').next().unwrap_or(lib).to_string();
    (lib, normalize_key(entry))
}

impl SolverRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Registers e.g. `("HeatSolve", "HeatSolver")`.
    pub fn register(&mut self, library: &str, entry: &str, f: SolverFactory) {
        self.map.insert(proc_key(library, entry), f);
    }

    pub fn get(&self, library: &str, entry: &str) -> Option<SolverFactory> {
        self.map.get(&proc_key(library, entry)).copied()
    }
}

/// Reads a `Procedure` keyword: a `Dependent` procedure value, or a string
/// holding the two (optionally quoted) names.
pub fn procedure_of(params: &ParamMap) -> Option<(String, String)> {
    match params.get("procedure")? {
        Value::Dependent(d) => match &d.kind {
            elmer_core::DependentKind::Procedure { library, entry } => Some((library.clone(), entry.clone())),
            _ => None,
        },
        Value::Str(s) | Value::File(s) => {
            let parts: Vec<&str> = s
                .split('"')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .flat_map(|t| t.split_whitespace())
                .collect();
            match parts.as_slice() {
                [lib, entry] => Some((lib.to_string(), entry.to_string())),
                _ => None,
            }
        }
        _ => None,
    }
}

/// What a running solver sees: the simulation (fields, other sections),
/// the mesh, the message sink and the communicator, focused on one solver
/// slot.
pub struct SolverCtx<'a> {
    pub sim: &'a mut Simulation,
    pub mesh: &'a dyn MeshView,
    pub log: &'a dyn Logger,
    pub comm: &'a dyn Communicator,
    pub(crate) slot: usize,
}

impl<'a> SolverCtx<'a> {
    pub fn new(
        sim: &'a mut Simulation,
        slot: usize,
        mesh: &'a dyn MeshView,
        log: &'a dyn Logger,
        comm: &'a dyn Communicator,
    ) -> Self {
        Self { sim, mesh, log, comm, slot }
    }

    /// 1-based Elmer solver id.
    pub fn solver_id(&self) -> usize {
        self.sim.slots[self.slot].id
    }

    /// The `Solver n` section.
    pub fn params(&self) -> &ParamMap {
        &self.sim.slots[self.slot].params
    }

    pub fn simulation_params(&self) -> &ParamMap {
        &self.sim.model.simulation
    }

    pub fn constants(&self) -> &ParamMap {
        &self.sim.model.constants
    }

    /// The first value of a global variable (`time`, `timestep size`, ...).
    pub fn scalar(&self, name: &str) -> Option<f64> {
        self.sim.vars.scalar(name)
    }

    /// The solver's primary variable name, if it has one.
    pub fn variable_name(&self) -> Option<&str> {
        let id = self.sim.slots[self.slot].variable?;
        Some(&self.sim.vars.by_id(id).name)
    }

    pub fn primary(&self) -> Option<FieldRef<'_>> {
        let id = self.sim.slots[self.slot].variable?;
        Some(self.sim.vars.by_id(id).field_ref())
    }

    pub fn primary_mut(&mut self) -> Option<FieldMut<'_>> {
        let id = self.sim.slots[self.slot].variable?;
        Some(self.sim.vars.by_id_mut(id).field_mut())
    }

    /// Any primary variable by name.
    pub fn variable(&self, name: &str) -> Option<FieldRef<'_>> {
        self.sim.vars.get(name).map(|v| v.field_ref())
    }

    pub fn variable_mut(&mut self, name: &str) -> Option<FieldMut<'_>> {
        self.sim.vars.get_mut(name).map(|v| v.field_mut())
    }

    /// A primary variable or a component alias (`"displacement 2"`).
    pub fn component(&self, name: &str) -> Option<Strided<'_>> {
        self.sim.vars.strided(name)
    }

    pub fn component_mut(&mut self, name: &str) -> Option<StridedMut<'_>> {
        self.sim.vars.strided_mut(name)
    }

    /// Bulk elements the solver is active on, in mesh order
    /// (`Solver % ActiveElements`).
    pub fn active_elements(&self) -> impl Iterator<Item = ElemId> + '_ {
        self.sim.slots[self.slot].active.iter().copied()
    }

    /// Boundary elements that carry a boundary condition, in mesh order.
    pub fn boundary_elements(&self) -> impl Iterator<Item = ElemId> + '_ {
        let nb = self.mesh.n_bulk() as u32;
        (nb..nb + self.mesh.n_boundary() as u32)
            .map(ElemId)
            .filter(|&e| self.bc_id(e).is_some())
    }

    /// Body, equation, material, body force and initial condition of a
    /// bulk element. `None` if the element's tag maps to no body.
    pub fn body_params(&self, e: ElemId) -> Option<BodyParams<'_>> {
        let tag = self.mesh.element(e).tag;
        let b = self.sim.model.body_of_tag(tag)?;
        self.sim.model.body_params(b)
    }

    /// 1-based boundary condition id of a boundary element.
    pub fn bc_id(&self, e: ElemId) -> Option<usize> {
        if (e.0 as usize) < self.mesh.n_bulk() {
            return None;
        }
        self.sim.model.bc_of_tag(self.mesh.element(e).tag)
    }

    pub fn bc_params(&self, e: ElemId) -> Option<&ParamMap> {
        self.bc_id(e).map(|i| &self.sim.model.boundary_conditions[i - 1])
    }

    /// Global equation indices of an element's dofs: `slot * dofs + k` for
    /// each node, [`NO_DOF`] where the node is not in the system (Elmer
    /// GetElementDOFs + Perm).
    pub fn element_dofs(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        let Some(id) = self.sim.slots[self.slot].variable else {
            return;
        };
        let v = self.sim.vars.by_id(id);
        let dofs = v.dofs as u32;
        for &n in self.mesh.element(e).nodes {
            match v.slot(n) {
                Some(s) => out.extend((0..dofs).map(|k| s as u32 * dofs + k)),
                None => out.extend(std::iter::repeat_n(NO_DOF, dofs as usize)),
            }
        }
    }

    /// Perm slot of each element node, [`NO_DOF`] where the node is not in
    /// the system (Elmer `Perm(Indexes)`), as the assembly glue takes it.
    pub fn element_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        let Some(id) = self.sim.slots[self.slot].variable else {
            return;
        };
        let v = self.sim.vars.by_id(id);
        out.extend(
            self.mesh
                .element(e)
                .nodes
                .iter()
                .map(|&n| v.slot(n).map_or(NO_DOF, |s| s as u32)),
        );
    }

    /// The simulation's fields as `elmer-expr` reads them.
    pub fn fields(&self) -> crate::expr::SimFields<'_> {
        crate::expr::SimFields {
            vars: &self.sim.vars,
            mesh: self.mesh,
        }
    }

    /// `ListGetReal(section, key, n, Indexes)`: a scalar keyword at the
    /// nodes of element `e`, `None` when the section lacks it.
    pub fn real_at_nodes(&self, section: &ParamMap, key: &str, e: ElemId) -> Result<Option<Vec<f64>>, SolverError> {
        let Some(v) = section.get(key) else { return Ok(None) };
        let nodes = self.mesh.element(e).nodes;
        crate::expr::eval_at_nodes(v, nodes, &self.fields(), self.sim.backends.matc.as_deref()).map(Some)
    }

    pub fn info(&self, caller: &str, msg: &str, level: u8) {
        self.log.info(caller, msg, level);
    }

    pub fn warn(&self, caller: &str, msg: &str) {
        self.log.warn(caller, msg);
    }

    /// `DefaultInitialize`: create the matrix on first use, then zero the
    /// matrix values and the right-hand side.
    pub fn default_initialize(&mut self) -> Result<(), SolverError> {
        let sim = &mut *self.sim;
        let slot = &mut sim.slots[self.slot];
        if slot.system.is_none() {
            let backend = sim
                .backends
                .assembly
                .as_ref()
                .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
            let id = slot
                .variable
                .ok_or_else(|| SolverError::fatal("DefaultInitialize", "solver has no variable"))?;
            let perm = slot
                .perm
                .clone()
                .ok_or_else(|| SolverError::fatal("DefaultInitialize", "solver variable has no permutation"))?;
            let v = sim.vars.by_id(id);
            let perm = Perm {
                perm: perm.to_vec(),
                n_slots: v.values.len() / v.dofs,
            };
            let matrix = backend.create_matrix(self.mesh, &slot.active, &perm, v.dofs);
            let n = matrix.nrows;
            slot.system = Some(LinearSystem { matrix, rhs: vec![0.0; n] });
        }
        let sys = slot.system.as_mut().expect("system created above");
        sys.matrix.values.iter_mut().for_each(|v| *v = 0.0);
        sys.rhs.iter_mut().for_each(|v| *v = 0.0);
        Ok(())
    }

    /// `Default1stOrderTime`: fold `local.mass` into the local system with
    /// the solver's time scheme, using this slot's committed histories at
    /// the element's dofs.
    pub fn default_first_order_time(&mut self, local: &mut LocalSystem, e: ElemId, dt: f64) -> Result<(), SolverError> {
        let mut idx = Vec::new();
        self.element_dofs(e, &mut idx);
        let id = self.sim.slots[self.slot]
            .variable
            .ok_or_else(|| SolverError::fatal("Default1stOrderTime", "solver has no variable"))?;
        let v = self.sim.vars.by_id(id);
        let prev: Vec<Vec<f64>> = v
            .prev_values
            .iter()
            .map(|h| idx.iter().map(|&i| if i == NO_DOF { 0.0 } else { h[i as usize] }).collect())
            .collect();
        let cur: Vec<f64> = idx
            .iter()
            .map(|&i| if i == NO_DOF { 0.0 } else { v.values[i as usize] })
            .collect();
        let backend = self
            .sim
            .backends
            .time
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("time backend (elmer-timeint)".into()))?;
        let slot = &self.sim.slots[self.slot];
        let terms = TimeTerms {
            dt,
            params: &slot.params,
            simulation: &self.sim.model.simulation,
            prev: &prev,
            cur: &cur,
            done_time: slot.done_time,
            order: slot.order,
            dts: &self.sim.prev_sizes,
        };
        backend.add_time_terms(local, &terms)
    }

    /// `DefaultUpdateEquations`: add a local system at element `e`'s dofs.
    pub fn default_update(&mut self, local: &LocalSystem, e: ElemId) -> Result<(), SolverError> {
        let mut idx = Vec::new();
        self.element_slots(e, &mut idx);
        let sim = &mut *self.sim;
        let backend = sim
            .backends
            .assembly
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
        let sys = sim.slots[self.slot]
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal("DefaultUpdateEquations", "call default_initialize first"))?;
        backend.add_local(&mut sys.matrix, &mut sys.rhs, &idx, local);
        Ok(())
    }

    /// The Dirichlet values this solver's variable gets from boundary
    /// conditions: for each boundary element whose BC sets the variable (or
    /// one of its `"<name> k"` components), the value at each node,
    /// evaluated with `elmer-expr` (`ListGetReal(BC, name, n, Indexes)`),
    /// in boundary-element order (a later element overwrites an earlier
    /// one). `<name> Condition` masks are refused.
    pub fn dirichlet_values(&self) -> Result<Vec<(u32, f64)>, SolverError> {
        let Some(id) = self.sim.slots[self.slot].variable else {
            return Ok(Vec::new());
        };
        let v = self.sim.vars.by_id(id);
        let keys: Vec<(usize, String)> = if v.dofs == 1 {
            vec![(0, v.name.clone())]
        } else {
            (0..v.dofs).map(|k| (k, format!("{} {}", v.name, k + 1))).collect()
        };
        let mut fixed: BTreeMap<u32, f64> = BTreeMap::new();
        for e in self.boundary_elements() {
            let bc = self.bc_params(e).expect("boundary_elements filters on BC");
            for (k, key) in &keys {
                if bc.get(&format!("{key} condition")).is_some() {
                    return Err(SolverError::Unsupported(format!("Dirichlet mask > {key} Condition <")));
                }
                let Some(vals) = self.real_at_nodes(bc, key, e)? else { continue };
                for (&n, &x) in self.mesh.element(e).nodes.iter().zip(&vals) {
                    if let Some(s) = v.slot(n) {
                        fixed.insert((s * v.dofs + k) as u32, x);
                    }
                }
            }
        }
        Ok(fixed.into_iter().collect())
    }

    /// `DefaultDirichletBCs` + `EnforceDirichletConditions`.
    pub fn default_dirichlet(&mut self) -> Result<(), SolverError> {
        let fixed = self.dirichlet_values()?;
        let sim = &mut *self.sim;
        let backend = sim
            .backends
            .assembly
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
        let slot = &mut sim.slots[self.slot];
        let sys = slot
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal("DefaultDirichletBCs", "call default_initialize first"))?;
        backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, &fixed, &slot.params)
    }

    /// `DefaultSolve`: solve the assembled system into the primary variable
    /// with the linear solver the section selects, then run the nonlinear
    /// `ComputeChange`. Returns the variable norm.
    pub fn default_solve(&mut self) -> Result<f64, SolverError> {
        let (mesh, log, comm, slot_ix) = (self.mesh, self.log, self.comm, self.slot);
        let sim = &mut *self.sim;
        let slot = &mut sim.slots[slot_ix];
        let id = slot
            .variable
            .ok_or_else(|| SolverError::fatal("DefaultSolve", "solver has no variable"))?;
        let sys = slot
            .system
            .as_ref()
            .ok_or_else(|| SolverError::fatal("DefaultSolve", "no assembled system"))?;
        let p = &slot.params;
        // SolveCore.F90:3158-3186: keep the previous iterate when relaxing or
        // measuring the solution change.
        let need_prev = p.real("nonlinear system relaxation factor").is_some_and(|r| r != 1.0)
            || matches!(
                p.string("nonlinear system convergence measure"),
                Some("residual") | Some("solution")
            );
        let var = sim.vars.by_id_mut(id);
        let n = sys.matrix.nrows;
        if need_prev {
            var.nonlin_values = Some(var.values[..n].to_vec());
        }
        match sim.backends.linear.select(p)? {
            LinearSolverChoice::Direct(mut d) => {
                d.factor(&sys.matrix)?;
                let mut x = sys.rhs.clone();
                d.solve_in_place(&mut x)?;
                var.values[..n].copy_from_slice(&x);
            }
            LinearSolverChoice::Iterative { mut solver, precond } => {
                let m: Box<dyn elmer_core::Preconditioner<f64> + Send> = match precond {
                    Some(f) => f(p, &sys.matrix)?,
                    None => Box::new(elmer_core::IdentityPrecond),
                };
                solver.solve(&sys.matrix, &*m, &sys.rhs, &mut var.values[..n], comm)?;
            }
        }
        var.values_changed = true;
        sim.compute_change(slot_ix, false, Some(n), true, mesh, log, comm)?;
        Ok(sim.vars.by_id(id).conv.norm)
    }

    /// The nonlinear `ComputeChange` alone, for a solver that updates its
    /// primary variable itself instead of through `default_solve`. Returns
    /// the variable norm.
    pub fn nonlinear_change(&mut self) -> Result<f64, SolverError> {
        let (mesh, log, comm) = (self.mesh, self.log, self.comm);
        self.sim.compute_change(self.slot, false, None, false, mesh, log, comm)?;
        Ok(self.sim.slots[self.slot]
            .variable
            .map_or(0.0, |id| self.sim.vars.by_id(id).conv.norm))
    }

    /// `DefaultConverged`: the last nonlinear `ComputeChange` converged.
    pub fn default_converged(&self) -> bool {
        self.sim.slots[self.slot]
            .variable
            .is_some_and(|id| self.sim.vars.by_id(id).conv.nonlin_converged > 0)
    }

    /// The assembled system, if any.
    pub fn system(&self) -> Option<&LinearSystem> {
        self.sim.slots[self.slot].system.as_ref()
    }

    pub fn system_mut(&mut self) -> Option<&mut LinearSystem> {
        self.sim.slots[self.slot].system.as_mut()
    }
}
