//! The `Solver` trait every physics crate implements, the registry that maps
//! a `Procedure = "Lib" "Entry"` to one, and [`SolverCtx`], the view a solver
//! gets of the running simulation (Elmer's `Solver_t` plus DefUtils.F90).

use std::collections::BTreeMap;
use std::sync::Arc;

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

/// An element assembly routine named by `Bulk Assembly Procedure` or
/// `Boundary Assembly Procedure` (Elmer's 11-argument `BulkAssembly(Model,
/// Solver, dt, Transient, Mass, Damp, Stiff, Force, Element, Nrow, Ncol)`,
/// called through `ExecLocalAssembly`).
pub trait ElementAssembly: Send {
    /// Adds element `e`'s contributions to `local`, which arrives zeroed and
    /// sized for the element's dofs with `mass` and `damp` present.
    /// `ctx` is the solver the routine assembles for
    /// (`CurrentModel % CurrentElement` is `e`).
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        dt: f64,
        transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError>;
}

/// Builds an assembly routine from its solver section. A closure, so the
/// composition can hand it the basis and quadrature it needs.
pub type AssemblyFactory = Arc<dyn Fn(&ParamMap) -> Result<Box<dyn ElementAssembly>, SolverError> + Send + Sync>;

/// `(library, entry)` -> factory, for solver procedures and for element
/// assembly routines. Names are matched case-insensitively and a library
/// suffix such as `.so` and any directory part are ignored.
#[derive(Default)]
pub struct SolverRegistry {
    map: BTreeMap<(String, String), SolverFactory>,
    assembly: BTreeMap<(String, String), AssemblyFactory>,
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

    /// Registers an element assembly routine, e.g.
    /// `("PoissonAssembly", "BulkAssembly")`.
    pub fn register_assembly(&mut self, library: &str, entry: &str, f: AssemblyFactory) {
        self.assembly.insert(proc_key(library, entry), f);
    }

    pub fn assembly(&self, library: &str, entry: &str) -> Option<AssemblyFactory> {
        self.assembly.get(&proc_key(library, entry)).cloned()
    }
}

/// Reads a `Procedure` keyword: a `Dependent` procedure value, or a string
/// holding the two (optionally quoted) names.
pub fn procedure_of(params: &ParamMap) -> Option<(String, String)> {
    procedure_value(params, "procedure")
}

/// [`procedure_of`] for any procedure-valued keyword (`Bulk Assembly
/// Procedure`, ...).
pub fn procedure_value(params: &ParamMap, key: &str) -> Option<(String, String)> {
    match params.get(key)? {
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
        let b = &self.sim.backends;
        crate::expr::eval_at_nodes_with(v, nodes, &self.fields(), b.matc.as_deref(), &b.udfs).map(Some)
    }

    /// `ListGetCReal(section, key)`: a scalar keyword at node 1, so it may
    /// depend on `time` or a field (MATC, tables, procedures). `None` when
    /// the section lacks it.
    pub fn creal(&self, section: &ParamMap, key: &str) -> Result<Option<f64>, SolverError> {
        let b = &self.sim.backends;
        crate::expr::eval_creal(section, key, &self.fields(), b.matc.as_deref(), &b.udfs)
    }

    /// Publishes the current field values to the [`crate::LiveModel`], for a
    /// solver that changed values itself before calling a property that reads
    /// them through a user function. The `default_*` steps do it already.
    pub fn sync_model(&self) {
        self.sim.publish_live();
    }

    pub fn info(&self, caller: &str, msg: &str, level: u8) {
        self.log.info(caller, msg, level);
    }

    pub fn warn(&self, caller: &str, msg: &str) {
        self.log.warn(caller, msg);
    }

    /// `DefaultInitialize`: run the `Slave Solvers` and `Nonlinear Pre
    /// Solvers`, create the matrix on first use, then zero the matrix values
    /// and the right-hand side.
    pub fn default_initialize(&mut self) -> Result<(), SolverError> {
        let (mesh, log, comm) = (self.mesh, self.log, self.comm);
        self.sim.run_slave_solvers(self.slot, "Slave Solvers", mesh, log, comm)?;
        self.sim.run_slave_solvers(self.slot, "Nonlinear Pre Solvers", mesh, log, comm)?;
        self.sim.publish_live();
        self.sim.slots[self.slot].asm = Default::default();
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
        self.time_terms(local, e, dt, false)
    }

    fn time_terms(&mut self, local: &mut LocalSystem, e: ElemId, dt: f64, second: bool) -> Result<(), SolverError> {
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
        if second {
            backend.add_2nd_order_time_terms(local, &terms)
        } else {
            backend.add_time_terms(local, &terms)
        }
    }

    /// `Default2ndOrderTime`: fold `local.mass` and `local.damp` into the
    /// local system with the solver's second-order scheme (Bossak or
    /// generalized-alpha), using the seven `PrevValues` columns at the
    /// element's dofs.
    pub fn default_second_order_time(&mut self, local: &mut LocalSystem, e: ElemId, dt: f64) -> Result<(), SolverError> {
        self.time_terms(local, e, dt, true)
    }

    /// `ListAddInteger(Params, 'Time Derivative Order', k)` from a solver's
    /// `init`, as Elmer's `_Init` routines do (ElasticSolver_Init sets 2).
    /// Those run before the variable's histories are allocated and `init`
    /// runs after, so in a transient run this also sets `Solver % TimeOrder`
    /// and reallocates the `PrevValues` columns (seven for order 2).
    pub fn set_time_derivative_order(&mut self, k: u8) -> Result<(), SolverError> {
        let k = k.clamp(1, 2);
        let slot = &mut self.sim.slots[self.slot];
        if slot.done_time > 0 {
            return Err(SolverError::fatal(
                "SetTimeDerivativeOrder",
                "Time Derivative Order changed after time stepping started",
            ));
        }
        slot.params.insert("Time Derivative Order", Value::Integer(vec![i32::from(k)]));
        if slot.time_order == 0 {
            return Ok(());
        }
        slot.time_order = k;
        let m = if k == 2 { 7 } else { slot.order };
        if let Some(id) = slot.variable {
            let v = self.sim.vars.by_id_mut(id);
            v.prev_values = vec![vec![0.0; v.values.len()]; m];
        }
        Ok(())
    }

    /// The solver variable's convergence record (`Norm`, `NonlinChange`, ...).
    pub fn convergence(&self) -> Option<&crate::variable::Convergence> {
        let id = self.sim.slots[self.slot].variable?;
        Some(&self.sim.vars.by_id(id).conv)
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
        if e.0 as usize >= self.mesh.n_bulk() && crate::loads::is_nonzero(local) {
            sim.slots[self.slot].asm.boundary_touched = true;
        }
        Ok(())
    }

    /// `DefaultFinishBulkAssembly`: keeps the bulk system for `Calculate
    /// Loads`. Call it between the bulk and the boundary element loops; when
    /// a solver does not, [`default_dirichlet`](Self::default_dirichlet)
    /// keeps it there instead, which is the same system unless boundary
    /// elements added something (then it refuses).
    pub fn default_finish_bulk_assembly(&mut self) -> Result<(), SolverError> {
        let slot = &mut self.sim.slots[self.slot];
        if !slot.asm.bulk_done && crate::loads::keeps_bulk(&slot.params)? {
            if let Some(sys) = slot.system.as_ref() {
                slot.asm.bulk = Some((sys.matrix.values.clone(), sys.rhs.clone()));
            }
        }
        slot.asm.bulk_done = true;
        Ok(())
    }

    /// `DefaultFinishBoundaryAssembly`: adds the `<var> Load` nodal loads of
    /// boundary conditions and body forces to the right-hand side
    /// (`SetNodalLoads`). [`default_dirichlet`](Self::default_dirichlet)
    /// does it when the solver has not.
    pub fn default_finish_boundary_assembly(&mut self) -> Result<(), SolverError> {
        if !self.sim.slots[self.slot].asm.bulk_done {
            if self.sim.slots[self.slot].asm.boundary_touched && crate::loads::keeps_bulk(self.params())? {
                return Err(SolverError::Unsupported(
                    "Calculate Loads after boundary contributions: call default_finish_bulk_assembly after the bulk loop"
                        .into(),
                ));
            }
            self.default_finish_bulk_assembly()?;
        }
        if !self.sim.slots[self.slot].asm.boundary_done {
            self.sim.set_nodal_loads(self.slot, self.mesh)?;
            self.sim.slots[self.slot].asm.boundary_done = true;
        }
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
        self.dirichlet_values_of(id, 0)
    }

    /// [`dirichlet_values`](Self::dirichlet_values) of any nodal variable,
    /// with `offset` added to every row (`DefaultDirichletBCs(Ux=Var,
    /// UOffset=offset)`).
    pub fn dirichlet_values_of(&self, id: crate::VarId, offset: usize) -> Result<Vec<(u32, f64)>, SolverError> {
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
                        fixed.insert((offset + s * v.dofs + k) as u32, x);
                    }
                }
            }
        }
        Ok(fixed.into_iter().collect())
    }

    /// `DefaultDirichletBCs` + `EnforceDirichletConditions`.
    pub fn default_dirichlet(&mut self) -> Result<(), SolverError> {
        self.default_finish_boundary_assembly()?;
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
        backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, &fixed, &slot.params)?;
        // `A % ConstrainedDOF`, which source control reads.
        slot.asm.constrained = fixed;
        Ok(())
    }

    /// `DefaultSolve` (DefUtils.F90:4470): returns early when `Linear System
    /// Solver Disabled`; with `Monolithic Slave` the system is merged with
    /// the slave solver's first (`MergeSlaveSolvers`, ContactUtils.F90:2619),
    /// solved and split again; otherwise see [`solve_system`](Self::solve_system).
    /// Returns the variable norm.
    pub fn default_solve(&mut self) -> Result<f64, SolverError> {
        crate::constraints::generate(self)?;
        let slot = &self.sim.slots[self.slot];
        if (slot.constraints.is_some() || slot.own_constraints.is_some())
            && ["apply source control", "monolithic slave"]
                .iter()
                .any(|k| slot.params.logical(k) == Some(true))
        {
            return Err(SolverError::Unsupported(
                "linear constraints with Apply Source Control or Monolithic Slave".into(),
            ));
        }
        if self.params().logical("linear system solver disabled") == Some(true) {
            self.info("DefaultSolve", "Solver disabled, exiting early!", 10);
            return Ok(self.convergence().map_or(0.0, |c| c.norm));
        }
        if self.params().logical("apply source control") == Some(true) {
            return crate::control::solve_controlled(self);
        }
        if self.params().logical("monolithic slave") == Some(true) {
            return crate::monolithic::solve_merged(self);
        }
        self.solve_system()
    }

    /// `SolveSystem` -> `SolveLinearSystem` (SolveCore.F90:1186): keep the
    /// previous iterate when needed, the zero right-hand-side shortcut,
    /// `Linear System Scaling` (default on) around the solve, the
    /// nullify-guess options, the solve with the linear solver the section
    /// selects, then the nonlinear `ComputeChange`. Returns the variable norm.
    pub(crate) fn solve_system(&mut self) -> Result<f64, SolverError> {
        let (mesh, log, comm, slot_ix) = (self.mesh, self.log, self.comm, self.slot);
        let sim = &mut *self.sim;
        let nonlin_iter = sim.vars.scalar("nonlin iter");
        // SolveCore.F90:3150-3163 (ListGetCReal).
        let relaxation = crate::expr::eval_creal(
            &sim.slots[slot_ix].params,
            "nonlinear system relaxation factor",
            &crate::expr::SimFields { vars: &sim.vars, mesh },
            sim.backends.matc.as_deref(),
            &sim.backends.udfs,
        )?;
        let slot = &mut sim.slots[slot_ix];
        let id = slot
            .variable
            .ok_or_else(|| SolverError::fatal("DefaultSolve", "solver has no variable"))?;
        let p = &slot.params;
        for k in [
            "nonlinear system acceleration",
            "nonlinear system compute change in scaled system",
            "linear system normalize guess",
            "apply limiter",
        ] {
            if p.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("SolveLinearSystem option {k}")));
            }
        }
        let sys = slot
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal("DefaultSolve", "no assembled system"))?;
        let var = sim.vars.by_id_mut(id);
        let n = sys.matrix.nrows;
        // SolveSystem's RestrictionMode: SolveWithLinearRestriction solves
        // the collection system [A C^T; C 0] instead.
        let restricted = match slot.constraints.as_ref().or(slot.own_constraints.as_ref()) {
            None => None,
            Some(c) => {
                crate::constraints::refuse_options(p)?;
                if var.values.len() != n {
                    return Err(SolverError::Unsupported(
                        "linear constraints on a variable longer than its system".into(),
                    ));
                }
                let (a, b) = crate::constraints::collection(&sys.matrix, &sys.rhs, c, &slot.asm.constrained, p, log);
                let mut x = var.values.clone();
                x.resize(a.nrows, 0.0);
                Some((a, b, x))
            }
        };
        let (mut ca, mut cb, mut cx) = match restricted {
            Some((a, b, x)) => (Some(a), Some(b), Some(x)),
            None => (None, None, None),
        };

        // SolveCore.F90:3150-3186: keep the previous iterate in residual
        // mode, when relaxing or when measuring the solution change.
        let residual_mode = p.logical("linear system residual mode") == Some(true);
        if residual_mode
            && matches!(
                p.string("nonlinear system convergence measure"),
                Some("residual") | Some("linear system residual")
            )
        {
            return Err(SolverError::Unsupported(
                "Linear System Residual Mode with a residual convergence measure".into(),
            ));
        }
        let need_prev = residual_mode
            || relaxation.is_some_and(|r| r != 1.0)
            || matches!(
                p.string("nonlinear system convergence measure"),
                Some("residual") | Some("solution")
            );
        if need_prev {
            var.nonlin_values = Some(var.values[..n].to_vec());
        }

        // SolveCore.F90:3177-3189: in residual mode solve A dx = b - A x0
        // from dx = 0; ComputeChange adds x0 back.
        let mut res = if residual_mode {
            let mut r = vec![0.0; n];
            elmer_core::LinearOperator::apply(&sys.matrix, &var.values[..n], &mut r);
            for (ri, bi) in r.iter_mut().zip(&sys.rhs[..n]) {
                *ri = bi - *ri;
            }
            var.values[..n].iter_mut().for_each(|x| *x = 0.0);
            Some(r)
        } else {
            None
        };
        let b: &mut [f64] = match (&mut res, &mut cb) {
            (Some(r), _) => r,
            (None, Some(bc)) => bc,
            (None, None) => &mut sys.rhs[..n],
        };

        // SolveCore.F90:1502-1527: b = 0 has only the trivial solution. In
        // residual mode that zeroes the solution too, as in Elmer: the
        // early return skips ComputeChange, which would add x0 back.
        let mut bb = [b.iter().fold(0.0, |s, b| s + b * b)];
        comm.sum_f64(&mut bb);
        let bnorm = bb[0].sqrt();
        if bnorm <= f64::MIN_POSITIVE && p.logical("skip zero rhs test") != Some(true) {
            log.info("SolveLinearSystem", "Solution trivially zero!", 5);
            var.values.iter_mut().for_each(|x| *x = 0.0);
            var.values_changed = true;
            if let Some(it) = nonlin_iter {
                var.conv.nonlin_iter = it.round() as i32;
                sim.vars.set_scalar("nonlin iter", it + 1.0);
            }
            let var = sim.vars.by_id_mut(id);
            var.conv.norm = 0.0;
            var.conv.nonlin_converged = 1;
            sim.update_2nd_order(slot_ix, n)?;
            sim.publish_live();
            sim.run_slave_solvers(slot_ix, "Nonlinear Post Solvers", mesh, log, comm)?;
            return Ok(0.0);
        }

        // SolveCore.F90:1254-1261: the interleaved real system of a complex
        // one is solved in complex arithmetic. EnforceDirichletConditions
        // scaled the Dirichlet rows by the real diagonal, where Elmer takes
        // the complex one's modulus; that changes their scaling, not x.
        let nullify = p.logical("linear system nullify guess") == Some(true)
            || (p.logical("linear system nullify first guess") == Some(true) && nonlin_iter.is_some_and(|i| i.round() <= 1.0));
        if let (Some(a), Some(x)) = (&mut ca, &mut cx) {
            log.info("SolveWithLinearRestriction", "Now solving the linear system with constraints!", 10);
            solve_linear(&sim.backends.linear, p, a, b, x, bnorm, nullify, log, comm)?;
            // ComputeChange over the collection solution, multipliers
            // included; then 'Picking solution from collection solution'.
            // Without NonlinValues the Fortran's x0 is the variable's own
            // (previous, n long) values.
            let nc = x.len();
            let own_x0 = var.nonlin_values.is_none();
            let prev = std::mem::replace(&mut var.values, std::mem::take(x));
            if own_x0 {
                var.nonlin_values = Some(prev);
            }
            var.values_changed = true;
            sim.compute_change(slot_ix, false, Some(nc), true, mesh, log, comm)?;
            log.info("SolveWithLinearRestriction", "Picking solution from collection solution", 10);
            let var = sim.vars.by_id_mut(id);
            var.values.truncate(n);
            if own_x0 {
                var.nonlin_values = None;
            }
            log.info("SolveWithLinearRestriction", "All done", 10);
        } else {
            let x = &mut var.values[..n];
            solve_linear(&sim.backends.linear, p, &mut sys.matrix, b, x, bnorm, nullify, log, comm)?;
            var.values_changed = true;
            // SolveCore.F90:1727-1741, then the change, then :1763-1771.
            sim.calculate_loads(slot_ix, n, true, mesh, log)?;
            sim.compute_change(slot_ix, false, Some(n), true, mesh, log, comm)?;
            sim.calculate_loads(slot_ix, n, false, mesh, log)?;
        }
        sim.update_2nd_order(slot_ix, n)?;
        sim.publish_live();
        sim.run_slave_solvers(slot_ix, "Nonlinear Post Solvers", mesh, log, comm)?;
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

/// `SolveLinearSystem`'s scaling, solve and back-scaling of `A x = b`
/// (SolveCore.F90:1254-1700): `Linear System Scaling` (with the complex
/// form for `Linear System Complex`), the nullify-guess flag, the solver
/// the section selects. `bnorm` is `|b|`, which turns the rhs scaling on.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_linear(
    linear: &crate::LinearSolverRegistry,
    p: &ParamMap,
    matrix: &mut elmer_core::CsrMatrix<f64>,
    b: &mut [f64],
    x: &mut [f64],
    bnorm: f64,
    nullify: bool,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    let complex = p.logical("linear system complex") == Some(true)
        && p.logical("linear system skip complex") != Some(true);
    if complex {
        log.info("SolveLinearSystem", "Assuming complex valued linear system", 6);
    }
    let scale = p.logical("linear system skip scaling") != Some(true)
        && p.logical("linear system scaling").unwrap_or(true);
    let scale_opts = || -> Result<elmer_sparse::ScaleOptions, SolverError> {
        let method = if p.logical("linear system row equilibration") == Some(true) {
            elmer_sparse::ScalingMethod::RowEquilibration
        } else {
            elmer_sparse::ScalingMethod::from_keyword(p.string("linear system scaling method").unwrap_or("diagonal"))?
        };
        Ok(elmer_sparse::ScaleOptions {
            method,
            apply: true,
            rhs_scaling: bnorm != 0.0,
        })
    };
    if complex {
        let mut ac = crate::linsys::to_complex(matrix)?;
        let mut bc = crate::linsys::pair_up(b);
        let mut xc = crate::linsys::pair_up(x);
        let scaling = if scale {
            Some(elmer_sparse::scale_system(&mut ac, Some(&mut bc[..]), Some(&mut xc[..]), &scale_opts()?))
        } else {
            None
        };
        if nullify {
            xc.iter_mut().for_each(|v| *v = num_complex::Complex64::new(0.0, 0.0));
        }
        let r = crate::linsys::solve_complex(p, &ac, &bc, &mut xc, log, comm);
        if let Some(rs) = scaling {
            rs.unscale(&mut ac, Some(&mut bc[..]), Some(&mut xc[..]));
        }
        for (k, z) in xc.iter().enumerate() {
            x[2 * k] = z.re;
            x[2 * k + 1] = z.im;
        }
        r
    } else {
        // SolveCore.F90:1263-1268, 1540 and MatrixScaling.F90 method choice.
        let scaling = if scale {
            Some(elmer_sparse::scale_system(matrix, Some(&mut *b), Some(&mut *x), &scale_opts()?))
        } else {
            None
        };
        if nullify {
            x.iter_mut().for_each(|v| *v = 0.0);
        }
        let solved = match linear.select_for(p, matrix, log) {
            Err(e) => Err(e),
            Ok(LinearSolverChoice::Direct(mut d)) => d
                .factor(matrix)
                .and_then(|_| {
                    let mut y = b.to_vec();
                    d.solve_in_place(&mut y).map(|_| x.copy_from_slice(&y))
                })
                .map_err(SolverError::from),
            Ok(LinearSolverChoice::Iterative { mut solver, precond }) => {
                let m: Result<Box<dyn elmer_core::Preconditioner<f64> + Send>, SolverError> = match precond {
                    Some(f) => f(p, matrix),
                    None => Ok(Box::new(elmer_core::IdentityPrecond)),
                };
                m.and_then(|m| {
                    solver
                        .solve(matrix, &*m, b, x, comm)
                        .map(|_| ())
                        .map_err(SolverError::from)
                })
            }
        };
        // BackScaleLinearSystem before anything else reads the system.
        if let Some(rs) = scaling {
            rs.unscale(matrix, Some(b), Some(x));
        }
        solved
    }
}
