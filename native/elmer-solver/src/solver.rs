//! The `Solver` trait every physics crate implements, the registry that maps
//! a `Procedure = "Lib" "Entry"` to one, and [`SolverCtx`], the view a solver
//! gets of the running simulation (Elmer's `Solver_t` plus DefUtils.F90).

use std::collections::BTreeMap;
use std::sync::Arc;

use elmer_core::{
    normalize_key, Communicator, DependentKind, ElemId, FieldRef, LocalSystem, Logger, MeshView, ParamMap, Params,
    Value, NO_DOF,
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
        if let Some(d) = &self.sim.slots[self.slot].dofs {
            let mut slots = Vec::new();
            d.element_slots(e, &mut slots);
            for s in slots {
                match s {
                    NO_DOF => out.extend(std::iter::repeat_n(NO_DOF, dofs as usize)),
                    s => out.extend((0..dofs).map(|k| s * dofs + k)),
                }
            }
            return;
        }
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
        if let Some(d) = &self.sim.slots[self.slot].dofs {
            d.element_slots(e, out);
            return;
        }
        let v = self.sim.vars.by_id(id);
        out.extend(
            self.mesh
                .element(e)
                .nodes
                .iter()
                .map(|&n| v.slot(n).map_or(NO_DOF, |s| s as u32)),
        );
    }

    /// The solver's non-nodal dofs (`Element = "p:..."`), when a
    /// [`crate::DofsBackend`] numbered them: kernels take the element basis,
    /// rule, node order and bubble count from it.
    pub fn solver_dofs(&self) -> Option<&dyn crate::SolverDofs> {
        self.sim.slots[self.slot].dofs.as_deref()
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

    /// `DefaultStart` (DefUtils.F90:3983), for a solver to call when its
    /// visit begins: `Nonlinear System Reset Newton`, `Nonlinear System
    /// Nullify Guess`, then the `Pre Solvers`. `CutFEM` and `Harmonic Mode`
    /// are refused.
    pub fn default_start(&mut self) -> Result<(), SolverError> {
        for (key, what) in [("cutfem", "CutFEM"), ("harmonic mode", "Harmonic Mode")] {
            if self.params().logical(key) == Some(true) {
                return Err(SolverError::Unsupported(what.into()));
            }
        }
        let eq = self.params().string("equation").unwrap_or_default().to_string();
        self.info("DefaultStart", &format!("Starting solver: {eq}"), 10);
        let slot = &mut self.sim.slots[self.slot];
        if slot.state.newton_active && slot.params.logical("nonlinear system reset newton") == Some(true) {
            slot.state.newton_active = false;
        }
        if slot.params.logical("nonlinear system nullify guess") == Some(true) {
            if let Some(id) = slot.variable {
                self.sim.vars.by_id_mut(id).values.fill(0.0);
            }
        }
        self.sim.slots[self.slot].started = true;
        let (mesh, log, comm) = (self.mesh, self.log, self.comm);
        self.sim.run_slave_solvers(self.slot, "Pre Solvers", mesh, log, comm)?;
        self.sim.publish_live();
        Ok(())
    }

    /// `DefaultFinish` (DefUtils.F90:4127), for a solver to call when its
    /// visit is done: the `Post Solvers` run, then `Apply Explicit Control`
    /// updates `cpar` from the new field. `Solve Adjoint Equation` is
    /// refused.
    pub fn default_finish(&mut self) -> Result<(), SolverError> {
        if self.params().logical("solve adjoint equation") == Some(true) {
            return Err(SolverError::Unsupported("Solve Adjoint Equation".into()));
        }
        self.sim.slots[self.slot].finished = true;
        let (mesh, log, comm) = (self.mesh, self.log, self.comm);
        self.sim.run_slave_solvers(self.slot, "Post Solvers", mesh, log, comm)?;
        self.sim.publish_live();
        if self.params().logical("apply explicit control") == Some(true) {
            crate::control::apply_explicit_control(self)?;
        }
        Ok(())
    }

    /// `DefaultInitialize`: run the `Slave Solvers` and `Nonlinear Pre
    /// Solvers`, create the matrix on first use, then zero the matrix values
    /// and the right-hand side.
    pub fn default_initialize(&mut self) -> Result<(), SolverError> {
        // DefUtils.F90:3897, before the slave solvers.
        if self.params().logical("apply explicit control") == Some(true) {
            crate::control::apply_explicit_control(self)?;
        }
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
            let matrix = match &slot.dofs {
                Some(d) => d.create_matrix(v.dofs),
                None => backend.create_matrix(self.mesh, &slot.active, &perm, v.dofs),
            };
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
        // UpdateGlobalEquations: a passive element adds nothing.
        if self.is_passive(e)? {
            return Ok(());
        }
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
    /// one of its components, by its declared name or `"<name> k"`), the
    /// value at each node,
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
        // DefaultDirichletBCs: `ComponentName(name, DOF)`, the declared
        // component names or `"<name> k"`.
        let names = self.sim.vars.component_names(id);
        let keys: Vec<(usize, String)> = if v.dofs == 1 {
            vec![(0, v.name.clone())]
        } else {
            (0..v.dofs)
                .map(|k| (k, names.get(k).cloned().unwrap_or_else(|| format!("{} {}", v.name, k + 1))))
                .collect()
        };
        let mut fixed: BTreeMap<u32, f64> = BTreeMap::new();
        // Non-nodal dofs (DefUtils.F90:6355): of this solver's own variable,
        // or of a variable numbered over the same dofs (a block variable),
        // through the dof index each layout slot stands for.
        let slot = &self.sim.slots[self.slot];
        let own = slot.variable == Some(id);
        let layout = slot.dofs.as_deref().filter(|d| own || v.perm.as_ref().is_some_and(|p| p.len() == d.perm().perm.len()));
        let raw: Vec<u32> = match layout {
            Some(d) if !own => {
                let perm = d.perm();
                let mut raw = vec![NO_DOF; perm.n_slots];
                for (i, &s) in perm.perm.iter().enumerate() {
                    if s != NO_DOF {
                        raw[s as usize] = i as u32;
                    }
                }
                raw
            }
            _ => Vec::new(),
        };
        let extra_slot = |s: u32| -> Option<usize> {
            if own {
                Some(s as usize)
            } else {
                raw.get(s as usize).copied().filter(|&i| i != NO_DOF).and_then(|i| v.slot(i))
            }
        };
        let mut extra = Vec::new();
        for e in self.dirichlet_element_order()? {
            let bc = self.bc_params(e).expect("boundary_elements filters on BC");
            if let Some(d) = layout {
                d.boundary_extra_slots(e, &mut extra);
            }
            for (k, key) in &keys {
                if bc.get(&format!("{key} condition")).is_some() {
                    return Err(SolverError::Unsupported(format!("Dirichlet mask > {key} Condition <")));
                }
                match bc.get(key) {
                    Some(Value::Dependent(dep))
                        if !extra.is_empty() && matches!(dep.kind, DependentKind::Procedure { .. }) =>
                    {
                        return Err(SolverError::Unsupported(format!(
                            "p-element Dirichlet value > {key} < given by a procedure (LocalBcBDOFs)"
                        )));
                    }
                    // Any other value fixes the extra (edge) dofs to zero.
                    Some(_) => {
                        for &s in &extra {
                            if let Some(t) = extra_slot(s) {
                                fixed.insert((offset + t * v.dofs + k) as u32, 0.0);
                            }
                        }
                    }
                    None => {}
                }
                let Some(vals) = self.real_at_nodes(bc, key, e)? else { continue };
                for (&n, &x) in self.mesh.element(e).nodes.iter().zip(&vals) {
                    if let Some(s) = v.slot(n) {
                        fixed.insert((offset + s * v.dofs + k) as u32, x);
                    }
                }
            }
        }
        // SetDirichletBoundaries' point values (BoundaryConditionUtils.F90:803):
        // a condition with `Target Nodes` (or `Target Coordinates`) sets the
        // variable at those nodes.
        for bc in &self.sim.model.boundary_conditions {
            if !keys.iter().any(|(_, key)| bc.get(key).is_some()) {
                continue;
            }
            let Some(nodes) = self.target_nodes(bc) else { continue };
            for (k, key) in &keys {
                let Some(value) = bc.get(key) else { continue };
                if bc.get(&format!("{key} condition")).is_some() {
                    return Err(SolverError::Unsupported(format!("Dirichlet mask > {key} Condition < on target nodes")));
                }
                let b = &self.sim.backends;
                let vals = crate::expr::eval_at_nodes_with(value, &nodes, &self.fields(), b.matc.as_deref(), &b.udfs)?;
                for (&n, &x) in nodes.iter().zip(&vals) {
                    if let Some(s) = v.slot(n) {
                        fixed.insert((offset + s * v.dofs + k) as u32, x);
                    }
                }
            }
        }
        // SetDirichletBoundaries: `<name> Single Node` after the others.
        for (k, key) in &keys {
            for (n, x) in self.single_nodes(id, key)? {
                if let Some(s) = v.slot(n) {
                    fixed.insert((offset + s * v.dofs + k) as u32, x);
                }
            }
        }
        Ok(fixed.into_iter().collect())
    }

    /// `CheckPassiveElement` (ElementGeometry.F90:341): a bulk element whose
    /// body force gives `<var> Passive` is passive when the value is
    /// positive at all of its nodes; at some of them, when `Passive Element
    /// Min Nodes` or `Active Element Min Nodes` (solver, else body force)
    /// says so, else when at more than half.
    pub fn is_passive(&self, e: ElemId) -> Result<bool, SolverError> {
        let Some(id) = self.sim.slots[self.slot].variable else { return Ok(false) };
        if (e.0 as usize) >= self.mesh.n_bulk() {
            return Ok(false);
        }
        let key = format!("{} passive", self.sim.vars.by_id(id).name);
        if !self.sim.model.body_forces.iter().any(|b| b.get(&key).is_some()) {
            return Ok(false);
        }
        let Some(bf) = self.body_params(e).and_then(|b| b.body_force) else { return Ok(false) };
        let Some(vals) = self.real_at_nodes(bf, &key, e)? else { return Ok(false) };
        let n = vals.len();
        let pass = vals.iter().filter(|&&v| v > 0.0).count();
        let min = |k: &str| self.params().integer(k).or_else(|| bf.integer(k)).filter(|&m| m > 0);
        Ok(if pass == 0 {
            false
        } else if pass == n {
            true
        } else if let Some(m) = min("passive element min nodes") {
            pass >= m as usize
        } else if let Some(m) = min("active element min nodes") {
            pass > n - (m as usize).min(n)
        } else {
            2 * pass > n
        })
    }

    /// `DefaultDirichletBCs`' passive part (BoundaryConditionUtils.F90:941):
    /// a node of a passive element that no other active element has keeps
    /// its current value, unless a condition already fixes it.
    fn passive_values(&self, fixed: &mut BTreeMap<u32, f64>) -> Result<(), SolverError> {
        let Some(id) = self.sim.slots[self.slot].variable else { return Ok(()) };
        let v = self.sim.vars.by_id(id);
        let key = format!("{} passive", v.name);
        if !self.sim.model.body_forces.iter().any(|b| b.get(&key).is_some()) {
            return Ok(());
        }
        if self.sim.model.boundary_conditions.iter().any(|b| b.logical("passive target") == Some(true)) {
            return Err(SolverError::Unsupported("Passive Target (GetPassiveBoundary)".into()));
        }
        let active: Vec<ElemId> = self.active_elements().collect();
        let mut passive = Vec::with_capacity(active.len());
        let mut used = vec![false; self.mesh.n_nodes()];
        for &e in &active {
            let p = self.is_passive(e)?;
            if !p {
                for &k in self.mesh.element(e).nodes {
                    used[k as usize] = true;
                }
            }
            passive.push(p);
        }
        for (&e, &p) in active.iter().zip(&passive) {
            if !p {
                continue;
            }
            for &k in self.mesh.element(e).nodes {
                if used[k as usize] {
                    continue;
                }
                let Some(s) = v.slot(k) else { continue };
                for l in 0..v.dofs {
                    let d = (s * v.dofs + l) as u32;
                    fixed.entry(d).or_insert(v.values[d as usize]);
                }
            }
        }
        Ok(())
    }

    /// The boundary elements in the order `SetDirichletBoundaries` visits
    /// them (BoundaryConditionUtils.F90:487-640): element order, or with
    /// `Set Dirichlet BCs by BC Numbering` condition by condition (in the
    /// solver's `Dirichlet BC Order`, if given), so that at a node two
    /// conditions share the later one wins.
    fn dirichlet_element_order(&self) -> Result<Vec<ElemId>, SolverError> {
        const CALLER: &str = "SetDirichletBoundaries";
        let by_bc = self.sim.model.simulation.logical("set dirichlet bcs by bc numbering") == Some(true);
        let nbc = self.sim.model.boundary_conditions.len();
        let order = match self.params().get("dirichlet bc order") {
            Some(Value::Integer(v)) if by_bc => {
                if v.len() != nbc {
                    return Err(SolverError::fatal(CALLER, "'Dirichlet BC Order' is the wrong length!"));
                }
                Some(v.clone())
            }
            Some(_) if !by_bc => {
                self.warn(
                    CALLER,
                    "Requested 'Dirichlet BC Order' but not 'Set Dirichlet BCs by BC Numbering', ignoring...",
                );
                None
            }
            _ => None,
        };
        let elements: Vec<ElemId> = self.boundary_elements().collect();
        if !by_bc {
            return Ok(elements);
        }
        let mut out = Vec::with_capacity(elements.len());
        for i in 1..=nbc {
            let bc = match &order {
                Some(o) => o[i - 1] as usize,
                None => i,
            };
            out.extend(elements.iter().copied().filter(|&e| self.bc_id(e) == Some(bc)));
        }
        Ok(out)
    }

    /// The nodes of a condition's `Target Nodes` (1-based in the input), or
    /// of its `Target Coordinates`: for each row the closest node within
    /// `Target Coordinates Eps` over as many coordinates as given
    /// (`TargetCoordinatesToTargetNodes`, `FindClosestNode`). `None` without
    /// either.
    pub(crate) fn target_nodes(&self, bc: &ParamMap) -> Option<Vec<u32>> {
        let n = self.mesh.n_nodes();
        if let Some(Value::Integer(ids)) = bc.get("target nodes") {
            return Some(ids.iter().filter(|&&i| i >= 1 && i as usize <= n).map(|&i| i as u32 - 1).collect());
        }
        let Some(Value::Real { rows, cols, data }) = bc.get("target coordinates") else { return None };
        let eps = bc.real("target coordinates eps").map_or(f64::MAX, |e| e * e);
        let mut out = Vec::new();
        for r in 0..*rows {
            let c = &data[r * cols..(r + 1) * cols];
            let mut best: Option<(f64, u32)> = None;
            for node in 0..n as u32 {
                let x = self.mesh.coord(node);
                let d: f64 = c.iter().take(3).enumerate().map(|(k, ck)| (x[k] - ck).powi(2)).sum();
                if d > eps {
                    continue;
                }
                if best.is_none_or(|(m, _)| d < m) {
                    best = Some((d, node));
                }
            }
            out.extend(best.map(|(_, node)| node));
        }
        Some(out)
    }

    /// `DefaultDirichletBCs` + `EnforceDirichletConditions`.
    pub fn default_dirichlet(&mut self) -> Result<(), SolverError> {
        self.default_finish_boundary_assembly()?;
        let mut fixed: BTreeMap<u32, f64> = self.dirichlet_values()?.into_iter().collect();
        self.passive_values(&mut fixed)?;
        crate::constraint_modes::set_boundaries(self, &mut fixed)?;
        let fixed: Vec<(u32, f64)> = fixed.into_iter().collect();
        self.lumped_rows()?;
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

    /// `DefaultDirichletBCs` with the dof values a solver gathered itself:
    /// applies `fixed` to the system and records them as the constrained
    /// dofs (`A % ConstrainedDOF`), as [`default_dirichlet`](Self::default_dirichlet)
    /// does, so that linear constraints, source control and loads see them.
    pub fn apply_dirichlet_values(&mut self, fixed: &[(u32, f64)]) -> Result<(), SolverError> {
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
        backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, fixed, &slot.params)?;
        slot.asm.constrained = fixed.to_vec();
        Ok(())
    }

    /// `DefaultSolve` (DefUtils.F90:4470): returns early when `Linear System
    /// Solver Disabled`; with `Monolithic Slave` the system is merged with
    /// the slave solver's first (`MergeSlaveSolvers`, ContactUtils.F90:2619),
    /// solved and split again; otherwise see [`solve_system`](Self::solve_system).
    /// Returns the variable norm.
    pub fn default_solve(&mut self) -> Result<f64, SolverError> {
        // DefUtils.F90:4502-4508: `Linear System Namespace Number` n > 0
        // pushes the namespace `linsys<n>:`, so that `linsys<n>: <key>`
        // stands for `<key>` until the solve returns.
        let p = self.params();
        let ns = self.creal(p, "linear system namespace number")?.map_or(0, |x| x.round() as i64);
        if p.logical("linear system trialing") == Some(true) {
            return self.default_solve_trialing(ns.max(1));
        }
        if ns <= 0 {
            return self.default_solve_in_namespace();
        }
        self.info("DefaultSolve", &format!("Linear system namespace number: {}", crate::fortran_fmt::i2s(ns)), 7);
        crate::block_system::with_namespace(self, &format!("linsys{ns}:"), |ctx| ctx.default_solve_in_namespace())
    }

    /// `Linear System Trialing` (DefUtils.F90:4579-4604): the solve starts
    /// in namespace `linsys<n>:` and, while the linear solver does not
    /// converge, `SolveSystem` is repeated in the next one; none left is
    /// fatal. With `Linear System Trialing Conserve` the converged namespace
    /// becomes the next solve's `Linear System Namespace Number`, back to the
    /// first after `Linear System Trialing Conserve Rounds` reuses. Only
    /// `SolveSystem` is repeated, so the source control, monolithic slave,
    /// constraint, complex and nonlinear post-solver paths are refused.
    fn default_solve_trialing(&mut self, mut ns: i64) -> Result<f64, SolverError> {
        let p = self.params();
        for k in ["apply source control", "monolithic slave", "linear system complex", "harmonic mode"] {
            if p.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("Linear System Trialing with {k}")));
            }
        }
        if p.get("nonlinear post solvers").is_some() {
            return Err(SolverError::Unsupported("Linear System Trialing with Nonlinear Post Solvers".into()));
        }
        crate::constraints::generate(self)?;
        let slot = &self.sim.slots[self.slot];
        if slot.constraints.is_some() || slot.own_constraints.is_some() {
            return Err(SolverError::Unsupported("Linear System Trialing with linear constraints".into()));
        }
        self.sim.slots[self.slot].lin_converged = true;
        let mut norm;
        let mut first = true;
        loop {
            let prefix = format!("linsys{ns}");
            if !first && !self.params().0.keys().any(|k| k.starts_with(&prefix)) {
                return Err(SolverError::fatal("DefaultSolve", "Exhausted all linear system strategies!"));
            }
            self.info("DefaultSolve", &format!("Linear system namespace number: {}", crate::fortran_fmt::i2s(ns)), 7);
            norm = crate::block_system::with_namespace(self, &format!("{prefix}:"), |ctx| {
                if ctx.params().logical("linear system solver disabled") == Some(true) {
                    ctx.info("DefaultSolve", "Solver disabled, exiting early!", 10);
                    return Ok(ctx.convergence().map_or(0.0, |c| c.norm));
                }
                ctx.solve_system()
            })?;
            if self.sim.slots[self.slot].lin_converged {
                break;
            }
            ns += 1;
            first = false;
        }
        let p = &mut self.sim.slots[self.slot].params;
        if p.logical("linear system trialing conserve") == Some(true) {
            if let Some(max) = p.integer("linear system trialing conserve rounds") {
                let prev = p.real("linear system namespace number").map(|x| x.round() as i64);
                let mut count = if prev == Some(ns) {
                    1 + p.integer("linear system namespace conserve count").unwrap_or(0)
                } else {
                    0
                };
                if count > max {
                    ns = 0;
                    count = 0;
                }
                p.insert("Linear System Namespace Conserve Count", Value::Integer(vec![count]));
            }
            p.insert("Linear System Namespace Number", Value::real(ns as f64));
        }
        Ok(norm)
    }

    fn default_solve_in_namespace(&mut self) -> Result<f64, SolverError> {
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
        if let Some(modes) = self.sim.slots[self.slot].cmodes.clone() {
            return crate::constraint_modes::solve_modes(self, modes);
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
        let schur = schur_matrix(sim, &sim.slots[slot_ix].params);
        let mg = crate::multigrid::transfers(sim, slot_ix, mesh)?;
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
        // SolveWithLinearRestriction: `Export Lagrange Multiplier` keeps the
        // multipliers in `Lagrange Multiplier <equation>` (or `Lagrange
        // Multiplier Name`), the initial guess of the next solve.
        let mult: Option<(crate::VarId, Vec<f64>)> = match slot.constraints.as_ref().or(slot.own_constraints.as_ref()) {
            Some(c) if slot.params.logical("export lagrange multiplier") == Some(true) => {
                let name = match slot.params.string("lagrange multiplier name") {
                    Some(n) => normalize_key(n),
                    None => normalize_key(&format!(
                        "lagrange multiplier {}",
                        slot.params.string("equation").unwrap_or("")
                    )),
                };
                let m = c.matrix.nrows;
                let mid = match sim.vars.id(&name) {
                    Some(mid) => mid,
                    None => {
                        log.info("SolveWithLinearRestriction", "Creating variable for Lagrange multiplier", 8);
                        sim.vars.add(crate::Variable::new(
                            &name,
                            elmer_core::FieldKind::Global,
                            1,
                            None,
                            vec![0.0; m],
                        ))?
                    }
                };
                let mv = sim.vars.by_id_mut(mid);
                if mv.values.len() < m {
                    log.info(
                        "SolveWithLinearRestriction",
                        &format!("Increasing Lagrange multiplier size to: {}", crate::fortran_fmt::i2s(m as i64)),
                        8,
                    );
                    mv.values.resize(m, 0.0);
                }
                Some((mid, mv.values[..m].to_vec()))
            }
            _ => None,
        };
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
                if let Some((_, init)) = &mult {
                    x[n..n + init.len()].copy_from_slice(init);
                }
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
            let env = BlockEnv {
                dofs: var.dofs,
                simulation: &sim.model.simulation,
                schur: schur.as_ref(),
                multigrid: mg.as_deref(),
            };
            sim.slots[slot_ix].lin_converged = solve_linear(&sim.backends.linear, &env, p, a, b, x, bnorm, nullify, log, comm)?;
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
            // CalculateLoads with the original matrix (BulkMatrix =
            // StiffMatrix) on the first n entries.
            sim.calculate_loads(slot_ix, n, true, mesh, log)?;
            sim.compute_change(slot_ix, false, Some(nc), true, mesh, log, comm)?;
            sim.calculate_loads(slot_ix, n, false, mesh, log)?;
            log.info("SolveWithLinearRestriction", "Picking solution from collection solution", 10);
            let var = sim.vars.by_id_mut(id);
            let tail = var.values[n..].to_vec();
            var.values.truncate(n);
            if own_x0 {
                var.nonlin_values = None;
            }
            if let Some((mid, _)) = mult {
                log.info("SolveWithLinearRestriction", "Separating Lagrange multiplier from collection solution", 10);
                let relax = crate::expr::eval_creal(
                    &sim.slots[slot_ix].params,
                    "lagrange multiplier relaxation factor",
                    &crate::expr::SimFields { vars: &sim.vars, mesh },
                    sim.backends.matc.as_deref(),
                    &sim.backends.udfs,
                )?;
                let mv = sim.vars.by_id_mut(mid);
                for (old, new) in mv.values.iter_mut().zip(&tail) {
                    *old = match relax {
                        Some(r) => (1.0 - r) * *old + r * new,
                        None => *new,
                    };
                }
                mv.values_changed = true;
            }
            log.info("SolveWithLinearRestriction", "All done", 10);
        } else {
            let x = &mut var.values[..n];
            let env = BlockEnv {
                dofs: var.dofs,
                simulation: &sim.model.simulation,
                schur: schur.as_ref(),
                multigrid: mg.as_deref(),
            };
            sim.slots[slot_ix].lin_converged =
                solve_linear(&sim.backends.linear, &env, p, &mut sys.matrix, b, x, bnorm, nullify, log, comm)?;
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

/// What `Linear System Solver = Block` needs beyond the section
/// (`BlockSolveExt`): the variable's components, the Simulation section and
/// the matrix of the `Block Matrix Schur Variable` solver.
pub(crate) struct BlockEnv<'a> {
    pub dofs: usize,
    pub simulation: &'a ParamMap,
    pub schur: Option<&'a elmer_core::CsrMatrix<f64>>,
    /// The interpolation operators of geometric multigrid, level 2 up.
    pub multigrid: Option<&'a [elmer_core::CsrMatrix<f64>]>,
}

/// The matrix of the solver whose variable `Block Matrix Schur Variable`
/// names (`BlockPrecMatrix`, BlockSolve.F90:2523), if any.
pub(crate) fn schur_matrix(sim: &Simulation, p: &ParamMap) -> Option<elmer_core::CsrMatrix<f64>> {
    let name = p.string("block matrix schur variable")?;
    let id = sim.vars.id(name)?;
    let slot = sim.slots.iter().find(|s| s.variable == Some(id))?;
    slot.system.as_ref().map(|s| s.matrix.clone())
}

/// `BlockError` as the framework's error, one to one.
pub(crate) fn block_error(e: elmer_block::BlockError) -> SolverError {
    match e {
        elmer_block::BlockError::Fatal { caller, msg } => SolverError::Fatal { caller, msg },
        elmer_block::BlockError::Numerical { caller, msg } => SolverError::Numerical { caller, msg },
        elmer_block::BlockError::Unsupported(m) => SolverError::Unsupported(m),
        elmer_block::BlockError::Linalg(e) => SolverError::Linalg(e),
    }
}

/// `SolveLinearSystem`'s scaling, solve and back-scaling of `A x = b`
/// (SolveCore.F90:1254-1700): `Linear System Scaling` (with the complex
/// form for `Linear System Complex`), the nullify-guess flag, the solver
/// the section selects (`block` through `elmer-block`, after the scaling
/// as in Elmer). `bnorm` is `|b|`, which turns the rhs scaling on.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_linear(
    linear: &crate::LinearSolverRegistry,
    block: &BlockEnv<'_>,
    p: &ParamMap,
    matrix: &mut elmer_core::CsrMatrix<f64>,
    b: &mut [f64],
    x: &mut [f64],
    bnorm: f64,
    nullify: bool,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<bool, SolverError> {
    let complex = p.logical("linear system complex") == Some(true)
        && p.logical("linear system skip complex") != Some(true);
    if complex {
        log.info("SolveLinearSystem", "Assuming complex valued linear system", 6);
        if p.string("linear system solver").is_some_and(|s| elmer_core::normalize_key(s) == "block") {
            return Err(SolverError::Unsupported("Linear System Solver = Block with Linear System Complex".into()));
        }
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
        let r = crate::linsys::solve_complex(p, &ac, &bc, &mut xc, log, comm).map(|_| true);
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
        let is_block = p
            .string("linear system solver")
            .is_some_and(|s| elmer_core::normalize_key(s) == "block");
        let is_multigrid = p
            .string("linear system solver")
            .is_some_and(|s| elmer_core::normalize_key(s) == "multigrid");
        let solved = if is_multigrid {
            let interps = block
                .multigrid
                .ok_or_else(|| SolverError::Unsupported("Linear System Solver = Multigrid without its mesh hierarchy".into()))?;
            crate::multigrid::solve(p, matrix, b, x, interps, log).map(|_| true)
        } else if is_block {
            let opts = elmer_block::BlockOptions {
                additive_namespaces: block.simulation.logical("additive namespaces") == Some(true),
                outer_scaling: scaling.as_ref().map(|rs| elmer_block::OuterScaling {
                    diag: rs.diag(),
                    rhs_scaling: rs.rhs_scaling(),
                }),
                global_abort: block.simulation.logical("global abort not converged").unwrap_or(true),
                schur: block.schur,
            };
            elmer_block::block_solve(p, block.dofs, matrix, b, x, &opts, log, comm)
                .map(|_| true)
                .map_err(block_error)
        } else {
            // `Linear System Preconditioning = multigrid` with `MG Method =
            // algebraic`: the Krylov method with the AMG cycle as its
            // preconditioner (`MultigridPrec`), built here as the registry's
            // factories do not know the dofs.
            let amg_prec = crate::amg::wants_preconditioner(p)?;
            if amg_prec && comm.size() > 1 {
                return Err(SolverError::fatal("AMGSolve", "This algebraic multigrid is not parallel"));
            }
            // `Linear System Preconditioning = BILU<n>`: ILU(n) of the
            // block diagonal, which needs the dofs too.
            let bilu = crate::linsys::bilu_level(p);
            let bilu_prec = bilu.is_some() && block.dofs > 1;
            if bilu.is_some() && block.dofs == 1 {
                log.warn("IterSolver", "BILU for one dofs is equal to ILU!");
            }
            let p_sel;
            let p = if let (Some(k), false) = (bilu, bilu_prec) {
                let mut q = p.clone();
                q.insert("Linear System Preconditioning", elmer_core::Value::Str(format!("ilu{k}")));
                p_sel = q;
                &p_sel
            } else if amg_prec || bilu_prec {
                let mut q = p.clone();
                q.insert("Linear System Preconditioning", elmer_core::Value::Str("none".into()));
                p_sel = q;
                &p_sel
            } else {
                p
            };
            match linear.select_for(p, matrix, log) {
            Err(e) => Err(e),
            Ok(LinearSolverChoice::Direct(mut d)) => d
                .factor(matrix)
                .and_then(|_| {
                    let mut y = b.to_vec();
                    d.solve_in_place(&mut y).map(|_| x.copy_from_slice(&y))
                })
                .map(|_| true)
                .map_err(SolverError::from),
            Ok(LinearSolverChoice::Iterative { mut solver, precond }) => {
                let m: Result<Box<dyn elmer_core::Preconditioner<f64> + Send>, SolverError> = match precond {
                    _ if amg_prec => crate::amg::preconditioner(p, matrix, block.dofs),
                    _ if bilu_prec => crate::linsys::bilu(matrix, block.dofs, bilu.unwrap_or(0)),
                    Some(f) => f(p, matrix),
                    None => Ok(Box::new(elmer_core::IdentityPrecond)),
                };
                m.and_then(|m| {
                    // `Variable % LinConverged` (IterSolve.F90:1646, 1664).
                    solver
                        .solve(matrix, &*m, b, x, comm)
                        .map(|s| matches!(s.status, elmer_core::SolveStatus::Converged))
                        .map_err(SolverError::from)
                })
            }
            }
        };
        // BackScaleLinearSystem before anything else reads the system.
        if let Some(rs) = scaling {
            rs.unscale(matrix, Some(b), Some(x));
        }
        solved
    }
}
