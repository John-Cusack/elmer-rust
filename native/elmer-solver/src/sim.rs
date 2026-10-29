//! [`Simulation`]: the model, the field registry, the solver slots and the
//! scheduling that drives them (ElmerSolver.F90 `ExecSimulation`,
//! MainUtils.F90 `SolveEquations`/`SolverActivate`/`SingleSolver`).

use std::sync::Arc;

use elmer_core::{
    normalize_key, Communicator, CsrMatrix, ElemId, FieldKind, Logger, MeshView, ParamMap, Params, Value, NO_DOF,
};

use crate::backend::{Backends, Perm};
use crate::change::{compute_change, ChangeEnv, ChangeReport, SolverState};
use crate::error::SolverError;
use crate::fortran_fmt::{es, i2s, list_int};
use crate::model::Model;
use crate::reference::ReferenceTally;
use crate::solver::{procedure_of, Solver, SolverCtx, SolverRegistry};
use crate::variable::{VarId, Variable, VariableStore};

/// `Exec Solver` (MainUtils.F90:966-1010 AddExecWhenFlag).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ExecWhen {
    Never,
    Always,
    /// `before simulation` / `before all`.
    BeforeAll,
    /// `after simulation` / `after all`.
    AfterAll,
    BeforeTimestep,
    AfterTimestep,
    BeforeSaving,
    AfterSaving,
}

impl ExecWhen {
    pub fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        if let Some(s) = p.string("exec solver") {
            return Ok(match normalize_key(s).as_str() {
                "never" => ExecWhen::Never,
                "always" => ExecWhen::Always,
                "after simulation" | "after all" => ExecWhen::AfterAll,
                "before simulation" | "before all" => ExecWhen::BeforeAll,
                "before timestep" => ExecWhen::BeforeTimestep,
                "after timestep" => ExecWhen::AfterTimestep,
                "before saving" => ExecWhen::BeforeSaving,
                "after saving" => ExecWhen::AfterSaving,
                "predictor-corrector" | "when created" | "after control" => {
                    return Err(SolverError::Unsupported(format!("Exec Solver = {s}")))
                }
                _ => return Err(SolverError::fatal("AddExecWhenFlag", format!("Unknown \"exec solver\" flag: {s}"))),
            });
        }
        let flag = |k: &str| p.logical(k) == Some(true);
        Ok(if flag("before all") || flag("before simulation") {
            ExecWhen::BeforeAll
        } else if flag("after all") || flag("after simulation") {
            ExecWhen::AfterAll
        } else if flag("before timestep") {
            ExecWhen::BeforeTimestep
        } else if flag("after timestep") {
            ExecWhen::AfterTimestep
        } else if flag("before saving") {
            ExecWhen::BeforeSaving
        } else if flag("after saving") {
            ExecWhen::AfterSaving
        } else {
            ExecWhen::Always
        })
    }
}

/// An assembled system: the matrix (pattern from the assembly backend) and
/// right-hand side.
#[derive(Clone, Debug)]
pub struct LinearSystem {
    pub matrix: CsrMatrix<f64>,
    pub rhs: Vec<f64>,
}

/// One `Solver n` section at run time (Elmer `Solver_t`).
pub struct SolverSlot {
    /// 1-based Elmer solver id.
    pub id: usize,
    pub params: ParamMap,
    pub exec_when: ExecWhen,
    pub variable: Option<VarId>,
    /// Shared with the variable.
    pub perm: Option<Arc<[u32]>>,
    /// Bulk elements on bodies whose equation activates this solver.
    pub active: Vec<ElemId>,
    /// `Solver % TimeOrder` (0 steady, 1 first order).
    pub time_order: u8,
    /// `Solver % Order` (BDF order, 1 otherwise).
    pub order: usize,
    /// `Solver % DoneTime`: timesteps initialised so far.
    pub done_time: usize,
    pub times_visited: usize,
    pub state: SolverState,
    pub system: Option<LinearSystem>,
    pub(crate) instance: Option<Box<dyn Solver>>,
}

impl SolverSlot {
    /// A slot with no variable or procedure; tests and the CLI use it to
    /// build slots by hand.
    pub fn new(id: usize, params: ParamMap) -> Result<Self, SolverError> {
        let exec_when = ExecWhen::from_params(&params)?;
        Ok(Self {
            id,
            params,
            exec_when,
            variable: None,
            perm: None,
            active: Vec::new(),
            time_order: 0,
            order: 1,
            done_time: 0,
            times_visited: 0,
            state: SolverState::default(),
            system: None,
            instance: None,
        })
    }

    pub fn has_procedure(&self) -> bool {
        self.instance.is_some()
    }

    pub fn set_instance(&mut self, s: Box<dyn Solver>) {
        self.instance = Some(s);
    }
}

/// A parsed `Variable` keyword.
#[derive(Clone, Debug, PartialEq)]
pub struct VarSpec {
    pub name: String,
    pub dofs: usize,
    /// Per-dof alias names (empty for a scalar).
    pub components: Vec<String>,
    pub kind: FieldKind,
    pub output: bool,
}

/// Parses `Variable = [-dofs n] [-nooutput] [-global|-elem|-dg|-ip] Name[A:k B:m]`
/// with `Variable DOFs` as the fallback dof count.
pub fn parse_variable_spec(s: &str, dofs_kw: Option<i32>) -> Result<VarSpec, SolverError> {
    let mut dofs: Option<usize> = None;
    let mut kind = FieldKind::Nodal;
    let mut output = true;
    let mut rest: Vec<&str> = Vec::new();
    let mut toks = s.split_whitespace().peekable();
    while let Some(t) = toks.next() {
        match t.to_ascii_lowercase().as_str() {
            "-dofs" => {
                let n = toks
                    .next()
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| SolverError::fatal("VariableSpec", format!("bad -dofs in: {s}")))?;
                dofs = Some(n);
            }
            "-nooutput" => output = false,
            "-global" => kind = FieldKind::Global,
            "-elem" => kind = FieldKind::Elemental,
            "-dg" => kind = FieldKind::Dg,
            "-ip" => kind = FieldKind::Ip,
            "-nodal" => kind = FieldKind::Nodal,
            o if o.starts_with('-') => return Err(SolverError::Unsupported(format!("variable option {t}"))),
            _ => rest.push(t),
        }
    }
    let full = rest.join(" ");
    let (name, comps) = match full.split_once('[') {
        Some((n, c)) => (n.trim(), Some(c.trim_end_matches(']'))),
        None => (full.as_str(), None),
    };
    let name = normalize_key(name);
    let mut components = Vec::new();
    if let Some(c) = comps {
        for item in c.split_whitespace() {
            let (cn, k) = item
                .split_once(':')
                .ok_or_else(|| SolverError::fatal("VariableSpec", format!("bad component {item}")))?;
            let k: usize = k
                .parse()
                .map_err(|_| SolverError::fatal("VariableSpec", format!("bad component {item}")))?;
            let cn = normalize_key(cn);
            if k == 1 {
                components.push(cn);
            } else {
                components.extend((1..=k).map(|j| format!("{cn} {j}")));
            }
        }
        let n = components.len();
        if dofs.is_some_and(|d| d != n) {
            return Err(SolverError::fatal("VariableSpec", format!("dofs mismatch in: {s}")));
        }
        dofs = Some(n);
    }
    let dofs = dofs.or(dofs_kw.map(|d| d.max(1) as usize)).unwrap_or(1);
    if components.is_empty() && dofs > 1 {
        components = (1..=dofs).map(|k| format!("{name} {k}")).collect();
    }
    Ok(VarSpec {
        name,
        dofs,
        components,
        kind,
        output,
    })
}

/// The simulation state.
pub struct Simulation {
    pub model: Model,
    pub vars: VariableStore,
    /// In solver id order.
    pub slots: Vec<SolverSlot>,
    pub backends: Backends,
    pub reference: ReferenceTally,
}

/// The run outcome of [`Simulation::run`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RunReport {
    /// Timesteps (steady iterations for steady state) taken.
    pub steps: usize,
    pub steady_state_reached: bool,
}

/// The globals `ExecSimulation` keeps (ElmerSolver.F90:1091-1107).
const GLOBALS: [&str; 6] = ["time", "timestep", "timestep size", "timestep interval", "nonlin iter", "coupled iter"];

impl Simulation {
    pub fn new(model: Model) -> Self {
        let mut vars = VariableStore::new();
        for g in GLOBALS {
            vars.set_scalar(g, 0.0);
        }
        Self {
            model,
            vars,
            slots: Vec::new(),
            backends: Backends::default(),
            reference: ReferenceTally::default(),
        }
    }

    pub fn is_transient(&self) -> bool {
        self.model
            .simulation
            .string("simulation type")
            .is_some_and(|s| normalize_key(s) == "transient")
    }

    pub fn is_scanning(&self) -> bool {
        self.model
            .simulation
            .string("simulation type")
            .is_some_and(|s| normalize_key(s) == "scanning")
    }

    /// Builds a slot for every `Solver n` section: procedure from the
    /// registry, active elements, primary variable (perm from the assembly
    /// backend, or node order of the active elements when none is
    /// installed), constant initial conditions, time order and history
    /// storage; then calls each solver's `init`.
    pub fn setup(
        &mut self,
        mesh: &dyn MeshView,
        registry: &SolverRegistry,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let transient = self.is_transient();
        for i in 0..self.model.solvers.len() {
            let params = self.model.solvers[i].clone();
            if params.0.is_empty() {
                continue;
            }
            let mut slot = SolverSlot::new(i + 1, params)?;
            if let Some((lib, entry)) = procedure_of(&slot.params) {
                if let Some(f) = registry.get(&lib, &entry) {
                    slot.instance = Some(f(&slot.params)?);
                }
            }
            slot.active = (0..mesh.n_bulk() as u32)
                .map(ElemId)
                .filter(|&e| {
                    self.model
                        .body_of_tag(mesh.element(e).tag)
                        .is_some_and(|b| self.model.solver_active_on_body(i + 1, b))
                })
                .collect();
            if transient {
                self.set_time_order(&mut slot)?;
            }
            if let Some(vs) = slot.params.string("variable").map(str::to_string) {
                let spec = parse_variable_spec(&vs, slot.params.integer("variable dofs"))?;
                self.add_solver_variable(&mut slot, &spec, mesh)?;
            }
            self.slots.push(slot);
        }
        for ix in 0..self.slots.len() {
            if let Some(mut inst) = self.slots[ix].instance.take() {
                let r = inst.init(&mut SolverCtx::new(self, ix, mesh, log, comm));
                self.slots[ix].instance = Some(inst);
                r?;
            }
        }
        Ok(())
    }

    /// `Solver % Order` and `TimeOrder` (MainUtils.F90:1429-1485, 2501-2520).
    fn set_time_order(&self, slot: &mut SolverSlot) -> Result<(), SolverError> {
        let p = &slot.params;
        let method = p
            .string("timestepping method")
            .or_else(|| self.model.simulation.string("timestepping method"))
            .map(normalize_key)
            .unwrap_or_else(|| "implicit euler".into());
        slot.order = match method.as_str() {
            "bdf" => p
                .integer("bdf order")
                .or_else(|| self.model.simulation.integer("bdf order"))
                .unwrap_or(2)
                .clamp(1, 5) as usize,
            "runge-kutta" | "adams-bashforth" | "adams-moulton" => {
                return Err(SolverError::Unsupported(format!("Timestepping Method = {method}")))
            }
            _ => 1,
        };
        if p.logical("eigen analysis") == Some(true) || p.logical("harmonic analysis") == Some(true) {
            slot.time_order = 0;
            return Ok(());
        }
        let k = p.integer("time derivative order").unwrap_or(1).clamp(1, 2);
        if k == 2 {
            return Err(SolverError::Unsupported("Time Derivative Order = 2".into()));
        }
        slot.time_order = 1;
        Ok(())
    }

    fn add_solver_variable(&mut self, slot: &mut SolverSlot, spec: &VarSpec, mesh: &dyn MeshView) -> Result<(), SolverError> {
        if spec.kind == FieldKind::Global {
            let mut v = Variable::new(&spec.name, FieldKind::Global, spec.dofs, None, vec![0.0; spec.dofs]);
            v.output = spec.output;
            slot.variable = Some(self.vars.add_with_components(v, &spec.components)?);
            return Ok(());
        }
        if spec.kind != FieldKind::Nodal {
            return Err(SolverError::Unsupported(format!("{:?} solver variable", spec.kind)));
        }
        let opt_bw = slot.params.logical("optimize bandwidth").unwrap_or(true);
        let perm = match &self.backends.assembly {
            Some(b) => b.build_perm(mesh, &slot.active, opt_bw),
            None => node_order_perm(mesh, &slot.active),
        };
        let perm_arc: Arc<[u32]> = perm.perm.into();
        let mut values = vec![0.0; perm.n_slots * spec.dofs];
        self.apply_initial_conditions(slot, spec, mesh, &perm_arc, &mut values)?;
        let mut v = Variable::new(&spec.name, FieldKind::Nodal, spec.dofs, Some(perm_arc.clone()), values);
        v.output = spec.output;
        if slot.time_order > 0 {
            let m = slot.order.max(slot.time_order as usize);
            v.prev_values = vec![vec![0.0; v.values.len()]; m];
        }
        slot.perm = Some(perm_arc);
        slot.variable = Some(self.vars.add_with_components(v, &spec.components)?);
        Ok(())
    }

    /// Constant `Initial Condition` values for the variable (or its
    /// components), node by node over the active elements.
    fn apply_initial_conditions(
        &self,
        slot: &SolverSlot,
        spec: &VarSpec,
        mesh: &dyn MeshView,
        perm: &[u32],
        values: &mut [f64],
    ) -> Result<(), SolverError> {
        let keys: Vec<(usize, &str)> = if spec.dofs == 1 {
            vec![(0, spec.name.as_str())]
        } else {
            spec.components.iter().enumerate().map(|(k, c)| (k, c.as_str())).collect()
        };
        for &e in &slot.active {
            let el = mesh.element(e);
            let Some(ic) = self
                .model
                .body_of_tag(el.tag)
                .and_then(|b| self.model.body_params(b))
                .and_then(|b| b.initial_condition)
            else {
                continue;
            };
            for &(k, key) in &keys {
                let Some(val) = ic.get(key) else { continue };
                let x = match val {
                    Value::Real { data, .. } if data.len() == 1 => data[0],
                    _ => return Err(SolverError::NotWired(format!("non-constant initial condition {key} (elmer-expr)"))),
                };
                for &n in el.nodes {
                    let s = perm[n as usize];
                    if s != NO_DOF {
                        values[s as usize * spec.dofs + k] = x;
                    }
                }
            }
        }
        Ok(())
    }

    fn slot_ix(&self, id: usize) -> Option<usize> {
        self.slots.iter().position(|s| s.id == id)
    }

    /// `InitializeTimestep` history commit for a first-order solver
    /// (SolverBasics.F90:2758-2770, 2786-2797): shift
    /// `PrevValues(:,Order..2)` down by one, then `PrevValues(:,1) = Values`,
    /// with `Order = min(DoneTime, Solver % Order)`.
    pub fn initialize_timestep(&mut self, ix: usize) -> Result<(), SolverError> {
        let slot = &mut self.slots[ix];
        slot.done_time += 1;
        let Some(id) = slot.variable else { return Ok(()) };
        if slot.time_order == 0 {
            return Ok(());
        }
        if slot
            .params
            .string("timestepping method")
            .is_some_and(|m| normalize_key(m) == "none")
        {
            return Ok(());
        }
        if self.model.simulation.logical("timestep extrapolation") == Some(true) {
            return Err(SolverError::Unsupported("Timestep Extrapolation".into()));
        }
        let order = slot.done_time.min(slot.order);
        let v = self.vars.by_id_mut(id);
        if v.prev_values.is_empty() {
            return Ok(());
        }
        for i in (1..order).rev() {
            let (lo, hi) = v.prev_values.split_at_mut(i);
            hi[0].copy_from_slice(&lo[i - 1]);
        }
        v.prev_values[0].copy_from_slice(&v.values);
        Ok(())
    }

    /// `ComputeChange` for slot `ix` with the iteration bookkeeping around
    /// it: the nonlinear level reads and advances `nonlin iter` (unless
    /// `Skip Advance Nonlinear iter`); the steady level reads `coupled iter`.
    /// `use_system` passes the slot's assembled matrix/rhs (residual
    /// measures, Aitken).
    #[allow(clippy::too_many_arguments)]
    pub fn compute_change(
        &mut self,
        ix: usize,
        steady: bool,
        nsize: Option<usize>,
        use_system: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<Option<ChangeReport>, SolverError> {
        let slot = &mut self.slots[ix];
        let Some(id) = slot.variable else {
            log.info(
                "ComputeChange",
                &format!("Solver variable not found for: {}", slot.params.string("equation").unwrap_or("")),
                10,
            );
            return Ok(None);
        };
        let p = &slot.params;
        for k in [
            "calculate derivative",
            "nonlinear calculate derivative",
            "nonlinear calculate velocity",
            "nonlinear update exported variables",
        ] {
            if p.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("ComputeChange option {k}")));
            }
        }
        if self.model.simulation.logical("convergence monitor") == Some(true) {
            return Err(SolverError::Unsupported("Convergence Monitor".into()));
        }
        let iter_no = if steady {
            self.vars.scalar("coupled iter").map_or(0, |v| v.round() as i32)
        } else {
            let it = self.vars.scalar("nonlin iter").map_or(0, |v| v.round() as i32);
            self.vars.by_id_mut(id).conv.nonlin_iter = it;
            if p.logical("skip advance nonlinear iter") != Some(true) {
                self.vars.set_scalar("nonlin iter", (it + 1) as f64);
            }
            it
        };
        let sys = if use_system { slot.system.as_ref() } else { None };
        let env = ChangeEnv {
            params: &slot.params,
            simulation: &self.model.simulation,
            steady,
            iter_no,
            nsize,
            matrix: sys.map(|s| &s.matrix),
            rhs: sys.map(|s| s.rhs.as_slice()),
            n_nodes: mesh.n_nodes(),
            log,
            comm,
        };
        let var = self.vars.by_id_mut(id);
        compute_change(&env, var, &mut slot.state).map(Some)
    }

    /// `SolverActivate` + `SingleSolver`: exec conditions, `nonlin iter`
    /// reset, then one call of the solver's procedure.
    pub fn solver_activate(
        &mut self,
        ix: usize,
        dt: f64,
        transient: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let p = &self.slots[ix].params;
        if p.get("start time").is_some() {
            return Err(SolverError::fatal("SolverActivate", "Use > Exec Condition = Real < instead of > Start Time <"));
        }
        if p.get("stop time").is_some() {
            return Err(SolverError::fatal("SolverActivate", "Use > Exec Condition = Real < instead of > Stop Time <"));
        }
        match p.get("exec condition") {
            Some(Value::Real { data, .. }) if data.len() == 1 => {
                if data[0] < 0.0 {
                    return Ok(());
                }
            }
            Some(_) => return Err(SolverError::NotWired("non-constant Exec Condition (elmer-expr)".into())),
            None => {}
        }
        let intervals = match p.get("exec intervals").or_else(|| p.get("exec interval")) {
            Some(Value::Integer(v)) => Some(v.clone()),
            _ => None,
        };
        let mut execi = 1;
        if let Some(iv) = intervals {
            let timei = self.vars.scalar("timestep interval").unwrap_or(1.0).round() as usize;
            let e = *iv
                .get(timei.max(1) - 1)
                .ok_or_else(|| SolverError::fatal("SolverActivate", "Exec Intervals shorter than Timestep Intervals"))?;
            if e == 0 {
                return Ok(());
            }
            let timestep = self.vars.scalar("timestep").unwrap_or(1.0).round() as i32;
            let off = match p.get("exec intervals offset") {
                Some(Value::Integer(v)) => v.get(timei.max(1) - 1).copied().unwrap_or(0),
                _ => 0,
            };
            execi = e;
            if (timestep - 1 - off).rem_euclid(e) != 0 {
                return Ok(());
            }
        }
        let auxiliary = p.logical("auxiliary solver") == Some(true);
        let mut dts = 1.0;
        let mut time_derivative_active = transient;
        if !auxiliary {
            dts = p.real("timestep scale").unwrap_or(1.0);
            if p.logical("timestep over intervals") == Some(true) {
                dts = execi as f64;
            }
            if transient {
                time_derivative_active = match p.logical("time derivative active") {
                    Some(b) => b,
                    None => p.real("time derivative condition").is_none_or(|c| c > 0.0),
                };
            }
        }
        if p.get("scanning loops").is_some() {
            return Err(SolverError::Unsupported("Scanning Loops".into()));
        }
        if !auxiliary {
            self.vars.set_scalar("nonlin iter", 1.0);
        }
        if let Some(mut inst) = self.slots[ix].instance.take() {
            let r = inst.solve(
                &mut SolverCtx::new(self, ix, mesh, log, comm),
                dts * dt,
                time_derivative_active,
            );
            self.slots[ix].instance = Some(inst);
            r?;
        }
        self.slots[ix].times_visited += 1;
        Ok(())
    }

    /// `SolveEquations` (MainUtils.F90:2987) without Runge-Kutta or
    /// predictor-corrector: before-timestep slots, the coupled loop over
    /// `Always` slots, after-timestep slots. Returns `SteadyStateReached`
    /// (steady runs only).
    #[allow(clippy::too_many_arguments)]
    pub fn solve_equations(
        &mut self,
        dt: f64,
        transient: bool,
        coupled_min: i32,
        coupled_max: i32,
        real_timestep: Option<usize>,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<bool, SolverError> {
        let scanning = self.is_scanning();
        let (mut cmin, mut cmax) = (coupled_min, coupled_max);
        if transient {
            cmin = self.model.simulation.integer("steady state min iterations").unwrap_or(0);
            cmax = self.model.simulation.integer("steady state max iterations").unwrap_or(1).max(1);
            for ix in 0..self.slots.len() {
                if self.slots[ix].has_procedure() {
                    self.initialize_timestep(ix)?;
                }
            }
        }
        if self
            .model
            .simulation
            .string("timestepping method")
            .is_some_and(|m| normalize_key(m) == "runge-kutta")
        {
            return Err(SolverError::Unsupported("Runge-Kutta timestepping".into()));
        }
        let skip_first = |s: &SolverSlot| real_timestep == Some(1) && s.params.logical("skip first timestep") == Some(true);

        log.info("SolveEquations", "Solvers before timestep", 12);
        for ix in 0..self.slots.len() {
            let s = &self.slots[ix];
            if s.has_procedure() && s.exec_when == ExecWhen::BeforeTimestep && !skip_first(s) {
                self.solver_activate(ix, dt, transient, mesh, log, comm)?;
                comm.barrier();
            }
        }

        let test_divergence = self.model.simulation.logical("convergence control within iterations") == Some(true);
        let after_converged: Vec<bool> = self
            .slots
            .iter()
            .map(|s| s.params.logical("coupled system after others converged") == Some(true))
            .collect();
        log.info("SolveEquations", "Solvers in main iteration loop", 12);
        let done = self.solve_coupled(
            dt,
            transient,
            scanning,
            cmin,
            cmax,
            test_divergence,
            &after_converged,
            real_timestep,
            mesh,
            log,
            comm,
        )?;
        let steady_reached = !transient && done.iter().all(|&d| d);

        log.info("SolveEquations", "Solvers after timestep", 12);
        for ix in 0..self.slots.len() {
            let s = &self.slots[ix];
            if s.has_procedure() && s.exec_when == ExecWhen::AfterTimestep && !skip_first(s) {
                self.solver_activate(ix, dt, transient, mesh, log, comm)?;
                comm.barrier();
            }
        }
        Ok(steady_reached)
    }

    /// `SolveCoupled` (MainUtils.F90:3321-3538). Returns `DoneThis`.
    #[allow(clippy::too_many_arguments)]
    fn solve_coupled(
        &mut self,
        dt: f64,
        transient: bool,
        scanning: bool,
        cmin: i32,
        cmax: i32,
        test_divergence: bool,
        after_converged: &[bool],
        real_timestep: Option<usize>,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<Vec<bool>, SolverError> {
        let ns = self.slots.len();
        let mut done = vec![false; ns];
        let mut test_convergence = false;
        let mut divergence_exit = false;
        let mut i = 1;
        log.info("SolveEquations", "Performing set of solvers in sequence", 12);
        while i <= cmax {
            if transient || scanning {
                if cmax > 1 {
                    log.info("SolveEquations", "-------------------------------------", 3);
                    log.info("SolveEquations", &format!(" Coupled system iteration: {}", list_int(i)), 3);
                    log.info("SolveEquations", "-------------------------------------", 3);
                }
                self.vars.set_scalar("coupled iter", i as f64);
            }
            done.iter_mut().for_each(|d| *d = false);
            for k in 0..ns {
                let s = &self.slots[k];
                if !s.has_procedure() {
                    log.warn("SolveEquations", "No routine related to solver!");
                    done[k] = true;
                    continue;
                }
                if real_timestep == Some(1) && s.params.logical("skip first timestep") == Some(true) {
                    done[k] = true;
                    continue;
                }
                if s.exec_when != ExecWhen::Always {
                    done[k] = true;
                    continue;
                }
                if after_converged[k] && !after_converged.iter().zip(&done).all(|(&a, &d)| a || d) {
                    continue;
                }
                let p = &s.params;
                let calc_derivative = p.logical("calculate derivative") == Some(true);
                if calc_derivative && !scanning {
                    return Err(SolverError::fatal(
                        "SolveEquations",
                        "> Calculate Derivative < should be active only for scanning!",
                    ));
                }
                let need_sol = calc_derivative
                    || p.string("steady state convergence measure").is_some_and(|m| m != "norm")
                    || p.get("steady state relaxation factor").is_some();
                if let Some(id) = s.variable {
                    let v = self.vars.by_id_mut(id);
                    v.conv.prev_norm = v.conv.norm;
                    if need_sol && !v.values.is_empty() {
                        v.steady_values = Some(v.values.clone());
                    }
                }
                let mut transient_solver = transient;
                if transient {
                    if let Some(c) = p.real("steady state condition") {
                        if c > 0.0 {
                            transient_solver = false;
                            log.info("SolveEquations", "Running solver in steady state", 6);
                        }
                    }
                }
                self.solver_activate(k, dt, transient_solver, mesh, log, comm)?;

                let s = &self.slots[k];
                let nonlin_converged = s.variable.map_or(-1, |id| self.vars.by_id(id).conv.nonlin_converged);
                if test_divergence {
                    divergence_exit = nonlin_converged > 1;
                    break;
                }
                if scanning || transient {
                    if cmax == 1 {
                        test_convergence = false;
                        if nonlin_converged < 0 {
                            test_convergence = s.params.get("reference norm").is_some();
                        }
                        if !test_convergence {
                            test_convergence = s.params.get("steady state relaxation factor").is_some();
                        }
                    } else {
                        test_convergence = i >= cmin;
                    }
                } else {
                    test_convergence = true;
                }
                if test_convergence || calc_derivative {
                    self.compute_change(k, true, None, false, mesh, log, comm)?;
                }
                if test_convergence {
                    done[k] = self.slots[k]
                        .variable
                        .is_some_and(|id| self.vars.by_id(id).conv.steady_converged != 0);
                }
                if done.iter().all(|&d| d) {
                    break;
                }
            }
            if divergence_exit || done.iter().all(|&d| d) {
                break;
            }
            i += 1;
        }
        if test_convergence && cmax > 1 && transient && !done.iter().all(|&d| d) {
            log.info(
                "SolveEquations",
                &format!("Coupled system iteration: {}", i2s(i.min(cmax) as i64)),
                4,
            );
            let abort = self.model.simulation.logical("coupled system abort not converged");
            crate::change::numerical_error(
                &self.model.simulation,
                "SolveEquations",
                "Coupled system did not converge",
                Some(abort.unwrap_or(false)),
                log,
            )?;
        }
        Ok(done)
    }

    fn exec_slots(
        &mut self,
        when: ExecWhen,
        dt: f64,
        transient: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        for ix in 0..self.slots.len() {
            if self.slots[ix].has_procedure() && self.slots[ix].exec_when == when {
                self.solver_activate(ix, dt, transient, mesh, log, comm)?;
            }
        }
        Ok(())
    }

    /// `ExecSimulation` (ElmerSolver.F90:2703) for steady state and
    /// fixed-size transient stepping, then `CompareToReferenceSolution`.
    /// Variable or adaptive steps, periodic time and scanning are refused.
    pub fn run(&mut self, mesh: &dyn MeshView, log: &dyn Logger, comm: &dyn Communicator) -> Result<RunReport, SolverError> {
        let simp = self.model.simulation.clone();
        let transient = self.is_transient();
        if self.is_scanning() {
            return Err(SolverError::Unsupported("Simulation Type = Scanning".into()));
        }
        for k in [
            "adaptive timestepping",
            "timestep size",
            "timestep function",
            "timestep ratios",
            "periodic time",
            "timestep start zero",
        ] {
            if simp.get(k).is_some_and(|v| *v != Value::Logical(false)) {
                return Err(SolverError::Unsupported(format!("Simulation keyword {k}")));
            }
        }
        let (intervals, sizes, coupled_max) = if transient {
            let iv = match simp.get("timestep intervals") {
                Some(Value::Integer(v)) => v.clone(),
                _ => {
                    return Err(SolverError::fatal(
                        "ElmerSolver",
                        "Keyword > Timestep Intervals < MUST be defined for transient and scanning simulations",
                    ))
                }
            };
            let sz = match simp.get("timestep sizes") {
                Some(Value::Real { data, .. }) => data.clone(),
                _ => {
                    return Err(SolverError::fatal(
                        "MAIN",
                        "Keyword [Timestep Sizes] MUST be defined for time dependent simulations",
                    ))
                }
            };
            if sz.len() < iv.len() {
                return Err(SolverError::fatal("MAIN", "Timestep Sizes shorter than Timestep Intervals"));
            }
            let cmax = simp.integer("steady state max iterations").unwrap_or(1).max(1);
            (iv, sz, cmax)
        } else {
            (vec![simp.integer("steady state max iterations").unwrap_or(1).max(1)], vec![1.0], 1)
        };
        let coupled_min = simp.integer("steady state min iterations").unwrap_or(1);
        let output_intervals = match simp.get("output intervals") {
            Some(Value::Integer(v)) => {
                if v.len() != intervals.len() {
                    return Err(SolverError::fatal(
                        "MAIN",
                        "> Output Intervals < should have the same size as > Timestep Intervals < !",
                    ));
                }
                v.clone()
            }
            _ => vec![1; intervals.len()],
        };
        let stepcount: i32 = intervals.iter().sum();

        self.exec_slots(ExecWhen::BeforeAll, 0.0, transient, mesh, log, comm)?;

        let mut report = RunReport::default();
        let mut cum = 0;
        let mut dt = 0.0;
        'outer: for (interval, &nsteps) in intervals.iter().enumerate() {
            // Without adaptive stepping RealTimestep equals the step in the interval.
            for timestep in 1..=nsteps {
                cum += 1;
                self.vars.set_scalar("timestep", cum as f64);
                dt = sizes[interval];
                let t = self.vars.scalar("time").unwrap_or(0.0) + dt;
                self.vars.set_scalar("time", t);
                self.vars.set_scalar("timestep size", dt);
                self.vars.set_scalar("timestep interval", (interval + 1) as f64);
                if !transient {
                    let c = self.vars.scalar("coupled iter").unwrap_or(0.0);
                    self.vars.set_scalar("coupled iter", c + 1.0);
                }
                if comm.rank() == 0 {
                    log.info("MAIN", " ", 3);
                    log.info("MAIN", "-------------------------------------", 3);
                    if transient {
                        log.info("MAIN", &format!("Time: {}/{}:{}", cum, stepcount, es(t, 12, 3)), 3);
                    } else {
                        // WRITE(Message,*) 'Steady state iteration: ',cum_Timestep
                        log.info("MAIN", &format!(" Steady state iteration: {}", list_int(cum)), 3);
                    }
                    log.info("MAIN", "-------------------------------------", 3);
                    log.info("MAIN", " ", 3);
                }
                let reached = self.solve_equations(
                    dt,
                    transient,
                    coupled_min,
                    coupled_max,
                    Some(timestep as usize),
                    mesh,
                    log,
                    comm,
                )?;
                report.steps = cum as usize;
                report.steady_state_reached = reached;
                let oi = output_intervals[interval];
                if oi != 0 && ((timestep - 1) % oi == 0 || reached) {
                    self.exec_slots(ExecWhen::BeforeSaving, dt, transient, mesh, log, comm)?;
                    self.exec_slots(ExecWhen::AfterSaving, dt, transient, mesh, log, comm)?;
                }
                if simp.real("exit condition").is_some_and(|c| c > 0.0) {
                    log.info("MAIN", "Found a positive exit condition, exiting...", 3);
                    break 'outer;
                }
                if reached && !transient && timestep >= coupled_min {
                    break 'outer;
                }
            }
        }

        self.exec_slots(ExecWhen::AfterAll, dt, transient, mesh, log, comm)?;
        crate::reference::compare_to_reference(self, mesh, log, comm)?;
        Ok(report)
    }

    /// Looks a slot up by 1-based solver id.
    pub fn slot(&self, id: usize) -> Option<&SolverSlot> {
        self.slot_ix(id).map(|i| &self.slots[i])
    }
}

/// Stand-in permutation until `elmer-assembly` is installed: nodes of the
/// active elements numbered in first-touch element order.
pub fn node_order_perm(mesh: &dyn MeshView, active: &[ElemId]) -> Perm {
    let mut perm = vec![NO_DOF; mesh.n_nodes()];
    let mut n = 0u32;
    for &e in active {
        for &node in mesh.element(e).nodes {
            if perm[node as usize] == NO_DOF {
                perm[node as usize] = n;
                n += 1;
            }
        }
    }
    Perm {
        perm,
        n_slots: n as usize,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variable_specs() {
        let s = parse_variable_spec("TempA", None).unwrap();
        assert_eq!((s.name.as_str(), s.dofs, s.components.len()), ("tempa", 1, 0));
        let s = parse_variable_spec("-dofs 2 Displacement", None).unwrap();
        assert_eq!(s.components, vec!["displacement 1", "displacement 2"]);
        let s = parse_variable_spec("Flow Solution[Velocity:2 Pressure:1]", None).unwrap();
        assert_eq!(s.name, "flow solution");
        assert_eq!(s.dofs, 3);
        assert_eq!(s.components, vec!["velocity 1", "velocity 2", "pressure"]);
        let s = parse_variable_spec("-nooutput -global Energy", None).unwrap();
        assert!(!s.output && s.kind == FieldKind::Global);
        assert_eq!(parse_variable_spec("U", Some(3)).unwrap().dofs, 3);
        assert!(parse_variable_spec("-weird U", None).is_err());
    }

    #[test]
    fn exec_when() {
        let mut p = ParamMap::new();
        assert_eq!(ExecWhen::from_params(&p).unwrap(), ExecWhen::Always);
        p.insert("Exec Solver", Value::Str("After Simulation".into()));
        assert_eq!(ExecWhen::from_params(&p).unwrap(), ExecWhen::AfterAll);
        p.insert("Exec Solver", Value::Str("bogus".into()));
        assert!(ExecWhen::from_params(&p).is_err());
        let mut p = ParamMap::new();
        p.insert("Before Timestep", Value::Logical(true));
        assert_eq!(ExecWhen::from_params(&p).unwrap(), ExecWhen::BeforeTimestep);
    }

    #[test]
    fn history_shift_follows_done_time() {
        let mut sim = Simulation::new(Model::default());
        let mut slot = SolverSlot::new(1, ParamMap::new()).unwrap();
        slot.time_order = 1;
        slot.order = 2;
        let mut v = Variable::new("u", FieldKind::Nodal, 1, None, vec![1.0]);
        v.prev_values = vec![vec![0.0], vec![0.0]];
        slot.variable = Some(sim.vars.add(v).unwrap());
        sim.slots.push(slot);
        let prev = |sim: &Simulation| sim.vars.get("u").unwrap().prev_values.clone();
        sim.initialize_timestep(0).unwrap();
        // DoneTime=1: order 1, no shift.
        assert_eq!(prev(&sim), vec![vec![1.0], vec![0.0]]);
        sim.vars.get_mut("u").unwrap().values[0] = 2.0;
        sim.initialize_timestep(0).unwrap();
        assert_eq!(prev(&sim), vec![vec![2.0], vec![1.0]]);
        sim.vars.get_mut("u").unwrap().values[0] = 3.0;
        sim.initialize_timestep(0).unwrap();
        assert_eq!(prev(&sim), vec![vec![3.0], vec![2.0]]);
    }
}
