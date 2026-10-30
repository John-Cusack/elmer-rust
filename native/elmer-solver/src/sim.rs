//! [`Simulation`]: the model, the field registry, the solver slots and the
//! scheduling that drives them (ElmerSolver.F90 `ExecSimulation`,
//! MainUtils.F90 `SolveEquations`/`SolverActivate`/`SingleSolver`).

use std::rc::Rc;
use std::sync::Arc;

use elmer_core::{
    normalize_key, Communicator, CsrMatrix, ElemId, FieldKind, Logger, MeshView, ParamMap, Params, Value, NO_DOF,
};
use elmer_expr::MatcEngine;

use crate::backend::{Backends, Perm};
use crate::change::{compute_change, ChangeEnv, ChangeReport, SolverState};
use crate::error::SolverError;
use crate::expr::{NoMatc, SharedMatc, SimFields};
use crate::fortran_fmt::{es, i2s, list_int};
use crate::live::LiveModel;
use crate::model::Model;
use crate::reference::ReferenceTally;
use crate::solver::{procedure_of, Solver, SolverCtx, SolverRegistry};
use crate::variable::{VarHandle, VarId, Variable, VariableStore};

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

/// `ddt` and `StepControl` of the adaptive scheme, kept across steps of
/// one run (`SAVE`d in `ExecSimulation`, reset on entry).
#[derive(Clone, Debug)]
struct AdaptiveState {
    ddt: f64,
    step_control: i32,
}

impl Default for AdaptiveState {
    fn default() -> Self {
        Self { ddt: -1.0, step_control: 0 }
    }
}

/// Fortran `Ew.d`: `0.ddd...dE+xx`, right-aligned in `w`.
fn e_fmt(x: f64, w: usize, d: usize) -> String {
    let body = if x == 0.0 {
        format!("0.{}E+00", "0".repeat(d))
    } else {
        let t = format!("{:.*e}", d - 1, x.abs());
        let (m, e) = t.split_once('e').expect("exponent form");
        let digits: String = m.chars().filter(|c| c.is_ascii_digit()).collect();
        let e: i32 = e.parse::<i32>().expect("exponent") + 1;
        let sign = if x < 0.0 { "-" } else { "" };
        format!("{sign}0.{digits}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    };
    format!("{body:>w$}")
}

/// Which solver groups a `SolveEquations` call runs (`BeforeTime`,
/// `AtTime`, `AfterTime`).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Phases {
    pub before: bool,
    pub at: bool,
    pub after: bool,
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
    /// `Solver % TimeOrder` (0 steady, 1 first order, 2 second order).
    pub time_order: u8,
    /// `Solver % Order` (BDF order, 1 otherwise).
    pub order: usize,
    /// `Solver % DoneTime`: timesteps initialised so far.
    pub done_time: usize,
    pub times_visited: usize,
    /// `Solver % dt`: the step size of the current visit (`Timestep Scale`
    /// applied).
    pub dt: f64,
    pub state: SolverState,
    pub system: Option<LinearSystem>,
    pub(crate) instance: Option<Box<dyn Solver>>,
    pub(crate) asm: crate::loads::AssemblyState,
    /// The constraint rows `Apply Integral BCs` generated for this solve.
    pub(crate) constraints: Option<crate::constraints::ConstraintMatrix>,
    /// Constraint rows the solver set itself.
    pub(crate) own_constraints: Option<crate::constraints::ConstraintMatrix>,
    /// Non-nodal dofs from [`crate::DofsBackend`]; `None` for a nodal solver.
    pub(crate) dofs: Option<Arc<dyn crate::SolverDofs>>,
    /// Boundaries (body forces past them) with a Robin projector so far.
    pub(crate) robin_projectors: std::collections::BTreeSet<usize>,
    /// Whether this visit called `DefaultStart` / `DefaultFinish`.
    pub(crate) started: bool,
    pub(crate) finished: bool,
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
            dt: 0.0,
            state: SolverState::default(),
            system: None,
            instance: None,
            asm: Default::default(),
            constraints: None,
            own_constraints: None,
            dofs: None,
            robin_projectors: Default::default(),
            started: false,
            finished: false,
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
        // ComponentNameStr (GeneralUtils.F90:1478): a component name runs
        // up to its ':' and may hold spaces (`Potential Re:1`); the count
        // follows the colon.
        let bad = || SolverError::fatal("ComponentName", format!("Syntax error in variable definition: {s}"));
        let mut rest = c;
        while !rest.trim().is_empty() {
            let colon = rest.find(':').ok_or_else(bad)?;
            let cn = normalize_key(rest[..colon].trim());
            let after = &rest[colon + 1..];
            let digits = after.len() - after.trim_start_matches(|ch: char| ch.is_ascii_digit()).len();
            let k: usize = after[..digits].parse().map_err(|_| bad())?;
            if cn.is_empty() {
                return Err(bad());
            }
            if k == 1 {
                components.push(cn);
            } else {
                components.extend((1..=k).map(|j| format!("{cn} {j}")));
            }
            rest = &after[digits..];
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
    /// Previous step sizes, newest first (`sPrevSizes`, ElmerSolver.F90:3064-3070).
    pub prev_sizes: Vec<f64>,
    /// The model view handed to user functions, once asked for.
    live: Option<Arc<LiveModel>>,
    /// `SharedWarned` of `CalculateLoads`: warn about shared boundary dofs once.
    pub(crate) flux_shared_warned: bool,
    /// `SetNormalizedKeywords` has run (`Mesh % EntityWeightsComputed`).
    pub(crate) keywords_normalized: bool,
}

/// The run outcome of [`Simulation::run`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RunReport {
    /// Timesteps (steady iterations for steady state) taken.
    pub steps: usize,
    pub steady_state_reached: bool,
}

/// Real keywords `ComputeChange` reads (SolverBasics.F90:3571-4100), most
/// with `ListGetCReal`: evaluated once per call when they are not constants.
const CHANGE_CREALS: &[&str] = &[
    "steady state relaxation factor",
    "nonlinear system relaxation factor",
    "steady state max norm",
    "nonlinear system max norm",
    "steady state convergence tolerance",
    "steady state divergence limit",
    "steady state exit condition",
    "nonlinear system convergence tolerance",
    "nonlinear system divergence limit",
    "nonlinear system exit condition",
    "nonlinear system newton after tolerance",
    "average solution value",
];

/// The globals `ExecSimulation` keeps (ElmerSolver.F90:1091-1107).
const GLOBALS: [&str; 6] = ["time", "timestep", "timestep size", "timestep interval", "nonlin iter", "coupled iter"];

impl Simulation {
    pub fn new(model: Model) -> Self {
        let mut vars = VariableStore::new();
        for g in GLOBALS {
            vars.set_scalar(g, 0.0);
        }
        // ElmerSolver.F90:1766: `scan` when any solver has `Scanning Loops`.
        if model.solvers.iter().any(|p| p.get("scanning loops").is_some()) {
            vars.set_scalar("scan", 0.0);
        }
        Self {
            model,
            vars,
            slots: Vec::new(),
            backends: Backends::default(),
            reference: ReferenceTally::default(),
            prev_sizes: vec![0.0; 5],
            live: None,
            flux_shared_warned: false,
            keywords_normalized: false,
        }
    }

    /// The live model view for user functions that read node coordinates,
    /// `Time` or other fields (Elmer passes them `Model`). Created on the
    /// first call and kept current from then on; see [`LiveModel`].
    pub fn live_model(&mut self, mesh: &dyn MeshView) -> Arc<LiveModel> {
        let live = self
            .live
            .get_or_insert_with(|| Arc::new(LiveModel::new(mesh, &self.vars)));
        live.publish(&self.vars);
        Arc::clone(live)
    }

    /// Copies the current field values into the live model view, if any.
    pub(crate) fn publish_live(&self) {
        if let Some(l) = &self.live {
            l.publish(&self.vars);
        }
    }

    /// Installs the MATC evaluator the front end parsed the `.sif` with, so
    /// keyword values see the variables and functions it defined.
    pub fn set_matc(&mut self, matc: Rc<dyn MatcEngine>) {
        self.backends.matc = Some(Box::new(SharedMatc(matc)));
    }

    /// The fields as `elmer-expr` reads them.
    pub fn fields<'a>(&'a self, mesh: &'a dyn MeshView) -> SimFields<'a> {
        SimFields { vars: &self.vars, mesh }
    }

    /// `ListGetCReal(section, key)` against this simulation's fields, MATC
    /// engine and function registry.
    pub fn creal(&self, section: &ParamMap, key: &str, mesh: &dyn MeshView) -> Result<Option<f64>, SolverError> {
        let b = &self.backends;
        crate::expr::eval_creal(section, key, &self.fields(mesh), b.matc.as_deref(), &b.udfs)
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
        // MeshLoad's MapInternalExternalBCs.
        let bc = |k: &str| self.model.boundary_conditions.iter().rposition(|b| b.logical(k) == Some(true));
        let (int_bc, ext_bc) = (bc("default internal bc"), bc("default external bc"));
        let (n_int, n_ext) = self.model.map_default_bcs(mesh)?;
        const MAP: &str = "MapInternalExternalBCs";
        if let Some(i) = int_bc {
            log.info(MAP, &format!("Default Internal BC: {}", crate::fortran_fmt::i2s(i as i64 + 1)), 8);
        }
        if let Some(i) = ext_bc {
            log.info(MAP, &format!("Default External BC: {}", crate::fortran_fmt::i2s(i as i64 + 1)), 8);
        }
        if n_int > 0 {
            log.info(MAP, &format!("\"Default Internal BC\" count: {}", crate::fortran_fmt::i2s(n_int as i64)), 6);
        }
        if n_ext > 0 {
            log.info(MAP, &format!("\"Default External BC\" count: {}", crate::fortran_fmt::i2s(n_ext as i64)), 6);
        }
        for i in 0..self.model.solvers.len() {
            if let Some(slot) = self.build_slot(i, mesh, registry, true)? {
                self.slots.push(slot);
            }
        }
        for ix in 0..self.slots.len() {
            self.init_slot(ix, mesh, log, comm)?;
        }
        // ElmerSolver.F90:578-585, after the initial conditions.
        for ix in 0..self.slots.len() {
            if self.slots[ix].params.logical("initialize exported variables") == Some(true) {
                self.update_exported_variables(ix, mesh, log)?;
            }
        }
        self.initialize_dirichlet(mesh)
    }

    /// The slot of `Solver i+1` (`AddEquationBasics`): procedure from the
    /// registry or the solver mode, active elements, time order and primary
    /// variable, with initial conditions when `apply_ics`. `None` for an
    /// empty section.
    fn build_slot(
        &mut self,
        i: usize,
        mesh: &dyn MeshView,
        registry: &SolverRegistry,
        apply_ics: bool,
    ) -> Result<Option<SolverSlot>, SolverError> {
        let params = self.model.solvers[i].clone();
        if params.0.is_empty() {
            return Ok(None);
        }
        let mut slot = SolverSlot::new(i + 1, params)?;
        if let Some((lib, entry)) = procedure_of(&slot.params) {
            if let Some(f) = registry.get(&lib, &entry) {
                slot.instance = Some(f(&slot.params)?);
            }
        }
        if let Some(s) = crate::coupled::solver_for_mode(&slot.params, registry)? {
            slot.instance = Some(s);
        }
        slot.active = (0..mesh.n_bulk() as u32)
            .map(ElemId)
            .filter(|&e| {
                self.model
                    .body_of_tag(mesh.element(e).tag)
                    .is_some_and(|b| self.model.solver_active_on_body(i + 1, b))
            })
            .collect();
        if self.is_transient() {
            self.set_time_order(&mut slot)?;
        }
        // AddEquationBasics (MainUtils.F90:1514-1520): a coupled or block
        // solver without a procedure makes its variables itself.
        let coupled = crate::coupled::is_coupled_system(&slot.params)
            || (slot.params.logical("block solver") == Some(true)
                && slot.params.get("procedure").is_none()
                && slot.params.logical("assembly solver") != Some(true));
        if let Some(vs) = slot.params.string("variable").filter(|_| !coupled).map(str::to_string) {
            let spec = parse_variable_spec(&vs, slot.params.integer("variable dofs"))?;
            self.publish_live();
            self.add_solver_variable(&mut slot, &spec, mesh, apply_ics)?;
            if let Some(id) = slot.variable {
                self.add_load_variables(&mut slot.params, &spec, id)?;
            }
        }
        self.add_exported_variables(&slot, mesh, apply_ics)?;
        Ok(Some(slot))
    }

    /// `Exported Variable n` of a solver section (AddEquationBasics,
    /// MainUtils.F90:1876-2160): nodal fields on the solver's permutation
    /// (all nodes of its active elements, or every node when it has none),
    /// or `-global` ones; `Exported Variable n Output` and `-nooutput`
    /// control output. Initial conditions are applied when `apply_ics`.
    /// `-elem` fields live on the active elements, `-ip` ones on their
    /// Gauss points (`exported`), both optionally masked and initialised.
    /// Masks of other fields, `-dg` and time derivatives are refused.
    fn add_exported_variables(&mut self, slot: &SolverSlot, mesh: &dyn MeshView, apply_ics: bool) -> Result<(), SolverError> {
        for l in 1.. {
            let key = format!("exported variable {l}");
            let Some(text) = slot.params.string(&key).map(str::to_string) else { break };
            for suffix in ["transient", "calculate velocity", "calculate acceleration"] {
                if slot.params.get(&format!("{key} {suffix}")).is_some() {
                    return Err(SolverError::Unsupported(format!("Exported Variable {l} {suffix}")));
                }
            }
            let mut spec = parse_variable_spec(&text, slot.params.integer(&format!("{key} dofs")))?;
            if let Some(false) = slot.params.logical(&format!("{key} output")) {
                spec.output = false;
            }
            if self.vars.id(&spec.name).is_some() {
                continue;
            }
            let mask = crate::exported::mask_of(&slot.params, &key);
            if mask.is_some() && !matches!(spec.kind, FieldKind::Elemental | FieldKind::Ip) {
                return Err(SolverError::Unsupported(format!(
                    "Exported Variable {l} mask on a {:?} field (MakePermUsingMask / CreateMaskedPerm)",
                    spec.kind
                )));
            }
            let mask = mask.as_ref().map(|(s, k)| (s.as_str(), k.as_str()));
            match spec.kind {
                FieldKind::Global => {
                    let mut v = Variable::new(&spec.name, FieldKind::Global, spec.dofs, None, vec![0.0; spec.dofs]);
                    v.output = spec.output;
                    self.vars.add_with_components(v, &spec.components)?;
                }
                FieldKind::Nodal => {
                    let perm: Arc<[u32]> = match &slot.perm {
                        Some(p) => p.clone(),
                        None => node_order_perm(mesh, &slot.active).perm.into(),
                    };
                    let n = perm.iter().filter(|&&s| s != NO_DOF).count();
                    let mut values = vec![0.0; n * spec.dofs];
                    if apply_ics {
                        self.apply_initial_conditions(slot, &spec, mesh, &perm, &mut values)?;
                    }
                    let mut v = Variable::new(&spec.name, FieldKind::Nodal, spec.dofs, Some(perm), values);
                    v.output = spec.output;
                    self.vars.add_with_components(v, &spec.components)?;
                }
                FieldKind::Elemental | FieldKind::Ip => {
                    let (perm, n) = match (spec.kind, mask) {
                        (FieldKind::Ip, _) => {
                            let offsets = crate::exported::ip_perm(self, mesh, slot, mask)?;
                            let n = *offsets.last().expect("one entry") as usize;
                            (offsets, n)
                        }
                        (_, Some(m)) => crate::exported::element_perm_masked(self, mesh, slot, m)?,
                        (_, None) => crate::exported::element_perm(mesh, slot),
                    };
                    let mut values = vec![0.0; n * spec.dofs];
                    if apply_ics {
                        self.exported_initial_values(mesh, &spec, &perm, &mut values)?;
                    }
                    let mut v = Variable::new(&spec.name, spec.kind, spec.dofs, Some(perm.into()), values);
                    v.output = spec.output;
                    self.vars.add_with_components(v, &spec.components)?;
                }
                k => return Err(SolverError::Unsupported(format!("{k:?} exported variable {}", spec.name))),
            }
        }
        Ok(())
    }

    /// `<Proc>_init` of slot `ix`, if it has a procedure.
    fn init_slot(&mut self, ix: usize, mesh: &dyn MeshView, log: &dyn Logger, comm: &dyn Communicator) -> Result<(), SolverError> {
        if let Some(mut inst) = self.slots[ix].instance.take() {
            let r = inst.init(&mut SolverCtx::new(self, ix, mesh, log, comm));
            self.slots[ix].instance = Some(inst);
            r?;
        }
        Ok(())
    }

    /// The next `RUN` block of the `.sif` (ElmerSolver.F90:396-598): `block`
    /// holds only that block's entries (`elmer-sif`'s `Sif::runs[k]`).
    ///
    /// - `ReloadInputFile`: every entry is added to the model's sections and
    ///   to the running solvers' sections, replacing a same-name entry;
    ///   sections the block introduces are appended.
    /// - `AddSolvers`: `Exec Solver` is read again for every solver; a
    ///   solver that had no section before is set up and its `init` called.
    ///   Initial conditions are not applied again (`SetInitialConditions`
    ///   runs on the first load only), so fields keep their values and
    ///   serve as the next run's initial guess; matrices are kept.
    /// - `InitializeIntervals`: `Time`, `Timestep`, `Timestep Size`,
    ///   `Timestep Interval`, `Coupled Iter`, `Nonlin Iter` and the previous
    ///   step sizes start from zero again.
    ///
    /// Then [`run`](Self::run) runs the block; the reference tally carries
    /// on across blocks until [`crate::finalize_reference`].
    pub fn reload(
        &mut self,
        block: &Model,
        mesh: &dyn MeshView,
        registry: &SolverRegistry,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        fn merge(into: &mut ParamMap, from: &ParamMap) {
            for (k, v) in &from.0 {
                into.0.insert(k.clone(), v.clone());
            }
        }
        fn merge_vec(into: &mut Vec<ParamMap>, from: &[ParamMap]) {
            if into.len() < from.len() {
                into.resize_with(from.len(), ParamMap::new);
            }
            for (a, b) in into.iter_mut().zip(from) {
                merge(a, b);
            }
        }
        let m = &mut self.model;
        merge(&mut m.header, &block.header);
        merge(&mut m.simulation, &block.simulation);
        merge(&mut m.constants, &block.constants);
        merge_vec(&mut m.solvers, &block.solvers);
        merge_vec(&mut m.equations, &block.equations);
        merge_vec(&mut m.bodies, &block.bodies);
        merge_vec(&mut m.materials, &block.materials);
        merge_vec(&mut m.body_forces, &block.body_forces);
        merge_vec(&mut m.initial_conditions, &block.initial_conditions);
        merge_vec(&mut m.boundary_conditions, &block.boundary_conditions);

        let mut added = Vec::new();
        for i in 0..self.model.solvers.len() {
            match self.slot_ix(i + 1) {
                Some(ix) => {
                    let slot = &mut self.slots[ix];
                    if let Some(b) = block.solvers.get(i) {
                        merge(&mut slot.params, b);
                    }
                    slot.exec_when = ExecWhen::from_params(&slot.params)?;
                }
                None => {
                    if let Some(slot) = self.build_slot(i, mesh, registry, false)? {
                        let at = self.slots.partition_point(|s| s.id < slot.id);
                        self.slots.insert(at, slot);
                        added.push(i + 1);
                    }
                }
            }
        }
        for id in added {
            let ix = self.slot_ix(id).expect("inserted above");
            self.init_slot(ix, mesh, log, comm)?;
        }

        for g in ["time", "timestep", "timestep size", "timestep interval", "coupled iter", "nonlin iter"] {
            self.vars.set_scalar(g, 0.0);
        }
        // ElmerSolver.F90:1766: `scan` when any solver has `Scanning Loops`.
        if self.slots.iter().any(|s| s.params.get("scanning loops").is_some()) {
            self.vars.set_scalar("scan", 0.0);
        }
        self.prev_sizes.iter_mut().for_each(|s| *s = 0.0);
        self.publish_live();
        Ok(())
    }

    /// The end of `SetInitialConditions` (ElmerSolver.F90:1862-2052): unless
    /// `Initialize Dirichlet Conditions = False`, every boundary element with
    /// a boundary condition writes the condition's value of each variable
    /// (and each component alias) into the variable at its nodes, overriding
    /// the initial condition there; for a transient second-order solver
    /// `<name> Velocity` and `<name> Acceleration` seed its `PrevValues`
    /// columns 1 and 2 (and 6). Elements go in mesh order, variables in
    /// list order, so a later element wins at a shared node. Normal-tangential
    /// conditions are refused; global variables are skipped.
    pub fn initialize_dirichlet(&mut self, mesh: &dyn MeshView) -> Result<(), SolverError> {
        if self.model.simulation.logical("initialize dirichlet conditions") == Some(false) {
            return Ok(());
        }
        let transient = self.is_transient();
        let names = self.vars.names_in_order();
        let nb = mesh.n_bulk() as u32;
        for e in (nb..nb + mesh.n_boundary() as u32).map(ElemId) {
            let el = mesh.element(e);
            let Some(bc_id) = self.model.bc_of_tag(el.tag) else { continue };
            let bc = &self.model.boundary_conditions[bc_id - 1];
            for (name, handle) in &names {
                let (id, comp) = match *handle {
                    VarHandle::Primary(id) => (id, None),
                    VarHandle::Component { parent, component } => (parent, Some(component)),
                };
                let (dofs, kind) = {
                    let v = self.vars.by_id(id);
                    (v.dofs, v.kind)
                };
                if kind == FieldKind::Global {
                    continue;
                }
                let Some(val) = bc.get(name) else { continue };
                let b = &self.backends;
                let fields = SimFields { vars: &self.vars, mesh };
                if dofs > 1 && comp.is_none() {
                    // ListGetRealArray: the first MIN(rows, dofs) entries of column 1.
                    let c = crate::expr::compile_value(val, &b.udfs);
                    let mut out = vec![0.0; el.nodes.len() * c.ncomp()];
                    c.eval_nodes(el.nodes, &fields, b.matc.as_deref().unwrap_or(&NoMatc), &mut out)?;
                    let (rows, cols, nc) = (c.rows(), c.cols(), c.ncomp());
                    let v = self.vars.by_id_mut(id);
                    for (j, &n) in el.nodes.iter().enumerate() {
                        if let Some(s) = v.slot(n) {
                            for l in 0..rows.min(dofs) {
                                v.values[s * dofs + l] = out[j * nc + l * cols];
                            }
                        }
                    }
                    continue;
                }
                if el.etype.family() != 1 && comp.is_some() {
                    let base = name.rsplit_once(' ').map_or(name.as_str(), |(b, _)| b);
                    if bc.logical(&format!("normal-tangential {base}")) == Some(true) {
                        return Err(SolverError::Unsupported(format!(
                            "Initialize Dirichlet Conditions with Normal-Tangential {base}"
                        )));
                    }
                }
                let xs = crate::expr::eval_at_nodes_with(val, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                let second = transient
                    && self
                        .slots
                        .iter()
                        .any(|s| s.variable == Some(id) && s.time_order == 2);
                let rates = if second {
                    let mut r = [None, None];
                    for (k, suffix) in ["velocity", "acceleration"].into_iter().enumerate() {
                        if let Some(v) = bc.get(&format!("{name} {suffix}")) {
                            r[k] = Some(crate::expr::eval_at_nodes_with(v, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?);
                        }
                    }
                    r
                } else {
                    [None, None]
                };
                let v = self.vars.by_id_mut(id);
                let (stride, offset) = match comp {
                    Some(c) => (dofs, c),
                    None => (1, 0),
                };
                for (j, &n) in el.nodes.iter().enumerate() {
                    let Some(s) = v.slot(n) else { continue };
                    let k = s * stride + offset;
                    v.values[k] = xs[j];
                    if let Some(vel) = &rates[0] {
                        v.prev_values[0][k] = vel[j];
                    }
                    if let Some(acc) = &rates[1] {
                        v.prev_values[1][k] = acc[j];
                        v.prev_values[5][k] = acc[j];
                    }
                }
            }
        }
        self.publish_live();
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
        slot.time_order = p.integer("time derivative order").unwrap_or(1).clamp(1, 2) as u8;
        Ok(())
    }

    fn add_solver_variable(
        &mut self,
        slot: &mut SolverSlot,
        spec: &VarSpec,
        mesh: &dyn MeshView,
        apply_ics: bool,
    ) -> Result<(), SolverError> {
        if spec.kind == FieldKind::Global {
            let mut v = Variable::new(&spec.name, FieldKind::Global, spec.dofs, None, vec![0.0; spec.dofs]);
            v.output = spec.output;
            slot.variable = Some(self.vars.add_with_components(v, &spec.components)?);
            return Ok(());
        }
        if spec.kind != FieldKind::Nodal {
            return Err(SolverError::Unsupported(format!("{:?} solver variable", spec.kind)));
        }
        let layout: Option<Arc<dyn crate::SolverDofs>> = match &self.backends.dofs {
            Some(b) => b
                .for_solver_id(mesh, &self.model, slot.id, &slot.params, &slot.active, spec.dofs)
                .transpose()?
                .map(Arc::from),
            None => None,
        };
        let perm = match (&layout, &self.backends.assembly) {
            (Some(d), _) => d.perm(),
            (None, Some(b)) => b.build_perm_for(mesh, &slot.active, &slot.params)?,
            (None, None) => node_order_perm(mesh, &slot.active),
        };
        slot.dofs = layout;
        let perm_arc: Arc<[u32]> = perm.perm.into();
        let mut values = vec![0.0; perm.n_slots * spec.dofs];
        if apply_ics {
            self.apply_initial_conditions(slot, spec, mesh, &perm_arc, &mut values)?;
        }
        let mut v = Variable::new(&spec.name, FieldKind::Nodal, spec.dofs, Some(perm_arc.clone()), values);
        v.output = spec.output;
        if slot.time_order > 0 {
            // MainUtils.F90:2510-2517: seven columns for a second-order equation.
            let m = if slot.time_order == 2 {
                7
            } else {
                slot.order.max(slot.time_order as usize)
            };
            v.prev_values = vec![vec![0.0; v.values.len()]; m];
        }
        slot.perm = Some(perm_arc);
        slot.variable = Some(self.vars.add_with_components(v, &spec.components)?);
        Ok(())
    }

    /// `Initial Condition` values for the variable (or its components),
    /// evaluated at each active element's nodes.
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
            let fields = self.fields(mesh);
            for &(k, key) in &keys {
                let Some(val) = ic.get(key) else { continue };
                let b = &self.backends;
                let xs = crate::expr::eval_at_nodes_with(val, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                for (&n, &x) in el.nodes.iter().zip(&xs) {
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
        let second = slot.time_order == 2;
        let v = self.vars.by_id_mut(id);
        if v.prev_values.is_empty() {
            return Ok(());
        }
        if second {
            // SolverBasics.F90:2788-2792: x, v, a of the step from the values
            // and the predicted v, a (and the generalized-alpha column).
            let pv = &mut v.prev_values;
            pv[2].copy_from_slice(&v.values);
            let (lo, hi) = pv.split_at_mut(3);
            hi[0].copy_from_slice(&lo[0]);
            hi[1].copy_from_slice(&lo[1]);
            let (lo, hi) = pv.split_at_mut(6);
            hi[0].copy_from_slice(&lo[5]);
            return Ok(());
        }
        for i in (1..order).rev() {
            let (lo, hi) = v.prev_values.split_at_mut(i);
            hi[0].copy_from_slice(&lo[i - 1]);
        }
        v.prev_values[0].copy_from_slice(&v.values);
        Ok(())
    }

    /// `Update2ndOrder` (SolveCore.F90:3317) after a solve of a second-order
    /// equation: the predicted velocity and acceleration columns from the
    /// first `n` new values, with the step's `Solver % dt`.
    pub fn update_2nd_order(&mut self, ix: usize, n: usize) -> Result<(), SolverError> {
        let slot = &self.slots[ix];
        let Some(id) = slot.variable else { return Ok(()) };
        if slot.time_order != 2 {
            return Ok(());
        }
        let scheme = elmer_timeint::Scheme::from_params(&slot.params, &self.model.simulation, 2)?;
        let Some((a, b)) = scheme.second_order_params() else {
            return Ok(());
        };
        let dt = slot.dt;
        let v = self.vars.by_id_mut(id);
        if v.prev_values.len() < 7 {
            return Ok(());
        }
        let mut state = elmer_timeint::SecondOrderState {
            cols: std::array::from_fn(|c| std::mem::take(&mut v.prev_values[c])),
        };
        elmer_timeint::update_2nd_order(dt, &v.values[..n], &mut state, a, b);
        for (c, col) in state.cols.into_iter().enumerate() {
            v.prev_values[c] = col;
        }
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
        let fields = SimFields { vars: &self.vars, mesh };
        let params = crate::expr::resolve_creals(
            &slot.params,
            CHANGE_CREALS,
            &fields,
            self.backends.matc.as_deref(),
            &self.backends.udfs,
        )?;
        let sys = if use_system { slot.system.as_ref() } else { None };
        let env = ChangeEnv {
            params: &params,
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
        let report = compute_change(&env, var, &mut slot.state)?;
        // UpdateDependentObjects(Solver, .FALSE.) (SolverBasics.F90:4242).
        if !steady && self.slots[ix].params.logical("nonlinear update exported variables") == Some(true) {
            self.update_exported_variables(ix, mesh, log)?;
        }
        Ok(Some(report))
    }

    /// `SolverActivate` + `SingleSolver`: exec conditions, then per
    /// `Scanning Loops` round (one without) the `nonlin iter` reset and one
    /// call of the solver's procedure.
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
        // SingleSolver (MainUtils.F90:5709), before the exec condition.
        self.normalize_keywords(mesh, log)?;
        let p = &self.slots[ix].params;
        // ListGetCReal (MainUtils.F90:5722): may follow `time` through MATC.
        if self.creal(p, "exec condition", mesh)?.is_some_and(|st| st < 0.0) {
            return Ok(());
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
            // ListGetConstReal (MainUtils.F90:5796).
            if let Some(v) = p.get("timestep scale") {
                dts = crate::expr::eval_const_with(v, self.backends.matc.as_deref(), &self.backends.udfs)?;
            }
            if p.logical("timestep over intervals") == Some(true) {
                dts = execi as f64;
            }
            if transient {
                time_derivative_active = match p.logical("time derivative active") {
                    Some(b) => b,
                    None => self.creal(p, "time derivative condition", mesh)?.is_none_or(|c| c > 0.0),
                };
            }
        }
        // MainUtils.F90:5840: `Scanning Loops` runs the solver that many
        // times with the global `scan` at 1, 2, ...
        let loops = match p.get("scanning loops") {
            None => None,
            Some(_) => Some(p.integer("scanning loops").ok_or_else(|| {
                SolverError::fatal("SolverActivate", "> Scanning Loops < should be an integer")
            })?),
        };
        if loops.is_some() {
            if p.logical("save scanning modes") == Some(true) {
                return Err(SolverError::Unsupported("Save Scanning Modes".into()));
            }
            if self.vars.scalar("scan").is_none() {
                return Err(SolverError::fatal("SolverActivate", "For scanning we should have scanning variable!"));
            }
        }
        if p.get("assembly solvers").is_some_and(|v| *v != Value::Logical(false)) {
            return Err(SolverError::Unsupported(format!("Solver {}: assembly solvers", self.slots[ix].id)));
        }
        self.slots[ix].dt = dts * dt;
        // MainUtils.F90:1392, 1518: `Block Solver` runs BlockSolver around
        // the procedure; without one the blocks are assembled one by one.
        // Without one, the slot's instance is the block assembly
        // (`block_system`).
        let p = &self.slots[ix].params;
        let with_procedure = p.get("procedure").is_some();
        // The solver mode, the last of these set (MainUtils.F90:1386-1395).
        let mode = ["coupled solver", "block solver", "assembly solver"]
            .into_iter()
            .rfind(|k| p.logical(k) == Some(true));
        let block = with_procedure && mode == Some("block solver");
        let coupled = with_procedure && matches!(mode, Some("coupled solver") | Some("assembly solver"));
        for scan in 1..=loops.unwrap_or(1) {
            if !auxiliary {
                self.vars.set_scalar("nonlin iter", 1.0);
            }
            if loops.is_some() {
                self.vars.set_scalar("scan", f64::from(scan));
            }
            self.publish_live();
            let slot = &mut self.slots[ix];
            (slot.started, slot.finished) = (false, false);
            if let Some(mut inst) = self.slots[ix].instance.take() {
                let r = if block {
                    crate::block_driver::block_solver(self, ix, &mut inst, dts * dt, time_derivative_active, mesh, log, comm)
                } else if coupled {
                    crate::block_driver::coupled_with_procedure(
                        self,
                        ix,
                        &mut inst,
                        dts * dt,
                        time_derivative_active,
                        mesh,
                        log,
                        comm,
                    )
                } else {
                    inst.solve(
                        &mut SolverCtx::new(self, ix, mesh, log, comm),
                        dts * dt,
                        time_derivative_active,
                    )
                };
                self.slots[ix].instance = Some(inst);
                r?;
            }
            // `Pre Solvers` and `Post Solvers` run from DefaultStart and
            // DefaultFinish; a procedure that calls neither would skip them.
            let slot = &self.slots[ix];
            for (key, ran, routine) in [
                ("pre solvers", slot.started, "DefaultStart"),
                ("post solvers", slot.finished, "DefaultFinish"),
            ] {
                if !ran && slot.params.get(key).is_some_and(|v| *v != Value::Logical(false)) {
                    return Err(SolverError::Unsupported(format!(
                        "Solver {}: {key} with a procedure that does not call {routine}",
                        slot.id
                    )));
                }
            }
            // UpdateDependentObjects(Solver, .TRUE.) (MainUtils.F90:5625).
            if self.slots[ix].params.logical("update exported variables") == Some(true) {
                self.update_exported_variables(ix, mesh, log)?;
            }
            self.slots[ix].times_visited += 1;
        }
        Ok(())
    }

    /// `DefaultSlaveSolvers(Solver, name)` (DefUtils.F90:3750): every solver
    /// listed in slot `master`'s `key` (`Slave Solvers`, `Nonlinear Pre
    /// Solvers`, `Nonlinear Post Solvers`) is activated with the current
    /// `Timestep Size` and the simulation's transient flag; `nonlin iter` is
    /// restored afterwards. `name` is the keyword as Elmer prints it.
    pub(crate) fn run_slave_solvers(
        &mut self,
        master: usize,
        name: &str,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let ids = match self.slots[master].params.get(&normalize_key(name)) {
            Some(Value::Integer(v)) => v.clone(),
            _ => return Ok(()),
        };
        log.info("DefaultSlaveSolvers", &format!("Executing slave solvers: {name}"), 6);
        let dt = self.vars.scalar("timestep size").unwrap_or(0.0);
        let transient = self.is_transient();
        let iter = self.vars.scalar("nonlin iter");
        for k in ids {
            let ix = usize::try_from(k)
                .ok()
                .and_then(|k| self.slot_ix(k))
                .ok_or_else(|| SolverError::fatal("DefaultSlaveSolvers", format!("no solver {k} to call as a slave")))?;
            if ix == master {
                return Err(SolverError::fatal("DefaultSlaveSolvers", format!("solver {k} is its own slave")));
            }
            log.info("DefaultSlaveSolvers", &format!("Calling slave solver: {k}"), 8);
            if self.slots[master].params.logical("monolithic slave") == Some(true)
                && self.slots[ix].params.get("linear system solver disabled").is_none()
            {
                log.info("DefaultSlaveSolvers", &format!("Disabling linear system solver for slave: {k}"), 6);
                self.slots[ix].params.insert("Linear System Solver Disabled", Value::Logical(true));
            }
            self.solver_activate(ix, dt, transient, mesh, log, comm)?;
        }
        if let Some(it) = iter {
            self.vars.set_scalar("nonlin iter", it.round());
        }
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
        let all = Phases { before: true, at: true, after: true };
        self.solve_equations_in(dt, transient, coupled_min, coupled_max, real_timestep, all, mesh, log, comm)
    }

    /// [`solve_equations`](Self::solve_equations) with `BeforeTime`,
    /// `AtTime` and `AfterTime` selecting which solver groups run. The
    /// histories advance (`InitializeTimestep`) on every call of a transient
    /// run, whatever the phases, as in the Fortran.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_equations_in(
        &mut self,
        dt: f64,
        transient: bool,
        coupled_min: i32,
        coupled_max: i32,
        real_timestep: Option<usize>,
        phases: Phases,
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

        if phases.before {
            log.info("SolveEquations", "Solvers before timestep", 12);
            for ix in 0..self.slots.len() {
                let s = &self.slots[ix];
                if s.has_procedure() && s.exec_when == ExecWhen::BeforeTimestep && !skip_first(s) {
                    self.solver_activate(ix, dt, transient, mesh, log, comm)?;
                    comm.barrier();
                }
            }
        }
        let mut steady_reached = false;
        if phases.at {
            steady_reached = self.solve_at_time(dt, transient, scanning, cmin, cmax, real_timestep, mesh, log, comm)?;
        }
        if phases.after {
            log.info("SolveEquations", "Solvers after timestep", 12);
            for ix in 0..self.slots.len() {
                let s = &self.slots[ix];
                if s.has_procedure() && s.exec_when == ExecWhen::AfterTimestep && !skip_first(s) {
                    self.solver_activate(ix, dt, transient, mesh, log, comm)?;
                    comm.barrier();
                }
            }
        }
        Ok(steady_reached)
    }

    /// The `AtTime` part of `SolveEquations`: the coupled loop over the
    /// `Always` solvers. Returns `SteadyStateReached` (steady runs only).
    #[allow(clippy::too_many_arguments)]
    fn solve_at_time(
        &mut self,
        dt: f64,
        transient: bool,
        scanning: bool,
        cmin: i32,
        cmax: i32,
        real_timestep: Option<usize>,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<bool, SolverError> {

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
        Ok(!transient && done.iter().all(|&d| d))
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
                    if let Some(c) = self.creal(p, "steady state condition", mesh)? {
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
    /// A scanning run steps like a transient one (`Timestep Intervals`,
    /// unit sizes by default, `time` advancing, no steady-state exit) but
    /// its solvers run steady. Periodic time is refused.
    pub fn run(&mut self, mesh: &dyn MeshView, log: &dyn Logger, comm: &dyn Communicator) -> Result<RunReport, SolverError> {
        let simp = self.model.simulation.clone();
        let transient = self.is_transient();
        let scanning = self.is_scanning();
        let timed = transient || scanning;
        for k in ["periodic time", "predictor-corrector control", "parallel timestepping"] {
            if simp.get(k).is_some_and(|v| *v != Value::Logical(false)) {
                return Err(SolverError::Unsupported(format!("Simulation keyword {k}")));
            }
        }
        let real_array = |k: &str| match simp.get(k) {
            Some(Value::Real { data, .. }) => Some(data.clone()),
            _ => None,
        };
        let (intervals, sizes, coupled_max) = if timed {
            let iv = match simp.get("timestep intervals") {
                Some(Value::Integer(v)) => v.clone(),
                _ => {
                    return Err(SolverError::fatal(
                        "ElmerSolver",
                        "Keyword > Timestep Intervals < MUST be defined for transient and scanning simulations",
                    ))
                }
            };
            // InitializeIntervals: unit sizes when a Timestep Size function gives them.
            let sz = match real_array("timestep sizes") {
                Some(d) => d,
                None if scanning || simp.get("timestep size").is_some() => vec![1.0; iv.len()],
                None => {
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
        let ratios = real_array("timestep ratios");
        let adaptive = transient && simp.logical("adaptive timestepping") == Some(true);
        let divergence_control = transient && !adaptive && simp.logical("convergence control") == Some(true);
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
        let mut dt = 1.0;
        let mut adapt = AdaptiveState::default();
        'outer: for (interval, &nsteps) in intervals.iter().enumerate() {
            let mut real_timestep = 1usize;
            for timestep in 1..=nsteps {
                cum += 1;
                self.vars.set_scalar("timestep", cum as f64);
                dt = self.step_size(&simp, timed, cum, interval, timestep, &sizes, ratios.as_deref(), dt, mesh, log)?;
                let mut t = self.vars.scalar("time").unwrap_or(0.0);
                if cum == 1 && simp.logical("timestep start zero") == Some(true) {
                    log.info("ExecSimulation", "Not advancing the 1st timestep!", 4);
                } else {
                    t += dt;
                }
                self.vars.set_scalar("time", t);
                if timestep > 1 || interval > 0 {
                    self.prev_sizes.pop();
                    let old = self.vars.scalar("timestep size").unwrap_or(0.0);
                    self.prev_sizes.insert(0, old);
                }
                self.vars.set_scalar("timestep size", dt);
                self.vars.set_scalar("timestep interval", (interval + 1) as f64);
                if !transient {
                    let c = self.vars.scalar("coupled iter").unwrap_or(0.0);
                    self.vars.set_scalar("coupled iter", c + 1.0);
                }
                if comm.rank() == 0 {
                    log.info("MAIN", " ", 3);
                    log.info("MAIN", "-------------------------------------", 3);
                    if timed {
                        log.info("MAIN", &format!("Time: {}/{}:{}", cum, stepcount, es(t, 12, 3)), 3);
                    } else {
                        // WRITE(Message,*) 'Steady state iteration: ',cum_Timestep
                        log.info("MAIN", &format!(" Steady state iteration: {}", list_int(cum)), 3);
                    }
                    log.info("MAIN", "-------------------------------------", 3);
                    log.info("MAIN", " ", 3);
                }
                let reached = if adaptive {
                    let rt = &mut real_timestep;
                    self.adaptive_step(&simp, dt, coupled_min, coupled_max, rt, &mut adapt, mesh, log, comm)?;
                    false
                } else if divergence_control {
                    let rt = &mut real_timestep;
                    self.divergence_step(&simp, dt, coupled_min, coupled_max, rt, &mut adapt, mesh, log, comm)?;
                    false
                } else {
                    let r = self.solve_equations(
                        dt,
                        transient,
                        coupled_min,
                        coupled_max,
                        Some(real_timestep),
                        mesh,
                        log,
                        comm,
                    )?;
                    real_timestep += 1;
                    r
                };
                report.steps = cum as usize;
                report.steady_state_reached = reached;
                let oi = output_intervals[interval];
                if oi != 0 && ((timestep - 1) % oi == 0 || reached) {
                    self.exec_slots(ExecWhen::BeforeSaving, dt, transient, mesh, log, comm)?;
                    self.exec_slots(ExecWhen::AfterSaving, dt, transient, mesh, log, comm)?;
                }
                if self.creal(&simp, "exit condition", mesh)?.is_some_and(|c| c > 0.0) {
                    log.info("MAIN", "Found a positive exit condition, exiting...", 3);
                    break 'outer;
                }
                if reached && !timed && timestep >= coupled_min {
                    break 'outer;
                }
            }
        }

        self.exec_slots(ExecWhen::AfterAll, dt, transient, mesh, log, comm)?;
        crate::reference::compare_to_reference(self, mesh, log, comm)?;
        Ok(report)
    }

    /// The step size of step `timestep` of `interval` (ElmerSolver.F90:2934-2985):
    /// `First Timestep Size` on the first step; then for a transient or scanning run
    /// `Timestep Function` (obsolete) or `Timestep Size` (`ListGetCReal`,
    /// so it may follow `time`), re-evaluated `Timestep Size Iterations`
    /// times with `time` moved ahead by the last value; else `Timestep
    /// Sizes(interval)`, times `Timestep Ratios(interval)^(timestep-1)`.
    /// A steady run keeps `dt`.
    #[allow(clippy::too_many_arguments)]
    fn step_size(
        &mut self,
        simp: &ParamMap,
        transient: bool,
        cum: i32,
        interval: usize,
        timestep: i32,
        sizes: &[f64],
        ratios: Option<&[f64]>,
        dt: f64,
        mesh: &dyn MeshView,
        log: &dyn Logger,
    ) -> Result<f64, SolverError> {
        if cum == 1 {
            if let Some(x) = self.creal(simp, "first timestep size", mesh)? {
                return Ok(x);
            }
        }
        if !transient {
            return Ok(dt);
        }
        let mut f = self.creal(simp, "timestep function", mesh)?;
        if f.is_some() {
            log.warn("ExecSimulation", "Obsolete keyword > Timestep Function < , use > Timestep Size < instead");
        } else {
            f = self.creal(simp, "timestep size", mesh)?;
        }
        let Some(mut f) = f else {
            let mut dt = sizes[interval];
            if let Some(q) = ratios.and_then(|r| r.get(interval)) {
                if (1.0 - q).abs() > f64::EPSILON {
                    dt *= q.powi(timestep - 1);
                }
            }
            return Ok(dt);
        };
        let iters = simp.integer("timestep size iterations").unwrap_or(0);
        if iters > 0 {
            let t0 = self.vars.scalar("time").unwrap_or(0.0);
            for _ in 0..iters {
                self.vars.set_scalar("time", t0 + f);
                f = self
                    .creal(simp, "timestep size", mesh)?
                    .ok_or_else(|| SolverError::fatal("ExecSimulation", "Timestep Size Iterations needs > Timestep Size <"))?;
            }
            self.vars.set_scalar("time", t0);
        }
        if f < f64::EPSILON {
            return Err(SolverError::fatal(
                "ExecSimulation",
                format!("Timestep smaller than epsilon: {}", es(f, 12, 3)),
            ));
        }
        Ok(f)
    }

    /// One step of `Adaptive Timestepping` (ElmerSolver.F90:3181-3330):
    /// substeps of size `ddt` cover `dt`. Each substep is solved once with
    /// `ddt` and again as two halves from the same start; the relative
    /// change of the variable norms (or of `Adaptive Error Measure`)
    /// between the two decides: below `Adaptive Time Error` (or at
    /// `Adaptive Min Timestep`) the halves are kept, the after-timestep
    /// solvers run and `ddt` may grow by `Adaptive Increase Coefficient`
    /// after `Adaptive Keep Smallest` accepted substeps; otherwise the
    /// values, norms and histories go back and `ddt` shrinks by `Adaptive
    /// Decrease Coefficient`. As in the Fortran every `SolveEquations` call
    /// advances the histories, and `Timestep` holds `ddt/2` from the halves
    /// on.
    #[allow(clippy::too_many_arguments)]
    fn adaptive_step(
        &mut self,
        simp: &ParamMap,
        dt: f64,
        cmin: i32,
        cmax: i32,
        real_timestep: &mut usize,
        st: &mut AdaptiveState,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let limit = simp
            .real("adaptive time error")
            .ok_or_else(|| SolverError::fatal("MAIN", "Adaptive Time Error must be given for adaptive stepping scheme."))?;
        let keep_smallest = simp.integer("adaptive keep smallest").unwrap_or(0).max(0);
        let max_dt = simp.real("adaptive max timestep").map_or(dt, |m| m.min(dt));
        let min_dt = simp.real("adaptive min timestep").unwrap_or(1.0e-8 * max_dt);
        let increase = simp.real("adaptive increase coefficient").unwrap_or(2.0);
        let decrease = simp.real("adaptive decrease coefficient").unwrap_or(0.5);
        let rough = simp.logical("adaptive rough timestep") == Some(true);
        let smart = simp.logical("adaptive smart timestep") == Some(true);
        let measure = |s: &Self| s.model.simulation.real("adaptive error measure");
        let solvers: Vec<VarId> = self.slots.iter().filter_map(|s| s.variable).collect();
        let at = Phases { before: true, at: true, after: false };
        let after = Phases { before: false, at: false, after: true };

        let mut cum_time = 0.0;
        if st.ddt < 0.0 || st.ddt > max_dt {
            st.ddt = max_dt;
        }
        let s0 = self.vars.scalar("time").unwrap_or(0.0) - dt;
        let mut smallest = 0;
        while cum_time < dt - 1.0e-12 {
            if !rough {
                st.ddt = st.ddt.min(dt - cum_time);
                if smart && dt - cum_time - st.ddt > 1.0e-12 {
                    log.info("ExecSimulation", "Splitted timestep into two equal parts", 12);
                    st.ddt = st.ddt.min((dt - cum_time) / 2.0);
                }
            }
            let ddt = st.ddt;
            let saved: Vec<(Vec<f64>, f64, Vec<Vec<f64>>)> = solvers
                .iter()
                .map(|&id| {
                    let v = self.vars.by_id(id);
                    (v.values.clone(), v.conv.norm, v.prev_values.clone())
                })
                .collect();
            self.vars.set_scalar("time", s0 + cum_time + ddt);
            self.vars.set_scalar("timestep size", ddt);
            let rt = Some(*real_timestep);
            self.solve_equations_in(ddt, true, cmin, cmax, rt, at, mesh, log, comm)?;
            let measure_full = measure(self);
            let mut full_norms = Vec::with_capacity(solvers.len());
            for (&id, (vals, _, prev)) in solvers.iter().zip(&saved) {
                let v = self.vars.by_id_mut(id);
                v.values.clone_from(vals);
                full_norms.push(v.conv.norm);
                v.prev_values.clone_from(prev);
            }
            self.vars.set_scalar("timestep", ddt / 2.0);
            self.vars.set_scalar("time", s0 + cum_time + ddt / 2.0);
            self.solve_equations_in(ddt / 2.0, true, cmin, cmax, rt, at, mesh, log, comm)?;
            self.vars.set_scalar("time", s0 + cum_time + ddt);
            self.solve_equations_in(ddt / 2.0, true, cmin, cmax, rt, at, mesh, log, comm)?;
            let max_err = match (measure_full, measure(self)) {
                (Some(a), Some(b)) => (a - b).abs(),
                _ => {
                    let mut e = 0.0_f64;
                    for (&id, &full) in solvers.iter().zip(&full_norms) {
                        let n = self.vars.by_id(id).conv.norm;
                        if full != n {
                            e = e.max((full - n).abs() / full);
                        }
                    }
                    e
                }
            };
            if max_err < limit || ddt <= min_dt {
                cum_time += ddt;
                *real_timestep += 1;
                if smallest >= keep_smallest || st.step_control > 0 {
                    st.ddt = (increase * ddt).min(max_dt);
                    st.step_control = 1;
                    smallest = 0;
                } else {
                    st.step_control = 0;
                    smallest += 1;
                }
                let rt = Some(*real_timestep);
                self.solve_equations_in(st.ddt / 2.0, true, cmin, cmax, rt, after, mesh, log, comm)?;
            } else {
                for (&id, (vals, norm, prev)) in solvers.iter().zip(&saved) {
                    let v = self.vars.by_id_mut(id);
                    v.conv.norm = *norm;
                    v.values.clone_from(vals);
                    v.prev_values.clone_from(prev);
                }
                st.ddt = decrease * ddt;
                st.step_control = -1;
            }
            // WRITE(*,'(a,3e20.12)') 'Adaptive(cum,ddt,err): ', cumtime, ddt, maxerr
            log.info(
                "ExecSimulation",
                &format!(
                    "Adaptive(cum,ddt,err): {}{}{}",
                    e_fmt(cum_time, 20, 12),
                    e_fmt(st.ddt, 20, 12),
                    e_fmt(max_err, 20, 12)
                ),
                4,
            );
        }
        self.vars.set_scalar("timestep size", dt);
        self.vars.set_scalar("time", s0 + dt);
        Ok(())
    }

    /// One step of `Convergence Control` (ElmerSolver.F90:3333-3395):
    /// substeps of size `ddt` (grown by `Adaptive Increase Coefficient`
    /// on entry, capped by `Adaptive Max Timestep`) cover `dt`. A substep
    /// whose nonlinear iteration of some solver diverged (`NonlinConverged >
    /// 1`, `Nonlinear System Divergence Limit`) is repeated with `ddt` shrunk
    /// by `Adaptive Decrease Coefficient` (not below `Adaptive Min Timestep`,
    /// where the run stops), from each variable's last committed values;
    /// otherwise the after-timestep solvers run and `ddt` grows again.
    #[allow(clippy::too_many_arguments)]
    fn divergence_step(
        &mut self,
        simp: &ParamMap,
        dt: f64,
        cmin: i32,
        cmax: i32,
        real_timestep: &mut usize,
        st: &mut AdaptiveState,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let max_dt = simp.real("adaptive max timestep").map_or(dt, |m| m.min(dt));
        let min_dt = simp.real("adaptive min timestep").unwrap_or(1.0e-8 * max_dt);
        let increase = simp.real("adaptive increase coefficient").unwrap_or(2.0);
        let decrease = simp.real("adaptive decrease coefficient").unwrap_or(0.5);
        let rough = simp.logical("adaptive rough timestep") == Some(true);
        let smart = simp.logical("adaptive smart timestep") == Some(true);
        let at = Phases { before: true, at: true, after: false };
        let after = Phases { before: false, at: false, after: true };

        log.info("ExecSimulation", "Solving equations with divergence control", 6);
        let mut cum_time = 0.0;
        st.ddt *= increase;
        if st.ddt < 0.0 || st.ddt > max_dt {
            st.ddt = max_dt;
        }
        let s0 = self.vars.scalar("time").unwrap_or(0.0) - dt;
        while cum_time < dt - 1.0e-12 {
            if !rough {
                st.ddt = st.ddt.min(dt - cum_time);
                if smart && dt - cum_time - st.ddt > 1.0e-12 {
                    log.info("ExecSimulation", "Splitted timestep into two equal parts", 12);
                    st.ddt = st.ddt.min((dt - cum_time) / 2.0);
                }
            }
            let ddt = st.ddt;
            self.vars.set_scalar("time", s0 + cum_time + ddt);
            self.vars.set_scalar("timestep size", ddt);
            self.solve_equations_in(ddt, true, cmin, cmax, Some(*real_timestep), at, mesh, log, comm)?;
            let diverged = self.slots.iter().enumerate().find(|(_, s)| {
                s.variable.is_some_and(|id| self.vars.by_id(id).conv.nonlin_converged > 1)
            });
            match diverged {
                None => {
                    log.info("ExecSimulation", "No solver has diverged", 8);
                    self.solve_equations_in(ddt, true, cmin, cmax, Some(*real_timestep), after, mesh, log, comm)?;
                    cum_time += ddt;
                    *real_timestep += 1;
                    st.ddt = increase * ddt;
                    if st.ddt > max_dt {
                        st.ddt = max_dt;
                        st.step_control = 0;
                    }
                }
                Some((ix, _)) => {
                    log.info("ExecSimulation", &format!("Solver {} has diverged", ix + 1), 8);
                    if ddt < min_dt * (1.0 + 1.0e-8) {
                        return Err(SolverError::fatal(
                            "ExecSimulation",
                            "Could not find stable timestep above given minimum",
                        ));
                    }
                    log.info("ExecSimulation", "Reducing timestep due to divergence problems!", 6);
                    st.ddt = (decrease * ddt).max(min_dt);
                    st.step_control = 1;
                    log.info("ExecSimulation", "Reverting to previous timestep as initial guess", 8);
                    let ids: Vec<VarId> = self.slots.iter().filter_map(|s| s.variable).collect();
                    for id in ids {
                        let v = self.vars.by_id_mut(id);
                        if let Some(p) = v.prev_values.first() {
                            let p = p.clone();
                            v.values.copy_from_slice(&p);
                            v.conv.norm = v.conv.prev_norm;
                        }
                    }
                }
            }
        }
        Ok(())
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
        // Component names with spaces (MagnetoDynamics2DHarmonic's default).
        let s = parse_variable_spec("Potential[Potential Re:1 Potential Im:1]", None).unwrap();
        assert_eq!((s.name.as_str(), s.dofs), ("potential", 2));
        assert_eq!(s.components, vec!["potential re", "potential im"]);
        let s = parse_variable_spec("-dofs 4 E[E Field Re:2  E Field Im:2]", None).unwrap();
        assert_eq!(s.components, vec!["e field re 1", "e field re 2", "e field im 1", "e field im 2"]);
        assert!(parse_variable_spec("P[A:1 B]", None).is_err());
        assert!(parse_variable_spec("P[:1]", None).is_err());
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
