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
    /// `predictor-corrector` (`SOLVER_EXEC_PREDCORR`): the predictor with
    /// the timestep's ahead-of-time solvers, the corrector with the
    /// after-time ones.
    PredCorr,
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
                "predictor-corrector" => ExecWhen::PredCorr,
                "when created" | "after control" => {
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
        } else if flag("predictor-corrector") {
            ExecWhen::PredCorr
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
    /// `Matrix % MassValues` on the matrix's pattern, kept for `Eigen
    /// Analysis` (`DefaultUpdateMass`); `None` otherwise.
    pub mass: Option<CsrMatrix<f64>>,
    /// `Matrix % DampValues` on the matrix's pattern: the global damping
    /// of a second-order `Ode Matrix` system; `None` otherwise.
    pub damp: Option<CsrMatrix<f64>>,
}

/// One `Solver n` section at run time (Elmer `Solver_t`).
/// `Solver % LocalSystem` (`Local Matrix Storage`): the last local system
/// of each active bulk element, as the solver assembled it and before
/// `Matrix Multiplier Name` scales it (`UseLocalMatrixStorage`).
#[derive(Clone, Debug, Default)]
pub struct LocalStore {
    /// Element index -> active index (`Solver % InvActiveElements`).
    pub inv: std::collections::HashMap<u32, usize>,
    /// `(K, F)` by active index, `K` row-major.
    pub systems: Vec<Option<(Vec<f64>, Vec<f64>)>>,
}

impl LocalStore {
    pub fn new(active: &[ElemId]) -> Self {
        LocalStore {
            inv: active.iter().enumerate().map(|(i, e)| (e.0, i)).collect(),
            systems: vec![None; active.len()],
        }
    }

    /// `GetLocalMatrixStorage`: the stored `(K, F)` of element `e`.
    pub fn get(&self, e: ElemId) -> Option<(&[f64], &[f64])> {
        let i = *self.inv.get(&e.0)?;
        self.systems[i].as_ref().map(|(k, f)| (k.as_slice(), f.as_slice()))
    }
}

pub struct SolverSlot {
    /// 1-based Elmer solver id.
    pub id: usize,
    pub params: ParamMap,
    pub exec_when: ExecWhen,
    pub variable: Option<VarId>,
    /// `Norm` of the nameless null variable `AddEquationBasics` gives a
    /// solver without `Variable` (MainUtils.F90:1504): what
    /// `CompareToReferenceSolution` reads for it. A procedure may set it, as
    /// `ResultOutputSolver` does with its `Reference Values` pseudonorm.
    pub null_norm: f64,
    /// `Solver % NOFEigenValues`: `Eigen System Values` of an `Eigen
    /// Analysis` solver (MainUtils.F90:2575), 0 otherwise.
    pub nof_eigen: usize,
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
    /// `Element % PropertyData` of the Adams-Bashforth / Adams-Moulton
    /// local schemes: per element, the residual and force of the last
    /// predictor.
    pub adams: std::collections::HashMap<u32, Vec<f64>>,
    pub times_visited: usize,
    /// Linear solves of the solver's matrix so far (`Matrix % SolveCount`).
    pub(crate) linear_solves: usize,
    /// `Solver % dt`: the step size of the current visit (`Timestep Scale`
    /// applied).
    pub dt: f64,
    pub state: SolverState,
    pub system: Option<LinearSystem>,
    pub(crate) instance: Option<Box<dyn Solver>>,
    pub(crate) asm: crate::loads::AssemblyState,
    /// `A % BulkValues`, `A % BulkRHS` as the matrix keeps them across
    /// assemblies: the last `DefaultFinishBulkAssembly` snapshot (`Constant
    /// Bulk System` and friends reuse it).
    pub(crate) saved_bulk: Option<(Vec<f64>, Vec<f64>)>,
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
    /// `Variable % LinConverged` of the last linear solve.
    pub(crate) lin_converged: bool,
    /// The coarser meshes of geometric multigrid, coarsest first
    /// ([`Simulation::set_multigrid_meshes`]).
    pub(crate) mg_meshes: Option<Vec<Box<dyn MeshView>>>,
    /// The constraint modes `DefaultDirichletBCs` set up.
    pub(crate) cmodes: Option<crate::constraint_modes::Modes>,
    /// The CutFEM state of a `CutFEM` solver (`CutFEMUtils` module data).
    pub(crate) cutfem: Option<Box<crate::cutfem::CutFemRun>>,
    /// `Var % Upper/LowerLimitActive` of `Apply Limiter`.
    pub(crate) limiter: Option<crate::limiter::LimiterState>,
    /// `Solver % LocalSystem` of `Local Matrix Storage`.
    pub local_store: Option<LocalStore>,
    /// The child solvers' systems (`CreateChildSolver`).
    pub(crate) children: Vec<crate::child::ChildSystem>,
    /// `Solver % NormalTangential`, set up by `DefaultInitialize`.
    pub(crate) nt: Option<crate::nt::NtSystem>,
    /// `ASSOCIATED(Solver % BlockMatrix)`: a block solve of this solver ran
    /// ([`crate::block_alias`]).
    pub(crate) block_matrix: bool,
    /// `Solver % PeriodicFlipActive`: the solver numbers its dofs by the
    /// conforming map and some of its nodes are flipped.
    pub(crate) conforming_flip: bool,
    /// Ownership of the system's rows when the run gives it per row (a DG
    /// solver's, owned by the element's partition), for the consistent norm.
    pub owned_rows: Option<Arc<[bool]>>,
    /// The slot's mesh: 0 for the run's first mesh, else the solver mesh
    /// of its `Mesh` keyword ([`crate::multimesh`]).
    pub mesh_ix: usize,
    /// `Calculate Weights` / `Calculate Boundary Weights` done on this
    /// mesh (SolverActivate, once per active-element table).
    pub(crate) weights_done: bool,
}

impl SolverSlot {
    /// The solver's non-nodal dofs (`None` for a nodal solver): what
    /// `Var % Solver` gives `EvaluateVariableAtGivenPoint`.
    pub fn solver_dofs(&self) -> Option<&dyn crate::SolverDofs> {
        self.dofs.as_deref()
    }

    /// A slot with no variable or procedure; tests and the CLI use it to
    /// build slots by hand.
    pub fn new(id: usize, params: ParamMap) -> Result<Self, SolverError> {
        let exec_when = ExecWhen::from_params(&params)?;
        Ok(Self {
            id,
            params,
            exec_when,
            variable: None,
            null_norm: 0.0,
            nof_eigen: 0,
            perm: None,
            active: Vec::new(),
            time_order: 0,
            order: 1,
            done_time: 0,
            adams: Default::default(),
            times_visited: 0,
            linear_solves: 0,
            dt: 0.0,
            state: SolverState::default(),
            system: None,
            instance: None,
            asm: Default::default(),
            saved_bulk: None,
            constraints: None,
            own_constraints: None,
            dofs: None,
            robin_projectors: Default::default(),
            started: false,
            finished: false,
            lin_converged: true,
            mg_meshes: None,
            cmodes: None,
            cutfem: None,
            limiter: None,
            local_store: None,
            children: Vec::new(),
            nt: None,
            block_matrix: false,
            conforming_flip: false,
            mesh_ix: 0,
            weights_done: false,
            owned_rows: None,
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
    /// `PredictorCorrectorControl`'s SAVEd `dtOld`, `etaOld`, `zeta`.
    pub pred_corr: (f64, f64, f64),
    /// The model view handed to user functions, once asked for.
    live: Option<Arc<LiveModel>>,
    /// `SharedWarned` of `CalculateLoads`: warn about shared boundary dofs once.
    pub(crate) flux_shared_warned: bool,
    /// `SetNormalizedKeywords` has run (`Mesh % EntityWeightsComputed`).
    pub(crate) keywords_normalized: bool,
    /// The `Output File` writer ([`Simulation::set_save_current`]).
    pub(crate) save_current: Option<Box<dyn SaveCurrent>>,
    /// Whether this rank owns each mesh node ([`Simulation::set_node_owned`]).
    node_owned: Option<Arc<[bool]>>,
    /// Elements of another partition (`Element % PartIndex /= myPE`): the
    /// halo of an `ElmerGrid -halo` partitioning ([`Simulation::set_halo_elements`]).
    halo: Option<Arc<[bool]>>,
    /// Boundary elements of other partitions imported for the radiation
    /// factors ([`Simulation::set_foreign_boundary`]).
    foreign: Option<Arc<[bool]>>,
    /// `GElementIndex` of every element ([`Simulation::set_element_gids`]).
    element_gids: Option<Arc<[u64]>>,
    /// `GlobalDOFs` of every node ([`Simulation::set_node_gids`]).
    node_gids: Option<Arc<[u64]>>,
    /// ARPACK's saved state (`dgetv0`'s seed), kept across eigen solves as
    /// in one ElmerSolver process.
    pub(crate) arpack: elmer_eigen::arpack::saupd::Arpack,
    /// The `SAVE`d state of ARPACK's complex routines, for the process.
    pub(crate) zarpack: elmer_eigen::arpack::znaupd::ZArpack,
    /// Meshes a CutFEM solver adds to the output for the step (the
    /// isoline mesh, the cut mesh), saved after the run's own.
    pub extra_meshes: Vec<crate::cutfem::ExtraMesh>,
    /// The solver BlockSolve.F90's `TotMatrix` points at
    /// ([`crate::block_alias`]).
    pub(crate) block_tot: Option<crate::block_alias::BlockOwner>,
    /// `Element % DGIndexes` of the mesh, once a DG variable asked for them
    /// ([`crate::exported::DgIndexes`]).
    pub dg_indexes: Option<crate::exported::DgIndexes>,
    /// `Mesh % Nodes` of a p mesh past `NumberOfNodes`: the coordinates of
    /// the p dofs by their global index (`mGetElementDOFs`), zero unless
    /// set (`EnlargeCoordinates`). `Follow P Curvature` fits them to a
    /// curved boundary ([`crate::follow_p::follow_p_curvature`]) and a p
    /// MeshSolve adds its displacement to them (`DisplaceMesh`).
    pub p_coords: Option<Vec<[f64; 3]>>,
    /// The mesh-wide p dofs of the p solvers ([`SolverDofs::ndeg_parts`],
    /// their maxima): what makes a nodal solver's `Ndeg` exceed the nodes.
    pub(crate) p_ndeg_parts: Option<(usize, usize)>,
    /// Every element's nodes in a p layout's order (p tetrahedra are
    /// reordered in `Mesh % Elements` itself), when it differs from the
    /// mesh's: the order a nodal solver's first-touch numbering sees.
    pub(crate) p_element_nodes: Option<Arc<Vec<Vec<u32>>>>,
    /// Adaptive mesh refinement: the pending mesh change and the current
    /// mesh's state ([`crate::remesh`]).
    pub(crate) adaptive: crate::remesh::AdaptiveMeshes,
    /// Solver-specific meshes and their variable lists ([`crate::multimesh`]).
    pub(crate) meshes: crate::multimesh::Meshes,
}

/// `mesh` with its elements' nodes in a p layout's order.
struct PNodeOrder<'a> {
    mesh: &'a dyn MeshView,
    nodes: &'a [Vec<u32>],
}

impl MeshView for PNodeOrder<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: ElemId) -> elmer_core::ElementRef<'_> {
        let mut el = self.mesh.element(e);
        el.nodes = &self.nodes[e.0 as usize];
        el
    }
}

/// `SaveCurrent` (ElmerSolver.F90:3700): the `Output File` writer, called
/// at each save point between the `before saving` and `after saving`
/// solvers. It lives outside the solver crate, like the VTU writer.
pub trait SaveCurrent {
    fn save(&mut self, sim: &Simulation, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError>;
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
    "derivative eps",
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
            dg_indexes: None,
            p_coords: None,
            prev_sizes: vec![0.0; 5],
            pred_corr: (0.0, 0.0, 1.0),
            live: None,
            flux_shared_warned: false,
            keywords_normalized: false,
            block_tot: None,
            save_current: None,
            node_owned: None,
            halo: None,
            foreign: None,
            element_gids: None,
            node_gids: None,
            arpack: Default::default(),
            zarpack: Default::default(),
            extra_meshes: Vec::new(),
            p_ndeg_parts: None,
            p_element_nodes: None,
            adaptive: Default::default(),
            meshes: Default::default(),
        }
    }

    /// `Ndeg` of a nodal solver (AddEquationBasics, MainUtils.F90:1657-1712):
    /// the nodes, the mesh's bubbles when the solver keeps bubbles in the
    /// global system (`SetGlobalBubblesFlag`), and the edge and face dofs.
    /// The p layouts of the other solvers are built once to find them.
    fn nodal_ndeg(&mut self, mesh: &dyn MeshView, slot: &SolverSlot) -> Result<usize, SolverError> {
        if self.p_ndeg_parts.is_none() {
            let mut parts = (0, 0);
            if let Some(b) = &self.backends.dofs {
                let nb = mesh.n_bulk();
                for (j, params) in self.model.solvers.iter().enumerate() {
                    if j + 1 == slot.id || params.logical("discontinuous galerkin") == Some(true) {
                        continue;
                    }
                    let active: Vec<ElemId> = (0..(nb + mesh.n_boundary()) as u32)
                        .map(ElemId)
                        .filter(|&e| {
                            let tag = mesh.element(e).tag;
                            let body = if (e.0 as usize) < nb { self.model.body_of_tag(tag) } else { self.model.boundary_body(tag) };
                            body.is_some_and(|bd| self.model.solver_active_on_body(j + 1, bd))
                        })
                        .collect();
                    let dofs = params.integer("variable dofs").unwrap_or(1).max(1) as usize;
                    if let Some(d) = b.for_solver_id(mesh, &self.model, j + 1, params, &active, dofs) {
                        let d = d?;
                        let (bub, ef) = d.ndeg_parts();
                        parts = (parts.0.max(bub), parts.1.max(ef));
                        if self.p_element_nodes.is_none() {
                            let all: Vec<Vec<u32>> = (0..(nb + mesh.n_boundary()) as u32)
                                .map(|e| {
                                    let own = mesh.element(ElemId(e)).nodes;
                                    let p = d.element_nodes(ElemId(e));
                                    if p.len() == own.len() { p } else { own.to_vec() }
                                })
                                .collect();
                            let differs = all.iter().enumerate().any(|(e, n)| n[..] != *mesh.element(ElemId(e as u32)).nodes);
                            if differs {
                                self.p_element_nodes = Some(Arc::new(all));
                            }
                        }
                    }
                }
            }
            self.p_ndeg_parts = Some(parts);
        }
        let (bub, ef) = self.p_ndeg_parts.unwrap_or((0, 0));
        let p = &slot.params;
        let gb = p
            .logical("bubbles in global system")
            .unwrap_or_else(|| p.string("element").is_some_and(|s| s.contains("p:") || s.contains("b:")));
        Ok(mesh.n_nodes() + if gb { bub } else { 0 } + ef)
    }

    /// Installs the `Output File` writer for this and later runs.
    pub(crate) fn clear_node_owned(&mut self) {
        self.node_owned = None;
    }

    pub fn set_save_current(&mut self, writer: Box<dyn SaveCurrent>) {
        self.save_current = Some(writer);
    }

    /// Node ownership in a partitioned run: `owned[i]` when this rank owns
    /// mesh node `i` (`ParallelInfo % NeighbourList(i) % Neighbours(1) ==
    /// MyPE`). `Nonlinear System Consistent Norm` then counts each shared
    /// dof once, on its owner (`ComputeNorm`, SolverBasics.F90:3216).
    /// The non-nodal dofs of slot `ix` (`None` for a nodal solver).
    pub fn slot_dofs(&self, ix: usize) -> Option<Arc<dyn crate::SolverDofs>> {
        self.slots[ix].dofs.clone()
    }

    pub fn set_node_owned(&mut self, owned: Arc<[bool]>) {
        self.node_owned = Some(owned);
    }

    /// Marks the elements (bulk then boundary, mesh order) that belong to
    /// another partition, before [`setup`](Simulation::setup). They are
    /// left out of every solver's active elements (`SetActiveElementsTable`
    /// skips `PartIndex /= myPE`), but a DG solver numbers their DG nodes
    /// (`InitialPermutation`'s DG branch tests only the equation), and a face
    /// between partitions is assembled at half weight on both
    /// ([`SolverCtx::partition_share`](crate::SolverCtx::partition_share)).
    pub fn set_halo_elements(&mut self, halo: Arc<[bool]>) {
        self.halo = Some(halo);
    }

    /// Marks the boundary elements (by element index, bulk then boundary)
    /// that `RadiationParallelMeshDistribute` imported from other
    /// partitions: `ActiveBoundaryElement` is false for them
    /// (`PartIndex /= myPE`), so no solver assembles them, but the
    /// radiation factors see them
    /// ([`SolverCtx::foreign_boundary_elements`](crate::SolverCtx::foreign_boundary_elements)).
    pub fn set_foreign_boundary(&mut self, foreign: Arc<[bool]>) {
        self.foreign = Some(foreign);
    }

    /// [`Simulation::set_foreign_boundary`]'s flags.
    pub fn foreign_boundary(&self) -> Option<&[bool]> {
        self.foreign.as_deref()
    }

    /// The global element ids (`GElementIndex`, bulk then boundary) of a
    /// partitioned mesh.
    pub fn set_element_gids(&mut self, gids: Arc<[u64]>) {
        self.element_gids = Some(gids);
    }

    /// [`Simulation::set_element_gids`]'s ids.
    pub fn element_gids(&self) -> Option<&[u64]> {
        self.element_gids.as_deref()
    }

    /// The global node ids (`GlobalDOFs`, 1-based) of a partitioned mesh.
    pub fn set_node_gids(&mut self, gids: Arc<[u64]>) {
        self.node_gids = Some(gids);
    }

    /// [`Simulation::set_node_gids`]'s ids.
    pub fn node_gids(&self) -> Option<&[u64]> {
        self.node_gids.as_deref()
    }

    /// [`Simulation::set_halo_elements`]'s flags, `None` without halo.
    pub fn halo_elements(&self) -> Option<&[bool]> {
        self.halo.as_deref()
    }

    /// The solver's active elements and, in a run with halo elements, the
    /// halo bulk elements whose body has the solver: the elements
    /// `InitialPermutation` numbers (it tests the equation, not
    /// `PartIndex`), in mesh order.
    pub(crate) fn active_with_halo(&self, slot: &SolverSlot, mesh: &dyn MeshView) -> Vec<ElemId> {
        let Some(h) = self.halo.as_ref() else { return slot.active.clone() };
        let mut v = slot.active.clone();
        v.extend((0..mesh.n_bulk() as u32).map(ElemId).filter(|e| {
            h[e.0 as usize]
                && self
                    .model
                    .body_of_tag(mesh.element(*e).tag)
                    .is_some_and(|b| self.model.solver_active_on_body(slot.id, b))
        }));
        v.sort_unstable_by_key(|e| e.0);
        v
    }

    /// [`Simulation::set_node_owned`]'s flags, `None` in a serial run.
    pub(crate) fn node_owned(&self) -> Option<&[bool]> {
        self.node_owned.as_deref()
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

    /// The coarser meshes of `MG Equal Split` geometric multigrid for solver
    /// `solver_id` (`Linear System Solver = Multigrid`), coarsest first: the
    /// meshes `SplitMeshEqual` refined into the simulation's mesh, which is
    /// the finest level (MainUtils.F90:2730-2770).
    pub fn set_multigrid_meshes(&mut self, solver_id: usize, coarse: Vec<Box<dyn MeshView>>) -> Result<(), SolverError> {
        let ix = self
            .slot_ix(solver_id)
            .ok_or_else(|| SolverError::fatal("MultiGrid", format!("no solver {solver_id}")))?;
        self.slots[ix].mg_meshes = Some(coarse);
        Ok(())
    }

    /// The fields as `elmer-expr` reads them.
    pub fn fields<'a>(&'a self, mesh: &'a dyn MeshView) -> SimFields<'a> {
        SimFields { vars: &self.vars, mesh, dg: self.dg_indexes.as_ref() }
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
        // InitializeRandomSeed (ElmerSolver.F90:874): the RANDOM_NUMBER
        // stream IDR(s) draws its shadow spaces from.
        let seed = self.model.simulation.integer("random number seed").unwrap_or(314_159_265);
        elmer_itermethods::initialize_random_seed(seed);
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
        // AddSolvers: each solver on its mesh (`Solver % Mesh`), with that
        // mesh's variable list.
        for i in 0..self.model.solvers.len() {
            let m = self.mesh_index_of(&self.model.solvers[i].clone())?;
            let mr = self.mesh_ref(m, mesh);
            self.enter_mesh(m);
            let built = self.build_slot(i, mr.get(), Some(registry), true);
            self.enter_mesh(0);
            if let Some(mut slot) = built? {
                slot.mesh_ix = m;
                self.slots.push(slot);
            }
        }
        for ix in 0..self.slots.len() {
            let m = self.slots[ix].mesh_ix;
            let mr = self.mesh_ref(m, mesh);
            self.enter_mesh(m);
            let r = self.init_slot(ix, mr.get(), log, comm);
            self.enter_mesh(0);
            r?;
        }
        // ElmerSolver.F90:578-585, after the initial conditions.
        for ix in 0..self.slots.len() {
            if self.slots[ix].params.logical("initialize exported variables") == Some(true) {
                let m = self.slots[ix].mesh_ix;
                let mr = self.mesh_ref(m, mesh);
                self.enter_mesh(m);
                let r = self.update_exported_variables(ix, mr.get(), log);
                self.enter_mesh(0);
                r?;
            }
        }
        if self.meshes.is_empty() {
            return self.initialize_dirichlet(mesh);
        }
        // SetInitialConditions on every mesh, in `Model % Meshes` order; the
        // Dirichlet values read the other meshes' fields through copies.
        for m in 0..=self.solver_meshes().len() {
            let mr = self.mesh_ref(m, mesh);
            self.enter_mesh(m);
            let r = self
                .make_mesh_copies(m, mr.get())
                .and_then(|_| self.refresh_mesh_copies(mr.get(), mesh, log))
                .and_then(|_| self.initialize_dirichlet(mr.get()));
            self.enter_mesh(0);
            r?;
        }
        Ok(())
    }

    /// The slot of `Solver i+1` (`AddEquationBasics`): procedure from the
    /// registry or the solver mode, active elements, time order and primary
    /// variable, with initial conditions when `apply_ics`. `None` for an
    /// empty section.
    pub(crate) fn build_slot(
        &mut self,
        i: usize,
        mesh: &dyn MeshView,
        registry: Option<&SolverRegistry>,
        apply_ics: bool,
    ) -> Result<Option<SolverSlot>, SolverError> {
        let params = self.model.solvers[i].clone();
        if params.0.is_empty() {
            return Ok(None);
        }
        let mut slot = SolverSlot::new(i + 1, params)?;
        // MainUtils.F90:2575: `Eigen System Values` eigenpairs for an `Eigen
        // Analysis` solver (`Solver % NOFEigenValues`). The harmonic path
        // keeps its own bookkeeping.
        if slot.params.logical("eigen analysis") == Some(true) {
            slot.nof_eigen = slot.params.integer("eigen system values").filter(|&n| n > 0).map_or(0, |n| n as usize);
        }
        // MainUtils.F90:2631: a `Harmonic Analysis` solver keeps `Harmonic
        // System Values` (1 by default) solutions as its eigenvectors, and
        // gets a mass matrix.
        else if slot.params.logical("harmonic analysis") == Some(true) {
            slot.nof_eigen = slot.params.integer("harmonic system values").filter(|&n| n > 1).map_or(1, |n| n as usize);
        }
        // Without a registry (a mesh change) the caller moves the instance.
        if let Some(registry) = registry {
            if let Some((lib, entry)) = procedure_of(&slot.params) {
                if let Some(f) = registry.get(&lib, &entry) {
                    slot.instance = Some(f(&slot.params)?);
                }
            }
            if let Some(s) = crate::coupled::solver_for_mode(&slot.params, registry)? {
                slot.instance = Some(s);
            }
        }
        // SetActiveElementsTable: bulk and boundary elements whose body runs
        // the solver; a boundary element has a body only through its BC's
        // `Body Id`.
        let nb = mesh.n_bulk();
        slot.active = (0..(nb + mesh.n_boundary()) as u32)
            .map(ElemId)
            .filter(|&e| {
                let tag = mesh.element(e).tag;
                let body = if (e.0 as usize) < nb { self.model.body_of_tag(tag) } else { self.model.boundary_body(tag) };
                body.is_some_and(|b| self.model.solver_active_on_body(i + 1, b))
            })
            .filter(|e| self.halo.as_ref().is_none_or(|h| !h[e.0 as usize]))
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
                if !self.is_transient() {
                    self.add_derivative_variable(&slot.params, id)?;
                }
            }
        }
        self.add_exported_variables(&slot, mesh, apply_ics)?;
        Ok(Some(slot))
    }

    /// `<var> Derivative` for `Calculate Derivative` (or `Nonlinear
    /// Calculate Derivative`) outside a transient run (AddEquationBasics,
    /// MainUtils.F90:2566-2577): the variable's dofs and permutation, zero.
    fn add_derivative_variable(&mut self, params: &ParamMap, primary: VarId) -> Result<(), SolverError> {
        let on = params
            .logical("calculate derivative")
            .or_else(|| params.logical("nonlinear calculate derivative"));
        if on != Some(true) {
            return Ok(());
        }
        let p = self.vars.by_id(primary);
        let name = format!("{} Derivative", p.name);
        if self.vars.id(&name).is_some() {
            return Ok(());
        }
        let v = Variable::new(&name, p.kind, p.dofs, p.perm.clone(), vec![0.0; p.values.len()]);
        self.vars.add(v)?;
        Ok(())
    }

    /// `Exported Variable n` of a solver section (AddEquationBasics,
    /// MainUtils.F90:1876-2160): nodal fields on the solver's permutation
    /// (all nodes of its active elements, or every node when it has none),
    /// or `-global` ones; `Exported Variable n Output` and `-nooutput`
    /// control output. Initial conditions are applied when `apply_ics`.
    /// `-elem` fields live on the active elements, `-ip` ones on their
    /// Gauss points (`exported`), both optionally masked and initialised.
    /// `-dg` fields live on the DG indexes of the active elements
    /// (`CreateDGPerm`), optionally masked. Masks of nodal fields, initial
    /// conditions of DG fields and time derivatives are refused.
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
            if mask.is_some() && !matches!(spec.kind, FieldKind::Elemental | FieldKind::Ip | FieldKind::Dg) {
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
                        // `InitialPermutation` tests the equation only, so
                        // the nodes of halo elements are numbered too.
                        None => node_order_perm(mesh, &self.active_with_halo(slot, mesh)).perm.into(),
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
                FieldKind::Dg => {
                    if apply_ics && self.model.initial_conditions.iter().any(|ic| {
                        ic.get(&spec.name).is_some() || spec.components.iter().any(|c| ic.get(c).is_some())
                    }) {
                        return Err(SolverError::Unsupported(format!(
                            "initial condition of the DG exported variable {}",
                            spec.name
                        )));
                    }
                    let (perm, n) = crate::exported::dg_perm(self, mesh, slot, mask)?;
                    let mut v = Variable::new(&spec.name, FieldKind::Dg, spec.dofs, Some(perm.into()), vec![0.0; n * spec.dofs]);
                    v.output = spec.output;
                    self.vars.add_with_components(v, &spec.components)?;
                }
                k => return Err(SolverError::Unsupported(format!("{k:?} exported variable {}", spec.name))),
            }
            self.add_exported_velocity(slot, &spec.name)?;
        }
        Ok(())
    }

    /// `<name> Transient` (else `<name> Calculate Velocity`, `<name>
    /// Calculate Acceleration`) of an exported variable (MainUtils.F90:2124):
    /// two previous values, and with `<name> Calculate Velocity`
    /// (CreateTimeDerivativeVariables) a `<name> Velocity` of the same
    /// kind, size and permutation. The previous values only change at a new
    /// time step, so a steady run keeps them zero; transient runs and the
    /// acceleration are refused.
    fn add_exported_velocity(&mut self, slot: &SolverSlot, name: &str) -> Result<(), SolverError> {
        let key = |k: &str| format!("{name} {k}");
        let p = &slot.params;
        let transient_var = match p.logical(&key("transient")) {
            Some(t) => t,
            None => match p.logical(&key("calculate velocity")) {
                Some(t) => t,
                None => p.logical(&key("calculate acceleration")).unwrap_or(false),
            },
        };
        if !transient_var {
            return Ok(());
        }
        if p.logical(&key("calculate acceleration")) == Some(true) {
            return Err(SolverError::Unsupported(format!("{name} Calculate Acceleration")));
        }
        if self.model.simulation.string("simulation type").is_some_and(|t| normalize_key(t) == "transient") {
            return Err(SolverError::Unsupported(format!("time derivatives of the exported variable {name} in a transient run")));
        }
        if p.logical(&key("calculate velocity")) != Some(true) {
            return Ok(());
        }
        let Some(id) = self.vars.id(name) else { return Ok(()) };
        let v = self.vars.by_id(id);
        let vel = format!("{name} velocity");
        if self.vars.id(&vel).is_some() {
            return Ok(());
        }
        let comps: Vec<String> = if v.dofs > 1 { (1..=v.dofs).map(|k| format!("{vel} {k}")).collect() } else { Vec::new() };
        let nv = Variable::new(&vel, v.kind, v.dofs, v.perm.clone(), vec![0.0; v.values.len()]);
        self.vars.add_with_components(nv, &comps)?;
        Ok(())
    }

    /// `DerivateExportedVariables` (SolverBasics.F90:6909) at the end of a
    /// steady visit: `<name> Velocity = (Values - PrevValues(:,1)) / dt`
    /// for the exported variables that keep previous values (zero in a
    /// steady run) and ask for their velocity.
    fn derivate_exported_variables(&mut self, ix: usize) {
        let dt = self.vars.scalar("timestep size").unwrap_or(1.0);
        let p = self.slots[ix].params.clone();
        for l in 1.. {
            let Some(text) = p.string(&format!("exported variable {l}")) else { break };
            let Ok(spec) = parse_variable_spec(text, p.integer(&format!("exported variable {l} dofs"))) else { continue };
            let name = spec.name;
            let transient_var = match p.logical(&format!("{name} transient")) {
                Some(t) => t,
                None => p.logical(&format!("{name} calculate velocity")).unwrap_or(false),
            };
            if !transient_var || p.logical(&format!("{name} calculate velocity")) != Some(true) {
                continue;
            }
            let (Some(src), Some(dst)) = (self.vars.id(&name), self.vars.id(&format!("{name} velocity"))) else { continue };
            let vals: Vec<f64> = self.vars.by_id(src).values.iter().map(|&x| (x - 0.0) / dt).collect();
            let d = self.vars.by_id_mut(dst);
            if d.values.len() == vals.len() {
                d.values.copy_from_slice(&vals);
            }
        }
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
                    if let Some(slot) = self.build_slot(i, mesh, Some(registry), false)? {
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
    /// list order, so a later element wins at a shared node. A component of a
    /// vector on a normal-tangential boundary sets the vector's normal
    /// (tangential) projection; global variables are skipped.
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
                let fields = SimFields { vars: &self.vars, mesh, dg: self.dg_indexes.as_ref() };
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
                // ElmerSolver.F90:1899-1993: component `vect_dof` of a vector
                // field on a normal-tangential boundary sets the projection of
                // the nodal vector on the element's normal (tangent) instead.
                let nt_dir = match (el.etype.family() != 1, comp, name.rsplit_once(' ')) {
                    (true, Some(c), Some((base, digit))) => {
                        let vect_dof: usize = digit.parse().unwrap_or(0);
                        let dim = mesh.dim();
                        if (1..=3).contains(&vect_dof)
                            && bc.logical(&format!("normal-tangential {base}")) == Some(true)
                            && dofs >= dim
                            && c + 1 >= vect_dof
                        {
                            let nrm = crate::nt::boundary_normal(mesh, &self.model, e, [0.0, 0.0])?;
                            let vec = match (vect_dof, el.etype.family()) {
                                (1, _) => nrm,
                                (2, 2) => [nrm[1], -nrm[0], 0.0],
                                (k, _) => {
                                    let (t1, t2) = crate::nt::tangent_directions(nrm);
                                    if k == 2 { t1 } else { t2 }
                                }
                            };
                            Some((vec, c + 1 - vect_dof, dim))
                        } else {
                            None
                        }
                    }
                    _ => None,
                };
                let xs = crate::expr::eval_at_nodes_with(val, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                if let Some((vec, first, dim)) = nt_dir {
                    let v = self.vars.by_id_mut(id);
                    for (j, &n) in el.nodes.iter().enumerate() {
                        let Some(s) = v.slot(n) else { continue };
                        let mut tmp = [0.0; 3];
                        for (l, t) in tmp.iter_mut().enumerate().take(dim) {
                            *t = v.values[dofs * s + l + first];
                        }
                        let udot = (0..dim).fold(0.0, |acc, l| acc + vec[l] * tmp[l]);
                        for l in 0..dim {
                            tmp[l] += (xs[j] - udot) * vec[l];
                        }
                        for (l, t) in tmp.iter().enumerate().take(dim) {
                            v.values[dofs * s + l + first] = *t;
                        }
                    }
                    continue;
                }
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
                // ElmerSolver.F90:1958-1979: a DG variable takes the node's
                // DG index in the left parent, else in the right one.
                let dg_slots: Option<Vec<Option<usize>>> = if self.vars.by_id(id).kind == FieldKind::Dg {
                    let ix = self.dg_indexes.as_ref().ok_or_else(|| {
                        SolverError::fatal("SetInitialConditions", "a DG variable without the mesh's DG indexes")
                    })?;
                    let v = self.vars.by_id(id);
                    let through = |p: Option<ElemId>, n: u32| -> Option<usize> {
                        let p = p?;
                        let pos = mesh.element(p).nodes.iter().position(|&m| m == n)?;
                        let d = *ix.index.get(p.0 as usize)?.get(pos)?;
                        (d != NO_DOF).then(|| v.slot(d)).flatten()
                    };
                    Some(el.nodes.iter().map(|&n| through(el.parents[0], n).or_else(|| through(el.parents[1], n))).collect())
                } else {
                    None
                };
                let v = self.vars.by_id_mut(id);
                let (stride, offset) = match comp {
                    Some(c) => (dofs, c),
                    None => (1, 0),
                };
                for (j, &n) in el.nodes.iter().enumerate() {
                    let s = match &dg_slots {
                        Some(d) => d[j],
                        None => v.slot(n),
                    };
                    let Some(s) = s else { continue };
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
        // MainUtils.F90:1437-1456: `Predictor Method` / `Corrector Method`
        // from the Simulation when the solver lacks them; with a corrector
        // the predictor-corrector defaults are read and the corrector becomes
        // the `Timestepping Method`.
        for key in ["predictor method", "corrector method"] {
            if slot.params.get(key).is_none() {
                if let Some(s) = self.model.simulation.string(key) {
                    let s = s.to_string();
                    slot.params.insert(key, Value::Str(s));
                }
            }
        }
        if let Some(c) = slot.params.string("corrector method").map(str::to_string) {
            read_pred_corr_params(&self.model.simulation, &mut slot.params);
            slot.params.insert("Timestepping Method", Value::Str(c));
        }
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
            "runge-kutta" => return Err(SolverError::Unsupported(format!("Timestepping Method = {method}"))),
            m if m.starts_with("adams") => 2,
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
            if slot.time_order > 0 {
                // The previous solutions, as for a field (MainUtils.F90:2510-2517).
                let m = if slot.time_order == 2 { 7 } else { slot.order.max(slot.time_order as usize) };
                v.prev_values = vec![vec![0.0; v.values.len()]; m];
            }
            slot.variable = Some(self.vars.add_with_components(v, &spec.components)?);
            return Ok(());
        }
        if spec.kind != FieldKind::Nodal {
            return Err(SolverError::Unsupported(format!("{:?} solver variable", spec.kind)));
        }
        // A DG solver also numbers the DG nodes of its halo elements.
        let dg_halo = slot.params.logical("discontinuous galerkin") == Some(true) && self.halo.is_some();
        let with_halo = if dg_halo { self.active_with_halo(slot, mesh) } else { Vec::new() };
        let layout_active: &[ElemId] = if dg_halo { &with_halo } else { &slot.active };
        let layout: Option<Arc<dyn crate::SolverDofs>> = match &self.backends.dofs {
            Some(b) => b
                .for_solver_id(mesh, &self.model, slot.id, &slot.params, layout_active, spec.dofs)
                .transpose()?
                .map(Arc::from),
            None => None,
        };
        let conforming = slot.params.logical("apply conforming bcs") == Some(true);
        let perm = match (&layout, &self.backends.assembly) {
            (Some(_), _) if conforming => {
                return Err(SolverError::Unsupported(format!(
                    "Solver {}: Apply Conforming BCs with p element dofs",
                    slot.id
                )))
            }
            (Some(d), _) => d.perm(),
            // CreateMatrix: no bandwidth optimisation for conforming BCs.
            (None, Some(b)) if conforming => b.build_perm(mesh, &slot.active, false),
            (None, Some(_)) => {
                let ndeg = self.nodal_ndeg(mesh, slot)?;
                let b = self.backends.assembly.as_ref().expect("matched above");
                if ndeg > mesh.n_nodes() {
                    match &self.p_element_nodes {
                        Some(nodes) => {
                            let view = PNodeOrder { mesh, nodes };
                            b.build_perm_ndeg(&view, &slot.active, &slot.params, ndeg)?
                        }
                        None => b.build_perm_ndeg(mesh, &slot.active, &slot.params, ndeg)?,
                    }
                } else {
                    b.build_perm_for(mesh, &slot.active, &slot.params)?
                }
            }
            (None, None) => node_order_perm(mesh, &slot.active),
        };
        let mut perm = perm;
        if conforming {
            // MainUtils: conforming BCs enforce `Use Global Mass Matrix`.
            slot.params.insert("Use Global Mass Matrix", Value::Logical(true));
            // Without a map (no `Conforming BC`) Elmer warns that the
            // conforming BC is requested but not generated, and numbers the
            // nodes as usual. A model with `Conforming BC` sections whose map
            // the front end did not build (the partitioned runner) would run
            // without them.
            match &self.backends.conforming {
                Some(map) => {
                    slot.conforming_flip = crate::conforming::apply_conforming(&mut perm, map);
                }
                None if self.model.boundary_conditions.iter().any(|b| b.get("conforming bc").is_some()) => {
                    return Err(SolverError::Unsupported(format!(
                        "Solver {}: Apply Conforming BCs without the conforming map (partitioned runs)",
                        slot.id
                    )));
                }
                None => {}
            }
        }
        // A DG solver: its variable lives on the mesh's DG indexes
        // (`Variable_on_nodes_on_elements`), made at mesh load.
        let dg = layout.as_ref().and_then(|d| d.dg_layout().cloned());
        let kind = if dg.is_some() { FieldKind::Dg } else { FieldKind::Nodal };
        if let Some(ix) = &dg {
            if self.dg_indexes.is_none() {
                self.dg_indexes = Some(ix.clone());
            }
        }
        // AddEquationBasics (MainUtils.F90:1764-1822): CreateMatrix returns
        // no matrix for a `No Matrix` solver, so its variable gets `DOFs *
        // Ndeg` values (every node's) whatever its permutation covers; the
        // norm runs over all of them.
        let n_values = if slot.dofs.is_none() && layout.is_none() && kind == FieldKind::Nodal && slot.params.logical("no matrix") == Some(true) {
            self.nodal_ndeg(mesh, slot)?.max(perm.n_slots)
        } else {
            perm.n_slots
        };
        slot.dofs = layout;
        let perm_arc: Arc<[u32]> = perm.perm.into();
        let mut values = vec![0.0; n_values * spec.dofs];
        if apply_ics {
            if kind == FieldKind::Dg {
                let ix = dg.as_ref().expect("a DG kind has a layout");
                self.apply_dg_initial_conditions(slot, spec, mesh, ix, &perm_arc, &mut values)?;
            } else {
                self.apply_initial_conditions(slot, spec, mesh, &perm_arc, &mut values)?;
            }
        }
        let mut v = Variable::new(&spec.name, kind, spec.dofs, Some(perm_arc.clone()), values);
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
            // A boundary element of a BC with `Body Id` is that body's.
            let body = if (e.0 as usize) < mesh.n_bulk() { self.model.body_of_tag(el.tag) } else { self.model.boundary_body(el.tag) };
            let Some(ic) = body.and_then(|b| self.model.body_params(b)).and_then(|b| b.initial_condition) else {
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

    /// `InitCond` (ElmerSolver.F90) for a `Variable_on_nodes_on_elements`:
    /// each bulk element with a body `Initial Condition` sets the values at
    /// its DG indexes, `Var % Values(Var % Perm(DGIndexes(k))) = Work(k)`,
    /// in element order (a node shared by elements of one reduced-basis
    /// group takes the last one). Initial conditions on boundary elements
    /// (a BC `Initial Condition`) are refused.
    fn apply_dg_initial_conditions(
        &self,
        slot: &SolverSlot,
        spec: &VarSpec,
        mesh: &dyn MeshView,
        ix: &crate::exported::DgIndexes,
        perm: &[u32],
        values: &mut [f64],
    ) -> Result<(), SolverError> {
        if spec.dofs != 1 {
            return Err(SolverError::Unsupported(format!(
                "Solver {}: initial conditions of the {}-component DG variable {}",
                slot.id, spec.dofs, spec.name
            )));
        }
        let key = spec.name.as_str();
        if !self.model.initial_conditions.iter().any(|ic| ic.get(key).is_some()) {
            return Ok(());
        }
        if self.model.boundary_conditions.iter().any(|bc| bc.get("initial condition").is_some()) {
            return Err(SolverError::Unsupported(format!(
                "Solver {}: boundary initial conditions of the DG variable {}",
                slot.id, spec.name
            )));
        }
        let fields = self.fields(mesh);
        let b = &self.backends;
        for (t, dgi) in ix.index.iter().enumerate().take(mesh.n_bulk()) {
            let el = mesh.element(ElemId(t as u32));
            let Some(ic) = self
                .model
                .body_of_tag(el.tag)
                .and_then(|b| self.model.body_params(b))
                .and_then(|b| b.initial_condition)
            else {
                continue;
            };
            let Some(val) = ic.get(key) else { continue };
            let xs = crate::expr::eval_at_nodes_with(val, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
            for (&d, &x) in dgi.iter().zip(&xs) {
                if d == NO_DOF {
                    continue;
                }
                let s = perm.get(d as usize).copied().unwrap_or(NO_DOF);
                if s != NO_DOF {
                    values[s as usize] = x;
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
        if p.logical("nonlinear calculate velocity") == Some(true) {
            return Err(SolverError::Unsupported("ComputeChange option nonlinear calculate velocity".into()));
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
        let fields = SimFields { vars: &self.vars, mesh, dg: self.dg_indexes.as_ref() };
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
            node_owned: self.node_owned.as_deref(),
            row_owned: slot.owned_rows.as_deref(),
            log,
            comm,
        };
        let var = self.vars.by_id_mut(id);
        let report = compute_change(&env, var, &mut slot.state)?;
        if let Some(dx) = &report.derivative {
            let eps = params.real("derivative eps");
            self.store_derivative(id, steady, eps, dx, log);
        }
        // UpdateDependentObjects(Solver, .TRUE.): the exported velocities.
        if steady {
            self.derivate_exported_variables(ix);
        }
        // UpdateDependentObjects(Solver, .FALSE.) (SolverBasics.F90:4242).
        if !steady && self.slots[ix].params.logical("nonlinear update exported variables") == Some(true) {
            self.update_exported_variables(ix, mesh, log)?;
        }
        if !steady && self.slots[ix].params.get("nonlinear update components").is_some() {
            return Err(SolverError::Unsupported("Nonlinear Update Components".into()));
        }
        Ok(Some(report))
    }

    /// The sensitivity `<var> Derivative = (x - x0) / eps` (SolverBasics.F90:
    /// 4198-4238): `eps` is 1 on the first coupled (steady) or nonlinear
    /// iteration, then the global `derivative eps` or the section's
    /// `Derivative Eps`, else 1.
    fn store_derivative(&mut self, id: VarId, steady: bool, eps_keyword: Option<f64>, dx: &[f64], log: &dyn Logger) {
        const CALLER: &str = "ComputeChange";
        let iter = if steady { "coupled iter" } else { "nonlin iter" };
        let iter_no = self.vars.scalar(iter).map_or(0, |v| v.round() as i64);
        let mut eps = 1.0;
        if iter_no > 1 {
            match self.vars.scalar("derivative eps").or(eps_keyword) {
                Some(e) => eps = e,
                None => log.info(CALLER, "Derivative Eps not given, using one", 7),
            }
        }
        let name = format!("{} Derivative", self.vars.by_id(id).name);
        match self.vars.get_mut(&name) {
            Some(v) => {
                log.info(CALLER, &format!("Computing variable:{name}"), 7);
                for (d, x) in v.values.iter_mut().zip(dx) {
                    *d = x / eps;
                }
                v.values_changed = true;
            }
            None => log.warn(CALLER, "Derivative variable not present"),
        }
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
        if self.meshes.is_empty() {
            return self.solver_activate_on(ix, dt, transient, mesh, log, comm);
        }
        // A solver of its own mesh (`mesh` is the run's first mesh): that
        // mesh's variables, with the other meshes' fields it reads
        // interpolated again where they changed.
        let m = self.slots[ix].mesh_ix;
        let prev = self.meshes.active;
        let mr = self.mesh_ref(m, mesh);
        self.enter_mesh(m);
        let r = self
            .refresh_mesh_copies(mr.get(), mesh, log)
            .and_then(|_| self.solver_activate_on(ix, dt, transient, mr.get(), log, comm));
        if r.is_ok() {
            self.invalidate_mesh_copies(ix);
        }
        self.enter_mesh(prev);
        r
    }

    #[allow(clippy::too_many_arguments)]
    fn solver_activate_on(
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
        // MainUtils.F90:5816: `Namespace` is pushed for the activation
        // (`ListPushNamespace`), so `<namespace> <key>` stands for `<key>` in
        // every list until it is popped.
        let namespace = if auxiliary { None } else { p.string("namespace").map(str::to_string) };
        let overlay = namespace.as_deref().map(|ns| self.push_namespace(ns, log));
        let r = self.activate_rounds(ix, dt, dts, time_derivative_active, auxiliary, mesh, log, comm);
        if let Some(o) = overlay {
            self.pop_namespace(o);
        }
        r
    }

    /// `ListPushNamespace(ns)` over the model's sections and the solvers'
    /// lists: each `<ns> <key>` entry (lower-cased, as `ListFind` builds it)
    /// takes the place of `<key>`. Returns what it replaced.
    fn push_namespace(&mut self, ns: &str, log: &dyn Logger) -> Vec<(usize, String, Option<Value>)> {
        let ns = ns.trim().to_lowercase();
        log.info("ListPushNameSpace", &format!("Adding name space: {ns}"), 12);
        let mut saved = Vec::new();
        for (i, sec) in self.namespace_lists().into_iter().enumerate() {
            let scoped: Vec<(String, Value)> = sec
                .0
                .iter()
                .filter_map(|(k, v)| Some((k.strip_prefix(&ns)?.trim_start().to_string(), v.clone())))
                .filter(|(k, _)| !k.is_empty())
                .collect();
            for (k, v) in scoped {
                let old = sec.0.insert(k.clone(), v);
                saved.push((i, k, old));
            }
        }
        saved
    }

    /// `ListPopNamespace`: puts back what [`push_namespace`](Self::push_namespace) replaced.
    fn pop_namespace(&mut self, saved: Vec<(usize, String, Option<Value>)>) {
        let mut lists = self.namespace_lists();
        for (i, k, old) in saved.into_iter().rev() {
            match old {
                Some(v) => lists[i].0.insert(k, v),
                None => lists[i].0.remove(&k),
            };
        }
    }

    /// Every keyword list a namespace applies to, in a fixed order.
    fn namespace_lists(&mut self) -> Vec<&mut ParamMap> {
        let m = &mut self.model;
        let mut v: Vec<&mut ParamMap> = vec![&mut m.header, &mut m.simulation, &mut m.constants];
        v.extend(m.solvers.iter_mut());
        v.extend(m.equations.iter_mut());
        v.extend(m.bodies.iter_mut());
        v.extend(m.materials.iter_mut());
        v.extend(m.body_forces.iter_mut());
        v.extend(m.initial_conditions.iter_mut());
        v.extend(m.boundary_conditions.iter_mut());
        v.extend(self.slots.iter_mut().map(|s| &mut s.params));
        v
    }

    /// The solver calls of one activation: one per `Scanning Loops` round.
    #[allow(clippy::too_many_arguments)]
    fn activate_rounds(
        &mut self,
        ix: usize,
        dt: f64,
        dts: f64,
        time_derivative_active: bool,
        auxiliary: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        let p = &self.slots[ix].params;
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
        // MainUtils.F90:5287: the integration weights of the solver's nodes
        // when the active elements are set up (once per mesh).
        if !self.slots[ix].weights_done {
            self.slots[ix].weights_done = true;
            let p = &self.slots[ix].params;
            let (bulk, bc) = (p.logical("calculate weights") == Some(true), p.logical("calculate boundary weights") == Some(true));
            if bulk || bc {
                SolverCtx::new(self, ix, mesh, log, comm).library_nodal_weights(bulk, bc)?;
            }
        }
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
            // ... and the components that depend on the solution.
            if let Some(Value::Integer(list)) = self.slots[ix].params.get("update components").cloned() {
                self.update_dependent_components(&list, mesh, log, comm)?;
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
                let pc = s.exec_when == ExecWhen::PredCorr;
                if s.has_procedure() && (s.exec_when == ExecWhen::BeforeTimestep || pc) && !skip_first(s) {
                    if pc {
                        log.info("SolveEquations", "Switching time-stepping method to predictor method", 7);
                        let m = self.slots[ix].params.string("predictor method").map(str::to_string).ok_or_else(|| {
                            SolverError::fatal("SolveEquations", "Predictor-corrector schemes require > Predictor Method <")
                        })?;
                        self.slots[ix].params.insert("Timestepping Method", Value::Str(m));
                        self.slots[ix].params.insert("Predictor Phase", Value::Logical(true));
                    }
                    self.solver_activate(ix, dt, transient, mesh, log, comm)?;
                    comm.barrier();
                    if pc {
                        log.info("SolveEquations", "Switching time-stepping method to corrector method", 7);
                        let m = self.slots[ix].params.string("corrector method").map(str::to_string).ok_or_else(|| {
                            SolverError::fatal("SolveEquations", "Predictor-corrector schemes require > Corrector Method <")
                        })?;
                        self.slots[ix].params.insert("Timestepping Method", Value::Str(m));
                    }
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
                let pc = s.exec_when == ExecWhen::PredCorr;
                if s.has_procedure() && (s.exec_when == ExecWhen::AfterTimestep || pc) && !skip_first(s) {
                    if pc {
                        self.initialize_timestep(ix)?;
                        self.slots[ix].params.insert("Predictor Phase", Value::Logical(false));
                    }
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
                // The slot's own mesh and variable list (solver meshes).
                let mr = self.mesh_ref(self.slots[k].mesh_ix, mesh);
                self.enter_mesh(self.slots[k].mesh_ix);
                let slot_mesh = mr.get();
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
                    if let Some(c) = self.creal(p, "steady state condition", slot_mesh)? {
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
                    self.compute_change(k, true, None, false, slot_mesh, log, comm)?;
                }
                if test_convergence {
                    done[k] = self.slots[k]
                        .variable
                        .is_some_and(|id| self.vars.by_id(id).conv.steady_converged != 0);
                }
                if self.adaptive.finished && !done[k] {
                    log.info("SolveEquations", "Overriding convergence due to Adaptive Meshing Finished!", 4);
                    done[k] = true;
                }
                if done.iter().all(|&d| d) {
                    break;
                }
            }
            self.enter_mesh(0);
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
    pub fn run(&mut self, first_mesh: &dyn MeshView, log: &dyn Logger, comm: &dyn Communicator) -> Result<RunReport, SolverError> {
        // The mesh a solver moved the run to (adaptive mesh refinement).
        let mut current: Option<Arc<dyn MeshView + Send + Sync>> = None;
        let mesh = first_mesh;
        let simp = self.model.simulation.clone();
        let transient = self.is_transient();
        let scanning = self.is_scanning();
        let timed = transient || scanning;
        for k in ["periodic time", "parallel timestepping"] {
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
        let mut last_saved = false;
        'outer: for (interval, &nsteps) in intervals.iter().enumerate() {
            let mut real_timestep = 1usize;
            for timestep in 1..=nsteps {
                let held = current.clone();
                let mesh: &dyn MeshView = match &held {
                    Some(m) => &**m,
                    None => first_mesh,
                };
                cum += 1;
                self.vars.set_scalar("timestep", cum as f64);
                dt = self.step_size(&simp, timed, cum, interval, timestep, &sizes, ratios.as_deref(), dt, mesh, log)?;
                if simp.logical("predictor-corrector control") == Some(true) {
                    dt = self.predictor_corrector_control(dt, timestep as usize, log)?;
                }
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
                last_saved = oi != 0 && ((timestep - 1) % oi == 0 || reached);
                if last_saved {
                    self.save_step(dt, transient, mesh, log, comm)?;
                }
                // RefineMesh's new mesh: the results above were the old
                // mesh's; the run goes on with the new one.
                if let Some(change) = self.adaptive.pending.take() {
                    // The refining slot's parameters as its solve left them
                    // (FlowSolve adds `Nonlinear System Norm Dofs`).
                    let refining = self.adaptive.pending_var.clone().and_then(|name| {
                        let s = self.slots.iter().find(|s| s.variable.is_some_and(|id| self.vars.by_id(id).name.eq_ignore_ascii_case(&name)))?;
                        Some((name, s.params.clone()))
                    });
                    let new_mesh = self.wrap_mesh(change.mesh.clone());
                    match self.adaptive.pending_coords.take() {
                        Some(c) => {
                            let old = crate::remesh::WithCoords { mesh, coords: &c };
                            self.swap_mesh(&old, &*new_mesh, log, comm)?;
                        }
                        None => self.swap_mesh(mesh, &*new_mesh, log, comm)?,
                    }
                    // Elmer's RefineMesh hands the solver its new variable
                    // before the steady-state ComputeChange, whose norm of the
                    // interpolated values the next nonlinear change compares
                    // against.
                    if let Some((name, params)) = refining {
                        if let Some(ix) = self.slots.iter().position(|s| {
                            s.variable.is_some_and(|id| self.vars.by_id(id).name.eq_ignore_ascii_case(&name))
                        }) {
                            let id = self.slots[ix].variable.expect("found by its variable");
                            let n_nodes = new_mesh.n_nodes();
                            let norm = crate::change::steady_norm(&params, self.vars.by_id(id), n_nodes, comm, log);
                            self.vars.by_id_mut(id).conv.norm = norm;
                        }
                    }
                    self.adaptive.depth += 1;
                    self.adaptive.finished = false;
                    self.adaptive.name = change.name.clone();
                    self.mesh_writer(&change.name)?;
                    current = Some(new_mesh);
                }
                let held2 = current.clone();
                let mesh: &dyn MeshView = match &held2 {
                    Some(m) => &**m,
                    None => first_mesh,
                };
                if self.creal(&simp, "exit condition", mesh)?.is_some_and(|c| c > 0.0) {
                    log.info("MAIN", "Found a positive exit condition, exiting...", 3);
                    break 'outer;
                }
                if reached && !timed && timestep >= coupled_min {
                    break 'outer;
                }
            }
        }

        let held = current.clone();
        let mesh: &dyn MeshView = match &held {
            Some(m) => &**m,
            None => first_mesh,
        };
        self.exec_slots(ExecWhen::AfterAll, dt, transient, mesh, log, comm)?;
        // "Always save the last step to output" (ElmerSolver.F90:3557).
        if !last_saved {
            self.save_step(dt, transient, mesh, log, comm)?;
        }
        crate::reference::compare_to_reference(self, mesh, log, comm)?;
        Ok(report)
    }

    /// One save point: the `before saving` solvers, `SaveCurrent`, then
    /// the `after saving` solvers.
    fn save_step(
        &mut self,
        dt: f64,
        transient: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        self.exec_slots(ExecWhen::BeforeSaving, dt, transient, mesh, log, comm)?;
        self.save_meshes(mesh, log)?;
        self.exec_slots(ExecWhen::AfterSaving, dt, transient, mesh, log, comm)
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

    /// `PredictorCorrectorControl(Model, dt, RealTimestep)` (MainUtils.F90:
    /// 5958-6052): the first step is taken with `dt = 0` (`dtOld` kept), the
    /// second with `dtOld`; from the third the local error estimate of the
    /// predictor-corrector solver's change sets the next size
    /// (`TimeStepController`, `TimeStepLimiter`) and `Adams Zeta`.
    fn predictor_corrector_control(&mut self, dt: f64, real_timestep: usize, log: &dyn Logger) -> Result<f64, SolverError> {
        let ix = self
            .slots
            .iter()
            .rposition(|s| s.exec_when == ExecWhen::PredCorr)
            .ok_or_else(|| SolverError::fatal("Predictor-Corrector Control", "Predictor-Corrector Solver is not found!"))?;
        let (mut dt_old, mut eta_old, mut zeta) = self.pred_corr;
        let mut dt = dt;
        match real_timestep {
            1 => {
                dt_old = dt;
                dt = 0.0;
                zeta = 1.0;
            }
            2 => {
                dt = dt_old;
                zeta = 1.0;
            }
            _ => {
                let simulation = self.model.simulation.clone();
                let (mut order, epsilon, beta1, beta2) = read_pred_corr_params(&simulation, &mut self.slots[ix].params);
                let id = self.slots[ix].variable.ok_or_else(|| SolverError::fatal("Predictor-Corrector Control", "no variable"))?;
                let v = self.vars.by_id(id);
                let prev = v.prev_values.first().ok_or_else(|| SolverError::fatal("Predictor-Corrector Control", "no history"))?;
                let time_error = v.values.iter().zip(prev).fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
                if real_timestep == 3 {
                    order = 1;
                }
                // PredCorrErrorEstimate.
                let eta = if dt_old > 0.0 {
                    if order == 2 {
                        time_error * zeta / dt_old / (zeta + 1.0) / 3.0
                    } else {
                        time_error / dt_old / 2.0
                    }
                } else {
                    log.warn("Predictor-Corrector", "Time Step is 0 in Local error estimate!");
                    0.0
                };
                if real_timestep == 3 {
                    eta_old = eta;
                }
                // TimeStepController / TimeStepLimiter.
                if eta != 0.0 && eta_old != 0.0 {
                    let g = (epsilon / eta).powf(beta1) * (epsilon / eta_old).powf(beta2);
                    let x = 1.0 + 2.0 * ((g - 1.0) * 0.5).atan();
                    dt = dt_old * x;
                } else {
                    dt = dt_old;
                }
                zeta = dt / dt_old;
                self.slots[ix].params.insert("Adams Zeta", Value::real(zeta));
                eta_old = eta;
                if self.slots[ix].params.logical("predictor-corrector save error") == Some(true) {
                    return Err(SolverError::Unsupported("Predictor-Corrector Save Error (ErrorPredictorCorrector.dat)".into()));
                }
                log.info("Predictor-Corrector", " ---------------- Predictor-Corrector Control ----------------------", 4);
                log.info(
                    "Predictor-Corrector",
                    &format!(" current dt={dt_old:>25.17} next dt={dt:>25.17}"),
                    4,
                );
                dt_old = dt;
            }
        }
        self.pred_corr = (dt_old, eta_old, zeta);
        Ok(dt)
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

/// `ReadPredCorrParams(Model, SolverParams, Order, Epsilon, Beta1, Beta2)`
/// (MainUtils.F90): each from the solver, else the Simulation, else its
/// default (order 2, tolerance 1e-6, `0.6/(order+1)`, `-0.2/(order+1)`),
/// written back to the solver section.
pub(crate) fn read_pred_corr_params(simulation: &ParamMap, p: &mut ParamMap) -> (i32, f64, f64, f64) {
    let order = match p.integer("predictor-corrector scheme order") {
        Some(o) => o,
        None => {
            let o = simulation.integer("predictor-corrector scheme order").unwrap_or(2);
            p.insert("Predictor-Corrector Scheme Order", Value::Integer(vec![o]));
            o
        }
    };
    let mut real = |key: &str, name: &str, default: f64| -> f64 {
        match p.real(key) {
            Some(x) => x,
            None => {
                let x = simulation.real(key).unwrap_or(default);
                p.insert(name, Value::real(x));
                x
            }
        }
    };
    let eps = real("predictor-corrector control tolerance", "Predictor-Corrector Control Tolerance", 1.0e-6);
    let b1 = real("predictor-corrector control beta 1", "Predictor-Corrector Control Beta 1", 0.6 / (f64::from(order) + 1.0));
    let b2 = real("predictor-corrector control beta 2", "Predictor-Corrector Control Beta 2", -0.2 / (f64::from(order) + 1.0));
    (order, eps, b1, b2)
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

