//! One case on every rank: front end, parallel preflight, `Simulation` with
//! the parallel backends, then `CompareToReferenceSolution`.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::Mutex;

use elmer_cli::driver::{Engine, EngineOutcome, RunRequest};
use elmer_cli::{ElmerLogger, ExitCode};
use elmer_core::{normalize_key, Communicator, LinalgError, Logger, MeshView, ParamMap, Params, Value};
use elmer_dist::Exchange;
use elmer_run::{Frontend, SifFrontend, SifMatc};
use elmer_solver::{finalize_reference, solver::procedure_of, Model, RunReport, Simulation, SolverError};

use crate::backend::{dist_options, linear_registry, method_of, prepare_slot, ParAssembly, ITERMETHODS, PARALLEL_METHODS};
use crate::context::RankContext;
use crate::mesh::ParMesh;

/// What a finished partitioned run produced, on every rank.
#[derive(Clone, Debug, PartialEq)]
pub struct ParOutcome {
    pub report: RunReport,
    /// `(solver id, ComputeNorm)` of every solver with a primary variable;
    /// the same on all ranks.
    pub norms: Vec<(usize, f64)>,
    /// The `TEST.PASSED_<np>` verdict, on rank 0 only.
    pub passed: Option<bool>,
}

/// The serial front end reads the `.sif`; the mesh is read per rank later.
struct NoMesh;

impl MeshView for NoMesh {
    fn dim(&self) -> usize {
        0
    }
    fn n_nodes(&self) -> usize {
        0
    }
    fn coord(&self, _: u32) -> [f64; 3] {
        unreachable!("placeholder mesh")
    }
    fn n_bulk(&self) -> usize {
        0
    }
    fn n_boundary(&self) -> usize {
        0
    }
    fn element(&self, _: elmer_core::ElemId) -> elmer_core::ElementRef<'_> {
        unreachable!("placeholder mesh")
    }
}

fn no_mesh(_: &Path, _: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    Ok(Box::new(NoMesh))
}

fn unsupported(msg: impl Into<String>) -> SolverError {
    SolverError::Unsupported(msg.into())
}

fn runs_never(params: &ParamMap) -> bool {
    params.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
}

/// What this composition cannot run in parallel, refused before any rank
/// reads its partition. The serial checks of `elmer-run` that are public
/// ([`elmer_run::unported::check`]) apply too.
fn preflight(model: &Model, registry: &elmer_solver::SolverRegistry) -> Result<(), SolverError> {
    let sim = &model.simulation;
    // The solver sections go to the Simulation as given; the parallel
    // switches are applied once the variables exist (`prepare_slot`).
    crate::levels::splits(sim, &model.boundary_conditions)?;
    for key in ["post file", "output file", "coordinate transformation", "extruded mesh levels"] {
        if sim.get(key).is_some() {
            return Err(unsupported(format!("Simulation {key} in parallel")));
        }
    }
    match sim.get("coordinate mapping") {
        None => {}
        Some(Value::Integer(v)) if v[..] == [1, 2, 3] => {}
        Some(v) => return Err(unsupported(format!("Simulation: Coordinate Mapping {v:?}"))),
    }
    if sim.get("coordinate scaling").is_some() && sim.real("coordinate scaling").is_none() {
        return Err(unsupported("Coordinate Scaling that is not a constant"));
    }
    for (i, p) in model.equations.iter().enumerate() {
        if p.string("element").is_some_and(|e| normalize_key(e) != "n:1") {
            return Err(unsupported(format!("Equation {}: Element (p-, bubble, edge or DG elements)", i + 1)));
        }
    }
    for (i, p) in model.solvers.iter().enumerate() {
        let id = i + 1;
        if p.0.is_empty() || runs_never(p) {
            continue;
        }
        let Some((lib, entry)) = procedure_of(p) else {
            return Err(unsupported(format!("Solver {id}: no Procedure")));
        };
        if elmer_run::output::is_result_output(p) {
            return Err(unsupported(format!("Solver {id}: ResultOutputSolver in parallel")));
        }
        if registry.get(&lib, &entry).is_none() {
            return Err(unsupported(format!("Solver {id}: Procedure \"{lib}\" \"{entry}\"")));
        }
        if p.string("element").is_some_and(|e| normalize_key(e) != "n:1") {
            return Err(unsupported(format!("Solver {id}: Element (p-, bubble, edge or DG elements)")));
        }
        for key in ["linear system complex", "coupled solver", "block solver", "monolithic slave"] {
            if p.logical(key) == Some(true) {
                return Err(unsupported(format!("Solver {id}: {key} in parallel")));
            }
        }
        let Some(kind) = p.string("linear system solver").map(normalize_key) else { continue };
        if kind != "iterative" {
            return Err(unsupported(format!("Solver {id}: Linear System Solver = {kind} in parallel")));
        }
        let method = method_of(p);
        if !PARALLEL_METHODS.contains(&method.as_str()) && !ITERMETHODS.contains(&method.as_str()) {
            return Err(unsupported(format!("Solver {id}: Linear System Iterative Method = {method} in parallel")));
        }
        dist_options(p)?;
        if ITERMETHODS.contains(&method.as_str()) {
            elmer_itermethods::IterMethodsConfig::from_params(p)?;
        }
    }
    elmer_run::unported::check(model)
}

/// Rank-local errors become an error on every rank, so no rank waits in a
/// collective another one never reaches.
fn agree<T>(comm: &dyn Communicator, local: Result<T, SolverError>) -> Result<T, SolverError> {
    let mut flag = [if local.is_err() { 1.0 } else { 0.0 }];
    comm.max_f64(&mut flag);
    match local {
        Err(e) => Err(e),
        Ok(_) if flag[0] != 0.0 => Err(SolverError::fatal("LoadMesh", "another rank could not load its partition")),
        Ok(v) => Ok(v),
    }
}

/// This rank's mesh, refusing element types `elmer-fe` does not define.
fn load_mesh(mesh_dir: &Path, rank: usize, np: usize, scale: f64) -> Result<ParMesh, SolverError> {
    let mesh = if np == 1 { ParMesh::read_serial(mesh_dir, scale) } else { ParMesh::read(mesh_dir, rank, np, scale) }
        .map_err(|e| SolverError::fatal("LoadMesh", e))?;
    if let Some(t) = mesh.element_types().into_iter().find(|&t| elmer_fe::element_type(t).is_err()) {
        return Err(unsupported(format!("element type {t} is not an elmer-fe Lagrange type")));
    }
    Ok(mesh)
}

/// Runs the case `sif` on the ranks of `ex`, each reading
/// `<Mesh DB>/partitioning.<np>/part.<rank+1>.*` (the mesh DB itself at
/// np=1). Collective. Unsupported input comes back as
/// [`SolverError::Unsupported`] before any file is written. Setup and run
/// errors are not yet agreed between ranks: one that only some ranks hit
/// leaves the others waiting, and the binary's `mpirun` then needs a kill.
pub fn run_partitioned(sif: &Path, ex: Box<dyn Exchange>, log: &dyn Logger) -> Result<ParOutcome, SolverError> {
    let (rank, np) = (ex.rank(), ex.size());
    let loaded = SifFrontend { mesh: no_mesh }.load(sif).map_err(SolverError::Unsupported)?;
    if !loaded.reloads.is_empty() {
        return Err(unsupported("RUN blocks in parallel"));
    }
    let mut model = loaded.model;
    let case_dir: PathBuf = loaded.case_dir;
    let registry = elmer_run::solver_registry();
    preflight(&model, &registry)?;
    // A solver that never runs and owns no variable is left out, as a gap in
    // the numbering: its procedure would still be built at setup.
    for p in &mut model.solvers {
        if runs_never(p) && p.get("variable").is_none() {
            p.0.clear();
        }
    }

    let scale = model.simulation.real("coordinate scaling").unwrap_or(1.0);
    let mut mesh = agree(&*ex, load_mesh(&loaded.mesh_dir, rank, np, scale))?;
    for _ in 0..crate::levels::splits(&model.simulation, &model.boundary_conditions)? {
        mesh = crate::levels::split(&mesh, &*ex)?;
    }
    let coords = elmer_run::coordinate_system(&model, &mesh)?;
    elmer_run::install_fe(coords.dim)?;
    let n_ghost = mesh.n_ghost();
    let (ctx, _guard) = RankContext::install(ex, mesh.node_gids().to_vec(), mesh.node_owners().to_vec());
    let comm: &dyn Communicator = &*ctx.ex;

    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(SifMatc));
    sim.backends.assembly = Some(Box::new(ParAssembly));
    sim.backends.linear = linear_registry();
    sim.setup(&mesh, &registry, log, comm)?;
    for slot in &mut sim.slots {
        let Some(v) = slot.variable else { continue };
        let dofs = sim.vars.by_id(v).dofs;
        prepare_slot(&mut slot.params, dofs, n_ghost)?;
    }
    let report = sim.run(&mesh, log, comm)?;
    let passed = finalize_reference(&sim.reference, &case_dir, log, comm)?;
    let norms = sim
        .slots
        .iter()
        .filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm)))
        .collect();
    Ok(ParOutcome { report, norms, passed })
}

/// The partitioned run behind `elmer-cli`'s [`Engine`]: one run per process
/// over the transport it was built with.
///
/// `elmer-cli`'s driver still refuses native runs under more than one MPI
/// rank before it asks the engine; until it gains a way to ask the engine
/// first, the `elmer-mpi-run` binary is the entry point.
pub struct MpiEngine {
    ex: Mutex<Option<Box<dyn Exchange>>>,
}

impl MpiEngine {
    pub fn new(ex: Box<dyn Exchange>) -> Self {
        MpiEngine { ex: Mutex::new(Some(ex)) }
    }

    /// Runs `sif` and maps the result to the process outcome, printing a
    /// fatal error the way `elmer-cli` does.
    pub fn run_sif(&self, sif: &Path, log: &dyn Logger, fatal: &dyn Fn(&str, &str)) -> EngineOutcome {
        let Some(ex) = self.ex.lock().unwrap_or_else(|e| e.into_inner()).take() else {
            return EngineOutcome::Refused(vec!["one partitioned run per process".into()]);
        };
        match run_partitioned(sif, ex, log) {
            Ok(_) => EngineOutcome::Done(ExitCode::Pass),
            Err(
                e @ (SolverError::Unsupported(_)
                | SolverError::NotWired(_)
                | SolverError::Linalg(LinalgError::Unsupported(_))),
            ) => EngineOutcome::Refused(vec![e.to_string()]),
            Err(e) => {
                let (caller, msg) = match &e {
                    SolverError::Fatal { caller, msg } | SolverError::Numerical { caller, msg } => (caller.as_str(), msg.clone()),
                    other => ("MAIN", other.to_string()),
                };
                fatal(caller, &msg);
                EngineOutcome::Done(ExitCode::Fail)
            }
        }
    }
}

impl Engine for MpiEngine {
    fn run(&self, req: &RunRequest<'_>, log: &ElmerLogger) -> EngineOutcome {
        self.run_sif(req.sif, log, &|caller, msg| log.fatal(caller, msg))
    }
}
