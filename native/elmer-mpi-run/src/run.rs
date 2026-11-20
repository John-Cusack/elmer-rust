//! One case on every rank: front end, parallel preflight, `Simulation` with
//! the parallel backends, then `CompareToReferenceSolution`.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_cli::driver::{Engine, EngineOutcome, RunRequest};
use elmer_cli::{ElmerLogger, ExitCode};
use elmer_core::{normalize_key, Communicator, ElemId, LinalgError, Logger, MeshView, ParamMap, Params, Value, NO_DOF};
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
    if !model.distributed.is_empty() {
        return Err(unsupported("-distribute keywords in parallel (the entity weights' reduction)"));
    }
    // `Post File` is a ResultOutputSolver section by now (`prepare_vtu_output`).
    for key in ["post file", "coordinate transformation", "extruded mesh levels"] {
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
        // VTU output, one piece per partition: checked by `prepare_vtu_output`.
        if elmer_run::output::is_result_output(p) {
            continue;
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
/// `MatrixParallelInfoDG` (ParallelUtils.F90:341): every DG solver's matrix
/// with the global dof `8 (GElementIndex - 1) + j` of its element's node `j`
/// (`MaxDGNodes = 8`), owned by the element's partition.
/// [`elmer_solver::NodeInterface`] over the thread's [`RankContext`].
struct ParNodes;

impl elmer_solver::NodeInterface for ParNodes {
    fn sum_shared(&self, v: &mut [f64]) -> Result<(), SolverError> {
        let ctx = RankContext::current().ok_or_else(|| unsupported("a nodal exchange outside a partitioned run"))?;
        ctx.sum_shared_nodes(v).map_err(|e| SolverError::Unsupported(e.to_string()))
    }

    fn min_shared(&self, v: &mut [f64]) -> Result<(), SolverError> {
        let ctx = RankContext::current().ok_or_else(|| unsupported("a nodal exchange outside a partitioned run"))?;
        ctx.pick_shared_nodes(v, f64::min).map_err(|e| SolverError::Unsupported(e.to_string()))
    }

    fn max_shared(&self, v: &mut [f64]) -> Result<(), SolverError> {
        let ctx = RankContext::current().ok_or_else(|| unsupported("a nodal exchange outside a partitioned run"))?;
        ctx.pick_shared_nodes(v, f64::max).map_err(|e| SolverError::Unsupported(e.to_string()))
    }
}

/// `MatrixParallelInfoDG` (ParallelUtils.F90:341-511) for every DG solver:
/// the global id and owner of each row. Without `DG Reduced Basis` index
/// `j` of an element is `8 (GElementIndex - 1) + j`; with it, a row is the
/// node's `GlobalDOFs` offset by its body group times the largest node id.
/// A row of an element whose nodes are all interface nodes is shared with
/// the partitions of its node (the element's own first when not reduced);
/// any other row is private to this rank. Collective.
fn register_dg_matrices(sim: &mut Simulation, mesh: &ParMesh, ctx: &RankContext) -> Result<(), SolverError> {
    if mesh.element_gids().is_empty() {
        return Ok(());
    }
    let mut slots = Vec::new();
    let rank = ctx.ex.rank();
    let parts = mesh.node_parts();
    let node_gid = mesh.node_gids();
    let mut maxnode = [node_gid.iter().copied().max().unwrap_or(0) as f64];
    ctx.ex.max_f64(&mut maxnode);
    let maxnode = maxnode[0] as u64;
    for ix in 0..sim.slots.len() {
        let slot = &sim.slots[ix];
        if slot.params.logical("discontinuous galerkin") != Some(true) {
            continue;
        }
        let (Some(d), Some(id)) = (slot.solver_dofs(), slot.variable) else { continue };
        let reduced = slot.params.logical("dg reduced basis") == Some(true);
        let ints = |key: &str| match slot.params.get(key) {
            Some(elmer_core::Value::Integer(v)) => Some(v.clone()),
            _ => None,
        };
        let (map, master) = (ints("dg reduced basis mapping"), ints("dg reduced basis master bodies"));
        if reduced && master.is_some() && ints("dg reduced basis slave bodies").is_some() {
            return Err(unsupported("DG Reduced Basis Master and Slave Bodies in a partitioned run"));
        }
        let dofs = sim.vars.by_id(id).dofs;
        let a = d.create_matrix(dofs);
        let mut gid = vec![0u64; a.nrows];
        let mut owner = vec![rank; a.nrows];
        for t in 0..mesh.n_bulk() {
            let el = mesh.element(ElemId(t as u32));
            let all_interface = el.nodes.iter().all(|&n| parts[n as usize].len() > 1);
            let group = if reduced {
                let body = sim.model.body_of_tag(el.tag).unwrap_or(0) as i32;
                (match (&master, &map) {
                    (Some(m), _) => if m.contains(&body) { 1 } else { 2 },
                    (None, Some(m)) => m.get((body - 1).max(0) as usize).copied().unwrap_or(0),
                    (None, None) => body,
                } - 1) as u64
            } else {
                0
            };
            d.element_slots(ElemId(t as u32), &mut slots);
            for (j, &s) in slots.iter().enumerate() {
                if s == NO_DOF {
                    continue;
                }
                let node = el.nodes[j] as usize;
                let base = if reduced {
                    group * maxnode + node_gid[node]
                } else {
                    8 * (mesh.element_gids()[t] - 1) + j as u64 + 1
                };
                for k in 0..dofs {
                    let row = s as usize * dofs + k;
                    if row >= a.nrows {
                        continue;
                    }
                    gid[row] = dofs as u64 * (base - 1) + k as u64;
                    if all_interface {
                        let nl = &parts[node];
                        let own = mesh.element_owners()[t];
                        owner[row] = if !reduced && nl.contains(&own) { own } else { nl[0] };
                    }
                }
            }
        }
        let owned: std::sync::Arc<[bool]> = owner.iter().map(|&o| o == rank).collect();
        ctx.register_rows(gid, owner, &a);
        sim.slots[ix].owned_rows = Some(owned);
    }
    Ok(())
}

pub fn run_partitioned(sif: &Path, ex: Box<dyn Exchange>, log: &dyn Logger) -> Result<ParOutcome, SolverError> {
    let (rank, np) = (ex.rank(), ex.size());
    let loaded = SifFrontend { mesh: no_mesh }.load(sif).map_err(SolverError::Unsupported)?;
    if loaded.model.simulation.logical("single mesh") == Some(true) {
        return run_single_mesh(loaded, &*ex, log);
    }
    if !loaded.reloads.is_empty() {
        return Err(unsupported("RUN blocks in parallel"));
    }
    let mut model = loaded.model;
    let case_dir: PathBuf = loaded.case_dir;
    // The case's test-local plugin ports, as elmer-run's serial path.
    let fixture_case = elmer_run::fixture_case(&case_dir);
    let registry = elmer_run::solver_registry_for(&fixture_case);
    // VTU output (`Post File` and ResultOutputSolver sections): each rank
    // writes the piece of its partition, `<base>_<np>np<rank+1>`, and rank 0
    // the `.pvtu` (VtuOutputSolver.F90:622-633), into the mesh directory.
    let mut outputs = elmer_run::output::prepare_vtu_output(&mut model, &case_dir, &loaded.mesh_dir, np > 1)?;
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
    elmer_run::case::legacy_solvers(&mut model, coords);
    elmer_run::install_fe(coords.dim)?;
    let (ctx, _guard) = RankContext::install(ex, mesh.node_gids().to_vec(), mesh.node_owners().to_vec());
    let comm: &dyn Communicator = &*ctx.ex;
    if np > 1 && model.solvers.iter().any(|p| p.logical("discontinuous galerkin") == Some(true)) {
        if mesh.dim() != 2 {
            return Err(unsupported("Discontinuous Galerkin on a partitioned 3D mesh (SParFaceNumbering)"));
        }
        let edges = elmer_solver::dg_edges(&mesh)?;
        let keys = crate::edges::edge_numbers(&*ctx.ex, &mesh, &edges).map_err(SolverError::Unsupported)?;
        mesh.set_edge_order(keys);
    }

    // The partition's node coordinates, movable by the solvers that assign
    // `Mesh % Nodes` (mesh mappers, SunAngle's heights), as the serial run's.
    let moving = elmer_run::MovableMesh::new(&mesh);
    let _motion = moving.install();
    elmer_run::elmer_fixtures::registry::init_sections(&fixture_case, &mut model);
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(SifMatc));
    sim.backends.assembly = Some(Box::new(ParAssembly));
    sim.backends.linear = linear_registry();
    // Basis and quadrature for the simulation's own integrals (component
    // reductions), as the serial run installs them.
    sim.backends.fe = Some(elmer_run::fe_backend(coords.dim));
    // The p and DG dofs, as the serial composition installs them.
    let degree: Arc<elmer_run::elmer_pdofs::solver::MatcFn> = Arc::new(|f: &str, cx: [f64; 4]| {
        use elmer_solver::elmer_expr::MatcEngine;
        elmer_run::SifMatc.eval_with("cx", f, &cx)?.first().copied().ok_or_else(|| format!("MATC {f}: no value"))
    });
    sim.backends.dofs = Some(Box::new(elmer_run::elmer_pdofs::PDofsBackend::new().with_matc(degree)));
    sim.backends.nodes = Some(Box::new(ParNodes));
    elmer_run::elmer_fixtures::registry::install_live(&fixture_case, &mut sim, &moving);
    if let Some(h) = mesh.halo_elements() {
        sim.set_halo_elements(h.into());
    }
    // `Apply Conforming BCs`: each rank's ConformingNodePerm /
    // ConformingEdgePerm of its own partition, as PeriodicPermutation
    // makes them.
    sim.backends.conforming = elmer_run::conforming_map(&sim.model, &mesh, coords.dim, log)?.map(Arc::new);
    sim.setup(&moving, &registry, log, comm)?;
    // `InitialPermutation` numbers the nodes of halo elements as well; only
    // the DG numbering follows it here, so a solver numbered by nodes is
    // refused.
    if mesh.halo_elements().is_some() {
        if let Some(s) = sim.slots.iter().find(|s| s.perm.is_some() && s.params.logical("discontinuous galerkin") != Some(true)) {
            let name = procedure_of(&s.params).map_or_else(|| format!("solver {}", s.id), |(_, entry)| entry);
            return Err(unsupported(format!("{name} without Discontinuous Galerkin on a halo partitioning")));
        }
    }
    register_dg_matrices(&mut sim, &mesh, &ctx)?;
    if let Some(writer) = elmer_run::result_file::writer(&sim.model, &loaded.mesh_dir, rank, np, log) {
        sim.set_save_current(Box::new(writer));
    }
    sim.set_node_owned(mesh.owned_nodes().into());
    for out in &mut outputs {
        out.select_fields(&sim, mesh.dim())?;
    }
    let written = Arc::new(Mutex::new(Vec::new()));
    for out in outputs {
        let id = out.index + 1;
        let slot = sim
            .slots
            .iter_mut()
            .find(|s| s.id == id)
            .ok_or_else(|| SolverError::fatal("ResultOutputSolver", format!("no slot for Solver {id}")))?;
        slot.set_instance(Box::new(out.into_solver(Arc::clone(&written)).in_job(rank, np)));
    }
    for slot in &mut sim.slots {
        if slot.variable.is_some() {
            prepare_slot(&mut slot.params);
        }
    }
    let report = sim.run(&moving, log, comm)?;
    let passed = finalize_reference(&sim.reference, &case_dir, log, comm)?;
    let norms = sim
        .slots
        .iter()
        .filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm)))
        .collect();
    Ok(ParOutcome { report, norms, passed })
}

/// `Single Mesh = True`: every rank reads the whole mesh DB and runs the
/// case serially through `elmer-run` (`Mesh % SingleMesh` turns Elmer's
/// parallel paths off), as one rank of the job: rank 0 writes
/// `TEST.PASSED_<np>`, VTU pieces are named per rank, and `Parallel Slices`
/// gives each rank its own slice. Not reproduced: the maximum of the
/// nonlinear change over slices (`ComputeChange` with `Number Of Slices >
/// 1`), which only matters when a slice would otherwise stop iterating
/// earlier than another.
fn run_single_mesh(loaded: elmer_run::LoadedCase, job: &dyn Exchange, log: &dyn Logger) -> Result<ParOutcome, SolverError> {
    let scale = loaded.model.simulation.real("coordinate scaling").unwrap_or(1.0);
    let mesh = agree(job, ParMesh::read_serial(&loaded.mesh_dir, scale).map_err(|e| SolverError::fatal("LoadMesh", e)))?;
    let case = elmer_run::Case {
        model: loaded.model,
        reloads: loaded.reloads,
        mesh: &mesh,
        case_dir: loaded.case_dir,
        mesh_dir: loaded.mesh_dir,
    };
    let job: &dyn Communicator = job;
    let out = elmer_run::run_case_with(case, log, &elmer_run::RunEnv { job, globals: Vec::new() })?;
    Ok(ParOutcome { report: out.report, norms: out.norms, passed: out.passed })
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
