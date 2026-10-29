//! One case through the native path: preflight, `Simulation` setup and run,
//! `Post File` output, then `CompareToReferenceSolution`'s `TEST.PASSED`.

use std::path::PathBuf;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_core::{normalize_key, Communicator, Logger, MeshView, ParamMap, Params, Value};
use elmer_fe::{CoordinateKind, CoordinateSystem};
use elmer_solver::{
    finalize_reference, solver::procedure_of, Model, RunReport, Simulation, SolverError, SolverRegistry,
};

use crate::compose::{coordinate_system, fe_backend, install_fe, solver_registry};
use crate::matc::SifMatc;
use crate::motion::MovableMesh;
use crate::output::{is_result_output, prepare_vtu_output, MappedTags};
use elmer_mesh_levels::{extrude, mesh_levels, ExtrudeOptions, Mesh as LevelsMesh, MeshLevelsError};

/// A parsed case and its mesh.
pub struct Case<'a> {
    pub model: Model,
    /// Further `RUN` blocks (`Simulation::reload`), each with only its own
    /// entries; empty for a single simulation.
    pub reloads: Vec<Model>,
    pub mesh: &'a dyn MeshView,
    /// Where `TEST.PASSED` goes and a relative `Output Directory` starts.
    pub case_dir: PathBuf,
    /// The `Mesh DB` directory, the default VTU output directory.
    pub mesh_dir: PathBuf,
}

/// What a finished run produced.
#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    pub report: RunReport,
    /// `(solver id, ComputeNorm)` of every solver with a primary variable.
    pub norms: Vec<(usize, f64)>,
    /// The `TEST.PASSED` verdict; `None` when no solver has a reference.
    pub passed: Option<bool>,
    /// Output files written.
    pub written: Vec<PathBuf>,
}

/// The procedure `AddEquationBasics` gives a solver section whose
/// `Equation` is the historical solver `equation`, for the ones this build
/// completes (see [`legacy_solvers`]).
pub fn legacy_procedure(equation: &str) -> Option<(&'static str, &'static str)> {
    match equation.trim().to_ascii_lowercase().as_str() {
        "heat equation" => Some(("HeatSolve", "HeatSolver")),
        "stress analysis" => Some(("StressSolve", "StressSolver")),
        "navier-stokes" => Some(("FlowSolve", "FlowSolver")),
        _ => None,
    }
}

/// The historical solvers `AddEquationBasics` (MainUtils.F90:1221-1305)
/// completes from the solver's `Equation` name when it has no `Procedure`,
/// for the ones whose procedure this build has. Others stay without a
/// procedure and are refused by [`preflight`].
fn legacy_solvers(model: &mut Model, coords: CoordinateSystem) {
    let cdim = coords.dim;
    for params in &mut model.solvers {
        if params.get("procedure").is_some() {
            continue;
        }
        let name = params.string("equation").map(|s| s.trim().to_ascii_lowercase());
        let (var, dofs, procedure) = match name.as_deref() {
            Some("heat equation") => {
                if params.get("radiation solver").is_none() {
                    params.insert("Radiation Solver", Value::Logical(true));
                }
                ("Temperature", 1, "HeatSolve HeatSolver")
            }
            Some("stress analysis") => ("Displacement", cdim, "StressSolve StressSolver"),
            // `dofs` is reset to 1: the components are in the name. The 1e-6
            // start value is applied by elmer-physics-flow.
            Some("navier-stokes") => {
                let velocity = cdim + usize::from(coords.kind == CoordinateKind::CylindricSymmetric);
                if velocity == 3 {
                    ("Flow Solution[Velocity:3 Pressure:1]", 1, "FlowSolve FlowSolver")
                } else {
                    ("Flow Solution[Velocity:2 Pressure:1]", 1, "FlowSolve FlowSolver")
                }
            }
            _ => continue,
        };
        if params.get("variable").is_none() {
            params.insert("Variable", Value::Str(var.into()));
            if dofs > 1 {
                params.insert("Variable Dofs", Value::Integer(vec![dofs as i32]));
            }
        }
        params.insert("Procedure", Value::Str(procedure.into()));
    }
}

/// Refuses solver sections this build cannot run: `Simulation` skips a
/// procedure missing from the registry, which would drop a solver silently.
/// `Exec Solver = Never` sections are exempt. Then refuses the model
/// features no provider implements ([`crate::unported`]).
fn preflight(model: &Model, registry: &SolverRegistry) -> Result<(), SolverError> {
    for (i, params) in model.solvers.iter().enumerate() {
        // A section that never runs needs no implementation (heateq_bdf2's SaveScalars).
        if params.0.is_empty() || params.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never")) {
            continue;
        }
        // Checked and served by [`crate::output`].
        if is_result_output(params) {
            continue;
        }
        let Some((lib, entry)) = procedure_of(params) else {
            return Err(SolverError::Unsupported(format!("Solver {}: no Procedure", i + 1)));
        };
        if registry.get(&lib, &entry).is_none() {
            return Err(SolverError::Unsupported(format!("Solver {}: Procedure \"{lib}\" \"{entry}\"", i + 1)));
        }
    }
    crate::unported::check(model)
}

/// HeatSolve convects relative to the mesh (`U - MU` with the `Mesh
/// Velocity` variable, HeatSolve.F90); `elmer-physics-heat` uses the
/// material velocity alone, so a convected heat equation on a moving mesh
/// is refused once the variables exist.
fn check_mesh_velocity(sim: &Simulation) -> Result<(), SolverError> {
    if sim.vars.id("mesh velocity").is_none() {
        return Ok(());
    }
    let heat = sim.model.solvers.iter().any(|p| {
        procedure_of(p).is_some_and(|(lib, _)| lib.eq_ignore_ascii_case("heatsolve"))
            && !p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
    });
    let convected =
        sim.model.equations.iter().any(|p| p.string("convection").is_some_and(|c| !c.eq_ignore_ascii_case("none")));
    if heat && convected {
        return Err(SolverError::Unsupported(
            "HeatSolve: convection relative to the Mesh Velocity (moving mesh) is not implemented".into(),
        ));
    }
    Ok(())
}

/// Section-wise merge of `block` into `model`, as `Simulation::reload` does
/// (`ReloadInputFile`): entries replace same-name ones, new sections append.
fn merged_with(model: &Model, block: &Model) -> Model {
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
    let mut m = model.clone();
    merge(&mut m.simulation, &block.simulation);
    merge(&mut m.constants, &block.constants);
    merge_vec(&mut m.solvers, &block.solvers);
    merge_vec(&mut m.equations, &block.equations);
    merge_vec(&mut m.bodies, &block.bodies);
    merge_vec(&mut m.materials, &block.materials);
    merge_vec(&mut m.body_forces, &block.body_forces);
    merge_vec(&mut m.initial_conditions, &block.initial_conditions);
    merge_vec(&mut m.boundary_conditions, &block.boundary_conditions);
    m
}

/// Refuses, before anything runs, what a later `RUN` block would make
/// unsupported. Fortran solvers read their keywords on every visit, while
/// several native solvers read their section once when built; until the
/// slots are rebuilt on reload, a block may only change `Reference Norm` /
/// `Reference Norm Tolerance` of existing solvers.
fn check_reloads(model: &Model, reloads: &[Model], registry: &SolverRegistry) -> Result<(), SolverError> {
    const ALLOWED: [&str; 2] = ["reference norm", "reference norm tolerance"];
    let mut merged = model.clone();
    for (k, block) in reloads.iter().enumerate() {
        let run = k + 2;
        let refuse = |what: String| {
            Err(SolverError::Unsupported(format!(
                "RUN block {run} changes {what} (re-reading changed keywords after a reload is not composed)"
            )))
        };
        if let Some(key) = block.simulation.0.keys().find(|k| *k != "solver input file") {
            return refuse(format!("Simulation {key}"));
        }
        let others = [
            ("Constants", std::slice::from_ref(&block.constants)),
            ("Equation", &block.equations[..]),
            ("Body", &block.bodies[..]),
            ("Material", &block.materials[..]),
            ("Body Force", &block.body_forces[..]),
            ("Initial Condition", &block.initial_conditions[..]),
            ("Boundary Condition", &block.boundary_conditions[..]),
        ];
        for (name, sections) in others {
            if let Some((i, p)) = sections.iter().enumerate().find(|(_, p)| !p.0.is_empty()) {
                let key = p.0.keys().next().cloned().unwrap_or_default();
                return refuse(format!("{name} {} {key}", i + 1));
            }
        }
        for (i, params) in block.solvers.iter().enumerate() {
            if i >= model.solvers.len() && !params.0.is_empty() {
                return refuse(format!("Solver {} (a new solver)", i + 1));
            }
            if let Some(key) = params.0.keys().find(|k| !ALLOWED.contains(&k.as_str())) {
                return refuse(format!("Solver {} {key}", i + 1));
            }
        }
        merged = merged_with(&merged, block);
        preflight(&merged, registry)?;
    }
    Ok(())
}

/// Inputs no provider applies yet: restart files, mesh-load transforms
/// beyond `Coordinate Scaling`, `Mesh Levels` and extrusion, and element
/// definitions other than the Lagrange (`n:1`) basis `elmer-fe` evaluates.
fn check_mesh_and_elements(model: &Model) -> Result<(), SolverError> {
    let sim = &model.simulation;
    // `Restart File` sets initial fields from a .result file; nothing reads
    // it yet, and running without it is a different case.
    if sim.get("restart file").is_some() {
        return Err(SolverError::Unsupported("Simulation: Restart File (restart from a .result file)".into()));
    }
    if sim.get("coordinate transformation").is_some() {
        return Err(SolverError::Unsupported("Simulation: Coordinate Transformation".into()));
    }
    match sim.get("coordinate mapping") {
        None => {}
        Some(Value::Integer(v)) if v[..] == [1, 2, 3] => {}
        Some(v) => return Err(SolverError::Unsupported(format!("Simulation: Coordinate Mapping {v:?}"))),
    }
    let sections = model.solvers.iter().map(|p| ("Solver", p)).chain(model.equations.iter().map(|p| ("Equation", p)));
    for (i, (name, p)) in sections.enumerate() {
        if let Some(e) = p.string("element") {
            if normalize_key(e) != "n:1" {
                return Err(SolverError::Unsupported(format!(
                    "{name} section {}: Element = {e} (p-, bubble, edge or DG elements)",
                    i + 1
                )));
            }
        }
    }
    Ok(())
}

/// Elmer's mesh-load transforms, in Elmer's order: `Mesh Levels` splits
/// (`LoadModel`, ModelDescription.F90:3049), then `Extruded Mesh *`
/// extrusion (`CreateExtrudedMesh`, ElmerSolver.F90:897). Extrusion works on
/// the `MapBodiesAndBCs` tags and leaves body and BC indices as the new
/// mesh's tags, so the model's `Target Bodies` / `Target Boundaries` become
/// the identity. `None` when the mesh is used as loaded.
fn transform_mesh(model: &mut Model, mesh: &dyn MeshView) -> Result<Option<LevelsMesh>, SolverError> {
    let err = |e: MeshLevelsError| match e {
        MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshLevelsError::Invalid(m) => SolverError::fatal("LoadModel", m),
    };
    let split = mesh_levels(mesh, &model.simulation, &model.boundary_conditions).map_err(err)?;
    let current: &dyn MeshView = match &split {
        Some(m) => m,
        None => mesh,
    };
    let Some(opts) = ExtrudeOptions::from_params(&model.simulation, current.dim(), None).map_err(err)? else {
        return Ok(split);
    };
    let extruded = extrude(&MappedTags::new(current, model), &opts).map_err(err)?;
    for (i, body) in model.bodies.iter_mut().enumerate() {
        body.insert("Target Bodies", Value::Integer(vec![i as i32 + 1]));
    }
    for (i, bc) in model.boundary_conditions.iter_mut().enumerate() {
        bc.insert("Target Boundaries", Value::Integer(vec![i as i32 + 1]));
    }
    Ok(Some(extruded))
}

/// Runs `case` serially. Unsupported input comes back as
/// [`SolverError::Unsupported`] or [`SolverError::NotWired`] before any
/// file is written.
pub fn run_case(case: Case<'_>, log: &dyn Logger, comm: &dyn Communicator) -> Result<Outcome, SolverError> {
    if comm.size() > 1 {
        return Err(SolverError::Unsupported("more than one MPI rank".into()));
    }
    let Case { mut model, reloads, mesh, case_dir, mesh_dir } = case;
    let registry = solver_registry();
    let transformed = transform_mesh(&mut model, mesh)?;
    let mesh: &dyn MeshView = match &transformed {
        Some(m) => m,
        None => mesh,
    };
    let mesh = MovableMesh::new(mesh);
    let _motion = mesh.install();
    let mesh: &dyn MeshView = &mesh;
    let mut outputs = prepare_vtu_output(&mut model, &case_dir, &mesh_dir)?;
    let coords = coordinate_system(&model, mesh)?;
    legacy_solvers(&mut model, coords);
    preflight(&model, &registry)?;
    check_mesh_and_elements(&model)?;
    check_reloads(&model, &reloads, &registry)?;
    install_fe(coords.dim)?;

    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(fe_backend(coords.dim));
    sim.set_matc(Rc::new(SifMatc));
    sim.setup(mesh, &registry, log, comm)?;
    check_mesh_velocity(&sim)?;
    // Field lists are fixed once every variable exists, before the run, so a
    // refusal still comes before any file.
    for out in &mut outputs {
        out.select_fields(&sim)?;
    }
    let written = Arc::new(Mutex::new(Vec::new()));
    for out in outputs {
        let id = out.index + 1;
        let slot = sim
            .slots
            .iter_mut()
            .find(|s| s.id == id)
            .ok_or_else(|| SolverError::fatal("ResultOutputSolver", format!("no slot for Solver {id}")))?;
        slot.set_instance(Box::new(out.into_solver(Arc::clone(&written))));
    }
    let mut report = sim.run(mesh, log, comm)?;
    for block in &reloads {
        sim.reload(block, mesh, &registry, log, comm)?;
        report = sim.run(mesh, log, comm)?;
    }
    let written = std::mem::take(&mut *written.lock().unwrap());
    let passed = finalize_reference(&sim.reference, &case_dir, log, comm)?;
    let norms = sim
        .slots
        .iter()
        .filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm)))
        .collect();
    Ok(Outcome { report, norms, passed, written })
}
