//! One case through the native path: preflight, `Simulation` setup and run,
//! `Post File` output, then `CompareToReferenceSolution`'s `TEST.PASSED`.

use std::path::PathBuf;
use std::rc::Rc;

use elmer_core::{Communicator, Logger, MeshView, Params, Value};
use elmer_fe::{CoordinateKind, CoordinateSystem};
use elmer_solver::{
    finalize_reference, solver::procedure_of, Model, RunReport, Simulation, SolverError, SolverRegistry,
};

use crate::compose::{coordinate_system, install_fe, solver_registry};
use crate::matc::SifMatc;
use crate::output::take_vtu_output;

/// A parsed case and its mesh.
pub struct Case<'a> {
    pub model: Model,
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
        let Some((lib, entry)) = procedure_of(params) else {
            return Err(SolverError::Unsupported(format!("Solver {}: no Procedure", i + 1)));
        };
        if registry.get(&lib, &entry).is_none() {
            return Err(SolverError::Unsupported(format!("Solver {}: Procedure \"{lib}\" \"{entry}\"", i + 1)));
        }
    }
    crate::unported::check(model)
}

/// Runs `case` serially. Unsupported input comes back as
/// [`SolverError::Unsupported`] or [`SolverError::NotWired`] before any
/// file is written.
pub fn run_case(case: Case<'_>, log: &dyn Logger, comm: &dyn Communicator) -> Result<Outcome, SolverError> {
    if comm.size() > 1 {
        return Err(SolverError::Unsupported("more than one MPI rank".into()));
    }
    let Case { mut model, mesh, case_dir, mesh_dir } = case;
    let registry = solver_registry();
    let mut post = take_vtu_output(&mut model, &case_dir, &mesh_dir)?;
    let coords = coordinate_system(&model, mesh)?;
    legacy_solvers(&mut model, coords);
    preflight(&model, &registry)?;
    install_fe(coords.dim)?;

    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(SifMatc));
    sim.setup(mesh, &registry, log, comm)?;
    if let Some(post) = &mut post {
        post.select_fields(&sim)?;
    }
    let report = sim.run(mesh, log, comm)?;

    let mut written = Vec::new();
    if let Some(post) = post {
        written.push(post.write(&sim, mesh)?);
    }
    let passed = finalize_reference(&sim.reference, &case_dir, log, comm)?;
    let norms = sim
        .slots
        .iter()
        .filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm)))
        .collect();
    Ok(Outcome { report, norms, passed, written })
}
