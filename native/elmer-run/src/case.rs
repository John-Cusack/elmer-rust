//! One case through the native path: preflight, `Simulation` setup and run,
//! `Post File` output, then `CompareToReferenceSolution`'s `TEST.PASSED`.

use std::path::PathBuf;

use elmer_core::{Communicator, Logger, MeshView, Params};
use elmer_solver::{
    finalize_reference, solver::procedure_of, Model, RunReport, Simulation, SolverError, SolverRegistry,
};

use crate::compose::solver_registry;
use crate::output::take_post_file;

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

/// Refuses solver sections this build cannot run: `Simulation` skips a
/// procedure missing from the registry, which would drop a solver silently.
fn preflight(model: &Model, registry: &SolverRegistry) -> Result<(), SolverError> {
    for (i, params) in model.solvers.iter().enumerate() {
        if params.0.is_empty() {
            continue;
        }
        let Some((lib, entry)) = procedure_of(params) else {
            return Err(SolverError::Unsupported(format!("Solver {}: no Procedure", i + 1)));
        };
        if registry.get(&lib, &entry).is_none() {
            return Err(SolverError::Unsupported(format!("Solver {}: Procedure \"{lib}\" \"{entry}\"", i + 1)));
        }
    }
    // MapBodiesAndBCs fallbacks `Model::bc_of_tag` does not apply.
    for (i, bc) in model.boundary_conditions.iter().enumerate() {
        for k in ["default target", "default bc"] {
            if bc.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("Boundary Condition {}: {k}", i + 1)));
            }
        }
    }
    Ok(())
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
    preflight(&model, &registry)?;
    let mut post = take_post_file(&mut model, &case_dir, &mesh_dir)?;

    let mut sim = Simulation::new(model);
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
