//! [`NativeEngine`]: the native solver behind `elmer-cli`'s [`Engine`] trait.

use std::path::{Path, PathBuf};

use elmer_cli::driver::{Engine, EngineOutcome, RunRequest};
use elmer_cli::{ElmerLogger, ExitCode};
use elmer_core::{LinalgError, MeshView, SerialComm};
use elmer_solver::{Model, SolverError};

use crate::case::{run_case, Case};

/// A case read from disk.
pub struct LoadedCase {
    pub model: Model,
    pub mesh: Box<dyn MeshView>,
    pub case_dir: PathBuf,
    pub mesh_dir: PathBuf,
}

/// The `.sif` parser and mesh reader (S1's `elmer-sif` and `elmer-mesh`).
pub trait Frontend {
    /// Parses `sif` and loads its `Mesh DB`. An `Err` is a refusal reason.
    fn load(&self, sif: &Path) -> Result<LoadedCase, String>;
}

/// Runs cases natively when a [`Frontend`] is installed and refuses them
/// otherwise.
#[derive(Default)]
pub struct NativeEngine {
    frontend: Option<Box<dyn Frontend>>,
}

impl NativeEngine {
    pub fn new(frontend: Box<dyn Frontend>) -> Self {
        NativeEngine { frontend: Some(frontend) }
    }
}

impl Engine for NativeEngine {
    fn run(&self, req: &RunRequest<'_>, log: &ElmerLogger) -> EngineOutcome {
        let Some(frontend) = &self.frontend else {
            return EngineOutcome::Refused(vec!["no .sif and mesh front end (elmer-sif, elmer-mesh) in this build".into()]);
        };
        let loaded = match frontend.load(req.sif) {
            Ok(c) => c,
            Err(e) => return EngineOutcome::Refused(vec![e]),
        };
        let case = Case {
            model: loaded.model,
            mesh: loaded.mesh.as_ref(),
            case_dir: loaded.case_dir,
            mesh_dir: loaded.mesh_dir,
        };
        match run_case(case, log, &SerialComm) {
            // Fortran exits 0 whatever the reference verdict; TEST.PASSED carries it.
            Ok(_) => EngineOutcome::Done(ExitCode::Pass),
            Err(
                e @ (SolverError::Unsupported(_)
                | SolverError::NotWired(_)
                | SolverError::Linalg(LinalgError::Unsupported(_))),
            ) => {
                EngineOutcome::Refused(vec![e.to_string()])
            }
            Err(e) => {
                let (caller, msg) = match &e {
                    SolverError::Fatal { caller, msg } | SolverError::Numerical { caller, msg } => {
                        (caller.as_str(), msg.clone())
                    }
                    other => ("MAIN", other.to_string()),
                };
                log.fatal(caller, &msg);
                EngineOutcome::Done(ExitCode::Fail)
            }
        }
    }
}
