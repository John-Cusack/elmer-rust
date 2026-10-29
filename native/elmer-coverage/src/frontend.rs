//! The front end of the harness: `elmer-run`'s [`SifFrontend`] (`elmer-sif`
//! parse and completion) with the test-local mesh loader and basis, plus the
//! same model read without a mesh for the static screen.

use std::path::{Path, PathBuf};

use elmer_run::output::take_vtu_output;
use elmer_run::sif::to_model;
use elmer_run::SifFrontend;
use elmer_solver::{Model, SolverError};

use crate::standin;

/// The front end `NativeEngine` runs the quick tests through.
pub fn frontend() -> SifFrontend {
    SifFrontend { mesh: standin::mesh::load }
}

/// The model `elmer-run` would run `sif` with, before the mesh is read and
/// after the `Post File` auto-solver has been taken out. The second value is
/// a refusal reason when that step already declines the case.
pub fn read_model(sif: &Path) -> Result<(Model, Option<String>), String> {
    let path = sif.to_str().ok_or("non-UTF-8 case path")?;
    let mut parsed = elmer_sif::parse_file(path).map_err(|e| e.to_string())?;
    let blocks = parsed.sif.runs.iter().skip(1).filter(|r| r.header.present || r.simulation.is_some() || !r.solvers.is_empty()).count();
    elmer_sif::complete(&mut parsed.sif, path)?;
    let merged = parsed.sif.merged();
    let mut model = to_model(&merged)?;
    let case_dir = sif.parent().filter(|p| !p.as_os_str().is_empty()).map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from("."));
    let mesh_dir = match &merged.header.mesh_db {
        Some((db, name)) => case_dir.join(db).join(name),
        None => case_dir.clone(),
    };
    let mut refusal = match take_vtu_output(&mut model, &case_dir, &mesh_dir) {
        Ok(_) => None,
        Err(SolverError::Unsupported(m) | SolverError::NotWired(m)) => Some(m),
        Err(e) => Some(e.to_string()),
    };
    if blocks > 0 {
        refusal.get_or_insert_with(|| format!("{} RUN blocks (one Simulation per process is composed)", blocks + 1));
    }
    Ok((model, refusal))
}
