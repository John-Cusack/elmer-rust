//! [`SifFrontend`]: `elmer-sif` as the case reader of [`crate::NativeEngine`].
//!
//! The `.sif` is parsed and completed by `elmer-sif` (`LoadInputFile`,
//! `CompleteModelKeywords`, the `Post File` auto-solver) and converted into
//! `elmer-solver`'s [`Model`]. The mesh reader is pluggable until
//! `elmer-mesh` merges.

use std::path::{Path, PathBuf};

use elmer_core::{MeshView, ParamMap};
use elmer_sif::{section_to_param_map, Section};
use elmer_solver::Model;

use crate::engine::{Frontend, LoadedCase};

/// Loads the `Mesh DB` directory after Elmer's mesh-load transforms
/// (`Coordinate Scaling`, ...), given the Simulation section.
pub type MeshLoader = fn(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String>;

/// `elmer-sif` plus a mesh loader.
pub struct SifFrontend {
    pub mesh: MeshLoader,
}

/// Sections by their 1-based id (`Body 3` at index 2); missing ids are
/// empty, as an unset section is in Elmer. A section without a number
/// (`elmer-sif` gives it id 0) has no place and is refused.
fn numbered(kind: &str, sections: &[Section]) -> Result<Vec<ParamMap>, String> {
    let n = sections.iter().map(|s| s.id as usize).max().unwrap_or(0);
    let mut out = vec![ParamMap::new(); n];
    for s in sections {
        let Some(i) = (s.id as usize).checked_sub(1) else {
            return Err(format!("{kind} section without a number"));
        };
        out[i] = section_to_param_map(s);
    }
    Ok(out)
}

fn single(section: &Option<Section>) -> ParamMap {
    section.as_ref().map(section_to_param_map).unwrap_or_default()
}

/// `elmer-solver`'s view of a parsed and completed model. Sections it has
/// no place for are refused.
pub fn to_model(m: &elmer_sif::Model) -> Result<Model, String> {
    if !m.boundaries.is_empty() {
        return Err("Boundary <n> sections".into());
    }
    if !m.components.is_empty() {
        return Err("Component sections".into());
    }
    if m.control.is_some() {
        return Err("Control section".into());
    }
    Ok(Model {
        header: ParamMap::new(),
        simulation: single(&m.simulation),
        constants: single(&m.constants),
        solvers: numbered("Solver", &m.solvers)?,
        equations: numbered("Equation", &m.equations)?,
        bodies: numbered("Body", &m.bodies)?,
        materials: numbered("Material", &m.materials)?,
        body_forces: numbered("Body Force", &m.body_forces)?,
        initial_conditions: numbered("Initial Condition", &m.ics)?,
        boundary_conditions: numbered("Boundary Condition", &m.bcs)?,
    })
}

fn is_empty(m: &elmer_sif::Model) -> bool {
    !m.header.present
        && m.control.is_none()
        && m.simulation.is_none()
        && m.constants.is_none()
        && m.boundaries.is_empty()
        && [&m.bodies, &m.bcs, &m.ics, &m.materials, &m.body_forces, &m.equations, &m.solvers, &m.components]
            .iter()
            .all(|v| v.is_empty())
}

impl Frontend for SifFrontend {
    fn load(&self, sif: &Path) -> Result<LoadedCase, String> {
        let path = sif.to_str().ok_or("non-UTF-8 case path")?;
        let mut parsed = elmer_sif::parse_file(path).map_err(|e| e.to_string())?;
        // A trailing `RUN` leaves an empty block after it; Elmer then stops at
        // end of file, so only blocks with content are further simulations.
        let runs = 1 + parsed.sif.runs.iter().skip(1).filter(|r| !is_empty(r)).count();
        if runs > 1 {
            return Err(format!("{runs} RUN blocks (one Simulation per process is composed)"));
        }
        elmer_sif::complete(&mut parsed.sif, path)?;
        let merged = parsed.sif.merged();
        let model = to_model(&merged)?;
        let case_dir = sif.parent().map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from("."));
        let (db, name) = merged.header.mesh_db.clone().ok_or("no Mesh DB in Header")?;
        let mesh_dir = case_dir.join(db).join(name);
        let mesh = (self.mesh)(&mesh_dir, &model.simulation)?;
        Ok(LoadedCase { model, mesh, case_dir, mesh_dir })
    }
}
