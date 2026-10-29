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
/// empty, as an unset section is in Elmer. `elmer-sif` leaves an empty
/// id-0 section for each gap in the numbering (`Solver 1`, `Solver 3`);
/// those are skipped, and an id-0 section with content is refused.
fn numbered(kind: &str, sections: &[Section]) -> Result<Vec<ParamMap>, String> {
    let n = sections.iter().map(|s| s.id as usize).max().unwrap_or(0);
    let mut out = vec![ParamMap::new(); n];
    for s in sections {
        match (s.id as usize).checked_sub(1) {
            Some(i) => out[i] = section_to_param_map(s),
            None if s.entries.is_empty() => {}
            None => return Err(format!("{kind} section without a number")),
        }
    }
    Ok(out)
}

fn single(section: &Option<Section>) -> ParamMap {
    section.as_ref().map(section_to_param_map).unwrap_or_default()
}

/// `-distribute` (SetNormalizedKeywords, MainUtils.F90:224) divides the
/// value by the area or volume it acts on. The tag lives only on
/// `elmer-sif`'s entries, so it is refused here.
fn refuse_distribute(m: &elmer_sif::Model) -> Result<(), String> {
    let numbered: [(&str, &[Section]); 7] = [
        ("Solver", &m.solvers),
        ("Equation", &m.equations),
        ("Body", &m.bodies),
        ("Material", &m.materials),
        ("Body Force", &m.body_forces),
        ("Initial Condition", &m.ics),
        ("Boundary Condition", &m.bcs),
    ];
    let singles = [("Simulation", &m.simulation), ("Constants", &m.constants)];
    let all = numbered
        .iter()
        .flat_map(|(kind, v)| v.iter().map(move |s| (*kind, s)))
        .chain(singles.iter().filter_map(|(kind, s)| s.as_ref().map(|s| (*kind, s))));
    for (kind, s) in all {
        if let Some(e) = s.entries.iter().find(|e| e.distribute) {
            return Err(format!(
                "{kind} {}: {} = -distribute: keyword normalisation by area or volume is not implemented",
                s.id, e.name
            ));
        }
    }
    Ok(())
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
    refuse_distribute(m)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Params;

    #[test]
    fn numbering_gaps_are_empty_sections() {
        let dir = std::env::temp_dir().join(format!("elmer-run-sif-gaps-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sif = dir.join("case.sif");
        std::fs::write(
            &sif,
            "Header\n  Mesh DB \".\" \"m\"\nEnd\nSolver 1\n  Equation = \"a\"\nEnd\nSolver 3\n  Equation = \"b\"\nEnd\n",
        )
        .unwrap();
        let parsed = elmer_sif::parse_file(sif.to_str().unwrap()).unwrap();
        let model = to_model(&parsed.sif.merged()).unwrap();
        assert_eq!(model.solvers.len(), 3);
        assert!(model.solvers[1].0.is_empty());
        assert_eq!(model.solvers[2].string("equation"), Some("b"));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
