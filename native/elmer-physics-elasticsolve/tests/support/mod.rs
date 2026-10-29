//! Test-local stand-ins for S1's `elmer-mesh` and `elmer-fe`, and the driver
//! that runs a pristine `.sif` through `Simulation` with `ElasticSolver`.
#![allow(dead_code)]

pub mod basis;
pub mod mesh;

use std::path::{Path, PathBuf};

use elmer_core::{MeshView, NullLogger, SerialComm};
use elmer_physics_elasticsolve::{install_fe, register, FeServices};
use elmer_solver::solver::procedure_of;
use elmer_solver::{finalize_reference, Simulation, SolverRegistry};

pub use basis::{Gauss, P1Basis};
use mesh::FileMesh;

pub fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case)
}

/// A case run to its end.
pub struct Run {
    pub sim: Simulation,
    /// `TEST.PASSED`: every `Reference Norm` held. `None` without references.
    pub passed: Option<bool>,
}

impl Run {
    /// `ComputeNorm` of the first solver's variable.
    pub fn norm(&self) -> f64 {
        let id = self.sim.slots[0].variable.expect("solver 1 has a variable");
        self.sim.vars.by_id(id).conv.norm
    }

    /// The solution of a variable in the order of its permutation slots.
    pub fn values(&self, name: &str) -> Vec<f64> {
        self.sim.vars.get(name).unwrap_or_else(|| panic!("no variable {name}")).values.clone()
    }
}

/// Runs `data/<case>/<sif>` (mesh in the directory its `Mesh DB` names) with
/// every solver other than `ElasticSolve` switched off, and `edit` applied to
/// the model first.
pub fn run_sif(case: &str, sif: &str, edit: impl FnOnce(&mut elmer_solver::Model)) -> Result<Run, elmer_solver::SolverError> {
    run_sif_with(case, sif, edit, |_| {})
}

/// [`run_sif`] with `after_setup` applied to the simulation once its variables exist.
pub fn run_sif_with(
    case: &str,
    sif: &str,
    edit: impl FnOnce(&mut elmer_solver::Model),
    after_setup: impl FnOnce(&mut Simulation),
) -> Result<Run, elmer_solver::SolverError> {
    let src = data(case);
    let scratch = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("elasticsolve-{case}-{}", sif.replace('.', "_")));
    let _ = std::fs::remove_dir_all(&scratch);
    copy_dir(&src, &scratch);
    let sif_path = scratch.join(sif);
    let path = sif_path.to_str().unwrap();
    let mut parsed = elmer_sif::parse_file(path).unwrap();
    elmer_sif::complete(&mut parsed.sif, path).unwrap();
    let merged = parsed.sif.merged();
    let mut model = elmer_run::sif::to_model(&merged).unwrap();
    let (db, name) = merged.header.mesh_db.clone().expect("Mesh DB");
    let mesh_dir = scratch.join(db).join(name);
    // `CoordinateSystemDimension()`: the mesh is 2D when it has no z extent.
    let mesh = FileMesh::read(&mesh_dir, 3);
    let flat = (0..mesh.n_nodes() as u32).all(|n| mesh.coord(n)[2] == 0.0);
    let mesh = if flat { FileMesh::read(&mesh_dir, 2) } else { mesh };

    let _ = install_fe(FeServices { basis: Box::new(P1Basis), quad: Box::new(Gauss) });
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    for p in &mut model.solvers {
        let known = procedure_of(p).is_some_and(|(l, e)| reg.get(&l, &e).is_some());
        if !known {
            p.insert("Exec Solver", elmer_core::Value::Str("never".into()));
        }
    }
    edit(&mut model);
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)?;
    after_setup(&mut sim);
    sim.run(&mesh, &NullLogger, &SerialComm)?;
    let passed = finalize_reference(&sim.reference, &scratch, &NullLogger, &SerialComm)?;
    Ok(Run { sim, passed })
}

fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for e in std::fs::read_dir(src).unwrap() {
        let e = e.unwrap();
        let to = dst.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &to);
        } else {
            std::fs::copy(e.path(), to).unwrap();
        }
    }
}
