//! `MinimalExample/case.sif` built by hand, standing in for `elmer-sif`.

use std::path::{Path, PathBuf};

use elmer_core::{ParamMap, Value};
use elmer_physics_heat::FeProvider;
use elmer_run::{Frontend, LoadedCase, SifFrontend};
use elmer_solver::Model;

use super::mesh::FileMesh;
use super::{basis::TestBasis, quad::TestQuadrature};

pub fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case)
}

/// A fresh copy of a case directory under the test scratch area.
pub fn scratch_case(case: &str, tag: &str) -> PathBuf {
    let dst = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("elmer-run-{tag}"));
    let _ = std::fs::remove_dir_all(&dst);
    copy_dir(&data(case), &dst);
    dst
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

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn ints(v: &[i32]) -> Value {
    Value::Integer(v.to_vec())
}

fn st(s: &str) -> Value {
    Value::Str(s.into())
}

/// The sections of `tests/data/minimal_example/case.sif`, with the string
/// values lowercased as Elmer's parser does outside quotes.
pub fn minimal_example() -> Model {
    Model {
        header: sec(&[("Check Keywords", st("warn")), ("Mesh DB", st("\".\" \"triangle\""))]),
        simulation: sec(&[("Simulation Type", st("steady state")), ("Post File", st("case.vtu"))]),
        solvers: vec![sec(&[
            ("Equation", st("heat equation")),
            ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
            ("Variable", st("temperature")),
            ("Linear System Solver", st("direct")),
            ("Reference Norm", Value::real(2.88675135e-01)),
        ])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Target Bodies", ints(&[1])), ("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(&[("Heat Conductivity", Value::real(1.0)), ("Density", Value::real(1.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[2])), ("Temperature", Value::real(0.0))]),
            sec(&[("Target Boundaries", ints(&[3])), ("Heat Flux", Value::real(1.0))]),
        ],
        ..Model::default()
    }
}

pub fn fe() -> FeProvider {
    || (Box::new(TestBasis), Box::new(TestQuadrature))
}

/// `elmer-sif` with the test-local mesh reader and basis.
pub fn sif_frontend() -> SifFrontend {
    SifFrontend { mesh: super::mesh::load, fe: fe() }
}

/// Serves the hand-built model for any `.sif`, with the mesh from
/// `<sif dir>/triangle`, after applying `edit`.
pub struct HandBuilt(pub fn(&mut Model));

impl Frontend for HandBuilt {
    fn fe(&self) -> FeProvider {
        fe()
    }

    fn load(&self, sif: &Path) -> Result<LoadedCase, String> {
        let case_dir = sif.parent().unwrap().to_path_buf();
        let mesh_dir = case_dir.join("triangle");
        let mut model = minimal_example();
        (self.0)(&mut model);
        Ok(LoadedCase { model, mesh: Box::new(FileMesh::read(&mesh_dir, 2)), case_dir, mesh_dir })
    }
}

