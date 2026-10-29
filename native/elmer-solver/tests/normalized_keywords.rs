//! `<keyword> Normalize by Volume`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) runs of `fem/tests/NormalizeBF` (a body force
//! `Heat Source`) and `NormalizeMaterial` (a material `Heat Conductivity`)
//! from their `case.sif` files on their ElmerGrid meshes
//! (`tests/data/normalize/`): the normalisation messages, `ComputeChange`
//! and the reference norm. HeatSolve is the test-local stand-in; the
//! SaveScalars solver of the cases is left out.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

/// The messages up to `Max Output Level = 7`, as the pristine runs.
struct Level7(Capture);

impl elmer_core::Logger for Level7 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 7 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

fn check(case: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/normalize").join(case);
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.solvers.truncate(1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Level7(Capture::default());
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("{case}: {e}"));

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.0.lines(&[
        "ListTagKeywords",
        "ListTagCount",
        "CalculateEntityWeights: Computing",
        "ListSetParameters",
        "ComputeChange",
        "CompareToReferenceSolution: Solver",
    ]);
    assert_eq!(got.len(), want.len(), "{case}: {got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "{case} line {i}");
            assert!((gr - wr).abs() <= 1e-6 * wr, "{case} line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "{case} line {i}");
        }
    }
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "{case}");
}

#[test]
fn body_force_normalized_by_volume() {
    check("bf");
}

#[test]
fn material_normalized_by_volume() {
    check("material");
}
