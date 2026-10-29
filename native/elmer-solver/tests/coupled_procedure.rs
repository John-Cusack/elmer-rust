//! `Coupled Solver = True` with a `Procedure` (CoupledSolver inheriting the
//! matrix), checked against the pristine ElmerSolver (26.2-devel, Rev
//! 9f6af2f85) run of `fem/tests/CoupledPoisson1b` from its `case.sif` on its
//! ElmerGrid mesh (`tests/data/angle`, `tests/data/coupled_poisson1b/`):
//! each coupled iteration runs HeatSolve with its linear solve disabled and
//! then solves the system itself; the CoupledSolver messages, `ComputeChange`
//! lines and the reference norm. HeatSolve is the test-local stand-in; the
//! case's SaveScalars solver never runs.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

/// The messages up to `Max Output Level = 10` (the default; the case comments its own out), as the pristine run.
struct Level10(Capture);

impl elmer_core::Logger for Level10 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 10 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

#[test]
fn coupled_solver_inherits_the_procedures_matrix() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("angle"), 2);
    let mut model = sif_model(&data.join("coupled_poisson1b/case.sif"));
    model.solvers.truncate(1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Level10(Capture::default());
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(data.join("coupled_poisson1b/stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.0.lines(&[
        "CoupledSolver",
        "DefaultSolve: Solver disabled",
        "ComputeChange",
        "CompareToReferenceSolution: Solver",
    ]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            let tol = if wr < 1e-6 { 1e-8 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else if w.starts_with("CompareTo") {
            let cut = w.find("Norm =").unwrap() + 14;
            assert_eq!(g[..cut], w[..cut], "line {i}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
}
