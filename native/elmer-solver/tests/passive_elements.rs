//! Passive elements (`Temperature Passive` in a body force), checked
//! against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/passive` from its `temp.sif` on its ElmerGrid mesh
//! (`tests/data/passive/`): body 2 is passive while the time table is
//! positive (t <= 5) and active after, BDF2 over six steps. Passive elements
//! add nothing to the system and their own nodes keep their values. Every
//! `ComputeChange` line and the reference norm. HeatSolve is the test-local
//! stand-in.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn passive_body_until_time_five() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/passive");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange", "CompareToReferenceSolution: Solver"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            // The second iterations change by the linear solver's tolerance.
            let tol = if wr < 1e-6 { 1e-8 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else {
            let cut = w.find("Norm =").unwrap() + 14;
            assert_eq!(g[..cut], w[..cut], "line {i}");
        }
    }
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
}

#[test]
fn without_passivity_the_first_step_differs() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/passive");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.body_forces[1].0.remove("temperature passive");
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    let first = log.lines(&["ComputeChange: NS"]).into_iter().next().unwrap();
    // The passive run's first step: 6.2087250.
    assert!(!first.contains("6.2087250"), "{first}");
}
