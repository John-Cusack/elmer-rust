//! `Convergence Control` (divergence-driven substeps), checked against
//! the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/ConvergenceControl`: its real `case.sif` (read with elmer-sif,
//! whose MATC also evaluates `Heat Source = Variable Temperature, time; Real
//! MATC "1.0-alpha*tx(0)*tx(1)"` here) on the `angle` mesh, all 143
//! `ComputeChange` lines (`tests/data/convctl/stdout.txt`). The source grows
//! with time until the nonlinear iteration diverges (`Nonlinear System
//! Divergence Limit = 1e-3`) and the step is split. HeatSolve is the
//! test-local nonlinear diffusion stand-in.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn convergence_control_matches_elmer() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("angle"), 2);
    let model = sif_model(&dir.join("convctl/case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let keep = |l: &&str| l.contains("(NRM,RELC)") || l.starts_with("MAIN: Time") || l.starts_with("CompareTo");
    let got = log.lines(&["ComputeChange", "MAIN", "CompareTo"]);
    let got: Vec<&str> = got.iter().map(String::as_str).filter(keep).filter(|l| !l.contains("Relative Error")).collect();
    // `PASSED all 1 tests!` is `finalize_reference`'s.
    let want: Vec<&str> = include_str!("data/convctl/stdout.txt")
        .lines()
        .filter(|l| !l.contains("Relative Error") && !l.contains("all 1 tests"))
        .collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            assert!((gr - wr).abs() <= 1e-6 * wr, "line {i}\n{g}\n{w}");
        } else if w.starts_with("CompareTo") && w.contains("Norm =") {
            assert_eq!(g[..w.find("Norm =").unwrap()], w[..w.find("Norm =").unwrap()], "line {i}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "{got:#?}");
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
    // The late steps diverge and are split (Info level 6, not in Elmer's
    // level-5 stdout).
    let reduced = log.lines(&["ExecSimulation: Reducing timestep due to divergence problems!"]).len();
    assert!(reduced > 0);
}
