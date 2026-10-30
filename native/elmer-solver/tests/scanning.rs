//! `Simulation Type = Scanning`, checked against a pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `tests/data/scanning/case.sif`
//! (HeatSolve on the LinearSystemNamespace angle mesh, `tests/data/namespace/`,
//! with `Heat Source = Equals "time"` over three scanning steps): the steps
//! advance `time` by the default unit size, the solvers run steady and, with
//! one coupled iteration, no steady-state change is printed. HeatSolve is the
//! test-local stand-in.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn scanning_steps_time_and_solves_steady() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("namespace"), 2);
    let model = sif_model(&data.join("scanning/case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let report = sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!(report.steps, 3);
    assert_eq!(sim.vars.scalar("time"), Some(3.0));

    let stdout = std::fs::read_to_string(data.join("scanning/stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange", "MAIN: Time"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("MAIN") {
            assert_eq!(g, w);
            continue;
        }
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
    }
}
