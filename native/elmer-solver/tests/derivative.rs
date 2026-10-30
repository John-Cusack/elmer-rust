//! `Calculate Derivative` (the sensitivity `<var> Derivative`), checked
//! against a pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `tests/data/derivative/case.sif` plus a SaveScalars solver
//! (`scalars.dat`: sum and max of `Temperature Derivative`, sum of
//! `Temperature`): a scanning step with two coupled iterations, the heat
//! source `1 + 0.1 (i - 1)` and `Derivative Eps = 0.1 (i - 1)` in the
//! coupled iteration `i`, on the LinearSystemNamespace angle mesh
//! (`tests/data/namespace/`). HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Simulation, SolverError, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

fn run(steady: bool) -> (Result<(), SolverError>, Simulation, Capture) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("namespace"), 2);
    let mut model = sif_model(&data.join("derivative/case.sif"));
    if steady {
        model.simulation.insert("Simulation Type", Value::Str("Steady State".into()));
    }
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let r = sim.run(&mesh, &log, &SerialComm).map(|_| ());
    (r, sim, log)
}

#[test]
fn derivative_is_the_coupled_change_over_eps() {
    let (r, sim, log) = run(false);
    r.unwrap();
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/derivative");
    let stdout = std::fs::read_to_string(data.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange: NS", "ComputeChange: SS"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
    }

    let scalars: Vec<f64> = std::fs::read_to_string(data.join("scalars.dat"))
        .unwrap()
        .split_whitespace()
        .map(|x| x.parse().unwrap())
        .collect();
    let d = sim.vars.get("temperature derivative").unwrap();
    let t = sim.vars.get("temperature").unwrap();
    let got = [
        d.values.iter().sum::<f64>(),
        d.values.iter().copied().fold(f64::MIN, f64::max),
        t.values.iter().sum::<f64>(),
    ];
    for (g, w) in got.iter().zip(&scalars) {
        assert!((g - w).abs() <= 1e-9 * w.abs(), "{got:?} vs pristine {scalars:?}");
    }
}

#[test]
fn derivative_outside_scanning_is_fatal() {
    let (r, _, _) = run(true);
    match r {
        Err(SolverError::Fatal { msg, .. }) => assert_eq!(msg, "> Calculate Derivative < should be active only for scanning!"),
        other => panic!("expected the steady-state refusal, got {other:?}"),
    }
}
