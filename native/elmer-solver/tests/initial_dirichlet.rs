//! `SetInitialConditions`' Dirichlet initialisation (ElmerSolver.F90:1862),
//! checked against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85):
//! HeatSolver, two BDF(1) steps of 0.1 on the 1dtests mesh, `T = 0`
//! initially, `T = 1` and `T = 3` at the ends; once as is
//! (`tests/data/initdir_stdout.txt`) and once with `Initialize Dirichlet
//! Conditions = False` (`initdir_off_stdout.txt`). The boundary values reach
//! the first step through the mass matrix, so the norms differ.

mod support;

use std::path::Path;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::poisson::heat;
use support::{i, r, s, sec, split_change, Capture};

fn model(init_dirichlet: Option<bool>) -> Model {
    let mut simulation = sec(&[
        ("Coordinate System", s("cartesian 1d")),
        ("Simulation Type", s("transient")),
        ("Timestepping Method", s("bdf")),
        ("BDF Order", i(&[1])),
        ("Timestep Intervals", i(&[2])),
        ("Timestep Sizes", r(0.1)),
        ("Steady State Max Iterations", i(&[1])),
    ]);
    if let Some(b) = init_dirichlet {
        simulation.insert("Initialize Dirichlet Conditions", Value::Logical(b));
    }
    Model {
        simulation,
        solvers: vec![sec(&[
            ("Equation", s("heat equation")),
            ("Variable", s("Temperature")),
            ("Procedure", s("\"HeatSolve\" \"HeatSolver\"")),
            ("Linear System Solver", s("direct")),
            ("Nonlinear System Max Iterations", i(&[1])),
        ])],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![sec(&[("Equation", i(&[1])), ("Initial Condition", i(&[1]))])],
        initial_conditions: vec![sec(&[("Temperature", r(0.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", i(&[1])), ("Temperature", r(1.0))]),
            sec(&[("Target Boundaries", i(&[2])), ("Temperature", r(3.0))]),
        ],
        ..Model::default()
    }
}

fn check(init_dirichlet: Option<bool>, stdout: &str) -> Simulation {
    let mesh = support::mesh::FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/1d"), 1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let mut sim = Simulation::new(model(init_dirichlet));
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let after_setup = sim.vars.get("temperature").unwrap().clone();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    let got = log.lines(&["ComputeChange"]);
    let want: Vec<&str> = stdout.lines().collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-7 * wr, "{g}\n{w}");
    }
    // Node 1 is x = 0 (BC 1), node 61 is x = pi (BC 2).
    let at = |n: u32| after_setup.values[after_setup.slot(n).unwrap()];
    let ends = if init_dirichlet == Some(false) { (0.0, 0.0) } else { (1.0, 3.0) };
    assert_eq!((at(0), at(60)), ends);
    assert_eq!(at(30), 0.0);
    sim
}

#[test]
fn boundary_values_start_the_transient() {
    check(None, include_str!("data/initdir_stdout.txt"));
    check(Some(true), include_str!("data/initdir_stdout.txt"));
}

#[test]
fn initialize_dirichlet_conditions_false_keeps_the_initial_condition() {
    check(Some(false), include_str!("data/initdir_off_stdout.txt"));
}
