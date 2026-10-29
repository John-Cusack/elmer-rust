//! SaveData after the real `HeatSolver` in `Simulation::run`: a steady
//! two-cell case with the hooks scheduled as solvers 2 and 3 of the same
//! coupled loop. The steady temperature (0 / 6 / 10 across the three node
//! columns) is the pristine ElmerSolver's, `tests/data/ops2`.

mod support;

use std::path::PathBuf;

use elmer_core::SerialComm;
use elmer_solver::Model;
use support::backends::SimMesh;
use support::harness::*;
use support::meshio::read_elmer_mesh;

fn model(out: &std::path::Path) -> Model {
    let o = || ("Output Directory", dir_str(out));
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ]),
        solvers: vec![
            sec(&[
                ("Equation", st("heat equation")),
                ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
                ("Variable", st("Temperature")),
                ("Linear System Solver", st("Direct")),
            ]),
            sec(&[
                ("Equation", st("savescalars")),
                ("Procedure", st("\"SaveData\" \"SaveScalars\"")),
                ("Filename", st("scalars.dat")),
                ("Variable 1", st("Temperature")),
                ("Operator 1", st("max")),
                ("Operator 2", st("mean")),
                ("Operator 3", st("int")),
                ("Operator 4", st("body int")),
                ("Operator 5", st("boundary sum")),
                ("Operator 6", st("norm")),
                ("Mask Name 4", st("Save This")),
                o(),
            ]),
            sec(&[
                ("Equation", st("saveline")),
                ("Procedure", st("\"SaveData\" \"SaveLine\"")),
                ("Filename", st("line.dat")),
                ("Variable 1", st("Temperature")),
                o(),
            ]),
        ],
        equations: vec![sec(&[("Active Solvers", ints(&[1, 2, 3]))])],
        bodies: vec![
            sec(&[("Target Bodies", ints(&[1])), ("Equation", ints(&[1])), ("Material", ints(&[1]))]),
            sec(&[
                ("Target Bodies", ints(&[2])),
                ("Equation", ints(&[1])),
                ("Material", ints(&[1])),
                ("Body Force", ints(&[1])),
                ("Save This", yes()),
            ]),
        ],
        materials: vec![sec(&[("Heat Conductivity", real(1.0)), ("Density", real(1.0))])],
        body_forces: vec![sec(&[("Heat Source", real(4.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Heat Flux", real(0.0)), ("Save This", yes()), ("Save Line", yes())]),
            sec(&[("Target Boundaries", ints(&[3])), ("Heat Flux", real(0.0)), ("Save Line", yes())]),
            sec(&[("Target Boundaries", ints(&[4])), ("Temperature", real(0.0))]),
            sec(&[("Target Boundaries", ints(&[2])), ("Temperature", real(10.0))]),
        ],
        ..Model::default()
    }
}

#[test]
fn hooks_run_inside_the_coupled_loop_after_the_heat_solver() {
    let out = tmp("heat-steady");
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/two_cell");
    let mesh = SimMesh::new(read_elmer_mesh(&dir, 2));
    let log = RecLogger::default();
    let mut sim = heat_simulation(model(&out));
    sim.setup(&mesh, &full_registry(), &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    // Heat before the hooks: the solution they saw is the final one.
    let t = sim.vars.get("temperature").unwrap();
    for (n, want) in [0.0, 6.0, 10.0, 0.0, 6.0, 10.0].iter().enumerate() {
        assert!((t.values[t.slot(n as u32).unwrap()] - want).abs() < 1e-10, "node {n}");
    }

    let text = std::fs::read_to_string(out.join("scalars.dat")).unwrap();
    let r: Vec<f64> = text.split_whitespace().map(|t| t.parse().unwrap()).collect();
    assert_eq!(r.len(), 6);
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-9 * b.abs().max(1.0), "{a} vs {b}");
    near(r[0], 10.0);
    near(r[1], 16.0 / 3.0);
    // int T: cell 1 average (0 + 6) / 2 * 1, cell 2 (6 + 10) / 2 * 2.
    near(r[2], 19.0);
    // Body 2 only.
    near(r[3], 16.0);
    // Bottom nodes 1, 2, 3.
    near(r[4], 16.0);
    near(r[5], sim.vars.get("temperature").unwrap().conv.norm);

    let line = std::fs::read_to_string(out.join("line.dat")).unwrap();
    assert_eq!(line.lines().count(), 6);
    let temps: Vec<f64> = line
        .lines()
        .map(|l| l.split_whitespace().last().unwrap().parse().unwrap())
        .collect();
    for t in &temps {
        assert!([0.0, 6.0, 10.0].iter().any(|v| (t - v).abs() < 1e-9));
    }
    // Every file exists once, written during the single steady iteration.
    for f in ["scalars.dat.names", "line.dat.names"] {
        assert!(out.join(f).exists(), "{f}");
    }
}
