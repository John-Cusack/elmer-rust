//! `Legacy Assembly = True` end to end: `fem/tests/HeatControl` (smart
//! heater control) and `HeatControl2` (smart heater with idealized radiation
//! and a Newton switch), compared with the `Reference Norm` of the pristine
//! ElmerSolver's tests. The mesh is the tests' `square.grd` through
//! `ElmerGrid`; `data/heat_control/case*.sif` are the unchanged upstream
//! sifs, kept for the keyword values.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{data, direct, heat_solver, ints, norm, run, sec, st};
use support::meshio::read_elmer_mesh;

fn model(radiation: bool) -> Model {
    let mut solver = direct();
    solver.extend([
        ("Legacy Assembly", Value::Logical(true)),
        ("Steady State Convergence Tolerance", Value::real(1e-5)),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
    ]);
    if radiation {
        solver.extend([
            ("Nonlinear System Convergence Tolerance", Value::real(1e-6)),
            ("Nonlinear System Max Iterations", ints(&[30])),
            ("Nonlinear System Newton After Tolerance", Value::real(1e-3)),
            ("Nonlinear System Relaxation Factor", Value::real(0.7)),
        ]);
    } else {
        solver.extend([
            ("Nonlinear System Convergence Tolerance", Value::real(1e-4)),
            ("Nonlinear System Max Iterations", ints(&[1])),
            ("Nonlinear System Newton After Tolerance", Value::real(1e-5)),
            ("Nonlinear System Convergence Measure", st("solution")),
        ]);
    }
    let heater = Value::Real {
        rows: 3,
        cols: 1,
        data: vec![0.5, 0.5, 0.0],
    };
    let mut body_force = sec(&[
        ("Heat Source", Value::real(1.0)),
        ("Smart Heater Control", Value::Logical(true)),
        ("Smart Heater Control Point", heater),
        (
            "Smart Heater Temperature",
            Value::real(if radiation { 500.0 } else { 100.0 }),
        ),
    ]);
    let mut material = sec(&[
        ("Heat Conductivity", Value::real(1.0)),
        ("Density", Value::real(1.0)),
    ]);
    let mut body = sec(&[
        ("Target Bodies", ints(&[1])),
        ("Equation", ints(&[1])),
        ("Material", ints(&[1])),
        ("Body Force", ints(&[1])),
    ]);
    let mut bc = sec(&[("Target Boundaries", ints(&[1]))]);
    let mut model = Model::default();
    if radiation {
        material.insert("Emissivity", Value::real(0.7));
        body.insert("Initial Condition", ints(&[1]));
        bc.insert("Radiation", st("Idealized"));
        bc.insert("External Temperature", Value::real(300.0));
        model.initial_conditions = vec![sec(&[("Temperature", Value::real(300.0))])];
        model.constants = sec(&[("Stefan Boltzmann", Value::real(5.67e-8))]);
    } else {
        bc.insert("Temperature", Value::real(0.0));
    }
    body_force.insert("Name", st("Heater"));
    model.simulation = sec(&[
        ("Simulation Type", st("steady")),
        ("Steady State Max Iterations", ints(&[1])),
        ("Output Intervals", ints(&[0])),
    ]);
    model.solvers = vec![heat_solver("Temperature", &solver)];
    model.equations = vec![sec(&[("Active Solvers", ints(&[1]))])];
    model.bodies = vec![body];
    model.materials = vec![material];
    model.body_forces = vec![body_force];
    model.boundary_conditions = vec![bc];
    model
}

fn check(radiation: bool, want: f64) {
    let mesh = SimMesh::new(read_elmer_mesh(&data("heat_control"), 2));
    let sim = run(model(radiation), &mesh).unwrap();
    let got = norm(&sim, "temperature");
    assert!(
        (got - want).abs() <= 1e-5 * want.abs(),
        "norm {got} vs reference {want}"
    );
}

#[test]
fn heat_control_matches_elmer() {
    check(false, 53.345_329_537_558_88);
}

#[test]
fn heat_control_with_idealized_radiation_matches_elmer() {
    check(true, 428.791_041_504_157_5);
}
