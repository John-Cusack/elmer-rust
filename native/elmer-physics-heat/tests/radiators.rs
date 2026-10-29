//! `Radiator BC` end to end: `fem/tests/radiator3d`, two point radiators
//! inside a cube whose walls lose the heat to 0 K by a heat transfer
//! coefficient, against the sif's `Reference Norm` (53.9315397, 1e-4). The
//! cube mesh is the one `elmer-viewfactors` compares the radiator factors on.

mod support;

use std::path::PathBuf;

use elmer_core::Value;
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{bicgstab_ilu0, heat_solver, ints, norm, run, sec, st};
use support::meshio::read_elmer_mesh;

fn model() -> Model {
    let mut solver = bicgstab_ilu0();
    solver.extend([
        ("Nonlinear System Newton After Tolerance", Value::real(1e-9)),
        ("Nonlinear System Newton After Iterations", ints(&[1])),
        ("Nonlinear System Max Iterations", ints(&[50])),
        ("Nonlinear System Convergence Tolerance", Value::real(1e-8)),
        ("Viewfactor Area Tolerance", Value::real(0.001)),
        ("Viewfactor Factor Tolerance", Value::real(0.001)),
        ("Viewfactor Number Of Rays", ints(&[8])),
        ("Calculate Loads", Value::Logical(true)),
    ]);
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("steady state")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ]),
        solvers: vec![heat_solver("Temperature", &solver)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Initial Condition", ints(&[1])),
        ])],
        initial_conditions: vec![sec(&[("Temperature", Value::real(250.0))])],
        materials: vec![sec(&[
            ("Density", Value::real(1.0)),
            ("Heat Conductivity", Value::real(1.0)),
        ])],
        body_forces: vec![sec(&[
            (
                "Radiator Coordinates",
                Value::Real {
                    rows: 2,
                    cols: 3,
                    data: vec![0.15, 0.15, 0.15, 0.85, 0.85, 0.85],
                },
            ),
            ("Radiator Power 1", Value::real(100.0)),
            ("Radiator Power 2", Value::real(200.0)),
        ])],
        boundary_conditions: vec![sec(&[
            ("Target Boundaries", ints(&[1, 2, 3])),
            ("Radiator BC", Value::Logical(true)),
            ("Radiation Target Body", ints(&[1])),
            ("External Temperature", Value::real(0.0)),
            ("Heat Transfer Coefficient", Value::real(1.0)),
        ])],
        ..Model::default()
    }
}

#[test]
fn cube_with_two_radiators_matches_elmer() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../elmer-viewfactors/tests/data/radiator3d");
    let mesh = SimMesh::new(read_elmer_mesh(&dir, 3));
    let sim = run(model(), &mesh).unwrap();
    let got = norm(&sim, "temperature");
    assert!(
        (got - 53.9315397).abs() <= 1e-4 * 53.9315397,
        "norm {got} vs 53.9315397"
    );
}
