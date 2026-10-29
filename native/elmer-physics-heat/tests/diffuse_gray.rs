//! `Radiation = Diffuse Gray` end to end: `fem/tests/radiation`, an
//! axisymmetric cavity whose inner walls radiate to each other, against the
//! `Reference Norm` of the pristine test (245.276478, tolerance 1e-4). The
//! view factors, their normalisation and the Gebhart factors are computed in
//! process; `elmer-viewfactors` compares them with the `ViewFactors` program.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{bicgstab_ilu0, data, heat_solver, ints, norm, run, sec, st};
use support::meshio::read_elmer_mesh;

fn model() -> Model {
    let mut solver = bicgstab_ilu0();
    solver.extend([
        ("Stabilize", Value::Logical(true)),
        ("Nonlinear System Newton After Iterations", ints(&[1])),
        ("Nonlinear System Newton After Tolerance", Value::real(1e-4)),
        ("Nonlinear System Max Iterations", ints(&[50])),
        ("Nonlinear System Consistent Norm", Value::Logical(true)),
        ("Nonlinear System Convergence Tolerance", Value::real(1e-8)),
        ("Nonlinear System Relaxation Factor", Value::real(0.7)),
    ]);
    let inner = sec(&[
        ("Target Boundaries", ints(&[1])),
        ("Heat Flux BC", Value::Logical(true)),
        ("Radiation", st("Diffuse Gray")),
        ("Radiation Target Body", ints(&[-1])),
        ("Emissivity", Value::real(0.6)),
    ]);
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Axi Symmetric")),
            ("Simulation Type", st("steady state")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ]),
        constants: sec(&[("Stefan Boltzmann", Value::real(5.67e-8))]),
        solvers: vec![heat_solver("Temperature", &solver)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![
            sec(&[
                ("Target Bodies", ints(&[1])),
                ("Equation", ints(&[1])),
                ("Body Force", ints(&[1])),
                ("Material", ints(&[1])),
                ("Initial Condition", ints(&[1])),
            ]),
            sec(&[
                ("Target Bodies", ints(&[2])),
                ("Equation", ints(&[1])),
                ("Material", ints(&[2])),
                ("Initial Condition", ints(&[1])),
            ]),
        ],
        initial_conditions: vec![sec(&[("Temperature", Value::real(250.0))])],
        body_forces: vec![sec(&[("Heat Source", Value::real(10000.0))])],
        materials: vec![
            sec(&[
                ("Density", Value::real(1.0)),
                ("Heat Conductivity", Value::real(10.0)),
                ("Heat Capacity", Value::real(1.0)),
            ]),
            sec(&[
                ("Density", Value::real(1.0)),
                ("Heat Conductivity", Value::real(1.0)),
                ("Heat Capacity", Value::real(1.0)),
            ]),
        ],
        boundary_conditions: vec![
            inner,
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Heat Flux BC", Value::Logical(true)),
                ("Radiation", st("Diffuse Gray")),
                ("Radiation Target Body", ints(&[-1])),
                ("Emissivity", Value::real(0.1)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[3])),
                ("Temperature", Value::real(100.0)),
            ]),
        ],
        ..Model::default()
    }
}

#[test]
fn axisymmetric_cavity_matches_elmer() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("radiation"), 2));
    let sim = run(model(), &mesh).unwrap();
    let got = norm(&sim, "temperature");
    assert!(
        (got - 245.276478).abs() <= 1e-4 * 245.276478,
        "norm {got} vs 245.276478"
    );
}

#[test]
fn cartesian_geometry_is_refused() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("radiation"), 2));
    let mut m = model();
    m.simulation.insert("Coordinate System", st("Cartesian"));
    let err = run(m, &mesh).err().expect("refused").to_string();
    assert!(err.contains("view"), "{err}");
}
