//! Options that are not ported are refused, never approximated; and a mapper
//! that has to move nodes says so when the runner installed no way to do it.

mod support;

use elmer_core::{NullLogger, SerialComm, Value};
use elmer_solver::{Model, SolverError};
use support::meshio::read_elmer_mesh;
use support::run::*;
use support::simmesh::SimMesh;

fn model(solver: &[(&str, Value)], body_force: &[(&str, Value)]) -> Model {
    Model {
        simulation: sec(&[
            ("Simulation Type", st("Steady")),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        solvers: vec![sec(solver)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Body Force", ints(&[1]))])],
        body_forces: vec![sec(body_force)],
        ..Model::default()
    }
}

fn mesh() -> SimMesh {
    SimMesh::new(read_elmer_mesh(&data("rigid1"), 2))
}

fn rigid<'a>(extra: &[(&'a str, Value)]) -> Vec<(&'a str, Value)> {
    let mut v = vec![
        ("Exec Solver", st("Before all")),
        ("Procedure", st("\"RigidMeshMapper\" \"RigidMeshMapper\"")),
    ];
    v.extend(extra.iter().map(|(k, x)| (*k, x.clone())));
    v
}

#[test]
fn rigid_mapper_refuses_rotor_mode_and_mesh_velocity() {
    for key in ["Rotor Mode", "Calculate Mesh Velocity", "mmg remesh"] {
        let m = model(
            &rigid(&[(key, lg(true))]),
            &[("Mesh Translate 1", Value::real(1.0))],
        );
        match run(m, &mesh()) {
            Err(SolverError::Unsupported(msg)) => assert!(msg.contains("RigidMeshMapper"), "{msg}"),
            other => panic!("{key}: {:?}", other.map(|_| ())),
        }
    }
}

#[test]
fn structured_mapper_refuses_masks() {
    let m = model(
        &[
            (
                "Procedure",
                st("\"StructuredMeshMapper\" \"StructuredMeshMapper\""),
            ),
            ("Active Coordinate", ints(&[2])),
            ("Mapping Mask Variable", st("TempB")),
        ],
        &[("Heat Source", Value::real(0.0))],
    );
    assert!(matches!(run(m, &mesh()), Err(SolverError::Unsupported(_))));
}

#[test]
fn moving_nodes_without_a_mesh_motion_is_not_wired() {
    let (mut sim, reg) = simulation(model(
        &rigid(&[]),
        &[("Mesh Translate 1", Value::real(1.0))],
    ));
    let mesh = mesh();
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    let r = sim.run(&mesh, &NullLogger, &SerialComm);
    assert!(
        matches!(r, Err(SolverError::NotWired(_))),
        "{:?}",
        r.map(|_| ())
    );
}

#[test]
fn a_mapper_that_moves_nothing_needs_no_mesh_motion() {
    let (mut sim, reg) = simulation(model(&rigid(&[]), &[("Heat Source", Value::real(1.0))]));
    let mesh = mesh();
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();
}
