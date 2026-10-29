//! `StructuredMeshMapper` on the pristine quick test `structmap`: the mapped
//! nodes against the coordinates of the pristine run's VTU file, then the
//! HeatSolve on the mapped mesh against the `.sif`'s `Reference Norm`.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::meshio::read_elmer_mesh;
use support::run::*;
use support::simmesh::SimMesh;

fn structmap_model(with_heat: bool) -> Model {
    let mut solvers = vec![];
    if with_heat {
        solvers.push(sec(&[
            ("Equation", st("Heat Equation")),
            ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
            ("Variable", st("Temperature")),
            ("Variable Dofs", ints(&[1])),
            ("Linear System Solver", st("Iterative")),
            ("Linear System Iterative Method", st("BiCGStab")),
            ("Linear System Max Iterations", ints(&[350])),
            ("Linear System Convergence Tolerance", Value::real(1.0e-12)),
            ("Linear System Abort Not Converged", lg(true)),
            ("Linear System Preconditioning", st("ILU0")),
            ("Linear System Residual Output", ints(&[1])),
            ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
            ("Stabilize", lg(true)),
            ("Nonlinear System Convergence Tolerance", Value::real(1.0e-5)),
            ("Nonlinear System Max Iterations", ints(&[1])),
            ("Nonlinear System Newton After Iterations", ints(&[3])),
            ("Nonlinear System Newton After Tolerance", Value::real(1.0e-2)),
            ("Nonlinear System Relaxation Factor", Value::real(1.0)),
            ("Reference Norm", Value::real(0.07352857410943443)),
        ]));
    }
    solvers.push(sec(&[
        ("Exec Solver", st("before all")),
        ("Equation", st("MapCoordinate")),
        ("Procedure", st("\"StructuredMeshMapper\" \"StructuredMeshMapper\"")),
        ("Active Coordinate", ints(&[2])),
    ]));
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[1])),
        ]),
        solvers,
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Body Force", ints(&[1])), ("Material", ints(&[1]))])],
        body_forces: vec![sec(&[("Heat Source", Value::real(1.0))])],
        materials: vec![sec(&[("Heat Conductivity", Value::real(1.0)), ("Density", Value::real(1.0))])],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("Bottom Surface", matc(&["coordinate 1"], "0.1*cos(5*tx)")),
                ("Temperature", Value::real(0.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[2])), ("Temperature", Value::real(0.0))]),
            sec(&[
                ("Target Boundaries", ints(&[3])),
                ("Top Surface", matc(&["coordinate 1"], "1+0.1*sin(5*tx)")),
                ("Temperature", Value::real(0.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[4])), ("Temperature", Value::real(0.0))]),
        ],
        ..Model::default()
    }
}

#[test]
fn structmap_nodes_match_pristine() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("structmap"), 2));
    let mut m = structmap_model(false);
    m.equations = vec![sec(&[("Active Solvers", ints(&[1]))])];
    run(m, &mesh).unwrap();
    let r = Reference::read("structmap.txt");
    for (k, c) in ["x", "y", "z"].into_iter().enumerate() {
        let d = max_diff(&coords_of(&mesh, k), &r.col(c));
        assert!(d < 1e-12, "{c}: max |diff| {d:e}");
    }
}

#[test]
fn structmap_heat_reaches_the_reference_norm() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("structmap"), 2));
    let sim = run(structmap_model(true), &mesh).unwrap();
    let r = Reference::read("structmap.txt");
    let d = max_diff(&nodal(&sim, "temperature", mesh.mesh.coords.len()), &r.col("temperature"));
    assert!(d < 1e-8, "temperature max |diff| {d:e}");
    let norm = sim.vars.get("temperature").unwrap().conv.norm;
    assert!(((norm - 0.07352857410943443) / norm).abs() < 1e-8, "{norm}");
}
