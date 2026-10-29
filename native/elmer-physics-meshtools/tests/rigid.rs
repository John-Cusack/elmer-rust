//! `RigidMeshMapper` on the pristine quick tests `RigidMeshMapper1` and
//! `RigidMeshMapper2`: the mapped nodes against the coordinates the pristine
//! run wrote to its VTU file, and the HeatSolve behind the mapper against the
//! `.sif`'s `Reference Norm`.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::meshio::read_elmer_mesh;
use support::run::*;
use support::simmesh::SimMesh;

fn mapper1_model(with_heat: bool) -> Model {
    let mut solvers = vec![sec(&[
        ("Exec Solver", st("Before all")),
        ("Equation", st("MeshDeform")),
        ("Procedure", st("\"RigidMeshMapper\" \"RigidMeshMapper\"")),
    ])];
    if with_heat {
        solvers.push(sec(&[
            ("Equation", st("HeatSolver")),
            ("Variable", st("Temperature")),
            ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
            ("Steady State Convergence Tolerance", Value::real(1e-5)),
            ("Nonlinear System Convergence Tolerance", Value::real(1e-4)),
            ("Nonlinear System Max Iterations", ints(&[1])),
            ("Linear System Solver", st("Iterative")),
            ("Linear System Iterative Method", st("BiCGStab")),
            ("Linear System Max Iterations", ints(&[500])),
            ("Linear System Convergence Tolerance", Value::real(1e-8)),
            ("Linear System Preconditioning", st("ILU1")),
            ("Linear System Abort Not Converged", lg(false)),
            ("Linear System Residual Output", ints(&[10])),
            ("Reference Norm", Value::real(3.657_515_699_823_906_5e-2)),
        ]));
    }
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[1])),
        ]),
        solvers,
        equations: vec![sec(&[("Active Solvers", ints(&[2]))])],
        bodies: vec![sec(&[
            ("Equation", ints(&[1])),
            ("Body Force", ints(&[1])),
            ("Material", ints(&[1])),
        ])],
        body_forces: vec![sec(&[
            (
                "Mesh Translate 2",
                matc(&["coordinate 1"], "0.2*sin(2*pi*tx)"),
            ),
            ("Mesh Relax", matc(&["coordinate 2"], "tx")),
            ("Heat Source", Value::real(1.0)),
        ])],
        materials: vec![sec(&[
            ("Heat Conductivity", Value::real(1.0)),
            ("Density", Value::real(1.0)),
        ])],
        boundary_conditions: vec![sec(&[
            ("Target Boundaries", ints(&[1])),
            ("Temperature", Value::real(0.0)),
        ])],
        ..Model::default()
    }
}

fn mapper2_model() -> Model {
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[1])),
        ]),
        solvers: vec![
            sec(&[
                ("Exec Solver", st("Before all")),
                ("Equation", st("MeshDeform")),
                ("Procedure", st("\"RigidMeshMapper\" \"RigidMeshMapper\"")),
                ("Cumulative Displacements", lg(true)),
                ("Variable", st("MeshRelax")),
                ("Nonlinear Conductivity Coefficient", Value::real(1.0)),
                ("Nonlinear System Relaxation Factor", Value::real(0.5)),
                ("Nonlinear System Relaxation After", ints(&[1])),
                ("Nonlinear System Max Iterations", ints(&[30])),
                (
                    "Nonlinear System Convergence Tolerance",
                    Value::real(1.4e-4),
                ),
                ("Linear System Solver", st("Iterative")),
                ("Linear System Iterative Method", st("BiCGStab")),
                ("Linear System Max Iterations", ints(&[500])),
                ("Linear System Convergence Tolerance", Value::real(1e-12)),
                ("Linear System Preconditioning", st("ILU1")),
                ("Linear System Abort Not Converged", lg(false)),
                ("Linear System Residual Output", ints(&[20])),
                ("Optimize Bandwidth", lg(false)),
            ]),
            sec(&[
                ("Equation", st("HeatSolver")),
                ("Variable", st("Temperature")),
                ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
                ("Steady State Convergence Tolerance", Value::real(1e-5)),
                ("Nonlinear System Convergence Tolerance", Value::real(1e-4)),
                ("Nonlinear System Max Iterations", ints(&[1])),
                ("Linear System Solver", st("Iterative")),
                ("Linear System Iterative Method", st("BiCGStab")),
                ("Linear System Max Iterations", ints(&[500])),
                ("Linear System Convergence Tolerance", Value::real(1e-8)),
                ("Linear System Preconditioning", st("ILU1")),
                ("Linear System Abort Not Converged", lg(false)),
                ("Linear System Residual Output", ints(&[10])),
                ("Reference Norm", Value::real(4.41553009e-2)),
            ]),
        ],
        equations: vec![
            sec(&[("Active Solvers", ints(&[1, 2]))]),
            sec(&[("Active Solvers", ints(&[2]))]),
        ],
        bodies: vec![
            sec(&[
                ("Equation", ints(&[1])),
                ("Body Force", ints(&[1])),
                ("Material", ints(&[1])),
            ]),
            sec(&[
                ("Equation", ints(&[2])),
                ("Body Force", ints(&[2])),
                ("Material", ints(&[2])),
            ]),
        ],
        body_forces: vec![
            sec(&[
                ("Mesh Rotate 3", Value::real(30.0)),
                ("Mesh Origin", reals(&[1.5, 1.25, 0.0])),
                ("Mesh Translate 2", Value::real(0.1)),
            ]),
            sec(&[
                ("Heat Source", Value::real(1.0)),
                ("Mesh Rotate 3", Value::real(30.0)),
                ("Mesh Origin", reals(&[1.5, 1.25, 0.0])),
                ("Mesh Translate 2", Value::real(0.1)),
            ]),
        ],
        materials: vec![
            sec(&[
                ("Heat Conductivity", Value::real(1.0)),
                ("Density", Value::real(1.0)),
            ]),
            sec(&[
                ("Heat Conductivity", Value::real(100.0)),
                ("Density", Value::real(1.0)),
            ]),
        ],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("Fixed Boundary", lg(true)),
                ("Temperature", Value::real(0.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Moving Boundary", lg(true)),
            ]),
        ],
        ..Model::default()
    }
}

#[test]
fn rigid1_nodes_match_pristine() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("rigid1"), 2));
    let mut m = mapper1_model(false);
    m.equations = vec![sec(&[("Active Solvers", ints(&[1]))])];
    run(m, &mesh).unwrap();
    let r = Reference::read("rigid1.txt");
    for (k, c) in ["x", "y", "z"].into_iter().enumerate() {
        let d = max_diff(&coords_of(&mesh, k), &r.col(c));
        assert!(d < 1e-12, "{c}: max |diff| {d:e}");
    }
}

#[test]
fn rigid1_heat_reaches_the_reference_norm() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("rigid1"), 2));
    let sim = run(mapper1_model(true), &mesh).unwrap();
    let r = Reference::read("rigid1.txt");
    let d = max_diff(
        &nodal(&sim, "temperature", mesh.mesh.coords.len()),
        &r.col("temperature"),
    );
    assert!(d < 1e-7, "temperature max |diff| {d:e}");
    let norm = sim.vars.get("temperature").unwrap().conv.norm;
    assert!(
        ((norm - 3.657_515_699_823_906_5e-2) / norm).abs() < 1e-6,
        "{norm}"
    );
}

#[test]
fn rigid2_relaxed_rotation_matches_pristine() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("rigid2"), 2));
    let n = mesh.mesh.coords.len();
    let sim = run(mapper2_model(), &mesh).unwrap();
    let r = Reference::read("rigid2.txt");
    let relax = max_diff(&nodal(&sim, "meshrelax", n), &r.col("meshrelax"));
    assert!(relax < 1e-6, "meshrelax max |diff| {relax:e}");
    for (k, c) in ["x", "y", "z"].into_iter().enumerate() {
        let d = max_diff(&coords_of(&mesh, k), &r.col(c));
        assert!(d < 1e-6, "{c}: max |diff| {d:e}");
    }
    let norm = sim.vars.get("temperature").unwrap().conv.norm;
    assert!(((norm - 4.41553009e-2) / norm).abs() < 1e-6, "{norm}");
}
