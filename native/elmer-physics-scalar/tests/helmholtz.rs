//! HelmholtzSolve end to end against the pristine ElmerSolver (`data/README.md`
//! records the runs): `fem/tests/HelmholtzFEM` (first pass, direct solve) and
//! `fem/tests/HelmholtzPlaneWaves`.

mod support;

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::Model;
use support::cases::*;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn helmholtz_solver(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("Helmholtz")),
        ("Variable", st("Pressure")),
        ("Variable Dofs", ints(&[2])),
        ("Procedure", st("\"HelmholtzSolve\" \"HelmholtzSolver\"")),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
        ("Nonlinear System Max Iterations", ints(&[1])),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

fn model(simulation: ParamMap, solver: ParamMap, bcs: Vec<ParamMap>) -> Model {
    Model {
        simulation,
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(&[
            ("Density", real(1.224)),
            ("Sound Speed", real(340.0)),
        ])],
        boundary_conditions: bcs,
        ..Model::default()
    }
}

fn check_pressure(sim: &elmer_solver::Simulation, mesh: &SimMesh, case: &str, tol: f64) {
    let n = mesh.mesh.coords.len();
    let cols = read_columns(&data(case).join("solution.dat"));
    for (k, want) in cols.iter().enumerate().take(2) {
        let d = rel_max_diff(&nodal(sim, "pressure", k, n), want);
        assert!(d < tol, "pressure {}: {d:e}", k + 1);
    }
}

#[test]
fn helmholtz_fem_circle_in_box_matches_elmer() {
    // Wave flux 1 on the circle (boundary 2), impedance 340 on the box
    // (boundary 1), 15000 rad/s.
    let mesh = SimMesh::new(read_elmer_mesh(&data("helmholtz_fem"), 2));
    let m = model(
        sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Angular Frequency", real(15000.0)),
        ]),
        helmholtz_solver(&[]),
        vec![
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Wave Flux 1", real(1.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("Wave Impedance 1", real(340.0)),
            ]),
        ],
    );
    let sim = run(m, &mesh).unwrap();
    check_pressure(&sim, &mesh, "helmholtz_fem", 1e-10);
    // The sif's Reference Norm is 0.128318191837E-01 (BiCGStab to 1e-8).
    assert_rel(norm(&sim, "pressure"), 0.01283181930, 1e-7);
    assert_rel(
        sim.model.simulation.real("res: frequency").unwrap(),
        15000.0 / (2.0 * std::f64::consts::PI),
        1e-14,
    );
}

#[test]
fn helmholtz_plane_waves_match_elmer() {
    // A plane wave from a pressure condition on the left, absorbed by the
    // impedance on the right; the frequency is in the solver section.
    let mesh = SimMesh::new(read_elmer_mesh(&data("helmholtz_plane_waves"), 2));
    let m = model(
        sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        helmholtz_solver(&[("Frequency", real(5000.0))]),
        vec![
            sec(&[("Target Boundaries", ints(&[1]))]),
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Wave Impedance 1", real(340.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[3]))]),
            sec(&[
                ("Target Boundaries", ints(&[4])),
                ("Pressure 1", real(1.0)),
                ("Pressure 2", real(0.0)),
            ]),
        ],
    );
    let sim = run(m, &mesh).unwrap();
    check_pressure(&sim, &mesh, "helmholtz_plane_waves", 1e-9);
    // min |p| over the nodes: the sif's Reference Norm 0.99967001.
    let n = mesh.mesh.coords.len();
    let (re, im) = (nodal(&sim, "pressure", 0, n), nodal(&sim, "pressure", 1, n));
    let min = re
        .iter()
        .zip(&im)
        .map(|(a, b)| a.hypot(*b))
        .fold(f64::INFINITY, f64::min);
    assert_rel(min, 0.9996700116143, 1e-8);
}
