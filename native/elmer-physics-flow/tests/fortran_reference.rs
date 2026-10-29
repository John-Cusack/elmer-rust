//! `FlowSolver` against the pristine ElmerSolver on the cases in
//! `tests/data` (see `tests/data/README.md`): the final norm Elmer prints
//! (`FlowSolve: Result Norm`, velocity components only, `Nonlinear System
//! Norm Dofs = dim`) and every nodal velocity and pressure value of the last
//! time level of its `.result`.

mod support;

use elmer_core::{ParamMap, Value};
use elmer_solver::Simulation;
use support::cases::{data, flow_solver, ints, logical, model, nodal, norm, read_solution, real, run, sec, st, steady};
use support::mesh::TestMesh;

fn compare(case: &str, sim: &Simulation, mesh: &TestMesh, ref_norm: f64, tol: f64) {
    let got = nodal(sim, mesh.coords.len());
    let want = read_solution(&data(case).join("solution.dat"));
    assert_eq!(got.len(), want.len());
    let dofs = want[0].len();
    let mut err = vec![0.0f64; dofs];
    let mut scale = vec![0.0f64; dofs];
    for (g, w) in got.iter().zip(&want) {
        for c in 0..dofs {
            err[c] = err[c].max((g[c] - w[c]).abs());
            scale[c] = scale[c].max(w[c].abs());
        }
    }
    let rel_norm = (norm(sim) - ref_norm).abs() / ref_norm;
    let errs: Vec<String> = err.iter().map(|e| format!("{e:.1e}")).collect();
    eprintln!(
        "{case}: norm {:.17e} (ref {ref_norm:.17e}, rel {rel_norm:.1e}), max abs error per component [{}]",
        norm(sim),
        errs.join(", ")
    );
    assert!(rel_norm < tol, "{case}: norm {} vs {ref_norm}", norm(sim));
    for c in 0..dofs {
        assert!(err[c] <= tol * scale[c].max(1.0), "{case}: component {c} max error {}", err[c]);
    }
}

fn bc(tags: &[i32], entries: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![("Target Boundaries", ints(tags))];
    v.extend(entries.iter().cloned());
    sec(&v)
}

fn material(density: f64, viscosity: f64) -> ParamMap {
    sec(&[("Density", real(density)), ("Viscosity", real(viscosity))])
}

/// Stabilized Stokes (`Flow Model = Stokes`) on 12 x 4 Q1 elements: a
/// body-force-driven channel with no-slip walls and an `External Pressure`
/// outlet.
#[test]
fn stokes_channel_matches_elmer() {
    let mesh = TestMesh::read(&data("stokes_channel"), 2);
    let solver = flow_solver(
        2,
        &[
            ("Flow Model", st("Stokes")),
            ("Stabilize", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[1])),
            ("Nonlinear System Convergence Tolerance", real(1.0e-8)),
        ],
    );
    let bf = sec(&[("Flow BodyForce 1", real(1.0))]);
    let bcs = vec![
        bc(&[1, 3], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
        bc(&[2], &[("Velocity 2", real(0.0)), ("External Pressure", real(0.25))]),
        bc(&[4], &[("Velocity 2", real(0.0))]),
    ];
    let sim = run(model(steady(), solver, material(1.0, 0.5), Some(bf), bcs), &mesh).unwrap();
    compare("stokes_channel", &sim, &mesh, 0.12283444526335803, 1e-9);
}

/// Stabilized Navier-Stokes on 12 x 4 x 2 P1 triangles: plug inflow, slip
/// walls (`Slip Coefficient 1`), Picard then Newton after 3 iterations, run
/// to a 1e-10 relative change (8 iterations in Elmer).
#[test]
fn ns_channel_tri_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_channel_tri"), 2);
    let solver = flow_solver(
        2,
        &[
            ("Stabilize", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[30])),
            ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
            ("Nonlinear System Newton After Iterations", ints(&[3])),
            ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
        ],
    );
    let bcs = vec![
        bc(&[1, 3], &[("Velocity 2", real(0.0)), ("Slip Coefficient 1", real(10.0))]),
        bc(&[4], &[("Velocity 1", real(1.0)), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(model(steady(), solver, material(1.0, 0.02), None, bcs), &mesh).unwrap();
    compare("ns_channel_tri", &sim, &mesh, 0.705_377_133_197_478_1, 1e-9);
}

/// Q1/Q1 Navier-Stokes with Newton from the first iteration (no `Newton
/// After` keywords), `Div Discretization`, `Noslip Wall BC` walls, a
/// hydrostatic body force plus `Flow BodyForce 2`, and a `Pressure 1`
/// outlet traction (17 iterations in Elmer).
#[test]
fn ns_quad_newton_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_quad_newton"), 2);
    let solver = flow_solver(
        2,
        &[
            ("Stabilize", logical(true)),
            ("Div Discretization", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[20])),
            ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ],
    );
    let bf = sec(&[("Hydrostatic Pressure", logical(true)), ("Flow BodyForce 2", real(0.5))]);
    let bcs = vec![
        bc(&[1, 3], &[("Noslip Wall BC", logical(true))]),
        bc(&[2], &[("Pressure 1", real(-0.1))]),
        bc(&[4], &[("Velocity 1", real(1.0)), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(model(steady(), solver, material(2.0, 0.05), Some(bf), bcs), &mesh).unwrap();
    compare("ns_quad_newton", &sim, &mesh, 0.795_557_415_930_786_5, 1e-9);
}

/// `Flow Model = No Convection`, `Stabilization Method = Stabilized` (the
/// string spelling that turns on `Mesh % Stabilize`), `Gradp` and `Laplace
/// Discretization` on
/// P1 triangles: a `Normal Pressure` and `Mass Flux` inlet and a fixed
/// `Pressure` outlet.
#[test]
fn nc_tri_gradp_matches_elmer() {
    let mesh = TestMesh::read(&data("nc_tri_gradp"), 2);
    let solver = flow_solver(
        2,
        &[
            ("Flow Model", st("No Convection")),
            ("Stabilization Method", st("Stabilized")),
            ("Gradp Discretization", logical(true)),
            ("Laplace Discretization", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[1])),
        ],
    );
    let bcs = vec![
        bc(&[1, 3], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
        bc(&[4], &[("Normal Pressure", real(-1.0)), ("Mass Flux", real(0.1)), ("Velocity 2", real(0.0))]),
        bc(&[2], &[("Pressure", real(0.0)), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(model(steady(), solver, material(1.5, 1.0), None, bcs), &mesh).unwrap();
    compare("nc_tri_gradp", &sim, &mesh, 4.9176098168166885E-002, 1e-9);
}

/// Transient Navier-Stokes, BDF2, four steps of 0.1 from rest on P1
/// triangles: plug inflow, slip walls, Newton after 2 iterations each step.
#[test]
fn transient_tri_matches_elmer() {
    let mesh = TestMesh::read(&data("transient_tri"), 2);
    let simulation = sec(&[
        ("Coordinate System", st("Cartesian")),
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[2])),
        ("Timestep Sizes", real(0.1)),
        ("Timestep Intervals", ints(&[4])),
    ]);
    let solver = flow_solver(
        2,
        &[
            ("Stabilize", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[10])),
            ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
            ("Nonlinear System Newton After Iterations", ints(&[2])),
            ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
        ],
    );
    let bcs = vec![
        bc(&[1, 3], &[("Velocity 2", real(0.0)), ("Slip Coefficient 1", real(5.0))]),
        bc(&[4], &[("Velocity 1", real(1.0)), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(model(simulation, solver, material(1.0, 0.05), None, bcs), &mesh).unwrap();
    compare("transient_tri", &sim, &mesh, 0.694_674_261_311_076, 1e-9);
}

/// Q1/Q1 Navier-Stokes in 3D on 4 x 2 x 2 bricks: plug inflow, slip walls,
/// `External Pressure` outlet along the 3D face normal, Newton after 2
/// iterations.
#[test]
fn ns_hex3d_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_hex3d"), 3);
    let solver = flow_solver(
        3,
        &[
            ("Stabilize", logical(true)),
            ("Nonlinear System Max Iterations", ints(&[20])),
            ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
            ("Nonlinear System Newton After Iterations", ints(&[2])),
            ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
        ],
    );
    let bcs = vec![
        bc(
            &[1, 3],
            &[("Velocity 2", real(0.0)), ("Slip Coefficient 1", real(2.0)), ("Slip Coefficient 3", real(2.0))],
        ),
        bc(
            &[5, 6],
            &[("Velocity 3", real(0.0)), ("Slip Coefficient 1", real(2.0)), ("Slip Coefficient 2", real(2.0))],
        ),
        bc(&[4], &[("Velocity 1", real(1.0)), ("Velocity 2", real(0.0)), ("Velocity 3", real(0.0))]),
        bc(&[2], &[("External Pressure", real(0.2))]),
    ];
    let sim = run(model(steady(), solver, material(1.0, 0.1), None, bcs), &mesh).unwrap();
    compare("ns_hex3d", &sim, &mesh, 0.521_889_885_631_366_3, 1e-9);
}
