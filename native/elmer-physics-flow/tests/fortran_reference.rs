//! `FlowSolver` against the pristine ElmerSolver on the cases in
//! `tests/data` (see `tests/data/README.md`): the final norm Elmer prints
//! (`ComputeChange ... NRM`, velocity components only, `Nonlinear System Norm
//! Dofs = 2`) and every nodal velocity and pressure value from its `.result`.

mod support;

use support::cases::{data, flow_solver, ints, nodal, norm, read_solution, real, run, sec, st, steady_model};
use support::mesh::TestMesh;

fn compare(case: &str, sim: &elmer_solver::Simulation, mesh: &TestMesh, ref_norm: f64, tol: f64) {
    let got = nodal(sim, mesh.coords.len());
    let want = read_solution(&data(case).join("solution.dat"));
    let mut err = [0.0f64; 3];
    let mut scale = [0.0f64; 3];
    for (g, w) in got.iter().zip(&want) {
        for c in 0..3 {
            err[c] = err[c].max((g[c] - w[c]).abs());
            scale[c] = scale[c].max(w[c].abs());
        }
    }
    let rel_norm = (norm(sim) - ref_norm).abs() / ref_norm;
    eprintln!("{case}: norm {:.17e} (ref {ref_norm:.17e}, rel {rel_norm:.2e}), max |du| {:.2e} |dv| {:.2e} |dp| {:.2e}", norm(sim), err[0], err[1], err[2]);
    assert!(rel_norm < tol, "{case}: norm {} vs {ref_norm}", norm(sim));
    for c in 0..3 {
        assert!(err[c] <= tol * scale[c].max(1.0), "{case}: component {c} max error {}", err[c]);
    }
}

/// Stabilized Stokes (`Flow Model = Stokes`) on 12 x 4 Q1 elements: a
/// body-force-driven channel with no-slip walls and an `External Pressure`
/// outlet.
#[test]
fn stokes_channel_matches_elmer() {
    let mesh = TestMesh::read(&data("stokes_channel"), 2);
    let solver = flow_solver(&[
        ("Flow Model", st("Stokes")),
        ("Stabilize", support::cases::logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-8)),
    ]);
    let material = sec(&[("Density", real(1.0)), ("Viscosity", real(0.5))]);
    let bf = sec(&[("Flow BodyForce 1", real(1.0))]);
    let bcs = vec![
        sec(&[("Target Boundaries", ints(&[1, 3])), ("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
        sec(&[("Target Boundaries", ints(&[2])), ("Velocity 2", real(0.0)), ("External Pressure", real(0.25))]),
        sec(&[("Target Boundaries", ints(&[4])), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(steady_model(solver, material, Some(bf), bcs), &mesh).unwrap();
    compare("stokes_channel", &sim, &mesh, 0.12283444526335803, 1e-9);
}

/// Stabilized Navier-Stokes on 12 x 4 x 2 P1 triangles: plug inflow, slip
/// walls (`Slip Coefficient 1`), Picard then Newton after 3 iterations, run
/// to a 1e-10 relative change (8 iterations in Elmer).
#[test]
fn ns_channel_tri_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_channel_tri"), 2);
    let solver = flow_solver(&[
        ("Stabilize", support::cases::logical(true)),
        ("Nonlinear System Max Iterations", ints(&[30])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
    ]);
    let material = sec(&[("Density", real(1.0)), ("Viscosity", real(0.02))]);
    let bcs = vec![
        sec(&[("Target Boundaries", ints(&[1, 3])), ("Velocity 2", real(0.0)), ("Slip Coefficient 1", real(10.0))]),
        sec(&[("Target Boundaries", ints(&[4])), ("Velocity 1", real(1.0)), ("Velocity 2", real(0.0))]),
    ];
    let sim = run(steady_model(solver, material, None, bcs), &mesh).unwrap();
    compare("ns_channel_tri", &sim, &mesh, 0.70537713319747808, 1e-9);
}
