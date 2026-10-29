//! `IncompressibleNsSolver` against the pristine ElmerSolver on the cases in
//! `tests/data` (see `tests/data/README.md`): every nodal velocity and
//! pressure value of the last time level of its `.result`, and the norm
//! Elmer prints (velocity components only, `Nonlinear System Norm Dofs =
//! dim`).

mod support;

use elmer_core::{ParamMap, Value};
use elmer_solver::Simulation;
use support::cases::{
    data, logical, model, nodal, norm, ns_solver, profile, read_solution, real, run, sec, steady, ints,
};
use support::mesh::TestMesh;

fn compare(case: &str, sim: &Simulation, mesh: &TestMesh, ref_norm: f64, norm_tol: f64, tol: f64) {
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
    let errs: Vec<String> = err.iter().zip(&scale).map(|(e, s)| format!("{:.1e}", e / s.max(1e-300))).collect();
    eprintln!(
        "{case}: norm {:.17e} (ref {ref_norm:.17e}, rel {rel_norm:.1e}), max relative error per component [{}]",
        norm(sim),
        errs.join(", ")
    );
    assert!(rel_norm < norm_tol, "{case}: norm {} vs {ref_norm}", norm(sim));
    for c in 0..dofs {
        assert!(err[c] <= tol * scale[c].max(1.0e-300), "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
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

fn step_bcs(mesh: &TestMesh) -> Vec<ParamMap> {
    vec![
        bc(
            &[1],
            &[("Velocity 1", profile(mesh, 1, 2, |y| (2.0 - y) * (y - 1.0))), ("Velocity 2", real(0.0))],
        ),
        bc(&[2], &[("Velocity 2", real(0.0))]),
        bc(&[3], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
    ]
}

fn stokes_step(case: &str, ref_norm: f64) {
    let mesh = TestMesh::read(&data(case), 2);
    let solver = ns_solver(&[
        ("Pressure Stabilization", logical(true)),
        ("Stokes Flow", logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-3)),
    ]);
    let sim = run(model(steady("Cartesian 2D"), solver, material(1.0, 1.0e3), None, step_bcs(&mesh)), &mesh).unwrap();
    compare(case, &sim, &mesh, ref_norm, 1e-7, 1e-9);
}

/// `StokesStabilized/tri.sif`: Stokes on the 1062 triangles of the backward
/// facing step, the Brezzi-Pitkaranta pressure block alone holding the
/// pressure mode.
#[test]
fn stokes_tri_matches_elmer() {
    stokes_step("stokes_tri", 7.58780146e-2);
}

/// `StokesStabilized/quad.sif`: the same on quadrilaterals.
#[test]
fn stokes_quad_matches_elmer() {
    stokes_step("stokes_quad", 7.55401834e-2);
}

/// `StokesStabilized/tet.sif`: a lid-driven unit cube of tetrahedra.
#[test]
fn stokes_tet_matches_elmer() {
    let mesh = TestMesh::read(&data("stokes_tet"), 3);
    let solver = ns_solver(&[
        ("Pressure Stabilization", logical(true)),
        ("Stokes Flow", logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-6)),
    ]);
    let bcs = vec![
        bc(
            &[5],
            &[
                ("Velocity 1", profile(&mesh, 5, 2, |y| y * (1.0 - y))),
                ("Velocity 2", real(0.0)),
                ("Velocity 3", real(0.0)),
            ],
        ),
        bc(&[6], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0)), ("Velocity 3", real(0.0))]),
    ];
    let sim = run(model(steady("Cartesian 3D"), solver, material(1.0, 1.0), None, bcs), &mesh).unwrap();
    compare("stokes_tet", &sim, &mesh, 5.66103337e-2, 1e-7, 1e-9);
}
