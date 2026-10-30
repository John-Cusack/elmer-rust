//! `IncompressibleNsSolver` against the pristine ElmerSolver on the cases in
//! `tests/data` (see `tests/data/README.md`): every nodal velocity and
//! pressure value of the last time level of its `.result`, and the norm
//! Elmer prints (velocity components only, `Nonlinear System Norm Dofs =
//! dim`).

mod support;

use elmer_core::{ParamMap, Value};
use elmer_solver::Simulation;
use support::cases::{
    data, logical, model, nodal, norm, ns_solver, profile, read_solution, real, run, sec, st, steady, ints,
};
use support::mesh::TestMesh;

fn compare(case: &str, sim: &Simulation, mesh: &TestMesh, ref_norm: f64, norm_tol: f64, tol: f64) {
    compare_with(case, sim, mesh, ref_norm, norm_tol, tol, false)
}

/// `floating_pressure`: the pressure level is a null space of the system
/// (the `GradP` form has no pressure in the traction), so only the pressure
/// up to a constant is compared.
fn compare_with(case: &str, sim: &Simulation, mesh: &TestMesh, ref_norm: f64, norm_tol: f64, tol: f64, floating_pressure: bool) {
    let mut got = nodal(sim, mesh.coords.len());
    let want = read_solution(&data(case).join("solution.dat"));
    if floating_pressure {
        let last = want[0].len() - 1;
        let shift = got.iter().zip(&want).map(|(g, w)| g[last] - w[last]).sum::<f64>() / got.len() as f64;
        for g in got.iter_mut() {
            g[last] -= shift;
        }
    }
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
    if std::env::var("DUMP").is_ok() {
        for (i, (g, w)) in got.iter().zip(&want).enumerate().step_by(37) {
            eprintln!("  node {} got {:?} want {:?} dp {:.6e}", i + 1, g, w, g[dofs - 1] - w[dofs - 1]);
        }
    }
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

/// `StokesStabilized/quad.sif` with `Slip Coefficient 1 = 200` on the
/// outflow and `Save Slip` (`tests/data/stokes_quad_saveslip`): the solution,
/// and the exported `Slip Coefficient` (normalized by the hidden `Slip
/// Weight`) and `Slip Speed` (zero without friction) at every node.
#[test]
fn save_slip_quad_matches_elmer() {
    let case = "stokes_quad_saveslip";
    let mesh = TestMesh::read(&data(case), 2);
    let solver = ns_solver(&[
        ("Pressure Stabilization", logical(true)),
        ("Save Slip", logical(true)),
        ("Stokes Flow", logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-3)),
    ]);
    let mut bcs = step_bcs(&mesh);
    bcs[1].insert("Slip Coefficient 1", real(200.0));
    let sim = run(model(steady("Cartesian 2D"), solver, material(1.0, 1.0e3), None, bcs), &mesh).unwrap();
    compare(case, &sim, &mesh, 0.07553691751285312, 1e-10, 1e-10);
    let rows: Vec<Vec<f64>> = std::fs::read_to_string(data(case).join("fields.dat"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(rows.len(), mesh.coords.len());
    for (c, name) in ["slip coefficient", "slip speed"].into_iter().enumerate() {
        let v = sim.vars.get(name).unwrap_or_else(|| panic!("no {name}"));
        for r in &rows {
            let node = r[0] as u32 - 1;
            let ours = v.values[v.slot(node).unwrap()];
            assert!((ours - r[1 + c]).abs() <= 1e-10 * 200.0, "{name} at node {}: {ours} vs {}", node + 1, r[1 + c]);
        }
    }
    assert!(sim.vars.get("slip weight").is_some_and(|w| !w.output), "-nooutput Slip Weight");
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

fn newton_options() -> Vec<(&'static str, Value)> {
    vec![
        ("Pressure Stabilization", logical(true)),
        ("Nonlinear System Max Iterations", ints(&[30])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
    ]
}

/// Navier-Stokes on the triangle step: Picard, then Newton after 3
/// iterations, `Pressure Stabilization Element Size = longest edge`; five
/// iterations in Elmer.
#[test]
fn ns_tri_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_tri"), 2);
    let mut opts = newton_options();
    opts.extend([
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
        ("Pressure Stabilization Element Size", st("longest edge")),
    ]);
    let sim = run(model(steady("Cartesian 2D"), ns_solver(&opts), material(1.0, 0.1), None, step_bcs(&mesh)), &mesh).unwrap();
    compare("ns_tri", &sim, &mesh, 0.75919202e-1, 1e-7, 1e-9);
}

/// Q1 Navier-Stokes with Newton after the first iteration, `GradP` and
/// `Div-Curl` discretizations, `Pressure Stabilization Element Size =
/// volume` with coefficient 1.5, a body force, a `Normal Surface Traction`
/// with a `Surface Traction 2` outlet and a slip wall; six iterations. The
/// pressure is compared up to its (undetermined) level.
#[test]
fn ns_quad_matches_elmer() {
    let mesh = TestMesh::read(&data("ns_quad"), 2);
    let mut opts = newton_options();
    opts.extend([
        ("Nonlinear System Newton After Iterations", ints(&[1])),
        ("GradP Discretization", logical(true)),
        ("Div-Curl Discretization", logical(true)),
        ("Pressure Stabilization Element Size", st("volume")),
        ("Pressure Stabilization Coefficient", real(1.5)),
    ]);
    let bf = sec(&[("Flow Bodyforce 2", real(-0.1))]);
    let bcs = vec![
        bc(
            &[1],
            &[("Velocity 1", profile(&mesh, 1, 2, |y| (2.0 - y) * (y - 1.0))), ("Velocity 2", real(0.0))],
        ),
        bc(
            &[2],
            &[("Normal Surface Traction", real(0.05)), ("Surface Traction 2", real(0.01)), ("Velocity 2", real(0.0))],
        ),
        bc(&[3], &[("Velocity 2", real(0.0)), ("Slip Coefficient 1", real(5.0))]),
    ];
    let sim = run(model(steady("Cartesian 2D"), ns_solver(&opts), material(2.0, 0.2), Some(bf), bcs), &mesh).unwrap();
    compare_with("ns_quad", &sim, &mesh, 0.14057138e-1, 1e-7, 1e-9, true);
}

/// Transient Navier-Stokes, BDF2, four steps of 0.1 from rest on triangles,
/// Newton after 2 iterations each step (the mass matrix carries the PSPG and
/// SUPG time terms).
#[test]
fn transient_tri_matches_elmer() {
    let mesh = TestMesh::read(&data("transient_tri"), 2);
    let simulation = sec(&[
        ("Coordinate System", st("Cartesian 2D")),
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[2])),
        ("Timestep Sizes", real(0.1)),
        ("Timestep Intervals", ints(&[4])),
        ("Steady State Max Iterations", ints(&[1])),
    ]);
    let solver = ns_solver(&[
        ("Pressure Stabilization", logical(true)),
        ("Nonlinear System Max Iterations", ints(&[3])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[2])),
    ]);
    let sim = run(model(simulation, solver, material(1.0, 0.1), None, step_bcs(&mesh)), &mesh).unwrap();
    compare("transient_tri", &sim, &mesh, 0.74655001e-1, 1e-7, 1e-8);
}
