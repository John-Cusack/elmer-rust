//! Expression-valued properties, properties at integration points, damping and
//! first-order transient against `SaveScalars` of pristine ElmerSolver
//! (`tests/data/a2d/{expr,ipt,tdamp,tray}.sif` and their `*_scalars.dat`).
//! The mesh is `a2d`, plane strain triangles.

mod support;

use elmer_core::{Dependent, DependentKind, ParamMap, Value};
use elmer_solver::{Model, Simulation};
use elmer_solver::SolverError;
use support::cases::{a2d, run, try_run};

fn table(var: &str, x: &[f64], y: &[f64]) -> Value {
    Value::Dependent(Dependent {
        variables: vec![var.into()],
        kind: DependentKind::Table { x: x.to_vec(), y: y.to_vec(), cubic: false },
        rows: 1,
        cols: 1,
    })
}

fn reference(file: &str) -> Vec<Vec<f64>> {
    let p: std::path::PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "data", "a2d", file].iter().collect();
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect()
}

fn sum_and_max(xs: impl Iterator<Item = f64>) -> (f64, f64) {
    xs.fold((0.0, 0.0_f64), |(s, m), x| (s + x, if x.abs() > m.abs() { x } else { m }))
}

/// Displacement then (optionally) Stress_xx, Stress_yy, Stress_xy, vonMises, each
/// as `sum, max abs`: the column order of the `SaveScalars` solver in the sifs.
fn columns(sim: &Simulation, stresses: bool) -> Vec<f64> {
    let mut out = Vec::new();
    let u = sim.vars.get("displacement").unwrap();
    for c in 0..2 {
        let (s, m) = sum_and_max(u.values.iter().skip(c).step_by(2).copied());
        out.extend([s, m]);
    }
    if stresses {
        let st = sim.vars.get("stress").unwrap();
        for c in [0, 1, 3] {
            let (s, m) = sum_and_max(st.values.iter().skip(c).step_by(st.dofs).copied());
            out.extend([s, m]);
        }
        let (s, m) = sum_and_max(sim.vars.get("vonmises").unwrap().values.iter().copied());
        out.extend([s, m]);
    }
    out
}

fn assert_close(case: &str, got: &[f64], want: &[f64], tol: f64) {
    assert_eq!(got.len(), want.len(), "{case}");
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!((g - w).abs() <= tol * w.abs(), "{case} col {i}: got {g:e}, pristine {w:e}");
    }
}

fn graded(model: &mut Model, at_ip: bool) {
    // Youngs Modulus(x) and Poisson Ratio(y) as in expr.sif / ipt.sif.
    let m = &mut model.materials[0];
    m.insert("Youngs Modulus", table("coordinate 1", &[0.0, 0.1, 0.3], &[7.0e10, 1.0e11, 2.0e11]));
    m.insert("Poisson Ratio", table("coordinate 2", &[0.0, 0.05, 0.1], &[0.25, 0.3, 0.4]));
    if at_ip {
        m.insert("Youngs Modulus at IP", Value::Logical(true));
        m.insert("Poisson Ratio at IP", Value::Logical(true));
    }
    model.solvers[0].insert("Calculate Stresses", Value::Logical(true));
    model.solvers[0].insert("Calculate Strains", Value::Logical(false));
    model.solvers[0].insert("Calculate Principal", Value::Logical(false));
}

#[test]
fn table_valued_properties_match_pristine() {
    let (mut m, mesh) = a2d();
    graded(&mut m, false);
    let sim = run(m, &mesh);
    assert_close("expr", &columns(&sim, true), &reference("expr_scalars.dat")[0], 1e-8);
}

#[test]
fn properties_at_integration_points_match_pristine() {
    let (mut m, mesh) = a2d();
    graded(&mut m, true);
    let sim = run(m, &mesh);
    assert_close("ipt", &columns(&sim, true), &reference("ipt_scalars.dat")[0], 1e-8);
}

#[test]
fn ip_and_nodal_evaluation_differ_for_nonlinear_tables() {
    let (mut a, mesh) = a2d();
    graded(&mut a, false);
    let (mut b, _) = a2d();
    graded(&mut b, true);
    let (ca, cb) = (columns(&run(a, &mesh), false), columns(&run(b, &mesh), false));
    assert!((ca[0] - cb[0]).abs() > 1e-4 * ca[0].abs(), "{ca:?} {cb:?}");
}

fn transient(model: &mut Model, steps: i32, mat: &[(&str, Value)]) {
    let sim: ParamMap = support::sim::sec(&[
        ("Coordinate System", Value::Str("Cartesian 2D".into())),
        ("Simulation Type", Value::Str("Transient".into())),
        ("Timestepping Method", Value::Str("BDF".into())),
        ("BDF Order", Value::Integer(vec![1])),
        ("Timestep Intervals", Value::Integer(vec![steps])),
        ("Timestep Sizes", Value::real(0.05)),
        ("Output Intervals", Value::Integer(vec![1])),
    ]);
    model.simulation = sim;
    model.solvers[0].insert("Time Derivative Order", Value::Integer(vec![1]));
    model.solvers[0].insert("Calculate Stresses", Value::Logical(false));
    model.solvers[0].insert("Calculate Strains", Value::Logical(false));
    model.solvers[0].insert("Calculate Principal", Value::Logical(false));
    model.materials[0] = support::sim::sec(mat);
}

fn transient_matches(file: &str, mat: &[(&str, Value)]) {
    let want = reference(file);
    for (k, row) in want.iter().enumerate() {
        let (mut m, mesh) = a2d();
        transient(&mut m, k as i32 + 1, mat);
        let sim = run(m, &mesh);
        assert_close(&format!("{file} step {}", k + 1), &columns(&sim, false), row, 1e-8);
    }
}

#[test]
fn first_order_damping_matches_pristine() {
    transient_matches(
        "tdamp_scalars.dat",
        &[
            ("Density", Value::real(1.0)),
            ("Youngs Modulus", Value::real(70e9)),
            ("Poisson Ratio", Value::real(0.3)),
            ("Damping", Value::real(5.0e9)),
        ],
    );
}

#[test]
fn rayleigh_damping_matches_pristine() {
    transient_matches(
        "tray_scalars.dat",
        &[
            ("Density", Value::real(7800.0)),
            ("Youngs Modulus", Value::real(70e9)),
            ("Poisson Ratio", Value::real(0.3)),
            ("Rayleigh Damping alpha", Value::real(2.0e5)),
            ("Rayleigh Damping beta", Value::real(0.02)),
        ],
    );
}

#[test]
fn second_order_transient_is_refused_with_the_missing_hook_named() {
    let (mut m, mesh) = a2d();
    transient(&mut m, 1, &[("Youngs Modulus", Value::real(70e9)), ("Poisson Ratio", Value::real(0.3))]);
    m.solvers[0].0.remove("time derivative order");
    match try_run(m, &mesh) {
        Err(SolverError::NotWired(what)) => assert!(what.contains("Default2ndOrderTime"), "{what}"),
        other => panic!("expected NotWired, got {:?}", other.err()),
    }
}
