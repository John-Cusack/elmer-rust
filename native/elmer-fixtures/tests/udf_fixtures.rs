mod support;

use std::sync::Arc;

use elmer_core::{Dependent, DependentKind, Value};
use elmer_expr::{compile_with, CompileOptions, UdfRegistry};
use elmer_fixtures::registry;
use elmer_fixtures::udf::{init_dist, my_source, OneDSource, OneDTolerances, INIT_DIST_RADIUS_SQ};
use elmer_fixtures::FixtureError;
use support::{NoMatc, TestModel};

fn procedure_value(vars: &[&str], lib: &str, entry: &str) -> Value {
    Value::Dependent(Dependent {
        variables: vars.iter().map(|s| s.to_string()).collect(),
        kind: DependentKind::Procedure { library: lib.into(), entry: entry.into() },
        rows: 1,
        cols: 1,
    })
}

#[test]
fn my_source_is_sin_of_time_plus_parameter() {
    assert_eq!(my_source([0.0, 0.0]), 0.0);
    assert_eq!(my_source([0.5, 1.0]), 1.5_f64.sin());
    assert_eq!(my_source([0.2, 2.01020304]), (0.2_f64 + 2.01020304).sin());
}

#[test]
fn my_source_through_expr_with_constant_parameter() {
    // ConstantParamFunc case.sif: `Heat Source = Variable "time, 4.01"`.
    let model = Arc::new(TestModel::line(&[0.0, 1.0, 2.0]));
    let reg = UdfRegistry::new();
    assert_eq!(registry::install("ConstantParamFunc", model.clone(), &reg), 1);
    let prop = compile_with(
        &procedure_value(&["time, 4.01"], "MySource", "MySource"),
        &reg,
        &CompileOptions::default(),
    )
    .unwrap();
    model.set_time(0.6);
    let mut out = [0.0; 3];
    prop.eval_nodes(&[0, 1, 2], &*model, &NoMatc, &mut out).unwrap();
    assert_eq!(out, [(0.6_f64 + 4.01).sin(); 3]);
}

#[test]
fn init_dist_radius_is_the_single_precision_square() {
    // Value read from the pristine InitDist.so (.rodata 0x2018).
    assert_eq!(INIT_DIST_RADIUS_SQ.to_bits(), 0x3FB7_0A3D_8000_0000);
    assert_ne!(INIT_DIST_RADIUS_SQ, 0.3_f64 * 0.3_f64);
}

#[test]
fn init_dist_inequality_is_strict_at_the_rim() {
    // Centre in, far out, and a point at distance exactly 0.3 in double
    // precision: inside, because the Fortran threshold is slightly larger.
    assert_eq!(init_dist(0.5, 0.0), 1.0);
    assert_eq!(init_dist(0.0, 0.0), 0.0);
    assert_eq!(init_dist(0.8, 0.0), 1.0);
    assert_eq!(init_dist(0.5, 0.3), 1.0);
    // Walk the doubles around the threshold along y: the last one below is
    // 1, the first one with y*y >= radius^2 is 0.
    let mut bits = INIT_DIST_RADIUS_SQ.sqrt().to_bits();
    while f64::from_bits(bits).powi(2) >= INIT_DIST_RADIUS_SQ {
        bits -= 1;
    }
    while f64::from_bits(bits + 1).powi(2) < INIT_DIST_RADIUS_SQ {
        bits += 1;
    }
    assert_eq!(init_dist(0.5, f64::from_bits(bits)), 1.0);
    assert_eq!(init_dist(0.5, f64::from_bits(bits + 1)), 0.0);
}

#[test]
fn init_dist_reads_node_coordinates_from_the_model() {
    let model = Arc::new(TestModel::points(&[(0.5, 0.0), (0.0, 0.0), (0.5, 0.29), (0.9, 0.0)]));
    let reg = UdfRegistry::new();
    assert_eq!(registry::install("AdvDiffFCT", model.clone(), &reg), 1);
    let prop = compile_with(
        &procedure_value(&["time"], "InitDist", "InitDist"),
        &reg,
        &CompileOptions::default(),
    )
    .unwrap();
    let mut out = [9.0; 4];
    prop.eval_nodes(&[0, 1, 2, 3], &*model, &NoMatc, &mut out).unwrap();
    assert_eq!(out, [1.0, 0.0, 1.0, 0.0]);
}

const XS: [f64; 4] = [0.0, 1.0, 2.0, 3.0];

fn one_d(model: &Arc<TestModel>, tol: OneDTolerances) -> OneDSource {
    OneDSource::new(model.clone(), tol)
}

fn fatal_text(r: Result<f64, FixtureError>) -> String {
    match r {
        Err(FixtureError::Fatal { caller, message }) => format!("{caller}: {message}"),
        other => panic!("expected Fatal, got {other:?}"),
    }
}

#[test]
fn source_walks_the_five_steps() {
    use std::f64::consts::PI;
    let m = Arc::new(TestModel::line(&XS));
    let src = one_d(&m, OneDTolerances::ONE_D_TESTS);
    assert_eq!(src.old_step(), -1);

    m.set_time(1.0);
    assert_eq!(src.eval(2, 7.0).unwrap(), 0.0);
    assert_eq!(src.old_step(), -1, "step 1 runs no check");

    m.set_time(2.0);
    m.set_potential(|_| 0.0);
    assert_eq!(src.eval(2, 7.0).unwrap(), 1.0);
    assert_eq!(src.old_step(), 2);

    m.set_time(3.0);
    m.set_potential(|x| -(x * x - PI * x) / 2.0);
    assert_eq!(src.eval(3, 7.0).unwrap(), 3.0, "step 3 returns x of the node");

    m.set_time(4.0);
    m.set_potential(|x| -(x * x * x - PI * PI * x) / 6.0);
    assert_eq!(src.eval(1, 7.5).unwrap(), 7.5, "step 4 returns the argument");

    m.set_time(5.0);
    m.set_potential(|x| 2.0 * x.sin() + 5.0e-5);
    assert_eq!(src.eval(1, -2.0).unwrap(), -2.0);
    assert_eq!(src.old_step(), 5);
}

#[test]
fn source_failed_checks_use_the_fortran_messages() {
    let m = Arc::new(TestModel::line(&XS));
    let src = one_d(&m, OneDTolerances::ONE_D_TESTS);
    type Field = fn(f64) -> f64;
    let cases: [(f64, &str, Field); 4] = [
        (2.0, "Not able to integrate a zero field ?", |_x: f64| 1.0e-300),
        (3.0, "Not able to integrate parabola ?", |x: f64| x),
        (4.0, "Not able to integrate 3d deg polynomial?", |x: f64| x),
        (5.0, "Not able to integrate sine?", |x: f64| 2.0 * x.sin() + 2.0e-4),
    ];
    for (t, text, field) in cases {
        m.set_time(t);
        m.set_potential(field);
        assert_eq!(fatal_text(src.eval(0, 0.0)), format!("Source: {text}"));
        assert_eq!(fatal_text(src.eval(0, 0.0)), format!("Source: {text}"), "OldStep stays put");
    }
}

#[test]
fn source_checks_each_step_once_per_instance() {
    let m = Arc::new(TestModel::line(&XS));
    let src = one_d(&m, OneDTolerances::ONE_D_TESTS);
    m.set_time(2.0);
    m.set_potential(|_| 0.0);
    src.eval(0, 0.0).unwrap();
    // A bad field after the check is not looked at again in the same step.
    m.set_potential(|_| 1.0);
    assert_eq!(src.eval(1, 0.0).unwrap(), 1.0);
    assert_eq!(src.old_step(), 2);
}

#[test]
fn source_instances_do_not_share_old_step() {
    let m = Arc::new(TestModel::line(&XS));
    let a = one_d(&m, OneDTolerances::ONE_D_TESTS);
    let b = one_d(&m, OneDTolerances::ONE_D_TESTS);
    m.set_time(2.0);
    m.set_potential(|_| 0.0);
    a.eval(0, 0.0).unwrap();
    assert_eq!((a.old_step(), b.old_step()), (2, -1));
    m.set_potential(|_| 1.0);
    a.eval(0, 0.0).unwrap(); // a already checked step 2
    assert!(b.eval(0, 0.0).is_err(), "b still has to check step 2");
}

#[test]
fn source_tolerance_differs_between_the_two_cases() {
    use std::f64::consts::PI;
    let m = Arc::new(TestModel::line(&XS));
    m.set_time(3.0);
    m.set_potential(|x| -(x * x - PI * x) / 2.0 + 1.0e-9);
    assert!(one_d(&m, OneDTolerances::ONE_D_TESTS).eval(0, 0.0).is_err());
    assert!(one_d(&m, OneDTolerances::MESH_LEVELS).eval(0, 0.0).is_ok());
}

#[test]
fn source_steps_outside_the_range_are_undefined() {
    let m = Arc::new(TestModel::line(&XS));
    let src = one_d(&m, OneDTolerances::ONE_D_TESTS);
    m.set_time(0.0);
    assert!(matches!(src.eval(0, 0.0), Err(FixtureError::Undefined { .. })));
    m.set_time(6.0);
    m.set_potential(|_| 0.0);
    assert!(matches!(src.eval(0, 0.0), Err(FixtureError::Undefined { .. })));
    // NINT: 4.6 is step 5.
    m.set_time(4.6);
    m.set_potential(|x| 2.0 * x.sin());
    assert_eq!(src.eval(0, 3.0).unwrap(), 3.0);
}

#[test]
fn source_missing_potential_is_an_error_not_a_skip() {
    let m = Arc::new(TestModel::line(&XS));
    let src = one_d(&m, OneDTolerances::ONE_D_TESTS);
    m.set_time(2.0);
    m.clear_potential_at(2);
    assert!(matches!(src.eval(0, 0.0), Err(FixtureError::Undefined { .. })));
}

#[test]
fn registered_source_uses_the_case_tolerance_and_is_a_fresh_instance() {
    use std::f64::consts::PI;
    let m = Arc::new(TestModel::line(&XS));
    let strict = UdfRegistry::new();
    let loose = UdfRegistry::new();
    assert_eq!(registry::install("1dtests", m.clone(), &strict), 1);
    assert_eq!(registry::install("1dtests_MeshLevels", m.clone(), &loose), 1);
    let prop = |reg: &UdfRegistry| {
        compile_with(&procedure_value(&["potential"], "Source", "Source"), reg, &CompileOptions::default())
            .unwrap()
    };
    m.set_time(3.0);
    m.set_potential(|x| -(x * x - PI * x) / 2.0 + 1.0e-9);
    let mut out = [0.0];
    assert!(prop(&loose).eval_nodes(&[1], &*m, &NoMatc, &mut out).is_ok());
    assert_eq!(out, [1.0], "step 3 gives x");
    let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut out = [0.0];
        let _ = prop(&strict).eval_nodes(&[1], &*m, &NoMatc, &mut out);
    }));
    assert!(caught.is_err(), "Fatal halts the strict case");
}

#[test]
fn unported_cases_and_plugins_are_not_found() {
    let m = Arc::new(TestModel::line(&XS));
    let reg = UdfRegistry::new();
    assert_eq!(registry::install("BDM2D", m.clone(), &reg), 0);
    assert!(reg.is_empty());
    assert!(registry::lookup("1dtests", "Poisson", "PoissonSolver", m.clone()).is_none());
    assert!(registry::lookup("1dtests", "Source", "source", m).is_some(), "entry names are case-insensitive");
    assert_eq!(registry::procedures("AdvDiffFCT"), vec![("InitDist", "InitDist")]);
}
