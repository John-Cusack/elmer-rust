//! Keyword handling as in IterSolve.F90.

use elmer_core::{ParamMap, Value};
use elmer_itermethods::*;

fn base() -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Linear System Max Iterations", Value::Integer(vec![500]));
    p.insert("Linear System Convergence Tolerance", Value::real(1e-8));
    p
}

#[test]
fn defaults_follow_itersolve() {
    let c = IterMethodsConfig::from_params(&base()).unwrap();
    assert_eq!(c.method_name, "bicgstab");
    assert_eq!(c.max_iterations, 500);
    assert_eq!(c.min_iterations, 0);
    assert_eq!(c.divergence_limit, 1.0e20);
    assert_eq!(c.output_interval, 1);
    assert_eq!(c.bicgstabl_degree, 2);
    assert_eq!(c.idrs_s, 4);
    assert_eq!(c.gcr_restart, 200);
    assert_eq!(c.sgs_omega, 1.8);
    assert!(c.robust.is_none());
    assert_eq!(c.pseudo_complex, PseudoComplex::Off);
    assert_eq!(c.random_seed, 314_159_265);
    // bicgstab is a HUTI method: not ours
    assert_eq!(c.method(false), None);
}

#[test]
fn gcr_restart_defaults_to_max_iterations_capped_at_200() {
    let mut p = base();
    p.insert("Linear System Max Iterations", Value::Integer(vec![80]));
    assert_eq!(IterMethodsConfig::from_params(&p).unwrap().gcr_restart, 80);
    p.insert("Linear System Max Iterations", Value::Integer(vec![3500]));
    assert_eq!(IterMethodsConfig::from_params(&p).unwrap().gcr_restart, 200);
    p.insert("Linear System GCR Restart", Value::Integer(vec![100]));
    assert_eq!(IterMethodsConfig::from_params(&p).unwrap().gcr_restart, 100);
}

#[test]
fn method_names_and_complex_bicgstab2_fallback() {
    let mut p = base();
    for (name, want) in [
        ("BiCGStabL", IterMethod::BiCgStabL),
        ("idrs", IterMethod::Idrs),
        ("GCR", IterMethod::Gcr),
        ("Jacobi", IterMethod::Jacobi),
        ("SGS", IterMethod::Sgs),
        ("Richardson", IterMethod::Richardson),
    ] {
        p.insert("Linear System Iterative Method", Value::Str(name.into()));
        let c = IterMethodsConfig::from_params(&p).unwrap();
        assert_eq!(c.method(false), Some(want), "{name}");
    }
    p.insert(
        "Linear System Iterative Method",
        Value::Str("BiCGStab2".into()),
    );
    let c = IterMethodsConfig::from_params(&p).unwrap();
    assert_eq!(c.method(false), None);
    assert_eq!(c.method(true), Some(IterMethod::BiCgStabL));
}

#[test]
fn optional_keywords_are_read() {
    let mut p = base();
    p.insert("BiCGstabl polynomial degree", Value::Integer(vec![4]));
    p.insert("IDRS parameter", Value::Integer(vec![2]));
    p.insert("IDRS Smoothing", Value::Logical(true));
    p.insert("SGS Overrelaxation Factor", Value::real(1.2));
    p.insert("Linear System Min Iterations", Value::Integer(vec![3]));
    p.insert("Linear System Divergence Limit", Value::real(1e6));
    p.insert("Linear System Residual Output", Value::Integer(vec![0]));
    p.insert("Linear System Pseudo Complex", Value::Logical(true));
    p.insert("Block Split Complex", Value::Logical(true));
    p.insert("Random Number Seed", Value::Integer(vec![42]));
    let c = IterMethodsConfig::from_params(&p).unwrap();
    assert_eq!((c.bicgstabl_degree, c.idrs_s, c.min_iterations), (4, 2, 3));
    assert!(c.idrs_smoothing);
    assert_eq!(c.sgs_omega, 1.2);
    assert_eq!(c.divergence_limit, 1e6);
    assert_eq!(c.output_interval, usize::MAX);
    assert_eq!(c.pseudo_complex, PseudoComplex::Split);
    assert_eq!(c.random_seed, 42);
}

#[test]
fn robust_defaults_derive_from_the_tolerance() {
    let mut p = base();
    p.insert("Linear System Robust", Value::Logical(true));
    let r = IterMethodsConfig::from_params(&p).unwrap().robust.unwrap();
    assert!((r.tolerance - 1e-8f64.powf(2.0 / 3.0)).abs() < 1e-20);
    assert!((r.max_tolerance - 1e-4).abs() < 1e-18);
    assert_eq!(r.step_size, 1.1);
    assert_eq!(r.max_bad_iterations, 250);
    assert_eq!(r.start_iteration, 1);
}

#[test]
fn missing_or_out_of_range_keywords_are_errors() {
    let mut p = ParamMap::new();
    p.insert("Linear System Convergence Tolerance", Value::real(1e-8));
    assert!(IterMethodsConfig::from_params(&p).is_err());
    let mut p = ParamMap::new();
    p.insert("Linear System Max Iterations", Value::Integer(vec![10]));
    assert!(IterMethodsConfig::from_params(&p).is_err());

    let mut p = base();
    p.insert("BiCGstabl polynomial degree", Value::Integer(vec![1]));
    assert!(IterMethodsConfig::from_params(&p).is_err());
    let mut p = base();
    p.insert("IDRS parameter", Value::Integer(vec![0]));
    assert!(IterMethodsConfig::from_params(&p).is_err());
    let mut p = base();
    p.insert("SGS Overrelaxation Factor", Value::real(2.5));
    assert!(IterMethodsConfig::from_params(&p).is_err());
}

#[test]
fn build_rejects_methods_from_other_crates() {
    let a = elmer_core::CsrMatrix::<f64> {
        nrows: 1,
        ncols: 1,
        row_ptr: vec![0, 1],
        col_idx: vec![0],
        values: vec![1.0],
        diag: Some(vec![0]),
    };
    let c = IterMethodsConfig {
        method_name: "cg".into(),
        ..Default::default()
    };
    assert!(c.build_real(&a).is_err());
    let c = IterMethodsConfig {
        method_name: "jacobi".into(),
        ..Default::default()
    };
    assert!(c.build_real(&a).is_ok());
    // Jacobi has no complex variant.
    assert!(c.build_complex().is_err());
}
