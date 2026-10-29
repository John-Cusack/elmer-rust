//! IterSolve keyword mapping (IterSolve.F90:724-1011).

use std::sync::Mutex;

use elmer_core::{LinalgError, Logger, ParamMap, Value};
use elmer_krylov::{HutiConfig, HutiMethod, HutiStop};

fn with_method(name: &str) -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Linear System Iterative Method", Value::Str(name.into()));
    p
}

#[derive(Default)]
struct Capture(Mutex<Vec<String>>);

impl Logger for Capture {
    fn info(&self, _: &str, msg: &str, level: u8) {
        self.0.lock().unwrap().push(format!("info{level}: {msg}"));
    }
    fn warn(&self, _: &str, msg: &str) {
        self.0.lock().unwrap().push(format!("warn: {msg}"));
    }
}

#[test]
fn method_names_case_insensitive() {
    for (name, m) in [
        ("CG", HutiMethod::Cg),
        ("cgs", HutiMethod::Cgs),
        ("BiCGStab", HutiMethod::BiCGStab),
        ("bicgstab2", HutiMethod::BiCGStab2),
        ("GMRES", HutiMethod::Gmres),
        (" tfqmr ", HutiMethod::Tfqmr),
    ] {
        assert_eq!(HutiConfig::from_params(&with_method(name)).unwrap().method, m);
    }
}

#[test]
fn unknown_names_take_the_bicgstab_default() {
    // CASE DEFAULT, IterSolve.F90:787; "qmr" is not an IterSolve name either.
    for name in ["qmr", "nonsense"] {
        let c = HutiConfig::from_params(&with_method(name)).unwrap();
        assert_eq!(c.method, HutiMethod::BiCGStab);
    }
}

#[test]
fn missing_method_warns_and_uses_bicgstab() {
    let log = Capture::default();
    let c = HutiConfig::from_params_with(&ParamMap::new(), false, &log).unwrap();
    assert_eq!(c.method, HutiMethod::BiCGStab);
    assert_eq!(
        log.0.lock().unwrap().as_slice(),
        ["warn: > Linear System Iterative Method < not found, using BiCGstab"]
    );
}

#[test]
fn itermethods_names_are_unsupported_until_l03() {
    for name in ["bicgstabl", "gcr", "idrs", "richardson", "jacobi", "sgs", "mprgp"] {
        let e = HutiConfig::from_params(&with_method(name)).unwrap_err();
        assert!(matches!(e, LinalgError::Unsupported(_)), "{name}");
    }
    // Complex BiCGStab(2) is remapped to BiCGStab(l) (IterSolve.F90:760-764).
    let e = HutiConfig::from_params_with(&with_method("bicgstab2"), true, &Capture::default())
        .unwrap_err();
    assert!(matches!(e, LinalgError::Unsupported(_)));
}

#[test]
fn defaults_when_only_the_method_is_given() {
    let c = HutiConfig::from_params(&with_method("gmres")).unwrap();
    assert_eq!(c.restart, 10);
    assert_eq!(c.verbosity, 1);
    assert_eq!(c.max_iter, 0, "ListGetInteger returns 0 when absent");
    assert_eq!(c.tol, 0.0, "ListGetCReal returns 0 when absent");
    assert_eq!(c.max_tol, 1.0e20);
    assert_eq!(c.stop_crit, HutiStop::TrueResidualScaledByB);
    assert!(c.abort_not_converged);
}

#[test]
fn values_are_read() {
    let mut p = with_method("gmres");
    p.insert("Linear System GMRES Restart", Value::Integer(vec![30]));
    p.insert("Linear System Max Iterations", Value::Integer(vec![700]));
    p.insert("Linear System Min Iterations", Value::Integer(vec![3]));
    p.insert("Linear System Convergence Tolerance", Value::real(1e-9));
    p.insert("Linear System Divergence Limit", Value::real(1e10));
    p.insert("Linear System Residual Output", Value::Integer(vec![25]));
    p.insert("Linear System Abort Not Converged", Value::Logical(false));
    let c = HutiConfig::from_params(&p).unwrap();
    assert_eq!(
        (c.restart, c.max_iter, c.min_iter, c.verbosity),
        (30, 700, 3, 25)
    );
    assert_eq!((c.tol, c.max_tol), (1e-9, 1e10));
    assert!(!c.abort_not_converged);
}

#[test]
fn preconditioning_side() {
    // BiCGStab(2), GMRES and TFQMR are always left-preconditioned.
    for name in ["bicgstab2", "gmres", "tfqmr"] {
        let mut p = with_method(name);
        p.insert("Linear System Left Preconditioning", Value::Logical(false));
        assert!(HutiConfig::from_params(&p).unwrap().left_precond, "{name}");
    }
    for name in ["cg", "cgs", "bicgstab"] {
        let mut p = with_method(name);
        assert!(!HutiConfig::from_params(&p).unwrap().left_precond, "{name}");
        p.insert("Linear System Left Preconditioning", Value::Logical(true));
        assert!(HutiConfig::from_params(&p).unwrap().left_precond, "{name}");
    }
}

#[test]
fn zero_guess_perturbation_only_for_the_bicgstab_family() {
    for (name, on) in [("bicgstab", true), ("bicgstab2", true), ("cg", false), ("gmres", false)] {
        assert_eq!(
            HutiConfig::from_params(&with_method(name)).unwrap().perturb_zero_guess,
            on,
            "{name}"
        );
    }
}

#[test]
fn invalid_values_are_refused() {
    let mut p = with_method("cg");
    p.insert("Linear System Max Iterations", Value::Integer(vec![0]));
    assert!(HutiConfig::from_params(&p).is_err(), "minv=1");
    let mut p = with_method("gmres");
    p.insert("Linear System GMRES Restart", Value::Integer(vec![0]));
    assert!(HutiConfig::from_params(&p).is_err());
}

#[test]
fn backward_error_criteria_are_unsupported() {
    for key in [
        "Linear System Componentwise Backward Error",
        "Linear System Normwise Backward Error",
    ] {
        let mut p = with_method("cg");
        p.insert(key, Value::Logical(true));
        assert!(matches!(
            HutiConfig::from_params(&p),
            Err(LinalgError::Unsupported(_))
        ));
    }
}
