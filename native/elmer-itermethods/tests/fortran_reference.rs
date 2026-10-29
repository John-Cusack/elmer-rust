//! Iteration counts and residual histories against the pristine ElmerSolver on a
//! dumped quick-test system (tests/fixtures/README.md).

mod common;

use common::*;
use elmer_core::{IdentityPrecond, IterativeSolver, SerialComm, SolveStatus};
use elmer_itermethods::*;

fn run(
    mut s: Box<dyn IterativeSolver<f64> + '_>,
    a: &elmer_core::CsrMatrix<f64>,
    b: &[f64],
    x0: &[f64],
) -> (Vec<f64>, elmer_core::SolveStats) {
    let mut x = x0.to_vec();
    let st = s
        .solve(a, &IdentityPrecond, b, &mut x, &SerialComm)
        .unwrap();
    (x, st)
}

fn cfg_tol(method: &str) -> IterMethodsConfig {
    IterMethodsConfig {
        max_iterations: 3500,
        tolerance: 1e-8,
        ..cfg(method)
    }
}

/// Compares histories to the 4 digits Elmer prints. Rounding differences grow in
/// the BiCG-type recurrences, so `tail` entries at the end get a looser bound.
fn assert_history(got: &[f64], want: &str, rel: f64) {
    let want = fortran_history(want);
    assert_eq!(got.len(), want.len(), "history length");
    for (k, (g, w)) in got.iter().zip(&want).enumerate() {
        assert!(
            (g - w).abs() <= rel * w.abs() + 1e-12,
            "iteration {}: {g:e} vs {w:e}",
            k + 1
        );
    }
}

#[test]
fn fixture_loads() {
    let (a, b, x0) = fortran_system();
    assert_eq!((a.nrows, b.len(), a.values.len()), (408, 408, 2690));
    assert_eq!(x0.iter().filter(|v| **v != 0.0).count(), 80);
}

#[test]
fn gcr_matches_fortran() {
    let (a, b, x0) = fortran_system();
    for (restart, count, file) in [
        (100, 53, include_str!("fixtures/ref_gcr.txt")),
        (10, 88, include_str!("fixtures/ref_gcr10.txt")),
    ] {
        let mut s = Gcr::<f64>::new(&IterMethodsConfig {
            gcr_restart: restart,
            ..cfg_tol("gcr")
        });
        let mut x = x0.clone();
        let st = s
            .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
            .unwrap();
        assert_eq!(st.status, SolveStatus::Converged);
        assert_eq!(st.iterations, count, "restart {restart}");
        assert_history(s.history(), file, 2e-3);
        assert!((x[0] - 1.0).abs() < 1e-6);
    }
}

#[test]
fn bicgstabl_matches_fortran() {
    let (a, b, x0) = fortran_system();
    for (l, count, file) in [
        (2, 20, include_str!("fixtures/ref_bicgstabl.txt")),
        (4, 9, include_str!("fixtures/ref_bicgstabl4.txt")),
    ] {
        let mut s = BiCgStabL::<f64>::new(&IterMethodsConfig {
            bicgstabl_degree: l,
            ..cfg_tol("bicgstabl")
        });
        let mut x = x0.clone();
        let st = s
            .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
            .unwrap();
        assert_eq!(st.status, SolveStatus::Converged);
        assert_eq!(st.iterations, count, "l = {l}");
        // The converging round exits early without a printed line.
        assert_history(s.history(), file, 2e-2);
        assert!((x[0] - 1.0).abs() < 1e-6);
    }
}

#[test]
fn jacobi_and_sgs_match_fortran() {
    let (a, b, x0) = fortran_system();
    let mut j = Jacobi::new(&cfg_tol("jacobi"), &a).unwrap();
    let mut x = x0.clone();
    let st = j
        .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    assert_eq!(st.iterations, 687);
    assert_eq!(st.status, SolveStatus::Converged);
    assert_history(j.history(), include_str!("fixtures/ref_jacobi.txt"), 2e-3);

    let mut s = Sgs::new(&cfg_tol("sgs"), &a);
    let mut x = x0.clone();
    let st = s
        .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    assert_eq!(st.iterations, 130);
    assert_history(s.history(), include_str!("fixtures/ref_sgs.txt"), 2e-3);
}

#[test]
fn idrs_count_is_near_fortran() {
    // RANDOM_NUMBER shadow vectors differ, so only the magnitude is comparable.
    let (a, b, x0) = fortran_system();
    for (s, fortran) in [(4, 59usize), (1, 66)] {
        let (x, st) = run(
            Box::new(Idrs::<f64>::new(&IterMethodsConfig {
                idrs_s: s,
                ..cfg_tol("idrs")
            })),
            &a,
            &b,
            &x0,
        );
        assert_eq!(st.status, SolveStatus::Converged);
        let (lo, hi) = (fortran * 7 / 10, fortran * 14 / 10);
        assert!(
            (lo..=hi).contains(&st.iterations),
            "s={s}: {} vs Fortran {fortran}",
            st.iterations
        );
        assert!((x[0] - 1.0).abs() < 1e-6);
    }
}

fn run_complex(
    mut s: Box<dyn IterativeSolver<num_complex::Complex64>>,
) -> (Vec<num_complex::Complex64>, elmer_core::SolveStats) {
    let (a, b) = fortran_helmholtz();
    let mut x = vec![num_complex::Complex64::new(0.0, 0.0); b.len()];
    let st = s
        .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    (x, st)
}

#[test]
fn complex_gcr_matches_fortran() {
    let (a, b) = fortran_helmholtz();
    assert_eq!((a.nrows, a.values.len()), (2110, 14258));
    for (restart, count, file) in [
        (100, 393, include_str!("fixtures/ref_zgcr.txt")),
        (10, 466, include_str!("fixtures/ref_zgcr10.txt")),
    ] {
        let mut s = Gcr::<num_complex::Complex64>::new(&IterMethodsConfig {
            gcr_restart: restart,
            ..cfg_tol("gcr")
        });
        let mut x = vec![num_complex::Complex64::new(0.0, 0.0); b.len()];
        let st = s
            .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
            .unwrap();
        assert_eq!(st.status, SolveStatus::Converged);
        assert_eq!(st.iterations, count, "restart {restart}");
        assert_history(s.history(), file, 2e-2);
        assert!(rel_resid(&a, &x, &b) < 2e-8);
    }
}

#[test]
fn complex_bicgstabl_matches_fortran() {
    let (a, b) = fortran_helmholtz();
    let mut s = BiCgStabL::<num_complex::Complex64>::new(&cfg_tol("bicgstabl"));
    let mut x = vec![num_complex::Complex64::new(0.0, 0.0); b.len()];
    let st = s
        .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    assert_eq!(st.status, SolveStatus::Converged);
    // Unpreconditioned Helmholtz makes BiCG-type recurrences chaotic: a 1e-13
    // relative change of the initial guess already moves round 22 by 1% and later
    // rounds completely, so rounding differences keep the count from matching
    // Fortran's 119 exactly. The first 20 rounds agree to the printed digits.
    assert!((100..=150).contains(&st.iterations), "{}", st.iterations);
    let want = fortran_history(include_str!("fixtures/ref_zbic.txt"));
    for (k, (g, w)) in s.history().iter().zip(&want).take(20).enumerate() {
        assert!((g - w).abs() <= 5e-3 * w, "round {}: {g:e} vs {w:e}", k + 1);
    }
    assert!(rel_resid(&a, &x, &b) < 2e-8);
}

#[test]
fn complex_idrs_count_is_near_fortran() {
    let (_x, st) = run_complex(Box::new(Idrs::<num_complex::Complex64>::new(&cfg_tol(
        "idrs",
    ))));
    assert_eq!(st.status, SolveStatus::Converged);
    // Fortran: 500 with its random shadow space.
    assert!((250..=1000).contains(&st.iterations), "{}", st.iterations);
}
