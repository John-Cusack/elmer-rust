//! Behaviour of each method on small hand-checkable systems.

mod common;

use common::*;
use elmer_core::{IdentityPrecond, IterativeSolver, SerialComm, SolveStatus};
use elmer_itermethods::*;
use num_complex::Complex64;

const N: usize = 40;

fn solve_real(
    solver: &mut dyn IterativeSolver<f64>,
    a: &elmer_core::CsrMatrix<f64>,
    b: &[f64],
) -> (Vec<f64>, elmer_core::SolveStats) {
    let mut x = vec![0.0; b.len()];
    let st = solver
        .solve(a, &IdentityPrecond, b, &mut x, &SerialComm)
        .unwrap();
    (x, st)
}

#[test]
fn jacobi_matches_hand_iteration() {
    // 2x2: A = [[4,1],[2,5]], b = [1,2]. x0 = 0 -> x1 = b/d = [0.25, 0.4].
    let a = tridiag(2, 2.0, 4.0, 1.0);
    let mut a = a;
    a.values = vec![4.0, 1.0, 2.0, 5.0];
    let b = [1.0, 2.0];
    let mut cfg = cfg("jacobi");
    cfg.max_iterations = 1;
    cfg.tolerance = 1e-30;
    let mut j = Jacobi::new(&cfg, &a).unwrap();
    let mut x = vec![0.0; 2];
    let st = j.solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
    assert_eq!(x, vec![0.25, 0.4]);
    assert_eq!(st.iterations, 1);
    assert_eq!(st.status, SolveStatus::MaxIterations);
    // relative residual: r = b - A x = [1-1.4, 2-(0.5+2)] = [-0.4, -0.5]
    let want = (0.16f64 + 0.25).sqrt() / 5.0f64.sqrt();
    assert!((st.residual - want).abs() < 1e-15);
}

#[test]
fn stationary_methods_converge() {
    let a = tridiag(N, -1.0, 4.0, -1.0);
    let b = rhs_real(N);
    let c = cfg("jacobi");
    let (x, st) = solve_real(&mut Jacobi::new(&c, &a).unwrap(), &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-9);

    let c = IterMethodsConfig {
        sgs_omega: 1.0,
        ..cfg("sgs")
    };
    let (x, st) = solve_real(&mut Sgs::new(&c, &a), &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-9);
    // SGS is far faster than Jacobi on this system.
    let (_, sj) = solve_real(&mut Jacobi::new(&cfg("jacobi"), &a).unwrap(), &a, &b);
    assert!(st.iterations < sj.iterations);
}

#[test]
fn richardson_first_sweep_discards_guess() {
    // Mass-like matrix: lumping by row sums is a good preconditioner.
    let a = tridiag(N, 1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0);
    let b = rhs_real(N);
    let mut r = Richardson::new(&cfg("richardson"), &a);
    let mut x = vec![123.0; N];
    let st = r.solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-9);

    // One iteration leaves exactly x = b / rowsum, whatever the guess was.
    let mut c = cfg("richardson");
    c.max_iterations = 1;
    let mut r = Richardson::new(&c, &a);
    let mut x = vec![123.0; N];
    r.solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
    assert!((x[5] - b[5] / 1.0).abs() < 1e-14);
    assert!((x[0] - b[0] / (5.0 / 6.0)).abs() < 1e-14);
}

#[test]
fn zero_rhs_is_divergence_for_stationary_and_converged_for_bicgstabl() {
    // Jacobi: 0/0 residual is NaN, which Elmer treats as diverged.
    let a = tridiag(N, -1.0, 4.0, -1.0);
    let b = vec![0.0; N];
    let (_, st) = solve_real(&mut Jacobi::new(&cfg("jacobi"), &a).unwrap(), &a, &b);
    assert_eq!(st.status, SolveStatus::Diverged);
    // BiCGStab(l) special-cases b = 0 and returns x = 0.
    let mut x = vec![1.0; N];
    let st = BiCgStabL::<f64>::new(&cfg("bicgstabl"))
        .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(x.iter().all(|v| *v == 0.0));
}

#[test]
fn bicgstabl_solves_nonsymmetric_for_several_degrees() {
    let a = nonsym(N);
    let b = rhs_real(N);
    for l in [2, 3, 4] {
        let c = IterMethodsConfig {
            bicgstabl_degree: l,
            ..cfg("bicgstabl")
        };
        let (x, st) = solve_real(&mut BiCgStabL::<f64>::new(&c), &a, &b);
        assert_eq!(st.status, SolveStatus::Converged, "l={l}");
        assert!(rel_resid(&a, &x, &b) < 1e-9, "l={l}");
    }
}

#[test]
fn bicgstabl_degree_below_two_is_rejected() {
    let a = nonsym(4);
    let c = IterMethodsConfig {
        bicgstabl_degree: 1,
        ..cfg("bicgstabl")
    };
    let mut x = vec![0.0; 4];
    let r = BiCgStabL::<f64>::new(&c).solve(&a, &IdentityPrecond, &[1.0; 4], &mut x, &SerialComm);
    assert!(r.is_err());
}

#[test]
fn bicgstabl_counts_rounds_like_elmer() {
    // A diagonal system converges in the first round for any l.
    let a = tridiag(6, 0.0, 2.0, 0.0);
    let b = vec![2.0; 6];
    let (x, st) = solve_real(&mut BiCgStabL::<f64>::new(&cfg("bicgstabl")), &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(x.iter().all(|v| (v - 1.0).abs() < 1e-8));
    assert!(st.iterations <= 1, "{}", st.iterations);
}

#[test]
fn bicgstabl_honours_max_iterations() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let c = IterMethodsConfig {
        max_iterations: 2,
        tolerance: 1e-14,
        ..cfg("bicgstabl")
    };
    let (_, st) = solve_real(&mut BiCgStabL::<f64>::new(&c), &a, &b);
    assert_eq!(st.status, SolveStatus::MaxIterations);
    assert_eq!(st.iterations, 2);
}

#[test]
fn gcr_solves_and_restart_changes_the_history() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let full = IterMethodsConfig {
        gcr_restart: 500,
        ..cfg("gcr")
    };
    let mut g = Gcr::<f64>::new(&full);
    let (x, st) = solve_real(&mut g, &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-9);
    // GCR is optimal: with an unrestarted space the residual decreases monotonically.
    let h = g.history();
    assert_eq!(h.len(), st.iterations);
    assert!(h.windows(2).all(|w| w[1] <= w[0] * (1.0 + 1e-12)));

    let rst = IterMethodsConfig {
        gcr_restart: 3,
        ..cfg("gcr")
    };
    let (x2, st2) = solve_real(&mut Gcr::<f64>::new(&rst), &a, &b);
    assert_eq!(st2.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x2, &b) < 1e-9);
    assert!(st2.iterations >= st.iterations);
}

#[test]
fn gcr_min_iterations_is_real_only() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let c = IterMethodsConfig {
        min_iterations: 100,
        gcr_restart: 500,
        ..cfg("gcr")
    };
    let (_, st) = solve_real(&mut Gcr::<f64>::new(&c), &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(st.iterations >= 100);

    let ah = helmholtz(N);
    let bh = rhs_complex(N);
    let mut x = vec![Complex64::new(0.0, 0.0); N];
    let st = Gcr::<Complex64>::new(&c)
        .solve(&ah, &IdentityPrecond, &bh, &mut x, &SerialComm)
        .unwrap();
    assert!(st.iterations < 100);
}

#[test]
fn gcr_custom_stop_criterion_replaces_the_residual() {
    let a = nonsym(N);
    let b = rhs_real(N);
    // Always "converged": stops after the first iteration.
    let mut g = Gcr::<f64>::new(&cfg("gcr")).with_stop_criterion(Box::new(|_, _, _| 0.0));
    let (_, st) = solve_real(&mut g, &a, &b);
    assert_eq!(st.iterations, 0);
    assert_eq!(st.status, SolveStatus::Converged);
}

#[test]
fn idrs_solves_for_several_s() {
    let a = nonsym(N);
    let b = rhs_real(N);
    for s in [1, 2, 4, 6] {
        let c = IterMethodsConfig {
            idrs_s: s,
            ..cfg("idrs")
        };
        let (x, st) = solve_real(&mut Idrs::<f64>::new(&c), &a, &b);
        assert_eq!(st.status, SolveStatus::Converged, "s={s}");
        assert!(rel_resid(&a, &x, &b) < 1e-8, "s={s}");
    }
}

#[test]
fn idrs_smoothing_and_robust_modes_converge() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let c = IterMethodsConfig {
        idrs_smoothing: true,
        ..cfg("idrs")
    };
    let mut s = Idrs::<f64>::new(&c);
    let (x, st) = solve_real(&mut s, &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-8);
    // Smoothed residuals never increase.
    assert!(s.history().windows(2).all(|w| w[1] <= w[0] * (1.0 + 1e-10)));

    let c = IterMethodsConfig {
        robust: Some(RobustConfig {
            tolerance: 1e-6,
            max_tolerance: 1.0,
            step_size: 1.1,
            max_bad_iterations: 50,
            start_iteration: 1,
        }),
        ..cfg("idrs")
    };
    let (x, st) = solve_real(&mut Idrs::<f64>::new(&c), &a, &b);
    assert_eq!(st.status, SolveStatus::Converged);
    assert!(rel_resid(&a, &x, &b) < 1e-6);
}

#[test]
fn idrs_is_deterministic_for_a_seed() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let run = |seed| {
        let c = IterMethodsConfig {
            random_seed: seed,
            ..cfg("idrs")
        };
        solve_real(&mut Idrs::<f64>::new(&c), &a, &b)
    };
    let (x1, s1) = run(7);
    let (x2, s2) = run(7);
    assert_eq!(x1, x2);
    assert_eq!(s1.iterations, s2.iterations);
}

#[test]
fn complex_methods_solve_helmholtz() {
    let a = helmholtz(N);
    let b = rhs_complex(N);
    let solvers: Vec<(&str, Box<dyn IterativeSolver<Complex64>>)> = vec![
        ("bicgstabl", cfg("bicgstabl").build_complex().unwrap()),
        ("gcr", cfg("gcr").build_complex().unwrap()),
        ("idrs", cfg("idrs").build_complex().unwrap()),
        // bicgstab2 on a complex system maps to BiCGStab(l)
        ("bicgstab2", cfg("bicgstab2").build_complex().unwrap()),
    ];
    for (name, mut s) in solvers {
        let mut x = vec![Complex64::new(0.0, 0.0); N];
        let st = s
            .solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm)
            .unwrap();
        assert_eq!(st.status, SolveStatus::Converged, "{name}");
        assert!(rel_resid(&a, &x, &b) < 1e-8, "{name}");
    }
    assert!(cfg("cg").build_complex().is_err());
}

#[test]
fn pseudo_complex_gcr_reproduces_complex_gcr() {
    let n = 24;
    let ac = helmholtz(n);
    let bc = rhs_complex(n);
    // Real 2n x 2n form [[A, -B], [B, A]] over interleaved (re, im) unknowns.
    let mut real = elmer_core::CsrMatrix::<f64> {
        nrows: 2 * n,
        ncols: 2 * n,
        row_ptr: vec![0],
        col_idx: vec![],
        values: vec![],
        diag: None,
    };
    for i in 0..n {
        let (s, e) = (ac.row_ptr[i], ac.row_ptr[i + 1]);
        for part in 0..2 {
            for k in s..e {
                let z = ac.values[k];
                let j = ac.col_idx[k];
                // row 2i: re = A_re*x_re - A_im*x_im; row 2i+1: im = A_im*x_re + A_re*x_im
                if part == 0 {
                    real.col_idx.push(2 * j);
                    real.values.push(z.re);
                    real.col_idx.push(2 * j + 1);
                    real.values.push(-z.im);
                } else {
                    real.col_idx.push(2 * j);
                    real.values.push(z.im);
                    real.col_idx.push(2 * j + 1);
                    real.values.push(z.re);
                }
            }
            real.row_ptr.push(real.col_idx.len());
        }
    }
    let br: Vec<f64> = bc.iter().flat_map(|z| [z.re, z.im]).collect();

    let c = IterMethodsConfig {
        gcr_restart: 500,
        ..cfg("gcr")
    };
    let mut xc = vec![Complex64::new(0.0, 0.0); n];
    let mut gz = Gcr::<Complex64>::new(&c);
    let sz = gz
        .solve(&ac, &IdentityPrecond, &bc, &mut xc, &SerialComm)
        .unwrap();

    for mode in [PseudoComplex::Interleaved] {
        let cp = IterMethodsConfig {
            pseudo_complex: mode,
            ..c.clone()
        };
        let mut gp = Gcr::<f64>::new(&cp);
        let (xr, sp) = solve_real(&mut gp, &real, &br);
        assert_eq!(sp.status, SolveStatus::Converged);
        // A pseudo-complex GCR step is the complex GCR step: same count, same solution.
        assert_eq!(sp.iterations, sz.iterations);
        for i in 0..n {
            assert!((xr[2 * i] - xc[i].re).abs() < 1e-8);
            assert!((xr[2 * i + 1] - xc[i].im).abs() < 1e-8);
        }
    }
}

#[test]
fn pseudo_complex_is_refused_by_other_methods() {
    let a = nonsym(4);
    let c = IterMethodsConfig {
        pseudo_complex: PseudoComplex::Interleaved,
        ..cfg("idrs")
    };
    let mut x = vec![0.0; 4];
    assert!(Idrs::<f64>::new(&c)
        .solve(&a, &IdentityPrecond, &[1.0; 4], &mut x, &SerialComm)
        .is_err());
}

#[test]
fn initial_guess_already_converged_returns_zero_iterations() {
    let a = nonsym(N);
    let b = rhs_real(N);
    // Exact solution via a long GCR run, then restart from it.
    let (x0, _) = solve_real(
        &mut Gcr::<f64>::new(&IterMethodsConfig {
            gcr_restart: 500,
            tolerance: 1e-14,
            ..cfg("gcr")
        }),
        &a,
        &b,
    );
    for m in ["gcr", "idrs"] {
        let mut x = x0.clone();
        let mut s = cfg(m).build_real(&a).unwrap();
        let st = s.solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
        assert_eq!(st.iterations, 0, "{m}");
        assert_eq!(st.status, SolveStatus::Converged, "{m}");
        assert_eq!(x, x0, "{m}");
    }
}

#[test]
fn divergence_limit_flags_divergence() {
    let a = nonsym(N);
    let b = rhs_real(N);
    let c = IterMethodsConfig {
        divergence_limit: 1e-3, // below the initial relative residual of 1
        ..cfg("gcr")
    };
    let (_, st) = solve_real(&mut Gcr::<f64>::new(&c), &a, &b);
    assert_eq!(st.status, SolveStatus::Diverged);
}

#[test]
fn dimension_mismatch_is_an_error() {
    let a = nonsym(5);
    let mut x = vec![0.0; 4];
    assert!(Gcr::<f64>::new(&cfg("gcr"))
        .solve(&a, &IdentityPrecond, &[1.0; 4], &mut x, &SerialComm)
        .is_err());
}
