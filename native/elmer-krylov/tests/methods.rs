//! Behaviour of each HUTI method on small manufactured systems: convergence,
//! the complex variants, stopping criteria and the breakdown/limit exits.
//! These check local correctness only; Elmer parity is in fortran_parity.rs.

mod common;

use std::sync::{Arc, Mutex};

use common::csr;
use elmer_core::{
    Communicator, CsrMatrix, IdentityPrecond, IterativeSolver, LinalgError, Preconditioner,
    SerialComm, SolveStatus,
};
use elmer_krylov::{Huti, HutiConfig, HutiInfo, HutiMethod, HutiStop};
use elmer_sparse::{CsrOp, ElmerScalar, Ilu};
use num_complex::Complex64;

const ALL: [HutiMethod; 7] = [
    HutiMethod::Cg,
    HutiMethod::Cgs,
    HutiMethod::BiCGStab,
    HutiMethod::BiCGStab2,
    HutiMethod::Gmres,
    HutiMethod::Qmr,
    HutiMethod::Tfqmr,
];

/// 5-point Laplacian on an `nx` x `ny` grid (1D when `ny == 1`), plus `shift` on
/// the diagonal and `conv` times a forward difference in x (nonsymmetric).
fn grid<T: elmer_core::Scalar>(nx: usize, ny: usize, shift: T, conv: f64) -> CsrMatrix<T> {
    let n = nx * ny;
    let f = T::from_f64;
    let mut t = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let r = j * nx + i;
            let diag = f(if ny == 1 { 2.0 } else { 4.0 }) + shift + f(conv);
            if j > 0 {
                t.push((r, (r - nx) as u32, f(-1.0)));
            }
            if i > 0 {
                t.push((r, (r - 1) as u32, f(-1.0 - conv)));
            }
            t.push((r, r as u32, diag));
            if i + 1 < nx {
                t.push((r, (r + 1) as u32, f(-1.0)));
            }
            if j + 1 < ny {
                t.push((r, (r + nx) as u32, f(-1.0)));
            }
        }
    }
    csr(n, &t)
}

fn cfg(method: HutiMethod, tol: f64, max_iter: usize) -> HutiConfig {
    let mut c = HutiConfig::new(method);
    c.tol = tol;
    c.max_iter = max_iter;
    c.restart = 20;
    c
}

fn residual<T: elmer_krylov::HutiScalar>(a: &CsrMatrix<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::zero(); b.len()];
    elmer_sparse::ops::spmv(a, x, &mut ax);
    let r: Vec<T> = ax.iter().zip(b).map(|(p, q)| *p - *q).collect();
    T::nrm2(&r) / T::nrm2(b)
}

fn solve<T: elmer_krylov::HutiScalar + ElmerScalar>(
    c: HutiConfig,
    a: &CsrMatrix<T>,
    m: &dyn Preconditioner<T>,
    b: &[T],
    x: &mut [T],
) -> (elmer_core::SolveStats, Huti<T>) {
    let mut s = Huti::new(c);
    s.set_output(Box::new(|_| {}));
    let st = s.solve(&CsrOp(a), m, b, x, &SerialComm).unwrap();
    (st, s)
}

#[test]
fn every_method_solves_poisson_1d_and_2d() {
    for (nx, ny) in [(40, 1), (12, 12)] {
        let a = grid(nx, ny, 0.0, 0.0);
        let b: Vec<f64> = (0..nx * ny).map(|i| 1.0 + (i % 7) as f64).collect();
        for m in ALL {
            let mut x = vec![0.0; b.len()];
            let (st, _) = solve(cfg(m, 1e-10, 2000), &a, &IdentityPrecond, &b, &mut x);
            assert_eq!(st.status, SolveStatus::Converged, "{m:?} {nx}x{ny}: {st:?}");
            assert!(residual(&a, &x, &b) < 1e-9, "{m:?} {nx}x{ny}");
        }
    }
}

#[test]
fn every_method_solves_a_nonsymmetric_system_with_both_preconditioner_sides() {
    let a = grid(10, 10, 0.0, 0.8);
    let b: Vec<f64> = (0..100).map(|i| (i as f64 * 0.37).sin()).collect();
    let m = Ilu::diagonal(&a);
    for method in ALL.into_iter().filter(|m| *m != HutiMethod::Cg) {
        for left in [false, true] {
            if method == HutiMethod::BiCGStab2 && !left {
                continue; // see bicgstab2_right_preconditioning_returns_the_inner_variable
            }
            let mut c = cfg(method, 1e-10, 2000);
            c.left_precond = left;
            let mut x = vec![0.0; 100];
            let (st, _) = solve(c, &a, &m, &b, &mut x);
            assert_eq!(st.status, SolveStatus::Converged, "{method:?} left={left}");
            assert!(residual(&a, &x, &b) < 1e-9, "{method:?} left={left}");
        }
    }
}

#[test]
fn complex_variants_solve_a_complex_symmetric_system() {
    // Helmholtz-like: Laplacian plus i*0.5 on the diagonal, complex symmetric
    // (A = A^T, not Hermitian), as Elmer's harmonic systems are.
    let a = grid(8, 8, Complex64::new(-0.3, 0.5), 0.0);
    let b: Vec<Complex64> = (0..64)
        .map(|i| Complex64::new(1.0 + (i % 3) as f64, (i % 5) as f64 - 2.0))
        .collect();
    for m in ALL {
        let mut x = vec![Complex64::new(0.0, 0.0); 64];
        let (st, _) = solve(cfg(m, 1e-10, 3000), &a, &IdentityPrecond, &b, &mut x);
        assert_eq!(st.status, SolveStatus::Converged, "{m:?}: {st:?}");
        assert!(residual(&a, &x, &b) < 1e-9, "{m:?}");
    }
}

#[test]
fn complex_gmres_prints_the_gmresz_tag() {
    let a = grid(6, 1, Complex64::new(0.0, 1.0), 0.0);
    let b = vec![Complex64::new(1.0, 0.0); 6];
    let mut c = cfg(HutiMethod::Gmres, 1e-12, 50);
    c.verbosity = 1;
    let lines = Arc::new(Mutex::new(Vec::<String>::new()));
    let sink = lines.clone();
    let mut s: Huti<Complex64> = Huti::new(c);
    s.set_output(Box::new(move |l| sink.lock().unwrap().push(l.into())));
    let mut x = vec![Complex64::new(0.0, 0.0); 6];
    s.solve(&CsrOp(&a), &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    let lines = lines.lock().unwrap();
    assert!(lines[0].starts_with("   gmresz:       1 "), "{lines:?}");
}

#[test]
fn zero_rhs_is_a_cg_rho_breakdown_with_x_untouched() {
    let a = grid(5, 1, 0.0, 0.0);
    let mut x = vec![0.0; 5];
    let (st, s) = solve(
        cfg(HutiMethod::Cg, 1e-8, 10),
        &a,
        &IdentityPrecond,
        &[0.0; 5],
        &mut x,
    );
    assert!(matches!(st.status, SolveStatus::Breakdown(_)));
    assert_eq!(s.info(), Some(HutiInfo::CgRho));
    assert_eq!(st.iterations, 1);
    assert_eq!(x, vec![0.0; 5]);
}

#[test]
fn exact_initial_guess_breaks_down_immediately() {
    let a = grid(5, 1, 0.0, 0.0);
    let x0 = [1.0, 2.0, 3.0, 2.0, 1.0];
    let mut b = [0.0; 5];
    elmer_sparse::ops::spmv(&a, &x0, &mut b);
    for (m, info) in [
        (HutiMethod::Cg, HutiInfo::CgRho),
        (HutiMethod::Cgs, HutiInfo::CgsRho),
        (HutiMethod::Gmres, HutiInfo::GmresAlpha),
        (HutiMethod::Tfqmr, HutiInfo::TfqmrRho),
        (HutiMethod::Qmr, HutiInfo::QmrRhoPsi),
    ] {
        let mut x = x0;
        let (_, s) = solve(cfg(m, 1e-8, 10), &a, &IdentityPrecond, &b, &mut x);
        assert_eq!(s.info(), Some(info), "{m:?}");
        assert_eq!(x, x0, "{m:?}");
    }
}

#[test]
fn bicgstab_perturbs_an_all_zero_guess() {
    let a = grid(5, 1, 0.0, 0.0);
    let b = [1.0; 5];
    // With the perturbation off, max_iter = 0 leaves one iteration from x = 0;
    // with it on, the same single iteration starts from 1e-8.
    let mut c = cfg(HutiMethod::BiCGStab, 1e-30, 0);
    c.perturb_zero_guess = false;
    let mut x0 = [0.0; 5];
    solve(c.clone(), &a, &IdentityPrecond, &b, &mut x0);
    c.perturb_zero_guess = true;
    let mut x1 = [0.0; 5];
    solve(c, &a, &IdentityPrecond, &b, &mut x1);
    assert_ne!(x0, x1);
}

#[test]
fn max_iterations_is_not_success() {
    let a = grid(30, 1, 0.0, 0.0);
    let b = [1.0; 30];
    for m in ALL {
        let mut x = [0.0; 30];
        let (st, _) = solve(cfg(m, 1e-14, 2), &a, &IdentityPrecond, &b, &mut x);
        assert_eq!(st.status, SolveStatus::MaxIterations, "{m:?}");
        assert_eq!(st.iterations, 3, "{m:?}: HUTI_ITERS is one past the limit");
    }
}

#[test]
fn divergence_limit_stops_the_iteration() {
    let a = grid(30, 1, 0.0, 0.0);
    let b = [1.0; 30];
    let mut c = cfg(HutiMethod::Cg, 1e-14, 100);
    c.max_tol = 1e-3; // below the first residual
    let mut x = [0.0; 30];
    let (st, _) = solve(c, &a, &IdentityPrecond, &b, &mut x);
    assert_eq!(st.status, SolveStatus::Diverged);
    assert_eq!(st.iterations, 1);
}

#[test]
fn b_scaled_and_preconditioned_b_scaled_criteria_differ() {
    // Every other row scaled by 1e4: ||b|| is ~1e4 times ||M^-1 b|| for the
    // diagonal M, so the same tolerance stops the two criteria apart.
    let mut a = grid(15, 15, 0.0, 0.0);
    let scale: Vec<f64> = (0..225)
        .map(|i| if i % 2 == 0 { 1e4 } else { 1.0 })
        .collect();
    for i in 0..225 {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            a.values[j] *= scale[i];
        }
    }
    let b = scale.clone();
    let m = Ilu::diagonal(&a);
    let mut iters = Vec::new();
    for stop in [
        HutiStop::PseudoResidualScaledByB,
        HutiStop::PseudoResidualScaledByPrecB,
    ] {
        let mut c = cfg(HutiMethod::BiCGStab, 1e-6, 500);
        c.stop_crit = stop;
        c.left_precond = true; // M^-1 b is formed with the left slot
        let mut x = vec![0.0; 225];
        let (st, _) = solve(c, &a, &m, &b, &mut x);
        assert_eq!(st.status, SolveStatus::Converged, "{stop:?}");
        iters.push(st.iterations);
    }
    assert!(iters[0] < iters[1], "{iters:?}");
}

#[test]
fn bicgstab2_right_preconditioning_returns_the_inner_variable() {
    // huti_bicgstab_2 iterates on y with x = M^-1 y and never maps back, which
    // is why IterSolve always runs it left-preconditioned (IterSolve.F90:1005).
    let a = grid(10, 10, 0.0, 0.8);
    let b: Vec<f64> = (0..100).map(|i| (i as f64 * 0.37).sin()).collect();
    let m = Ilu::diagonal(&a);
    let mut c = cfg(HutiMethod::BiCGStab2, 1e-10, 2000);
    c.left_precond = false;
    let mut y = vec![0.0; 100];
    let (st, _) = solve(c, &a, &m, &b, &mut y);
    assert_eq!(st.status, SolveStatus::Converged);
    let mut x = vec![0.0; 100];
    m.apply(&y, &mut x);
    assert!(residual(&a, &x, &b) < 1e-9);
    assert!(residual(&a, &y, &b) > 1e-3);
}

#[test]
fn every_stop_criterion_converges() {
    let a = grid(15, 1, 0.0, 0.3);
    let b: Vec<f64> = (0..15).map(|i| i as f64).collect();
    let m = Ilu::diagonal(&a);
    let stops = [
        HutiStop::TrueResidual,
        HutiStop::TrueResidualScaledByB,
        HutiStop::PseudoResidual,
        HutiStop::PseudoResidualScaledByB,
        HutiStop::PseudoResidualScaledByPrecB,
        HutiStop::XDiffNorm,
        HutiStop::UpperBound,
    ];
    for method in ALL.into_iter().filter(|m| *m != HutiMethod::Cg) {
        for stop in stops {
            if method == HutiMethod::Gmres && stop == HutiStop::XDiffNorm {
                // HUTI forms the update from Y(1:m) even after an early inner
                // exit at i < m, so stale Y entries keep it from shrinking.
                continue;
            }
            let mut c = cfg(method, 1e-9, 500);
            c.stop_crit = stop;
            let mut x = vec![0.0; 15];
            let (st, _) = solve(c, &a, &m, &b, &mut x);
            assert_eq!(st.status, SolveStatus::Converged, "{method:?} {stop:?}");
            assert!(residual(&a, &x, &b) < 1e-5, "{method:?} {stop:?}");
        }
    }
}

#[test]
fn gmres_restart_shorter_than_the_system_and_longer() {
    let a = grid(25, 1, 0.1, 0.0);
    let b = [1.0; 25];
    for restart in [1, 3, 40] {
        let mut c = cfg(HutiMethod::Gmres, 1e-10, 5000);
        c.restart = restart;
        let mut x = [0.0; 25];
        let (st, _) = solve(c, &a, &IdentityPrecond, &b, &mut x);
        assert_eq!(st.status, SolveStatus::Converged, "restart {restart}");
        assert!(residual(&a, &x, &b) < 1e-9);
    }
}

#[test]
fn qmr_needs_a_transpose() {
    struct NoTranspose(CsrMatrix<f64>);
    impl elmer_core::LinearOperator<f64> for NoTranspose {
        fn nrows(&self) -> usize {
            self.0.nrows
        }
        fn ncols(&self) -> usize {
            self.0.ncols
        }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            elmer_sparse::ops::spmv(&self.0, x, y)
        }
    }
    let a = NoTranspose(grid(5, 1, 0.0, 0.0));
    let mut s: Huti = Huti::new(cfg(HutiMethod::Qmr, 1e-8, 10));
    let err = s
        .solve(&a, &IdentityPrecond, &[1.0; 5], &mut [0.0; 5], &SerialComm)
        .unwrap_err();
    assert!(matches!(err, LinalgError::Unsupported(_)));
}

#[test]
fn dimension_mismatch_is_an_error() {
    let a = grid(5, 1, 0.0, 0.0);
    let mut s: Huti = Huti::new(cfg(HutiMethod::Cg, 1e-8, 10));
    let err = s
        .solve(
            &CsrOp(&a),
            &IdentityPrecond,
            &[1.0; 4],
            &mut [0.0; 4],
            &SerialComm,
        )
        .unwrap_err();
    assert!(matches!(err, LinalgError::Dim(_)));
}

#[test]
fn output_interval_and_final_line() {
    let a = grid(30, 1, 0.0, 0.0);
    let b = [1.0; 30];
    let mut c = cfg(HutiMethod::Cg, 1e-12, 100);
    c.verbosity = 4;
    let lines = Arc::new(Mutex::new(Vec::<String>::new()));
    let sink = lines.clone();
    let mut s: Huti = Huti::new(c);
    s.set_output(Box::new(move |l| sink.lock().unwrap().push(l.into())));
    let mut x = [0.0; 30];
    let st = s
        .solve(&CsrOp(&a), &IdentityPrecond, &b, &mut x, &SerialComm)
        .unwrap();
    let lines = lines.lock().unwrap();
    // Every 4th iteration, then the exit line.
    assert_eq!(lines.len(), st.iterations / 4 + 1);
    assert!(lines[0].starts_with("       4 "));
    assert_eq!(
        lines.last().unwrap(),
        &elmer_krylov::iteration_line(st.iterations, st.residual)
    );
    assert_eq!(s.history().len(), st.iterations);
}

#[test]
fn only_rank_zero_prints() {
    struct Rank1;
    impl Communicator for Rank1 {
        fn rank(&self) -> usize {
            1
        }
        fn size(&self) -> usize {
            1
        }
        fn sum_f64(&self, _: &mut [f64]) {}
        fn max_f64(&self, _: &mut [f64]) {}
        fn barrier(&self) {}
    }
    let a = grid(10, 1, 0.0, 0.0);
    let mut c = cfg(HutiMethod::Cg, 1e-10, 100);
    c.verbosity = 1;
    let lines = Arc::new(Mutex::new(0usize));
    let sink = lines.clone();
    let mut s: Huti = Huti::new(c);
    s.set_output(Box::new(move |_| *sink.lock().unwrap() += 1));
    s.solve(
        &CsrOp(&a),
        &IdentityPrecond,
        &[1.0; 10],
        &mut [0.0; 10],
        &Rank1,
    )
    .unwrap();
    assert_eq!(*lines.lock().unwrap(), 0);
}
