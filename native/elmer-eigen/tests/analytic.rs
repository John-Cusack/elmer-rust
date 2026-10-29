//! Problems with closed-form spectra: 1-D finite differences and linear
//! finite elements, symmetric and not.
mod common;
use common::*;

use elmer_core::{CsrMatrix, LinalgError, LinearOperator};
use elmer_direct::FaerLu;
use elmer_eigen::{eigs, EigenConfig, Mass, Normalize, Which};
use num_complex::Complex64;
use std::f64::consts::PI;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

fn sorted_by_magnitude(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
    v
}

/// max_i ||K x_i - lambda_i M x_i|| / (|lambda_i| ||M x_i||), complex vectors.
fn max_residual(
    k: &CsrMatrix<f64>,
    mass: Mass<'_>,
    values: &[Complex64],
    vectors: &[Vec<Complex64>],
) -> f64 {
    let n = k.nrows;
    let mut worst = 0.0_f64;
    for (lam, x) in values.iter().zip(vectors) {
        let (xr, xi): (Vec<f64>, Vec<f64>) = x.iter().map(|z| (z.re, z.im)).unzip();
        let mut r = vec![Complex64::new(0.0, 0.0); n];
        let mut mxn = 0.0;
        for (part, unit) in [
            (&xr, Complex64::new(1.0, 0.0)),
            (&xi, Complex64::new(0.0, 1.0)),
        ] {
            let (mut kx, mut mx) = (vec![0.0; n], vec![0.0; n]);
            k.apply(part, &mut kx);
            match mass {
                Mass::Identity => mx.copy_from_slice(part),
                Mass::Lumped(d) => mx
                    .iter_mut()
                    .zip(part.iter().zip(d))
                    .for_each(|(o, (a, b))| *o = a * b),
                Mass::Consistent(m) => m.apply(part, &mut mx),
            }
            for i in 0..n {
                r[i] += unit * (Complex64::new(kx[i], 0.0) - lam * mx[i]);
            }
            mxn += mx.iter().map(|v| v * v).sum::<f64>();
        }
        let rn = r.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        worst = worst.max(rn / (lam.norm() * mxn.sqrt()));
    }
    worst
}

#[test]
fn laplacian_standard_problem() {
    let n = 60;
    let k = tridiag(n, -1.0, 2.0, -1.0);
    let mut cfg = EigenConfig::new(6, true);
    cfg.tol = 1e-12;
    let res = eigs(&k, Mass::Identity, &mut FaerLu::new(), &cfg).unwrap();
    for (j, v) in res.values.iter().enumerate() {
        let exact = 2.0 - 2.0 * ((j + 1) as f64 * PI / (n + 1) as f64).cos();
        assert!(rel(v.re, exact) < 1e-12, "{j}: {} vs {exact}", v.re);
        assert_eq!(v.im, 0.0);
    }
    assert!(max_residual(&k, Mass::Identity, &res.values, &res.vectors) < 1e-10);
    // unit 2-norm vectors (B = I)
    for x in &res.vectors {
        let nrm: f64 = x.iter().map(|z| z.norm_sqr()).sum();
        assert!((nrm - 1.0).abs() < 1e-10);
    }
}

/// Linear elements on (0,1) with Elmer-style Dirichlet rows at both ends: the
/// row is cut to its diagonal and the mass row zeroed, the columns kept.
fn fe_system(ne: usize) -> (CsrMatrix<f64>, CsrMatrix<f64>, f64) {
    let n = ne + 1;
    let h = 1.0 / ne as f64;
    let mut k = tridiag(n, -1.0 / h, 2.0 / h, -1.0 / h);
    let mut m = tridiag(n, h / 6.0, 4.0 * h / 6.0, h / 6.0);
    for bc in [0, n - 1] {
        for p in k.row_ptr[bc]..k.row_ptr[bc + 1] {
            let diag = k.col_idx[p] as usize == bc;
            k.values[p] = if diag { 1.0 } else { 0.0 };
            m.values[p] = 0.0;
        }
    }
    (k, m, h)
}

fn fe_exact(ne: usize, j: usize) -> f64 {
    let h = 1.0 / ne as f64;
    let c = (j as f64 * PI / ne as f64).cos();
    6.0 / (h * h) * (1.0 - c) / (2.0 + c)
}

#[test]
fn consistent_mass_with_dirichlet_rows() {
    let ne = 40;
    let (k, m, _) = fe_system(ne);
    for symmetric in [true, false] {
        let mut cfg = EigenConfig::new(5, symmetric);
        cfg.tol = 1e-10;
        let res = eigs(&k, Mass::Consistent(&m), &mut FaerLu::new(), &cfg).unwrap();
        for (j, v) in res.values.iter().enumerate() {
            let exact = fe_exact(ne, j + 1);
            assert!(
                rel(v.re, exact) < 1e-11,
                "sym {symmetric} {j}: {} vs {exact}",
                v.re
            );
            assert!(v.im.abs() < 1e-9 * exact);
        }
        for x in &res.vectors {
            // Dirichlet components vanish; x^H M x = 1.
            assert!(x[0].norm() < 1e-12 && x[ne].norm() < 1e-12);
            let xr = re(x);
            let mut mx = vec![0.0; ne + 1];
            m.apply(&xr, &mut mx);
            let xm: f64 = xr.iter().zip(&mx).map(|(a, b)| a * b).sum();
            let xi: f64 = x.iter().map(|z| z.im * z.im).sum();
            assert!((xm - 1.0).abs() < 1e-8, "x^T M x = {xm}, imag {xi}");
        }
        assert!(max_residual(&k, Mass::Consistent(&m), &res.values, &res.vectors) < 1e-9);
    }
}

#[test]
fn shift_selects_nearest_eigenvalues() {
    let ne = 40;
    let (k, m, _) = fe_system(ne);
    let mut cfg = EigenConfig::new(3, true);
    cfg.tol = 1e-12;
    let target = fe_exact(ne, 10);
    cfg.shift = target * 1.001;
    let res = eigs(&k, Mass::Consistent(&m), &mut FaerLu::new(), &cfg).unwrap();
    let mut want = vec![fe_exact(ne, 9), fe_exact(ne, 10), fe_exact(ne, 11)];
    want = sorted_by_magnitude(want);
    for (v, w) in res.values.iter().zip(&want) {
        assert!(rel(v.re, *w) < 1e-11, "{} vs {w}", v.re);
    }
    // Eigen System Sorting = largest magnitude reverses the order.
    cfg.sorting = Which::LargestMagnitude;
    let res = eigs(&k, Mass::Consistent(&m), &mut FaerLu::new(), &cfg).unwrap();
    assert!(res.values[0].re > res.values[1].re && res.values[1].re > res.values[2].re);
}

#[test]
fn lumped_mass_regular_mode() {
    // Lumped linear elements, no Dirichlet rows (mode 2 needs M nonsingular):
    // K = tridiag(-1,2,-1)/h with M = h I has lambda = (2 - 2 cos)/h^2.
    let n = 24;
    let h = 0.1;
    let k = tridiag(n, -1.0 / h, 2.0 / h, -1.0 / h);
    let d = vec![h; n];
    let mut cfg = EigenConfig::new(3, true);
    cfg.tol = 1e-12;
    cfg.shift = 123.0; // ignored for a lumped mass, as in Elmer
    let res = eigs(&k, Mass::Lumped(&d), &mut FaerLu::new(), &cfg).unwrap();
    for (j, v) in res.values.iter().enumerate() {
        let exact = (2.0 - 2.0 * ((j + 1) as f64 * PI / (n + 1) as f64).cos()) / (h * h);
        assert!(rel(v.re, exact) < 1e-10, "{j}: {} vs {exact}", v.re);
    }
    assert!(max_residual(&k, Mass::Lumped(&d), &res.values, &res.vectors) < 1e-9);
}

#[test]
fn nonsymmetric_real_spectrum() {
    // Toeplitz tridiag(a, b, c): lambda_j = b + 2 sqrt(ac) cos(j pi / (n+1)).
    // The eigenvector condition number is |a/c|^((n-1)/2), about 140 here.
    let n = 50;
    let g = 0.1;
    let (a, b, c) = (-1.0 - g, 2.0, -1.0 + g);
    let k = tridiag(n, a, b, c);
    let mut cfg = EigenConfig::new(4, false);
    cfg.tol = 1e-12;
    let res = eigs(&k, Mass::Identity, &mut FaerLu::new(), &cfg).unwrap();
    let s = (a * c).sqrt();
    let mut exact: Vec<f64> = (1..=n)
        .map(|j| b + 2.0 * s * (j as f64 * PI / (n + 1) as f64).cos())
        .collect();
    exact = sorted_by_magnitude(exact);
    for (v, e) in res.values.iter().zip(&exact) {
        assert!(rel(v.re, *e) < 1e-9, "{} vs {e}", v.re);
        assert!(v.im.abs() < 1e-9);
    }
    assert!(max_residual(&k, Mass::Identity, &res.values, &res.vectors) < 1e-8);
}

#[test]
fn nonsymmetric_complex_pairs() {
    // ac < 0: lambda_j = b + 2 i sqrt(-ac) cos(j pi/(n+1)), conjugate pairs.
    // a = -c keeps the matrix normal (2 I plus a skew part), so the pairs are
    // well conditioned.
    let n = 41;
    let (a, b, c) = (-1.0, 2.0, 1.0);
    let k = tridiag(n, a, b, c);
    let mut cfg = EigenConfig::new(5, false);
    cfg.tol = 1e-12;
    cfg.sorting = Which::LargestImag;
    let res = eigs(&k, Mass::Identity, &mut FaerLu::new(), &cfg).unwrap();
    let s = (-(a * c)).sqrt();
    // Smallest magnitude (LM on 1/lambda): lambda = 2 and the two pairs next to
    // it. (Selecting on the imaginary part, 'SI' on 1/lambda, stalls here as it
    // can in ARPACK: unresolved real Ritz values all have imaginary part 0.)
    let exact: Vec<Complex64> = (1..=n)
        .map(|j| Complex64::new(b, 2.0 * s * (j as f64 * PI / (n + 1) as f64).cos()))
        .collect();
    let mut wanted = exact.clone();
    wanted.sort_by(|p, q| p.norm().total_cmp(&q.norm()));
    for v in &res.values {
        assert!(
            wanted[..5].iter().any(|e| (v - e).norm() < 1e-9 * e.norm()),
            "{v} is not among the 5 wanted {:?}",
            &wanted[..5]
        );
    }
    for v in &res.values {
        let best = exact
            .iter()
            .map(|e| (v - e).norm() / e.norm())
            .fold(f64::MAX, f64::min);
        assert!(best < 1e-10, "{v} not in the spectrum ({best:e})");
    }
    // Sorted largest imaginary part first; a pair's members both present.
    for w in res.values.windows(2) {
        assert!(w[0].im >= w[1].im);
    }
    let pos = res.values.iter().filter(|z| z.im > 1e-8).count();
    let neg = res.values.iter().filter(|z| z.im < -1e-8).count();
    assert!(pos >= 2 && neg >= 2, "{:?}", res.values);
    assert!(max_residual(&k, Mass::Identity, &res.values, &res.vectors) < 1e-9);
}

#[test]
fn normalize_to_unity_and_residual_report() {
    let (k, m, _) = fe_system(30);
    let mut cfg = EigenConfig::new(3, true);
    cfg.tol = 1e-12;
    cfg.normalize = Normalize::Unity;
    cfg.compute_residuals = true;
    let res = eigs(&k, Mass::Consistent(&m), &mut FaerLu::new(), &cfg).unwrap();
    for x in &res.vectors {
        let big = x.iter().fold(0.0_f64, |a, z| a.max(z.norm()));
        assert!((big - 1.0).abs() < 1e-14);
    }
    let r = res.residuals.unwrap();
    assert_eq!(r.len(), 3);
    // computed before the scaling, on unit M-norm vectors
    assert!(
        r.iter().zip(&res.values).all(|(ri, v)| *ri < 1e-8 * v.re),
        "{r:?}"
    );
}

#[test]
fn argument_errors() {
    let k = tridiag(20, -1.0, 2.0, -1.0);
    let mut lu = FaerLu::new();
    let mut cfg = EigenConfig::new(5, true);
    cfg.ncv = 5;
    assert!(matches!(
        eigs(&k, Mass::Identity, &mut lu, &cfg),
        Err(LinalgError::Dim(_))
    ));
    cfg.ncv = 21; // > n
    assert!(matches!(
        eigs(&k, Mass::Identity, &mut lu, &cfg),
        Err(LinalgError::Dim(_))
    ));
    let mut cfg = EigenConfig::new(2, true);
    cfg.select = Which::SmallestReal; // DSAUPD info = -5
    assert!(matches!(
        eigs(&k, Mass::Identity, &mut lu, &cfg),
        Err(LinalgError::Unsupported(_))
    ));
    let mut cfg = EigenConfig::new(2, false);
    cfg.ncv = 3; // DNAUPD needs NCV - NEV >= 2
    assert!(matches!(
        eigs(&k, Mass::Identity, &mut lu, &cfg),
        Err(LinalgError::Dim(_))
    ));
    let d = vec![1.0, 0.0];
    let k2 = tridiag(2, -1.0, 2.0, -1.0);
    let cfg = EigenConfig::new(1, true);
    assert!(matches!(
        eigs(&k2, Mass::Lumped(&d), &mut lu, &cfg),
        Err(LinalgError::Singular(1))
    ));
}

#[test]
fn max_iterations_is_an_error() {
    // Smallest-magnitude without shift-invert (lumped, mode 2) converges slowly.
    let n = 200;
    let k = tridiag(n, -1.0, 2.0, -1.0);
    let d = vec![1.0; n];
    let mut cfg = EigenConfig::new(4, true);
    cfg.ncv = 6;
    cfg.max_iter = 1;
    cfg.tol = 1e-14;
    match eigs(&k, Mass::Lumped(&d), &mut FaerLu::new(), &cfg) {
        Err(LinalgError::Breakdown { detail, .. }) => {
            assert!(detail.contains("Maximum number of iterations"))
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn small_invariant_subspace_restarts_cleanly() {
    // Diagonal matrix: a random start sees every eigenvalue, ncv = n.
    let n = 8;
    let mut k = tridiag(n, 0.0, 1.0, 0.0);
    for i in 0..n {
        for p in k.row_ptr[i]..k.row_ptr[i + 1] {
            if k.col_idx[p] as usize == i {
                k.values[p] = (i + 1) as f64;
            }
        }
    }
    let mut cfg = EigenConfig::new(3, true);
    cfg.ncv = n;
    let res = eigs(&k, Mass::Identity, &mut FaerLu::new(), &cfg).unwrap();
    let got: Vec<f64> = res.values.iter().map(|z| z.re).collect();
    for (g, e) in got.iter().zip([1.0, 2.0, 3.0]) {
        assert!((g - e).abs() < 1e-12, "{got:?}");
    }
    let _ = Lcg(1).next_f64();
}
