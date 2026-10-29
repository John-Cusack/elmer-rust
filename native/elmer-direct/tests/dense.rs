mod common;
use common::*;

use elmer_direct::dense::{geev, gesv, gesvd, getrf, getrs, potrf, potrs, syev, sygv};
use num_complex::Complex64;

fn rand_mat(m: usize, n: usize, seed: u64) -> Vec<f64> {
    let mut g = Lcg(seed);
    (0..m * n).map(|_| g.next_f64()).collect()
}

fn matmul(m: usize, k: usize, n: usize, a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0; m * n];
    for j in 0..n {
        for p in 0..k {
            for i in 0..m {
                c[j * m + i] += a[p * m + i] * b[j * k + p];
            }
        }
    }
    c
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn getrf_reconstructs_pa_equals_lu() {
    for (m, n) in [(6, 6), (7, 4), (4, 7)] {
        let a0 = rand_mat(m, n, 31 + m as u64);
        let mut a = a0.clone();
        let ipiv = getrf(m, n, &mut a).unwrap();
        let k = m.min(n);
        assert_eq!(ipiv.len(), k);
        let mut l = vec![0.0; m * k];
        let mut u = vec![0.0; k * n];
        for j in 0..k {
            for i in 0..m {
                l[j * m + i] = if i == j { 1.0 } else if i > j { a[j * m + i] } else { 0.0 };
            }
        }
        for j in 0..n {
            for i in 0..k.min(j + 1) {
                u[j * k + i] = a[j * m + i];
            }
        }
        let mut lu = matmul(m, k, n, &l, &u);
        // Undo the interchanges (apply in reverse) to get A back.
        for i in (0..k).rev() {
            for c in 0..n {
                lu.swap(c * m + i, c * m + ipiv[i]);
            }
        }
        assert!(max_diff(&lu, &a0) < 1e-12, "{m}x{n}");
    }
}

#[test]
fn getrs_gesv_solve() {
    let n = 9;
    let a0 = rand_mat(n, n, 5);
    let x0 = rand_mat(n, 2, 6);
    let mut b = matmul(n, n, 2, &a0, &x0);
    let mut a = a0.clone();
    gesv(n, 2, &mut a, &mut b).unwrap();
    assert!(max_diff(&b, &x0) < 1e-10);

    let mut a = a0.clone();
    let ipiv = getrf(n, n, &mut a).unwrap();
    let mut b = matmul(n, n, 2, &a0, &x0);
    getrs(n, 2, &a, &ipiv, &mut b).unwrap();
    assert!(max_diff(&b, &x0) < 1e-10);
}

#[test]
fn getrf_complex_and_singular() {
    let n = 5;
    let a0: Vec<Complex64> = rand_mat(n, n, 1).iter().zip(rand_mat(n, n, 2)).map(|(&r, i)| Complex64::new(r, i)).collect();
    let x0: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 1.0)).collect();
    let mut b = vec![Complex64::new(0.0, 0.0); n];
    for j in 0..n {
        for i in 0..n {
            b[i] += a0[j * n + i] * x0[j];
        }
    }
    let mut a = a0.clone();
    gesv(n, 1, &mut a, &mut b).unwrap();
    for (p, q) in b.iter().zip(&x0) {
        assert!((p - q).norm() < 1e-11);
    }
    let mut z = vec![0.0; 9];
    assert!(getrf(3, 3, &mut z).is_err());
}

#[test]
fn potrf_potrs() {
    let n = 8;
    let g = rand_mat(n, n, 12);
    // A = G G^T + n I
    let mut a = matmul(n, n, n, &g, &{
        let mut t = vec![0.0; n * n];
        for i in 0..n { for j in 0..n { t[j * n + i] = g[i * n + j]; } }
        t
    });
    for i in 0..n { a[i * n + i] += n as f64; }
    let x0 = rand_mat(n, 3, 13);
    let mut b = matmul(n, n, 3, &a, &x0);
    let mut l = a.clone();
    potrf(n, &mut l).unwrap();
    potrs(n, 3, &l, &mut b).unwrap();
    assert!(max_diff(&b, &x0) < 1e-11);
    let mut bad = vec![1.0, 2.0, 2.0, 1.0];
    assert!(potrf(2, &mut bad).is_err());
}

fn sym(n: usize, seed: u64) -> Vec<f64> {
    let g = rand_mat(n, n, seed);
    let mut a = vec![0.0; n * n];
    for i in 0..n { for j in 0..n { a[j * n + i] = 0.5 * (g[j * n + i] + g[i * n + j]); } }
    a
}

#[test]
fn syev_residuals_and_order() {
    let n = 10;
    let a = sym(n, 20);
    let e = syev(n, &a).unwrap();
    assert!(e.values.windows(2).all(|w| w[0] <= w[1]));
    let av = matmul(n, n, n, &a, &e.vectors);
    for k in 0..n {
        for i in 0..n {
            assert!((av[k * n + i] - e.values[k] * e.vectors[k * n + i]).abs() < 1e-12);
        }
    }
}

#[test]
fn sygv_b_orthonormal() {
    let n = 7;
    let a = sym(n, 21);
    let g = rand_mat(n, n, 22);
    let mut b = vec![0.0; n * n];
    for i in 0..n { for j in 0..n { for p in 0..n { b[j * n + i] += g[p * n + i] * g[p * n + j]; } } }
    for i in 0..n { b[i * n + i] += 1.0; }
    let e = sygv(n, &a, &b).unwrap();
    let ax = matmul(n, n, n, &a, &e.vectors);
    let bx = matmul(n, n, n, &b, &e.vectors);
    for k in 0..n {
        for i in 0..n {
            assert!((ax[k * n + i] - e.values[k] * bx[k * n + i]).abs() < 1e-9, "k={k}");
        }
    }
    for p in 0..n {
        for q in 0..n {
            let d: f64 = (0..n).map(|i| e.vectors[p * n + i] * bx[q * n + i]).sum();
            assert!((d - if p == q { 1.0 } else { 0.0 }).abs() < 1e-10);
        }
    }
}

#[test]
fn geev_residuals_and_normalisation() {
    let n = 8;
    let a = rand_mat(n, n, 40);
    let e = geev(n, &a).unwrap();
    for k in 0..n {
        let v = &e.vectors[k * n..(k + 1) * n];
        let norm: f64 = v.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-12);
        let big = v.iter().copied().fold(Complex64::new(0.0, 0.0), |m, z| if z.norm_sqr() > m.norm_sqr() { z } else { m });
        assert!(big.im.abs() < 1e-12 && big.re > 0.0);
        for i in 0..n {
            let av: Complex64 = (0..n).map(|j| a[j * n + i] * v[j]).sum();
            assert!((av - e.values[k] * v[i]).norm() < 1e-10);
        }
    }
    // A rotation has the pair +/- i, positive imaginary first.
    let r = geev(2, &[0.0, 1.0, -1.0, 0.0]).unwrap();
    assert!(r.values[0].im > 0.0 && r.values[1].im < 0.0);
}

#[test]
fn gesvd_reconstructs() {
    for (m, n) in [(6, 4), (4, 6), (5, 5)] {
        let a = rand_mat(m, n, 50 + m as u64);
        let s = gesvd(m, n, &a).unwrap();
        let k = m.min(n);
        assert_eq!(s.s.len(), k);
        assert!(s.s.windows(2).all(|w| w[0] >= w[1]));
        let mut rec = vec![0.0; m * n];
        for j in 0..n {
            for p in 0..k {
                for i in 0..m {
                    rec[j * m + i] += s.u[p * m + i] * s.s[p] * s.vt[j * n + p];
                }
            }
        }
        assert!(max_diff(&rec, &a) < 1e-12, "{m}x{n}");
    }
}
