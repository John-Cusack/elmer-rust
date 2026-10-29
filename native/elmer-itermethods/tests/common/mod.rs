//! Small fixtures shared by the integration tests.
#![allow(dead_code)]

use elmer_core::{CsrMatrix, LinearOperator, Scalar};
use num_complex::Complex64;

/// Tridiagonal `n x n` matrix with the given sub-, main and super-diagonals.
pub fn tridiag<T: Scalar>(n: usize, lo: T, d: T, up: T) -> CsrMatrix<T> {
    let mut row_ptr = vec![0];
    let mut col_idx = Vec::new();
    let mut values = Vec::new();
    let mut diag = Vec::new();
    for i in 0..n {
        if i > 0 {
            col_idx.push((i - 1) as u32);
            values.push(lo);
        }
        diag.push(values.len());
        col_idx.push(i as u32);
        values.push(d);
        if i + 1 < n {
            col_idx.push((i + 1) as u32);
            values.push(up);
        }
        row_ptr.push(values.len());
    }
    CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr,
        col_idx,
        values,
        diag: Some(diag),
    }
}

/// `b - A x` in the max norm relative to `|b|`.
pub fn rel_resid<T: elmer_core::Scalar>(a: &CsrMatrix<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::zero(); b.len()];
    a.apply(x, &mut ax);
    let r: f64 = ax
        .iter()
        .zip(b)
        .map(|(p, q)| (*q - *p).abs().powi(2))
        .sum::<f64>()
        .sqrt();
    let bn: f64 = b.iter().map(|v| v.abs().powi(2)).sum::<f64>().sqrt();
    r / bn
}

pub fn rhs_real(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1.0 + (i as f64 * 0.7).sin()).collect()
}

pub fn rhs_complex(n: usize) -> Vec<Complex64> {
    (0..n)
        .map(|i| Complex64::new(1.0 + (i as f64 * 0.7).sin(), (i as f64 * 0.3).cos()))
        .collect()
}

/// Convection-diffusion style nonsymmetric matrix, diagonally dominant.
pub fn nonsym(n: usize) -> CsrMatrix<f64> {
    tridiag(n, -1.6, 2.4, -0.4)
}

/// Complex-symmetric Helmholtz-like matrix with damping, as in the
/// `linearsolvers_cmplx` family.
pub fn helmholtz(n: usize) -> CsrMatrix<Complex64> {
    tridiag(
        n,
        Complex64::new(-1.0, 0.0),
        Complex64::new(2.0 - 0.6, 0.5),
        Complex64::new(-1.0, 0.0),
    )
}

pub fn cfg(method: &str) -> elmer_itermethods::IterMethodsConfig {
    elmer_itermethods::IterMethodsConfig {
        method_name: method.into(),
        max_iterations: 500,
        tolerance: 1e-10,
        ..Default::default()
    }
}
