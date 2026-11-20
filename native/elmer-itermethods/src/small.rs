//! Tiny dense helpers for the BiCGStab(l) polynomial step: LU with partial
//! pivoting (`dgetrf`/`dgetrs`, `zgetrf`/`zgetrs`) on matrices of order l-1.
//!
//! The rounding follows reference LAPACK (pristine links the system netlib
//! 3.12 `liblapack`/`libblas`). `dgetrf` hands small matrices to the
//! recursive `dgetrf2`, which for these orders is right-looking elimination
//! that scales the pivot column by `ONE / pivot` (`DSCAL`, division only
//! below `SFMIN`). `dgetrs` is `DLASWP` and two column-oriented `DTRSM`s.
//! Order l-1 <= 2 (polynomial degree up to 3) rounds the same either way
//! for the solve; the reciprocal matters from order 2 on.

use crate::scalar::IterScalar;

/// BLAS `idamax` / `izamax` pivot magnitude: `|x|` real, `|re| + |im|` complex.
fn cabs1<T: IterScalar>(x: T) -> f64 {
    x.re().abs() + x.im().abs()
}

/// In-place LU of the `m x m` matrix `a` (row-major, `a[i][j]`), returning the
/// pivot rows. A zero pivot is left in place (LAPACK reports `info > 0` and the
/// callers here ignore it, so the solve then divides by zero as in Fortran).
pub(crate) fn lu_factor<T: IterScalar>(a: &mut [Vec<T>]) -> Vec<usize> {
    let m = a.len();
    let mut piv = vec![0usize; m];
    for k in 0..m {
        let mut p = k;
        let mut best = cabs1(a[k][k]);
        for (i, row) in a.iter().enumerate().skip(k + 1) {
            let c = cabs1(row[k]);
            if c > best {
                best = c;
                p = i;
            }
        }
        piv[k] = p;
        if p != k {
            a.swap(p, k);
        }
        let d = a[k][k];
        if cabs1(d) != 0.0 {
            // dgetrf2: DSCAL by ONE / A(1,1) when ABS(A(1,1)) >= SFMIN.
            if d.abs() >= f64::MIN_POSITIVE {
                let rcp = T::one().fortran_div(d);
                for i in k + 1..m {
                    a[i][k] = a[i][k] * rcp;
                }
            } else {
                for i in k + 1..m {
                    a[i][k] = a[i][k].fortran_div(d);
                }
            }
        }
        for i in k + 1..m {
            let f = a[i][k];
            for j in k + 1..m {
                let v = a[k][j];
                a[i][j] = a[i][j] - f * v;
            }
        }
    }
    piv
}

/// Solves `A x = rhs` in place given [`lu_factor`] output.
pub(crate) fn lu_solve<T: IterScalar>(a: &[Vec<T>], piv: &[usize], rhs: &mut [T]) {
    let m = a.len();
    for k in 0..m {
        rhs.swap(k, piv[k]);
    }
    // DTRSM('L','L','N','U'): column k is subtracted from the rows below.
    for k in 0..m {
        let bk = rhs[k];
        if bk != T::zero() {
            for i in k + 1..m {
                rhs[i] = rhs[i] - bk * a[i][k];
            }
        }
    }
    // DTRSM('L','U','N','N'): from the last column back.
    for k in (0..m).rev() {
        if rhs[k] != T::zero() {
            rhs[k] = rhs[k].fortran_div(a[k][k]);
            let bk = rhs[k];
            for i in 0..k {
                rhs[i] = rhs[i] - bk * a[i][k];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn solves_with_pivoting() {
        // Needs a row swap: leading entry is zero.
        let mut a = vec![
            vec![0.0, 2.0, 1.0],
            vec![1.0, 1.0, 0.0],
            vec![3.0, 0.0, 1.0],
        ];
        let a0 = a.clone();
        let piv = lu_factor(&mut a);
        let mut x = vec![1.0, 2.0, 3.0];
        let rhs = x.clone();
        lu_solve(&a, &piv, &mut x);
        for i in 0..3 {
            let s: f64 = (0..3).map(|j| a0[i][j] * x[j]).sum();
            assert!((s - rhs[i]).abs() < 1e-12, "row {i}: {s} vs {}", rhs[i]);
        }
    }

    /// Reference LAPACK's rounding: `ONE / pivot` scaling and the
    /// column-oriented upper solve, written out for a 3x3 without pivoting.
    #[test]
    fn rounds_as_reference_lapack() {
        let a0: Vec<Vec<f64>> = vec![
            vec![3.9, -0.8, -0.6],
            vec![0.3, 0.9, 0.2],
            vec![-0.2, 1.0, -0.9],
        ];
        // Division and row-oriented substitution give different last bits.
        let b0 = [-1.26, -2.13, -2.29];
        // dgetrf2 / DTRSM by hand.
        let r = 1.0 / a0[0][0];
        let (l10, l20) = (a0[1][0] * r, a0[2][0] * r);
        let u11 = a0[1][1] - l10 * a0[0][1];
        let u12 = a0[1][2] - l10 * a0[0][2];
        let a21 = a0[2][1] - l20 * a0[0][1];
        let u22i = a0[2][2] - l20 * a0[0][2];
        assert!(u11.abs() >= a21.abs(), "no second pivot swap");
        let l21 = a21 * (1.0 / u11);
        let u22 = u22i - l21 * u12;
        let y1 = b0[1] - b0[0] * l10;
        let y2 = (b0[2] - b0[0] * l20) - y1 * l21;
        let x2 = y2 / u22;
        let x1 = (y1 - x2 * u12) / u11;
        let x0 = ((b0[0] - x2 * a0[0][2]) - x1 * a0[0][1]) / a0[0][0];

        let mut a = a0.clone();
        let piv = lu_factor(&mut a);
        assert_eq!(piv, vec![0, 1, 2]);
        let mut x = b0.to_vec();
        lu_solve(&a, &piv, &mut x);
        assert_eq!(x, vec![x0, x1, x2]);
    }

    #[test]
    fn solves_complex() {
        let c = |r, i| Complex64::new(r, i);
        let mut a = vec![
            vec![c(1.0, 1.0), c(2.0, 0.0)],
            vec![c(0.0, 1.0), c(3.0, -1.0)],
        ];
        let a0 = a.clone();
        let piv = lu_factor(&mut a);
        let rhs = vec![c(1.0, 0.0), c(0.0, 2.0)];
        let mut x = rhs.clone();
        lu_solve(&a, &piv, &mut x);
        for i in 0..2 {
            let s = a0[i][0] * x[0] + a0[i][1] * x[1];
            assert!((s - rhs[i]).norm() < 1e-12);
        }
    }
}
