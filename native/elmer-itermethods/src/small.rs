//! Tiny dense helpers for the BiCGStab(l) polynomial step: LU with partial
//! pivoting (`dgetrf`/`dgetrs`, `zgetrf`/`zgetrs`) on matrices of order l-1.

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
            for i in k + 1..m {
                a[i][k] = a[i][k] / d;
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
    for i in 0..m {
        for j in 0..i {
            rhs[i] = rhs[i] - a[i][j] * rhs[j];
        }
    }
    for i in (0..m).rev() {
        for j in i + 1..m {
            rhs[i] = rhs[i] - a[i][j] * rhs[j];
        }
        rhs[i] = rhs[i] / a[i][i];
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
