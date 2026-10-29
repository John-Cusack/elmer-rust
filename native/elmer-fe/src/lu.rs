//! `LUDecomp` and `InvertMatrix` (LinearAlgebra.F90:56,191), the dense
//! inverse Elmer uses to turn the elements.def node/term matrix into basis
//! coefficients. Column-major `a[i + j*n]` like the Fortran `A(i,j)`.

/// LUDecomp: in-place LU with column pivoting; the diagonal is replaced by
/// its reciprocal. Returns `None` when a pivot is exactly zero.
fn lu_decomp(a: &mut [f64], n: usize, pivot: &mut [usize]) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[ix(i, k)].abs() > a[ix(i, j)].abs() {
                j = k;
            }
        }
        if a[ix(i, j)].abs() == 0.0 {
            return None;
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a.swap(ix(k, j), ix(k, i));
            }
        }
        for k in i + 1..n {
            a[ix(i, k)] /= a[ix(i, i)];
        }
        for k in i + 1..n {
            if j != i {
                a.swap(ix(k, i), ix(k, j));
            }
            for l in i + 1..n {
                a[ix(k, l)] -= a[ix(k, i)] * a[ix(i, l)];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[ix(i, i)].abs() == 0.0 {
            return None;
        }
        a[ix(i, i)] = 1.0 / a[ix(i, i)];
    }
    Some(())
}

/// InvertMatrix: `a` (n x n, column-major) is replaced by its inverse.
pub(crate) fn invert_matrix(a: &mut [f64], n: usize) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    let mut pivot = vec![0usize; n];
    lu_decomp(a, n, &mut pivot)?;

    // INV(U)
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = -a[ix(i, j)];
            for k in i + 1..j {
                s -= a[ix(i, k)] * a[ix(k, j)];
            }
            a[ix(i, j)] = s;
        }
    }

    // INV(L)
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s -= a[ix(j, k)] * a[ix(k, i)];
            }
            a[ix(j, i)] = a[ix(i, i)] * s;
        }
    }

    // A = INV(AP)
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                if k != i {
                    s += a[ix(i, k)] * a[ix(k, j)];
                } else {
                    s += a[ix(k, j)];
                }
            }
            a[ix(i, j)] = s;
        }
    }

    // A = INV(A) (at last)
    for i in (0..n).rev() {
        if pivot[i] != i {
            for j in 0..n {
                a.swap(ix(i, j), ix(pivot[i], j));
            }
        }
    }
    Some(())
}
