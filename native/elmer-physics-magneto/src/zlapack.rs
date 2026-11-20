//! Reference LAPACK 3.12 / BLAS `ZGETRF` + `ZGETRS` (the `liblapack.so.3`
//! and `libblas.so.3` Elmer links), loop for loop, with gfortran's complex
//! arithmetic: products `(ac - bd, ad + bc)` and the range-reduced (Smith)
//! division of `-fcx-fortran-rules`. Column-major, `lda == n`, 0-based
//! pivots in LAPACK's sequential interchange form.

// Index loops mirror the reference Fortran.
#![allow(clippy::needless_range_loop)]

use elmer_itermethods::IterScalar;
use num_complex::Complex64 as C;

const ONE: C = C::new(1.0, 0.0);
const ZERO: C = C::new(0.0, 0.0);
/// `ILAENV(1, 'ZGETRF', ...)`.
const NB: usize = 64;

#[inline]
fn mul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// `DCABS1`.
fn cabs1(z: C) -> f64 {
    z.re.abs() + z.im.abs()
}

/// A column-major view `a[(r0 + i) + (c0 + j) * lda]`.
#[derive(Clone, Copy)]
struct At {
    r0: usize,
    c0: usize,
}

impl At {
    fn idx(self, lda: usize, i: usize, j: usize) -> usize {
        (self.r0 + i) + (self.c0 + j) * lda
    }
    fn sub(self, i: usize, j: usize) -> At {
        At { r0: self.r0 + i, c0: self.c0 + j }
    }
}

/// `ZLASWP(N, A, LDA, K1, K2, IPIV, 1)` on the `n` columns at `at`, rows
/// `k1..k2` (0-based, inclusive-exclusive) interchanged with `ipiv`.
fn zlaswp(a: &mut [C], lda: usize, at: At, n: usize, k1: usize, k2: usize, ipiv: &[usize]) {
    for i in k1..k2 {
        let ip = ipiv[i];
        if ip != i {
            for j in 0..n {
                a.swap(at.idx(lda, i, j), at.idx(lda, ip, j));
            }
        }
    }
}

/// `ZTRSM('L', 'L', 'N', 'U', M, N, ONE, A, LDA, B, LDB)`.
fn ztrsm_llnu(a: &mut [C], lda: usize, ta: At, m: usize, n: usize, tb: At) {
    for j in 0..n {
        for k in 0..m {
            let bkj = a[tb.idx(lda, k, j)];
            if bkj != ZERO {
                for i in k + 1..m {
                    let p = mul(bkj, a[ta.idx(lda, i, k)]);
                    let b = &mut a[tb.idx(lda, i, j)];
                    *b -= p;
                }
            }
        }
    }
}

/// `ZTRSM('L', 'U', 'N', 'N', M, N, ONE, A, LDA, B, LDB)` with `B` a
/// separate array.
fn ztrsm_lunn(a: &[C], lda: usize, m: usize, b: &mut [C]) {
    for k in (0..m).rev() {
        if b[k] != ZERO {
            b[k] = b[k].fortran_div(a[k + k * lda]);
            let bk = b[k];
            for i in 0..k {
                b[i] -= mul(bk, a[i + k * lda]);
            }
        }
    }
}

/// `ZGEMM('N', 'N', M, N, K, -ONE, A, LDA, B, LDB, ONE, C, LDC)`, all three
/// in `a`.
#[allow(clippy::too_many_arguments)]
fn zgemm_nn_minus(a: &mut [C], lda: usize, m: usize, n: usize, k: usize, ta: At, tb: At, tc: At) {
    let alpha = C::new(-1.0, 0.0);
    for j in 0..n {
        for l in 0..k {
            let temp = mul(alpha, a[tb.idx(lda, l, j)]);
            for i in 0..m {
                let p = mul(temp, a[ta.idx(lda, i, l)]);
                let c = &mut a[tc.idx(lda, i, j)];
                *c += p;
            }
        }
    }
}

/// `ZGETRF2` (recursive) on the `m x n` block at `at`; `ipiv[..min(m,n)]`
/// relative to the block.
fn zgetrf2(a: &mut [C], lda: usize, at: At, m: usize, n: usize, ipiv: &mut [usize]) -> Result<(), String> {
    if m == 0 || n == 0 {
        return Ok(());
    }
    if m == 1 {
        ipiv[0] = 0;
        if a[at.idx(lda, 0, 0)] == ZERO {
            return Err("ZGETRF2: singular".into());
        }
    } else if n == 1 {
        let sfmin = f64::MIN_POSITIVE;
        // IZAMAX.
        let mut imax = 0;
        let mut dmax = cabs1(a[at.idx(lda, 0, 0)]);
        for i in 1..m {
            let d = cabs1(a[at.idx(lda, i, 0)]);
            if d > dmax {
                imax = i;
                dmax = d;
            }
        }
        ipiv[0] = imax;
        if a[at.idx(lda, imax, 0)] != ZERO {
            if imax != 0 {
                a.swap(at.idx(lda, 0, 0), at.idx(lda, imax, 0));
            }
            let a11 = a[at.idx(lda, 0, 0)];
            if a11.norm() >= sfmin {
                let za = ONE.fortran_div(a11);
                // ZSCAL (returns at once for ZA = ONE).
                if za != ONE {
                    for i in 1..m {
                        let x = &mut a[at.idx(lda, i, 0)];
                        *x = mul(za, *x);
                    }
                }
            } else {
                for i in 1..m {
                    let x = &mut a[at.idx(lda, i, 0)];
                    *x = x.fortran_div(a11);
                }
            }
        } else {
            return Err("ZGETRF2: singular".into());
        }
    } else {
        let n1 = m.min(n) / 2;
        let n2 = n - n1;
        zgetrf2(a, lda, at, m, n1, &mut ipiv[..n1])?;
        zlaswp(a, lda, at.sub(0, n1), n2, 0, n1, ipiv);
        ztrsm_llnu(a, lda, at, n1, n2, at.sub(0, n1));
        zgemm_nn_minus(a, lda, m - n1, n2, n1, at.sub(n1, 0), at.sub(0, n1), at.sub(n1, n1));
        let k = m.min(n);
        zgetrf2(a, lda, at.sub(n1, n1), m - n1, n2, &mut ipiv[n1..k])?;
        for p in &mut ipiv[n1..k] {
            *p += n1;
        }
        zlaswp(a, lda, at, n1, n1, k, ipiv);
    }
    Ok(())
}

/// `ZGETRF(N, N, A, N, IPIV, INFO)`.
pub fn zgetrf(n: usize, a: &mut [C]) -> Result<Vec<usize>, String> {
    let mut ipiv = vec![0usize; n];
    let lda = n;
    let base = At { r0: 0, c0: 0 };
    if NB >= n {
        zgetrf2(a, lda, base, n, n, &mut ipiv)?;
        return Ok(ipiv);
    }
    let mut j = 0;
    while j < n {
        let jb = (n - j).min(NB);
        zgetrf2(a, lda, base.sub(j, j), n - j, jb, &mut ipiv[j..j + jb])?;
        for p in &mut ipiv[j..j + jb] {
            *p += j;
        }
        zlaswp(a, lda, base, j, j, j + jb, &ipiv);
        if j + jb < n {
            zlaswp(a, lda, base.sub(0, j + jb), n - j - jb, j, j + jb, &ipiv);
            ztrsm_llnu(a, lda, base.sub(j, j), jb, n - j - jb, base.sub(j, j + jb));
            zgemm_nn_minus(a, lda, n - j - jb, n - j - jb, jb, base.sub(j + jb, j), base.sub(j, j + jb), base.sub(j + jb, j + jb));
        }
        j += NB;
    }
    Ok(ipiv)
}

/// `ZGETRS('N', N, 1, A, N, IPIV, B, N, INFO)`.
pub fn zgetrs(n: usize, lu: &[C], ipiv: &[usize], b: &mut [C]) {
    for i in 0..n {
        let ip = ipiv[i];
        if ip != i {
            b.swap(i, ip);
        }
    }
    // ZTRSM('L', 'L', 'N', 'U') on the one column.
    for k in 0..n {
        let bk = b[k];
        if bk != ZERO {
            for i in k + 1..n {
                b[i] -= mul(bk, lu[i + k * n]);
            }
        }
    }
    ztrsm_lunn(lu, n, n, b);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_blocked_system() {
        let n = 70;
        let mut a = vec![ZERO; n * n];
        for j in 0..n {
            for i in 0..n {
                let x = ((i * 7 + j * 13) % 17) as f64 - 8.0;
                a[i + j * n] = C::new(x, ((i + 2 * j) % 5) as f64 - 2.0);
            }
            a[j + j * n] += C::new(40.0, 3.0);
        }
        let xs: Vec<C> = (0..n).map(|i| C::new(i as f64, 1.0 - i as f64)).collect();
        let mut b = vec![ZERO; n];
        for i in 0..n {
            for j in 0..n {
                b[i] += a[i + j * n] * xs[j];
            }
        }
        let mut lu = a.clone();
        let ipiv = zgetrf(n, &mut lu).unwrap();
        zgetrs(n, &lu, &ipiv, &mut b);
        for (x, want) in b.iter().zip(&xs) {
            assert!((x - want).norm() < 1e-9, "{x} {want}");
        }
    }
}
