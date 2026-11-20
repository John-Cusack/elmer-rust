//! The reference BLAS routines (Netlib, the `libblas.so.3` the pristine
//! Elmer build links) under ARPACK's complex driver (`znaupd` / `zneupd`),
//! in their operation order: unit strides only, column-major matrices with
//! leading dimension `ld`.
//!
//! Complex products are the plain `(ac - bd, ad + bc)` gfortran emits under
//! its Fortran rules, as `num_complex` computes them. Complex quotients are
//! not: gfortran divides with Smith's range reduction ([`cdiv`]), which
//! `num_complex`'s `Div` does not, so every Fortran `/` of two complex
//! numbers goes through [`cdiv`].

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use num_complex::Complex64 as C;

pub const ZERO: C = C::new(0.0, 0.0);
pub const ONE: C = C::new(1.0, 0.0);

/// A complex quotient as gfortran computes it (`expand_complex_div_wide`,
/// Smith's method without the NaN recovery).
#[inline]
pub fn cdiv(a: C, b: C) -> C {
    let (ar, ai, br, bi) = (a.re, a.im, b.re, b.im);
    if br.abs() < bi.abs() {
        let ratio = br / bi;
        let div = br * ratio + bi;
        C::new((ar * ratio + ai) / div, (ai * ratio - ar) / div)
    } else {
        let ratio = bi / br;
        let div = bi * ratio + br;
        C::new((ai * ratio + ar) / div, (ai - ar * ratio) / div)
    }
}

/// A real divided by a complex number: gfortran promotes the real.
#[inline]
pub fn rdiv(a: f64, b: C) -> C {
    cdiv(C::new(a, 0.0), b)
}

/// `ABS` of a complex number (`cabs`, the libm `hypot`).
#[inline]
pub fn cabs(z: C) -> f64 {
    z.re.hypot(z.im)
}

/// `DCABS1` / `CABS1`: `|re| + |im|`.
#[inline]
pub fn cabs1(z: C) -> f64 {
    z.re.abs() + z.im.abs()
}

/// `ZCOPY(n, x, 1, y, 1)`.
pub fn zcopy(n: usize, x: &[C], y: &mut [C]) {
    y[..n].copy_from_slice(&x[..n]);
}

/// `ZSWAP(n, x, 1, y, 1)`.
pub fn zswap(n: usize, x: &mut [C], y: &mut [C]) {
    for i in 0..n {
        std::mem::swap(&mut x[i], &mut y[i]);
    }
}

/// `ZSCAL(n, za, x, 1)` (3.12 returns at once when `za` is one, keeping
/// the signs of zeros).
pub fn zscal(n: usize, za: C, x: &mut [C]) {
    if za == ONE {
        return;
    }
    for v in &mut x[..n] {
        *v = za * *v;
    }
}

/// `ZDSCAL(n, da, x, 1)` (3.12: both parts scaled by the real `da`, and
/// nothing done when `da` is one).
pub fn zdscal(n: usize, da: f64, x: &mut [C]) {
    if da == 1.0 {
        return;
    }
    for v in &mut x[..n] {
        *v = C::new(da * v.re, da * v.im);
    }
}

/// `ZAXPY(n, za, x, 1, y, 1)`: `y = y + za x`.
pub fn zaxpy(n: usize, za: C, x: &[C], y: &mut [C]) {
    if n == 0 || cabs1(za) == 0.0 {
        return;
    }
    for i in 0..n {
        y[i] += za * x[i];
    }
}

/// `ZDOTC(n, x, 1, y, 1)`: `sum conj(x_i) y_i`.
pub fn zdotc(n: usize, x: &[C], y: &[C]) -> C {
    let mut t = ZERO;
    for i in 0..n {
        t += x[i].conj() * y[i];
    }
    t
}

/// `ZDOTU(n, x, 1, y, 1)`: `sum x_i y_i`.
pub fn zdotu(n: usize, x: &[C], y: &[C]) -> C {
    let mut t = ZERO;
    for i in 0..n {
        t += x[i] * y[i];
    }
    t
}

/// `DZNRM2(n, x, 1)`: Blue's algorithm (dznrm2.f90, reference BLAS since
/// 3.10), the real and the imaginary part of each entry in turn.
pub fn dznrm2(n: usize, x: &[C]) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let tsml = 2f64.powi(-511);
    let tbig = 2f64.powi(486);
    let ssml = 2f64.powi(537);
    let sbig = 2f64.powi(-538);
    let maxn = f64::MAX;
    let mut notbig = true;
    let (mut asml, mut amed, mut abig) = (0.0f64, 0.0f64, 0.0f64);
    for z in &x[..n] {
        for ax in [z.re.abs(), z.im.abs()] {
            if ax > tbig {
                abig += (ax * sbig) * (ax * sbig);
                notbig = false;
            } else if ax < tsml {
                if notbig {
                    asml += (ax * ssml) * (ax * ssml);
                }
            } else {
                amed += ax * ax;
            }
        }
    }
    let (scl, sumsq);
    if abig > 0.0 {
        if amed > 0.0 || amed > maxn || amed.is_nan() {
            abig += (amed * sbig) * sbig;
        }
        scl = 1.0 / sbig;
        sumsq = abig;
    } else if asml > 0.0 {
        if amed > 0.0 || amed > maxn || amed.is_nan() {
            let amed = amed.sqrt();
            let asml = asml.sqrt() / ssml;
            let (ymin, ymax) = if asml > amed { (amed, asml) } else { (asml, amed) };
            scl = 1.0;
            sumsq = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax));
        } else {
            scl = 1.0 / ssml;
            sumsq = asml;
        }
    } else {
        scl = 1.0;
        sumsq = amed;
    }
    scl * sumsq.sqrt()
}

/// `IZAMAX(n, x, 1)`: 1-based index of the first largest `|re| + |im|`.
pub fn izamax(n: usize, x: &[C]) -> usize {
    if n == 0 {
        return 0;
    }
    let mut imax = 1;
    let mut dmax = cabs1(x[0]);
    for i in 1..n {
        let v = cabs1(x[i]);
        if v > dmax {
            imax = i + 1;
            dmax = v;
        }
    }
    imax
}

/// `DZASUM(n, x, 1)`.
pub fn dzasum(n: usize, x: &[C]) -> f64 {
    let mut s = 0.0;
    for z in &x[..n] {
        s += cabs1(*z);
    }
    s
}

fn scale_y(len: usize, beta: C, y: &mut [C]) {
    if beta != ONE {
        if beta == ZERO {
            y[..len].iter_mut().for_each(|v| *v = ZERO);
        } else {
            y[..len].iter_mut().for_each(|v| *v = beta * *v);
        }
    }
}

/// `ZGEMV('N', m, n, alpha, a, lda, x, 1, beta, y, 1)`: `y = beta y + alpha A x`.
#[allow(clippy::too_many_arguments)]
pub fn zgemv_n(m: usize, n: usize, alpha: C, a: &[C], lda: usize, x: &[C], beta: C, y: &mut [C]) {
    if m == 0 || n == 0 || (alpha == ZERO && beta == ONE) {
        return;
    }
    scale_y(m, beta, y);
    if alpha == ZERO {
        return;
    }
    for j in 0..n {
        let temp = alpha * x[j];
        let col = &a[j * lda..j * lda + m];
        for i in 0..m {
            y[i] += temp * col[i];
        }
    }
}

/// `ZGEMV('C', m, n, alpha, a, lda, x, 1, beta, y, 1)`: `y = beta y + alpha A^H x`.
#[allow(clippy::too_many_arguments)]
pub fn zgemv_c(m: usize, n: usize, alpha: C, a: &[C], lda: usize, x: &[C], beta: C, y: &mut [C]) {
    if m == 0 || n == 0 || (alpha == ZERO && beta == ONE) {
        return;
    }
    scale_y(n, beta, y);
    if alpha == ZERO {
        return;
    }
    for j in 0..n {
        let col = &a[j * lda..j * lda + m];
        let mut temp = ZERO;
        for i in 0..m {
            temp += col[i].conj() * x[i];
        }
        y[j] += alpha * temp;
    }
}

/// `ZGERU(m, n, alpha, x, 1, y, 1, a, lda)`: `A = A + alpha x y^T`.
#[allow(clippy::too_many_arguments)]
pub fn zgeru(m: usize, n: usize, alpha: C, x: &[C], y: &[C], a: &mut [C], lda: usize) {
    if m == 0 || n == 0 || alpha == ZERO {
        return;
    }
    for j in 0..n {
        let temp = alpha * y[j];
        let col = &mut a[j * lda..j * lda + m];
        for i in 0..m {
            col[i] += x[i] * temp;
        }
    }
}

/// `ZTRMM('Right', 'Upper', 'No transpose', 'Non-unit', m, n, alpha, a,
/// lda, b, ldb)`: `B = alpha B A` with `A` upper triangular (`n x n`).
#[allow(clippy::too_many_arguments)]
pub fn ztrmm_runn(m: usize, n: usize, alpha: C, a: &[C], lda: usize, b: &mut [C], ldb: usize) {
    if m == 0 || n == 0 {
        return;
    }
    if alpha == ZERO {
        for j in 0..n {
            for i in 0..m {
                b[i + j * ldb] = ZERO;
            }
        }
        return;
    }
    for j in (0..n).rev() {
        let temp = alpha * a[j + j * lda];
        for i in 0..m {
            b[i + j * ldb] = temp * b[i + j * ldb];
        }
        for k in 0..j {
            let akj = a[k + j * lda];
            if akj != ZERO {
                let temp = alpha * akj;
                for i in 0..m {
                    let t = b[i + k * ldb];
                    b[i + j * ldb] += temp * t;
                }
            }
        }
    }
}
