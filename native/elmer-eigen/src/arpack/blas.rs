//! The reference BLAS routines ARPACK's symmetric driver calls, in the
//! operation order of Netlib's reference implementation (the `libblas.so.3`
//! the pristine Elmer build links): unit strides only, which is all ARPACK
//! uses. Matrices are column-major with leading dimension `ld`.

/// `DCOPY(n, x, 1, y, 1)`.
pub fn dcopy(n: usize, x: &[f64], y: &mut [f64]) {
    y[..n].copy_from_slice(&x[..n]);
}

/// `DSWAP(n, x, 1, y, 1)`.
pub fn dswap(n: usize, x: &mut [f64], y: &mut [f64]) {
    for i in 0..n {
        std::mem::swap(&mut x[i], &mut y[i]);
    }
}

/// `DSWAP` of two ranges of one array (`x(a..)` and `x(b..)`, may not overlap).
pub fn dswap_within(n: usize, x: &mut [f64], a: usize, b: usize) {
    for i in 0..n {
        x.swap(a + i, b + i);
    }
}

/// `DSCAL(n, a, x, 1)`.
pub fn dscal(n: usize, a: f64, x: &mut [f64]) {
    for v in &mut x[..n] {
        *v *= a;
    }
}

/// `DAXPY(n, a, x, 1, y, 1)`: `y = y + a x`; nothing when `a == 0`.
pub fn daxpy(n: usize, a: f64, x: &[f64], y: &mut [f64]) {
    if n == 0 || a == 0.0 {
        return;
    }
    for i in 0..n {
        y[i] += a * x[i];
    }
}

/// `DDOT(n, x, 1, y, 1)`: one accumulator, left to right (the unrolled
/// reference loop adds in the same order).
pub fn ddot(n: usize, x: &[f64], y: &[f64]) -> f64 {
    let mut t = 0.0;
    for i in 0..n {
        t += x[i] * y[i];
    }
    t
}

/// `DNRM2(n, x, 1)`: Blue's algorithm with three accumulators, as the
/// reference BLAS has it since 3.10 (dnrm2.f90, Anderson 2017).
pub fn dnrm2(n: usize, x: &[f64]) -> f64 {
    if n == 0 {
        return 0.0;
    }
    // radix 2, minexponent -1021, maxexponent 1024, digits 53.
    let tsml = 2f64.powi(-511); // 2**ceiling((-1021 - 1) * 0.5)
    let tbig = 2f64.powi(486); // 2**floor((1024 - 53 + 1) * 0.5)
    let ssml = 2f64.powi(537); // 2**(-floor((-1021 - 53) * 0.5))
    let sbig = 2f64.powi(-538); // 2**(-ceiling((1024 + 53 - 1) * 0.5))
    let maxn = f64::MAX;
    let mut notbig = true;
    let (mut asml, mut amed, mut abig) = (0.0f64, 0.0f64, 0.0f64);
    for &xi in &x[..n] {
        let ax = xi.abs();
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

/// `DGEMV('N', m, n, alpha, a, lda, x, 1, beta, y, 1)`: `y = beta y + alpha A x`.
#[allow(clippy::too_many_arguments)]
pub fn dgemv_n(m: usize, n: usize, alpha: f64, a: &[f64], lda: usize, x: &[f64], beta: f64, y: &mut [f64]) {
    if m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0) {
        return;
    }
    scale_y(m, beta, y);
    if alpha == 0.0 {
        return;
    }
    for j in 0..n {
        if x[j] != 0.0 {
            let temp = alpha * x[j];
            let col = &a[j * lda..j * lda + m];
            for i in 0..m {
                y[i] += temp * col[i];
            }
        }
    }
}

/// `DGEMV('T', m, n, alpha, a, lda, x, 1, beta, y, 1)`: `y = beta y + alpha A^T x`.
#[allow(clippy::too_many_arguments)]
pub fn dgemv_t(m: usize, n: usize, alpha: f64, a: &[f64], lda: usize, x: &[f64], beta: f64, y: &mut [f64]) {
    if m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0) {
        return;
    }
    scale_y(n, beta, y);
    if alpha == 0.0 {
        return;
    }
    for j in 0..n {
        let col = &a[j * lda..j * lda + m];
        let mut temp = 0.0;
        for i in 0..m {
            temp += col[i] * x[i];
        }
        y[j] += alpha * temp;
    }
}

fn scale_y(len: usize, beta: f64, y: &mut [f64]) {
    if beta != 1.0 {
        if beta == 0.0 {
            y[..len].iter_mut().for_each(|v| *v = 0.0);
        } else {
            y[..len].iter_mut().for_each(|v| *v *= beta);
        }
    }
}

/// `DGER(m, n, alpha, x, 1, y, 1, a, lda)`: `A = A + alpha x y^T`.
#[allow(clippy::too_many_arguments)]
pub fn dger(m: usize, n: usize, alpha: f64, x: &[f64], y: &[f64], a: &mut [f64], lda: usize) {
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    for j in 0..n {
        if y[j] != 0.0 {
            let temp = alpha * y[j];
            let col = &mut a[j * lda..j * lda + m];
            for i in 0..m {
                col[i] += x[i] * temp;
            }
        }
    }
}
