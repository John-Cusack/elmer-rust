//! Vector kernels shared by the methods. Inner products reduce through the
//! [`Communicator`], so a distributed run only has to supply `sum_f64`.

use crate::scalar::IterScalar;
use elmer_core::{Communicator, LinalgError, LinearOperator, Preconditioner, Scalar};

/// `x^H y` (`ddot` for real, `zdotc` for complex).
pub(crate) fn dot<T: IterScalar>(comm: &dyn Communicator, x: &[T], y: &[T]) -> T {
    let mut s = T::zero();
    for (a, b) in x.iter().zip(y) {
        s = s + a.conj() * *b;
    }
    let mut buf = [s.re(), s.im()];
    comm.sum_f64(&mut buf);
    T::from_parts(buf[0], buf[1])
}

/// Euclidean norm (`dnrm2` / `dznrm2`).
pub(crate) fn nrm2<T: IterScalar>(comm: &dyn Communicator, x: &[T]) -> f64 {
    let mut s = 0.0;
    for v in x {
        let a = v.abs();
        s += a * a;
    }
    let mut buf = [s];
    comm.sum_f64(&mut buf);
    buf[0].sqrt()
}

/// `y += a x`.
pub(crate) fn axpy<T: Scalar>(a: T, x: &[T], y: &mut [T]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi = *yi + a * *xi;
    }
}

/// `y -= a x`.
pub(crate) fn axmy<T: Scalar>(a: T, x: &[T], y: &mut [T]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi = *yi - a * *xi;
    }
}

/// `r = b - A x`, using `r` as the matvec output.
pub(crate) fn residual<T: Scalar>(
    a: &dyn LinearOperator<T>,
    x: &[T],
    b: &[T],
    r: &mut [T],
) {
    a.apply(x, r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri = *bi - *ri;
    }
}

pub(crate) fn check_dims<T: Scalar>(
    method: &str,
    a: &dyn LinearOperator<T>,
    b: &[T],
    x: &[T],
) -> Result<(), LinalgError> {
    let n = b.len();
    if x.len() != n || a.nrows() != n || a.ncols() != n {
        return Err(LinalgError::Dim(format!(
            "{method}: operator is {}x{}, b has {}, x has {}",
            a.nrows(),
            a.ncols(),
            n,
            x.len()
        )));
    }
    Ok(())
}

/// `z = M^-1 r` (the right preconditioner `C_rpcond`).
pub(crate) fn pcond<T: Scalar>(m: &dyn Preconditioner<T>, r: &[T], z: &mut [T]) {
    m.apply(r, z);
}

/// Elmer's divergence test: over the limit or NaN.
pub(crate) fn diverged(res: f64, max_tol: f64) -> bool {
    res > max_tol || res.is_nan()
}
