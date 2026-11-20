//! What PARPACK (ARPACK's MPI version, `p`-prefixed routines) adds to the
//! serial routines: each rank holds its rows of the Arnoldi vectors, and the
//! inner products, norms and `V^T x` products are reduced over the ranks
//! (`MPI_ALLREDUCE`, one call per reduction, as the Fortran makes them).
//!
//! The serial routines take a [`Par`] and run unchanged when it is `None`.

use elmer_core::Communicator;

use super::blas::{ddot, dnrm2};
use super::lapack::EPS;

/// The communicator of a PARPACK call, or `None` for serial ARPACK.
pub type Par<'a> = Option<&'a dyn Communicator>;

/// `ddot`, then `MPI_ALLREDUCE(SUM)` in PARPACK.
pub(crate) fn pddot(par: Par<'_>, n: usize, x: &[f64], y: &[f64]) -> f64 {
    let d = ddot(n, x, y);
    match par {
        None => d,
        Some(c) => {
            let mut b = [d];
            c.sum_f64(&mut b);
            b[0]
        }
    }
}

/// `dnrm2`, or PARPACK's `pdnorm2`: the local norms' maximum `m`, then
/// `m * sqrt(|sum over ranks of (local/m)**2|)`.
pub(crate) fn pdnorm2(par: Par<'_>, n: usize, x: &[f64]) -> f64 {
    let local = dnrm2(n, x);
    let Some(c) = par else { return local };
    let mut max = [local];
    c.max_f64(&mut max);
    if max[0] == 0.0 {
        return 0.0;
    }
    // `(pdnorm2/max)**2.0`, which gfortran evaluates as a product.
    let r = local / max[0];
    let mut s = [r * r];
    c.sum_f64(&mut s);
    max[0] * s[0].abs().sqrt()
}

/// `MPI_ALLREDUCE(SUM)` of `buf` in place under PARPACK; nothing serially.
pub(crate) fn psum(par: Par<'_>, buf: &mut [f64]) {
    if let Some(c) = par {
        c.sum_f64(buf);
    }
}

/// `eps23` of `_saup2`, `_naup2` and `_seupd`: `eps**(2/3)`, the exponent
/// a double-precision literal in ARPACK and a default-REAL one, `2.0/3.0`,
/// in PARPACK.
pub(crate) fn eps23(par: Par<'_>) -> f64 {
    if par.is_some() {
        EPS.powf(f64::from(2.0f32 / 3.0f32))
    } else {
        EPS.powf(2.0 / 3.0)
    }
}
