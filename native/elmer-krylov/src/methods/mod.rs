//! One module per HUTI method. Each is a statement-by-statement port of the
//! double (`d`) and double complex (`z`) routine in `fhutiter/src/huti_*.F90`:
//! array expressions become left-to-right element loops in the same operation
//! order, so the arithmetic matches the reference build bit for bit. That is
//! why the index loops stay index loops (`clippy::needless_range_loop`).
#![allow(clippy::needless_range_loop)]

pub(crate) mod bicgstab;
pub(crate) mod bicgstab2;
pub(crate) mod cg;
pub(crate) mod cgs;
pub(crate) mod gmres;
pub(crate) mod qmr;
pub(crate) mod tfqmr;

use crate::config::HutiStop;

/// The criteria that need `||b||` up front (`rhsnorm`).
pub(crate) fn needs_rhsnorm(stop: HutiStop) -> bool {
    matches!(
        stop,
        HutiStop::TrueResidualScaledByB | HutiStop::PseudoResidualScaledByB
    )
}

/// `HUTI_EPSILON` (huti_fdefs.h:16), a default-real literal compared in double.
pub(crate) const HUTI_EPSILON: f64 = f32::MIN_POSITIVE as f64;

/// Convergence then divergence, as every method checks them after printing.
pub(crate) fn converged(residual: f64, tol: f64) -> bool {
    residual < tol
}

pub(crate) fn diverged(residual: f64, max_tol: f64) -> bool {
    residual.is_nan() || residual > max_tol
}
