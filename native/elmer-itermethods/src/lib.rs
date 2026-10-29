//! Elmer's `IterativeMethods.F90` family: BiCGStab(l), GCR and IDR(s) for real
//! and complex systems, plus Jacobi, Richardson and symmetric Gauss-Seidel.
//!
//! Reference: `fem/src/IterativeMethods.F90` (Jacobi `:296`, Richardson `:404`,
//! BiCGStabl `:577/693`, GCR `:1188/1271`, IDRS `:1481/1590`, complex GCR,
//! BiCGStabl and IDRS `:2456-3371`) and the keyword handling and dispatch in
//! `fem/src/IterSolve.F90`.
//!
//! All methods apply the preconditioner on the right and stop on `||r|| / ||b||`
//! (or a caller-supplied criterion, see `with_stop_criterion`), with the
//! initial-guess, divergence and breakdown rules of the Fortran. The iteration
//! count in [`SolveStats`] is the number Elmer prints (`k` for GCR, `iter` for
//! IDR(s), `min(MaxRounds, Round)` for BiCGStab(l)).
//!
//! Not covered: the constraint-matrix extension of the real methods (the caller
//! passes a block operator and preconditioner over the extended unknowns), MPRGP,
//! and the HUTI methods (`elmer-krylov`).
#![forbid(unsafe_code)]
// Index loops mirror the Fortran, and `!(x > 0.0)` is the NaN-catching test it
// uses on purpose.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

mod bicgstabl;
mod config;
mod gcr;
mod idrs;
mod ops;
mod rng;
mod scalar;
mod small;
mod stationary;

pub use bicgstabl::BiCgStabL;
pub use config::{IterMethod, IterMethodsConfig, PseudoComplex, RobustConfig};
pub use gcr::Gcr;
pub use idrs::Idrs;
pub use scalar::IterScalar;
pub use stationary::{Jacobi, Richardson, Sgs};

use elmer_core::{SolveStats, SolveStatus};

/// Replaces the default `||r|| / ||b||` (the `stopcfun` of HUTI): called as
/// `f(x, b, r)` and returns the value compared against the tolerance.
pub type StopFn<T> = Box<dyn Fn(&[T], &[T], &[T]) -> f64 + Send + Sync>;

/// Per-iteration relative residuals of the last solve.
#[derive(Clone, Debug, Default)]
pub(crate) struct History(pub(crate) Vec<f64>);

impl History {
    pub(crate) fn push(&mut self, v: f64) {
        self.0.push(v);
    }
}

/// Maps the Fortran flags onto `HUTI_INFO` (convergence, divergence, halted,
/// otherwise max-iterations).
pub(crate) fn finish(
    iterations: usize,
    residual: f64,
    converged: bool,
    diverged: bool,
    halted: Option<String>,
) -> SolveStats {
    let status = if converged {
        SolveStatus::Converged
    } else if diverged {
        SolveStatus::Diverged
    } else if let Some(why) = halted {
        SolveStatus::Breakdown(why)
    } else {
        SolveStatus::MaxIterations
    };
    SolveStats {
        iterations,
        residual,
        status,
    }
}
