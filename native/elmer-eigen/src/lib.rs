//! `Eigen Analysis` for the Rust-native Elmer rewrite: the ARPACK-based
//! `ArpackEigenSolve` (EigenSolve.F90:140) without ARPACK.
//!
//! * [`eigs`]: `K x = lambda M x` with a consistent mass (shift-invert,
//!   ARPACK mode 3), a lumped mass (mode 2) or `M = I`, following Elmer's
//!   `Eigen System Select`, `Eigen System Sorting` and normalisation keywords.
//! * [`krylov_schur`]: the restarted iteration itself, over any
//!   [`EigenOperator`]: thick-restart Lanczos for symmetric problems
//!   (`DSAUPD`), Krylov-Schur style thick-restart Arnoldi otherwise (`DNAUPD`).
//! * [`EigenConfig::from_params`]: the keyword defaults of `ArpackEigenSolve`.
//!
//! Not covered (refused by [`EigenConfig::from_params`] or not offered):
//! `Eigen System Damped` (`ArpackDampedEigenSolve`), `Stability Analysis`
//! (`ArpackStabEigenSolve`), complex matrices (`ArpackEigenSolveComplex`) and
//! the parallel solvers (`ParallelEigenSolve.F90`).
//!
//! Eigenvectors come back with unit `M`-norm and an arbitrary sign, like
//! ARPACK's; Elmer's random start vector differs from this crate's, so
//! compare vectors sign-normalised.
#![forbid(unsafe_code)]

pub mod config;
pub mod krylov;
pub mod scale;
pub mod solve;
pub mod sort;

pub use config::{ArpackWhich, EigenConfig, Normalize, Which};
pub use krylov::{krylov_schur, EigenOperator, KrylovParams, RitzPairs};
pub use scale::{check_residuals, scale_eigen_vectors};
pub use solve::{eigs, shifted_matrix, EigenResult, Mass};
pub use sort::{eigen_system_sorting, sort_r};

use elmer_core::Logger;

/// Fortran `ESw.d`: one leading digit, `d` decimals, `E+xx` exponent,
/// right-justified in `w` columns.
pub fn fortran_es(x: f64, w: usize, d: usize) -> String {
    let s = format!("{:.*e}", d, x);
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let e: i32 = exp.parse().unwrap_or(0);
    let digits = if e.abs() > 99 { 3 } else { 2 };
    let body = format!(
        "{mant}E{}{:0digits$}",
        if e < 0 { '-' } else { '+' },
        e.abs()
    );
    format!("{body:>w$}")
}

/// The Info lines `ArpackEigenSolve` prints after a solve (caller
/// `EigenSolve`). `tol` is the configured tolerance as Elmer prints it.
pub fn log_summary(res: &EigenResult, tol: f64, log: &dyn Logger) {
    const CALLER: &str = "EigenSolve";
    log.info(CALLER, " ", 4);
    log.info(CALLER, "Eigen system solution complete: ", 4);
    log.info(CALLER, " ", 4);
    log.info(
        CALLER,
        &format!("Convergence criterion is: {}", fortran_es(tol, 12, 3)),
        7,
    );
    log.info(
        CALLER,
        &format!(
            "Number of eigensystem iterations is: {}",
            res.op_applications + 1
        ),
        4,
    );
    log.info(
        CALLER,
        &format!("Number of converged Ritz values is: {}", res.converged),
        4,
    );
    log.info(
        CALLER,
        &format!("Number of update iterations taken: {}", res.iterations),
        4,
    );
    log.info(
        CALLER,
        &format!("Computed {} Eigen Values", res.values.len()),
        4,
    );
    for (i, v) in res.values.iter().enumerate() {
        log.info(
            CALLER,
            &format!(
                "{}: {}{}",
                i + 1,
                fortran_es(v.re, 15, 6),
                fortran_es(v.im, 15, 6)
            ),
            4,
        );
    }
    log.info(CALLER, "--------------------------------", 4);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es_format_matches_gfortran() {
        assert_eq!(fortran_es(121834.3221170, 15, 6), "   1.218343E+05");
        assert_eq!(fortran_es(0.0, 15, 6), "   0.000000E+00");
        assert_eq!(fortran_es(1e-6, 12, 3), "   1.000E-06");
        assert_eq!(fortran_es(-2.5e-120, 15, 6), " -2.500000E-120");
    }
}
