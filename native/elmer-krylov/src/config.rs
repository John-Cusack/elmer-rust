//! HUTI parameters and the IterSolve keyword mapping that fills them.

use elmer_core::{LinalgError, Logger, NullLogger, Params};

/// The HUTI methods (fhutiter/src/huti_*.F90).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum HutiMethod {
    Cg,
    Cgs,
    BiCGStab,
    BiCGStab2,
    Gmres,
    /// Not reachable from IterSolve ("qmr" falls to its BiCGStab default);
    /// available by direct construction. Needs `apply_transpose` on both the
    /// operator and the preconditioner.
    Qmr,
    Tfqmr,
}

impl HutiMethod {
    pub fn name(self) -> &'static str {
        match self {
            HutiMethod::Cg => "cg",
            HutiMethod::Cgs => "cgs",
            HutiMethod::BiCGStab => "bicgstab",
            HutiMethod::BiCGStab2 => "bicgstab2",
            HutiMethod::Gmres => "gmres",
            HutiMethod::Qmr => "qmr",
            HutiMethod::Tfqmr => "tfqmr",
        }
    }
}

/// HUTI stopping criteria (`HUTI_STOPC`, huti_fdefs.h:86-95). The
/// user-supplied criterion (10) is not here: IterSolve only uses it for the
/// backward-error keywords, which [`HutiConfig::from_params`] refuses.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum HutiStop {
    /// `||A x - b||` (0).
    TrueResidual,
    /// `||A x - b|| / ||b||` (1); what IterSolve always selects.
    TrueResidualScaledByB,
    /// The method's recursive residual `||r||` (2).
    PseudoResidual,
    /// `||r|| / ||b||` (3).
    PseudoResidualScaledByB,
    /// `||r|| / ||M_L^-1 b||` (4).
    PseudoResidualScaledByPrecB,
    /// Norm of the last update of `x` (5).
    XDiffNorm,
    /// TFQMR's `sqrt(2 it) tau / ||b||` bound, checked against the true
    /// residual near convergence (6). Other methods treat it as their
    /// `case default`, the unscaled true residual.
    UpperBound,
}

/// Everything one HUTI call reads from `ipar`/`dpar` plus the IterSolve
/// decisions around it.
#[derive(Clone, Debug, PartialEq)]
pub struct HutiConfig {
    pub method: HutiMethod,
    /// `HUTI_MAXIT`. An iteration count above this ends with `MaxIterations`.
    pub max_iter: usize,
    /// `HUTI_MINIT`. IterSolve sets it, but no HUTI method reads it; kept so a
    /// caller handing the same config to IterativeMethods solvers loses nothing.
    pub min_iter: usize,
    /// `HUTI_TOLERANCE`: converged when the monitored residual is below it.
    pub tol: f64,
    /// `HUTI_MAXTOLERANCE`: diverged when the residual is NaN or above it.
    pub max_tol: f64,
    pub stop_crit: HutiStop,
    /// `HUTI_GMRES_RESTART` (m).
    pub restart: usize,
    /// `HUTI_BICGSTABL_L`, only for BiCGStab(l) in IterativeMethods (lane L03).
    pub ell: usize,
    /// `HUTI_DBUGLVL`: print every this many iterations; 0 prints nothing.
    pub verbosity: usize,
    /// Apply the preconditioner from the left (`pcondlsubr`) rather than the
    /// right (`pcondrsubr`), IterSolve.F90:1005-1011.
    pub left_precond: bool,
    /// IterSolve.F90:962-963: BiCGStab, BiCGStab(2) and BiCGStab(l) replace an
    /// all-zero initial guess by 1e-8 in every real component.
    pub perturb_zero_guess: bool,
    /// "Linear System Abort Not Converged" (default true): whether the caller
    /// should treat `MaxIterations` as a numerical error (IterSolve.F90:1652-1659).
    pub abort_not_converged: bool,
}

impl HutiConfig {
    /// HUTI defaults (huti_fdefs.h:11-14: 5000 iterations, tolerance 1e-6, no
    /// output) with IterSolve's divergence limit, stopping rule, GMRES restart
    /// and preconditioning side for `method`.
    pub fn new(method: HutiMethod) -> Self {
        HutiConfig {
            method,
            max_iter: 5000,
            min_iter: 0,
            tol: 1.0e-6,
            max_tol: 1.0e20,
            stop_crit: HutiStop::TrueResidualScaledByB,
            restart: 10,
            ell: 2,
            verbosity: 0,
            left_precond: matches!(
                method,
                HutiMethod::BiCGStab2 | HutiMethod::Gmres | HutiMethod::Tfqmr
            ),
            perturb_zero_guess: matches!(method, HutiMethod::BiCGStab | HutiMethod::BiCGStab2),
            abort_not_converged: true,
        }
    }

    /// The real-valued IterSolve mapping, logging nothing.
    pub fn from_params(params: &dyn Params) -> Result<Self, LinalgError> {
        Self::from_params_with(params, false, &NullLogger)
    }

    /// IterSolve's keyword mapping (IterSolve.F90:724-1011) for a real or
    /// complex system. Methods that IterSolve sends to IterativeMethods.F90
    /// (and complex BiCGStab(2), which it remaps to BiCGStab(l)) return
    /// `Unsupported` until lane L03 lands them.
    pub fn from_params_with(
        params: &dyn Params,
        complex: bool,
        log: &dyn Logger,
    ) -> Result<Self, LinalgError> {
        let name = match params.string("Linear System Iterative Method") {
            Some(s) => {
                let s = s.trim().to_ascii_lowercase();
                log.info("IterSolver", &format!("Using iterative method: {s}"), 9);
                s
            }
            None => {
                log.warn(
                    "IterSolver",
                    "> Linear System Iterative Method < not found, using BiCGstab",
                );
                "bicgstab".to_string()
            }
        };
        let method = match name.as_str() {
            "bicgstab2" if complex => return Err(itermethods("bicgstabl (complex bicgstab2)")),
            "bicgstab2" => HutiMethod::BiCGStab2,
            "tfqmr" => HutiMethod::Tfqmr,
            "cgs" => HutiMethod::Cgs,
            "cg" => HutiMethod::Cg,
            "gmres" => HutiMethod::Gmres,
            "bicgstabl" | "sgs" | "jacobi" | "richardson" | "gcr" | "idrs" | "mprgp" => {
                return Err(itermethods(&name))
            }
            // CASE DEFAULT (IterSolve.F90:787): any other name, "bicgstab" included.
            _ => HutiMethod::BiCGStab,
        };
        let mut cfg = HutiConfig::new(method);

        if method == HutiMethod::Gmres {
            let m = params.integer("Linear System GMRES Restart").unwrap_or(10);
            // Fortran would index an empty Krylov basis; refuse instead.
            if m < 1 {
                return Err(LinalgError::Unsupported(format!(
                    "Linear System GMRES Restart = {m}: HUTI GMRES needs at least 1"
                )));
            }
            cfg.restart = m as usize;
        }
        if let Some(l) = params.integer("BiCGstabl polynomial degree") {
            cfg.ell = l.max(0) as usize;
        }

        if params.logical("Linear System Componentwise Backward Error") == Some(true) {
            return Err(LinalgError::Unsupported(
                "Linear System Componentwise Backward Error (BackwardError.F90)".into(),
            ));
        }
        if params.logical("Linear System Normwise Backward Error") == Some(true) {
            return Err(LinalgError::Unsupported(
                "Linear System Normwise Backward Error (BackwardError.F90)".into(),
            ));
        }
        cfg.stop_crit = HutiStop::TrueResidualScaledByB;

        // Ranks other than 0 print nothing; Huti::solve applies that from the
        // communicator (IterSolve.F90:932).
        cfg.verbosity = params
            .integer("Linear System Residual Output")
            .unwrap_or(1)
            .unsigned_abs() as usize;

        // ListGetInteger(..., minv=1): absent is 0, present and below 1 is fatal.
        cfg.max_iter = match params.integer("Linear System Max Iterations") {
            Some(v) if v < 1 => {
                return Err(LinalgError::Unsupported(format!(
                    "Linear System Max Iterations = {v} is smaller than the minimum 1"
                )))
            }
            Some(v) => v as usize,
            None => 0,
        };
        cfg.min_iter = params
            .integer("Linear System Min Iterations")
            .unwrap_or(0)
            .max(0) as usize;

        // ListGetCReal returns 0 when absent.
        cfg.tol = params
            .real("Linear System Convergence Tolerance")
            .unwrap_or(0.0);
        cfg.max_tol = params
            .real("Linear System Divergence Limit")
            .unwrap_or(1.0e20);

        cfg.left_precond = match method {
            HutiMethod::BiCGStab2 | HutiMethod::Gmres | HutiMethod::Tfqmr => true,
            _ => params
                .logical("Linear System Left Preconditioning")
                .unwrap_or(false),
        };
        cfg.abort_not_converged = params
            .logical("Linear System Abort Not Converged")
            .unwrap_or(true);
        Ok(cfg)
    }
}

fn itermethods(name: &str) -> LinalgError {
    LinalgError::Unsupported(format!(
        "iterative method '{name}' is an IterativeMethods.F90 solver (lane L03), not HUTI"
    ))
}
