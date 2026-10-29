//! Solver keywords, read as `IterSolve.F90` reads them (lines 718-1000).

use crate::{BiCgStabL, Gcr, Idrs, Jacobi, Richardson, Sgs};
use elmer_core::{CsrMatrix, IterativeSolver, LinalgError, Params};
use num_complex::Complex64;

/// The methods of `IterativeMethods.F90` this crate implements. HUTI methods
/// (CG, CGS, BiCGStab, BiCGStab2, TFQMR, GMRES) belong to `elmer-krylov`;
/// MPRGP (bound constrained) is not covered.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IterMethod {
    Jacobi,
    Richardson,
    Sgs,
    BiCgStabL,
    Gcr,
    Idrs,
}

impl IterMethod {
    /// Method for a `Linear System Iterative Method` value (case-insensitive).
    /// `complex` mirrors IterSolve: `bicgstab2` falls back to BiCGStab(l) on
    /// complex systems.
    pub fn from_name(name: &str, complex: bool) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "bicgstabl" => Some(Self::BiCgStabL),
            "bicgstab2" if complex => Some(Self::BiCgStabL),
            "gcr" => Some(Self::Gcr),
            "idrs" => Some(Self::Idrs),
            "jacobi" => Some(Self::Jacobi),
            "richardson" => Some(Self::Richardson),
            "sgs" => Some(Self::Sgs),
            _ => None,
        }
    }
}

/// Layout of a complex system stored in real arrays (`Linear System Pseudo
/// Complex`). Only GCR supports it, as in the Fortran.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum PseudoComplex {
    #[default]
    Off,
    /// Re/Im pairs adjacent (`HUTI_PSEUDOCOMPLEX = 1`).
    Interleaved,
    /// Re block then Im block (`Block Split Complex`, `= 2`).
    Split,
}

/// `Linear System Robust` family (BiCGStab(l) and real IDR(s) only).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct RobustConfig {
    pub tolerance: f64,
    pub max_tolerance: f64,
    pub step_size: f64,
    pub max_bad_iterations: usize,
    pub start_iteration: usize,
}

#[derive(Clone, Debug)]
pub struct IterMethodsConfig {
    /// `Linear System Iterative Method`, lower-cased; `bicgstab` if absent.
    pub method_name: String,
    /// `Linear System Max Iterations` (`HUTI_MAXIT`).
    pub max_iterations: usize,
    /// `Linear System Min Iterations` (`HUTI_MINIT`); real GCR only.
    pub min_iterations: usize,
    /// `Linear System Convergence Tolerance`.
    pub tolerance: f64,
    /// `Linear System Divergence Limit`, default 1e20.
    pub divergence_limit: f64,
    /// `Linear System Residual Output`, default 1. `usize::MAX` stands for
    /// Fortran's `HUGE`: the internal methods map 0 to it. Only a hint for
    /// callers that print the recorded history.
    pub output_interval: usize,
    /// `BiCGstabl polynomial degree`, default 2, at least 2.
    pub bicgstabl_degree: usize,
    /// `Linear System GCR Restart`; default `min(max_iterations, 200)`.
    pub gcr_restart: usize,
    /// `IDRS parameter`, default 4.
    pub idrs_s: usize,
    /// `IDRS Smoothing` (real only).
    pub idrs_smoothing: bool,
    /// `SGS Overrelaxation Factor`, default 1.8.
    pub sgs_omega: f64,
    pub robust: Option<RobustConfig>,
    pub pseudo_complex: PseudoComplex,
    /// IterSolve.F90:1000 sets `x = 1e-8` before BiCGStab-family calls when the
    /// initial guess is all zero. `None` skips it for callers that did it.
    pub zero_guess_fill: Option<f64>,
    /// IDR(s) shadow-space seed (`Simulation: Random Number Seed`, default 314159265).
    pub random_seed: u64,
}

impl Default for IterMethodsConfig {
    fn default() -> Self {
        IterMethodsConfig {
            method_name: "bicgstab".into(),
            max_iterations: 100,
            min_iterations: 0,
            tolerance: 1e-8,
            divergence_limit: 1e20,
            output_interval: 1,
            bicgstabl_degree: 2,
            gcr_restart: 100,
            idrs_s: 4,
            idrs_smoothing: false,
            sgs_omega: 1.8,
            robust: None,
            pseudo_complex: PseudoComplex::Off,
            zero_guess_fill: Some(1e-8),
            random_seed: 314_159_265,
        }
    }
}

fn bad(key: &str, why: &str) -> LinalgError {
    LinalgError::Unsupported(format!("keyword '{key}': {why}"))
}

fn positive(p: &dyn Params, key: &str) -> Result<Option<usize>, LinalgError> {
    match p.integer(key) {
        None => Ok(None),
        Some(v) if v >= 1 => Ok(Some(v as usize)),
        Some(v) => Err(bad(key, &format!("value {v} below minimum 1"))),
    }
}

impl IterMethodsConfig {
    /// Reads the keywords of a solver section. `Linear System Max Iterations`
    /// and `Linear System Convergence Tolerance` are mandatory, as in
    /// `IterSolver`.
    pub fn from_params(p: &dyn Params) -> Result<Self, LinalgError> {
        let method_name = p
            .string("Linear System Iterative Method")
            .map(|s| s.trim().to_ascii_lowercase())
            .unwrap_or_else(|| "bicgstab".into());

        let max_iterations = positive(p, "Linear System Max Iterations")?
            .ok_or_else(|| bad("Linear System Max Iterations", "not found"))?;
        let tolerance = p
            .real("Linear System Convergence Tolerance")
            .ok_or_else(|| bad("Linear System Convergence Tolerance", "not found"))?;

        let bicgstabl_degree = match p.integer("BiCGstabl polynomial degree") {
            None => 2,
            Some(v) if v >= 2 => v as usize,
            Some(v) => return Err(bad("BiCGstabl polynomial degree", &format!("{v} below 2"))),
        };
        let gcr_restart = match p.integer("Linear System GCR Restart") {
            Some(v) if v >= 1 => v as usize,
            Some(v) => return Err(bad("Linear System GCR Restart", &format!("{v} below 1"))),
            None => max_iterations.min(200),
        };
        let sgs_omega = p.real("SGS Overrelaxation Factor").unwrap_or(1.8);
        if !(0.0..=2.0).contains(&sgs_omega) {
            return Err(bad("SGS Overrelaxation Factor", "outside [0, 2]"));
        }

        let tolerance_robust = tolerance.powf(2.0 / 3.0);
        let robust = if p.logical("Linear System Robust").unwrap_or(false) {
            Some(RobustConfig {
                tolerance: p
                    .real("Linear System Robust Tolerance")
                    .unwrap_or(tolerance_robust),
                max_tolerance: p
                    .real("Linear System Robust Limit")
                    .unwrap_or_else(|| tolerance.sqrt()),
                step_size: p.real("Linear System Robust Margin").unwrap_or(1.1),
                max_bad_iterations: p
                    .integer("Linear System Robust Max Iterations")
                    .map(|v| v.max(0) as usize)
                    .unwrap_or(max_iterations / 2),
                start_iteration: p
                    .integer("Linear System Robust Start Iteration")
                    .map(|v| v.max(0) as usize)
                    .unwrap_or(1),
            })
        } else {
            None
        };

        let pseudo_complex = if p.logical("Linear System Pseudo Complex").unwrap_or(false) {
            if p.logical("Block Split Complex").unwrap_or(false) {
                PseudoComplex::Split
            } else {
                PseudoComplex::Interleaved
            }
        } else {
            PseudoComplex::Off
        };

        Ok(IterMethodsConfig {
            method_name,
            max_iterations,
            min_iterations: p
                .integer("Linear System Min Iterations")
                .map(|v| v.max(0) as usize)
                .unwrap_or(0),
            tolerance,
            divergence_limit: p
                .real("Linear System Divergence Limit")
                .unwrap_or(1.0e20),
            output_interval: match p.integer("Linear System Residual Output") {
                None => 1,
                Some(v) if v <= 0 => usize::MAX,
                Some(v) => v as usize,
            },
            bicgstabl_degree,
            gcr_restart,
            idrs_s: positive(p, "IDRS parameter")?.unwrap_or(4),
            idrs_smoothing: p.logical("IDRS Smoothing").unwrap_or(false),
            sgs_omega,
            robust,
            pseudo_complex,
            zero_guess_fill: Some(1e-8),
            random_seed: p
                .integer("Random Number Seed")
                .map(|v| v as u32 as u64)
                .unwrap_or(314_159_265),
        })
    }

    /// The method this section selects, or `None` if it is a HUTI method (or
    /// unknown) that another crate implements.
    pub fn method(&self, complex: bool) -> Option<IterMethod> {
        IterMethod::from_name(&self.method_name, complex)
    }

    /// Real-valued solver for the selected method. Jacobi, Richardson and SGS
    /// read matrix entries, so they take `a` (Jacobi and Richardson copy what
    /// they need; SGS borrows it).
    pub fn build_real<'a>(
        &self,
        a: &'a CsrMatrix<f64>,
    ) -> Result<Box<dyn IterativeSolver<f64> + 'a>, LinalgError> {
        Ok(match self.method(false) {
            Some(IterMethod::Jacobi) => Box::new(Jacobi::new(self, a)?),
            Some(IterMethod::Richardson) => Box::new(Richardson::new(self, a)),
            Some(IterMethod::Sgs) => Box::new(Sgs::new(self, a)),
            Some(IterMethod::BiCgStabL) => Box::new(BiCgStabL::<f64>::new(self)),
            Some(IterMethod::Gcr) => Box::new(Gcr::<f64>::new(self)),
            Some(IterMethod::Idrs) => Box::new(Idrs::<f64>::new(self)),
            None => {
                return Err(LinalgError::Unsupported(format!(
                    "iterative method '{}' is not in elmer-itermethods",
                    self.method_name
                )))
            }
        })
    }

    /// Complex-valued solver (`itermethod_z_*`): GCR, BiCGStab(l) and IDR(s).
    pub fn build_complex(&self) -> Result<Box<dyn IterativeSolver<Complex64>>, LinalgError> {
        Ok(match self.method(true) {
            Some(IterMethod::BiCgStabL) => Box::new(BiCgStabL::<Complex64>::new(self)),
            Some(IterMethod::Gcr) => Box::new(Gcr::<Complex64>::new(self)),
            Some(IterMethod::Idrs) => Box::new(Idrs::<Complex64>::new(self)),
            _ => {
                return Err(LinalgError::Unsupported(
                    "Complex arithmetic version of the given linear solver is not available".into(),
                ))
            }
        })
    }
}
