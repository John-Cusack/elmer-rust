//! The default [`LinearSolverRegistry`]: HUTI Krylov methods from
//! `elmer-krylov`, BiCGStab(l)/GCR/IDR(s) from `elmer-itermethods`, direct
//! methods from `elmer-direct` and the ILU family and
//! diagonal preconditioner from `elmer-sparse`, keyed as IterSolve.F90 and
//! DirectSolve.F90 read the `Linear System *` keywords.

use elmer_core::{normalize_key, CsrMatrix, DirectSolver, IterativeSolver, ParamMap, Params, Preconditioner};
use elmer_direct::{kind_for_method, BandCholesky, BandLu, DirectKind, FaerCholesky, FaerLu};
use elmer_itermethods::{BiCgStabL, Gcr, Idrs, IterMethodsConfig};
use elmer_krylov::{Huti, HutiConfig};
use elmer_sparse::Ilu;

use crate::backend::LinearSolverRegistry;
use crate::error::SolverError;

/// Every `Linear System Direct Method` value `elmer-direct` serves
/// (DirectSolve.F90:4266).
pub const DIRECT_METHODS: &[&str] = &[
    "banded",
    "symmetric banded",
    "umfpack",
    "big umfpack",
    "smumps",
    "cmumps",
    "mumps",
    "dmumps",
    "zmumps",
    "mumpslocal",
    "superlu",
    "pardiso",
    "cpardiso",
    "cholmod",
];

/// The registry the default [`crate::Backends`] carry.
pub fn elmer_registry() -> LinearSolverRegistry {
    let mut r = LinearSolverRegistry::new();
    // IterSolve's CASE DEFAULT sends any unlisted method name to BiCGStab;
    // HutiConfig::from_params does that mapping and refuses the
    // IterativeMethods.F90 methods (lane L03).
    r.set_iterative_default(huti);
    for m in ["bicgstabl", "gcr", "idrs"] {
        r.register_iterative(m, itermethods);
    }
    for m in ["jacobi", "richardson", "sgs"] {
        r.register_iterative(m, stationary);
    }
    for m in DIRECT_METHODS {
        r.register_direct(m, direct);
    }
    r.register_precond("diagonal", precond);
    r.register_precond("ilut", precond);
    r.register_precond("ilu", precond);
    for k in 0..=9 {
        r.register_precond(&format!("ilu{k}"), precond);
    }
    r
}

fn huti(p: &ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    Ok(Box::new(Huti::<f64>::new(HutiConfig::from_params(p)?)))
}

/// IterativeMethods.F90 Krylov methods; IterSolve.F90 routes these names there.
fn itermethods(p: &ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    let cfg = IterMethodsConfig::from_params(p)?;
    Ok(match cfg.method_name.as_str() {
        "bicgstabl" => Box::new(BiCgStabL::<f64>::new(&cfg)),
        "gcr" => Box::new(Gcr::<f64>::new(&cfg)),
        _ => Box::new(Idrs::<f64>::new(&cfg)),
    })
}

/// Jacobi, Richardson and SGS are built from the matrix, which the
/// registry's factories do not receive yet.
fn stationary(p: &ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    Err(SolverError::NotWired(format!(
        "iterative method {} (needs the matrix at construction)",
        p.string("linear system iterative method").unwrap_or("")
    )))
}

/// `Linear System Direct Method`, default `banded` (DirectSolve.F90).
fn direct(p: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError> {
    let kind = kind_for_method(p.string("linear system direct method").unwrap_or("banded"))?;
    Ok(match kind {
        DirectKind::BandLu => Box::new(BandLu::<f64>::new()),
        DirectKind::BandCholesky => Box::new(BandCholesky::<f64>::new()),
        DirectKind::SparseLu => Box::new(FaerLu::<f64>::new()),
        DirectKind::SparseCholesky => Box::new(FaerCholesky::<f64>::new()),
    })
}

/// IterSolve.F90:1050-1100 for the CRS preconditioners: `diagonal`,
/// `ilut` (`Linear System ILUT Tolerance`), `ilu[n]` (`Linear System ILU
/// Order` first, else the digit in the name, out of 0..9 -> 0), and incomplete
/// Cholesky with `Linear System Symmetric ILU`. Recomputed every solve, as
/// Elmer's default; reuse options and `Linear System ILU Factor` are refused.
fn precond(p: &ParamMap, a: &CsrMatrix<f64>) -> Result<Box<dyn Preconditioner<f64> + Send>, SolverError> {
    if p.logical("no precondition recompute") == Some(true)
        || p.integer("linear system precondition recompute").is_some_and(|n| n > 1)
        || p.logical("linear system refactorize") == Some(false)
    {
        return Err(SolverError::Unsupported("preconditioner reuse across solves".into()));
    }
    if p.real("linear system ilu factor").is_some_and(|f| f > f64::EPSILON) {
        return Err(SolverError::Unsupported("Linear System ILU Factor".into()));
    }
    let name = normalize_key(p.string("linear system preconditioning").unwrap_or("none"));
    if name == "diagonal" {
        return Ok(Box::new(Ilu::diagonal(a)));
    }
    if name == "ilut" {
        let tol = p.real("linear system ilut tolerance").unwrap_or(0.0);
        return Ok(Box::new(Ilu::ilut(a, tol)?));
    }
    let mut level = match p.integer("linear system ilu order") {
        Some(k) => k,
        None => name
            .as_bytes()
            .get(3)
            .map_or(-1, |&c| i32::from(c) - i32::from(b'0')),
    };
    if !(0..=9).contains(&level) {
        level = 0;
    }
    if p.logical("linear system symmetric ilu") == Some(true) {
        let retry = if p.logical("linear system cholesky retry") == Some(true) {
            Some(p.real("linear system cholesky factor").ok_or_else(|| {
                SolverError::Unsupported("Linear System Cholesky Retry without Linear System Cholesky Factor".into())
            })?)
        } else {
            None
        };
        return Ok(Box::new(Ilu::ic_n(a, level as usize, retry)?));
    }
    Ok(Box::new(Ilu::ilu_n(a, level as usize)?))
}
