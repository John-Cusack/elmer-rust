//! The default [`LinearSolverRegistry`]: HUTI Krylov methods from
//! `elmer-krylov`, BiCGStab(l)/GCR/IDR(s) from `elmer-itermethods`, direct
//! methods from `elmer-direct` and the ILU family and
//! diagonal preconditioner from `elmer-sparse`, keyed as IterSolve.F90 and
//! DirectSolve.F90 read the `Linear System *` keywords.

use elmer_core::{
    normalize_key, Communicator, CsrMatrix, DirectSolver, IterativeSolver, LinalgError, LinearOperator, ParamMap, Params,
    Preconditioner, SolveStats,
};
use elmer_direct::{kind_for_method, BandCholesky, BandLu, DirectKind, FaerCholesky, FaerLu};
use num_complex::Complex64;
use elmer_itermethods::{BiCgStabL, Gcr, Idrs, IterMethodsConfig, Jacobi, Richardson, Sgs};
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
        r.register_iterative_with_matrix(m, stationary_for);
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

/// Jacobi, Richardson and SGS outside [`LinearSolverRegistry::select_for`]:
/// they are built from the matrix.
fn stationary(p: &ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    Err(SolverError::NotWired(format!(
        "iterative method {} (needs the matrix at construction)",
        p.string("linear system iterative method").unwrap_or("")
    )))
}

/// Jacobi, Richardson and SGS (IterativeMethods.F90:296, :404, :178) from
/// the assembled matrix. They ignore the preconditioner, as the Fortran
/// does.
fn stationary_for(p: &ParamMap, a: &CsrMatrix<f64>) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    let cfg = IterMethodsConfig::from_params(p)?;
    Ok(match cfg.method_name.as_str() {
        "jacobi" => Box::new(Jacobi::new(&cfg, a)?),
        "richardson" => Box::new(Richardson::new(&cfg, a)),
        _ => Box::new(OwnedSgs { cfg, a: a.clone() }),
    })
}

/// [`Sgs`] borrows its matrix; the registry hands out owned solvers, so this
/// keeps a copy of the assembled matrix for the one solve it is built for.
struct OwnedSgs {
    cfg: IterMethodsConfig,
    a: CsrMatrix<f64>,
}

impl IterativeSolver<f64> for OwnedSgs {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<f64>,
        m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        Sgs::new(&self.cfg, &self.a).solve(a, m, b, x, comm)
    }
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

/// Complex systems (`Linear System Complex`): the solve on the native
/// complex form of the interleaved real system, with the same keywords as
/// [`elmer_registry`] reads. Direct: `banded` (`ZGBTF2`/`ZGBTRS`) and the
/// sparse LU methods; iterative: the HUTI methods in complex arithmetic and
/// BiCGStab(l), GCR and IDR(s); preconditioners `none`, `diagonal`, `ilu[n]`
/// and `ilut`. Cholesky methods, the stationary methods, complex
/// BiCGStab(2), incomplete Cholesky and preconditioner reuse are refused.
pub(crate) fn solve_complex(
    p: &ParamMap,
    a: &CsrMatrix<Complex64>,
    b: &[Complex64],
    x: &mut [Complex64],
    log: &dyn elmer_core::Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    let kind = p.string("linear system solver").map(normalize_key).ok_or_else(|| {
        SolverError::fatal("SolveLinearSystem", "Give \"Linear System Solver\", e.g. \"iterative\" or \"direct\"")
    })?;
    match kind.as_str() {
        "direct" => {
            let mut d: Box<dyn DirectSolver<Complex64>> =
                match kind_for_method(p.string("linear system direct method").unwrap_or("banded"))? {
                    DirectKind::BandLu => Box::new(BandLu::<Complex64>::new()),
                    DirectKind::SparseLu => Box::new(FaerLu::<Complex64>::new()),
                    DirectKind::BandCholesky | DirectKind::SparseCholesky => {
                        return Err(SolverError::Unsupported("complex Cholesky direct method".into()))
                    }
                };
            d.factor(a)?;
            let mut y = b.to_vec();
            d.solve_in_place(&mut y)?;
            x.copy_from_slice(&y);
            Ok(())
        }
        "iterative" => {
            let m = p
                .string("linear system iterative method")
                .map(normalize_key)
                .unwrap_or_else(|| "bicgstab".into());
            let mut solver: Box<dyn IterativeSolver<Complex64>> = match m.as_str() {
                "bicgstabl" | "gcr" | "idrs" => {
                    let cfg = IterMethodsConfig::from_params(p)?;
                    match m.as_str() {
                        "bicgstabl" => Box::new(BiCgStabL::<Complex64>::new(&cfg)),
                        "gcr" => Box::new(Gcr::<Complex64>::new(&cfg)),
                        _ => Box::new(Idrs::<Complex64>::new(&cfg)),
                    }
                }
                "jacobi" | "richardson" | "sgs" | "bicgstab2" => {
                    return Err(SolverError::Unsupported(format!("complex iterative method {m}")))
                }
                _ => Box::new(Huti::<Complex64>::new(HutiConfig::from_params_with(p, true, log)?)),
            };
            let precond = complex_precond(p, a, log)?;
            solver.solve(a, &*precond, b, x, comm)?;
            Ok(())
        }
        other => Err(SolverError::Unsupported(format!("Linear System Solver = {other}"))),
    }
}

fn complex_precond(
    p: &ParamMap,
    a: &CsrMatrix<Complex64>,
    log: &dyn elmer_core::Logger,
) -> Result<Box<dyn Preconditioner<Complex64>>, SolverError> {
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
    if name == "none" {
        return Ok(Box::new(elmer_core::IdentityPrecond));
    }
    if name == "diagonal" {
        return Ok(Box::new(Ilu::diagonal(a)));
    }
    if name == "ilut" {
        let tol = p.real("linear system ilut tolerance").unwrap_or(0.0);
        return Ok(Box::new(Ilu::ilut(a, tol)?));
    }
    if name.starts_with("ilu") {
        if p.logical("linear system symmetric ilu") == Some(true) {
            return Err(SolverError::Unsupported("complex incomplete Cholesky".into()));
        }
        let mut level = match p.integer("linear system ilu order") {
            Some(k) => k,
            None => name.as_bytes().get(3).map_or(-1, |&c| i32::from(c) - i32::from(b'0')),
        };
        if !(0..=9).contains(&level) {
            level = 0;
        }
        return Ok(Box::new(Ilu::ilu_n(a, level as usize)?));
    }
    if ["bilu", "multigrid", "vanka", "auxiliary space solver", "slave", "circuit"]
        .iter()
        .any(|k| name.starts_with(k))
    {
        return Err(SolverError::NotWired(format!("preconditioner {name}")));
    }
    log.warn("IterSolve", &format!("Unknown preconditioner type: {name}, feature disabled."));
    Ok(Box::new(elmer_core::IdentityPrecond))
}

/// The native complex form of an interleaved real system (row `2i` holds
/// `(re, -im)` of each complex entry, CRSMatrix.F90's complex routines):
/// entry `(i, j)` is `A(2i, 2j) - i A(2i, 2j+1)`; the vectors pair
/// `(v(2i), v(2i+1))`.
pub(crate) fn to_complex(a: &CsrMatrix<f64>) -> Result<CsrMatrix<Complex64>, SolverError> {
    if !a.nrows.is_multiple_of(2) {
        return Err(SolverError::fatal("SolveLinearSystem", "complex system with an odd number of rows"));
    }
    let n = a.nrows / 2;
    let mut c = CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(n);
    for i in 0..n {
        let mut row: std::collections::BTreeMap<u32, Complex64> = std::collections::BTreeMap::new();
        for k in a.row_ptr[2 * i]..a.row_ptr[2 * i + 1] {
            let col = a.col_idx[k];
            let e = row.entry(col / 2).or_insert(Complex64::new(0.0, 0.0));
            if col.is_multiple_of(2) {
                e.re += a.values[k];
            } else {
                e.im -= a.values[k];
            }
        }
        let mut d = None;
        for (j, v) in row {
            if j as usize == i {
                d = Some(c.values.len());
            }
            c.col_idx.push(j);
            c.values.push(v);
        }
        diag.push(d);
        c.row_ptr.push(c.values.len());
    }
    c.diag = diag.into_iter().collect();
    Ok(c)
}

pub(crate) fn pair_up(v: &[f64]) -> Vec<Complex64> {
    v.as_chunks::<2>().0.iter().map(|&[re, im]| Complex64::new(re, im)).collect()
}
