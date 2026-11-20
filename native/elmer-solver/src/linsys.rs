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
    for m in ["gcr", "idrs"] {
        r.register_iterative_with_matrix(m, itermethods_for);
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
    refuse_backward_error(p)?;
    let cfg = IterMethodsConfig::from_params(p)?;
    Ok(match cfg.method_name.as_str() {
        "bicgstabl" => Box::new(BiCgStabL::<f64>::new(&cfg)),
        "gcr" => Box::new(Gcr::<f64>::new(&cfg)),
        _ => Box::new(Idrs::<f64>::new(&cfg)),
    })
}

/// GCR and IDR(s) with `Linear System Componentwise Backward Error`: HUTI's
/// user stop criterion `ComponentwiseBackwardError` (BackwardError.F90:170)
/// on the assembled matrix. Without the keyword, `None` defers to
/// [`itermethods`].
fn itermethods_for(p: &ParamMap, a: &CsrMatrix<f64>) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    if p.logical("linear system componentwise backward error") != Some(true) {
        return itermethods(p);
    }
    let cfg = IterMethodsConfig::from_params(p)?;
    let a = a.clone();
    let stop: elmer_itermethods::StopFn<f64> =
        Box::new(move |x: &[f64], b: &[f64], _r: &[f64]| componentwise_backward_error(&a, x, b));
    Ok(match cfg.method_name.as_str() {
        "gcr" => Box::new(Gcr::<f64>::new(&cfg).with_stop_criterion(stop)),
        _ => Box::new(Idrs::<f64>::new(&cfg).with_stop_criterion(stop)),
    })
}

/// IterSolve.F90:894 `Linear System Componentwise Backward Error` makes
/// [`componentwise_backward_error`] the stop criterion; it needs the
/// assembled matrix, which only [`itermethods_for`] has, and GCR or IDR(s),
/// which take a user criterion. Anything else is refused rather than run to
/// the default `||r|| / ||b||`.
fn refuse_backward_error(p: &ParamMap) -> Result<(), SolverError> {
    if p.logical("linear system componentwise backward error") == Some(true) {
        return Err(SolverError::Unsupported(
            "Linear System Componentwise Backward Error (BackwardError.F90) for this iterative method".into(),
        ));
    }
    Ok(())
}

/// `ComponentwiseBackwardError( x, b, r )`: `max_i |Ax - b|_i / (|A||x| + |b|)_i`.
fn componentwise_backward_error(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> f64 {
    const AEPS: f64 = 10.0 * f64::EPSILON;
    let n = a.nrows;
    let mut res = vec![0.0f64; n];
    crate::matvec::crs_matrix_vector_multiply(a, 1, x, &mut res);
    let ax: Vec<f64> = x.iter().map(|v| v.abs()).collect();
    let mut err = 0.0f64;
    for i in 0..n {
        let r = res[i] - b[i];
        // CRS_ABSMatrixVectorMultiply: unrolled by four.
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        let k = hi - lo;
        let (mut r1, mut r2, mut r3, mut r4) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut j = lo;
        while j < lo + 4 * (k / 4) {
            r1 += ax[a.col_idx[j] as usize] * a.values[j].abs();
            r2 += ax[a.col_idx[j + 1] as usize] * a.values[j + 1].abs();
            r3 += ax[a.col_idx[j + 2] as usize] * a.values[j + 2].abs();
            r4 += ax[a.col_idx[j + 3] as usize] * a.values[j + 3].abs();
            j += 4;
        }
        while j < hi {
            r1 += ax[a.col_idx[j] as usize] * a.values[j].abs();
            j += 1;
        }
        let d = r1 + r2 + r3 + r4 + b[i].abs();
        if d.abs() < AEPS {
            if r.abs() > AEPS {
                return f64::MAX;
            }
            err = err.max(0.0);
        } else {
            err = err.max(r.abs() / d);
        }
    }
    err
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
/// `umfpack` / `big umfpack` run the bit-exact UMFPACK 4.4 port
/// (`elmer-umfpack`), the other sparse LU methods faer's LU.
fn direct(p: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError> {
    let method = p.string("linear system direct method").unwrap_or("banded");
    if elmer_umfpack::is_umfpack_method(method) {
        return Ok(Box::new(elmer_umfpack::UmfpackLu::new(|| Box::new(FaerLu::<f64>::new()))));
    }
    let kind = kind_for_method(method)?;
    Ok(match kind {
        DirectKind::BandLu => Box::new(BandLu::<f64>::new()),
        DirectKind::BandCholesky => Box::new(BandCholesky::<f64>::new()),
        DirectKind::SparseLu => Box::new(FaerLu::<f64>::new()),
        DirectKind::SparseCholesky => Box::new(FaerCholesky::<f64>::new()),
    })
}

/// `No Precondition Recompute`, `Linear System Precondition Recompute > 1`
/// or `Linear System Refactorize = False`: the preconditioner (or the
/// factorization) of a matrix is kept for its later solves.
pub fn precond_reuse(p: &ParamMap) -> bool {
    p.logical("no precondition recompute") == Some(true)
        || p.integer("linear system precondition recompute").is_some_and(|n| n > 1)
        || p.logical("linear system refactorize") == Some(false)
}

/// `p` without the [`precond_reuse`] keywords: for the first solve of a
/// matrix, whose preconditioner is computed whatever they say.
pub fn without_precond_reuse(p: &ParamMap) -> ParamMap {
    let mut q = p.clone();
    for k in ["no precondition recompute", "linear system precondition recompute", "linear system refactorize"] {
        q.0.remove(k);
    }
    q
}

/// IterSolve.F90:1050-1100 for the CRS preconditioners: `diagonal`,
/// `ilut` (`Linear System ILUT Tolerance`), `ilu[n]` (`Linear System ILU
/// Order` first, else the digit in the name, out of 0..9 -> 0), and incomplete
/// Cholesky with `Linear System Symmetric ILU`. Recomputed every solve, as
/// Elmer's default; reuse options and `Linear System ILU Factor` are refused.
fn precond(p: &ParamMap, a: &CsrMatrix<f64>) -> Result<Box<dyn Preconditioner<f64> + Send>, SolverError> {
    if precond_reuse(p) {
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
    solve_complex_op(p, a, false, b, x, log, comm)
}

/// [`solve_complex`] with the iterative methods applying `a` as Elmer's
/// IterSolver does to a complex child matrix
/// ([`crate::matvec::crs_complex_matrix_vector_multiply`]).
pub(crate) fn solve_complex_elmer_product(
    p: &ParamMap,
    a: &CsrMatrix<Complex64>,
    b: &[Complex64],
    x: &mut [Complex64],
    log: &dyn elmer_core::Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    solve_complex_op(p, a, true, b, x, log, comm)
}

fn solve_complex_op(
    p: &ParamMap,
    a: &CsrMatrix<Complex64>,
    elmer_product: bool,
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
            let method = p.string("linear system direct method").unwrap_or("banded");
            if elmer_umfpack::is_umfpack_method(method) {
                // Umfpack_SolveSystem (DirectSolve.F90) has only the real
                // umf4 routines: a complex matrix goes to UMFPACK as its
                // interleaved real form (after the complex scaling).
                let ar = to_interleaved(a);
                let mut d = direct(p)?;
                d.factor(&ar)?;
                let mut y: Vec<f64> = b.iter().flat_map(|z| [z.re, z.im]).collect();
                d.solve_in_place(&mut y)?;
                for (k, xk) in x.iter_mut().enumerate() {
                    *xk = Complex64::new(y[2 * k], y[2 * k + 1]);
                }
                return Ok(());
            }
            let mut d: Box<dyn DirectSolver<Complex64>> =
                match kind_for_method(method)? {
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
            if elmer_product {
                solver.solve(&crate::matvec::ElmerComplexMatvec(a), &*precond, b, x, comm)?;
            } else {
                solver.solve(a, &*precond, b, x, comm)?;
            }
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

/// The interleaved real form of a native complex matrix, the inverse of
/// [`to_complex`]: entry `(i, j) = re + i im` gives `A(2i, 2j) = re`,
/// `A(2i, 2j+1) = -im`, `A(2i+1, 2j) = im`, `A(2i+1, 2j+1) = re`.
pub(crate) fn to_interleaved(c: &CsrMatrix<Complex64>) -> CsrMatrix<f64> {
    let n = 2 * c.nrows;
    let mut a = CsrMatrix {
        nrows: n,
        ncols: 2 * c.ncols,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(n);
    for i in 0..c.nrows {
        let row = c.row_ptr[i]..c.row_ptr[i + 1];
        for half in 0..2u32 {
            let r = 2 * i as u32 + half;
            let mut d = None;
            for k in row.clone() {
                let (j, z) = (c.col_idx[k], c.values[k]);
                for (col, v) in [(2 * j, if half == 0 { z.re } else { z.im }), (2 * j + 1, if half == 0 { -z.im } else { z.re })] {
                    if col == r {
                        d = Some(a.values.len());
                    }
                    a.col_idx.push(col);
                    a.values.push(v);
                }
            }
            diag.push(d);
            a.row_ptr.push(a.values.len());
        }
    }
    a.diag = diag.into_iter().collect();
    a
}

pub(crate) fn pair_up(v: &[f64]) -> Vec<Complex64> {
    v.as_chunks::<2>().0.iter().map(|&[re, im]| Complex64::new(re, im)).collect()
}

/// The ILU level of `Linear System Preconditioning = BILU<n>` (`n` one
/// digit, 0 otherwise), or `None` for another preconditioner
/// (IterSolve.F90:1078-1087).
pub(crate) fn bilu_level(p: &ParamMap) -> Option<usize> {
    let name = normalize_key(p.string("linear system preconditioning")?);
    if !name.starts_with("bilu") {
        return None;
    }
    let k = name.as_bytes().get(4).map_or(-1, |&c| i32::from(c) - i32::from(b'0'));
    Some(if (0..=9).contains(&k) { k as usize } else { 0 })
}

/// `CRS_BlockDiagonal(A, B, Blocks)` (CRSMatrix.F90:3002): the entries of
/// `a` whose row and column belong to the same of `blocks` components.
pub(crate) fn block_diagonal(a: &CsrMatrix<f64>, blocks: usize) -> CsrMatrix<f64> {
    let mut b = CsrMatrix { nrows: a.nrows, ncols: a.ncols, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
    for i in 0..a.nrows {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let c = a.col_idx[k] as usize;
            if (i + 1) % blocks == (c + 1) % blocks {
                b.col_idx.push(a.col_idx[k]);
                b.values.push(a.values[k]);
            }
        }
        b.row_ptr.push(b.values.len());
    }
    b
}

/// `PRECOND_BILUn` (IterSolve.F90:1277-1290): ILU(`level`) of the block
/// diagonal of `a` with `dofs` blocks, applied as the preconditioner of `a`.
pub(crate) fn bilu(
    a: &CsrMatrix<f64>,
    dofs: usize,
    level: usize,
) -> Result<Box<dyn Preconditioner<f64> + Send>, SolverError> {
    Ok(Box::new(Ilu::ilu_n(&block_diagonal(a, dofs), level)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `to_interleaved` is the inverse of `to_complex`, with the diagonal.
    #[test]
    fn interleaved_round_trip() {
        let c = CsrMatrix::<Complex64> {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 2, 3],
            col_idx: vec![0, 1, 1],
            values: vec![Complex64::new(2.0, -3.0), Complex64::new(-1.0, 0.5), Complex64::new(4.0, 1.0)],
            diag: Some(vec![0, 2]),
        };
        let a = to_interleaved(&c);
        assert_eq!(a.row_ptr, vec![0, 4, 8, 10, 12]);
        assert_eq!(a.values[..4], [2.0, 3.0, -1.0, -0.5]);
        assert_eq!(a.values[4..8], [-3.0, 2.0, 0.5, -1.0]);
        assert_eq!(a.diag, Some(vec![0, 5, 8, 11]));
        let back = to_complex(&a).unwrap();
        assert_eq!(back.values, c.values);
        assert_eq!(back.col_idx, c.col_idx);
        assert_eq!(back.diag, c.diag);
    }

    /// `max_i |Ax - b|_i / (|A||x| + |b|)_i` on a 2 x 2 system; a row with
    /// `|A||x| + |b|` below AEPS and no residual counts zero.
    #[test]
    fn componentwise_backward_error_is_the_worst_row() {
        let a = CsrMatrix::<f64> {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 2, 3],
            col_idx: vec![0, 1, 1],
            values: vec![2.0, -1.0, 4.0],
            diag: None,
        };
        let x = [1.0, 1.0];
        // Ax = [1, 4]; b = [1, 3]: rows 0 / (3 + 1) and 1 / (4 + 3).
        let e = componentwise_backward_error(&a, &x, &[1.0, 3.0]);
        assert_eq!(e, 1.0 / 7.0);
        let z = CsrMatrix::<f64> {
            nrows: 1,
            ncols: 1,
            row_ptr: vec![0, 1],
            col_idx: vec![0],
            values: vec![0.0],
            diag: None,
        };
        assert_eq!(componentwise_backward_error(&z, &[1.0], &[0.0]), 0.0);
        assert_eq!(componentwise_backward_error(&z, &[1.0], &[1.0]), 1.0);
    }
}

/// CheckLinearSolverOptions (MainUtils.F90:164): a serial `Linear System
/// Solver = Direct` without `Linear System Direct Method` is UMFPACK in
/// Elmer (HAVE_UMFPACK), not `banded`. The section with the method added,
/// or `None` when nothing changes.
pub(crate) fn serial_direct_default(p: &ParamMap, comm: &dyn Communicator) -> Option<ParamMap> {
    let direct = p.string("linear system solver").is_some_and(|s| normalize_key(s) == "direct");
    if comm.size() > 1 || !direct || p.string("linear system direct method").is_some() {
        return None;
    }
    let mut q = p.clone();
    q.insert("Linear System Direct Method", elmer_core::Value::Str("umfpack".into()));
    Some(q)
}
