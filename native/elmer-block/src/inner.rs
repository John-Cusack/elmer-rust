//! The solve of one diagonal block: `SolveLinearSystem` (SolveCore.F90:1186)
//! as `BlockMatrixPrec` calls it, under the `block NN:` namespace.

use elmer_core::{
    normalize_key, Communicator, CsrMatrix, DirectSolver, IdentityPrecond, IterativeSolver, LinalgError, Logger,
    ParamMap, Params, Preconditioner, SolveStatus,
};
use elmer_direct::{kind_for_method, BandCholesky, BandLu, DirectKind, FaerCholesky, FaerLu};
use elmer_itermethods::{BiCgStabL, Gcr, Idrs, IterMethodsConfig, Jacobi, Richardson, Sgs};
use elmer_krylov::{Huti, HutiConfig};
use elmer_sparse::{CsrOp, Ilu};

use crate::BlockError;

/// What one inner solve did.
#[derive(Clone, Debug, PartialEq)]
pub struct InnerStats {
    /// Iterations of an iterative method; 0 for a direct solve or a zero
    /// right-hand side.
    pub iterations: usize,
    /// The last monitored residual (0 when not iterated).
    pub residual: f64,
}

/// Solves `A x = b` for one diagonal block. `params` is the solver section
/// as seen under the block's namespaces; `x` holds the initial guess. The
/// matrix of a block does not change during one [`crate::block_solve`], so
/// an implementation may keep factorizations and preconditioners per block.
pub trait InnerSolver: Send {
    fn solve(
        &mut self,
        block: usize,
        params: &ParamMap,
        a: &CsrMatrix<f64>,
        b: &[f64],
        x: &mut [f64],
    ) -> Result<InnerStats, BlockError>;
}

/// The serial `SolveLinearSystem` path for a block: `Linear System Skip
/// Scaling` / `Linear System Scaling` (default on) around the solve, the zero
/// right-hand side shortcut, `Linear System Solver = Iterative` with the
/// IterSolve methods and the `diagonal`, `ilu[n]` and `ilut`
/// preconditioners, or `Direct` with the DirectSolve methods. Factorizations
/// and preconditioners are built once per block and reused, which is what
/// `Linear System Refactorize = False` asks for and gives the same numbers
/// as refactorizing the unchanged block.
pub struct ElmerInner<'a> {
    log: &'a dyn Logger,
    comm: &'a dyn Communicator,
    /// `Global Abort Not Converged` (Simulation, default true).
    global_abort: bool,
    cache: Vec<Option<Cached>>,
}

enum Cached {
    Precond(Box<dyn Preconditioner<f64> + Send>),
    Direct(Box<dyn DirectSolver<f64> + Send>),
}

impl<'a> ElmerInner<'a> {
    pub fn new(log: &'a dyn Logger, comm: &'a dyn Communicator, global_abort: bool) -> Self {
        ElmerInner {
            log,
            comm,
            global_abort,
            cache: Vec::new(),
        }
    }

    fn solve_scaled(
        &mut self,
        block: usize,
        p: &ParamMap,
        a: &CsrMatrix<f64>,
        b: &[f64],
        x: &mut [f64],
    ) -> Result<InnerStats, BlockError> {
        let kind = p.string("linear system solver").map(normalize_key).ok_or_else(|| {
            BlockError::fatal(
                "SolveLinearSystem",
                "Give \"Linear System Solver\", e.g. \"iterative\" or \"direct\"",
            )
        })?;
        if self.cache.len() <= block {
            self.cache.resize_with(block + 1, || None);
        }
        match kind.as_str() {
            "iterative" => {
                // IterSolve.F90:1108: `No Precondition Recompute` skips the
                // factorization; an ILU that was never built then leaves
                // CRS_LUSolve dividing by the diagonal (CRSMatrix.F90:5372).
                let no_recompute = p.logical("no precondition recompute") == Some(true);
                let fallback;
                let m: &dyn Preconditioner<f64> = match &self.cache[block] {
                    Some(Cached::Precond(m)) => &**m,
                    Some(Cached::Direct(_)) => {
                        return Err(BlockError::fatal("SolveLinearSystem", "block changed its solver kind"))
                    }
                    None if no_recompute && is_ilu(p) => {
                        fallback = DiagonalFallback(a);
                        &fallback
                    }
                    None => {
                        self.cache[block] = Some(Cached::Precond(precond(p, a, self.log)?));
                        let Some(Cached::Precond(m)) = &self.cache[block] else {
                            unreachable!()
                        };
                        &**m
                    }
                };
                let mut solver = iterative(p, Some(a))?;
                let stats = solver.solve(&CsrOp(a), m, b, x, self.comm)?;
                check_converged(p, &stats.status, self.global_abort)?;
                Ok(InnerStats {
                    iterations: stats.iterations,
                    residual: stats.residual,
                })
            }
            "direct" => {
                if self.cache[block].is_none() {
                    let mut d = direct(p)?;
                    d.factor(a)?;
                    self.cache[block] = Some(Cached::Direct(d));
                }
                let Some(Cached::Direct(d)) = &self.cache[block] else {
                    return Err(BlockError::fatal("SolveLinearSystem", "block changed its solver kind"));
                };
                let mut y = b.to_vec();
                d.solve_in_place(&mut y)?;
                x.copy_from_slice(&y);
                Ok(InnerStats {
                    iterations: 0,
                    residual: 0.0,
                })
            }
            "block" => Err(BlockError::Unsupported("nested block linear systems".into())),
            other => Err(BlockError::Unsupported(format!(
                "block {}: Linear System Solver = {other}",
                11 * (block + 1)
            ))),
        }
    }
}

impl InnerSolver for ElmerInner<'_> {
    fn solve(
        &mut self,
        block: usize,
        p: &ParamMap,
        a: &CsrMatrix<f64>,
        b: &[f64],
        x: &mut [f64],
    ) -> Result<InnerStats, BlockError> {
        if p.logical("linear system skip complex") != Some(true) && p.logical("linear system complex") == Some(true) {
            return Err(BlockError::Unsupported("complex blocks (Linear System Complex)".into()));
        }
        for k in [
            "linear system normalize guess",
            "nonlinear system acceleration",
            "linear system skip mask",
            "linear system use hypre",
            "linear system use rocalution",
        ] {
            if p.logical(k) == Some(true) {
                return Err(BlockError::Unsupported(format!("block solve option {k}")));
            }
        }
        // SolveCore.F90:1263-1268.
        let scale =
            p.logical("linear system skip scaling") != Some(true) && p.logical("linear system scaling").unwrap_or(true);
        // SolveCore.F90:1502-1520: b = 0 has only the trivial solution.
        let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        if bnorm <= f64::MIN_POSITIVE && p.logical("skip zero rhs test") != Some(true) {
            self.log.info("SolveLinearSystem", "Solution trivially zero!", 5);
            x.fill(0.0);
            return Ok(InnerStats {
                iterations: 0,
                residual: 0.0,
            });
        }
        if p.logical("linear system nullify guess") == Some(true) {
            x.fill(0.0);
        }
        if !scale {
            return self.solve_scaled(block, p, a, b, x);
        }
        // ScaleLinearSystem / BackScaleLinearSystem (SolveCore.F90:1540, :1714)
        // on copies: the block matrix and the caller's vectors stay as they were.
        let method = if p.logical("linear system row equilibration") == Some(true) {
            elmer_sparse::ScalingMethod::RowEquilibration
        } else {
            elmer_sparse::ScalingMethod::from_keyword(p.string("linear system scaling method").unwrap_or("diagonal"))?
        };
        let opts = elmer_sparse::ScaleOptions {
            method,
            apply: true,
            rhs_scaling: bnorm != 0.0,
        };
        let mut sa = a.clone();
        let mut sb = b.to_vec();
        let scaling = elmer_sparse::scale_system(&mut sa, Some(&mut sb[..]), Some(&mut *x), &opts);
        let r = self.solve_scaled(block, p, &sa, &sb, x);
        scaling.unscale(&mut sa, Some(&mut sb[..]), Some(&mut *x));
        r
    }
}

/// IterSolve.F90:1640-1660: a solve that hit its iteration limit is a
/// numerical error unless `Linear System Abort Not Converged = False`, and
/// fatal only under `Global Abort Not Converged` (default true).
pub(crate) fn check_converged(p: &dyn Params, status: &SolveStatus, global_abort: bool) -> Result<(), BlockError> {
    if *status == SolveStatus::MaxIterations
        && p.logical("linear system abort not converged").unwrap_or(true)
        && global_abort
    {
        return Err(BlockError::Numerical {
            caller: "IterSolve".into(),
            msg: "Too many iterations were needed.".into(),
        });
    }
    Ok(())
}

/// IterSolve's method dispatch (IterSolve.F90:724-780): the
/// IterativeMethods.F90 names go to `elmer-itermethods`, every other name to
/// the HUTI methods of `elmer-krylov` (an unknown one to BiCGStab).
/// `a` is the matrix the stationary methods read; without it they are refused.
pub(crate) fn iterative(
    p: &ParamMap,
    a: Option<&CsrMatrix<f64>>,
) -> Result<Box<dyn IterativeSolver<f64> + Send>, BlockError> {
    let m = p
        .string("linear system iterative method")
        .map(normalize_key)
        .unwrap_or_else(|| "bicgstab".into());
    Ok(match m.as_str() {
        "bicgstabl" | "gcr" | "idrs" => {
            let cfg = IterMethodsConfig::from_params(p)?;
            match m.as_str() {
                "bicgstabl" => Box::new(BiCgStabL::<f64>::new(&cfg)),
                "gcr" => Box::new(Gcr::<f64>::new(&cfg)),
                _ => Box::new(Idrs::<f64>::new(&cfg)),
            }
        }
        "jacobi" | "richardson" | "sgs" => {
            let a =
                a.ok_or_else(|| BlockError::Unsupported(format!("outer iterative method {m} of a block system")))?;
            let cfg = IterMethodsConfig::from_params(p)?;
            match m.as_str() {
                "jacobi" => Box::new(Jacobi::new(&cfg, a)?),
                "richardson" => Box::new(Richardson::new(&cfg, a)),
                _ => Box::new(OwnedSgs { cfg, a: a.clone() }),
            }
        }
        _ => Box::new(Huti::<f64>::new(HutiConfig::from_params(p)?)),
    })
}

/// [`Sgs`] borrows its matrix; this keeps a copy for the solve.
struct OwnedSgs {
    cfg: IterMethodsConfig,
    a: CsrMatrix<f64>,
}

impl IterativeSolver<f64> for OwnedSgs {
    fn solve(
        &mut self,
        a: &dyn elmer_core::LinearOperator<f64>,
        m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<elmer_core::SolveStats, LinalgError> {
        Sgs::new(&self.cfg, &self.a).solve(a, m, b, x, comm)
    }
}

/// IterSolve.F90:1050-1100 for the CRS preconditioners: `none`, `diagonal`,
/// `ilut` (`Linear System ILUT Tolerance`), `ilu[n]` (`Linear System ILU
/// Order` first, else the digit in the name, out of 0..9 -> 0) and incomplete
/// Cholesky with `Linear System Symmetric ILU`. An unknown name is warned
/// about and treated as `none`.
fn precond(
    p: &ParamMap,
    a: &CsrMatrix<f64>,
    log: &dyn Logger,
) -> Result<Box<dyn Preconditioner<f64> + Send>, BlockError> {
    if p.real("linear system ilu factor").is_some_and(|f| f > f64::EPSILON) {
        return Err(BlockError::Unsupported("Linear System ILU Factor".into()));
    }
    let name = normalize_key(p.string("linear system preconditioning").unwrap_or("none"));
    if name == "none" {
        return Ok(Box::new(IdentityPrecond));
    }
    if name == "diagonal" {
        return Ok(Box::new(Ilu::diagonal(a)));
    }
    if name == "ilut" {
        let tol = p.real("linear system ilut tolerance").unwrap_or(0.0);
        return Ok(Box::new(Ilu::ilut(a, tol)?));
    }
    if name.starts_with("ilu") {
        let mut level = match p.integer("linear system ilu order") {
            Some(k) => k,
            None => name.as_bytes().get(3).map_or(-1, |&c| i32::from(c) - i32::from(b'0')),
        };
        if !(0..=9).contains(&level) {
            level = 0;
        }
        if p.logical("linear system symmetric ilu") == Some(true) {
            let retry = if p.logical("linear system cholesky retry") == Some(true) {
                Some(p.real("linear system cholesky factor").ok_or_else(|| {
                    BlockError::Unsupported("Linear System Cholesky Retry without Linear System Cholesky Factor".into())
                })?)
            } else {
                None
            };
            return Ok(Box::new(Ilu::ic_n(a, level as usize, retry)?));
        }
        return Ok(Box::new(Ilu::ilu_n(a, level as usize)?));
    }
    if [
        "bilu",
        "multigrid",
        "vanka",
        "auxiliary space solver",
        "slave",
        "circuit",
    ]
    .iter()
    .any(|k| name.starts_with(k))
    {
        return Err(BlockError::Unsupported(format!("block preconditioner {name}")));
    }
    log.warn(
        "IterSolve",
        &format!("Unknown preconditioner type: {name}, feature disabled."),
    );
    Ok(Box::new(IdentityPrecond))
}

/// The preconditioner names CRS_LUPrecondition serves.
fn is_ilu(p: &ParamMap) -> bool {
    normalize_key(p.string("linear system preconditioning").unwrap_or("none")).starts_with("ilu")
}

/// CRS_LUSolve without factors: `b(i) / A(i,i)` where the diagonal is nonzero.
struct DiagonalFallback<'a>(&'a CsrMatrix<f64>);

impl Preconditioner<f64> for DiagonalFallback<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let a = self.0;
        for i in 0..a.nrows {
            let d = (a.row_ptr[i]..a.row_ptr[i + 1])
                .find(|&k| a.col_idx[k] as usize == i)
                .map_or(0.0, |k| a.values[k]);
            z[i] = if d != 0.0 { r[i] / d } else { r[i] };
        }
    }
}

/// `Linear System Direct Method`, default `banded` (DirectSolve.F90).
fn direct(p: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, BlockError> {
    Ok(
        match kind_for_method(p.string("linear system direct method").unwrap_or("banded"))? {
            DirectKind::BandLu => Box::new(BandLu::<f64>::new()),
            DirectKind::BandCholesky => Box::new(BandCholesky::<f64>::new()),
            DirectKind::SparseLu => Box::new(FaerLu::<f64>::new()),
            DirectKind::SparseCholesky => Box::new(FaerCholesky::<f64>::new()),
        },
    )
}
