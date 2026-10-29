//! One parallel iterative solve as `SolveLinearSystem` runs it with
//! `Linear System Solver = Iterative`: zero-rhs shortcut, diagonal scaling
//! (`ScaleLinearSystem`, MatrixScaling.F90), the split matrix, a block-local
//! preconditioner on the inside matrix and HUTI on owned dofs
//! (`SolveHutiter`, SParIterSolver.F90:2511), then back-scaling and the ghost
//! update (`ExchangeResult`).

use elmer_core::{
    CsrMatrix, IdentityPrecond, IterativeSolver, LinalgError, Params, Preconditioner, SolveStats,
    SolveStatus,
};
use elmer_krylov::{Huti, HutiConfig, OutputSink};
use elmer_sparse::Ilu;

use crate::exchange::Exchange;
use crate::layout::DofLayout;
use crate::matrix::DistCsr;
use crate::{agree, DistError};

/// `Linear System Preconditioning`, applied to each rank's inside matrix only
/// (block Jacobi across ranks), as the parallel IterSolve does.
#[derive(Clone, Debug, PartialEq)]
pub enum DistPrecond {
    None,
    Diagonal,
    /// ILU(n), `n` in 0..=9.
    Ilu(usize),
    /// ILUT with `Linear System ILUT Tolerance`.
    Ilut(f64),
}

/// What one distributed solve reads from the solver section.
#[derive(Clone, Debug)]
pub struct DistSolveOptions {
    pub huti: HutiConfig,
    pub precond: DistPrecond,
    /// `Linear System Scaling` with the `diagonal` method (default on).
    pub scaling: bool,
}

impl DistSolveOptions {
    pub fn new(huti: HutiConfig) -> Self {
        DistSolveOptions {
            huti,
            precond: DistPrecond::None,
            scaling: true,
        }
    }

    /// IterSolve's preconditioner choice (IterSolve.F90:1050-1100) and
    /// SolveLinearSystem's scaling switch. Choices this crate does not
    /// implement are refused rather than replaced.
    pub fn from_params(p: &dyn Params) -> Result<Self, LinalgError> {
        let huti = HutiConfig::from_params(p)?;
        let name = p
            .string("Linear System Preconditioning")
            .map_or("none".to_string(), |s| s.trim().to_ascii_lowercase());
        if p.logical("Linear System Symmetric ILU") == Some(true) {
            return Err(LinalgError::Unsupported(
                "Linear System Symmetric ILU in parallel".into(),
            ));
        }
        let precond = match name.as_str() {
            "none" => DistPrecond::None,
            "diagonal" => DistPrecond::Diagonal,
            "ilut" => DistPrecond::Ilut(p.real("Linear System ILUT Tolerance").unwrap_or(0.0)),
            s if s.starts_with("ilu") => {
                let n = p.integer("Linear System ILU Order").unwrap_or_else(|| {
                    s.as_bytes()
                        .get(3)
                        .map_or(-1, |&c| i32::from(c) - i32::from(b'0'))
                });
                DistPrecond::Ilu(if (0..=9).contains(&n) { n as usize } else { 0 })
            }
            s if s.starts_with("bilu")
                || s.starts_with("vanka")
                || matches!(
                    s,
                    "multigrid" | "auxiliary space solver" | "slave" | "circuit"
                ) =>
            {
                return Err(LinalgError::Unsupported(format!(
                    "Linear System Preconditioning = {s} in parallel"
                )))
            }
            // IterSolve warns and runs without a preconditioner.
            _ => DistPrecond::None,
        };
        // ScaleLinearSystem (MatrixScaling.F90:84-108).
        let mut scaling = p.logical("Linear System Scaling").unwrap_or(true);
        if scaling {
            if p.logical("Linear System Row Equilibration") == Some(true) {
                return Err(LinalgError::Unsupported(
                    "Linear System Row Equilibration in parallel".into(),
                ));
            }
            let method = p
                .string("Linear System Scaling Method")
                .map_or("diagonal".to_string(), |s| s.trim().to_ascii_lowercase());
            match method.as_str() {
                "diagonal" => {}
                "none" => scaling = false,
                s => {
                    return Err(LinalgError::Unsupported(format!(
                        "Linear System Scaling Method = {s} in parallel"
                    )))
                }
            }
        }
        Ok(DistSolveOptions {
            huti,
            precond,
            scaling,
        })
    }
}

/// Solve the distributed system whose rank-local pieces are `a` and `b`
/// (local assembled matrix and rhs, partial sums at shared dofs), starting
/// from `x` (consistent across copies). On return `x` holds the solution at
/// every local dof, ghosts included. Collective.
///
/// HUTI prints its residual lines on rank 0 only, to `sink` when given.
/// `MaxIterations` and `Diverged` come back in the status; whether that is
/// fatal (`Linear System Abort Not Converged`) is the caller's decision.
pub fn dist_solve(
    comm: &dyn Exchange,
    layout: &DofLayout,
    a: &CsrMatrix<f64>,
    b: &[f64],
    x: &mut [f64],
    opts: &DistSolveOptions,
    sink: Option<OutputSink>,
) -> Result<SolveStats, DistError> {
    let n = layout.n_local();
    agree(
        comm,
        if b.len() == n && x.len() == n {
            Ok(())
        } else {
            Err(format!(
                "rhs {} and x {} for {n} local dofs",
                b.len(),
                x.len()
            ))
        },
    )?;

    // SolveCore.F90:1514-1520: the norm of the unsummed local pieces.
    let mut bb = [b.iter().map(|v| v * v).sum::<f64>()];
    comm.sum_f64(&mut bb);
    if bb[0].sqrt() <= f64::MIN_POSITIVE {
        x.fill(0.0);
        return Ok(SolveStats {
            iterations: 0,
            residual: 0.0,
            status: SolveStatus::Converged,
        });
    }

    let mut a = a.clone();
    let mut b = b.to_vec();
    let (diag, rhs_scaling) = if opts.scaling {
        scale_diagonal(comm, layout, &mut a, &mut b, x)
    } else {
        (vec![1.0; n], 1.0)
    };

    let m = DistCsr::from_local(comm, layout, &a)?;
    let b_own = layout.sum_to_owned(comm, &b);
    let mut x_own = layout.gather_owned(x);

    let prec: Box<dyn Preconditioner<f64>> = match &opts.precond {
        DistPrecond::None => Box::new(IdentityPrecond),
        DistPrecond::Diagonal => Box::new(Ilu::diagonal(m.inside())),
        DistPrecond::Ilu(k) => Box::new(agree(
            comm,
            Ilu::ilu_n(m.inside(), *k).map_err(|e| e.to_string()),
        )?),
        DistPrecond::Ilut(tol) => Box::new(agree(
            comm,
            Ilu::ilut(m.inside(), *tol).map_err(|e| e.to_string()),
        )?),
    };

    let mut huti: Huti = Huti::new(opts.huti.clone());
    if let Some(s) = sink {
        huti.set_output(s);
    }
    let op = m.operator(comm);
    let stats = huti.solve(&op, prec.as_ref(), &b_own, &mut x_own, comm)?;

    layout.owned_to_local(comm, &x_own, x);
    for (xi, d) in x.iter_mut().zip(&diag) {
        *xi *= d * rhs_scaling;
    }
    Ok(stats)
}

/// `ScaleLinearSystemDiagonal` with `RhsScaling`: `A <- D A D`, `b <- D b / |D b|`,
/// `x <- x / (D |D b|)` with `D = 1/sqrt(|diag A|)` from the summed diagonal.
/// Returns `D` per local dof and the rhs scaling.
fn scale_diagonal(
    comm: &dyn Exchange,
    layout: &DofLayout,
    a: &mut CsrMatrix<f64>,
    b: &mut [f64],
    x: &mut [f64],
) -> (Vec<f64>, f64) {
    let n = a.nrows;
    let diag_of = |a: &CsrMatrix<f64>, i: usize| {
        (a.row_ptr[i]..a.row_ptr[i + 1])
            .find(|&k| a.col_idx[k] as usize == i)
            .map(|k| a.values[k])
    };
    let mut d: Vec<f64> = (0..n).map(|i| diag_of(a, i).unwrap_or(0.0)).collect();
    layout.sum_shared(comm, &mut d);

    // A (numerically) zero diagonal anywhere: those rows use the local absolute
    // row sum instead, and the vector is summed again.
    let mut zero = [if d.iter().any(|v| v.abs() <= f64::MIN_POSITIVE) {
        1.0
    } else {
        0.0
    }];
    comm.max_f64(&mut zero);
    if zero[0] > f64::MIN_POSITIVE {
        for (i, di) in d.iter_mut().enumerate() {
            *di = if di.abs() <= f64::MIN_POSITIVE {
                a.values[a.row_ptr[i]..a.row_ptr[i + 1]]
                    .iter()
                    .map(|v| v.abs())
                    .sum()
            } else {
                diag_of(a, i).unwrap_or(0.0)
            };
        }
        layout.sum_shared(comm, &mut d);
    }
    for di in d.iter_mut() {
        *di = if di.abs() > f64::MIN_POSITIVE {
            1.0 / di.abs().sqrt()
        } else {
            1.0
        };
    }

    for i in 0..n {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            a.values[k] *= d[i] * d[a.col_idx[k] as usize];
        }
    }
    for (bi, di) in b.iter_mut().zip(&d) {
        *bi *= di;
    }
    let summed = layout.sum_to_owned(comm, b);
    let mut s = [summed.iter().map(|v| v * v).sum::<f64>()];
    comm.sum_f64(&mut s);
    let bnorm = s[0].sqrt();
    let rhs_scaling = if bnorm < f64::MIN_POSITIVE.sqrt() {
        1.0
    } else {
        for bi in b.iter_mut() {
            *bi /= bnorm;
        }
        bnorm
    };
    for (xi, di) in x.iter_mut().zip(&d) {
        *xi /= di * rhs_scaling;
    }
    (d, rhs_scaling)
}
