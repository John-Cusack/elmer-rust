//! `BlockStandardIter` (BlockSolve.F90:4239) with `BlockUpdateRhs` (:2894):
//! the block Jacobi or Gauss-Seidel outer iteration of `Block Preconditioner
//! = False`, each block solving for a correction of its residual.

use elmer_core::{ParamMap, Params, SolveStats, SolveStatus};

use crate::solve::comp_norm;
use crate::{BlockError, BlockSystem, InnerRecord, InnerSolver, PrecBlock};

/// What the outer iteration did.
pub(crate) struct StandardOutcome {
    pub stats: SolveStats,
    pub records: Vec<InnerRecord>,
}

/// Runs the iteration on the block-ordered `x`. `outer` is the section under
/// `outer:` (iteration count, tolerance, `Block Gauss-Seidel`, `Block Order`,
/// `Block Relaxation Factor`); `inner[0][i]` the one block `i` is solved
/// with in the first iteration (`block:` and `block NN:` pushed),
/// `inner[1][i]` in later ones (with `No Precondition Recompute`);
/// `prec[i]` a matrix solved in place of block `i`'s.
#[allow(clippy::too_many_arguments)]
pub(crate) fn standard_iter(
    outer: &ParamMap,
    order: &[usize],
    inner: [&[ParamMap]; 2],
    prec: &[Option<PrecBlock<'_>>],
    sys: &BlockSystem,
    b: &[f64],
    x: &mut [f64],
    solver: &mut dyn InnerSolver,
    log: &dyn elmer_core::Logger,
) -> Result<StandardOutcome, BlockError> {
    let gs = outer.logical("block gauss-seidel") == Some(true);
    let max_iter = outer.integer("linear system max iterations").unwrap_or(0).max(0) as usize;
    let min_iter = outer.integer("linear system min iterations").unwrap_or(0).max(0) as usize;
    let tol = outer.real("linear system convergence tolerance").unwrap_or(0.0);
    let relax = outer.real("block relaxation factor").unwrap_or(1.0);
    let nb = sys.nblocks();

    let mut records = Vec::new();
    let mut rhs: Vec<Vec<f64>> = (0..nb).map(|i| vec![0.0; sys.range(i).len()]).collect();
    let mut iterations = 0;
    let mut max_change = 0.0;
    let mut converged = false;
    for iter in 1..=max_iter {
        iterations = iter;
        if gs {
            log.info("BlockStandardIter", &format!("Block Gauss-Seidel iteration: {iter}"), 5);
        } else {
            log.info("BlockStandardIter", &format!("Block Jacobi iteration: {iter}"), 5);
            for (i, r) in rhs.iter_mut().enumerate() {
                update_rhs(sys, i, b, x, r);
            }
        }
        max_change = 0.0f64;
        let mut tot_norm = 0.0;
        for &i in order {
            if gs {
                update_rhs(sys, i, b, x, &mut rhs[i]);
            }
            let a = match &prec[i] {
                Some(p) => p.matrix,
                None => sys.block(i, i).expect("diagonal blocks are present"),
            };
            let mut dx = vec![0.0; a.nrows];
            let view = &inner[usize::from(iter > 1)][i];
            let stats = solver.solve(i, view, a, &rhs[i], &mut dx)?;
            records.push(InnerRecord {
                block: i,
                iterations: stats.iterations,
                residual: stats.residual,
            });
            let xi = &mut x[sys.range(i)];
            if iter > 1 {
                for (v, d) in xi.iter_mut().zip(&dx) {
                    *v += relax * d;
                }
            } else {
                for (v, d) in xi.iter_mut().zip(&dx) {
                    *v += d;
                }
            }
            let dxnorm = comp_norm(&dx);
            let xnorm = comp_norm(xi);
            let change = dxnorm / xnorm;
            log.info(
                "BlockStandardIter",
                &format!("Block {} norms: {xnorm:12.3E}{change:12.3E}", i + 1),
                5,
            );
            tot_norm += xnorm;
            // Fortran MAX keeps the first argument against a NaN, as f64::max does.
            max_change = max_change.max(change);
        }
        log.info(
            "BlockStandardIter",
            &format!("Sum of norms: {tot_norm:12.3E}{max_change:12.3E}"),
            4,
        );
        if max_change < tol && iter >= min_iter {
            log.info("BlockStandardIter", &format!("Converged after iterations: {iter}"), 5);
            converged = true;
            break;
        }
    }
    Ok(StandardOutcome {
        stats: SolveStats {
            iterations,
            residual: max_change,
            status: if converged {
                SolveStatus::Converged
            } else {
                SolveStatus::MaxIterations
            },
        },
        records,
    })
}

/// `BlockUpdateRhs` for row block `i`: `b(i) - sum_j A(i,j) x(j)`, the
/// off-diagonal blocks first in column order, then the diagonal one.
fn update_rhs(sys: &BlockSystem, i: usize, b: &[f64], x: &[f64], rhs: &mut [f64]) {
    rhs.copy_from_slice(&b[sys.range(i)]);
    let mut t = vec![0.0; rhs.len()];
    let cols = (0..sys.nblocks()).filter(|&j| j != i).chain(std::iter::once(i));
    for j in cols {
        let Some(a) = sys.block(i, j) else { continue };
        elmer_sparse::spmv(a, &x[sys.range(j)], &mut t);
        for (r, v) in rhs.iter_mut().zip(&t) {
            *r -= v;
        }
    }
}
