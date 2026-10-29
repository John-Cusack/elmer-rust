//! `BlockMatrixPrec` (BlockSolve.F90:3777): one block Jacobi or block
//! Gauss-Seidel sweep, each diagonal block solved from a zero initial guess.

use std::sync::Mutex;

use elmer_core::{CsrMatrix, ParamMap, Preconditioner};

use crate::{BlockError, BlockSystem, InnerSolver};

/// A tailored preconditioning matrix solved in place of a diagonal block
/// (`SubMatrix(i,i) % PrecMat`, here the `Block Matrix Schur Variable`
/// matrix of the last block).
#[derive(Clone, Debug)]
pub struct PrecBlock<'a> {
    pub matrix: &'a CsrMatrix<f64>,
    /// `DoPrecScaling` (BlockSolve.F90:3934-3961): the matrix belongs to the
    /// unscaled system while the outer iteration runs in the scaled one: the
    /// block's `DiagScaling` entries, per row, and `RhsScaling`.
    pub unscale: Option<(Vec<f64>, f64)>,
}

/// The keywords a sweep reads, resolved once per block solve.
#[derive(Clone, Debug)]
pub struct PrecSettings<'a> {
    /// Solution order of the blocks, 0-based (`Block Order`, default 1..n).
    pub order: Vec<usize>,
    /// `Block Gauss-Seidel`: later blocks see the updated right-hand side.
    pub gauss_seidel: bool,
    /// `Dummy block i`: the block passes its right-hand side through.
    pub dummy: Vec<bool>,
    /// The section as each block's inner solve sees it (`block:` then `block
    /// NN:` pushed, with BlockMatrixPrec's own `ListAddLogical` settings).
    pub inner: Vec<ParamMap>,
    /// `passive[k][i]`: `Block Gauss-Seidel Passive ki`, the update of block
    /// `k` from block `i` is skipped.
    pub passive: Vec<Vec<bool>>,
    /// Per block, the matrix solved instead of `A(i,i)`, if any.
    pub prec: Vec<Option<PrecBlock<'a>>>,
}

/// One inner solve of a sweep.
#[derive(Clone, Debug, PartialEq)]
pub struct InnerRecord {
    /// The block solved, 0-based.
    pub block: usize,
    pub iterations: usize,
    pub residual: f64,
}

/// `u = M^-1 v` for the block system: blocks in [`PrecSettings::order`],
/// each solving `A(i,i) u(i) = v(i)` (Gauss-Seidel: `v(i)` less the
/// already solved blocks' coupling). Inner failures are kept and reported
/// by [`take_error`](Self::take_error) since `apply` cannot fail.
pub struct BlockPreconditioner<'a> {
    sys: &'a BlockSystem,
    settings: PrecSettings<'a>,
    state: Mutex<State<'a>>,
}

struct State<'a> {
    inner: Box<dyn InnerSolver + 'a>,
    records: Vec<InnerRecord>,
    applications: usize,
    error: Option<BlockError>,
}

impl<'a> BlockPreconditioner<'a> {
    pub fn new(sys: &'a BlockSystem, settings: PrecSettings<'a>, inner: Box<dyn InnerSolver + 'a>) -> Self {
        BlockPreconditioner {
            sys,
            settings,
            state: Mutex::new(State {
                inner,
                records: Vec::new(),
                applications: 0,
                error: None,
            }),
        }
    }

    /// The first inner failure, if any.
    pub fn take_error(&self) -> Option<BlockError> {
        self.lock().error.take()
    }

    /// Every inner solve so far, in order.
    pub fn records(&self) -> Vec<InnerRecord> {
        self.lock().records.clone()
    }

    /// Number of sweeps so far (`TotMatrix % NoIters`).
    pub fn applications(&self) -> usize {
        self.lock().applications
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, State<'a>> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn sweep(&self, st: &mut State<'a>, v: &[f64], u: &mut [f64]) -> Result<(), BlockError> {
        let sys = self.sys;
        let s = &self.settings;
        let n = sys.size();
        u[..n].fill(0.0);
        let mut vtmp = if s.gauss_seidel { Some(v[..n].to_vec()) } else { None };
        let mut rtmp = Vec::new();
        for (j, &i) in s.order.iter().enumerate() {
            let ri = sys.range(i);
            if s.dummy[i] {
                u[ri.clone()].copy_from_slice(&v[ri]);
                continue;
            }
            let b = match &vtmp {
                Some(t) => &t[ri.clone()],
                None => &v[ri.clone()],
            };
            let stats = match &s.prec[i] {
                Some(PrecBlock {
                    matrix,
                    unscale: Some((d, rhs)),
                }) => {
                    // x and b to the unscaled system and the solution back.
                    let x = &mut u[ri.clone()];
                    for (xk, dk) in x.iter_mut().zip(d) {
                        *xk = *xk * dk * rhs;
                    }
                    let bt: Vec<f64> = b.iter().zip(d).map(|(bk, dk)| bk / dk * rhs).collect();
                    let stats = st.inner.solve(i, &s.inner[i], matrix, &bt, x)?;
                    for (xk, dk) in x.iter_mut().zip(d) {
                        *xk /= dk * rhs;
                    }
                    stats
                }
                Some(PrecBlock { matrix, unscale: None }) => {
                    st.inner.solve(i, &s.inner[i], matrix, b, &mut u[ri.clone()])?
                }
                None => {
                    let Some(a) = sys.block(i, i) else { continue };
                    st.inner.solve(i, &s.inner[i], a, b, &mut u[ri.clone()])?
                }
            };
            st.records.push(InnerRecord {
                block: i,
                iterations: stats.iterations,
                residual: stats.residual,
            });
            let Some(vt) = vtmp.as_mut() else { continue };
            for &k in &s.order[j + 1..] {
                if s.passive[k][i] {
                    continue;
                }
                let Some(aki) = sys.block(k, i) else { continue };
                rtmp.clear();
                rtmp.resize(aki.nrows, 0.0);
                elmer_sparse::spmv(aki, &u[ri.clone()], &mut rtmp);
                for (t, r) in vt[sys.range(k)].iter_mut().zip(&rtmp) {
                    *t -= r;
                }
            }
        }
        Ok(())
    }
}

impl Preconditioner<f64> for BlockPreconditioner<'_> {
    fn apply(&self, v: &[f64], u: &mut [f64]) {
        let mut st = self.lock();
        st.applications += 1;
        if st.error.is_some() {
            u[..self.sys.size()].copy_from_slice(&v[..self.sys.size()]);
            return;
        }
        if let Err(e) = self.sweep(&mut st, v, u) {
            st.error = Some(e);
        }
    }
}
