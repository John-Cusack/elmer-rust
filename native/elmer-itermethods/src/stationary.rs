//! Jacobi (`IterativeMethods.F90:296`), Richardson (`:404`) and symmetric
//! Gauss-Seidel (`:178`). All three ignore the preconditioner and work on the
//! plain `||r|| / ||b||` residual; they read matrix entries, so they are built
//! from the CSR matrix and are real-valued only (the Fortran has no complex
//! variant).

use crate::config::IterMethodsConfig;
use crate::ops::{check_dims, diverged, nrm2, residual};
use crate::{finish, History};
use elmer_core::{
    Communicator, CsrMatrix, IterativeSolver, LinalgError, LinearOperator, Preconditioner,
    SolveStats,
};

/// Shared loop of the three methods; `update(k, x, r)` is one sweep, `k` from 1.
fn run(
    cfg: &IterMethodsConfig,
    a: &dyn LinearOperator<f64>,
    b: &[f64],
    x: &mut [f64],
    comm: &dyn Communicator,
    history: &mut History,
    mut update: impl FnMut(usize, &mut [f64], &[f64]),
) -> SolveStats {
    let n = b.len();
    let mut r = vec![0.0; n];
    residual(a, x, b, &mut r);
    let bnorm = nrm2(comm, b);
    let rnorm = nrm2(comm, &r);

    let mut res = rnorm / bnorm;
    let mut conv = res < cfg.tolerance;
    let mut div = diverged(res, cfg.divergence_limit);
    let mut iters = 0;
    if !(conv || div) {
        for k in 1..=cfg.max_iterations {
            iters = k;
            update(k, x, &r);
            residual(a, x, b, &mut r);
            res = nrm2(comm, &r) / bnorm;
            history.push(res);
            conv = res < cfg.tolerance;
            div = diverged(res, cfg.divergence_limit);
            if conv || div {
                break;
            }
        }
    }
    finish(iters, res, conv, div, None)
}

/// Jacobi iteration `x += D^-1 (b - A x)`.
#[derive(Clone, Debug)]
pub struct Jacobi {
    cfg: IterMethodsConfig,
    diag: Vec<f64>,
    history: History,
}

impl Jacobi {
    /// Captures the diagonal of `a` (`A % Values(A % diag(j))`). Errors if a
    /// row has no stored diagonal entry. Call [`refresh`](Self::refresh) if the
    /// matrix values change between solves.
    pub fn new(cfg: &IterMethodsConfig, a: &CsrMatrix<f64>) -> Result<Self, LinalgError> {
        let mut j = Jacobi {
            cfg: cfg.clone(),
            diag: Vec::new(),
            history: History::default(),
        };
        j.refresh(a)?;
        Ok(j)
    }

    pub fn refresh(&mut self, a: &CsrMatrix<f64>) -> Result<(), LinalgError> {
        self.diag = diagonal(a)?;
        Ok(())
    }

    /// Relative residual after each iteration of the last solve.
    pub fn history(&self) -> &[f64] {
        &self.history.0
    }
}

impl IterativeSolver<f64> for Jacobi {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<f64>,
        _m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims("Jacobi", a, b, x)?;
        if self.diag.len() != b.len() {
            return Err(LinalgError::Dim("Jacobi: matrix changed size".into()));
        }
        self.history.0.clear();
        let diag = &self.diag;
        Ok(run(&self.cfg, a, b, x, comm, &mut self.history, |_, x, r| {
            for ((xi, ri), d) in x.iter_mut().zip(r).zip(diag) {
                *xi += ri / d;
            }
        }))
    }
}

/// Richardson iteration preconditioned with the lumped (row-sum) matrix. The
/// first sweep discards the initial guess (`x = b / M`), as the Fortran does.
#[derive(Clone, Debug)]
pub struct Richardson {
    cfg: IterMethodsConfig,
    lumped: Vec<f64>,
    history: History,
}

impl Richardson {
    pub fn new(cfg: &IterMethodsConfig, a: &CsrMatrix<f64>) -> Self {
        Richardson {
            cfg: cfg.clone(),
            lumped: lump(a),
            history: History::default(),
        }
    }

    pub fn refresh(&mut self, a: &CsrMatrix<f64>) {
        self.lumped = lump(a);
    }

    pub fn history(&self) -> &[f64] {
        &self.history.0
    }
}

impl IterativeSolver<f64> for Richardson {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<f64>,
        _m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims("Richardson", a, b, x)?;
        if self.lumped.len() != b.len() {
            return Err(LinalgError::Dim("Richardson: matrix changed size".into()));
        }
        self.history.0.clear();
        let lumped = &self.lumped;
        Ok(run(&self.cfg, a, b, x, comm, &mut self.history, |k, x, r| {
            for (i, xi) in x.iter_mut().enumerate() {
                if k == 1 {
                    *xi = b[i] / lumped[i];
                } else {
                    *xi += r[i] / lumped[i];
                }
            }
        }))
    }
}

/// Symmetric Gauss-Seidel with over-relaxation `omega`
/// (`SGS Overrelaxation Factor`, default 1.8).
#[derive(Debug)]
pub struct Sgs<'a> {
    cfg: IterMethodsConfig,
    a: &'a CsrMatrix<f64>,
    history: History,
}

impl<'a> Sgs<'a> {
    pub fn new(cfg: &IterMethodsConfig, a: &'a CsrMatrix<f64>) -> Self {
        Sgs {
            cfg: cfg.clone(),
            a,
            history: History::default(),
        }
    }

    pub fn history(&self) -> &[f64] {
        &self.history.0
    }
}

impl IterativeSolver<f64> for Sgs<'_> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<f64>,
        _m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims("SGS", a, b, x)?;
        let csr = self.a;
        if csr.nrows != b.len() {
            return Err(LinalgError::Dim("SGS: matrix changed size".into()));
        }
        let diag = diagonal(csr)?;
        let omega = self.cfg.sgs_omega;
        self.history.0.clear();
        let row_update = |i: usize, x: &mut [f64]| {
            let mut s = 0.0;
            for k in csr.row_ptr[i]..csr.row_ptr[i + 1] {
                s += x[csr.col_idx[k] as usize] * csr.values[k];
            }
            x[i] += omega * (b[i] - s) / diag[i];
        };
        Ok(run(&self.cfg, a, b, x, comm, &mut self.history, |_, x, _| {
            for i in 0..b.len() {
                row_update(i, x);
            }
            for i in (0..b.len()).rev() {
                row_update(i, x);
            }
        }))
    }
}

/// `A[i,i]` for every row, via `diag` if present, else by searching the row.
fn diagonal(a: &CsrMatrix<f64>) -> Result<Vec<f64>, LinalgError> {
    (0..a.nrows)
        .map(|i| match &a.diag {
            Some(d) => Ok(a.values[d[i]]),
            None => (a.row_ptr[i]..a.row_ptr[i + 1])
                .find(|&k| a.col_idx[k] as usize == i)
                .map(|k| a.values[k])
                .ok_or(LinalgError::Singular(i)),
        })
        .collect()
}

/// Row sums, the lumped mass `M(i)` of `itermethod_richardson`.
fn lump(a: &CsrMatrix<f64>) -> Vec<f64> {
    (0..a.nrows)
        .map(|i| a.values[a.row_ptr[i]..a.row_ptr[i + 1]].iter().sum())
        .collect()
}
