//! `ArpackEigenSolve` (EigenSolve.F90:140): the real, undamped eigenproblem
//! `K x = lambda M x` of one assembled system.

use elmer_core::{CsrMatrix, DirectSolver, LinalgError, LinearOperator};
use num_complex::Complex64;

use crate::config::{ArpackWhich, EigenConfig, Normalize};
use crate::krylov::{krylov_schur, EigenOperator, KrylovParams};
use crate::scale::{check_residuals, scale_eigen_vectors};
use crate::sort::eigen_system_sorting;

/// The mass side of the problem.
#[derive(Copy, Clone, Debug)]
pub enum Mass<'a> {
    /// `M = I`: the standard problem `K x = lambda x`, run like a consistent
    /// mass (shift-invert).
    Identity,
    /// `Matrix % Lumped`: the diagonal of `M`. ARPACK mode 2, `OP = M^-1 K`;
    /// `Eigen System Shift` is ignored, as in Elmer.
    Lumped(&'a [f64]),
    /// Consistent `M` (`Matrix % MassValues`). ARPACK mode 3,
    /// `OP = (K - sigma M)^-1 M`.
    Consistent(&'a CsrMatrix<f64>),
}

/// Eigenpairs in `Eigen System Sorting` order.
#[derive(Clone, Debug)]
pub struct EigenResult {
    pub values: Vec<Complex64>,
    /// `vectors[i]` belongs to `values[i]`; length `K.nrows`.
    pub vectors: Vec<Vec<Complex64>>,
    /// ARPACK `IPARAM(3)`: "Number of update iterations taken".
    pub iterations: usize,
    /// `OP` applications; Elmer's "Number of eigensystem iterations" is this plus one.
    pub op_applications: usize,
    /// ARPACK `IPARAM(5)`: "Number of converged Ritz values".
    pub converged: usize,
    /// `CheckResiduals`: `||(K - Re(lambda) M) Re(x)||_2`, when requested.
    pub residuals: Option<Vec<f64>>,
}

struct ShiftInvert<'a> {
    n: usize,
    mass: Mass<'a>,
    solver: &'a dyn DirectSolver<f64>,
}

impl EigenOperator for ShiftInvert<'_> {
    fn dim(&self) -> usize {
        self.n
    }
    fn generalized(&self) -> bool {
        !matches!(self.mass, Mass::Identity)
    }
    fn apply_b(&self, x: &[f64], y: &mut [f64]) {
        apply_mass(self.mass, x, y);
    }
    fn apply_op(&mut self, _x: &[f64], bx: &[f64], y: &mut [f64]) -> Result<(), LinalgError> {
        y.copy_from_slice(bx);
        self.solver.solve_in_place(y)
    }
}

struct LumpedRegular<'a> {
    k: &'a CsrMatrix<f64>,
    diag: &'a [f64],
}

impl EigenOperator for LumpedRegular<'_> {
    fn dim(&self) -> usize {
        self.k.nrows
    }
    fn generalized(&self) -> bool {
        true
    }
    fn apply_b(&self, x: &[f64], y: &mut [f64]) {
        for ((yi, xi), d) in y.iter_mut().zip(x).zip(self.diag) {
            *yi = xi * d;
        }
    }
    fn apply_op(&mut self, x: &[f64], _bx: &[f64], y: &mut [f64]) -> Result<(), LinalgError> {
        self.k.apply(x, y);
        for (yi, d) in y.iter_mut().zip(self.diag) {
            *yi /= d;
        }
        Ok(())
    }
}

pub(crate) fn apply_mass(mass: Mass<'_>, x: &[f64], y: &mut [f64]) {
    match mass {
        Mass::Identity => y.copy_from_slice(x),
        Mass::Lumped(d) => {
            for ((yi, xi), di) in y.iter_mut().zip(x).zip(d) {
                *yi = xi * di;
            }
        }
        Mass::Consistent(m) => m.apply(x, y),
    }
}

/// `K - sigma M` on the union of the two patterns (Elmer's
/// `Values - SigmaR * MassValues` when they share one).
pub fn shifted_matrix(
    k: &CsrMatrix<f64>,
    mass: Mass<'_>,
    sigma: f64,
) -> Result<CsrMatrix<f64>, LinalgError> {
    let n = k.nrows;
    let mut out = CsrMatrix {
        nrows: n,
        ncols: k.ncols,
        row_ptr: Vec::with_capacity(n + 1),
        col_idx: Vec::with_capacity(k.nnz()),
        values: Vec::with_capacity(k.nnz()),
        diag: None,
    };
    out.row_ptr.push(0);
    let mut diag = Vec::with_capacity(n);
    for i in 0..n {
        let krow = k.row_ptr[i]..k.row_ptr[i + 1];
        let mut entries: Vec<(u32, f64)> = k.col_idx[krow.clone()]
            .iter()
            .copied()
            .zip(k.values[krow].iter().copied())
            .collect();
        match mass {
            Mass::Identity => entries.push((i as u32, -sigma)),
            Mass::Lumped(d) => entries.push((i as u32, -sigma * d[i])),
            Mass::Consistent(m) => {
                let r = m.row_ptr[i]..m.row_ptr[i + 1];
                for (&c, &v) in m.col_idx[r.clone()].iter().zip(&m.values[r]) {
                    entries.push((c, -sigma * v));
                }
            }
        }
        entries.sort_by_key(|e| e.0);
        let mut d = None;
        for (c, v) in entries {
            if out.col_idx.len() > out.row_ptr[i] && *out.col_idx.last().unwrap() == c {
                *out.values.last_mut().unwrap() += v;
            } else {
                if c as usize == i {
                    d = Some(out.col_idx.len());
                }
                out.col_idx.push(c);
                out.values.push(v);
            }
        }
        diag.push(
            d.ok_or_else(|| LinalgError::Dim(format!("row {i} of K - sigma M has no diagonal")))?,
        );
        out.row_ptr.push(out.col_idx.len());
    }
    out.diag = Some(diag);
    Ok(out)
}

fn check_dims(k: &CsrMatrix<f64>, mass: Mass<'_>) -> Result<(), LinalgError> {
    if k.nrows != k.ncols {
        return Err(LinalgError::Dim(format!("K is {}x{}", k.nrows, k.ncols)));
    }
    let ok = match mass {
        Mass::Identity => true,
        Mass::Lumped(d) => d.len() == k.nrows,
        Mass::Consistent(m) => m.nrows == k.nrows && m.ncols == k.ncols,
    };
    if ok {
        Ok(())
    } else {
        Err(LinalgError::Dim("mass matrix does not match K".into()))
    }
}

/// Solve `K x = lambda M x` for `cfg.nev` eigenpairs as `ArpackEigenSolve`
/// does, then sort (`EigenSystemSorting`), scale (`ScaleEigenVectors`, when
/// `cfg.normalize` asks) and check residuals (`CheckResiduals`).
///
/// For a consistent or identity mass, `solver` factors `K - sigma M` once
/// (Elmer's `Linear System Refactorize = False` after the first solve). An
/// iterative inner solve plugs in through any [`DirectSolver`] impl. The
/// matrices are expected as Elmer leaves them for the eigensolver: Dirichlet
/// rows of `K` reduced to the diagonal and those rows of `M` zero; any
/// `ScaleLinearSystem` is the caller's, as in `SolveCore`.
pub fn eigs(
    k: &CsrMatrix<f64>,
    mass: Mass<'_>,
    solver: &mut dyn DirectSolver<f64>,
    cfg: &EigenConfig,
) -> Result<EigenResult, LinalgError> {
    check_dims(k, mass)?;
    if cfg.ncv <= cfg.nev {
        return Err(LinalgError::Dim(
            "Number of Lanczos vectors must exceed the number of eigenvalues.".into(),
        ));
    }
    let lumped = matches!(mass, Mass::Lumped(_));
    let mut params = KrylovParams {
        nev: cfg.nev,
        ncv: cfg.ncv,
        which: ArpackWhich::for_mode(cfg.select, !lumped),
        tol: cfg.tol,
        max_iter: cfg.max_iter,
        symmetric: cfg.symmetric,
        purify: !lumped,
        seed: 1357,
    };
    let (ritz, sigma) = match mass {
        Mass::Lumped(d) => {
            if let Some(i) = d.iter().position(|&v| v == 0.0) {
                return Err(LinalgError::Singular(i));
            }
            params.purify = false;
            (
                krylov_schur(&mut LumpedRegular { k, diag: d }, &params)?,
                None,
            )
        }
        _ => {
            let sigma = cfg.shift;
            if sigma != 0.0 {
                solver.factor(&shifted_matrix(k, mass, sigma)?)?;
            } else {
                solver.factor(k)?;
            }
            let mut op = ShiftInvert {
                n: k.nrows,
                mass,
                solver: &*solver,
            };
            (krylov_schur(&mut op, &params)?, Some(sigma))
        }
    };

    // dseupd/dneupd: back-transform the Ritz values of OP.
    let mut values: Vec<Complex64> = ritz
        .theta
        .iter()
        .map(|&t| match sigma {
            Some(s) => Complex64::new(s, 0.0) + Complex64::new(1.0, 0.0) / t,
            None => t,
        })
        .collect();
    let mut vectors = ritz.vectors;
    if cfg.symmetric {
        // dseupd returns the pairs in ascending algebraic order.
        let mut idx: Vec<usize> = (0..values.len()).collect();
        idx.sort_by(|&a, &b| values[a].re.total_cmp(&values[b].re));
        values = idx.iter().map(|&i| values[i]).collect();
        vectors = idx
            .iter()
            .map(|&i| std::mem::take(&mut vectors[i]))
            .collect();
    }

    let perm = eigen_system_sorting(cfg.sorting, &values);
    let values: Vec<Complex64> = perm.iter().map(|&i| values[i]).collect();
    let mut vectors: Vec<Vec<Complex64>> = perm
        .iter()
        .map(|&i| std::mem::take(&mut vectors[i]))
        .collect();

    let residuals = cfg
        .compute_residuals
        .then(|| check_residuals(k, mass, &values, &vectors));
    if cfg.normalize != Normalize::None {
        scale_eigen_vectors(mass, &mut vectors, cfg.normalize == Normalize::Unity);
    }
    Ok(EigenResult {
        values,
        vectors,
        iterations: ritz.iterations,
        op_applications: ritz.op_applications,
        converged: ritz.converged,
        residuals,
    })
}
