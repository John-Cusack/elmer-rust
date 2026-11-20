//! `Eigen Analysis` of one solver (SolveCore.F90:1355-1495, `SolveEigenSystem`
//! and EigenSolve.F90:140 `ArpackEigenSolve`) with `elmer-eigen` as the
//! eigensolver.
//!
//! The solver assembles its stiffness into the system matrix and its mass
//! with [`crate::SolverCtx::default_update_mass`] (`Matrix % MassValues`, on
//! the same pattern). `DefaultDirichletBCs` zeroes the Dirichlet rows of the
//! mass, and with `Linear System Symmetric` the columns too
//! ([`dirichlet_mass`]: `CRS_ZeroRow`, `CRS_SetSymmDirichlet`). The solve
//! then:
//!
//! 1. scales `K` and `M` with the same diagonal (`ScaleLinearSystem(Solver,
//!    A)`, `Linear System Scaling` default on);
//! 2. solves `K x = lambda M x` in shift-invert mode, the inner solves with
//!    the section's linear solver (a direct one is factored once, as
//!    `Linear System Refactorize = False` after the first solve in
//!    `ArpackEigenSolve`);
//! 3. scales the eigenvectors back (`BackScaleLinearSystem(EigenScaling)`),
//!    sets the variable's norm from its (unchanged) values (`ComputeNorm`)
//!    and normalises the modes when asked (`ScaleEigenVectors`, after the
//!    back scaling);
//! 4. returns without `ComputeChange`, as `SolveSystem` does for an eigen
//!    system.
//!
//! The modes are kept in [`crate::Variable::eigen_values`] and
//! [`crate::Variable::eigen_vectors`].
//!
//! The eigenpairs come from `elmer_eigen::arpack`, the port of ARPACK
//! answered as `ArpackEigenSolve` answers it ([`arpack_eigen_solve`]):
//! `DSAUPD`/`DSEUPD` for a symmetric system (`Linear System Symmetric`,
//! Elmer's `Matrix % Symmetric`), `DNAUPD`/`DNEUPD` otherwise, complex
//! pairs included. They equal the pristine run's when the inner linear
//! solves do.
//! A `Stability Analysis` goes through `ArpackStabEigenSolve`'s mode-2
//! iteration instead ([`arpack_stab_eigen_solve`]).
//! A `Harmonic Mode` eigen analysis goes through `ChangeToHarmonicSystem`'s
//! complex system `A = K + i D` and `ArpackEigenSolveComplex`
//! (`ZNAUPD`/`ZNEUPD`, [`arpack_eigen_solve_complex`]), the inner solves
//! by the complex `SolveLinearSystem` with Elmer's complex product; it
//! needs a solver that assembles its damping matrix (`DefaultUpdateDamp`).
//! A partitioned run hands the real, consistent-mass system with an
//! iterative inner solver to its [`crate::ParEigenBackend`]
//! (`ParallelArpackEigenSolve`); without one, parallel eigen systems are
//! refused.
//! Refused: `Linear System Complex` eigen systems, a shifted complex one,
//! lumped-mass eigen systems, block and multigrid inner solvers, and
//! scaling methods other than the diagonal one.

use num_complex::Complex64;
use elmer_core::{CsrMatrix, LinalgError, ParamMap, Params};
use elmer_eigen::{EigenConfig, Mass, Normalize};

use crate::backend::LinearSolverChoice;
use crate::error::SolverError;
use crate::solver::{solve_linear, BlockEnv, SolverCtx};

const CALLER: &str = "EigenSolve";

/// `CRS_ZeroRow` / `CRS_SetSymmDirichlet` on `Matrix % MassValues`: the
/// mass rows of the Dirichlet dofs are zero, and with `Linear System
/// Symmetric` their columns as well.
pub(crate) fn dirichlet_mass(m: &mut CsrMatrix<f64>, fixed: &[(u32, f64)], params: &ParamMap) -> Result<(), SolverError> {
    if fixed.is_empty() {
        return Ok(());
    }
    let symmetric = crate::assembly::dirichlet_mode(params)?.symmetric;
    let mut on = vec![false; m.nrows];
    for &(r, _) in fixed {
        if let Some(o) = on.get_mut(r as usize) {
            *o = true;
        }
    }
    for r in 0..m.nrows {
        let row = m.row_ptr[r]..m.row_ptr[r + 1];
        if on[r] {
            m.values[row].fill(0.0);
        } else if symmetric {
            for k in row {
                if on[m.col_idx[k] as usize] {
                    m.values[k] = 0.0;
                }
            }
        }
    }
    Ok(())
}

/// Whether this visit's solve is an eigen solve (`EigenAnalysis` of
/// `SolveSystem`: `NOFEigenValues > 0` and `Eigen Analysis` on now).
pub(crate) fn active(ctx: &SolverCtx<'_>) -> bool {
    let p = ctx.params();
    ctx.nof_eigen() > 0 && p.logical("eigen analysis") == Some(true) && p.logical("skip eigen analysis") != Some(true)
}

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("Eigen Analysis: {}", what.into()))
}

/// An inner solve of `ArpackEigenSolve`: `a x = b` from the guess in `x`.
type InnerSolve<'a> = dyn FnMut(&CsrMatrix<f64>, &[f64], &mut [f64]) -> Result<(), SolverError> + 'a;

/// `ArpackEigenSolve` (EigenSolve.F90:140) of a symmetric (`Matrix %
/// Symmetric`), consistent-mass system with the ported ARPACK (mode 3,
/// shift-invert): `K - sigma M` in place of `K` for the iteration
/// (`Values - SigmaR * MassValues`, put back after); the reverse
/// communication as Elmer answers it; `DSEUPD`; the pairs in `Eigen System
/// Sorting` order. `inner(a, b, x)` solves `a x = b` from the guess in `x`.
#[allow(clippy::too_many_arguments)]
fn arpack_eigen_solve(
    arp: &mut elmer_eigen::arpack::saupd::Arpack,
    k: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    dofs: usize,
    cfg: &EigenConfig,
    symmetric: bool,
    inner: &mut InnerSolve<'_>,
) -> Result<elmer_eigen::EigenResult, SolverError> {
    use elmer_eigen::arpack::{naupd, neupd, saupd, seupd};
    use elmer_eigen::ArpackWhich as W;
    use num_complex::Complex64;
    let n = k.nrows;
    let (nev, ncv) = (cfg.nev, cfg.ncv);
    if ncv <= nev {
        return Err(SolverError::fatal(CALLER, "Number of Lanczos vectors must exceed the number of eigenvalues."));
    }
    let wsel = elmer_eigen::ArpackWhich::for_mode(cfg.select, true);
    let which = match wsel {
        W::LM => saupd::Which::LM,
        W::SM => saupd::Which::SM,
        // DSAUPD refuses the real/imaginary-part selections (info = -5).
        _ if symmetric => return Err(SolverError::fatal(CALLER, " Error with DNAUPD, info =           -5")),
        _ => saupd::Which::LM,
    };
    let nwhich = match wsel {
        W::LM => naupd::NWhich::LM,
        W::SM => naupd::NWhich::SM,
        W::LR => naupd::NWhich::LR,
        W::SR => naupd::NWhich::SR,
        W::LI => naupd::NWhich::LI,
        W::SI => naupd::NWhich::SI,
    };
    let sigma = cfg.shift;
    if sigma != 0.0 {
        for (kv, mv) in k.values.iter_mut().zip(&mass.values) {
            *kv -= sigma * mv;
        }
    }
    let mut iter = 1usize;
    let mut resid = vec![0.0; n];
    let mut v = vec![0.0; n * ncv];
    let mut workd = vec![0.0; 3 * n];
    let mut workl = vec![0.0; 3 * ncv * ncv + 6 * ncv];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = cfg.max_iter as i32;
    iparam[6] = 3;
    let mut tol = cfg.tol;
    let mut failure: Option<SolverError> = None;
    let au = {
        let kk: &CsrMatrix<f64> = k;
        let mm = crate::matvec::ElmerMatvec::new(mass, dofs);
        let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
            let r = match ido {
                -1 | 1 => {
                    iter += 1;
                    let b: Vec<f64> = if ido == -1 {
                        // y = M x, then x := M x, the right-hand side.
                        let mut y = vec![0.0; n];
                        elmer_core::LinearOperator::apply(&mm, &wd[ip[0]..ip[0] + n], &mut y);
                        wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                        wd[ip[0]..ip[0] + n].copy_from_slice(&y);
                        y
                    } else {
                        wd[ip[2]..ip[2] + n].to_vec()
                    };
                    inner(kk, &b, &mut wd[ip[1]..ip[1] + n])
                }
                2 => {
                    let mut y = vec![0.0; n];
                    elmer_core::LinearOperator::apply(&mm, &wd[ip[0]..ip[0] + n], &mut y);
                    wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                    Ok(())
                }
                _ => Err(SolverError::fatal(CALLER, format!("unexpected ARPACK request ido = {ido}"))),
            };
            r.map_err(|e| {
                let msg = e.to_string();
                failure = Some(e);
                msg
            })
        };
        if symmetric {
            saupd::dsaupd(&mut op, arp, true, n, which, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
                .map(|a| (a.iparam, a.info, Some(a), None))
        } else {
            naupd::dnaupd(&mut op, arp, true, n, nwhich, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
                .map(|a| (a.iparam, a.info, None, Some(a)))
        }
    };
    let (iparam_out, info, sau, nau) = match au {
        Ok(a) => a,
        Err(m) => return Err(failure.unwrap_or_else(|| SolverError::fatal(CALLER, m))),
    };
    let (values, vectors);
    if let Some(au) = sau {
        if info != 0 {
            return Err(SolverError::fatal(CALLER, format!(" Error with DNAUPD, info = {info:>11}")));
        }
        let (d, ierr) = seupd::dseupd(sigma, true, n, which, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl);
        if ierr != 0 {
            return Err(SolverError::fatal(CALLER, format!("  Error with DNEUPD, info = {ierr:>11}")));
        }
        let vals: Vec<Complex64> = (0..nev).map(|i| Complex64::new(d.get(i).copied().unwrap_or(0.0), 0.0)).collect();
        let perm = elmer_eigen::eigen_system_sorting(cfg.sorting, &vals);
        values = perm.iter().map(|&p| vals[p]).collect::<Vec<_>>();
        vectors = perm.iter().map(|&p| v[p * n..p * n + n].iter().map(|&x| Complex64::new(x, 0.0)).collect()).collect::<Vec<Vec<_>>>();
    } else {
        let au = nau.expect("the non-symmetric path ran DNAUPD");
        // `IF (kinfo < 0)`: DNAUPD's error, before DNEUPD.
        if info < 0 {
            return Err(SolverError::fatal(CALLER, format!(" Error with DNAUPD, info = {info:>11}")));
        }
        let r = neupd::dneupd(sigma, 0.0, true, n, nwhich, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl);
        if r.info != 0 {
            return Err(SolverError::fatal(CALLER, format!(" Error with DNEUPD, info = {:>11}", r.info)));
        }
        if info == 1 {
            return Err(SolverError::fatal(CALLER, "Maximum number of iterations reached."));
        } else if info == 3 {
            return Err(SolverError::fatal(CALLER, "No shifts could be applied during implicit Arnoldi update, try increasing NCV."));
        }
        // `EigValues(i) = CMPLX(D(i,1), D(i,2))`, sorted by
        // EigenSystemSorting. The vector of value p is taken from column
        // k = 1 + sum over j < p of (1 if D(j,2) == 0 else 2) of V, a
        // complex one as V(:,k) + i V(:,k+1): ArpackEigenSolve's count,
        // in which both members of a conjugate pair add 2.
        let vals: Vec<Complex64> = (0..nev).map(|i| Complex64::new(r.dr[i], r.di[i])).collect();
        let perm = elmer_eigen::eigen_system_sorting(cfg.sorting, &vals);
        values = perm.iter().map(|&p| vals[p]).collect::<Vec<_>>();
        let col = |c: usize| -> Result<&[f64], SolverError> {
            v.get(c * n..c * n + n).ok_or_else(|| SolverError::fatal(CALLER, "eigenvector column outside the Arnoldi basis"))
        };
        let mut vecs = Vec::with_capacity(nev);
        for &p in &perm {
            let k0: usize = (0..p).map(|j| if r.di[j] == 0.0 { 1 } else { 2 }).sum();
            let re = col(k0)?;
            let x: Vec<Complex64> = if r.di[p] != 0.0 {
                let im = col(k0 + 1)?;
                re.iter().zip(im).map(|(&a, &b)| Complex64::new(a, b)).collect()
            } else {
                re.iter().map(|&a| Complex64::new(a, 0.0)).collect()
            };
            vecs.push(x);
        }
        vectors = vecs;
    }
    if sigma != 0.0 {
        for (kv, mv) in k.values.iter_mut().zip(&mass.values) {
            *kv += sigma * mv;
        }
    }
    let residuals = cfg
        .compute_residuals
        .then(|| elmer_eigen::check_residuals(k, Mass::Consistent(mass), &values, &vectors));
    Ok(elmer_eigen::EigenResult {
        values,
        vectors,
        iterations: iparam_out[2].max(0) as usize,
        op_applications: iter - 1,
        converged: iparam_out[4].max(0) as usize,
        residuals,
    })
}

/// `ArpackStabEigenSolve` (EigenSolve.F90:803), the eigen solve of a
/// `Stability Analysis`: `DSAUPD` in mode 2 with `B = K` and `OP = K^-1 M`,
/// `M` being what the solver put in the mass slot (minus the geometric
/// stiffness), no shift, `NCV = 3 NEIG + 1`. `ido = -1, 1`: `y = M x` copied
/// over `x`, then `K y = x` from the guess `M x`; `ido = 2`: `y = K x`. The
/// eigenvalues are `1/theta` in `Eigen System Sorting` order, and the
/// vector of a negative one is multiplied by `i`.
fn arpack_stab_eigen_solve(
    arp: &mut elmer_eigen::arpack::saupd::Arpack,
    k: &CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    dofs: usize,
    cfg: &EigenConfig,
    inner: &mut InnerSolve<'_>,
) -> Result<elmer_eigen::EigenResult, SolverError> {
    use elmer_eigen::arpack::{saupd, seupd};
    use elmer_eigen::ArpackWhich as W;
    use num_complex::Complex64;
    const STAB: &str = "StabEigenSolve";
    let n = k.nrows;
    let nev = cfg.nev;
    let ncv = 3 * nev + 1;
    if ncv > n {
        return Err(SolverError::fatal("StabEigenSolve", format!("Error with DSAUPD, info = {:>11}", -3)));
    }
    let which = match elmer_eigen::ArpackWhich::for_mode(cfg.select, true) {
        W::LM => saupd::Which::LM,
        W::SM => saupd::Which::SM,
        _ => return Err(SolverError::fatal(STAB, format!("Error with DSAUPD, info = {:>11}", -5))),
    };
    let mut iter = 1usize;
    let mut resid = vec![0.0; n];
    let mut v = vec![0.0; n * ncv];
    let mut workd = vec![0.0; 3 * n];
    let mut workl = vec![0.0; ncv * (ncv + 8)];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = cfg.max_iter as i32;
    iparam[6] = 2;
    let mut tol = cfg.tol;
    let mut failure: Option<SolverError> = None;
    let au = {
        let mm = crate::matvec::ElmerMatvec::new(mass, dofs);
        let kk = crate::matvec::ElmerMatvec::new(k, dofs);
        let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
            let r = match ido {
                -1 | 1 => {
                    iter += 1;
                    let mut y = vec![0.0; n];
                    elmer_core::LinearOperator::apply(&mm, &wd[ip[0]..ip[0] + n], &mut y);
                    wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                    wd[ip[0]..ip[0] + n].copy_from_slice(&y);
                    inner(k, &y, &mut wd[ip[1]..ip[1] + n])
                }
                2 => {
                    let mut y = vec![0.0; n];
                    elmer_core::LinearOperator::apply(&kk, &wd[ip[0]..ip[0] + n], &mut y);
                    wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                    Ok(())
                }
                _ => Err(SolverError::fatal(STAB, format!("unexpected ARPACK request ido = {ido}"))),
            };
            r.map_err(|e| {
                let msg = e.to_string();
                failure = Some(e);
                msg
            })
        };
        saupd::dsaupd(&mut op, arp, true, n, which, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
    };
    let au = match au {
        Ok(a) => a,
        Err(m) => return Err(failure.unwrap_or_else(|| SolverError::fatal(STAB, m))),
    };
    if au.info != 0 {
        return Err(SolverError::fatal(STAB, format!("Error with DSAUPD, info = {:>11}", au.info)));
    }
    let (d, ierr) = seupd::dseupd(0.0, true, n, which, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl);
    if ierr != 0 {
        return Err(SolverError::fatal(STAB, format!(" Error with DSEUPD, info = {ierr:>11}")));
    }
    // `EigValues(i) = CMPLX(1.0d0 / D(i,1), D(i,2))`; D(:,2) is zero here.
    let vals: Vec<Complex64> = (0..nev).map(|i| Complex64::new(1.0 / d.get(i).copied().unwrap_or(0.0), 0.0)).collect();
    let perm = elmer_eigen::eigen_system_sorting(cfg.sorting, &vals);
    let values: Vec<Complex64> = perm.iter().map(|&p| vals[p]).collect();
    let vectors: Vec<Vec<Complex64>> = perm
        .iter()
        .zip(&values)
        .map(|(&p, lam)| {
            v[p * n..p * n + n]
                .iter()
                .map(|&x| if lam.re < 0.0 { Complex64::new(x, 0.0) * Complex64::new(0.0, 1.0) } else { Complex64::new(x, 0.0) })
                .collect()
        })
        .collect();
    Ok(elmer_eigen::EigenResult {
        values,
        vectors,
        iterations: au.iparam[2].max(0) as usize,
        op_applications: iter - 1,
        converged: au.iparam[4].max(0) as usize,
        residuals: None,
    })
}

/// An inner solve `A x = b` of the complex eigen iteration, `x` holding the guess.
type ComplexInner<'a> = dyn FnMut(&[Complex64], &mut [Complex64]) -> Result<(), SolverError> + 'a;

/// `ArpackEigenSolveComplex` (EigenSolve.F90:1150): `ZNAUPD` in mode 3
/// with `BMAT = 'G'` on `A x = lambda M x` (complex `A`, `M`), `OP = A^-1
/// M` (no shift), `inner` solving `A x = b`; `ZNEUPD` for the Ritz
/// pairs, then `EigenSystemSorting`.
fn arpack_eigen_solve_complex(
    zarp: &mut elmer_eigen::arpack::znaupd::ZArpack,
    a: &CsrMatrix<Complex64>,
    mass: &CsrMatrix<Complex64>,
    cfg: &EigenConfig,
    sigma: Complex64,
    inner: &mut ComplexInner<'_>,
) -> Result<elmer_eigen::EigenResult, SolverError> {
    use elmer_eigen::arpack::naupd::NWhich;
    use elmer_eigen::arpack::{znaupd, zneupd};
    use elmer_eigen::ArpackWhich as W;
    let n = a.nrows;
    let (nev, ncv) = (cfg.nev, cfg.ncv);
    if ncv <= nev {
        return Err(SolverError::fatal(CALLER, "Number of Lanczos vectors must exceed the number of eigenvalues."));
    }
    let which = match elmer_eigen::ArpackWhich::for_mode(cfg.select, true) {
        W::LM => NWhich::LM,
        W::SM => NWhich::SM,
        W::LR => NWhich::LR,
        W::SR => NWhich::SR,
        W::LI => NWhich::LI,
        W::SI => NWhich::SI,
    };
    let z = Complex64::new(0.0, 0.0);
    let mut iter = 1usize;
    let mut resid = vec![z; n];
    let mut v = vec![z; n * ncv];
    let mut workd = vec![z; 3 * n];
    let mut workl = vec![z; 3 * ncv * ncv + 6 * ncv];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = cfg.max_iter as i32;
    iparam[6] = 3;
    let mut tol = cfg.tol;
    let mut failure: Option<SolverError> = None;
    let au = {
        let mut op = |ido: i32, ip: [usize; 3], wd: &mut [Complex64]| -> Result<(), String> {
            let r = match ido {
                -1 | 1 => {
                    iter += 1;
                    let b: Vec<Complex64> = if ido == -1 {
                        let mut y = vec![z; n];
                        crate::matvec::crs_complex_matrix_vector_multiply(mass, &wd[ip[0]..ip[0] + n], &mut y);
                        wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                        y
                    } else {
                        wd[ip[2]..ip[2] + n].to_vec()
                    };
                    // The guess is what workd(ipntr(2)) holds.
                    inner(&b, &mut wd[ip[1]..ip[1] + n])
                }
                2 => {
                    let mut y = vec![z; n];
                    crate::matvec::crs_complex_matrix_vector_multiply(mass, &wd[ip[0]..ip[0] + n], &mut y);
                    wd[ip[1]..ip[1] + n].copy_from_slice(&y);
                    Ok(())
                }
                _ => Err(SolverError::fatal(CALLER, format!("unexpected ARPACK request ido = {ido}"))),
            };
            r.map_err(|e| {
                let msg = e.to_string();
                failure = Some(e);
                msg
            })
        };
        znaupd::znaupd(&mut op, zarp, true, n, which, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
    };
    let au = match au {
        Ok(a) => a,
        Err(msg) => return Err(failure.unwrap_or_else(|| SolverError::fatal(CALLER, msg))),
    };
    match au.info {
        0 | 2 => {}
        1 => return Err(SolverError::fatal(CALLER, "Maximum number of iterations reached.")),
        3 => return Err(SolverError::fatal(CALLER, "No shifts could be applied during implicit Arnoldi update, try increasing NCV.")),
        i => return Err(SolverError::fatal(CALLER, format!(" Error with ZNAUPD, info = {i:>11}"))),
    }
    let r = zneupd::zneupd(true, sigma, true, n, which, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl);
    if r.info != 0 {
        return Err(SolverError::fatal(CALLER, format!(" Error with ZNEUPD, info = {:>11}", r.info)));
    }
    let vals: Vec<Complex64> = (0..nev).map(|i| r.d.get(i).copied().unwrap_or(z)).collect();
    let perm = elmer_eigen::eigen_system_sorting(cfg.sorting, &vals);
    let values: Vec<Complex64> = perm.iter().map(|&p| vals[p]).collect();
    let vectors: Vec<Vec<Complex64>> = perm.iter().map(|&p| v[p * n..p * n + n].to_vec()).collect();
    Ok(elmer_eigen::EigenResult {
        values,
        vectors,
        iterations: au.iparam[2].max(0) as usize,
        op_applications: iter - 1,
        converged: au.iparam[4].max(0) as usize,
        residuals: None,
    })
}

fn sim_linear_zarpack<'a>(
    linear: &'a mut crate::LinearSolverRegistry,
    zarp: &'a mut elmer_eigen::arpack::znaupd::ZArpack,
) -> (&'a mut crate::LinearSolverRegistry, &'a mut elmer_eigen::arpack::znaupd::ZArpack) {
    (linear, zarp)
}

/// `ArpackEigenSolveComplex` of a `Linear System Complex` system (the
/// interleaved real `Matrix % Values`, `MassValues`), in place as the
/// Fortran works on the solver's matrix:
/// `Values = Values - Sigma * MassValues` (the real part:
/// `V - Re(sigma) M`) for `Eigen System Shift` (+ `Im`);
/// `ScaleLinearSystem`'s complex diagonal `1/SQRT(ABS(CMPLX(A_ii, -A_i,i+1)))`
/// on `Values` and `MassValues` (`Values(j) * (Diag(i) * Diag(col))`);
/// ZNAUPD (mode 3, `BMAT = 'G'`) with the inner `DirectSolver` /
/// `IterSolver` on the interleaved system (`Linear System Nullify Guess`
/// starting from zero) and `CRS_ComplexMatrixVectorMultiply` with the
/// mass; ZNEUPD with `sigma`; `EigenSystemSorting`; the vectors times
/// `Diag(1:n:2)`, the values divided back (`/ (Diag(i) * Diag(col))`) and
/// the shift added back (`V + Re(sigma) M`).
fn complex_system_eigen(
    p: &ParamMap,
    k: &mut CsrMatrix<f64>,
    mass: &mut CsrMatrix<f64>,
    (linear, zarp): (&mut crate::LinearSolverRegistry, &mut elmer_eigen::arpack::znaupd::ZArpack),
    cfg: &EigenConfig,
    log: &dyn elmer_core::Logger,
) -> Result<elmer_eigen::EigenResult, SolverError> {
    let n = k.nrows;
    if n % 2 != 0 || mass.values.len() != k.values.len() {
        return Err(SolverError::fatal(CALLER, "a complex eigen system needs an interleaved matrix and mass on its pattern"));
    }
    let sigma = Complex64::new(cfg.shift, p.real("eigen system shift im").unwrap_or(0.0));
    if sigma != Complex64::new(0.0, 0.0) {
        for (v, m) in k.values.iter_mut().zip(&mass.values) {
            *v -= sigma.re * m;
        }
    }
    let scale = p.logical("linear system skip scaling") != Some(true) && p.logical("linear system scaling").unwrap_or(true);
    let mut diag = vec![1.0f64; n];
    if scale {
        for i in (0..n).step_by(2) {
            let j = (k.row_ptr[i]..k.row_ptr[i + 1]).find(|&j| k.col_idx[j] as usize == i);
            let (dr, di) = match j {
                Some(j) => (k.values[j], k.values[j + 1]),
                None => (0.0, 0.0),
            };
            let s = Complex64::new(dr, -di).norm().sqrt();
            let d = if s > f64::MIN_POSITIVE { 1.0 / s } else { 1.0 };
            diag[i] = d;
            diag[i + 1] = d;
        }
        for mat in [&mut *k, &mut *mass] {
            for i in 0..n {
                for j in mat.row_ptr[i]..mat.row_ptr[i + 1] {
                    mat.values[j] *= diag[i] * diag[mat.col_idx[j] as usize];
                }
            }
        }
    }
    let ac = crate::linsys::to_complex(k)?;
    let mc = crate::linsys::to_complex(mass)?;
    let kind = p.string("linear system solver").map(elmer_core::normalize_key).unwrap_or_default();
    let mut direct = match kind.as_str() {
        "direct" => match linear.select_for(p, k, log)? {
            LinearSolverChoice::Direct(d) => Some(d),
            LinearSolverChoice::Iterative { .. } => None,
        },
        _ => None,
    };
    if direct.is_none() {
        return Err(unsupported("an iterative inner solve of a complex eigen system"));
    }
    let mut factored = false;
    let kk: &CsrMatrix<f64> = k;
    let mut inner = |b: &[Complex64], x: &mut [Complex64]| -> Result<(), SolverError> {
        let d = direct.as_mut().expect("direct solver");
        if !factored {
            d.factor(kk)?;
            factored = true;
        }
        let mut y: Vec<f64> = b.iter().flat_map(|z| [z.re, z.im]).collect();
        d.solve_in_place(&mut y)?;
        for (i, xi) in x.iter_mut().enumerate() {
            *xi = Complex64::new(y[2 * i], y[2 * i + 1]);
        }
        Ok(())
    };
    let r = arpack_eigen_solve_complex(zarp, &ac, &mc, cfg, sigma, &mut inner);
    let mut result = r?;
    for v in &mut result.vectors {
        for (i, x) in v.iter_mut().enumerate() {
            *x *= diag[2 * i];
        }
    }
    if scale {
        for mat in [&mut *k, &mut *mass] {
            for i in 0..n {
                for j in mat.row_ptr[i]..mat.row_ptr[i + 1] {
                    mat.values[j] /= diag[i] * diag[mat.col_idx[j] as usize];
                }
            }
        }
    }
    if sigma != Complex64::new(0.0, 0.0) {
        for (v, m) in k.values.iter_mut().zip(&mass.values) {
            *v += sigma.re * m;
        }
    }
    Ok(result)
}

/// The complex system `ChangeToHarmonicSystem` makes of a second-order
/// eigen system (SolveCore.F90:4811): `A = K + i D` with `D` the damping,
/// `M` the mass, on `K`'s pattern, as `CRS_ComplexMatrixVectorMultiply`
/// reads the interleaved child matrix (`Values(kr+1) = -D`, then
/// `CMPLX(Values(j), -Values(j+1))`); the mass child's off-diagonal stays
/// zero, so `M` gets the imaginary part `-0`.
fn harmonic_system(k: &CsrMatrix<f64>, damp: &CsrMatrix<f64>, mass: &CsrMatrix<f64>) -> (CsrMatrix<Complex64>, CsrMatrix<Complex64>) {
    let mk = |vals: Vec<Complex64>| CsrMatrix {
        nrows: k.nrows,
        ncols: k.ncols,
        row_ptr: k.row_ptr.clone(),
        col_idx: k.col_idx.clone(),
        values: vals,
        diag: k.diag.clone(),
    };
    let a: Vec<Complex64> = (0..k.values.len())
        .map(|j| {
            Complex64::new(k.values[j], -(-damp.values[j]))
        })
        .collect();
    let m: Vec<Complex64> = mass.values.iter().map(|&mv| Complex64::new(mv, -0.0)).collect();
    (mk(a), mk(m))
}

/// The eigen solve of `SolveSystem`; returns the variable's norm.
pub(crate) fn solve(ctx: &mut SolverCtx<'_>) -> Result<f64, SolverError> {
    let (mesh, log, comm, ix) = (ctx.mesh, ctx.log, ctx.comm, ctx.slot);
    let p = ctx.params().clone();
    // A partitioned run solves through its backend (ParallelArpackEigenSolve)
    // the real, consistent-mass system with an iterative inner solver.
    let par_eigen = if comm.size() > 1 {
        let Some(b) = ctx.sim.backends.eigen.clone() else {
            return Err(unsupported("parallel eigen systems (ParallelArpackEigenSolve)"));
        };
        let refused = [
            ("linear system complex", "a complex system"),
            ("harmonic mode", "Harmonic Mode"),
            ("stability analysis", "Stability Analysis"),
            ("eigen system compute residuals", "Eigen System Compute Residuals"),
            ("eigen system normalize to unity", "eigenvector normalization"),
            ("eigen system mass normalize", "eigenvector normalization"),
            ("eigen system damped", "Eigen System Damped"),
        ];
        if let Some((_, what)) = refused.iter().find(|(k, _)| p.logical(k) == Some(true)) {
            return Err(unsupported(format!("{what} in parallel")));
        }
        let kind = p.string("linear system solver").map(elmer_core::normalize_key).unwrap_or_default();
        if kind != "iterative" {
            return Err(unsupported(format!("Linear System Solver = {kind} in parallel")));
        }
        Some(b)
    } else {
        None
    };
    // The inner direct solves use the section's method, UMFPACK by default.
    let p = crate::linsys::serial_direct_default(&p, comm).unwrap_or(p);
    let complex_system = p.logical("linear system complex") == Some(true) && p.logical("linear system skip complex") != Some(true);
    if complex_system && p.logical("harmonic mode") == Some(true) {
        return Err(unsupported("a Harmonic Mode of a complex system"));
    }
    if p.logical("lumped mass matrix") == Some(true) {
        return Err(unsupported("a lumped mass matrix"));
    }
    if p.logical("linear system row equilibration") == Some(true)
        || p.string("linear system scaling method").is_some_and(|m| elmer_core::normalize_key(m) != "diagonal")
    {
        return Err(unsupported("scaling other than the diagonal one"));
    }
    let kind = p.string("linear system solver").map(elmer_core::normalize_key).unwrap_or_default();
    if kind != "direct" && kind != "iterative" {
        return Err(unsupported(format!("Linear System Solver = {kind}")));
    }
    let symmetric = p.logical("linear system symmetric").unwrap_or(false); // Matrix % Symmetric
    let mut cfg = EigenConfig::from_params(&p, symmetric).map_err(|e| match e {
        LinalgError::Unsupported(m) => unsupported(m),
        e => SolverError::from(e),
    })?;
    // ScaleEigenVectors runs after the back scaling, on the original mass.
    let normalize = std::mem::replace(&mut cfg.normalize, Normalize::None);

    let sim = &mut *ctx.sim;
    let slot = &mut sim.slots[ix];
    let id = slot.variable.ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let sys = slot.system.as_mut().ok_or_else(|| SolverError::fatal(CALLER, "no assembled system"))?;
    let mut mass = sys.mass.take().ok_or_else(|| SolverError::fatal(CALLER, "no mass matrix (DefaultUpdateMass)"))?;
    let k = &mut sys.matrix;
    let n = k.nrows;
    if complex_system {
        let r = complex_system_eigen(&p, k, &mut mass, sim_linear_zarpack(&mut sim.backends.linear, &mut sim.zarpack), &cfg, log);
        sys.mass = Some(mass);
        let result = r?;
        elmer_eigen::log_summary(&result, cfg.tol, log);
        if normalize != Normalize::None {
            return Err(unsupported("eigenvector normalization of a complex system"));
        }
        let var = sim.vars.by_id_mut(id);
        let spec = crate::norm::NormSpec::from_params(&p, var.dofs);
        let (vals, m) = crate::norm::select_norm_entries(&spec, &var.values, n.min(var.values.len()), var.perm.as_deref(), mesh.n_nodes(), log);
        let norm = crate::norm::compute_norm(&spec, &vals, m, comm);
        var.conv.norm = norm;
        var.eigen_values = result.values;
        var.eigen_vectors = result.vectors;
        log.info("SolveLinearSystem", &format!("Repointing {} eigenvalue components for: {}", var.dofs, var.name), 3);
        return Ok(norm);
    }
    if p.logical("harmonic mode") == Some(true) {
        // ChangeToHarmonicSystem + ArpackEigenSolveComplex. The damping
        // matrix is what the solver gave DefaultUpdateDamp; a solver that
        // never called it may have dropped a damping term, so it is refused.
        let Some(damp) = sys.damp.as_ref() else {
            sys.mass = Some(mass);
            return Err(unsupported("Harmonic Mode for a solver that does not assemble its damping matrix (DefaultUpdateDamp)"));
        };
        let (mut a, mut mc) = harmonic_system(k, damp, &mass);
        sys.mass = Some(mass);
        // ScaleLinearSystem on the complex matrix: 1/sqrt|A_ii|.
        let scale = p.logical("linear system skip scaling") != Some(true) && p.logical("linear system scaling").unwrap_or(true);
        let mut dsc = vec![1.0f64; n];
        if scale {
            for i in 0..n {
                let mut di = Complex64::new(0.0, 0.0);
                for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                    if a.col_idx[j] as usize == i {
                        di = a.values[j];
                    }
                }
                let s = di.re.hypot(di.im).sqrt();
                dsc[i] = if s > f64::MIN_POSITIVE { 1.0 / s } else { 1.0 };
            }
            for mat in [&mut a, &mut mc] {
                for i in 0..n {
                    for j in mat.row_ptr[i]..mat.row_ptr[i + 1] {
                        let f = dsc[i] * dsc[mat.col_idx[j] as usize];
                        mat.values[j] = Complex64::new(mat.values[j].re * f, mat.values[j].im * f);
                    }
                }
            }
        }
        let mut q = p.clone();
        q.insert("Linear System Scaling", elmer_core::Value::Logical(false));
        let mut inner = |b: &[Complex64], x: &mut [Complex64]| -> Result<(), SolverError> {
            crate::linsys::solve_complex_elmer_product(&q, &a, b, x, log, comm)
        };
        if cfg.shift != 0.0 {
            return Err(unsupported("a shifted complex eigen system (Eigen System Shift)"));
        }
        let mut result = arpack_eigen_solve_complex(&mut sim.zarpack, &a, &mc, &cfg, Complex64::new(0.0, 0.0), &mut inner)?;
        // BackScaleLinearSystem(EigenScaling = .TRUE.).
        if scale {
            for v in &mut result.vectors {
                for (x, s) in v.iter_mut().zip(&dsc) {
                    *x *= *s;
                }
            }
        }
        elmer_eigen::log_summary(&result, cfg.tol, log);
        if normalize != Normalize::None {
            if let Some(mm) = sys.mass.as_ref() {
                elmer_eigen::scale_eigen_vectors(Mass::Consistent(mm), &mut result.vectors, normalize == Normalize::Unity);
            }
        }
        // ChangeToHarmonicSystem's imaginary part `<name> im` (and its
        // components), an output variable that keeps the zero values of the
        // harmonic system's (never solved) imaginary dofs.
        let (vname, vdofs, vperm, vlen) = {
            let v = sim.vars.by_id(id);
            (v.name.clone(), v.dofs, v.perm.clone(), v.values.len())
        };
        let im = format!("{vname} im");
        if sim.vars.get(&im).is_none() {
            sim.vars.add(crate::Variable::new(&im, elmer_core::FieldKind::Nodal, vdofs, vperm, vec![0.0; vlen]))?;
        }
        let var = sim.vars.by_id_mut(id);
        let spec = crate::norm::NormSpec::from_params(&p, var.dofs);
        let (vals, m) = crate::norm::select_norm_entries(&spec, &var.values, n.min(var.values.len()), var.perm.as_deref(), mesh.n_nodes(), log);
        let norm = crate::norm::compute_norm(&spec, &vals, m, comm);
        var.conv.norm = norm;
        var.eigen_values = result.values;
        var.eigen_vectors = result.vectors;
        log.info("SolveLinearSystem", &format!("Repointing {} eigenvalue components for: {} complex", 2 * var.dofs, var.name), 3);
        return Ok(norm);
    }

    // ScaleLinearSystem(Solver, A): K and MassValues with one diagonal.
    let scale = p.logical("linear system skip scaling") != Some(true) && p.logical("linear system scaling").unwrap_or(true);
    let scaling = scale.then(|| {
        let opts = elmer_sparse::ScaleOptions { method: elmer_sparse::ScalingMethod::Diagonal, apply: true, rhs_scaling: false };
        elmer_sparse::scale_system(k, None, None, &opts)
    });
    if let Some(rs) = &scaling {
        let mut v = std::mem::take(&mut mass.values);
        rs.scale_values(&mass, &mut v);
        mass.values = v;
    }

    let result: Result<elmer_eigen::EigenResult, SolverError> = {
        let dofs = sim.vars.by_id(id).dofs;
        let mut direct = match kind.as_str() {
            "direct" => match sim.backends.linear.select_for(&p, k, log)? {
                LinearSolverChoice::Direct(d) => Some(d),
                LinearSolverChoice::Iterative { .. } => None,
            },
            _ => None,
        };
        let mut q = p.clone();
        q.insert("Linear System Scaling", elmer_core::Value::Logical(false));
        // ArpackEigenSolve's inner solves: a direct solver factors the
        // shifted matrix once (`Linear System Refactorize = False` after
        // the first); IterSolver starts from what `workd(ipntr(2))` holds.
        let (linear, simulation) = (&sim.backends.linear, &sim.model.simulation);
        let mut factored = false;
        let mut inner = |a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]| -> Result<(), SolverError> {
            match direct.as_mut() {
                Some(d) => {
                    if !factored {
                        d.factor(a)?;
                        factored = true;
                    }
                    let mut y = b.to_vec();
                    d.solve_in_place(&mut y)?;
                    x.copy_from_slice(&y);
                    Ok(())
                }
                None => {
                    let mut a2 = a.clone();
                    let mut bb = b.to_vec();
                    let bnorm = bb.iter().map(|v| v * v).sum::<f64>().sqrt();
                    let env = BlockEnv { dofs, ndeg: dofs, simulation, schur: None, multigrid: None, alias: None, precond: None , prec_values: None};
                    solve_linear(linear, &env, &q, &mut a2, &mut bb, x, bnorm, false, log, comm).map(|_| ())
                }
            }
        };
        if let Some(b) = &par_eigen {
            // The partitioned run's slots have `Linear System Scaling`
            // off here (the backend scales with the summed diagonal).
            b.solve(&p, k, &mut mass, &cfg, symmetric, &mut sim.arpack, log)
        } else if p.logical("stability analysis") == Some(true) {
            arpack_stab_eigen_solve(&mut sim.arpack, k, &mass, dofs, &cfg, &mut inner)
        } else {
            arpack_eigen_solve(&mut sim.arpack, k, &mass, dofs, &cfg, symmetric, &mut inner)
        }
    };

    // BackScaleLinearSystem(EigenScaling = .TRUE.), then K and M unscaled.
    let mut result = match result {
        Ok(r) => r,
        Err(e) => {
            if let Some(rs) = &scaling {
                rs.unscale(k, None, None);
                let mut v = std::mem::take(&mut mass.values);
                rs.unscale_values(&mass, &mut v);
                mass.values = v;
            }
            sys.mass = Some(mass);
            return Err(match e {
                e @ (SolverError::Fatal { .. } | SolverError::Unsupported(_)) => e,
                e => SolverError::fatal(CALLER, e.to_string()),
            });
        }
    };
    if let Some(rs) = &scaling {
        let d = rs.diag();
        for v in &mut result.vectors {
            for (x, s) in v.iter_mut().zip(d) {
                *x *= *s;
            }
        }
        rs.unscale(k, None, None);
        let mut v = std::mem::take(&mut mass.values);
        rs.unscale_values(&mass, &mut v);
        mass.values = v;
    }
    elmer_eigen::log_summary(&result, cfg.tol, log);
    if normalize != Normalize::None {
        elmer_eigen::scale_eigen_vectors(Mass::Consistent(&mass), &mut result.vectors, normalize == Normalize::Unity);
    }
    sys.mass = Some(mass);

    // Norm = ComputeNorm(Solver, n, x) of the unchanged solution vector.
    let var = sim.vars.by_id_mut(id);
    let spec = crate::norm::NormSpec::from_params(&p, var.dofs);
    let (vals, m) = crate::norm::select_norm_entries(&spec, &var.values, n.min(var.values.len()), var.perm.as_deref(), mesh.n_nodes(), log);
    let norm = crate::norm::compute_norm(&spec, &vals, m, comm);
    var.conv.norm = norm;
    var.eigen_values = result.values;
    var.eigen_vectors = result.vectors;
    log.info(
        "SolveLinearSystem",
        &format!("Repointing {} eigenvalue components for: {}", var.dofs, var.name),
        3,
    );
    Ok(norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Value;

    fn full3() -> CsrMatrix<f64> {
        CsrMatrix {
            nrows: 3,
            ncols: 3,
            row_ptr: vec![0, 3, 6, 9],
            col_idx: vec![0, 1, 2, 0, 1, 2, 0, 1, 2],
            values: vec![1.0; 9],
            diag: None,
        }
    }

    /// The ported ARPACK path on `K = tridiag(-1, 2, -1)`, `M = I` (n = 30)
    /// with a shift: the 4 eigenvalues nearest `sigma` are the analytic
    /// `2 - 2 cos(k pi / (n+1))`, in `Eigen System Sorting` order, with
    /// unit-M-norm vectors, and `K` is the same after the shift is undone.
    #[test]
    fn arpack_path_finds_the_eigenpairs_nearest_the_shift() {
        let n = 30usize;
        let mut k = CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: vec![], values: vec![], diag: None };
        let mut m = CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: vec![], values: vec![], diag: None };
        for i in 0..n {
            for j in i.saturating_sub(1)..(i + 2).min(n) {
                k.col_idx.push(j as u32);
                k.values.push(if i == j { 2.0 } else { -1.0 });
                m.col_idx.push(j as u32);
                m.values.push(if i == j { 1.0 } else { 0.0 });
            }
            k.row_ptr.push(k.values.len());
            m.row_ptr.push(m.values.len());
        }
        let k0 = k.values.clone();
        let mut p = ParamMap::new();
        p.insert("Eigen System Values", Value::Integer(vec![4]));
        p.insert("Eigen System Shift", Value::real(0.5));
        p.insert("Eigen System Convergence Tolerance", Value::real(1e-12));
        let cfg = EigenConfig::from_params(&p, true).unwrap();
        // Dense Gaussian elimination as the inner solver.
        let mut inner = |a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]| -> Result<(), SolverError> {
            let mut d = vec![0.0; n * n];
            for i in 0..n {
                for q in a.row_ptr[i]..a.row_ptr[i + 1] {
                    d[i * n + a.col_idx[q] as usize] = a.values[q];
                }
            }
            let mut r = b.to_vec();
            for c in 0..n {
                for i in c + 1..n {
                    let f = d[i * n + c] / d[c * n + c];
                    for j in c..n {
                        d[i * n + j] -= f * d[c * n + j];
                    }
                    r[i] -= f * r[c];
                }
            }
            for i in (0..n).rev() {
                let mut t = r[i];
                for j in i + 1..n {
                    t -= d[i * n + j] * x[j];
                }
                x[i] = t / d[i * n + i];
            }
            Ok(())
        };
        let mut arp = elmer_eigen::arpack::saupd::Arpack::default();
        let r = arpack_eigen_solve(&mut arp, &mut k, &m, 1, &cfg, true, &mut inner).unwrap();
        let mut want: Vec<f64> = (1..=n).map(|j| 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos()).collect();
        want.sort_by(|a, b| (a - 0.5).abs().total_cmp(&(b - 0.5).abs()));
        let mut got: Vec<f64> = r.values.iter().map(|z| z.re).collect();
        let mut near = want[..4].to_vec();
        near.sort_by(f64::total_cmp);
        got.sort_by(f64::total_cmp);
        for (g, w) in got.iter().zip(&near) {
            assert!((g - w).abs() < 1e-10, "{g} vs {w}");
        }
        for v in &r.vectors {
            let nrm: f64 = v.iter().map(|z| z.re * z.re).sum();
            assert!((nrm - 1.0).abs() < 1e-10, "{nrm}");
        }
        for (a, b) in k.values.iter().zip(&k0) {
            assert!((a - b).abs() <= 1e-15 * b.abs().max(1.0));
        }
        // The seed moved on, as dgetv0's SAVE does.
        assert_ne!(arp.iseed, [1, 3, 5, 7]);
    }

    /// Without `Linear System Symmetric` the same solve runs DNAUPD /
    /// DNEUPD: a tridiagonal `K` with sub- and superdiagonals -1.3 and -0.7
    /// has the real eigenvalues `2 - 2 sqrt(0.91) cos(j pi / (n+1))`; the
    /// four nearest the shift come back with `K x = lambda x`.
    #[test]
    fn nonsymmetric_arpack_path_finds_the_eigenpairs_nearest_the_shift() {
        let n = 30;
        let mut k = CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: vec![], values: vec![], diag: None };
        let mut m = k.clone();
        for i in 0..n {
            for j in i.saturating_sub(1)..(i + 2).min(n) {
                k.col_idx.push(j as u32);
                k.values.push(if i == j { 2.0 } else if j < i { -1.3 } else { -0.7 });
                m.col_idx.push(j as u32);
                m.values.push(if i == j { 1.0 } else { 0.0 });
            }
            k.row_ptr.push(k.values.len());
            m.row_ptr.push(m.values.len());
        }
        let k0 = k.clone();
        let mut p = ParamMap::new();
        p.insert("Eigen System Values", Value::Integer(vec![4]));
        p.insert("Eigen System Shift", Value::real(0.5));
        p.insert("Eigen System Convergence Tolerance", Value::real(1e-12));
        let cfg = EigenConfig::from_params(&p, false).unwrap();
        let mut inner = |a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]| -> Result<(), SolverError> {
            let mut d = vec![0.0; n * n];
            for i in 0..n {
                for q in a.row_ptr[i]..a.row_ptr[i + 1] {
                    d[i * n + a.col_idx[q] as usize] = a.values[q];
                }
            }
            let mut r = b.to_vec();
            for c in 0..n {
                for i in c + 1..n {
                    let f = d[i * n + c] / d[c * n + c];
                    for j in c..n {
                        d[i * n + j] -= f * d[c * n + j];
                    }
                    r[i] -= f * r[c];
                }
            }
            for i in (0..n).rev() {
                let mut t = r[i];
                for j in i + 1..n {
                    t -= d[i * n + j] * x[j];
                }
                x[i] = t / d[i * n + i];
            }
            Ok(())
        };
        let mut arp = elmer_eigen::arpack::saupd::Arpack::default();
        let r = arpack_eigen_solve(&mut arp, &mut k, &m, 1, &cfg, false, &mut inner).unwrap();
        let c = 2.0 * 0.91f64.sqrt();
        let mut want: Vec<f64> = (1..=n).map(|j| 2.0 - c * (j as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos()).collect();
        want.sort_by(|a, b| (a - 0.5).abs().total_cmp(&(b - 0.5).abs()));
        let mut near = want[..4].to_vec();
        near.sort_by(f64::total_cmp);
        let mut got: Vec<f64> = r.values.iter().map(|z| z.re).collect();
        got.sort_by(f64::total_cmp);
        for (g, w) in got.iter().zip(&near) {
            assert!((g - w).abs() < 1e-9, "{g} vs {w}");
        }
        for (lam, x) in r.values.iter().zip(&r.vectors) {
            assert!(lam.im.abs() < 1e-12);
            let mut res = 0.0f64;
            for i in 0..n {
                let kx: f64 = (k0.row_ptr[i]..k0.row_ptr[i + 1]).map(|q| k0.values[q] * x[k0.col_idx[q] as usize].re).sum();
                res = res.max((kx - lam.re * x[i].re).abs());
            }
            assert!(res < 1e-8, "residual {res}");
        }
        // K is put back after the shift.
        assert_eq!(k.values, k0.values);
    }

    /// Dirichlet dof 1: its mass row is zero; with `Linear System
    /// Symmetric` its column too (`CRS_SetSymmDirichlet`).
    #[test]
    fn dirichlet_rows_and_symmetric_columns_of_the_mass_are_zero() {
        let mut m = full3();
        dirichlet_mass(&mut m, &[(1, 0.5)], &ParamMap::new()).unwrap();
        assert_eq!(m.values, [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let mut m = full3();
        let mut p = ParamMap::new();
        p.insert("Linear System Symmetric", Value::Logical(true));
        dirichlet_mass(&mut m, &[(1, 0.5)], &p).unwrap();
        assert_eq!(m.values, [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]);
    }
}
