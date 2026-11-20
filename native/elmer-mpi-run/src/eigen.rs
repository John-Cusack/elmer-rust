//! The eigen solve of a partitioned run: `ScaleLinearSystem` (SolveCore.F90:
//! 1420), `ParallelArpackEigenSolve` (ParallelEigenSolve.F90:82) and
//! `BackScaleLinearSystem(EigenScaling)`.
//!
//! The Arnoldi vectors hold this rank's owned dofs, in local order
//! (`ParallelVector`). PARPACK (`PDSAUPD`/`PDSEUPD` for a symmetric
//! system, `PDNAUPD`/`PDNEUPD` otherwise, `elmer_eigen::arpack`) asks for:
//!
//! * `y = M x`: the split mass matrix (`MGmv`, `ParallelMatrixVector` with
//!   the mass values);
//! * `y = OP x`: `(K - sigma M) y = b` by `ParallelIter`, with `b` and the
//!   starting guess `PartitionVector`ed from `workd` (zero at the dofs owned
//!   elsewhere), the solve unscaled (`ParallelIter` is called directly, not
//!   through `SolveLinearSystem`).
//!
//! The eigenvalues are ordered by `SortC` (increasing magnitude), not by
//! `Eigen System Sorting`, and each vector is spread back over the local
//! dofs (`ParallelUpdateResult`).

use elmer_core::{CsrMatrix, LinalgError, Logger, ParamMap, SolveStatus};
use elmer_dist::{dist_solve, DistCsr, DofLayout, Exchange};
use elmer_itermethods::IterMethodsConfig;
use elmer_solver::elmer_eigen::arpack::{naupd, neupd, saupd, seupd};
use elmer_solver::elmer_eigen::{self, ArpackWhich, Complex64, EigenConfig, EigenResult};
use elmer_solver::{ParEigenBackend, SolverError};

use crate::backend::{
    dist_err, dist_options, itermethod_solve, method_of, scaling_of, summed_diagonal_scaling,
    ITERMETHODS,
};
use crate::context::RankContext;

const CALLER: &str = "ParallelEigenSolve";

/// [`ParEigenBackend`] of `elmer-mpi-run`.
pub struct ParEigen;

impl ParEigenBackend for ParEigen {
    fn solve(
        &self,
        params: &ParamMap,
        k: &mut CsrMatrix<f64>,
        mass: &mut CsrMatrix<f64>,
        cfg: &EigenConfig,
        symmetric: bool,
        arp: &mut saupd::Arpack,
        _log: &dyn Logger,
    ) -> Result<EigenResult, SolverError> {
        let ctx = RankContext::current().ok_or_else(|| {
            SolverError::NotWired("parallel eigen solve outside a partitioned run".into())
        })?;
        let layout = ctx.layout_of(k).map_err(dist_err)?;
        let ex: &dyn Exchange = &*ctx.ex;
        // EnforceDirichletConditions zeroed the mass rows of every
        // constrained copy (after CommunicateDirichletBCs), and with
        // `Linear System Symmetric` their columns; the solver zeroed only
        // its own constrained rows.
        let rows = ctx.dirichlet_rows(k);
        if !rows.is_empty() {
            let mut on = vec![false; mass.nrows];
            for &r in &rows {
                if let Some(o) = on.get_mut(r as usize) {
                    *o = true;
                }
            }
            let columns = params.logical("linear system symmetric").unwrap_or(false);
            for r in 0..mass.nrows {
                let row = mass.row_ptr[r]..mass.row_ptr[r + 1];
                if on[r] {
                    mass.values[row].fill(0.0);
                } else if columns {
                    for j in row {
                        if on[mass.col_idx[j] as usize] {
                            mass.values[j] = 0.0;
                        }
                    }
                }
            }
        }
        // ScaleLinearSystem(Solver, A): K and MassValues by the diagonal
        // summed over the holders.
        let d = scaling_of(params).then(|| summed_diagonal_scaling(ex, &layout, k));
        if let Some(d) = &d {
            for m in [&mut *k, &mut *mass] {
                for i in 0..m.nrows {
                    for j in m.row_ptr[i]..m.row_ptr[i + 1] {
                        m.values[j] *= d[i] * d[m.col_idx[j] as usize];
                    }
                }
            }
        }
        let r = arpack_eigen_solve(ex, &layout, params, k, mass, cfg, symmetric, arp);
        // BackScaleLinearSystem(EigenScaling = .TRUE.).
        let Some(d) = &d else { return r };
        for m in [&mut *k, &mut *mass] {
            for i in 0..m.nrows {
                for j in m.row_ptr[i]..m.row_ptr[i + 1] {
                    m.values[j] /= d[i] * d[m.col_idx[j] as usize];
                }
            }
        }
        let mut r = r?;
        for v in &mut r.vectors {
            for (x, s) in v.iter_mut().zip(d) {
                *x *= *s;
            }
        }
        Ok(r)
    }
}

/// `ParallelArpackEigenSolve` on the (scaled) local pieces `k`, `mass`.
#[allow(clippy::too_many_arguments)]
fn arpack_eigen_solve(
    ex: &dyn Exchange,
    layout: &DofLayout,
    params: &ParamMap,
    k: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    cfg: &EigenConfig,
    symmetric: bool,
    arp: &mut saupd::Arpack,
) -> Result<EigenResult, SolverError> {
    use ArpackWhich as W;
    let n = k.nrows;
    let pn = layout.n_owned();
    let (nev, ncv) = (cfg.nev, cfg.ncv);
    if ncv <= nev {
        return Err(SolverError::fatal(
            CALLER,
            "Number of Lanczos vectors must exceed the number of eigenvalues.",
        ));
    }
    let wsel = ArpackWhich::for_mode(cfg.select, true);
    let which = match wsel {
        W::LM => saupd::Which::LM,
        W::SM => saupd::Which::SM,
        // PDSAUPD refuses the real/imaginary-part selections (info = -5).
        _ if symmetric => {
            return Err(SolverError::fatal(
                CALLER,
                " Error with DNAUPD, info =           -5",
            ))
        }
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
    // The split mass matrix (ParallelInitSolve), for MGmv.
    let m = DistCsr::from_local(ex, layout, mass).map_err(dist_err)?;
    // The inner solves: ParallelIter, no scaling.
    let mut opts = dist_options(params)?;
    opts.scaling = false;
    let itermethod = if ITERMETHODS.contains(&method_of(params).as_str()) {
        Some(IterMethodsConfig::from_params(params)?)
    } else {
        None
    };
    let sigma = cfg.shift;
    if sigma != 0.0 {
        for (kv, mv) in k.values.iter_mut().zip(&mass.values) {
            *kv -= sigma * mv;
        }
    }
    let mut iter = 1usize;
    let mut resid = vec![0.0; pn];
    let mut v = vec![0.0; pn * ncv];
    let mut workd = vec![0.0; 3 * pn];
    let mut workl = vec![0.0; 3 * ncv * ncv + 6 * ncv];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = cfg.max_iter as i32;
    iparam[6] = 3;
    let mut tol = cfg.tol;
    let mut failure: Option<SolverError> = None;
    let comm: &dyn elmer_core::Communicator = ex;
    let au = {
        let kk: &CsrMatrix<f64> = k;
        // PartitionVector: owned entries in place, zero elsewhere.
        let partition = |own: &[f64]| {
            let mut local = vec![0.0; n];
            for (&i, &x) in layout.owned().iter().zip(own) {
                local[i as usize] = x;
            }
            local
        };
        let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
            let r = (|| -> Result<(), SolverError> {
                match ido {
                    -1 | 1 => {
                        iter += 1;
                        if ido == -1 {
                            // b = M x, then x := b.
                            let x = wd[ip[0]..ip[0] + pn].to_vec();
                            let mut y = vec![0.0; pn];
                            m.matvec(ex, &x, &mut y);
                            wd[ip[1]..ip[1] + pn].copy_from_slice(&y);
                            wd[ip[0]..ip[0] + pn].copy_from_slice(&y);
                        }
                        let mut x = partition(&wd[ip[1]..ip[1] + pn]);
                        let src = if ido == 1 { ip[2] } else { ip[0] };
                        let b = partition(&wd[src..src + pn]);
                        let stats = match &itermethod {
                            None => dist_solve(ex, layout, kk, &b, &mut x, &opts, None)
                                .map_err(dist_err)?,
                            Some(c) => itermethod_solve(ex, layout, kk, &b, &mut x, &opts, c)?,
                        };
                        if opts.huti.abort_not_converged
                            && matches!(
                                stats.status,
                                SolveStatus::MaxIterations | SolveStatus::Diverged
                            )
                        {
                            return Err(SolverError::Linalg(LinalgError::Breakdown {
                                method: "IterSolve",
                                detail: format!(
                                    "failed convergence tolerances after {} iterations",
                                    stats.iterations
                                ),
                            }));
                        }
                        // ParallelVector.
                        wd[ip[1]..ip[1] + pn].copy_from_slice(&layout.gather_owned(&x));
                        Ok(())
                    }
                    2 => {
                        let x = wd[ip[0]..ip[0] + pn].to_vec();
                        let mut y = vec![0.0; pn];
                        m.matvec(ex, &x, &mut y);
                        wd[ip[1]..ip[1] + pn].copy_from_slice(&y);
                        Ok(())
                    }
                    _ => Err(SolverError::fatal(
                        CALLER,
                        format!("unexpected ARPACK request ido = {ido}"),
                    )),
                }
            })();
            r.map_err(|e| {
                let msg = e.to_string();
                failure = Some(e);
                msg
            })
        };
        if symmetric {
            saupd::pdsaupd(
                comm, &mut op, arp, true, pn, which, nev, &mut tol, &mut resid, ncv, &mut v, pn,
                iparam, &mut workd, &mut workl, 0,
            )
            .map(|a| (a.iparam, a.info, Some(a), None))
        } else {
            naupd::pdnaupd(
                comm, &mut op, arp, true, pn, nwhich, nev, &mut tol, &mut resid, ncv, &mut v, pn,
                iparam, &mut workd, &mut workl, 0,
            )
            .map(|a| (a.iparam, a.info, None, Some(a)))
        }
    };
    let (iparam_out, kinfo, sau, nau) = match au {
        Ok(a) => a,
        Err(msg) => return Err(failure.unwrap_or_else(|| SolverError::fatal(CALLER, msg))),
    };
    if kinfo != 0 {
        let name = if symmetric { "DSAUPD" } else { "DNAUPD" };
        return Err(SolverError::fatal(
            CALLER,
            format!(" Error with {name}, info = {kinfo:>11}"),
        ));
    }
    // D(:,1), D(:,2) of pDSEUPD / pDNEUPD.
    let (dr, di): (Vec<f64>, Vec<f64>) = if let Some(au) = sau {
        let (d, ierr) = seupd::pdseupd(
            comm, sigma, true, pn, which, nev, tol, &resid, ncv, &mut v, pn, &au, &mut workd,
            &mut workl,
        );
        if ierr != 0 {
            return Err(SolverError::fatal(
                CALLER,
                format!(" Error with DSEUPD, info = {ierr:>11}"),
            ));
        }
        let d: Vec<f64> = (0..nev).map(|i| d.get(i).copied().unwrap_or(0.0)).collect();
        (d, vec![0.0; nev])
    } else {
        let au = nau.expect("the non-symmetric path ran PDNAUPD");
        let r = neupd::pdneupd(
            comm, sigma, 0.0, true, pn, nwhich, nev, tol, &resid, ncv, &mut v, pn, &au, &mut workd,
            &mut workl,
        );
        if r.info != 0 {
            return Err(SolverError::fatal(
                CALLER,
                format!(" Error with DNEUPD, info = {:>11}", r.info),
            ));
        }
        (r.dr[..nev].to_vec(), r.di[..nev].to_vec())
    };
    // EigValues, ordered by SortC.
    let mut values: Vec<Complex64> = (0..nev).map(|i| Complex64::new(dr[i], di[i])).collect();
    let mut perm: Vec<usize> = (0..nev).collect();
    sort_c(&mut values, &mut perm);
    if sigma != 0.0 {
        for (kv, mv) in k.values.iter_mut().zip(&mass.values) {
            *kv += sigma * mv;
        }
    }
    // The vector of value p is column k = 1 + sum over j < p of (1 if
    // D(j,2) == 0 else 2) of V, a complex one V(:,k) + i V(:,k+1).
    let col = |c: usize| -> Result<&[f64], SolverError> {
        v.get(c * pn..c * pn + pn).ok_or_else(|| {
            SolverError::fatal(CALLER, "eigenvector column outside the Arnoldi basis")
        })
    };
    let mut vectors = Vec::with_capacity(nev);
    for &p in &perm {
        let k0: usize = (0..p).map(|j| if di[j] == 0.0 { 1 } else { 2 }).sum();
        let mut re = vec![0.0; n];
        layout.owned_to_local(ex, col(k0)?, &mut re);
        let mut im = vec![0.0; n];
        if di[p] != 0.0 {
            layout.owned_to_local(ex, col(k0 + 1)?, &mut im);
        } else {
            layout.owned_to_local(ex, &vec![0.0; pn], &mut im);
        }
        vectors.push(
            re.iter()
                .zip(&im)
                .map(|(&a, &b)| Complex64::new(a, b))
                .collect(),
        );
    }
    Ok(EigenResult {
        values,
        vectors,
        iterations: iparam_out[2].max(0) as usize,
        op_applications: iter - 1,
        converged: iparam_out[4].max(0) as usize,
        residuals: None,
    })
}

/// `SortC(n, a, b)` (GeneralUtils.F90:724): heap sort of `a` by increasing
/// magnitude, `b` along.
fn sort_c(a: &mut [Complex64], b: &mut [usize]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    // 1-based, as the Fortran.
    let mut l = n / 2 + 1;
    let mut ir = n;
    loop {
        let (ra, rb);
        if l > 1 {
            l -= 1;
            ra = a[l - 1];
            rb = b[l - 1];
        } else {
            ra = a[ir - 1];
            rb = b[ir - 1];
            a[ir - 1] = a[0];
            b[ir - 1] = b[0];
            ir -= 1;
            if ir == 1 {
                a[0] = ra;
                b[0] = rb;
                return;
            }
        }
        let mut i = l;
        let mut j = l + l;
        while j <= ir {
            if j < ir && a[j - 1].norm() < a[j].norm() {
                j += 1;
            }
            if ra.norm() < a[j - 1].norm() {
                a[i - 1] = a[j - 1];
                b[i - 1] = b[j - 1];
                i = j;
                j += i;
            } else {
                j = ir + 1;
            }
            a[i - 1] = ra;
            b[i - 1] = rb;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sort_c_orders_by_magnitude() {
        let mut a: Vec<Complex64> = [3.0, -1.0, 2.0, -5.0, 4.0]
            .iter()
            .map(|&x| Complex64::new(x, 0.0))
            .collect();
        let mut b: Vec<usize> = (0..5).collect();
        sort_c(&mut a, &mut b);
        let m: Vec<f64> = a.iter().map(|z| z.re).collect();
        assert_eq!(m, vec![-1.0, 2.0, 3.0, 4.0, -5.0]);
        assert_eq!(b, vec![1, 2, 0, 4, 3]);
    }
}
