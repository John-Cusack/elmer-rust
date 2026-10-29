//! BiCGStab(l) with right preconditioning: `RealBiCGStabl`
//! (`IterativeMethods.F90:693`) and `ComplexBiCGStabl` (`:2722`), after
//! Fokkema's zbistbl. One generic body; the places where the two Fortran
//! routines differ are branches on `T::IS_COMPLEX`.

use crate::config::{IterMethodsConfig, PseudoComplex, RobustConfig};
use crate::ops::{axpy, check_dims, diverged, dot, nrm2, pcond, residual};
use crate::scalar::IterScalar;
use crate::small::{lu_factor, lu_solve};
use crate::{finish, History};
use elmer_core::{
    Communicator, IterativeSolver, LinalgError, LinearOperator, Preconditioner, SolveStats,
};

pub struct BiCgStabL<T: IterScalar> {
    cfg: IterMethodsConfig,
    history: History,
    messages: Vec<String>,
    _t: std::marker::PhantomData<T>,
}

impl<T: IterScalar> BiCgStabL<T> {
    pub fn new(cfg: &IterMethodsConfig) -> Self {
        BiCgStabL {
            cfg: cfg.clone(),
            history: History::default(),
            messages: Vec::new(),
            _t: std::marker::PhantomData,
        }
    }

    /// Relative residual after each outer round of the last solve.
    pub fn history(&self) -> &[f64] {
        &self.history.0
    }

    /// Warnings the Fortran prints through `Warn` (halts and breakdowns).
    pub fn messages(&self) -> &[String] {
        &self.messages
    }
}

const NAME: &str = "BiCGStab(l)";

impl<T: IterScalar> IterativeSolver<T> for BiCgStabL<T> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<T>,
        m: &dyn Preconditioner<T>,
        b: &[T],
        x: &mut [T],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims(NAME, a, b, x)?;
        let cfg = &self.cfg;
        let l = cfg.bicgstabl_degree;
        if l < 2 {
            return Err(LinalgError::Unsupported(
                "BiCGStab(l): polynomial degree < 2".into(),
            ));
        }
        if cfg.pseudo_complex != PseudoComplex::Off && !T::IS_COMPLEX {
            return Err(LinalgError::Unsupported(
                "BiCGStab(l): pseudo-complex storage is only supported by GCR".into(),
            ));
        }
        self.history.0.clear();
        self.messages.clear();
        let cplx = T::IS_COMPLEX;
        let n = b.len();
        let tol = cfg.tolerance;
        let max_tol = cfg.divergence_limit;
        let breakdown = |what: &str| LinalgError::Breakdown {
            method: NAME,
            detail: what.to_string(),
        };

        // IterSolve: an all-zero guess becomes 1e-8 for the BiCGStab family.
        if let Some(fill) = cfg.zero_guess_fill {
            if x.iter().all(|v| *v == T::zero()) {
                x.fill(T::from_parts(fill, fill));
            }
        }
        if x.iter().all(|v| *v == T::zero()) {
            x.copy_from_slice(b);
        }

        // Column layout of `w` (the Fortran `work(n, 3+2*(l+1))`), 0-based.
        let (rr, r, u, xp, bp) = (0, 1, l + 2, 2 * l + 3, 2 * l + 4);
        let mut w = vec![vec![T::zero(); n]; 2 * l + 5];
        let mut t = vec![T::zero(); n];
        // Small Gram matrix `rw[col][row]` over the l+1 residuals, plus y0, yl, y.
        let mut gram = vec![vec![T::zero(); l + 1]; l + 1];
        let mut y0 = vec![T::zero(); l + 1];
        let mut yl = vec![T::zero(); l + 1];
        let mut yv = vec![T::zero(); l + 1];

        residual(a, x, b, &mut w[r]);
        let bnrm = nrm2(comm, b);
        let rnrm0 = nrm2(comm, &w[r]);

        if !cplx {
            if bnrm.is_nan() {
                return Err(breakdown("bnrm = NaN"));
            }
            if rnrm0.is_nan() {
                return Err(breakdown("nrm0 = NaN"));
            }
        }
        if bnrm == 0.0 {
            // Zero right-hand side: the solution is x = 0.
            x.fill(T::zero());
            return Ok(finish(0, 0.0, true, false, None));
        }
        let mut errorind = rnrm0 / bnrm;
        if !cplx && errorind.is_nan() {
            return Err(breakdown("errorind = NaN"));
        }
        let mut conv = errorind < tol;
        let mut div = if cplx {
            diverged(errorind, max_tol)
        } else {
            errorind > max_tol
        };
        if conv || div {
            return Ok(finish(0, errorind, conv, div, None));
        }

        let r0 = w[r].clone();
        w[rr].copy_from_slice(&r0);
        w[bp].copy_from_slice(&r0);
        w[xp].copy_from_slice(x);
        x.fill(T::zero());

        let mut rnrm: f64;
        let mut mxnrmx = rnrm0;
        let mut mxnrmr = rnrm0;
        let mut alpha = T::zero();
        let mut omega = T::one();
        let mut rho0 = T::one();
        let delta = 1.0e-2;

        let robust: Option<RobustConfig> = if cplx { None } else { cfg.robust };
        let mut best_norm = f64::MAX.sqrt();
        let mut bad_count = 0usize;
        let mut best_x = vec![T::zero(); if robust.is_some() { n } else { 0 }];

        let mut halted: Option<String> = None;
        let mut early_exit = false;
        let mut rounds_done = 0usize;

        'rounds: for rd in 1..=cfg.max_iterations {
            rounds_done = rd;
            // --- The BiCG part ---
            rho0 = -omega * rho0;
            for k in 1..=l {
                let rho1 = dot(comm, &w[rr], &w[r + k - 1]);
                if rho0 == T::zero() {
                    if cplx {
                        return Err(breakdown("Breakdown error 1."));
                    }
                    let msg = "Iteration halted: rho0 == zero.".to_string();
                    self.messages.push(msg.clone());
                    halted = Some(msg);
                    break 'rounds;
                }
                if !cplx && rho1.is_nan() {
                    return Err(breakdown("rho1 == NaN"));
                }
                let beta = alpha * (rho1 / rho0);
                rho0 = rho1;
                for j in 0..k {
                    for i in 0..n {
                        w[u + j][i] = w[r + j][i] - beta * w[u + j][i];
                    }
                }
                pcond(m, &w[u + k - 1], &mut t);
                a.apply(&t, &mut w[u + k]);
                let sigma = dot(comm, &w[rr], &w[u + k]);
                if sigma == T::zero() {
                    if cplx {
                        return Err(breakdown("Breakdown error 2."));
                    }
                    let msg = "Iteration halted: sigma == zero.".to_string();
                    self.messages.push(msg.clone());
                    halted = Some(msg);
                    break 'rounds;
                }
                if !cplx && sigma.is_nan() {
                    return Err(breakdown("sigma == NaN"));
                }
                alpha = rho1 / sigma;
                axpy(alpha, &w[u], x);
                for j in 0..k {
                    for i in 0..n {
                        w[r + j][i] = w[r + j][i] - alpha * w[u + j + 1][i];
                    }
                }
                pcond(m, &w[r + k - 1], &mut t);
                a.apply(&t, &mut w[r + k]);
                rnrm = nrm2(comm, &w[r]);
                if !cplx && rnrm.is_nan() {
                    return Err(breakdown("rnrm == NaN"));
                }
                mxnrmx = mxnrmx.max(rnrm);
                mxnrmr = mxnrmr.max(rnrm);

                // A few BiCG updates may already be enough.
                errorind = rnrm / bnrm;
                conv = errorind < tol;
                if !cplx {
                    div = errorind.is_nan();
                }
                if conv || div {
                    early_exit = true;
                    break;
                }
            }
            if early_exit {
                break 'rounds;
            }

            // --- The convex polynomial part ---
            for i in 0..=l {
                for j in 0..=i {
                    gram[j][i] = dot(comm, &w[r + i], &w[r + j]);
                }
            }
            for j in 1..=l {
                for i in 0..j {
                    gram[j][i] = gram[i][j].conj();
                }
            }
            // gram[col][row] holds rwork(row, col); the (l-1)x(l-1) block is
            // rows/cols 1..l, copied before factorisation.
            let mut lu: Vec<Vec<T>> = (1..l)
                .map(|i| (1..l).map(|j| gram[j][i]).collect())
                .collect();
            let piv = lu_factor(&mut lu);

            y0[0] = -T::one();
            y0[1..l].copy_from_slice(&gram[0][1..l]);
            lu_solve(&lu, &piv, &mut y0[1..l]);
            y0[l] = T::zero();

            yl[0] = T::zero();
            yl[1..l].copy_from_slice(&gram[l][1..l]);
            lu_solve(&lu, &piv, &mut yl[1..l]);
            yl[l] = -T::one();

            // --- Convex combination ---
            // y = Z v with Z the Hermitian (symmetric) Gram matrix.
            let zmul = |v: &[T], out: &mut [T]| {
                for (i, o) in out.iter_mut().enumerate() {
                    let mut s = T::zero();
                    for (j, vj) in v.iter().enumerate() {
                        s = s + gram[j][i] * *vj;
                    }
                    *o = s;
                }
            };
            let small_dot = |p: &[T], q: &[T]| {
                let mut s = T::zero();
                for (pi, qi) in p.iter().zip(q) {
                    s = s + pi.conj() * *qi;
                }
                s
            };

            zmul(&y0, &mut yv);
            let q0 = small_dot(&y0, &yv);
            let kappa0;
            let kappal;
            if cplx {
                kappa0 = q0.abs().sqrt();
                zmul(&yl, &mut yv);
                kappal = small_dot(&yl, &yv).abs().sqrt();
                // zero or NaN kappa: the Fortran warns and leaves the loop
                if !(kappa0 > 0.0) || !(kappal > 0.0) {
                    self.messages
                        .push("kappa^2 is zero or NaN, iteration halted".into());
                    break 'rounds;
                }
            } else {
                // `.NOT. (x > 0)` also catches NaN
                if !(q0.re() > 0.0) {
                    let msg = "kappa0^2 is non-positive, iteration halted".to_string();
                    self.messages.push(msg.clone());
                    halted = Some(msg);
                    break 'rounds;
                }
                kappa0 = q0.re().sqrt();
                zmul(&yl, &mut yv);
                let ql = small_dot(&yl, &yv).re();
                if !(ql > 0.0) {
                    let msg = "kappal^2 is non-positive, iteration halted".to_string();
                    self.messages.push(msg.clone());
                    halted = Some(msg);
                    break 'rounds;
                }
                kappal = ql.sqrt();
            }

            zmul(&y0, &mut yv);
            let varrho = small_dot(&yl, &yv).unscale(kappa0 * kappal);
            let avar = varrho.abs();
            // SIGN(1, varrho) real; unit direction of varrho (1 when zero) complex
            let hatgamma = if avar > 0.0 {
                varrho.unscale(avar).scale(avar.max(0.7) * kappa0 / kappal)
            } else {
                T::from_f64(0.7 * kappa0 / kappal)
            };
            for i in 0..=l {
                y0[i] = y0[i] - hatgamma * yl[i];
            }

            // --- Update ---
            omega = y0[l];
            for j in 1..=l {
                let c = y0[j];
                for i in 0..n {
                    w[u][i] = w[u][i] - c * w[u + j][i];
                }
                axpy(c, &w[r + j - 1], x);
                for i in 0..n {
                    w[r][i] = w[r][i] - c * w[r + j][i];
                }
            }

            zmul(&y0, &mut yv);
            let q = small_dot(&y0, &yv);
            if cplx {
                rnrm = q.abs().sqrt();
            } else {
                if q.re() < 0.0 {
                    let msg = "rnrm^2 is negative, iteration halted".to_string();
                    self.messages.push(msg.clone());
                    halted = Some(msg);
                    break 'rounds;
                }
                rnrm = q.re().sqrt();
            }

            // --- The reliable update part ---
            mxnrmx = mxnrmx.max(rnrm);
            mxnrmr = mxnrmr.max(rnrm);
            let xpdt = rnrm < delta * rnrm0 && rnrm0 < mxnrmx;
            let rcmp = (rnrm < delta * mxnrmr && rnrm0 < mxnrmr) || xpdt;
            if rcmp {
                pcond(m, x, &mut t);
                a.apply(&t, &mut w[r]);
                mxnrmr = rnrm;
                for i in 0..n {
                    w[r][i] = w[bp][i] - w[r][i];
                }
                if xpdt {
                    for i in 0..n {
                        w[xp][i] = w[xp][i] + t[i];
                    }
                    x.fill(T::zero());
                    let (lo, hi) = w.split_at_mut(bp);
                    hi[0].copy_from_slice(&lo[r]);
                    mxnrmx = rnrm;
                }
            }

            errorind = rnrm / bnrm;
            self.history.push(errorind);

            if let Some(rb) = robust {
                if rd >= rb.start_iteration {
                    if errorind < rb.step_size * best_norm {
                        best_norm = errorind;
                        best_x.copy_from_slice(x);
                        bad_count = 0;
                    } else {
                        bad_count += 1;
                    }
                    if best_norm < rb.tolerance
                        && (errorind > rb.max_tolerance || bad_count > rb.max_bad_iterations)
                    {
                        break 'rounds;
                    }
                }
            }

            conv = errorind < tol;
            div = diverged(errorind, max_tol);
            if conv || div {
                break 'rounds;
            }
        }
        if let Some(rb) = robust {
            if best_norm < rb.tolerance {
                conv = true;
            }
            if best_norm < errorind {
                x.copy_from_slice(&best_x);
            }
        }

        // We solved z = P x; recover the true unknown.
        t.copy_from_slice(x);
        pcond(m, &t, x);
        axpy(T::one(), &w[xp], x);

        Ok(finish(
            rounds_done.min(cfg.max_iterations),
            errorind,
            conv,
            div,
            halted,
        ))
    }
}
