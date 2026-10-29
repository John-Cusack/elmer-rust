//! IDR(s) with right preconditioning, after van Gijzen: `RealIDRS`
//! (`IterativeMethods.F90:1590`) and `ComplexIDRS` (`:3118`). The real variant
//! adds residual smoothing and the robust mode.

use crate::config::{IterMethodsConfig, PseudoComplex, RobustConfig};
use crate::ops::{axmy, axpy, check_dims, diverged, dot, nrm2, pcond, residual};
use crate::rng::Rng;
use crate::scalar::IterScalar;
use crate::{finish, History, StopFn};
use elmer_core::{
    Communicator, IterativeSolver, LinalgError, LinearOperator, Preconditioner, SolveStats,
};

pub struct Idrs<T: IterScalar> {
    cfg: IterMethodsConfig,
    stop: Option<StopFn<T>>,
    history: History,
}

impl<T: IterScalar> Idrs<T> {
    pub fn new(cfg: &IterMethodsConfig) -> Self {
        Idrs {
            cfg: cfg.clone(),
            stop: None,
            history: History::default(),
        }
    }

    /// Uses `f(x, b, r)` instead of `||r|| / ||b||`. With smoothing it sees the
    /// smoothed iterate and residual.
    pub fn with_stop_criterion(mut self, f: StopFn<T>) -> Self {
        self.stop = Some(f);
        self
    }

    /// Relative residual after each iteration of the last solve.
    pub fn history(&self) -> &[f64] {
        &self.history.0
    }
}

const NAME: &str = "IDR(s)";

impl<T: IterScalar> IterativeSolver<T> for Idrs<T> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<T>,
        mm: &dyn Preconditioner<T>,
        b: &[T],
        x: &mut [T],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims(NAME, a, b, x)?;
        let cfg = &self.cfg;
        let s = cfg.idrs_s;
        if s < 1 {
            return Err(LinalgError::Unsupported("IDR(s): s must be at least 1".into()));
        }
        if cfg.pseudo_complex != PseudoComplex::Off && !T::IS_COMPLEX {
            return Err(LinalgError::Unsupported(
                "IDR(s): pseudo-complex storage is only supported by GCR".into(),
            ));
        }
        self.history.0.clear();
        let cplx = T::IS_COMPLEX;
        let n = b.len();
        let tol = cfg.tolerance;
        let max_tol = cfg.divergence_limit;
        let max_rounds = cfg.max_iterations;
        let smoothing = cfg.idrs_smoothing && !cplx;
        let robust: Option<RobustConfig> = if cplx { None } else { cfg.robust };

        let mut r = vec![T::zero(); n];
        let mut v = vec![T::zero(); n];
        let mut t = vec![T::zero(); n];

        let normb = nrm2(comm, b);
        residual(a, x, b, &mut r);
        let stop_val = |stop: &Option<StopFn<T>>, x: &[T], r: &[T]| match stop {
            Some(f) => f(x, b, r),
            None => nrm2(comm, r) / normb,
        };
        let mut errorind = stop_val(&self.stop, x, &r);
        let mut conv = errorind < tol;
        let mut div = diverged(errorind, max_tol);
        if conv || div {
            return Ok(finish(0, errorind, conv, div, None));
        }

        let (mut r_s, mut x_s) = if smoothing {
            (r.clone(), x.to_vec())
        } else {
            (Vec::new(), Vec::new())
        };

        // Shadow space P: random, then orthonormalised.
        let mut rng = Rng::new(cfg.random_seed);
        let mut p = vec![vec![T::zero(); n]; s];
        for col in p.iter_mut() {
            for e in col.iter_mut() {
                *e = T::from_parts(rng.uniform(), 0.0);
            }
        }
        if cplx {
            // P = Pr + i*Pi, with Pi drawn after all of Pr.
            for col in p.iter_mut() {
                for e in col.iter_mut() {
                    *e = T::from_parts(e.re(), rng.uniform());
                }
            }
        }
        for j in 0..s {
            for k in 0..j {
                let alpha = dot(comm, &p[k], &p[j]);
                let (lo, hi) = p.split_at_mut(j);
                axmy(alpha, &lo[k], &mut hi[0]);
            }
            let nrm = nrm2(comm, &p[j]);
            for e in p[j].iter_mut() {
                *e = e.unscale(nrm);
            }
        }
        let kappa = 0.7;

        let mut g = vec![vec![T::zero(); n]; s];
        let mut u = vec![vec![T::zero(); n]; s];
        let mut mmat = vec![vec![T::zero(); s]; s];
        let mut f = vec![T::zero(); s];
        let mut mu = vec![T::zero(); s];
        let mut alpha = vec![T::zero(); s];
        let mut beta = vec![T::zero(); s];
        let mut gamma = vec![T::zero(); s];
        let mut om = T::one();
        let mut iter = 0usize;
        let mut jj = 0usize;

        let mut best_norm = f64::MAX.sqrt();
        let mut bad_count = 0usize;
        let mut best_x = vec![T::zero(); if robust.is_some() { n } else { 0 }];

        // Smoothing step shared by the inner and outer updates.
        let smooth = |r: &[T], x: &[T], r_s: &mut [T], x_s: &mut [T], t: &mut [T]| {
            for i in 0..n {
                t[i] = r_s[i] - r[i];
            }
            let tr_s = dot(comm, t, r_s);
            let tt = dot(comm, t, t);
            let theta = tr_s / tt;
            for i in 0..n {
                r_s[i] = r_s[i] - theta * t[i];
                x_s[i] = x_s[i] - theta * (x_s[i] - x[i]);
            }
        };

        'outer: while !conv && !div {
            // Generate s vectors in G_j.
            for k in 0..s {
                f[k] = dot(comm, &p[k], &r);
            }
            for k in 0..s {
                v.copy_from_slice(&r);
                if jj > 0 {
                    // Solve the small lower-triangular system and make v orthogonal to P.
                    for i in k..s {
                        gamma[i] = f[i];
                        for j in k..i {
                            gamma[i] = gamma[i] - mmat[i][j] * gamma[j];
                        }
                        gamma[i] = gamma[i] / mmat[i][i];
                        axmy(gamma[i], &g[i], &mut v);
                    }
                    pcond(mm, &v, &mut t);
                    for e in t.iter_mut() {
                        *e = om * *e;
                    }
                    for i in k..s {
                        axpy(gamma[i], &u[i], &mut t);
                    }
                    u[k].copy_from_slice(&t);
                } else {
                    pcond(mm, &v, &mut u[k]);
                }
                a.apply(&u[k], &mut g[k]);

                // Bi-orthogonalise the new basis vectors.
                for i in 0..s {
                    mu[i] = dot(comm, &p[i], &g[k]);
                }
                for i in 0..k {
                    alpha[i] = mu[i];
                    for j in 0..i {
                        alpha[i] = alpha[i] - mmat[i][j] * alpha[j];
                    }
                    alpha[i] = alpha[i] / mmat[i][i];
                    let (glo, ghi) = g.split_at_mut(k);
                    axmy(alpha[i], &glo[i], &mut ghi[0]);
                    let (ulo, uhi) = u.split_at_mut(k);
                    axmy(alpha[i], &ulo[i], &mut uhi[0]);
                    for q in k..s {
                        mu[q] = mu[q] - mmat[q][i] * alpha[i];
                    }
                }
                for q in k..s {
                    mmat[q][k] = mu[q];
                }

                // Breakdown?
                if mmat[k][k].abs() <= f64::MIN_POSITIVE {
                    div = true;
                    break;
                }

                // Make r orthogonal to p_1..p_k; update solution and residual.
                beta[k] = f[k] / mmat[k][k];
                axmy(beta[k], &g[k], &mut r);
                axpy(beta[k], &u[k], x);
                for q in k + 1..s {
                    f[q] = f[q] - beta[k] * mmat[q][k];
                }

                if smoothing {
                    smooth(&r, x, &mut r_s, &mut x_s, &mut t);
                }

                iter += 1;
                errorind = if smoothing {
                    stop_val(&self.stop, &x_s, &r_s)
                } else {
                    stop_val(&self.stop, x, &r)
                };
                self.history.push(errorind);
                conv = errorind < tol;
                div = diverged(errorind, max_tol);
                if conv || div || iter == max_rounds {
                    break;
                }
            }
            if conv || div || iter == max_rounds {
                break 'outer;
            }

            // First residual in G_{j+1}: r is already orthogonal to P, so v = r.
            jj += 1;
            pcond(mm, &r, &mut v);
            a.apply(&v, &mut t);

            // Choose omega, 'maintaining the convergence'.
            let nr = nrm2(comm, &r);
            let nt = nrm2(comm, &t);
            let tr = dot(comm, &t, &r);
            let rho = tr.unscale(nt * nr).abs();
            om = tr.unscale(nt * nt);
            if rho < kappa {
                om = om.scale(kappa).unscale(rho);
            }
            if om.abs() <= f64::EPSILON {
                div = true;
                break 'outer;
            }

            axmy(om, &t, &mut r);
            axpy(om, &v, x);
            if smoothing {
                smooth(&r, x, &mut r_s, &mut x_s, &mut t);
            }

            iter += 1;
            errorind = if smoothing {
                stop_val(&self.stop, &x_s, &r_s)
            } else {
                stop_val(&self.stop, x, &r)
            };
            self.history.push(errorind);

            if let Some(rb) = robust {
                // Keep the best iterate so far (with a small margin).
                if errorind < rb.step_size * best_norm {
                    best_norm = errorind;
                    best_x.copy_from_slice(if smoothing { &x_s } else { x });
                    bad_count = 0;
                } else {
                    bad_count += 1;
                }
                if best_norm < rb.tolerance
                    && (errorind > rb.max_tolerance || bad_count > rb.max_bad_iterations)
                {
                    break 'outer;
                }
            }

            conv = errorind < tol;
            div = diverged(errorind, max_tol);
            if iter == max_rounds {
                break 'outer;
            }
        }

        if smoothing {
            x.copy_from_slice(&x_s);
        }
        if let Some(rb) = robust {
            if best_norm < rb.tolerance {
                conv = true;
            }
            if best_norm < errorind {
                x.copy_from_slice(&best_x);
            }
        }
        Ok(finish(iter, errorind, conv, div, None))
    }
}
