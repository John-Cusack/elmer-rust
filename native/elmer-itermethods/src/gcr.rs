//! GCR with restart: `GCR` (`IterativeMethods.F90:1271`) and `GCR_Z` (`:2517`).
//! Right preconditioned; the real variant also supports `Linear System Min
//! Iterations` and the pseudo-complex (real storage of a complex system) mode.

use crate::config::{IterMethodsConfig, PseudoComplex};
use crate::ops::{axmy, axpy, check_dims, diverged, dot, nrm2, pcond, residual};
use crate::scalar::IterScalar;
use crate::{finish, History, StopFn};
use elmer_core::{
    Communicator, IterativeSolver, LinalgError, LinearOperator, Preconditioner, SolveStats,
};

pub struct Gcr<T: IterScalar> {
    cfg: IterMethodsConfig,
    stop: Option<StopFn<T>>,
    history: History,
    messages: Vec<String>,
}

impl<T: IterScalar> Gcr<T> {
    pub fn new(cfg: &IterMethodsConfig) -> Self {
        Gcr {
            cfg: cfg.clone(),
            stop: None,
            history: History::default(),
            messages: Vec::new(),
        }
    }

    /// Uses `f(x, b, r)` instead of `||r|| / ||b||` (HUTI's user stop criterion).
    pub fn with_stop_criterion(mut self, f: StopFn<T>) -> Self {
        self.stop = Some(f);
        self
    }

    /// Relative residual after each iteration of the last solve.
    pub fn history(&self) -> &[f64] {
        &self.history.0
    }

    /// The accuracy warning GCR prints when the iterated and true residuals disagree.
    pub fn messages(&self) -> &[String] {
        &self.messages
    }
}

/// `(Re, Im)` of `x^H y` for a complex vector stored in real arrays: the real
/// part is the plain dot product, the imaginary part `PseudoZDotProd[2]`'s second call.
fn pseudo_dot<T: IterScalar>(mode: PseudoComplex, x: &[T], y: &[T]) -> f64 {
    let mut b = 0.0;
    match mode {
        PseudoComplex::Interleaved => {
            for k in 0..x.len() / 2 {
                b += x[2 * k].re() * y[2 * k + 1].re() - x[2 * k + 1].re() * y[2 * k].re();
            }
        }
        PseudoComplex::Split => {
            let h = x.len() / 2;
            for k in 0..h {
                b += x[k].re() * y[h + k].re() - x[h + k].re() * y[k].re();
            }
        }
        PseudoComplex::Off => {}
    }
    b
}

/// `dst += s * rot(src)` where `rot` multiplies the stored complex vector by `i`:
/// `(re, im) -> (-im, re)`, in the layout `mode`.
fn add_rot<T: IterScalar>(mode: PseudoComplex, s: f64, src: &[T], dst: &mut [T]) {
    let n = src.len();
    match mode {
        PseudoComplex::Interleaved => {
            for k in 0..n / 2 {
                dst[2 * k] = dst[2 * k] - src[2 * k + 1].scale(s);
                dst[2 * k + 1] = dst[2 * k + 1] + src[2 * k].scale(s);
            }
        }
        PseudoComplex::Split => {
            let h = n / 2;
            for k in 0..h {
                dst[k] = dst[k] - src[h + k].scale(s);
                dst[h + k] = dst[h + k] + src[k].scale(s);
            }
        }
        PseudoComplex::Off => {}
    }
}

impl<T: IterScalar> IterativeSolver<T> for Gcr<T> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<T>,
        m: &dyn Preconditioner<T>,
        b: &[T],
        x: &mut [T],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        check_dims("GCR", a, b, x)?;
        let cfg = &self.cfg;
        let cplx = T::IS_COMPLEX;
        // The complex Fortran routine has no pseudo-complex or min-iteration logic.
        let pseudo = if cplx { PseudoComplex::Off } else { cfg.pseudo_complex };
        let min_iter = if cplx { 0 } else { cfg.min_iterations };
        let restart = cfg.gcr_restart;
        if restart < 1 {
            return Err(LinalgError::Unsupported("GCR restart must be at least 1".into()));
        }
        let n = b.len();
        if pseudo != PseudoComplex::Off && n & 1 == 1 {
            return Err(LinalgError::Dim("GCR pseudo-complex needs an even length".into()));
        }
        self.history.0.clear();
        self.messages.clear();
        let tol = cfg.tolerance;
        let max_tol = cfg.divergence_limit;

        let mut r = vec![T::zero(); n];
        let mut t1 = vec![T::zero(); n];
        let mut t2 = vec![T::zero(); n];
        let mut s = vec![vec![T::zero(); n]; restart - 1];
        let mut v = vec![vec![T::zero(); n]; restart - 1];

        residual(a, x, b, &mut r);
        let bnorm = nrm2(comm, b);
        let mut rnorm = nrm2(comm, &r);
        let mut res = match &self.stop {
            Some(f) => f(x, b, &r),
            None => rnorm / bnorm,
        };
        let mut conv = res < tol && min_iter == 0;
        let mut div = diverged(res, max_tol);
        if conv || div {
            return Ok(finish(0, res, conv, div, None));
        }

        let mut iters = 0;
        for k in 1..=cfg.max_iterations {
            iters = k;
            // Restart bookkeeping: j is the position inside the cycle.
            let j = if k % restart == 0 {
                restart
            } else {
                let j = k % restart;
                if j == 1 && k > 1 {
                    // true residual at a restart
                    residual(a, x, b, &mut r);
                }
                j
            };

            pcond(m, &r, &mut t1);
            a.apply(&t1, &mut t2);

            // Orthogonalise against the stored search directions.
            for i in 0..j - 1 {
                let beta = dot(comm, &v[i], &t2);
                if pseudo != PseudoComplex::Off {
                    let bim = pseudo_dot(pseudo, &v[i], &t2);
                    add_rot(pseudo, -bim, &s[i], &mut t1);
                    add_rot(pseudo, -bim, &v[i], &mut t2);
                }
                axmy(beta, &s[i], &mut t1);
                axmy(beta, &v[i], &mut t2);
            }

            let alpha = nrm2(comm, &t2);
            let inv = T::from_f64(1.0 / alpha);
            for e in t1.iter_mut() {
                *e = inv * *e;
            }
            for e in t2.iter_mut() {
                *e = inv * *e;
            }

            // Update the solution and save the search data.
            let beta = dot(comm, &t2, &r);
            if pseudo != PseudoComplex::Off {
                let bim = pseudo_dot(pseudo, &t2, &r);
                add_rot(pseudo, bim, &t1, x);
                add_rot(pseudo, -bim, &t2, &mut r);
            }
            axpy(beta, &t1, x);
            axmy(beta, &t2, &mut r);
            if j != restart {
                s[j - 1].copy_from_slice(&t1);
                v[j - 1].copy_from_slice(&t2);
            }

            rnorm = nrm2(comm, &r);
            res = match &self.stop {
                Some(f) => f(x, b, &r),
                None => rnorm / bnorm,
            };
            self.history.push(res);

            conv = res < tol && k >= min_iter;
            if conv {
                // The iterated residual may drift from the true one.
                let mut tr = vec![T::zero(); n];
                residual(a, x, b, &mut tr);
                let true_norm = nrm2(comm, &tr);
                let err = (true_norm - rnorm).abs() / true_norm;
                if err > 0.1 {
                    self.messages
                        .push("Iterated GCR solution may not be accurate".into());
                }
            }
            div = diverged(res, max_tol);
            if conv || div {
                break;
            }
        }
        Ok(finish(iters, res, conv, div, None))
    }
}
