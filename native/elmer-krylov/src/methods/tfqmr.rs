//! Transpose-free QMR, `huti_dtfqmrsolv` / `huti_ztfqmrsolv`
//! (huti_tfqmr.F90:455-803, 1171-1519). Each iteration takes two half steps,
//! each checked for convergence; the line is printed only after the second.
//! Like GMRES it maps `x` back through `pcondr` on exit.

use super::{converged, diverged};
use crate::config::{HutiConfig, HutiStop};
use crate::huti::{HutiInfo, Monitor, Ops, Outcome};
use crate::output;
use crate::scalar::HutiScalar;

/// `UPPERB_TOL_RATIO` (huti_tfqmr.F90:87).
const UPPERB_TOL_RATIO: f64 = 10.0;

struct Norms {
    rhs: f64,
    precrhs: f64,
}

pub(crate) fn solve<T: HutiScalar>(
    o: &Ops<T>,
    c: &HutiConfig,
    m: &mut Monitor,
    x: &mut [T],
    b: &[T],
) -> Outcome {
    let n = x.len();
    let zero = vec![T::zero(); n];
    let mut v = zero.clone();
    let mut y = zero.clone();
    let mut ynew = zero.clone();
    let mut rtld = zero.clone();
    let mut t1v = zero.clone();
    let mut t2v = zero.clone();
    let mut w = zero.clone();
    let mut d = zero.clone();
    let mut r = zero.clone();
    let mut trv = zero;

    let mut iter = 1;
    let mut nb = Norms {
        rhs: 0.0,
        precrhs: 0.0,
    };
    if matches!(
        c.stop_crit,
        HutiStop::TrueResidualScaledByB | HutiStop::PseudoResidualScaledByB | HutiStop::UpperBound
    ) {
        nb.rhs = o.norm(b);
    }
    if c.stop_crit == HutiStop::PseudoResidualScaledByPrecB {
        o.pl(b, &mut d);
        nb.precrhs = o.norm(&d);
    }

    o.pr(x, &mut d);
    o.mv(&d, &mut r);
    for i in 0..n {
        d[i] = b[i] - r[i];
    }
    o.pl(&d, &mut r);
    y.copy_from_slice(&r);
    w.copy_from_slice(&r);
    o.pr(&y, &mut v);
    o.mv(&v, &mut d);
    o.pl(&d, &mut v);
    t2v.copy_from_slice(&v);
    d.fill(T::zero());
    let mut tau = T::from_f64(o.norm(&r));
    let mut oldgamma = T::zero();
    let mut eta = T::zero();
    rtld.copy_from_slice(&r);
    let mut oldrho = o.dot(&rtld, &r);

    let mut residual = 0.0;
    let mut upperb = 0.0;
    let info = 'run: {
        if oldrho.is_zero() {
            break 'run HutiInfo::TfqmrRho;
        }
        loop {
            let alpha = oldrho.fdiv(o.dot(&rtld, &v));
            for i in 0..n {
                ynew[i] = y[i] - alpha * v[i];
            }
            for i in 0..n {
                w[i] = w[i] - alpha * t2v[i];
            }

            // First half step, from y.
            let gamma = T::from_f64(o.norm(&w)).fdiv(tau);
            let cc = T::one().fdiv((T::one() + gamma * gamma).fsqrt());
            tau = tau * gamma * cc;
            let f = (oldgamma * oldgamma * eta).fdiv(alpha);
            for i in 0..n {
                d[i] = y[i] + f * d[i];
            }
            eta = cc * cc * alpha;
            for i in 0..n {
                x[i] = x[i] + eta * d[i];
            }
            oldgamma = gamma;
            residual = stop(o, c, &nb, x, b, &d, eta, tau, iter, &mut r, &mut trv, &mut upperb);
            if converged(residual, c.tol) {
                break 'run HutiInfo::Convergence;
            }
            if diverged(residual, c.max_tol) {
                break 'run HutiInfo::Divergence;
            }

            // Second half step, from ynew.
            o.pr(&ynew, &mut t1v);
            o.mv(&t1v, &mut r);
            o.pl(&r, &mut t1v);
            for i in 0..n {
                w[i] = w[i] - alpha * t1v[i];
            }
            let gamma = T::from_f64(o.norm(&w)).fdiv(tau);
            let cc = T::one().fdiv((T::one() + gamma * gamma).fsqrt());
            tau = tau * gamma * cc;
            let f = (oldgamma * oldgamma * eta).fdiv(alpha);
            for i in 0..n {
                d[i] = ynew[i] + f * d[i];
            }
            eta = cc * cc * alpha;
            for i in 0..n {
                x[i] = x[i] + eta * d[i];
            }
            oldgamma = gamma;
            residual = stop(o, c, &nb, x, b, &d, eta, tau, iter, &mut r, &mut trv, &mut upperb);
            if converged(residual, c.tol) {
                break 'run HutiInfo::Convergence;
            }
            if diverged(residual, c.max_tol) {
                break 'run HutiInfo::Divergence;
            }

            if m.record(iter, residual) {
                m.emit(&line(c, iter, residual, upperb));
            }

            let rho = o.dot(&rtld, &w);
            let beta = rho.fdiv(oldrho);
            for i in 0..n {
                ynew[i] = w[i] + beta * ynew[i];
            }
            o.pr(&ynew, &mut t2v);
            o.mv(&t2v, &mut r);
            o.pl(&r, &mut t2v);
            for i in 0..n {
                v[i] = t2v[i] + beta * t1v[i] + beta * beta * v[i];
            }
            y.copy_from_slice(&ynew);
            oldrho = rho;
            iter += 1;
            if iter > c.max_iter {
                break 'run HutiInfo::MaxIter;
            }
        }
    };
    o.pr(x, &mut trv);
    x.copy_from_slice(&trv);
    if m.on() {
        m.emit(&line(c, iter, residual, upperb));
    }
    Outcome {
        info,
        iters: iter,
        residual,
    }
}

fn line(c: &HutiConfig, it: usize, res: f64, upperb: f64) -> String {
    let u = (c.stop_crit == HutiStop::UpperBound).then_some(upperb);
    output::tfqmr_line(it, res, u)
}

/// The `select case (HUTI_STOPC)` block, identical after both half steps.
#[allow(clippy::too_many_arguments)]
fn stop<T: HutiScalar>(
    o: &Ops<T>,
    c: &HutiConfig,
    nb: &Norms,
    x: &[T],
    b: &[T],
    d: &[T],
    eta: T,
    tau: T,
    iter: usize,
    r: &mut [T],
    trv: &mut [T],
    upperb: &mut f64,
) -> f64 {
    let n = x.len();
    let true_res = |r: &mut [T], trv: &mut [T]| {
        o.pr(x, trv);
        o.mv(trv, r);
        for i in 0..n {
            trv[i] = r[i] - b[i];
        }
        o.pl(trv, r);
        o.norm(r)
    };
    let pseudo_res = |r: &mut [T], trv: &mut [T]| {
        o.mv(x, r);
        for i in 0..n {
            r[i] = r[i] - b[i];
        }
        o.pl(r, trv);
        o.norm(trv)
    };
    match c.stop_crit {
        HutiStop::TrueResidual => true_res(r, trv),
        HutiStop::TrueResidualScaledByB => true_res(r, trv) / nb.rhs,
        HutiStop::PseudoResidual => pseudo_res(r, trv),
        HutiStop::PseudoResidualScaledByB => pseudo_res(r, trv) / nb.rhs,
        HutiStop::PseudoResidualScaledByPrecB => pseudo_res(r, trv) / nb.precrhs,
        HutiStop::XDiffNorm => {
            for i in 0..n {
                r[i] = eta * d[i];
            }
            o.norm(r)
        }
        HutiStop::UpperBound => {
            // `real( sqrt( 2.0 * iter_count ) * tau / rhsnorm )`: a single
            // precision square root, and REAL() of the quotient.
            let s = (2.0f32 * iter as f32).sqrt() as f64;
            *upperb = (T::from_f64(s) * tau).fdiv(T::from_f64(nb.rhs)).freal();
            if *upperb / c.tol < UPPERB_TOL_RATIO {
                true_res(r, trv) / nb.rhs
            } else {
                *upperb
            }
        }
    }
}
