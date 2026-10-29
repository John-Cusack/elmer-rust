//! Conjugate gradient squared, `huti_dcgssolv` / `huti_zcgssolv`
//! (huti_cgs.F90:283-470, 679-866).

use super::{converged, diverged, needs_rhsnorm};
use crate::config::{HutiConfig, HutiStop};
use crate::huti::{HutiInfo, Monitor, Ops, Outcome};
use crate::scalar::HutiScalar;

pub(crate) fn solve<T: HutiScalar>(
    o: &Ops<T>,
    c: &HutiConfig,
    m: &mut Monitor,
    x: &mut [T],
    b: &[T],
) -> Outcome {
    let n = x.len();
    let zero = vec![T::zero(); n];
    let mut rtld = zero.clone();
    let mut p = zero.clone();
    let mut q = zero.clone();
    let mut u = zero.clone();
    let mut t1v = zero.clone();
    let mut t2v = zero.clone();
    let mut r = zero;

    let mut iter = 1;
    let mut rhsnorm = 0.0;
    let mut precrhsnorm = 0.0;
    if needs_rhsnorm(c.stop_crit) {
        rhsnorm = o.norm(b);
    }
    if c.stop_crit == HutiStop::PseudoResidualScaledByPrecB {
        o.pl(b, &mut p);
        precrhsnorm = o.norm(&p);
    }

    o.mv(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    rtld.copy_from_slice(&r);

    let mut oldrho = T::zero();
    let mut residual = 0.0;
    let info = loop {
        let rho = o.dot(&rtld, &r);
        if rho.is_zero() {
            break HutiInfo::CgsRho;
        }
        if iter == 1 {
            u.copy_from_slice(&r);
            p.copy_from_slice(&u);
        } else {
            let beta = rho.fdiv(oldrho);
            for i in 0..n {
                u[i] = r[i] + beta * q[i];
            }
            for i in 0..n {
                p[i] = u[i] + beta * q[i] + beta * beta * p[i];
            }
        }
        o.pl(&p, &mut t2v);
        o.pr(&t2v, &mut t1v);
        o.mv(&t1v, &mut t2v);
        let alpha = rho.fdiv(o.dot(&rtld, &t2v));
        for i in 0..n {
            q[i] = u[i] - alpha * t2v[i];
        }
        for i in 0..n {
            t2v[i] = u[i] + q[i];
        }
        o.pl(&t2v, &mut u);
        o.pr(&u, &mut t1v);
        for i in 0..n {
            x[i] = x[i] + alpha * t1v[i];
        }
        o.mv(&t1v, &mut t2v);
        for i in 0..n {
            r[i] = r[i] - alpha * t2v[i];
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidualScaledByB => {
                o.mv(x, &mut t1v);
                for i in 0..n {
                    t1v[i] = t1v[i] - b[i];
                }
                o.norm(&t1v) / rhsnorm
            }
            HutiStop::PseudoResidual => o.norm(&r),
            HutiStop::PseudoResidualScaledByB => o.norm(&r) / rhsnorm,
            HutiStop::PseudoResidualScaledByPrecB => o.norm(&r) / precrhsnorm,
            HutiStop::XDiffNorm => {
                for i in 0..n {
                    t1v[i] = alpha * t1v[i];
                }
                o.norm(&t1v)
            }
            HutiStop::TrueResidual | HutiStop::UpperBound => {
                o.mv(x, &mut t1v);
                for i in 0..n {
                    t1v[i] = t1v[i] - b[i];
                }
                o.norm(&t1v)
            }
        };

        m.step(iter, residual);
        if converged(residual, c.tol) {
            break HutiInfo::Convergence;
        }
        if diverged(residual, c.max_tol) {
            break HutiInfo::Divergence;
        }
        oldrho = rho;
        iter += 1;
        if iter > c.max_iter {
            break HutiInfo::MaxIter;
        }
    };
    m.finish(iter, residual);
    Outcome {
        info,
        iters: iter,
        residual,
    }
}
