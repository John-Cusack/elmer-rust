//! BiCGStab(2), `huti_dbicgstab_2solv` / `huti_zbicgstab_2solv`
//! (huti_bicgstab_2.F90:339-578, 841-1080). Works in the right-preconditioned
//! variable, so the true residual applies `pcondr` to `x` first; IterSolve
//! always runs it left-preconditioned, where that is the identity.

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
    let mut u = zero.clone();
    let mut t1v = zero.clone();
    let mut v = zero.clone();
    let mut s = zero.clone();
    let mut w = zero.clone();
    let mut t = zero.clone();
    let mut r = zero;

    let mut iter = 1;
    let mut rhsnorm = 0.0;
    let mut precrhsnorm = 0.0;
    if needs_rhsnorm(c.stop_crit) {
        rhsnorm = o.norm(b);
    }
    if c.stop_crit == HutiStop::PseudoResidualScaledByPrecB {
        o.pl(b, &mut t1v);
        precrhsnorm = o.norm(&t1v);
    }

    o.pr(x, &mut u);
    o.mv(&u, &mut r);
    for i in 0..n {
        u[i] = b[i] - r[i];
    }
    o.pl(&u, &mut r);
    rtld.copy_from_slice(&r);
    u.fill(T::zero());
    let mut oldrho = T::one();
    let mut omega2 = T::one();
    let mut alpha = T::zero();

    let mut residual = 0.0;
    let info = loop {
        oldrho = -(omega2 * oldrho);

        // Even Bi-CG step.
        let rho = o.dot(&rtld, &r);
        if rho.is_zero() {
            break HutiInfo::BiCGStab2Rho;
        }
        let beta = (rho * alpha).fdiv(oldrho);
        oldrho = rho;
        for i in 0..n {
            u[i] = r[i] - beta * u[i];
        }
        o.pr(&u, &mut v);
        o.mv(&v, &mut t1v);
        o.pl(&t1v, &mut v);
        alpha = oldrho.fdiv(o.dot(&rtld, &v));
        for i in 0..n {
            r[i] = r[i] - alpha * v[i];
        }
        o.pr(&r, &mut s);
        o.mv(&s, &mut t1v);
        o.pl(&t1v, &mut s);
        for i in 0..n {
            x[i] = x[i] + alpha * u[i];
        }

        // Odd Bi-CG step.
        let rho = o.dot(&rtld, &s);
        if rho.is_zero() {
            break HutiInfo::BiCGStab2Rho;
        }
        let beta = (rho * alpha).fdiv(oldrho);
        oldrho = rho;
        for i in 0..n {
            v[i] = s[i] - beta * v[i];
        }
        o.pr(&v, &mut w);
        o.mv(&w, &mut t1v);
        o.pl(&t1v, &mut w);
        alpha = oldrho.fdiv(o.dot(&rtld, &w));
        for i in 0..n {
            u[i] = r[i] - beta * u[i];
        }
        for i in 0..n {
            r[i] = r[i] - alpha * v[i];
        }
        for i in 0..n {
            s[i] = s[i] - alpha * w[i];
        }
        o.pr(&s, &mut t);
        o.mv(&t, &mut t1v);
        o.pl(&t1v, &mut t);

        // GCR(2) part.
        let mut omega1 = o.dot(&r, &s);
        let myy = o.dot(&s, &s);
        let delta = o.dot(&s, &t);
        let mut tau = o.dot(&t, &t);
        omega2 = o.dot(&r, &t);
        tau = tau - (delta * delta).fdiv(myy);
        omega2 = (omega2 - (delta * omega1).fdiv(myy)).fdiv(tau);
        omega1 = (omega1 - delta * omega2).fdiv(myy);
        for i in 0..n {
            x[i] = x[i] + omega1 * r[i] + omega2 * s[i] + alpha * u[i];
        }
        for i in 0..n {
            r[i] = r[i] - omega1 * s[i] - omega2 * t[i];
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidual | HutiStop::TrueResidualScaledByB | HutiStop::UpperBound => {
                o.pr(x, &mut s);
                o.mv(&s, &mut t1v);
                for i in 0..n {
                    t1v[i] = t1v[i] - b[i];
                }
                o.pl(&t1v, &mut s);
                let res = o.norm(&s);
                if c.stop_crit == HutiStop::TrueResidualScaledByB {
                    res / rhsnorm
                } else {
                    res
                }
            }
            HutiStop::PseudoResidual => o.norm(&r),
            HutiStop::PseudoResidualScaledByB => o.norm(&r) / rhsnorm,
            HutiStop::PseudoResidualScaledByPrecB => o.norm(&r) / precrhsnorm,
            HutiStop::XDiffNorm => {
                for i in 0..n {
                    t1v[i] = omega1 * r[i] + omega2 * s[i] + alpha * u[i];
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
        for i in 0..n {
            u[i] = u[i] - omega1 * v[i] - omega2 * w[i];
        }
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
