//! BiCGStab, `huti_dbicgstabsolv` / `huti_zbicgstabsolv`
//! (huti_bicgstab.F90:279-566, 783-978).

use super::{converged, diverged, needs_rhsnorm, HUTI_EPSILON};
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
    let mut t1v = zero.clone();
    let mut v = zero.clone();
    let mut s = zero.clone();
    let mut t2v = zero.clone();
    let mut t = zero.clone();
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
    p.fill(T::zero());
    v.fill(T::zero());
    let mut oldrho = T::one();
    let mut omega = T::one();
    let mut alpha = T::zero();

    let mut residual = 0.0;
    let info = loop {
        let rho = o.dot(&rtld, &r);
        if rho.is_zero() {
            break HutiInfo::BiCGStabRho;
        }
        let beta = (rho * alpha).fdiv(oldrho * omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        o.pl(&p, &mut v);
        o.pr(&v, &mut t1v);
        o.mv(&t1v, &mut v);
        alpha = rho.fdiv(o.dot(&rtld, &v));
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        residual = o.norm(&s);
        if residual < HUTI_EPSILON {
            for i in 0..n {
                x[i] = x[i] + alpha * t1v[i];
            }
            break HutiInfo::Convergence;
        }
        o.pl(&s, &mut t);
        o.pr(&t, &mut t2v);
        o.mv(&t2v, &mut t);
        omega = o.dot(&t, &s).fdiv(o.dot(&t, &t));
        for i in 0..n {
            x[i] = x[i] + alpha * t1v[i] + omega * t2v[i];
        }
        for i in 0..n {
            r[i] = s[i] - omega * t[i];
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidualScaledByB => {
                o.mv(x, &mut t2v);
                for i in 0..n {
                    t1v[i] = t2v[i] - b[i];
                }
                o.norm(&t1v) / rhsnorm
            }
            HutiStop::PseudoResidual => o.norm(&r),
            HutiStop::PseudoResidualScaledByB => o.norm(&r) / rhsnorm,
            HutiStop::PseudoResidualScaledByPrecB => o.norm(&r) / precrhsnorm,
            HutiStop::XDiffNorm => {
                for i in 0..n {
                    t1v[i] = alpha * t1v[i] + omega * t2v[i];
                }
                o.norm(&t1v)
            }
            HutiStop::TrueResidual | HutiStop::UpperBound => {
                o.mv(x, &mut t2v);
                for i in 0..n {
                    t1v[i] = t2v[i] - b[i];
                }
                o.norm(&t1v)
            }
        };

        m.step(iter, residual);
        if converged(residual, c.tol) {
            break HutiInfo::Convergence;
        }
        // The real routine tests omega before divergence, the complex one after
        // (huti_bicgstab.F90:527-535 against 938-946).
        if !T::COMPLEX && omega.is_zero() {
            break HutiInfo::BiCGStabOmega;
        }
        if diverged(residual, c.max_tol) {
            break HutiInfo::Divergence;
        }
        if T::COMPLEX && omega.is_zero() {
            break HutiInfo::BiCGStabOmega;
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
