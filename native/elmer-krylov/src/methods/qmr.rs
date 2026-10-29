//! Quasi-minimal residual, `huti_dqmrsolv` / `huti_zqmrsolv`
//! (huti_qmr.F90:363-605, 869-1111). The `HUTI_MAT_TRPSED` calls map to the
//! operator's and preconditioner's `apply_transpose`.

use elmer_core::LinalgError;

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
) -> Result<Outcome, LinalgError> {
    let n = x.len();
    let zero = vec![T::zero(); n];
    let mut v = zero.clone();
    let mut vtld = zero.clone();
    let mut y = zero.clone();
    let mut ytld = zero.clone();
    let mut w = zero.clone();
    let mut wtld = zero.clone();
    let mut z = zero.clone();
    let mut ztld = zero.clone();
    let mut p = zero.clone();
    let mut ptld = zero.clone();
    let mut q = zero.clone();
    let mut d = zero.clone();
    let mut s = zero.clone();
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
    vtld.copy_from_slice(&r);
    o.pl(&vtld, &mut y);
    wtld.copy_from_slice(&r);
    o.pr(&wtld, &mut z);
    let mut rho = T::from_f64(o.norm(&y));
    let mut psi = T::from_f64(o.norm(&z));
    let mut oldgamma = T::one();
    let mut eta = -T::one();
    // Undefined in Fortran until the first iteration assigns them.
    let mut epsilon = T::zero();
    let mut oldtheta = T::zero();

    let mut residual = 0.0;
    let info = loop {
        if rho.is_zero() || psi.is_zero() {
            break HutiInfo::QmrRhoPsi;
        }
        for i in 0..n {
            v[i] = vtld[i].fdiv(rho);
        }
        for i in 0..n {
            y[i] = y[i].fdiv(rho);
        }
        for i in 0..n {
            w[i] = wtld[i].fdiv(psi);
        }
        for i in 0..n {
            z[i] = z[i].fdiv(psi);
        }
        let delta = o.dot(&z, &y);
        if delta.is_zero() {
            break HutiInfo::QmrDelta;
        }
        o.pr(&y, &mut ytld);
        o.pl_t(&z, &mut ztld)?;
        if iter == 1 {
            p.copy_from_slice(&ytld);
            q.copy_from_slice(&ztld);
        } else {
            let fp = (psi * delta).fdiv(epsilon);
            for i in 0..n {
                p[i] = ytld[i] - fp * p[i];
            }
            let fq = (rho * delta).fdiv(epsilon);
            for i in 0..n {
                q[i] = ztld[i] - fq * q[i];
            }
        }
        o.mv(&p, &mut ptld);
        epsilon = o.dot(&q, &ptld);
        if epsilon.is_zero() {
            break HutiInfo::QmrEpsilon;
        }
        let beta = epsilon.fdiv(delta);
        if beta.is_zero() {
            break HutiInfo::QmrBeta;
        }
        for i in 0..n {
            vtld[i] = ptld[i] - beta * v[i];
        }
        o.pl(&vtld, &mut y);
        let rhonext = T::from_f64(o.norm(&y));
        o.mv_t(&q, &mut wtld)?;
        for i in 0..n {
            wtld[i] = wtld[i] - beta * w[i];
        }
        o.pr_t(&wtld, &mut z)?;
        psi = T::from_f64(o.norm(&z));
        let theta = rhonext.fdiv(oldgamma * T::from_f64(beta.abs()));
        let gamma = T::one().fdiv((T::one() + theta * theta).fsqrt());
        if gamma.is_zero() {
            break HutiInfo::QmrGamma;
        }
        // `-1 * (...) / (...)` parses as `-((1 * ...) / ...)`.
        eta = -((eta * rho * gamma * gamma).fdiv(beta * oldgamma * oldgamma));
        if iter == 1 {
            for i in 0..n {
                d[i] = eta * p[i];
            }
            for i in 0..n {
                s[i] = eta * ptld[i];
            }
        } else {
            for i in 0..n {
                d[i] = eta * p[i] + (oldtheta * gamma) * (oldtheta * gamma) * d[i];
            }
            for i in 0..n {
                s[i] = eta * ptld[i] + (oldtheta * gamma) * (oldtheta * gamma) * s[i];
            }
        }
        for i in 0..n {
            x[i] = x[i] + d[i];
        }
        for i in 0..n {
            r[i] = r[i] - s[i];
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidualScaledByB => {
                o.mv(x, &mut ytld);
                for i in 0..n {
                    ytld[i] = ytld[i] - b[i];
                }
                o.norm(&ytld) / rhsnorm
            }
            HutiStop::PseudoResidual => o.norm(&r),
            HutiStop::PseudoResidualScaledByB => o.norm(&r) / rhsnorm,
            HutiStop::PseudoResidualScaledByPrecB => o.norm(&r) / precrhsnorm,
            HutiStop::XDiffNorm => o.norm(&d),
            HutiStop::TrueResidual | HutiStop::UpperBound => {
                o.mv(x, &mut ytld);
                for i in 0..n {
                    ytld[i] = ytld[i] - b[i];
                }
                o.norm(&ytld)
            }
        };

        m.step(iter, residual);
        if converged(residual, c.tol) {
            break HutiInfo::Convergence;
        }
        if diverged(residual, c.max_tol) {
            break HutiInfo::Divergence;
        }
        rho = rhonext;
        oldgamma = gamma;
        oldtheta = theta;
        iter += 1;
        if iter > c.max_iter {
            break HutiInfo::MaxIter;
        }
    };
    m.finish(iter, residual);
    Ok(Outcome {
        info,
        iters: iter,
        residual,
    })
}
