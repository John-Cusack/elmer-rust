//! Preconditioned conjugate gradient, `huti_dcgsolv` / `huti_zcgsolv`
//! (huti_cg.F90:267-518, 720-899).

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
    let (mut z, mut p, mut q, mut r) = (zero.clone(), zero.clone(), zero.clone(), zero);

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

    // Fortran leaves both undefined until first assigned.
    let mut oldrho = T::zero();
    let mut residual = 0.0;
    let info = loop {
        o.pl(&r, &mut q);
        o.pr(&q, &mut z);
        let rho = o.dot(&r, &z);
        if rho.is_zero() {
            break HutiInfo::CgRho;
        }
        if iter == 1 {
            p.copy_from_slice(&z);
        } else {
            let beta = rho.fdiv(oldrho);
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
        }
        o.mv(&p, &mut q);
        let alpha = rho.fdiv(o.dot(&p, &q));
        for i in 0..n {
            x[i] = x[i] + alpha * p[i];
        }
        for i in 0..n {
            r[i] = r[i] - alpha * q[i];
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidualScaledByB => {
                o.mv(x, &mut z);
                for i in 0..n {
                    z[i] = z[i] - b[i];
                }
                o.norm(&z) / rhsnorm
            }
            HutiStop::PseudoResidual => o.norm(&r),
            HutiStop::PseudoResidualScaledByB => o.norm(&r) / rhsnorm,
            HutiStop::PseudoResidualScaledByPrecB => o.norm(&r) / precrhsnorm,
            HutiStop::XDiffNorm => {
                for i in 0..n {
                    z[i] = alpha * p[i];
                }
                o.norm(&z)
            }
            HutiStop::TrueResidual | HutiStop::UpperBound => {
                o.mv(x, &mut z);
                for i in 0..n {
                    z[i] = z[i] - b[i];
                }
                o.norm(&z)
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
