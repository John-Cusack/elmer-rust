//! Restarted GMRES(m), `huti_dgmressolv` / `huti_zgmressolv`
//! (huti_gmres.F90:390-822, 1143-1443) with `huti_[dz]lusolve`
//! (huti_aux.F90:221-289, 483-551). One HUTI "iteration" is one restart
//! cycle; the solver works in the right-preconditioned variable and maps `x`
//! back through `pcondr` on exit.

use super::{converged, diverged, needs_rhsnorm};
use crate::config::{HutiConfig, HutiStop};
use crate::huti::{HutiInfo, Monitor, Ops, Outcome};
use crate::output;
use crate::scalar::HutiScalar;

pub(crate) fn solve<T: HutiScalar>(
    o: &Ops<T>,
    c: &HutiConfig,
    mon: &mut Monitor,
    x: &mut [T],
    b: &[T],
) -> Outcome {
    let n = x.len();
    let m = c.restart;
    let zero = vec![T::zero(); n];
    let mut w = zero.clone();
    let mut r = zero.clone();
    let mut t1v = zero;
    // S and VTMP are work columns of length n in Fortran, indexed up to m+1;
    // they are sized to hold both so n < m+1 does not run off the end.
    let glen = n.max(m + 1);
    let mut s = vec![T::zero(); glen];
    let mut vtmp = vec![T::zero(); glen];
    let mut v = vec![vec![T::zero(); n]; m + 1];
    let mut h = vec![vec![T::zero(); m + 1]; m + 1]; // h[row][col], 0-based
    let mut cs = vec![T::zero(); m + 1];
    let mut sn = vec![T::zero(); m + 1];
    let mut y = vec![T::zero(); m + 1];
    let mut hlu = vec![T::zero(); m * m];

    let mut iter = 1;
    let bnrm = o.norm(b);
    let mut rhsnorm = 0.0;
    let mut precrhsnorm = 0.0;
    if needs_rhsnorm(c.stop_crit) {
        rhsnorm = bnrm;
    }
    if c.stop_crit == HutiStop::PseudoResidualScaledByPrecB {
        o.pl(b, &mut t1v);
        precrhsnorm = o.norm(&t1v);
    }

    o.pr(x, &mut t1v);
    o.mv(&t1v, &mut r);
    for i in 0..n {
        t1v[i] = b[i] - r[i];
    }
    o.pl(&t1v, &mut r);
    vtmp[0] = T::one();

    let mut residual = 0.0;
    let mut error = T::zero();
    let info = 'outer: loop {
        o.pr(x, &mut t1v);
        o.mv(&t1v, &mut r);
        for i in 0..n {
            t1v[i] = b[i] - r[i];
        }
        o.pl(&t1v, &mut r);
        let alpha = o.norm(&r);
        if alpha == 0.0 {
            break HutiInfo::GmresAlpha;
        }
        let ta = T::from_f64(alpha);
        for i in 0..n {
            v[0][i] = r[i].fdiv(ta);
        }
        for i in 0..glen {
            s[i] = ta * vtmp[i];
        }

        let mut inner_done = false;
        // `i` is the Fortran column index, 1..=m; arrays below are 0-based.
        for i in 1..=m {
            o.pr(&v[i - 1], &mut w);
            o.mv(&w, &mut t1v);
            o.pl(&t1v, &mut w);
            for k in 1..=i {
                let hk = o.dot(&w, &v[k - 1]);
                h[k - 1][i - 1] = hk;
                for j in 0..n {
                    w[j] = w[j] - hk * v[k - 1][j];
                }
            }
            let beta = o.norm(&w);
            if beta == 0.0 {
                break 'outer HutiInfo::GmresBeta;
            }
            h[i][i - 1] = T::from_f64(beta);
            for j in 0..n {
                v[i][j] = w[j].fdiv(h[i][i - 1]);
            }
            for k in 1..i {
                let temp = cs[k - 1] * h[k - 1][i - 1] + sn[k - 1] * h[k][i - 1];
                h[k][i - 1] = -(sn[k - 1] * h[k - 1][i - 1]) + cs[k - 1] * h[k][i - 1];
                h[k - 1][i - 1] = temp;
            }
            let (hii, hi1) = (h[i - 1][i - 1], h[i][i - 1]);
            if hi1.is_zero() {
                cs[i - 1] = T::one();
                sn[i - 1] = T::zero();
            } else if hi1.abs() > hii.abs() {
                let temp2 = hii.fdiv(hi1);
                sn[i - 1] = T::one().fdiv((T::one() + temp2 * temp2).fsqrt());
                cs[i - 1] = temp2 * sn[i - 1];
            } else {
                let temp2 = hi1.fdiv(hii);
                cs[i - 1] = T::one().fdiv((T::one() + temp2 * temp2).fsqrt());
                sn[i - 1] = temp2 * cs[i - 1];
            }
            let temp = cs[i - 1] * s[i - 1];
            s[i] = -(sn[i - 1] * s[i - 1]);
            s[i - 1] = temp;
            h[i - 1][i - 1] = cs[i - 1] * h[i - 1][i - 1] + sn[i - 1] * h[i][i - 1];
            h[i][i - 1] = T::zero();
            error = T::from_f64(s[i].abs() / bnrm);
            if error.freal() < c.tol {
                lusolve(i, &h, &mut hlu, &mut y, &s, mon);
                update(x, &v, &y, i);
                inner_done = true;
                break;
            }
        }
        if !inner_done {
            lusolve(m, &h, &mut hlu, &mut y, &s, mon);
            update(x, &v, &y, m);
        }

        residual = match c.stop_crit {
            HutiStop::TrueResidual | HutiStop::TrueResidualScaledByB => {
                o.pr(x, &mut t1v);
                o.mv(&t1v, &mut r);
                for i in 0..n {
                    t1v[i] = b[i] - r[i];
                }
                o.pl(&t1v, &mut r);
                let res = o.norm(&r);
                if c.stop_crit == HutiStop::TrueResidualScaledByB {
                    res / rhsnorm
                } else {
                    res
                }
            }
            HutiStop::PseudoResidual
            | HutiStop::PseudoResidualScaledByB
            | HutiStop::PseudoResidualScaledByPrecB => {
                o.mv(x, &mut r);
                for i in 0..n {
                    r[i] = r[i] - b[i];
                }
                o.pl(&r, &mut t1v);
                let res = o.norm(&t1v);
                match c.stop_crit {
                    HutiStop::PseudoResidualScaledByB => res / rhsnorm,
                    HutiStop::PseudoResidualScaledByPrecB => res / precrhsnorm,
                    _ => res,
                }
            }
            HutiStop::XDiffNorm => {
                t1v.fill(T::zero());
                update(&mut t1v, &v, &y, m);
                o.norm(&t1v)
            }
            HutiStop::UpperBound => {
                o.pr(x, &mut t1v);
                o.mv(&t1v, &mut r);
                for i in 0..n {
                    t1v[i] = r[i] - b[i];
                }
                o.pl(&t1v, &mut r);
                o.norm(&r)
            }
        };
        s[m] = T::from_f64(o.norm(&r));

        if mon.record(iter, residual) {
            mon.emit(&output::gmres_line(T::COMPLEX, iter, residual));
        }
        if converged(residual, c.tol) {
            break HutiInfo::Convergence;
        }
        if diverged(residual, c.max_tol) {
            break HutiInfo::Divergence;
        }
        iter += 1;
        if iter > c.max_iter {
            break HutiInfo::MaxIter;
        }
    };
    let _ = error;
    if mon.on() {
        mon.emit(&output::gmres_line(T::COMPLEX, iter, residual));
    }
    o.pr(x, &mut t1v);
    x.copy_from_slice(&t1v);
    Outcome {
        info,
        iters: iter,
        residual,
    }
}

/// `X = X + MATMUL(V(:, 1:k), Y(1:k))`: libgfortran forms the product in a
/// zeroed temporary, accumulating column by column, then adds it to `X`.
fn update<T: HutiScalar>(x: &mut [T], v: &[Vec<T>], y: &[T], k: usize) {
    for (j, xj) in x.iter_mut().enumerate() {
        let mut acc = T::zero();
        for (col, yc) in v.iter().zip(y).take(k) {
            acc = col[j].mul_acc(*yc, acc);
        }
        *xj = *xj + acc;
    }
}

/// `huti_[dz]lusolve`: copy the leading `k x k` block of `h` (column-major,
/// as HLU), Doolittle LU without pivoting (Saad, Algorithm 10.4), then the two
/// triangular solves into `y`.
fn lusolve<T: HutiScalar>(
    k: usize,
    h: &[Vec<T>],
    hlu: &mut [T],
    y: &mut [T],
    rhs: &[T],
    mon: &mut Monitor,
) {
    // lu[(row, col)] = hlu[row + col * k]
    for col in 0..k {
        for row in 0..k {
            hlu[row + col * k] = h[row][col];
        }
    }
    let at = |r: usize, c: usize| r + c * k;
    // `1.0e-16` is a default-real literal.
    let small = 1.0e-16_f32 as f64;
    for i in 1..k {
        for kk in 0..i {
            if hlu[at(kk, kk)].abs() < small {
                // List-directed in Fortran; the spacing is not reproduced.
                let p = hlu[at(kk, kk)].parts();
                mon.emit(&format!(
                    " (libhuti.a) GMRES: small pivot {:e} {:e}",
                    p[0], p[1]
                ));
            }
            hlu[at(i, kk)] = hlu[at(i, kk)].fdiv(hlu[at(kk, kk)]);
            for j in kk + 1..k {
                hlu[at(i, j)] = hlu[at(i, j)] - hlu[at(i, kk)] * hlu[at(kk, j)];
            }
        }
    }
    for i in 0..k {
        y[i] = rhs[i];
        for kk in 0..i {
            y[i] = y[i] - hlu[at(i, kk)] * y[kk];
        }
    }
    for i in (0..k).rev() {
        for kk in i + 1..k {
            y[i] = y[i] - hlu[at(i, kk)] * y[kk];
        }
        y[i] = y[i].fdiv(hlu[at(i, i)]);
    }
}
