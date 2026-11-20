//! ARPACK's `zneupd` (the version Elmer bundles as `libarpack.so`): the
//! converged Ritz values and, with `rvec`, the Ritz vectors (`HOWMNY =
//! 'A'`) of the Arnoldi factorisation [`super::znaupd::znaupd`] left in
//! `v` and `workl`, in the Fortran's operation order.

// Index loops and argument lists mirror the Fortran's.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use num_complex::Complex64 as C;

use super::lapack::{dlapy2, EPS};
use super::naupd::NWhich;
use super::zblas::{cdiv, dznrm2, zcopy, zdotc, zdscal, zgeru, zscal, ztrmm_runn, ONE, ZERO};
use super::zlapack::{zgeqr2, zlahqr, zlaset, ztrevc_rs, ztrsen_nv, zunm2r};
use super::znaupd::Znaupd;

/// What [`zneupd`] returns: `D(1:nconv)`, `NCONV` and `INFO`.
#[derive(Clone, Debug, Default)]
pub struct Zneupd {
    pub d: Vec<C>,
    pub nconv: usize,
    pub info: i32,
}

/// `ZNEUPD(rvec, 'A', select, d, z, ldz, sigma, workev, bmat, n, which,
/// nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
/// rwork, info)` with `Z` the same array as `V` (Elmer passes `V` for
/// both), after [`super::znaupd::znaupd`] returned `au`.
pub fn zneupd(
    rvec: bool,
    sigma: C,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: usize,
    tol: f64,
    resid: &[C],
    ncv: usize,
    v: &mut [C],
    ldv: usize,
    au: &Znaupd,
    workd: &mut [C],
    workl: &mut [C],
) -> Zneupd {
    let mut out = Zneupd { d: vec![ZERO; nev.max(1)], ..Default::default() };
    let mode = au.iparam[6];
    let nconv = au.iparam[4].max(0) as usize;
    out.nconv = nconv;
    let eps23 = EPS.powf(2.0 / 3.0);
    let lworkl = 3 * ncv * ncv + 4 * ncv;
    let mut ierr = if nconv == 0 {
        -14
    } else if n == 0 {
        -1
    } else if nev == 0 {
        -2
    } else if ncv <= nev + 1 || ncv > n {
        -3
    } else if workl.len() < lworkl {
        -7
    } else {
        0
    };
    let shifti = match mode {
        1 | 2 => false,
        3 => true,
        _ => {
            ierr = -10;
            false
        }
    };
    if mode == 1 && bmat_g {
        ierr = -11;
    }
    if ierr != 0 {
        out.info = ierr;
        return out;
    }
    let _ = bmat_g;
    // 0-based offsets of the 1-based IPNTR pointers.
    let ih = au.ipntr[4] - 1;
    let ritz = au.ipntr[5] - 1;
    let iq = au.ipntr[6] - 1;
    let bounds = au.ipntr[7] - 1;
    let _ = iq;
    let (ldh, ldq) = (ncv, ncv);
    let iheig = bounds + ldh;
    let ihbds = iheig + ldh;
    let iuptri = ihbds + ldh;
    let invsub = iuptri + ldh * ncv;
    let irz = au.ipntr[13] - 1 + ncv * ncv;
    let ibd = irz + ncv;
    let mut workev = vec![ZERO; 3 * ncv];
    let rnorm = workl[ih + 2];
    workl[ih + 2] = ZERO;
    let mut select = vec![false; ncv];
    if rvec {
        let w0 = workl[ritz];
        let thres = match which {
            NWhich::LM | NWhich::SM => dlapy2(w0.re, w0.im),
            NWhich::LR | NWhich::SR => w0.re,
            NWhich::LI | NWhich::SI => w0.im,
        };
        let mut reord = false;
        for j in 0..ncv {
            let z = workl[irz + j];
            let passes = match which {
                NWhich::LM => dlapy2(z.re, z.im) >= thres,
                NWhich::SM => dlapy2(z.re, z.im) <= thres,
                NWhich::LR => z.re >= thres,
                NWhich::SR => z.re <= thres,
                NWhich::LI => z.im >= thres,
                NWhich::SI => z.im <= thres,
            };
            if passes {
                // ARPACK reads the Ritz value one place before (IRZ+J-1).
                let zp = workl[irz + j - 1];
                let rtemp = eps23.max(dlapy2(zp.re, zp.im));
                let b = workl[ibd + j];
                if dlapy2(b.re, b.im) <= tol * rtemp {
                    select[j] = true;
                }
            }
            if j + 1 > nconv {
                reord = select[j] || reord;
            }
        }
        // The Schur form of H and its Schur vectors.
        let hcopy: Vec<C> = workl[ih..ih + ldh * ncv].to_vec();
        zcopy(ldh * ncv, &hcopy, &mut workl[iuptri..]);
        {
            let (lo, hi) = workl.split_at_mut(invsub);
            zlaset(ncv, ncv, ZERO, ONE, hi, ldq);
            let (head, tri) = lo.split_at_mut(iuptri);
            let mut heig = vec![ZERO; ncv];
            let e = zlahqr(true, true, ncv, 1, ncv, tri, ldh, &mut heig, 1, ncv, hi, ldq);
            head[iheig..iheig + ncv].copy_from_slice(&heig);
            if e != 0 {
                out.info = -8;
                return out;
            }
        }
        for k in 0..ncv {
            workl[ihbds + k] = workl[invsub + ncv - 1 + k * ldq];
        }
        if reord {
            let (lo, hi) = workl.split_at_mut(invsub);
            let (head, tri) = lo.split_at_mut(iuptri);
            let mut heig = vec![ZERO; ncv];
            ztrsen_nv(&select, ncv, tri, ldh, hi, ldq, &mut heig);
            head[iheig..iheig + ncv].copy_from_slice(&heig);
        }
        for k in 0..ncv {
            workl[ihbds + k] = workl[invsub + ncv - 1 + k * ldq];
        }
        if !shifti {
            out.d[..nconv].copy_from_slice(&workl[iheig..iheig + nconv]);
        }
        // The Ritz vectors: the first NCONV columns of V Q.
        let tau = zgeqr2(ncv, nconv, &mut workl[invsub..], ldq);
        workev[..tau.len()].copy_from_slice(&tau);
        zunm2r(false, false, n, ncv, nconv, &mut workl[invsub..], ldq, &tau, v, ldv);
        // Z = V (the same array): ZLACPY is a no-op.
        for j in 1..=nconv {
            if workl[invsub + (j - 1) * ldq + j - 1].re < 0.0 {
                // ZSCAL(NCONV, -ONE, WORKL(IUPTRI+J-1), LDQ): row J of T.
                for c in 0..nconv {
                    let p = iuptri + j - 1 + c * ldq;
                    workl[p] = -ONE * workl[p];
                }
                zscal(nconv, -ONE, &mut workl[iuptri + (j - 1) * ldq..]);
            }
        }
        for j in 0..ncv {
            select[j] = j < nconv;
        }
        {
            let (lo, hi) = workl.split_at_mut(invsub);
            ztrevc_rs(&select, ncv, &mut lo[iuptri..], ldq, hi, ldq);
        }
        for j in 1..=nconv {
            let c0 = invsub + (j - 1) * ldq;
            let rtemp = 1.0 / dznrm2(ncv, &workl[c0..]);
            zdscal(ncv, rtemp, &mut workl[c0..]);
            workev[j - 1] = zdotc(j, &workl[ihbds..], &workl[c0..]);
        }
        let we: Vec<C> = workev[..nconv].to_vec();
        zcopy(nconv, &we, &mut workl[ihbds..]);
        let tri: Vec<C> = workl[invsub..invsub + ldq * ncv].to_vec();
        ztrmm_runn(n, nconv, ONE, &tri, ldq, v, ldv);
    } else {
        for k in 0..nconv {
            out.d[k] = workl[ritz + k];
            workl[iheig + k] = workl[ritz + k];
            workl[ihbds + k] = workl[bounds + k];
        }
    }
    if rvec {
        zscal(ncv, rnorm, &mut workl[ihbds..]);
    }
    if shifti {
        for k in 0..ncv {
            let temp = workl[iheig + k];
            workl[ihbds + k] = cdiv(cdiv(workl[ihbds + k], temp), temp);
        }
        for k in 0..nconv {
            out.d[k] = cdiv(ONE, workl[iheig + k]) + sigma;
        }
    }
    if rvec && shifti {
        for j in 1..=nconv {
            if workl[iheig + j - 1] != ZERO {
                workev[j - 1] = cdiv(workl[invsub + (j - 1) * ldq + ncv - 1], workl[iheig + j - 1]);
            }
        }
        zgeru(n, nconv, ONE, resid, &workev[..nconv], v, ldv);
    }
    let _ = workd;
    out
}
