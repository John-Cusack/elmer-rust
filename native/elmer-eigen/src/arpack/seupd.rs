//! `DSEUPD` (ARPACK 2.4): Ritz values and Ritz vectors from the final
//! Lanczos factorisation [`super::saupd::dsaupd`] left in `v` and `workl`,
//! back-transformed for the spectral transformation of `mode`.
//!
//! Only what Elmer's `ArpackEigenSolve` asks for is ported: `rvec = .true.`
//! with `howmny = 'A'`, and the Ritz vectors written over `v` (Elmer passes
//! `V` as both `v` and `z`).

use super::blas::{dcopy, dger, dscal};
use super::lapack::{dgeqr2, dorm2r, dsteqr_i};
use super::par::{eps23, pdnorm2, Par};
use super::saupd::{dsesrt, dsgets, dsortr, Saupd, Which};

#[derive(Copy, Clone, PartialEq, Eq)]
enum Kind {
    Regular,
    ShiftI,
    Buckle,
    Cayley,
}

/// `DSEUPD(.true., 'A', select, d, v, ldv, sigma, bmat, n, which, nev, tol,
/// resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)`. Returns
/// `D(1:ncv)` (the first `IPARAM(5)` are the Ritz values) and `info`; the
/// Ritz vectors are in the first `IPARAM(5)` columns of `v`.
#[allow(clippy::too_many_arguments)]
pub fn dseupd(
    sigma: f64,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Saupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> (Vec<f64>, i32) {
    seupd(None, sigma, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, au, workd, workl)
}

/// PARPACK's `PDSEUPD(comm, ...)` after [`super::saupd::pdsaupd`]: the B-norm
/// of the residual reduced over `comm`, and the wanted Ritz values picked
/// as PARPACK does, by sorting them with `_sgets` and taking the first
/// `IPARAM(5)` whose estimates pass (`info = -17` if fewer do), in place of
/// ARPACK 2.4's threshold test. Collective.
#[allow(clippy::too_many_arguments)]
pub fn pdseupd(
    comm: &dyn elmer_core::Communicator,
    sigma: f64,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Saupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> (Vec<f64>, i32) {
    seupd(Some(comm), sigma, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, au, workd, workl)
}

#[allow(clippy::too_many_arguments)]
fn seupd(
    par: Par<'_>,
    sigma: f64,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Saupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> (Vec<f64>, i32) {
    let mut d = vec![0.0; ncv];
    let mode = au.iparam[6];
    let nconv = au.iparam[4].max(0) as usize;
    if nconv == 0 {
        return (d, 0);
    }
    let mut ierr = 0;
    if n == 0 {
        ierr = -1;
    }
    if nev == 0 {
        ierr = -2;
    }
    if ncv <= nev || (ncv > n && par.is_none()) {
        ierr = -3;
    }
    if workl.len() < ncv * ncv + 8 * ncv {
        ierr = -7;
    }
    let kind = match mode {
        1 | 2 => Kind::Regular,
        3 => Kind::ShiftI,
        4 => Kind::Buckle,
        5 => Kind::Cayley,
        _ => {
            ierr = -10;
            Kind::Regular
        }
    };
    if mode == 1 && bmat_g {
        ierr = -11;
    }
    if nev == 1 && which == Which::BE {
        ierr = -12;
    }
    if ierr != 0 {
        return (d, ierr);
    }
    // 1-based offsets into workl, as the Fortran.
    let ih = au.ipntr[4];
    let ritz = au.ipntr[5];
    let bounds = au.ipntr[6];
    let ldh = ncv;
    let ldq = ncv;
    let ihd = bounds + ldh;
    let ihb = ihd + ldh;
    let iq = ihb + ldh;
    let iw = iq + ldh * ncv;
    let irz = au.ipntr[10] + ncv;
    let ibd = irz + ncv;
    macro_rules! w {
        ($i:expr) => {
            workl[$i - 1]
        };
    }
    let eps23 = eps23(par);
    let rnorm = w!(ih);
    let bnorm2 = if bmat_g { pdnorm2(par, n, workd) } else { rnorm };

    let mut select = vec![false; ncv];
    let mut reord = false;
    if par.is_some() {
        // pdseupd: sort the Ritz values (and their indices, in the bounds
        // slot) with _sgets; take the first nconv converged from the wanted
        // end.
        for j in 0..ncv {
            w!(bounds + j) = (j + 1) as f64;
        }
        {
            let (lo, hi) = workl.split_at_mut(irz - 1);
            dsgets(0, which, nev, ncv - nev, &mut hi[..ncv], &mut lo[bounds - 1..bounds - 1 + ncv], &mut []);
        }
        let mut numcnv = 0;
        for j in 1..=ncv {
            let temp1 = eps23.max(w!(irz + ncv - j).abs());
            let jj = w!(bounds + ncv - j) as usize;
            if numcnv < nconv && w!(ibd + jj - 1) <= tol * temp1 {
                select[jj - 1] = true;
                numcnv += 1;
                if jj > nev {
                    reord = true;
                }
            }
        }
        if numcnv != nconv {
            return (d, -17);
        }
    }
    // Select the wanted Ritz values (rvec = .true.).
    let (thres1, thres2) = if which == Which::BE {
        let ism = nev.max(nconv) / 2;
        (w!(ism), w!(ism + 1))
    } else {
        (w!(ritz), 0.0)
    };
    for j in (0..ncv).filter(|_| par.is_none()) {
        let r = w!(irz + j);
        let inside = match which {
            Which::LM => r.abs() >= thres1.abs(),
            Which::SM => r.abs() <= thres1.abs(),
            Which::LA => r >= thres1,
            Which::SA => r <= thres1,
            Which::BE => r <= thres1 || r >= thres2,
        };
        if inside {
            let tempbnd = eps23.max(r.abs());
            if w!(ibd + j) <= tol * tempbnd {
                select[j] = true;
            }
        }
        if j + 1 > nconv {
            reord = select[j] || reord;
        }
    }
    // The tridiagonal H: diagonal to ihd, subdiagonal to ihb; its Schur
    // vectors to iq.
    for k in 0..ncv - 1 {
        w!(ihb + k) = w!(ih + 1 + k);
    }
    for k in 0..ncv {
        w!(ihd + k) = w!(ih + ldh + k);
    }
    {
        let (_, rest) = workl.split_at_mut(ihd - 1);
        let (hd, rest) = rest.split_at_mut(ldh);
        let (hb, rest) = rest.split_at_mut(ldh);
        let (qm, rest) = rest.split_at_mut(ldh * ncv);
        let wk = &mut rest[..2 * ncv];
        if dsteqr_i(ncv, hd, hb, qm, ldq, wk) != 0 {
            return (d, -8);
        }
    }
    if reord && ncv > 1 {
        let (mut leftptr, mut rghtptr) = (1usize, ncv);
        loop {
            // 20
            if select[leftptr - 1] {
                leftptr += 1;
            } else if !select[rghtptr - 1] {
                rghtptr -= 1;
            } else {
                workl.swap(ihd + leftptr - 2, ihd + rghtptr - 2);
                let (l, r) = (iq - 1 + ncv * (leftptr - 1), iq - 1 + ncv * (rghtptr - 1));
                let tmp: Vec<f64> = workl[l..l + ncv].to_vec();
                workl.copy_within(r..r + ncv, l);
                workl[r..r + ncv].copy_from_slice(&tmp);
                // the Fortran passes the column through workl(iw).
                workl[iw - 1..iw - 1 + ncv].copy_from_slice(&tmp);
                leftptr += 1;
                rghtptr -= 1;
            }
            if leftptr >= rghtptr {
                break;
            }
        }
    }
    dcopy(nconv, &workl[ihd - 1..], &mut d);

    // Transform the Ritz values back and sort them with the vectors.
    if kind == Kind::Regular {
        let (_, rest) = workl.split_at_mut(iq - 1);
        dsesrt(Which::LA, true, nconv, &mut d, ncv, rest, ldq);
    } else {
        for k in 0..ncv {
            w!(iw + k) = w!(ihd + k);
        }
        for k in 0..ncv {
            let x = w!(ihd + k);
            w!(ihd + k) = match kind {
                Kind::ShiftI => 1.0 / x + sigma,
                Kind::Buckle => sigma * x / (x - 1.0),
                _ => sigma * (x + 1.0) / (x - 1.0),
            };
        }
        dcopy(nconv, &workl[ihd - 1..], &mut d);
        {
            let (lo, hi) = workl.split_at_mut(iw - 1);
            dsortr(Which::LA, true, nconv, &mut lo[ihd - 1..], hi);
        }
        let (_, rest) = workl.split_at_mut(iq - 1);
        dsesrt(Which::LA, true, nconv, &mut d, ncv, rest, ldq);
    }

    // Ritz vectors: v Q, through the QR factorisation of the selected Schur
    // vectors; then the last row of Q for the error bounds.
    {
        let (_, rest) = workl.split_at_mut(ihb - 1);
        let (hb, rest) = rest.split_at_mut(ldh);
        let (qm, rest) = rest.split_at_mut(ldh * ncv);
        let tau = &mut rest[ncv..2 * ncv];
        dgeqr2(ncv, nconv, qm, ldq, tau, hb);
        let (wd_lo, wd_hi) = workd.split_at_mut(n);
        let _ = wd_lo;
        dorm2r(false, false, n, ncv, nconv, qm, ldq, tau, v, ldv, wd_hi);
        hb[..ncv - 1].iter_mut().for_each(|x| *x = 0.0);
        hb[ncv - 1] = 1.0;
        let mut temp = [0.0f64; 1];
        dorm2r(true, true, ncv, 1, nconv, qm, ldq, tau, hb, ncv, &mut temp);
    }
    if kind == Kind::Regular {
        for j in 0..ncv {
            w!(ihb + j) = rnorm * w!(ihb + j).abs();
        }
    } else {
        dscal(ncv, bnorm2, &mut workl[ihb - 1..]);
        for k in 0..ncv {
            let hbk = w!(ihb + k);
            let wk = w!(iw + k);
            w!(ihb + k) = match kind {
                Kind::ShiftI => hbk.abs() / (wk * wk),
                Kind::Buckle => sigma * hbk.abs() / ((wk - 1.0) * (wk - 1.0)),
                _ => (hbk / wk * (wk - 1.0)).abs(),
            };
        }
    }
    match kind {
        Kind::ShiftI | Kind::Cayley => {
            for k in 0..nconv {
                w!(iw + k) = w!(iq + k * ldq + ncv - 1) / w!(iw + k);
            }
        }
        Kind::Buckle => {
            for k in 0..nconv {
                w!(iw + k) = w!(iq + k * ldq + ncv - 1) / (w!(iw + k) - 1.0);
            }
        }
        Kind::Regular => {}
    }
    if kind != Kind::Regular {
        dger(n, nconv, 1.0, resid, &workl[iw - 1..], v, ldv);
    }
    (d, 0)
}
