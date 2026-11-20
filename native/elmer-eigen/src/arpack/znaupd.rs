//! ARPACK's complex implicitly restarted Arnoldi driver (`znaupd`, the
//! version Elmer bundles and links as `libarpack.so`), with its helpers
//! `znaup2`, `znaitr`, `znapps`, `zneigh`, `zngets`, `zsortc` and `zgetv0`,
//! in the Fortran's operation order.
//!
//! Reverse communication is a callback, as in [`super::naupd`]:
//! `op(ido, ipntr, workd)` performs the request `ido` on `workd` (0-based
//! offsets in `ipntr`) and the routines continue where the Fortran
//! re-enters. The `SAVE`d state of the complex routines (the seed of
//! `zgetv0` and the `smlnum` of `znaitr` and `znapps`) is [`ZArpack`],
//! separate from the real routines' [`super::saupd::Arpack`].

// Index loops and argument lists mirror the Fortran's.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use num_complex::Complex64 as C;

use super::lapack::{dlapy2, EPS, PREC, SAFMIN};
use super::naupd::NWhich;
use super::saupd::{c0717, Stats};
use super::zblas::{cabs1, dznrm2, zaxpy, zcopy, zdotc, zdscal, zgemv_c, zgemv_n, zscal, ONE, ZERO};
use super::zlapack::{zlacpy, zlahqr, zlanhs_one, zlarnv2, zlartg, zlascl_g, zlaset, ztrevc_rb};

/// The complex reverse-communication callback: `op(ido, ipntr, workd)`.
pub type ZOp<'a> = dyn FnMut(i32, [usize; 3], &mut [C]) -> Result<(), String> + 'a;

/// The `SAVE`d state of ARPACK's complex routines: `zgetv0`'s seed
/// (starting at 1, 3, 5, 7) and the `smlnum` that `znaitr` and `znapps`
/// each fix on their first call.
#[derive(Clone, Debug)]
pub struct ZArpack {
    pub iseed: [i32; 4],
    pub naitr_smlnum: Option<f64>,
    pub napps_smlnum: Option<f64>,
}

impl Default for ZArpack {
    fn default() -> Self {
        ZArpack { iseed: [1, 3, 5, 7], naitr_smlnum: None, napps_smlnum: None }
    }
}

#[inline]
fn ix(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

/// `sqrt(dlapy2(re, im))` of a complex inner product: ARPACK's B-norm.
fn bnorm(c: C) -> f64 {
    dlapy2(c.re, c.im).sqrt()
}

/// `ZSORTC(which, apply, n, x, y)`: Shell sort of `x` (and `y` alongside
/// when `apply`) into the order that puts the wanted values last.
pub fn zsortc(which: NWhich, apply: bool, n: usize, x: &mut [C], y: &mut [C]) {
    let key_gt = |a: C, b: C| -> bool {
        match which {
            NWhich::LM => dlapy2(a.re, a.im) > dlapy2(b.re, b.im),
            NWhich::SM => dlapy2(a.re, a.im) < dlapy2(b.re, b.im),
            NWhich::LR => a.re > b.re,
            NWhich::SR => a.re < b.re,
            NWhich::LI => a.im > b.im,
            NWhich::SI => a.im < b.im,
        }
    };
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let ju = j as usize;
                if key_gt(x[ju], x[ju + igap]) {
                    x.swap(ju, ju + igap);
                    if apply {
                        y.swap(ju, ju + igap);
                    }
                } else {
                    break;
                }
                j -= igap as isize;
            }
        }
        igap /= 2;
    }
}

/// `ZNGETS(ishift, which, kev, np, ritz, bounds)`: sort the Ritz values so
/// the `np` unwanted ones come first, those with the largest error bounds
/// first among them (exact shifts).
pub fn zngets(ishift: i32, which: NWhich, kev: usize, np: usize, ritz: &mut [C], bounds: &mut [C]) {
    zsortc(which, true, kev + np, ritz, bounds);
    if ishift == 1 {
        zsortc(NWhich::SM, true, np, bounds, ritz);
    }
}

/// `ZGETV0(ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr,
/// workd, ierr)`: a (random) starting vector in the range of `OP`,
/// orthogonal to the first `j - 1` columns of `v`. Returns `(rnorm,
/// ierr)`.
pub fn zgetv0(
    op: &mut ZOp<'_>,
    st: &mut Stats,
    arp: &mut ZArpack,
    bmat_g: bool,
    initv: bool,
    n: usize,
    j: usize,
    v: &[C],
    ldv: usize,
    resid: &mut [C],
    workd: &mut [C],
) -> Result<(f64, i32), String> {
    let mut ierr = 0;
    let mut iter = 0;
    if !initv {
        zlarnv2(&mut arp.iseed, n, resid);
    }
    if bmat_g {
        st.nopx += 1;
        zcopy(n, resid, workd);
        op(-1, [0, n, 0], workd)?;
        st.nbx += 1;
        let hi: Vec<C> = workd[n..2 * n].to_vec();
        zcopy(n, &hi, resid);
        op(2, [n, 0, 0], workd)?;
    } else {
        zcopy(n, resid, workd);
    }
    let mut rnorm0 = if bmat_g { bnorm(zdotc(n, resid, workd)) } else { dznrm2(n, resid) };
    let mut rnorm = rnorm0;
    if j == 1 {
        return Ok((rnorm, ierr));
    }
    loop {
        // 30
        {
            let (lo, hi) = workd.split_at_mut(n);
            zgemv_c(n, j - 1, ONE, v, ldv, lo, ZERO, hi);
            zgemv_n(n, j - 1, -ONE, v, ldv, hi, ONE, resid);
        }
        if bmat_g {
            st.nbx += 1;
            zcopy(n, resid, &mut workd[n..]);
            op(2, [n, 0, 0], workd)?;
        } else {
            zcopy(n, resid, workd);
        }
        // 40
        rnorm = if bmat_g { bnorm(zdotc(n, resid, workd)) } else { dznrm2(n, resid) };
        if rnorm > c0717() * rnorm0 {
            break;
        }
        iter += 1;
        if iter <= 1 {
            rnorm0 = rnorm;
            continue;
        }
        resid[..n].iter_mut().for_each(|x| *x = ZERO);
        rnorm = 0.0;
        ierr = -1;
        break;
    }
    Ok((rnorm, ierr))
}

fn smlnum_once(cell: &mut Option<f64>, n: usize) -> f64 {
    *cell.get_or_insert(SAFMIN * (n as f64 / PREC))
}

/// `ZNAITR(ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, ipntr,
/// workd, info)`: extend a `k`-step Arnoldi factorisation by `np` steps.
/// Returns `info`.
fn znaitr(
    op: &mut ZOp<'_>,
    st: &mut Stats,
    arp: &mut ZArpack,
    bmat_g: bool,
    n: usize,
    k: usize,
    np: usize,
    resid: &mut [C],
    rnorm: &mut f64,
    v: &mut [C],
    ldv: usize,
    h: &mut [C],
    ldh: usize,
    workd: &mut [C],
) -> Result<usize, String> {
    let smlnum = smlnum_once(&mut arp.naitr_smlnum, n);
    let unfl = SAFMIN;
    let ulp = PREC;
    let (ipj, irj, ivj) = (0usize, n, 2 * n);
    let mut j = k + 1;
    loop {
        let mut betaj = *rnorm;
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(*rnorm > 0.0) {
            betaj = 0.0;
            st.nrstrt += 1;
            let mut itry = 1;
            loop {
                let (rn, ierr) = zgetv0(op, st, arp, bmat_g, false, n, j, v, ldv, resid, workd)?;
                *rnorm = rn;
                if ierr < 0 {
                    itry += 1;
                    if itry <= 3 {
                        continue;
                    }
                    return Ok(j - 1);
                }
                break;
            }
        }
        // 40
        let vj = (j - 1) * ldv;
        zcopy(n, resid, &mut v[vj..]);
        if *rnorm >= unfl {
            let temp1 = 1.0 / *rnorm;
            zdscal(n, temp1, &mut v[vj..]);
            zdscal(n, temp1, &mut workd[ipj..]);
        } else {
            zlascl_g(*rnorm, 1.0, n, 1, &mut v[vj..], n);
            zlascl_g(*rnorm, 1.0, n, 1, &mut workd[ipj..], n);
        }
        // STEP 3: r_j = OP v_j.
        st.nopx += 1;
        let vcol: Vec<C> = v[vj..vj + n].to_vec();
        zcopy(n, &vcol, &mut workd[ivj..]);
        op(1, [ivj, irj, ipj], workd)?;
        // 50
        let rj: Vec<C> = workd[irj..irj + n].to_vec();
        zcopy(n, &rj, resid);
        if bmat_g {
            st.nbx += 1;
            op(2, [irj, ipj, 0], workd)?;
        } else {
            zcopy(n, resid, &mut workd[ipj..]);
        }
        // 60
        let wnorm = if bmat_g { bnorm(zdotc(n, resid, &workd[ipj..])) } else { dznrm2(n, resid) };
        let hj = (j - 1) * ldh;
        zgemv_c(n, j, ONE, v, ldv, &workd[ipj..], ZERO, &mut h[hj..]);
        {
            let col: Vec<C> = h[hj..hj + j].to_vec();
            zgemv_n(n, j, -ONE, v, ldv, &col, ONE, resid);
        }
        if j > 1 {
            h[ix(j, j - 1, ldh)] = C::new(betaj, 0.0);
        }
        if bmat_g {
            st.nbx += 1;
            zcopy(n, resid, &mut workd[irj..]);
            op(2, [irj, ipj, 0], workd)?;
        } else {
            zcopy(n, resid, &mut workd[ipj..]);
        }
        // 70
        *rnorm = if bmat_g { bnorm(zdotc(n, resid, &workd[ipj..])) } else { dznrm2(n, resid) };
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(*rnorm > c0717() * wnorm) {
            let mut iter = 0;
            st.nrorth += 1;
            loop {
                // 80
                {
                    let (p, r) = workd.split_at_mut(irj);
                    zgemv_c(n, j, ONE, v, ldv, p, ZERO, r);
                }
                {
                    let r: Vec<C> = workd[irj..irj + j].to_vec();
                    zgemv_n(n, j, -ONE, v, ldv, &r, ONE, resid);
                    zaxpy(j, ONE, &r, &mut h[hj..]);
                }
                if bmat_g {
                    st.nbx += 1;
                    zcopy(n, resid, &mut workd[irj..]);
                    op(2, [irj, ipj, 0], workd)?;
                } else {
                    zcopy(n, resid, &mut workd[ipj..]);
                }
                // 90
                let rnorm1 = if bmat_g { bnorm(zdotc(n, resid, &workd[ipj..])) } else { dznrm2(n, resid) };
                if rnorm1 > c0717() * *rnorm {
                    *rnorm = rnorm1;
                    break;
                }
                st.nitref += 1;
                *rnorm = rnorm1;
                iter += 1;
                if iter <= 1 {
                    continue;
                }
                resid[..n].iter_mut().for_each(|x| *x = ZERO);
                *rnorm = 0.0;
                break;
            }
        }
        // 100
        j += 1;
        if j > k + np {
            for i in k.max(1)..k + np {
                let (hii, hi1) = (h[ix(i, i, ldh)], h[ix(i + 1, i + 1, ldh)]);
                let mut tst1 = dlapy2(hii.re, hii.im) + dlapy2(hi1.re, hi1.im);
                if tst1 == 0.0 {
                    tst1 = zlanhs_one(k + np, h, ldh);
                }
                let s = h[ix(i + 1, i, ldh)];
                if dlapy2(s.re, s.im) <= (ulp * tst1).max(smlnum) {
                    h[ix(i + 1, i, ldh)] = ZERO;
                }
            }
            return Ok(0);
        }
    }
}

/// `ZNAPPS(n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workl,
/// workd)`: apply the `np` shifts implicitly, keeping `kev` vectors.
fn znapps(
    arp: &mut ZArpack,
    n: usize,
    kev: usize,
    np: usize,
    shift: &[C],
    v: &mut [C],
    ldv: usize,
    h: &mut [C],
    ldh: usize,
    resid: &mut [C],
    q: &mut [C],
    ldq: usize,
    workd: &mut [C],
) {
    let smlnum = smlnum_once(&mut arp.napps_smlnum, n);
    let ulp = PREC;
    let kplusp = kev + np;
    zlaset(kplusp, kplusp, ZERO, ONE, q, ldq);
    if np == 0 {
        return;
    }
    for jj in 1..=np {
        let sigma = shift[jj - 1];
        let mut istart = 1;
        loop {
            // 20: look for a small subdiagonal element to split at.
            let mut iend = kplusp;
            for i in istart..kplusp {
                let mut tst1 = cabs1(h[ix(i, i, ldh)]) + cabs1(h[ix(i + 1, i + 1, ldh)]);
                if tst1 == 0.0 {
                    tst1 = zlanhs_one(kplusp - jj + 1, h, ldh);
                }
                if h[ix(i + 1, i, ldh)].re.abs() <= (ulp * tst1).max(smlnum) {
                    iend = i;
                    h[ix(i + 1, i, ldh)] = ZERO;
                    break;
                }
            }
            // 40
            if !(istart == iend || istart > kev) {
                let h11 = h[ix(istart, istart, ldh)];
                let h21 = h[ix(istart + 1, istart, ldh)];
                let mut f = h11 - sigma;
                let mut g = h21;
                for i in istart..iend {
                    let (c, s, r) = zlartg(f, g);
                    if i > istart {
                        h[ix(i, i - 1, ldh)] = r;
                        h[ix(i + 1, i - 1, ldh)] = ZERO;
                    }
                    for j in i..=kplusp {
                        let (a, b) = (h[ix(i, j, ldh)], h[ix(i + 1, j, ldh)]);
                        let t = c * a + s * b;
                        h[ix(i + 1, j, ldh)] = -(s.conj() * a) + c * b;
                        h[ix(i, j, ldh)] = t;
                    }
                    let jend = (i + 2).min(iend);
                    for j in 1..=jend {
                        let (a, b) = (h[ix(j, i, ldh)], h[ix(j, i + 1, ldh)]);
                        let t = c * a + s.conj() * b;
                        h[ix(j, i + 1, ldh)] = -(s * a) + c * b;
                        h[ix(j, i, ldh)] = t;
                    }
                    // `DO 70 J = 1, MIN(J+JJ, KPLUSP)` reads J as loop 60 left it.
                    let j60 = jend + 1;
                    for j in 1..=(j60 + jj).min(kplusp) {
                        let (a, b) = (q[ix(j, i, ldq)], q[ix(j, i + 1, ldq)]);
                        let t = c * a + s.conj() * b;
                        q[ix(j, i + 1, ldq)] = -(s * a) + c * b;
                        q[ix(j, i, ldq)] = t;
                    }
                    if i + 1 < iend {
                        f = h[ix(i + 1, i, ldh)];
                        g = h[ix(i + 2, i, ldh)];
                    }
                }
            }
            // 100
            istart = iend + 1;
            if iend >= kplusp {
                break;
            }
        }
    }
    // Make the subdiagonals of the leading KEV block real and
    // nonnegative.
    for j in 1..=kev {
        let hj = h[ix(j + 1, j, ldh)];
        if hj.re < 0.0 || hj.im != 0.0 {
            let t = hj / dlapy2(hj.re, hj.im);
            // ZSCAL(KPLUSP-J+1, CONJG(T), H(J+1,J), LDH).
            let ct = t.conj();
            if ct != ONE {
                for c in j..=kplusp {
                    let p = ix(j + 1, c, ldh);
                    h[p] = ct * h[p];
                }
            }
            zscal((j + 2).min(kplusp), t, &mut h[ix(1, j + 1, ldh)..]);
            zscal((j + np + 1).min(kplusp), t, &mut q[ix(1, j + 1, ldq)..]);
            let p = ix(j + 1, j, ldh);
            h[p] = C::new(h[p].re, 0.0);
        }
    }
    for i in 1..=kev {
        let mut tst1 = cabs1(h[ix(i, i, ldh)]) + cabs1(h[ix(i + 1, i + 1, ldh)]);
        if tst1 == 0.0 {
            tst1 = zlanhs_one(kev, h, ldh);
        }
        if h[ix(i + 1, i, ldh)].re <= (ulp * tst1).max(smlnum) {
            h[ix(i + 1, i, ldh)] = ZERO;
        }
    }
    // Compute the KEV+1 st column of (V*Q) and temporarily store it.
    let hk = h[ix(kev + 1, kev, ldh)];
    if hk.re > 0.0 {
        let (w0, w1) = workd.split_at_mut(n);
        let _ = w0;
        zgemv_n(n, kplusp, ONE, v, ldv, &q[ix(1, kev + 1, ldq)..], ZERO, w1);
    }
    // Compute column 1 to KEV of (V*Q) in backward order.
    for i in 1..=kev {
        zgemv_n(n, kplusp - i + 1, ONE, v, ldv, &q[ix(1, kev - i + 1, ldq)..], ZERO, workd);
        let w: Vec<C> = workd[..n].to_vec();
        zcopy(n, &w, &mut v[ix(1, kplusp - i + 1, ldv)..]);
    }
    // Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).
    let src: Vec<C> = v[ix(1, kplusp - kev + 1, ldv)..].to_vec();
    zlacpy(n, kev, &src, ldv, v, ldv);
    if hk.re > 0.0 {
        let w: Vec<C> = workd[n..2 * n].to_vec();
        zcopy(n, &w, &mut v[ix(1, kev + 1, ldv)..]);
    }
    // Update the residual vector: r <- sigmak*r + betak*v(:,kev+1).
    zscal(n, q[ix(kplusp, kev, ldq)], resid);
    if hk.re > 0.0 {
        let vk: Vec<C> = v[ix(1, kev + 1, ldv)..ix(1, kev + 1, ldv) + n].to_vec();
        zaxpy(n, hk, &vk, resid);
    }
}

/// `ZNEIGH(rnorm, n, h, ldh, ritz, bounds, q, ldq, workl, rwork, ierr)`:
/// the eigenvalues of `H` and their Ritz estimates. `workl` (at least `n *
/// n + 2 n`) gets the Schur form; `q` the eigenvectors. Returns `ierr`.
pub fn zneigh(rnorm: f64, n: usize, h: &[C], ldh: usize, ritz: &mut [C], bounds: &mut [C], q: &mut [C], ldq: usize, workl: &mut [C]) -> usize {
    zlacpy(n, n, h, ldh, workl, n);
    zlaset(n, n, ZERO, ONE, q, ldq);
    let ierr = zlahqr(true, true, n, 1, n, workl, ldh, ritz, 1, n, q, ldq);
    if ierr != 0 {
        return ierr;
    }
    for j in 0..n {
        bounds[j] = q[ix(n - 1, j + 1, ldq)];
    }
    ztrevc_rb(n, workl, n, q, ldq);
    for j in 1..=n {
        let temp = dznrm2(n, &q[ix(1, j, ldq)..]);
        zdscal(n, 1.0 / temp, &mut q[ix(1, j, ldq)..]);
    }
    // ZCOPY(N, Q(N,1), N, BOUNDS, 1): the last row (stride N = LDQ).
    for j in 0..n {
        bounds[j] = q[(n - 1) + j * n];
    }
    zdscal(n, rnorm, bounds);
    0
}

/// The parameters and results of [`znaupd`] (`iparam`, `ipntr` and `info`).
#[derive(Clone, Debug, Default)]
pub struct Znaupd {
    /// `IPARAM(1..11)` (index 0 is `IPARAM(1)`).
    pub iparam: [i32; 11],
    /// `IPNTR(1..14)`, 1-based offsets into `workl` (as the Fortran).
    pub ipntr: [usize; 14],
    pub info: i32,
    pub stats: Stats,
}

/// `ZNAUPD(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
/// ipntr, workd, workl, lworkl, rwork, info)` driven to completion.
/// `iparam` carries `ishift`, `mxiter` and `mode` in (`IPARAM(1)`, `(3)`,
/// `(7)`); `info_in` is 0 for a random starting vector, nonzero to use
/// `resid`. `workd` has `3 n`, `workl` at least `3 ncv**2 + 5 ncv` entries.
pub fn znaupd(
    op: &mut ZOp<'_>,
    arp: &mut ZArpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: usize,
    tol: &mut f64,
    resid: &mut [C],
    ncv: usize,
    v: &mut [C],
    ldv: usize,
    iparam: [i32; 11],
    workd: &mut [C],
    workl: &mut [C],
    info_in: i32,
) -> Result<Znaupd, String> {
    let mut out = Znaupd { iparam, ..Default::default() };
    let mut ishift = iparam[0];
    let mut mxiter = iparam[2];
    let mode = iparam[6];
    let lworkl = 3 * ncv * ncv + 5 * ncv;
    let ierr = if n == 0 {
        -1
    } else if nev == 0 {
        -2
    } else if ncv <= nev || ncv > n {
        -3
    } else if mxiter <= 0 {
        -4
    } else if workl.len() < lworkl {
        -7
    } else if !(1..=3).contains(&mode) {
        -10
    } else if mode == 1 && bmat_g {
        -11
    } else {
        0
    };
    if ierr != 0 {
        out.info = ierr;
        return Ok(out);
    }
    if *tol <= 0.0 {
        *tol = EPS;
    }
    if !(0..=2).contains(&ishift) {
        ishift = 1;
    }
    if ishift != 1 {
        return Err("ZNAUPD with user shifts (ishift /= 1) is not ported".into());
    }
    let np = ncv - nev;
    workl[..lworkl].iter_mut().for_each(|x| *x = ZERO);
    let (ldh, ldq) = (ncv, ncv);
    let ih = 1;
    let ritz = ih + ldh * ncv;
    let bounds = ritz + ncv;
    let iq = bounds + ncv;
    let iw = iq + ldq * ncv;
    let next = iw + ncv * ncv + 3 * ncv;
    out.ipntr[3] = next;
    out.ipntr[4] = ih;
    out.ipntr[5] = ritz;
    out.ipntr[6] = iq;
    out.ipntr[7] = bounds;
    out.ipntr[13] = iw;

    let mut st = Stats::default();
    let (h, rest) = workl.split_at_mut(ldh * ncv);
    let (rz, rest) = rest.split_at_mut(ncv);
    let (bd, rest) = rest.split_at_mut(ncv);
    let (q, rest) = rest.split_at_mut(ldq * ncv);
    let wl = &mut rest[..ncv * ncv + 3 * ncv];
    let mut nev0 = nev;
    let mut np_ = np;
    let info = znaup2(
        op, &mut st, arp, bmat_g, n, which, &mut nev0, &mut np_, *tol, resid, ishift, &mut mxiter, v, ldv, h, ldh, rz,
        bd, q, ldq, wl, workd, info_in,
    )?;
    out.iparam[2] = mxiter;
    out.iparam[4] = np_ as i32;
    out.iparam[8] = st.nopx;
    out.iparam[9] = st.nbx;
    out.iparam[10] = st.nrorth;
    out.info = if info == 2 { 3 } else { info };
    out.stats = st;
    Ok(out)
}

/// `ZNAUP2`: the main loop of the implicitly restarted Arnoldi method.
/// Returns `info`.
fn znaup2(
    op: &mut ZOp<'_>,
    st: &mut Stats,
    arp: &mut ZArpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: &mut usize,
    np: &mut usize,
    tol: f64,
    resid: &mut [C],
    ishift: i32,
    mxiter: &mut i32,
    v: &mut [C],
    ldv: usize,
    h: &mut [C],
    ldh: usize,
    ritz: &mut [C],
    bounds: &mut [C],
    q: &mut [C],
    ldq: usize,
    workl: &mut [C],
    workd: &mut [C],
    info_in: i32,
) -> Result<i32, String> {
    let eps23 = EPS.powf(2.0 / 3.0);
    let nev0 = *nev;
    let np0 = *np;
    let kplusp = *nev + *np;
    let mut nconv;
    let mut iter = 0i32;
    let initv = info_in != 0;
    let mut info = 0i32;

    let (mut rnorm, _) = zgetv0(op, st, arp, bmat_g, initv, n, 1, v, ldv, resid, workd)?;
    if rnorm == 0.0 {
        return Ok(-9);
    }
    let r = znaitr(op, st, arp, bmat_g, n, 0, *nev, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
    if r > 0 {
        *np = r;
        *mxiter = iter;
        return Ok(-9999);
    }
    loop {
        // 1000
        iter += 1;
        *np = kplusp - *nev;
        // 20
        let r = znaitr(op, st, arp, bmat_g, n, *nev, *np, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
        if r > 0 {
            *np = r;
            *mxiter = iter;
            return Ok(-9999);
        }
        if zneigh(rnorm, kplusp, h, ldh, ritz, bounds, q, ldq, workl) != 0 {
            return Ok(-8);
        }
        *nev = nev0;
        *np = np0;
        let k2 = kplusp * kplusp;
        zcopy(kplusp, ritz, &mut workl[k2..]);
        zcopy(kplusp, bounds, &mut workl[k2 + kplusp..]);
        zngets(ishift, which, *nev, *np, ritz, bounds);
        nconv = 0;
        for i in 1..=*nev {
            let rz = ritz[*np + i - 1];
            let rtemp = eps23.max(dlapy2(rz.re, rz.im));
            let b = bounds[*np + i - 1];
            if dlapy2(b.re, b.im) <= tol * rtemp {
                nconv += 1;
            }
        }
        let nptemp = *np;
        for j in 1..=nptemp {
            if bounds[j - 1] == ZERO {
                *np -= 1;
                *nev += 1;
            }
        }
        if nconv >= nev0 || iter > *mxiter || *np == 0 {
            h[ix(3, 1, ldh)] = C::new(rnorm, 0.0);
            let wprime = match which {
                NWhich::LM => NWhich::SM,
                NWhich::SM => NWhich::LM,
                NWhich::LR => NWhich::SR,
                NWhich::SR => NWhich::LR,
                NWhich::LI => NWhich::SI,
                NWhich::SI => NWhich::LI,
            };
            zsortc(wprime, true, kplusp, ritz, bounds);
            for j in 0..nev0 {
                let rtemp = eps23.max(dlapy2(ritz[j].re, ritz[j].im));
                bounds[j] /= rtemp;
            }
            zsortc(NWhich::LM, true, nev0, bounds, ritz);
            for j in 0..nev0 {
                let rtemp = eps23.max(dlapy2(ritz[j].re, ritz[j].im));
                bounds[j] *= rtemp;
            }
            zsortc(which, true, nconv, ritz, bounds);
            if iter > *mxiter && nconv < nev0 {
                info = 1;
            }
            if *np == 0 && nconv < nev0 {
                info = 2;
            }
            *np = nconv;
            break;
        } else if nconv < nev0 && ishift == 1 {
            let nevbef = *nev;
            *nev += nconv.min(*np / 2);
            if *nev == 1 && kplusp >= 6 {
                *nev = kplusp / 2;
            } else if *nev == 1 && kplusp > 3 {
                *nev = 2;
            }
            *np = kplusp - *nev;
            if nevbef < *nev {
                zngets(ishift, which, *nev, *np, ritz, bounds);
            }
        }
        // 50: exact shifts are the unwanted Ritz values.
        let shifts: Vec<C> = ritz[..*np].to_vec();
        znapps(arp, n, *nev, *np, &shifts, v, ldv, h, ldh, resid, q, ldq, workd);
        // 100: the B-norm of the updated residual.
        if bmat_g {
            st.nbx += 1;
            zcopy(n, resid, &mut workd[n..]);
            op(2, [n, 0, 0], workd)?;
            rnorm = bnorm(zdotc(n, resid, workd));
        } else {
            zcopy(n, resid, workd);
            rnorm = dznrm2(n, resid);
        }
    }
    // 1100
    *mxiter = iter;
    *nev = nconv;
    Ok(info)
}
