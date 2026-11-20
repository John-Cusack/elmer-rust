//! `DSAUPD` and the routines under it (ARPACK 2.4): the implicitly
//! restarted Lanczos iteration for symmetric `OP` (`dsaup2`, `dsaitr`,
//! `dgetv0`, `dsapps`, `dseigt`, `dsconv`, `dsgets`, `dsortr`, `dsesrt`).
//!
//! ARPACK's reverse communication becomes a callback: where the Fortran
//! sets `ido` and returns to its caller, [`dsaupd`] calls `op(ido, ipntr,
//! workd)` and carries on where the Fortran would resume. `ipntr` holds
//! 0-based offsets into `workd`:
//!
//! * `ido = -1`: `y = OP x`, `x` at `ipntr[0]`, `y` to `ipntr[1]` (the first
//!   call, `B x` not provided);
//! * `ido = 1`: `y = OP x` with `B x` already at `ipntr[2]` (mode 3);
//! * `ido = 2`: `y = B x`, `x` at `ipntr[0]`, `y` to `ipntr[1]`.
//!
//! The whole `workd` is handed over, as the Fortran caller sees it, so a
//! caller that uses its other parts (Elmer's `ArpackEigenSolve` overwrites
//! `x` and starts its inner solves from what `y` holds) does the same.
//!
//! Debug output (`msglvl > 0`) and timing are left out; the counters of
//! `stat.h` are kept in [`Stats`].

use super::blas::{daxpy, dcopy, dgemv_n, dgemv_t, dscal, dswap, dswap_within};
use super::lapack::{dlartg, dlascl_g, dlaset, dlarnv, dstqrb, EPS, SAFMIN};
use super::par::{eps23, pddot, pdnorm2, psum, Par};

/// The callback of the reverse communication: `(ido, ipntr, workd)`.
pub type Op<'a> = dyn FnMut(i32, [usize; 3], &mut [f64]) -> Result<(), String> + 'a;

/// ARPACK's process state: the seed `dgetv0` keeps (`SAVE`, starting at
/// 1, 3, 5, 7), so successive starting vectors continue one random stream,
/// and the `smlnum` that `dnaitr` and `dnapps` each `SAVE` on their first
/// call (computed with that call's `n`).
#[derive(Clone, Debug)]
pub struct Arpack {
    pub iseed: [i32; 4],
    pub naitr_smlnum: Option<f64>,
    pub napps_smlnum: Option<f64>,
}

impl Default for Arpack {
    fn default() -> Self {
        Arpack { iseed: [1, 3, 5, 7], naitr_smlnum: None, napps_smlnum: None }
    }
}

/// The counters of `stat.h` (`dstats` zeroes them at each `dsaupd` start).
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub nopx: i32,
    pub nbx: i32,
    pub nrorth: i32,
    pub nitref: i32,
    pub nrstrt: i32,
}

/// `0.717` is a default-REAL literal in ARPACK: single precision.
pub(crate) fn c0717() -> f64 {
    f64::from(0.717f32)
}

/// `which` as ARPACK spells it.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Which {
    LM,
    SM,
    LA,
    SA,
    BE,
}

/// `DSORTR(which, apply, n, x1, x2)`: shell sort of `x1` (and `x2` with it)
/// so that the wanted end comes last: `SA` decreasing, `SM` decreasing in
/// magnitude, `LA` increasing, `LM` increasing in magnitude.
pub fn dsortr(which: Which, apply: bool, n: usize, x1: &mut [f64], x2: &mut [f64]) {
    let swap_if = |a: f64, b: f64| match which {
        Which::SA => a < b,
        Which::SM => a.abs() < b.abs(),
        Which::LA => a > b,
        _ => a.abs() > b.abs(),
    };
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let ju = j as usize;
                if swap_if(x1[ju], x1[ju + igap]) {
                    x1.swap(ju, ju + igap);
                    if apply {
                        x2.swap(ju, ju + igap);
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

/// `DSESRT(which, apply, n, x, na, a, lda)`: [`dsortr`] carrying the columns
/// of `a` (`na` rows) along.
pub fn dsesrt(which: Which, apply: bool, n: usize, x: &mut [f64], na: usize, a: &mut [f64], lda: usize) {
    let swap_if = |p: f64, q: f64| match which {
        Which::SA => p < q,
        Which::SM => p.abs() < q.abs(),
        Which::LA => p > q,
        _ => p.abs() > q.abs(),
    };
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let ju = j as usize;
                if swap_if(x[ju], x[ju + igap]) {
                    x.swap(ju, ju + igap);
                    if apply {
                        let (lo, hi) = a.split_at_mut((ju + igap) * lda);
                        dswap(na, &mut lo[ju * lda..], hi);
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

/// `DSCONV(n, ritz, bounds, tol, nconv)`: Ritz values whose error bound is
/// within `tol` relative to `max(eps**(2/3), |ritz|)`.
pub fn dsconv(n: usize, ritz: &[f64], bounds: &[f64], tol: f64) -> usize {
    let eps23 = EPS.powf(2.0 / 3.0);
    (0..n).filter(|&i| bounds[i] <= tol * eps23.max(ritz[i].abs())).count()
}

/// `DSGETS(ishift, which, kev, np, ritz, bounds, shifts)`: sort the Ritz
/// values so the `np` unwanted ones come first; with exact shifts those,
/// ordered by error bound, are the shifts.
#[allow(clippy::too_many_arguments)]
pub fn dsgets(ishift: i32, which: Which, kev: usize, np: usize, ritz: &mut [f64], bounds: &mut [f64], shifts: &mut [f64]) {
    if which == Which::BE {
        dsortr(Which::LA, true, kev + np, ritz, bounds);
        let kevd2 = kev / 2;
        if kev > 1 {
            let m = kevd2.min(np);
            let b = kevd2.max(np);
            dswap_within(m, ritz, 0, b);
            dswap_within(m, bounds, 0, b);
        }
    } else {
        dsortr(which, true, kev + np, ritz, bounds);
    }
    if ishift == 1 && np > 0 {
        dsortr(Which::SM, true, np, bounds, ritz);
        dcopy(np, ritz, shifts);
    }
}

/// `DSEIGT(rnorm, n, h, ldh, eig, bounds, workl, ierr)`: eigenvalues of the
/// tridiagonal `H` and their error bounds `rnorm * |last eigenvector entry|`.
pub fn dseigt(rnorm: f64, n: usize, h: &[f64], ldh: usize, eig: &mut [f64], bounds: &mut [f64], workl: &mut [f64]) -> usize {
    dcopy(n, &h[ldh..], eig);
    dcopy(n - 1, &h[1..], workl);
    let (e, work) = workl.split_at_mut(n);
    let ierr = dstqrb(n, eig, e, bounds, work);
    if ierr != 0 {
        return ierr;
    }
    for b in &mut bounds[..n] {
        *b = rnorm * b.abs();
    }
    0
}

/// `DGETV0(ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd,
/// ierr)`: a (random, unless `initv`) starting vector in the range of `OP`,
/// B-orthogonal to the first `j - 1` columns of `v`. Returns `(rnorm, ierr)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn dgetv0(
    par: Par<'_>,
    op: &mut Op<'_>,
    st: &mut Stats,
    arp: &mut Arpack,
    bmat_g: bool,
    initv: bool,
    n: usize,
    j: usize,
    v: &[f64],
    ldv: usize,
    resid: &mut [f64],
    workd: &mut [f64],
) -> Result<(f64, i32), String> {
    let mut ierr = 0;
    let mut iter = 0;
    if !initv {
        dlarnv(2, &mut arp.iseed, n, resid);
    }
    if bmat_g {
        st.nopx += 1;
        dcopy(n, resid, workd);
        op(-1, [0, n, 0], workd)?;
        st.nbx += 1;
        let (lo, hi) = workd.split_at_mut(n);
        let _ = lo;
        dcopy(n, hi, resid);
        op(2, [n, 0, 0], workd)?;
    } else {
        dcopy(n, resid, workd);
    }
    // 20
    let mut rnorm0 = if bmat_g { pddot(par, n, resid, workd).abs().sqrt() } else { pdnorm2(par, n, resid) };
    let mut rnorm = rnorm0;
    if j == 1 {
        return Ok((rnorm, ierr));
    }
    loop {
        // 30
        if par.is_some() {
            // PARPACK reduces V^T B r in its own workspace, not workd(n+1).
            let mut t = vec![0.0; j - 1];
            dgemv_t(n, j - 1, 1.0, v, ldv, &workd[..n], 0.0, &mut t);
            psum(par, &mut t);
            dgemv_n(n, j - 1, -1.0, v, ldv, &t, 1.0, resid);
        } else {
            let (lo, hi) = workd.split_at_mut(n);
            dgemv_t(n, j - 1, 1.0, v, ldv, lo, 0.0, hi);
            dgemv_n(n, j - 1, -1.0, v, ldv, hi, 1.0, resid);
        }
        if bmat_g {
            st.nbx += 1;
            dcopy(n, resid, &mut workd[n..]);
            op(2, [n, 0, 0], workd)?;
        } else {
            dcopy(n, resid, workd);
        }
        // 40
        rnorm = if bmat_g { pddot(par, n, resid, workd).abs().sqrt() } else { pdnorm2(par, n, resid) };
        if rnorm > c0717() * rnorm0 {
            break;
        }
        iter += 1;
        if iter <= 1 {
            rnorm0 = rnorm;
            continue;
        }
        resid[..n].iter_mut().for_each(|x| *x = 0.0);
        rnorm = 0.0;
        ierr = -1;
        break;
    }
    Ok((rnorm, ierr))
}

/// `DSAITR(ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh, ipntr,
/// workd, info)`: extend a `k`-step Lanczos factorisation by `np` steps.
/// Returns `info` (0, or the step at which no new starting vector could be
/// found).
#[allow(clippy::too_many_arguments)]
fn dsaitr(
    par: Par<'_>,
    op: &mut Op<'_>,
    st: &mut Stats,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    k: usize,
    np: usize,
    mode: i32,
    resid: &mut [f64],
    rnorm: &mut f64,
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    workd: &mut [f64],
) -> Result<usize, String> {
    let (ipj, irj, ivj) = (0usize, n, 2 * n);
    let mut j = k + 1;
    let mut rstart = false;
    // h(i,1) = h[i-1], h(i,2) = h[ldh+i-1].
    loop {
        // 1000
        if *rnorm <= 0.0 {
            st.nrstrt += 1;
            let mut itry = 1;
            loop {
                // 20
                rstart = true;
                let (rn, ierr) = dgetv0(par, op, st, arp, bmat_g, false, n, j, v, ldv, resid, workd)?;
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
        dcopy(n, resid, &mut v[vj..]);
        if *rnorm >= SAFMIN {
            let temp1 = 1.0 / *rnorm;
            dscal(n, temp1, &mut v[vj..]);
            dscal(n, temp1, &mut workd[ipj..]);
        } else {
            dlascl_g(*rnorm, 1.0, n, 1, &mut v[vj..], n);
            dlascl_g(*rnorm, 1.0, n, 1, &mut workd[ipj..], n);
        }
        st.nopx += 1;
        dcopy(n, &v[vj..], &mut workd[ivj..]);
        op(1, [ivj, irj, ipj], workd)?;
        // 50
        dcopy(n, &workd[irj..], resid);
        if mode != 2 {
            if bmat_g {
                st.nbx += 1;
                op(2, [irj, ipj, 0], workd)?;
            } else {
                dcopy(n, resid, &mut workd[ipj..]);
            }
        }
        // 65
        let wnorm = if mode == 2 {
            pddot(par, n, resid, &workd[ivj..]).abs().sqrt()
        } else if bmat_g {
            pddot(par, n, resid, &workd[ipj..]).abs().sqrt()
        } else {
            pdnorm2(par, n, resid)
        };
        if par.is_some() {
            // PARPACK: V^T B r reduced in its workspace; workd(irj) untouched.
            let mut t = vec![0.0; j];
            let x = if mode != 2 { &workd[ipj..] } else { &workd[ivj..] };
            dgemv_t(n, j, 1.0, v, ldv, x, 0.0, &mut t);
            psum(par, &mut t);
            dgemv_n(n, j, -1.0, v, ldv, &t, 1.0, resid);
            h[ldh + j - 1] = t[j - 1];
        } else {
            {
                let (lo, hi) = workd.split_at_mut(ivj);
                if mode != 2 {
                    let (p, r) = lo.split_at_mut(irj);
                    dgemv_t(n, j, 1.0, v, ldv, p, 0.0, r);
                } else {
                    dgemv_t(n, j, 1.0, v, ldv, hi, 0.0, &mut lo[irj..]);
                }
            }
            dgemv_n(n, j, -1.0, v, ldv, &workd[irj..], 1.0, resid);
            h[ldh + j - 1] = workd[irj + j - 1];
        }
        h[j - 1] = if j == 1 || rstart { 0.0 } else { *rnorm };
        let mut iter = 0;
        if bmat_g {
            st.nbx += 1;
            dcopy(n, resid, &mut workd[irj..]);
            op(2, [irj, ipj, 0], workd)?;
        } else {
            dcopy(n, resid, &mut workd[ipj..]);
        }
        // 70
        *rnorm = if bmat_g { pddot(par, n, resid, &workd[ipj..]).abs().sqrt() } else { pdnorm2(par, n, resid) };
        // `IF (rnorm .GT. 0.717*wnorm) GO TO 100`: a NaN norm refines too.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(*rnorm > c0717() * wnorm) {
            st.nrorth += 1;
            loop {
                // 80
                if par.is_some() {
                    let mut t = vec![0.0; j];
                    dgemv_t(n, j, 1.0, v, ldv, &workd[ipj..], 0.0, &mut t);
                    psum(par, &mut t);
                    dgemv_n(n, j, -1.0, v, ldv, &t, 1.0, resid);
                    if j == 1 || rstart {
                        h[j - 1] = 0.0;
                    }
                    h[ldh + j - 1] += t[j - 1];
                } else {
                    {
                        let (p, r) = workd.split_at_mut(irj);
                        dgemv_t(n, j, 1.0, v, ldv, p, 0.0, r);
                    }
                    dgemv_n(n, j, -1.0, v, ldv, &workd[irj..], 1.0, resid);
                    if j == 1 || rstart {
                        h[j - 1] = 0.0;
                    }
                    h[ldh + j - 1] += workd[irj + j - 1];
                }
                if bmat_g {
                    st.nbx += 1;
                    dcopy(n, resid, &mut workd[irj..]);
                    op(2, [irj, ipj, 0], workd)?;
                } else {
                    dcopy(n, resid, &mut workd[ipj..]);
                }
                // 90
                let rnorm1 = if bmat_g { pddot(par, n, resid, &workd[ipj..]).abs().sqrt() } else { pdnorm2(par, n, resid) };
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
                resid[..n].iter_mut().for_each(|x| *x = 0.0);
                *rnorm = 0.0;
                break;
            }
        }
        // 100
        rstart = false;
        if h[j - 1] < 0.0 {
            h[j - 1] = -h[j - 1];
            if j < k + np {
                dscal(n, -1.0, &mut v[j * ldv..]);
            } else {
                dscal(n, -1.0, resid);
            }
        }
        j += 1;
        if j > k + np {
            return Ok(0);
        }
    }
}

/// `DSAPPS(n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd)`: apply
/// `np` implicit shifts to the factorisation, keeping `kev` vectors.
#[allow(clippy::too_many_arguments)]
fn dsapps(
    par: Par<'_>,
    n: usize,
    kev: usize,
    np: usize,
    shift: &[f64],
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    resid: &mut [f64],
    q: &mut [f64],
    ldq: usize,
    workd: &mut [f64],
) {
    let epsmch = EPS;
    let kplusp = kev + np;
    // 1-based accessors: h(i,1), h(i,2), q(i,j).
    macro_rules! h1 { ($i:expr) => { h[$i - 1] } }
    macro_rules! h2 { ($i:expr) => { h[ldh + $i - 1] } }
    macro_rules! qq { ($i:expr, $j:expr) => { q[($j - 1) * ldq + $i - 1] } }
    let mut itop = 1usize;
    dlaset(kplusp, kplusp, 0.0, 1.0, q, ldq);
    if np == 0 {
        return;
    }
    for jj in 1..=np {
        let mut istart = itop;
        loop {
            // 20
            let mut iend = kplusp;
            for i in istart..kplusp {
                let big = h2!(i).abs() + h2!(i + 1).abs();
                if h1!(i + 1) <= epsmch * big {
                    h1!(i + 1) = 0.0;
                    iend = i;
                    break;
                }
            }
            // 40
            if istart < iend {
                let f = h2!(istart) - shift[jj - 1];
                let g = h1!(istart + 1);
                let (mut c, mut s, _) = dlartg(f, g);
                let a1 = c * h2!(istart) + s * h1!(istart + 1);
                let a2 = c * h1!(istart + 1) + s * h2!(istart + 1);
                let a4 = c * h2!(istart + 1) - s * h1!(istart + 1);
                let a3 = c * h1!(istart + 1) - s * h2!(istart);
                h2!(istart) = c * a1 + s * a2;
                h2!(istart + 1) = c * a4 - s * a3;
                h1!(istart + 1) = c * a3 + s * a4;
                // DO 60 j = 1, min(istart+jj, kplusp); the loop leaves j at
                // its bound + 1, which the next loop's bound reads.
                let ub = (istart + jj).min(kplusp);
                for j in 1..=ub {
                    let a1 = c * qq!(j, istart) + s * qq!(j, istart + 1);
                    qq!(j, istart + 1) = -s * qq!(j, istart) + c * qq!(j, istart + 1);
                    qq!(j, istart) = a1;
                }
                let mut jlast = ub + 1;
                for i in istart + 1..iend {
                    let f = h1!(i);
                    let g = s * h1!(i + 1);
                    h1!(i + 1) *= c;
                    let (cc, ss, mut r) = dlartg(f, g);
                    c = cc;
                    s = ss;
                    if r < 0.0 {
                        r = -r;
                        c = -c;
                        s = -s;
                    }
                    h1!(i) = r;
                    let a1 = c * h2!(i) + s * h1!(i + 1);
                    let a2 = c * h1!(i + 1) + s * h2!(i + 1);
                    let a3 = c * h1!(i + 1) - s * h2!(i);
                    let a4 = c * h2!(i + 1) - s * h1!(i + 1);
                    h2!(i) = c * a1 + s * a2;
                    h2!(i + 1) = c * a4 - s * a3;
                    h1!(i + 1) = c * a3 + s * a4;
                    // DO 50 j = 1, min(j+jj, kplusp) with j from the loop
                    // before; PARPACK's pdsapps bounds it by min(i+jj, kplusp).
                    let ub = if par.is_some() { (i + jj).min(kplusp) } else { (jlast + jj).min(kplusp) };
                    for j in 1..=ub {
                        let a1 = c * qq!(j, i) + s * qq!(j, i + 1);
                        qq!(j, i + 1) = -s * qq!(j, i) + c * qq!(j, i + 1);
                        qq!(j, i) = a1;
                    }
                    jlast = ub + 1;
                }
            }
            istart = iend + 1;
            if h1!(iend) < 0.0 {
                h1!(iend) = -h1!(iend);
                dscal(kplusp, -1.0, &mut q[(iend - 1) * ldq..]);
            }
            if iend < kplusp {
                continue;
            }
            break;
        }
        // DO 80 i = itop, kplusp-1 (bounds fixed at entry, as in Fortran).
        let (lo, hi) = (itop, kplusp);
        for i in lo..hi {
            if h1!(i + 1) > 0.0 {
                break;
            }
            itop += 1;
        }
    }
    // 90
    for i in itop..kplusp {
        let big = h2!(i).abs() + h2!(i + 1).abs();
        if h1!(i + 1) <= epsmch * big {
            h1!(i + 1) = 0.0;
        }
    }
    if h1!(kev + 1) > 0.0 {
        dgemv_n(n, kplusp, 1.0, v, ldv, &q[kev * ldq..], 0.0, &mut workd[n..]);
    }
    for i in 1..=kev {
        dgemv_n(n, kplusp - i + 1, 1.0, v, ldv, &q[(kev - i) * ldq..], 0.0, workd);
        dcopy(n, workd, &mut v[(kplusp - i) * ldv..]);
    }
    // DLACPY('All', n, kev, v(1,np+1), ldv, v, ldv), column by column.
    for j in 0..kev {
        v.copy_within((np + j) * ldv..(np + j) * ldv + n, j * ldv);
    }
    if h1!(kev + 1) > 0.0 {
        dcopy(n, &workd[n..], &mut v[kev * ldv..]);
    }
    dscal(n, qq!(kplusp, kev), resid);
    if h1!(kev + 1) > 0.0 {
        daxpy(n, h1!(kev + 1), &v[kev * ldv..], resid);
    }
}

/// The parameters and results of [`dsaupd`] (`iparam`, `ipntr` and `info`).
#[derive(Clone, Debug, Default)]
pub struct Saupd {
    /// `IPARAM(1..11)` (1-based meaning; index 0 is `IPARAM(1)`).
    pub iparam: [i32; 11],
    /// `IPNTR(1..11)`, 1-based offsets into `workl` (as the Fortran).
    pub ipntr: [usize; 11],
    pub info: i32,
    pub stats: Stats,
}

/// `DSAUPD(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr,
/// workd, workl, lworkl, info)` driven to completion. `iparam` carries
/// `ishift`, `mxiter` and `mode` in (`IPARAM(1)`, `(3)`, `(7)`); `info` in
/// is 0 for a random starting vector, nonzero to use `resid`. `workd` has
/// `3 n`, `workl` at least `ncv**2 + 8 ncv` entries.
#[allow(clippy::too_many_arguments)]
pub fn dsaupd(
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: &mut f64,
    resid: &mut [f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    iparam: [i32; 11],
    workd: &mut [f64],
    workl: &mut [f64],
    info_in: i32,
) -> Result<Saupd, String> {
    saupd(None, op, arp, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, iparam, workd, workl, info_in)
}

/// PARPACK's `PDSAUPD(comm, ...)`: [`dsaupd`] on this rank's `n` rows of
/// the vectors, the reductions over `comm` (see [`super::par`]). Collective.
#[allow(clippy::too_many_arguments)]
pub fn pdsaupd(
    comm: &dyn elmer_core::Communicator,
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: &mut f64,
    resid: &mut [f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    iparam: [i32; 11],
    workd: &mut [f64],
    workl: &mut [f64],
    info_in: i32,
) -> Result<Saupd, String> {
    saupd(Some(comm), op, arp, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, iparam, workd, workl, info_in)
}

#[allow(clippy::too_many_arguments)]
fn saupd(
    par: Par<'_>,
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: usize,
    tol: &mut f64,
    resid: &mut [f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    iparam: [i32; 11],
    workd: &mut [f64],
    workl: &mut [f64],
    info_in: i32,
) -> Result<Saupd, String> {
    let mut out = Saupd { iparam, ..Default::default() };
    let ishift = iparam[0];
    let mut mxiter = iparam[2];
    let mode = iparam[6];
    let mut ierr = 0;
    if n == 0 && par.is_none() {
        ierr = -1;
    } else if nev == 0 {
        ierr = -2;
    } else if ncv <= nev || (ncv > n && par.is_none()) {
        ierr = -3;
    }
    if mxiter <= 0 {
        ierr = -4;
    }
    if workl.len() < ncv * ncv + 8 * ncv {
        ierr = -7;
    }
    if !(1..=5).contains(&mode) {
        ierr = -10;
    } else if mode == 1 && bmat_g {
        ierr = -11;
    } else if !(0..=1).contains(&ishift) {
        ierr = -12;
    } else if nev == 1 && which == Which::BE {
        ierr = -13;
    }
    if ishift == 0 {
        return Err("DSAUPD with user shifts (ishift = 0) is not ported".into());
    }
    if ierr != 0 {
        out.info = ierr;
        return Ok(out);
    }
    if *tol <= 0.0 {
        *tol = EPS;
    }
    let np = ncv - nev;
    workl[..ncv * ncv + 8 * ncv].iter_mut().for_each(|x| *x = 0.0);
    let ldh = ncv;
    let ldq = ncv;
    // 1-based, as IPNTR reports them.
    let ih = 1;
    let ritz = ih + 2 * ldh;
    let bounds = ritz + ncv;
    let iq = bounds + ncv;
    let iw = iq + ncv * ncv;
    let next = iw + 3 * ncv;
    out.ipntr[3] = next;
    out.ipntr[4] = ih;
    out.ipntr[5] = ritz;
    out.ipntr[6] = bounds;
    out.ipntr[10] = iw;

    let mut st = Stats::default();
    let (h, rest) = workl.split_at_mut(2 * ldh);
    let (ritz_s, rest) = rest.split_at_mut(ncv);
    let (bounds_s, rest) = rest.split_at_mut(ncv);
    let (q, rest) = rest.split_at_mut(ncv * ncv);
    let wl = &mut rest[..3 * ncv];
    let mut nev0 = nev;
    let mut np_ = np;
    let info = dsaup2(
        par, op, &mut st, arp, bmat_g, n, which, &mut nev0, &mut np_, *tol, resid, mode, ishift, &mut mxiter, v, ldv, h, ldh,
        ritz_s, bounds_s, q, ldq, wl, workd, info_in,
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

/// `DSAUP2`: the main loop of the implicitly restarted Lanczos method.
/// Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dsaup2(
    par: Par<'_>,
    op: &mut Op<'_>,
    st: &mut Stats,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: Which,
    nev: &mut usize,
    np: &mut usize,
    tol: f64,
    resid: &mut [f64],
    mode: i32,
    ishift: i32,
    mxiter: &mut i32,
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    ritz: &mut [f64],
    bounds: &mut [f64],
    q: &mut [f64],
    ldq: usize,
    workl: &mut [f64],
    workd: &mut [f64],
    info_in: i32,
) -> Result<i32, String> {
    let eps23 = eps23(par);
    let nev0 = *nev;
    let np0 = *np;
    let kplusp = nev0 + np0;
    let mut nconv;
    let mut iter = 0i32;
    let initv = info_in != 0;
    let mut info = 0i32;

    let (mut rnorm, _) = dgetv0(par, op, st, arp, bmat_g, initv, n, 1, v, ldv, resid, workd)?;
    if rnorm == 0.0 {
        return Ok(-9);
    }
    let r = dsaitr(par, op, st, arp, bmat_g, n, 0, nev0, mode, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
    if r > 0 {
        *np = r;
        *mxiter = iter;
        return Ok(-9999);
    }
    loop {
        // 1000
        iter += 1;
        // 20
        let r = dsaitr(par, op, st, arp, bmat_g, n, *nev, *np, mode, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
        if r > 0 {
            *np = r;
            *mxiter = iter;
            return Ok(-9999);
        }
        if dseigt(rnorm, kplusp, h, ldh, ritz, bounds, workl) != 0 {
            return Ok(-8);
        }
        dcopy(kplusp, ritz, &mut workl[kplusp..]);
        dcopy(kplusp, bounds, &mut workl[2 * kplusp..]);
        *nev = nev0;
        *np = np0;
        dsgets(ishift, which, *nev, *np, ritz, bounds, workl);
        {
            let np_ = *np;
            let nev_ = *nev;
            let b: Vec<f64> = bounds[np_..np_ + nev_].to_vec();
            workl[np_..np_ + nev_].copy_from_slice(&b);
            nconv = dsconv(nev_, &ritz[np_..], &workl[np_..], tol);
        }
        let nptemp = *np;
        for &b in &bounds[..nptemp] {
            if b == 0.0 {
                *np -= 1;
                *nev += 1;
            }
        }
        if nconv >= nev0 || iter > *mxiter || *np == 0 {
            if which == Which::BE {
                dsortr(Which::SA, true, kplusp, ritz, bounds);
                // PARPACK halves nev0, ARPACK the current nev.
                let nevd2 = if par.is_some() { nev0 / 2 } else { *nev / 2 };
                let nevm2 = if par.is_some() { nev0 - nevd2 } else { *nev - nevd2 };
                if *nev > 1 {
                    let m = nevd2.min(*np);
                    let a = nevm2;
                    let b = (kplusp - nevd2 + 1).max(kplusp - *np + 1) - 1;
                    dswap_within(m, ritz, a, b);
                    let b2 = if par.is_some() {
                        (kplusp - nevd2 + 1).max(kplusp - *np + 1) - 1
                    } else {
                        (kplusp - nevd2 + 1).max(kplusp - *np)
                    };
                    dswap_within(m, bounds, a, b2);
                }
            } else {
                let wprime = match which {
                    Which::LM => Which::SM,
                    Which::SM => Which::LM,
                    Which::LA => Which::SA,
                    _ => Which::LA,
                };
                dsortr(wprime, true, kplusp, ritz, bounds);
            }
            for j in 0..nev0 {
                let temp = eps23.max(ritz[j].abs());
                bounds[j] /= temp;
            }
            dsortr(Which::LA, true, nev0, bounds, ritz);
            for j in 0..nev0 {
                let temp = eps23.max(ritz[j].abs());
                bounds[j] *= temp;
            }
            if which == Which::BE {
                dsortr(Which::LA, true, nconv, ritz, bounds);
            } else {
                dsortr(which, true, nconv, ritz, bounds);
            }
            h[0] = rnorm;
            if iter > *mxiter && nconv < *nev {
                info = 1;
            }
            if *np == 0 && nconv < nev0 {
                info = 2;
            }
            *np = nconv;
            break;
        } else if nconv < *nev && ishift == 1 {
            let nevbef = *nev;
            *nev += nconv.min(*np / 2);
            if *nev == 1 && kplusp >= 6 {
                *nev = kplusp / 2;
            } else if *nev == 1 && kplusp > 2 {
                *nev = 2;
            }
            *np = kplusp - *nev;
            if nevbef < *nev {
                dsgets(ishift, which, *nev, *np, ritz, bounds, workl);
            }
        }
        // 50: the shifts are the first np Ritz values.
        let shifts: Vec<f64> = ritz[..*np].to_vec();
        dsapps(par, n, *nev, *np, &shifts, v, ldv, h, ldh, resid, q, ldq, workd);
        if bmat_g {
            st.nbx += 1;
            dcopy(n, resid, &mut workd[n..]);
            op(2, [n, 0, 0], workd)?;
        } else {
            dcopy(n, resid, workd);
        }
        // 100
        rnorm = if bmat_g { pddot(par, n, resid, workd).abs().sqrt() } else { pdnorm2(par, n, resid) };
    }
    // 1100
    *mxiter = iter;
    *nev = nconv;
    Ok(info)
}
