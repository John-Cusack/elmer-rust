//! `DNAUPD` and the routines under it (ARPACK 2.4, as Elmer's
//! `mathlibs/src/arpack` builds it): the implicitly restarted Arnoldi
//! iteration for non-symmetric `OP` (`dnaup2`, `dnaitr`, `dgetv0`, `dnapps`,
//! `dneigh`, `dlaqrb`, `dngets`, `dnconv`, `dsortc`).
//!
//! The reverse communication is the callback of [`super::saupd`]
//! (`op(ido, ipntr, workd)`, 0-based offsets): `ido = -1` and `1` ask for
//! `y = OP x` (`x` at `ipntr[0]`, `y` to `ipntr[1]`, and for `1` in mode 3
//! `B x` at `ipntr[2]`), `ido = 2` for `y = B x`. User shifts (`ishift =
//! 0`, `ido = 3`) are not ported.
//!
//! `dnaitr` and `dnapps` compute their `smlnum = unfl*(n/ulp)` once per
//! process (`SAVE`d after the first call, with that call's `n`); the value
//! is kept in [`Arpack`] with `dgetv0`'s seed.

use super::blas::{daxpy, dcopy, dgemv_n, dgemv_t, dnrm2, dscal};
use super::lapack::{dlacpy, dlapy2, dlarf, dlarfg, dlartg, dlascl_g, dlaset, EPS, PREC, SAFMIN};
use super::lapack_ns::{dlanhs_one, dlanv2, drot_within, dtrevc_right};
use super::par::{eps23, pddot, pdnorm2, psum, Par};
use super::saupd::{c0717, dgetv0, Arpack, Op, Stats};

/// `which` of the non-symmetric routines.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum NWhich {
    LM,
    SM,
    LR,
    SR,
    LI,
    SI,
}

#[inline]
fn ix(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

/// `DSORTC(which, apply, n, xreal, ximag, y)`: shell sort of the complex
/// numbers `xreal + i ximag` (and `y` with them, when `apply`) so that the
/// wanted end comes last: `LM` increasing magnitude, `SM` decreasing, `LR`
/// increasing real part, `SR` decreasing, `LI` increasing `|imag|`, `SI`
/// decreasing.
pub fn dsortc(which: NWhich, apply: bool, n: usize, xreal: &mut [f64], ximag: &mut [f64], y: &mut [f64]) {
    let swap = |xr: &[f64], xi: &[f64], j: usize, k: usize| -> bool {
        match which {
            NWhich::LM => dlapy2(xr[j], xi[j]) > dlapy2(xr[k], xi[k]),
            NWhich::SM => dlapy2(xr[j], xi[j]) < dlapy2(xr[k], xi[k]),
            NWhich::LR => xr[j] > xr[k],
            NWhich::SR => xr[j] < xr[k],
            NWhich::LI => xi[j].abs() > xi[k].abs(),
            NWhich::SI => xi[j].abs() < xi[k].abs(),
        }
    };
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let (ju, ku) = (j as usize, j as usize + igap);
                if !swap(xreal, ximag, ju, ku) {
                    break;
                }
                xreal.swap(ju, ku);
                ximag.swap(ju, ku);
                if apply {
                    y.swap(ju, ku);
                }
                j -= igap as isize;
            }
        }
        igap /= 2;
    }
}

/// `DNCONV(n, ritzr, ritzi, bounds, tol, nconv)`: the number of Ritz values
/// whose estimate is at most `tol * max(eps23, |ritz|)`.
pub fn dnconv(n: usize, ritzr: &[f64], ritzi: &[f64], bounds: &[f64], tol: f64) -> usize {
    let eps23 = EPS.powf(2.0 / 3.0);
    (0..n).filter(|&i| bounds[i] <= tol * eps23.max(dlapy2(ritzr[i], ritzi[i]))).count()
}

/// `DNGETS(ishift, which, kev, np, ritzr, ritzi, bounds, shiftr, shifti)`
/// with exact shifts: sorts the Ritz values so the `np` unwanted ones come
/// first, keeps a complex pair on one side (moving the split by one), and
/// sorts the shifts by their estimates (largest last, `SR` of the bounds).
pub fn dngets(
    ishift: i32,
    which: NWhich,
    kev: &mut usize,
    np: &mut usize,
    ritzr: &mut [f64],
    ritzi: &mut [f64],
    bounds: &mut [f64],
) -> Result<(), String> {
    let n = *kev + *np;
    let pre = match which {
        NWhich::LM => NWhich::LR,
        NWhich::SM => NWhich::SR,
        NWhich::LR => NWhich::LM,
        NWhich::SR => NWhich::SM,
        NWhich::LI => NWhich::LM,
        NWhich::SI => NWhich::SM,
    };
    dsortc(pre, true, n, ritzr, ritzi, bounds);
    dsortc(which, true, n, ritzr, ritzi, bounds);
    if *np == 0 {
        // `ritzr(np)` would be `ritzr(0)`, outside the array.
        return Err("DNGETS with np = 0".into());
    }
    let p = *np;
    if (ritzr[p] - ritzr[p - 1]) == 0.0 && (ritzi[p] + ritzi[p - 1]) == 0.0 {
        *np -= 1;
        *kev += 1;
    }
    if ishift == 1 {
        let p = *np;
        dsortc(NWhich::SR, true, p, bounds, ritzr, ritzi);
    }
    Ok(())
}

/// `DLAQRB(wantt, n, ilo, ihi, h, ldh, wr, wi, z, info)`: the eigenvalues
/// (and with `wantt` the Schur form) of the Hessenberg `h`, computing only
/// the last components of the Schur vectors into `z`. Returns `info`.
#[allow(clippy::too_many_arguments)]
pub fn dlaqrb(
    wantt: bool,
    n: usize,
    ilo: usize,
    ihi: usize,
    h: &mut [f64],
    ldh: usize,
    wr: &mut [f64],
    wi: &mut [f64],
    z: &mut [f64],
) -> usize {
    const DAT1: f64 = 7.5e-1;
    const DAT2: f64 = -4.375e-1;
    if n == 0 {
        return 0;
    }
    if ilo == ihi {
        wr[ilo - 1] = h[ix(ilo, ilo, ldh)];
        wi[ilo - 1] = 0.0;
        return 0;
    }
    // Initialize the vector of last components of the Schur vectors for
    // accumulation.
    for zj in z.iter_mut().take(n - 1) {
        *zj = 0.0;
    }
    z[n - 1] = 1.0;
    let nh = ihi - ilo + 1;
    let unfl = SAFMIN;
    let ulp = PREC;
    let smlnum = unfl * (nh as f64 / ulp);
    let (mut i1, mut i2) = (0usize, 0usize);
    if wantt {
        i1 = 1;
        i2 = n;
        for i in 1..=i2.saturating_sub(2) {
            h[ix(i1 + i + 1, i, ldh)] = 0.0;
        }
    } else {
        for i in 1..=(ihi - ilo).saturating_sub(1) {
            h[ix(ilo + i + 1, ilo + i - 1, ldh)] = 0.0;
        }
    }
    // ITN is the total number of QR iterations allowed.
    let mut itn = 30 * nh;
    let mut i = ihi as isize;
    let mut v: [f64; 3];
    loop {
        // 10
        let mut l = ilo;
        if i < ilo as isize {
            return 0;
        }
        let iu = i as usize;
        let mut done = None;
        // The loop bound is evaluated once, on entry.
        for its in 0..=itn {
            // Look for a single small subdiagonal element.
            let mut k = iu;
            while k > l {
                let mut tst1 = h[ix(k - 1, k - 1, ldh)].abs() + h[ix(k, k, ldh)].abs();
                if tst1 == 0.0 {
                    tst1 = dlanhs_one(iu - l + 1, &h[ix(l, l, ldh)..], ldh);
                }
                if h[ix(k, k - 1, ldh)].abs() <= (ulp * tst1).max(smlnum) {
                    break;
                }
                k -= 1;
            }
            l = k;
            if l > ilo {
                // H(L,L-1) is negligible.
                h[ix(l, l - 1, ldh)] = 0.0;
            }
            // Exit from loop if a submatrix of order 1 or 2 has split off.
            if l + 1 >= iu {
                done = Some(its);
                break;
            }
            // Now the active submatrix is in rows and columns L to I.
            if !wantt {
                i1 = l;
                i2 = iu;
            }
            let (h44, h33, h43h34);
            if its == 10 || its == 20 {
                // Exceptional shift.
                let s = h[ix(iu, iu - 1, ldh)].abs() + h[ix(iu - 1, iu - 2, ldh)].abs();
                h44 = DAT1 * s;
                h33 = h44;
                h43h34 = DAT2 * s * s;
            } else {
                // Prepare to use Wilkinson's double shift.
                h44 = h[ix(iu, iu, ldh)];
                h33 = h[ix(iu - 1, iu - 1, ldh)];
                h43h34 = h[ix(iu, iu - 1, ldh)] * h[ix(iu - 1, iu, ldh)];
            }
            // Look for two consecutive small subdiagonal elements.
            let mut m = iu - 2;
            loop {
                let h11 = h[ix(m, m, ldh)];
                let h22 = h[ix(m + 1, m + 1, ldh)];
                let h21 = h[ix(m + 1, m, ldh)];
                let h12 = h[ix(m, m + 1, ldh)];
                let h44s = h44 - h11;
                let h33s = h33 - h11;
                let mut v1 = (h33s * h44s - h43h34) / h21 + h12;
                let mut v2 = h22 - h11 - h33s - h44s;
                let mut v3 = h[ix(m + 2, m + 1, ldh)];
                let s = v1.abs() + v2.abs() + v3.abs();
                v1 /= s;
                v2 /= s;
                v3 /= s;
                v = [v1, v2, v3];
                if m == l {
                    break;
                }
                let h00 = h[ix(m - 1, m - 1, ldh)];
                let h10 = h[ix(m, m - 1, ldh)];
                let tst1 = v1.abs() * (h00.abs() + h11.abs() + h22.abs());
                if h10.abs() * (v2.abs() + v3.abs()) <= ulp * tst1 {
                    break;
                }
                m -= 1;
            }
            // Double-shift QR step.
            for k in m..iu {
                let nr = 3.min(iu - k + 1);
                if k > m {
                    for q in 0..nr {
                        v[q] = h[ix(k + q, k - 1, ldh)];
                    }
                }
                let (beta, t1) = dlarfg(nr, v[0], &mut v[1..nr]);
                v[0] = beta;
                if k > m {
                    h[ix(k, k - 1, ldh)] = v[0];
                    h[ix(k + 1, k - 1, ldh)] = 0.0;
                    if k + 1 < iu {
                        h[ix(k + 2, k - 1, ldh)] = 0.0;
                    }
                } else if m > l {
                    h[ix(k, k - 1, ldh)] = -h[ix(k, k - 1, ldh)];
                }
                let v2 = v[1];
                let t2 = t1 * v2;
                if nr == 3 {
                    let v3 = v[2];
                    let t3 = t1 * v3;
                    for j in k..=i2 {
                        let sum = h[ix(k, j, ldh)] + v2 * h[ix(k + 1, j, ldh)] + v3 * h[ix(k + 2, j, ldh)];
                        h[ix(k, j, ldh)] -= sum * t1;
                        h[ix(k + 1, j, ldh)] -= sum * t2;
                        h[ix(k + 2, j, ldh)] -= sum * t3;
                    }
                    for j in i1..=(k + 3).min(iu) {
                        let sum = h[ix(j, k, ldh)] + v2 * h[ix(j, k + 1, ldh)] + v3 * h[ix(j, k + 2, ldh)];
                        h[ix(j, k, ldh)] -= sum * t1;
                        h[ix(j, k + 1, ldh)] -= sum * t2;
                        h[ix(j, k + 2, ldh)] -= sum * t3;
                    }
                    // Accumulate transformations for Z.
                    let sum = z[k - 1] + v2 * z[k] + v3 * z[k + 1];
                    z[k - 1] -= sum * t1;
                    z[k] -= sum * t2;
                    z[k + 1] -= sum * t3;
                } else if nr == 2 {
                    for j in k..=i2 {
                        let sum = h[ix(k, j, ldh)] + v2 * h[ix(k + 1, j, ldh)];
                        h[ix(k, j, ldh)] -= sum * t1;
                        h[ix(k + 1, j, ldh)] -= sum * t2;
                    }
                    for j in i1..=iu {
                        let sum = h[ix(j, k, ldh)] + v2 * h[ix(j, k + 1, ldh)];
                        h[ix(j, k, ldh)] -= sum * t1;
                        h[ix(j, k + 1, ldh)] -= sum * t2;
                    }
                    let sum = z[k - 1] + v2 * z[k];
                    z[k - 1] -= sum * t1;
                    z[k] -= sum * t2;
                }
            }
        }
        let Some(its) = done else {
            // Failure to converge in remaining number of iterations.
            return iu;
        };
        // 140
        if l == iu {
            // H(I,I-1) is negligible: one eigenvalue has converged.
            wr[iu - 1] = h[ix(iu, iu, ldh)];
            wi[iu - 1] = 0.0;
        } else if l + 1 == iu {
            // H(I-1,I-2) is negligible; a pair of eigenvalues have
            // converged. Transform the 2-by-2 submatrix to standard Schur
            // form, and compute and store the eigenvalues.
            let (a, b, c, d, r1, i1v, r2, i2v, cs, sn) =
                dlanv2(h[ix(iu - 1, iu - 1, ldh)], h[ix(iu - 1, iu, ldh)], h[ix(iu, iu - 1, ldh)], h[ix(iu, iu, ldh)]);
            h[ix(iu - 1, iu - 1, ldh)] = a;
            h[ix(iu - 1, iu, ldh)] = b;
            h[ix(iu, iu - 1, ldh)] = c;
            h[ix(iu, iu, ldh)] = d;
            wr[iu - 2] = r1;
            wi[iu - 2] = i1v;
            wr[iu - 1] = r2;
            wi[iu - 1] = i2v;
            if wantt {
                // Apply the transformation to the rest of H and to Z, as
                // required.
                if i2 > iu {
                    drot_within(i2 - iu, h, ix(iu - 1, iu + 1, ldh), ldh, ix(iu, iu + 1, ldh), ldh, cs, sn);
                }
                drot_within(iu - i1 - 1, h, ix(i1, iu - 1, ldh), 1, ix(i1, iu, ldh), 1, cs, sn);
                let sum = cs * z[iu - 2] + sn * z[iu - 1];
                z[iu - 1] = cs * z[iu - 1] - sn * z[iu - 2];
                z[iu - 2] = sum;
            }
        }
        // Decrement number of remaining iterations, and return to start of
        // the main loop with new value of I.
        itn -= its;
        i = l as isize - 1;
    }
}

/// `DNEIGH(rnorm, n, h, ldh, ritzr, ritzi, bounds, q, ldq, workl, ierr)`:
/// the Ritz values of the Hessenberg `h` and their error bounds. `q`
/// receives the eigenvectors of the Schur form (`ldq` x `n`), `workl` has
/// `n*(n+3)` entries. Returns `ierr`.
#[allow(clippy::too_many_arguments)]
pub fn dneigh(
    rnorm: f64,
    n: usize,
    h: &[f64],
    ldh: usize,
    ritzr: &mut [f64],
    ritzi: &mut [f64],
    bounds: &mut [f64],
    q: &mut [f64],
    ldq: usize,
    workl: &mut [f64],
) -> usize {
    // 1. Compute the eigenvalues, the last components of the corresponding
    // Schur vectors and the full Schur form T of the current upper
    // Hessenberg matrix H. dlaqrb returns the full Schur form of H in WORKL
    // and the last components of the Schur vectors in BOUNDS.
    dlacpy(n, n, h, ldh, workl, n);
    let ierr = dlaqrb(true, n, 1, n, workl, n, ritzr, ritzi, bounds);
    if ierr != 0 {
        return ierr;
    }
    // 2. Compute the eigenvectors of the full Schur form T and apply the
    // last components of the Schur vectors to get the last components of
    // the corresponding eigenvectors.
    {
        let (t, work) = workl.split_at_mut(n * n);
        dtrevc_right(None, n, t, n, q, ldq, work);
    }
    // Scale the returning eigenvectors so that their Euclidean norms are
    // all one. A complex conjugate pair shares one scale.
    let mut iconj = false;
    for i in 1..=n {
        let c = ix(1, i, ldq);
        if ritzi[i - 1].abs() <= 0.0 {
            let temp = dnrm2(n, &q[c..]);
            dscal(n, 1.0 / temp, &mut q[c..]);
        } else if !iconj {
            let c2 = ix(1, i + 1, ldq);
            let temp = dlapy2(dnrm2(n, &q[c..]), dnrm2(n, &q[c2..]));
            dscal(n, 1.0 / temp, &mut q[c..]);
            dscal(n, 1.0 / temp, &mut q[c2..]);
            iconj = true;
        } else {
            iconj = false;
        }
    }
    dgemv_t(n, n, 1.0, q, ldq, bounds, 0.0, workl);
    // 3. Compute the Ritz estimates.
    let mut iconj = false;
    for i in 0..n {
        if ritzi[i].abs() <= 0.0 {
            bounds[i] = rnorm * workl[i].abs();
        } else if !iconj {
            // The Ritz estimate of a complex pair is the same for both.
            bounds[i] = rnorm * dlapy2(workl[i], workl[i + 1]);
            bounds[i + 1] = bounds[i];
            iconj = true;
        } else {
            iconj = false;
        }
    }
    0
}

/// `smlnum = unfl*(n/ulp)` as `dnaitr` and `dnapps` `SAVE` it on their
/// first call.
fn smlnum_once(cell: &mut Option<f64>, n: usize) -> f64 {
    *cell.get_or_insert(SAFMIN * (n as f64 / PREC))
}

/// `DNAITR(ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, ipntr,
/// workd, info)`: extend a `k`-step Arnoldi factorisation by `np` steps.
/// Returns `info` (0, or the step at which no new starting vector could be
/// found).
#[allow(clippy::too_many_arguments)]
fn dnaitr(
    par: Par<'_>,
    op: &mut Op<'_>,
    st: &mut Stats,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    k: usize,
    np: usize,
    resid: &mut [f64],
    rnorm: &mut f64,
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    workd: &mut [f64],
) -> Result<usize, String> {
    let smlnum = smlnum_once(&mut arp.naitr_smlnum, n);
    let unfl = SAFMIN;
    let ulp = PREC;
    let (ipj, irj, ivj) = (0usize, n, 2 * n);
    let mut j = k + 1;
    loop {
        // 1000
        let mut betaj = *rnorm;
        // `IF (rnorm .GT. zero) GO TO 40`: a NaN norm restarts too.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(*rnorm > 0.0) {
            // Invariant subspace found: generate a new starting vector
            // that is orthogonal to the current basis.
            betaj = 0.0;
            st.nrstrt += 1;
            let mut itry = 1;
            loop {
                // 20 / 30
                let (rn, ierr) = dgetv0(par, op, st, arp, bmat_g, false, n, j, v, ldv, resid, workd)?;
                *rnorm = rn;
                if ierr < 0 {
                    itry += 1;
                    if itry <= 3 {
                        continue;
                    }
                    // Give up after several restarts.
                    return Ok(j - 1);
                }
                break;
            }
        }
        // 40: STEP 2. v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm.
        let vj = (j - 1) * ldv;
        dcopy(n, resid, &mut v[vj..]);
        if *rnorm >= unfl {
            let temp1 = 1.0 / *rnorm;
            dscal(n, temp1, &mut v[vj..]);
            dscal(n, temp1, &mut workd[ipj..]);
        } else {
            dlascl_g(*rnorm, 1.0, n, 1, &mut v[vj..], n);
            dlascl_g(*rnorm, 1.0, n, 1, &mut workd[ipj..], n);
        }
        // STEP 3. r_{j} = OP*v_{j}.
        st.nopx += 1;
        dcopy(n, &v[vj..], &mut workd[ivj..]);
        op(1, [ivj, irj, ipj], workd)?;
        // 50: WORKD(IRJ:IRJ+N-1) := OP*v_{j}; if bmat = 'G', compute B*OP*v_{j}.
        dcopy(n, &workd[irj..], resid);
        if bmat_g {
            st.nbx += 1;
            op(2, [irj, ipj, 0], workd)?;
        } else {
            dcopy(n, resid, &mut workd[ipj..]);
        }
        // 60
        let wnorm = if bmat_g { pddot(par, n, resid, &workd[ipj..]).abs().sqrt() } else { pdnorm2(par, n, resid) };
        // STEP 5: h = V^T B r; r = r - V h.
        let hj = (j - 1) * ldh;
        dgemv_t(n, j, 1.0, v, ldv, &workd[ipj..], 0.0, &mut h[hj..]);
        // PARPACK: computed in its workspace, reduced into h(1,j).
        psum(par, &mut h[hj..hj + j]);
        {
            let col: Vec<f64> = h[hj..hj + j].to_vec();
            dgemv_n(n, j, -1.0, v, ldv, &col, 1.0, resid);
        }
        if j > 1 {
            h[ix(j, j - 1, ldh)] = betaj;
        }
        if bmat_g {
            st.nbx += 1;
            dcopy(n, resid, &mut workd[irj..]);
            op(2, [irj, ipj, 0], workd)?;
        } else {
            dcopy(n, resid, &mut workd[ipj..]);
        }
        // 70: compute the B-norm of r_{j}.
        *rnorm = if bmat_g { pddot(par, n, resid, &workd[ipj..]).abs().sqrt() } else { pdnorm2(par, n, resid) };
        // STEP 5: Re-orthogonalization / Iterative refinement phase.
        // `IF (rnorm .GT. 0.717*wnorm) GO TO 100`: a NaN norm refines too.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(*rnorm > c0717() * wnorm) {
            let mut iter = 0;
            st.nrorth += 1;
            loop {
                // 80
                if par.is_some() {
                    // PARPACK: V^T B r reduced in its workspace; workd(irj)
                    // untouched.
                    let mut t = vec![0.0; j];
                    dgemv_t(n, j, 1.0, v, ldv, &workd[ipj..], 0.0, &mut t);
                    psum(par, &mut t);
                    dgemv_n(n, j, -1.0, v, ldv, &t, 1.0, resid);
                    daxpy(j, 1.0, &t, &mut h[hj..]);
                } else {
                    {
                        let (p, r) = workd.split_at_mut(irj);
                        dgemv_t(n, j, 1.0, v, ldv, p, 0.0, r);
                    }
                    dgemv_n(n, j, -1.0, v, ldv, &workd[irj..], 1.0, resid);
                    daxpy(j, 1.0, &workd[irj..], &mut h[hj..]);
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
                // Otherwise RESID is numerically in the span of V.
                resid[..n].iter_mut().for_each(|x| *x = 0.0);
                *rnorm = 0.0;
                break;
            }
        }
        // 100
        j += 1;
        if j > k + np {
            // Check for splitting and deflation (DLAHQR's test).
            for i in k.max(1)..k + np {
                let mut tst1 = h[ix(i, i, ldh)].abs() + h[ix(i + 1, i + 1, ldh)].abs();
                if tst1 == 0.0 {
                    tst1 = dlanhs_one(k + np, h, ldh);
                }
                if h[ix(i + 1, i, ldh)].abs() <= (ulp * tst1).max(smlnum) {
                    h[ix(i + 1, i, ldh)] = 0.0;
                }
            }
            return Ok(0);
        }
    }
}

/// `DNAPPS(n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq,
/// workl, workd)`: apply the `np` shifts implicitly, keeping `kev` vectors.
/// A complex shift with its conjugate left out moves the split by one
/// (`kev` grows).
#[allow(clippy::too_many_arguments)]
fn dnapps(
    par: Par<'_>,
    arp: &mut Arpack,
    n: usize,
    kev: &mut usize,
    np: usize,
    shiftr: &[f64],
    shifti: &[f64],
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    resid: &mut [f64],
    q: &mut [f64],
    ldq: usize,
    workl: &mut [f64],
    workd: &mut [f64],
) {
    let smlnum = smlnum_once(&mut arp.napps_smlnum, n);
    let ulp = PREC;
    let kplusp = *kev + np;
    // Initialize Q to the identity to accumulate the rotations and
    // reflections.
    dlaset(kplusp, kplusp, 0.0, 1.0, q, ldq);
    if np == 0 {
        return;
    }
    let mut cconj = false;
    for jj in 1..=np {
        let sigmar = shiftr[jj - 1];
        let sigmai = shifti[jj - 1];
        if cconj {
            // The conjugate of the previous shift: already applied.
            cconj = false;
            continue;
        } else if jj < np && sigmai.abs() > 0.0 {
            cconj = true;
        } else if jj == np && sigmai.abs() > 0.0 {
            // The last shift has a nonzero imaginary part. Don't apply it;
            // thus the order of the compressed H is order KEV+1 since only
            // np-1 were applied.
            *kev += 1;
            continue;
        }
        let mut istart = 1;
        loop {
            // 20: Check for splitting and deflation (DLAHQR's test).
            let mut iend = kplusp;
            for i in istart..kplusp {
                let mut tst1 = h[ix(i, i, ldh)].abs() + h[ix(i + 1, i + 1, ldh)].abs();
                if tst1 == 0.0 {
                    tst1 = dlanhs_one(kplusp - jj + 1, h, ldh);
                }
                if h[ix(i + 1, i, ldh)].abs() <= (ulp * tst1).max(smlnum) {
                    iend = i;
                    h[ix(i + 1, i, ldh)] = 0.0;
                    break;
                }
            }
            // 40: no reason to apply a shift to a block of order 1, nor a
            // complex shift to a block of order 2.
            let skip = istart == iend || (istart + 1 == iend && sigmai.abs() > 0.0);
            if !skip {
                let h11 = h[ix(istart, istart, ldh)];
                let h21 = h[ix(istart + 1, istart, ldh)];
                if sigmai.abs() <= 0.0 {
                    // Real-valued shift ==> apply single shift QR.
                    let mut f = h11 - sigmar;
                    let mut g = h21;
                    for i in istart..iend {
                        // Construct the plane rotation G to zero out the
                        // bulge.
                        let (mut c, mut s, mut r) = dlartg(f, g);
                        if i > istart {
                            // The following ensures that h(1:iend-1,1), the
                            // first iend-2 off diagonal of elements H,
                            // remain non negative.
                            if r < 0.0 {
                                r = -r;
                                c = -c;
                                s = -s;
                            }
                            h[ix(i, i - 1, ldh)] = r;
                            h[ix(i + 1, i - 1, ldh)] = 0.0;
                        }
                        // Apply rotation to the left of H; H <- G'*H.
                        for j in i..=kplusp {
                            let t = c * h[ix(i, j, ldh)] + s * h[ix(i + 1, j, ldh)];
                            h[ix(i + 1, j, ldh)] = -s * h[ix(i, j, ldh)] + c * h[ix(i + 1, j, ldh)];
                            h[ix(i, j, ldh)] = t;
                        }
                        // Apply rotation to the right of H; H <- H*G.
                        let jend = (i + 2).min(iend);
                        for j in 1..=jend {
                            let t = c * h[ix(j, i, ldh)] + s * h[ix(j, i + 1, ldh)];
                            h[ix(j, i + 1, ldh)] = -s * h[ix(j, i, ldh)] + c * h[ix(j, i + 1, ldh)];
                            h[ix(j, i, ldh)] = t;
                        }
                        // Accumulate the rotation in the matrix Q; Q <- Q*G.
                        // The bound reads loop 60's index after it ended
                        // (`DO 70 J = 1, MIN(J+JJ, KPLUSP)`): jend + 1.
                        // PARPACK's pdnapps bounds it by min(i+jj, kplusp).
                        let ub = if par.is_some() { (i + jj).min(kplusp) } else { (jend + 1 + jj).min(kplusp) };
                        for j in 1..=ub {
                            let t = c * q[ix(j, i, ldq)] + s * q[ix(j, i + 1, ldq)];
                            q[ix(j, i + 1, ldq)] = -s * q[ix(j, i, ldq)] + c * q[ix(j, i + 1, ldq)];
                            q[ix(j, i, ldq)] = t;
                        }
                        // Prepare for next rotation.
                        if i + 1 < iend {
                            f = h[ix(i + 1, i, ldh)];
                            g = h[ix(i + 2, i, ldh)];
                        }
                    }
                } else {
                    // Complex conjugate shifts ==> apply double shift QR.
                    let h12 = h[ix(istart, istart + 1, ldh)];
                    let h22 = h[ix(istart + 1, istart + 1, ldh)];
                    let h32 = h[ix(istart + 2, istart + 1, ldh)];
                    // Compute 1st column of (H - shift*I)*(H - conj(shift)*I).
                    let s = 2.0 * sigmar;
                    let t = dlapy2(sigmar, sigmai);
                    let mut u = [(h11 * (h11 - s) + t * t) / h21 + h12, h11 + h22 - s, h32];
                    for i in istart..iend {
                        let nr = 3.min(iend - i + 1);
                        // Construct Householder reflector G to zero out u(1).
                        // G is of the form I - tau*( 1 u )' * ( 1 u' ).
                        let (beta, tau) = dlarfg(nr, u[0], &mut u[1..nr]);
                        u[0] = beta;
                        if i > istart {
                            h[ix(i, i - 1, ldh)] = u[0];
                            h[ix(i + 1, i - 1, ldh)] = 0.0;
                            if i + 1 < iend {
                                h[ix(i + 2, i - 1, ldh)] = 0.0;
                            }
                        }
                        u[0] = 1.0;
                        // Apply the reflector to the left of H.
                        let hc = ix(i, i, ldh);
                        dlarf(true, nr, kplusp - i + 1, &u[..nr], tau, &mut h[hc..], ldh, workl);
                        // Apply the reflector to the right of H.
                        let ir = (i + 3).min(iend);
                        let hc = ix(1, i, ldh);
                        dlarf(false, ir, nr, &u[..nr], tau, &mut h[hc..], ldh, workl);
                        // Accumulate the reflector in the matrix Q; Q <- Q*G.
                        let qc = ix(1, i, ldq);
                        dlarf(false, kplusp, nr, &u[..nr], tau, &mut q[qc..], ldq, workl);
                        // Prepare for next reflector.
                        if i + 1 < iend {
                            u[0] = h[ix(i + 1, i, ldh)];
                            u[1] = h[ix(i + 2, i, ldh)];
                            if i + 2 < iend {
                                u[2] = h[ix(i + 3, i, ldh)];
                            }
                        }
                    }
                }
            }
            // 100: Apply the same shift to the next block if there is any.
            istart = iend + 1;
            if iend >= kplusp {
                break;
            }
        }
    }
    let kev_ = *kev;
    // Perform a similarity transformation that makes the sub-diagonal
    // elements of H non negative.
    for j in 1..=kev_ {
        if h[ix(j + 1, j, ldh)] < 0.0 {
            let c = ix(j + 1, j, ldh);
            for t in 0..kplusp - j + 1 {
                h[c + t * ldh] = -h[c + t * ldh];
            }
            dscal((j + 2).min(kplusp), -1.0, &mut h[ix(1, j + 1, ldh)..]);
            dscal((j + np + 1).min(kplusp), -1.0, &mut q[ix(1, j + 1, ldq)..]);
        }
    }
    for i in 1..=kev_ {
        // Final check for splitting and deflation (no ABS on H(I+1,I)).
        let mut tst1 = h[ix(i, i, ldh)].abs() + h[ix(i + 1, i + 1, ldh)].abs();
        if tst1 == 0.0 {
            tst1 = dlanhs_one(kev_, h, ldh);
        }
        if h[ix(i + 1, i, ldh)] <= (ulp * tst1).max(smlnum) {
            h[ix(i + 1, i, ldh)] = 0.0;
        }
    }
    // Compute the (kev+1)-st column of (V*Q) and temporarily store the
    // result in WORKD(N+1:2*N). This is needed in the residual update since
    // we cannot GUARANTEE that the corresponding entry of H would be zero
    // as in exact arithmetic.
    let hk = h[ix(kev_ + 1, kev_, ldh)];
    if hk > 0.0 {
        let (lo, hi) = workd.split_at_mut(n);
        let _ = lo;
        dgemv_n(n, kplusp, 1.0, v, ldv, &q[ix(1, kev_ + 1, ldq)..], 0.0, hi);
    }
    // Compute column 1 to kev of (V*Q) in backward order taking advantage
    // of the upper Hessenberg structure of Q.
    for i in 1..=kev_ {
        dgemv_n(n, kplusp - i + 1, 1.0, v, ldv, &q[ix(1, kev_ - i + 1, ldq)..], 0.0, workd);
        let c = ix(1, kplusp - i + 1, ldv);
        dcopy(n, workd, &mut v[c..]);
    }
    // Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).
    {
        let src: Vec<f64> = v[ix(1, kplusp - kev_ + 1, ldv)..ix(1, kplusp - kev_ + 1, ldv) + ldv * (kev_ - 1) + n].to_vec();
        dlacpy(n, kev_, &src, ldv, v, ldv);
    }
    // Copy the (kev+1)-st column of (V*Q) in the appropriate place.
    if hk > 0.0 {
        let c = ix(1, kev_ + 1, ldv);
        dcopy(n, &workd[n..], &mut v[c..]);
    }
    // Update the residual vector:
    // r <- sigmak*r + betak*v(:,kev+1)
    // where sigmak = (e_{kplusp}'*Q)*e_{kev}, betak = e_{kev+1}'*H*e_{kev}.
    dscal(n, q[ix(kplusp, kev_, ldq)], resid);
    if hk > 0.0 {
        daxpy(n, hk, &v[ix(1, kev_ + 1, ldv)..], resid);
    }
}

/// The parameters and results of [`dnaupd`] (`iparam`, `ipntr` and `info`).
#[derive(Clone, Debug, Default)]
pub struct Naupd {
    /// `IPARAM(1..11)` (index 0 is `IPARAM(1)`).
    pub iparam: [i32; 11],
    /// `IPNTR(1..14)`, 1-based offsets into `workl` (as the Fortran).
    pub ipntr: [usize; 14],
    pub info: i32,
    pub stats: Stats,
}

/// `DNAUPD(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
/// ipntr, workd, workl, lworkl, info)` driven to completion. `iparam`
/// carries `ishift`, `mxiter` and `mode` in (`IPARAM(1)`, `(3)`, `(7)`);
/// `info_in` is 0 for a random starting vector, nonzero to use `resid`.
/// `workd` has `3 n`, `workl` at least `3 ncv**2 + 6 ncv` entries.
#[allow(clippy::too_many_arguments)]
pub fn dnaupd(
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
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
) -> Result<Naupd, String> {
    naupd(None, op, arp, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, iparam, workd, workl, info_in)
}

/// PARPACK's `PDNAUPD(comm, ...)`: [`dnaupd`] on this rank's `n` rows of
/// the vectors, the reductions over `comm` (see [`super::par`]). Collective.
#[allow(clippy::too_many_arguments)]
pub fn pdnaupd(
    comm: &dyn elmer_core::Communicator,
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
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
) -> Result<Naupd, String> {
    naupd(Some(comm), op, arp, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, iparam, workd, workl, info_in)
}

#[allow(clippy::too_many_arguments)]
fn naupd(
    par: Par<'_>,
    op: &mut Op<'_>,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
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
) -> Result<Naupd, String> {
    let mut out = Naupd { iparam, ..Default::default() };
    let ishift = iparam[0];
    let mut mxiter = iparam[2];
    let mode = iparam[6];
    let lworkl = 3 * ncv * ncv + 6 * ncv;
    let ierr = if n == 0 {
        -1
    } else if nev == 0 {
        -2
    } else if ncv <= nev + 1 || (ncv > n && par.is_none()) {
        -3
    } else if mxiter <= 0 {
        -4
    } else if workl.len() < lworkl {
        -7
    } else if !(1..=if par.is_some() { 4 } else { 5 }).contains(&mode) {
        -10
    } else if mode == 1 && bmat_g {
        -11
    } else if !(0..=1).contains(&ishift) {
        -12
    } else {
        0
    };
    if ierr != 0 {
        out.info = ierr;
        return Ok(out);
    }
    if ishift == 0 {
        return Err("DNAUPD with user shifts (ishift = 0) is not ported".into());
    }
    if *tol <= 0.0 {
        *tol = EPS;
    }
    let np = ncv - nev;
    workl[..lworkl].iter_mut().for_each(|x| *x = 0.0);
    let (ldh, ldq) = (ncv, ncv);
    // 1-based, as IPNTR reports them.
    let ih = 1;
    let ritzr = ih + ldh * ncv;
    let ritzi = ritzr + ncv;
    let bounds = ritzi + ncv;
    let iq = bounds + ncv;
    let iw = iq + ldq * ncv;
    let next = iw + ncv * ncv + 3 * ncv;
    out.ipntr[3] = next;
    out.ipntr[4] = ih;
    out.ipntr[5] = ritzr;
    out.ipntr[6] = ritzi;
    out.ipntr[7] = bounds;
    out.ipntr[13] = iw;

    let mut st = Stats::default();
    let (h, rest) = workl.split_at_mut(ldh * ncv);
    let (rr, rest) = rest.split_at_mut(ncv);
    let (ri, rest) = rest.split_at_mut(ncv);
    let (bd, rest) = rest.split_at_mut(ncv);
    let (q, rest) = rest.split_at_mut(ldq * ncv);
    let wl = &mut rest[..ncv * ncv + 3 * ncv];
    let mut nev0 = nev;
    let mut np_ = np;
    let info = dnaup2(
        par, op, &mut st, arp, bmat_g, n, which, &mut nev0, &mut np_, *tol, resid, ishift, &mut mxiter, v, ldv, h, ldh, rr,
        ri, bd, q, ldq, wl, workd, info_in,
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

/// `DNAUP2`: the main loop of the implicitly restarted Arnoldi method.
/// Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dnaup2(
    par: Par<'_>,
    op: &mut Op<'_>,
    st: &mut Stats,
    arp: &mut Arpack,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: &mut usize,
    np: &mut usize,
    tol: f64,
    resid: &mut [f64],
    ishift: i32,
    mxiter: &mut i32,
    v: &mut [f64],
    ldv: usize,
    h: &mut [f64],
    ldh: usize,
    ritzr: &mut [f64],
    ritzi: &mut [f64],
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
    let kplusp = *nev + *np;
    let mut nconv;
    let mut numcnv;
    let mut iter = 0i32;
    let initv = info_in != 0;
    let mut info = 0i32;

    // Get a possibly random starting vector and force it into the range
    // of the operator OP.
    let (mut rnorm, _) = dgetv0(par, op, st, arp, bmat_g, initv, n, 1, v, ldv, resid, workd)?;
    if rnorm == 0.0 {
        // The initial vector is zero. Error exit.
        return Ok(-9);
    }
    // Compute the first NEV steps of the Arnoldi factorization.
    let r = dnaitr(par, op, st, arp, bmat_g, n, 0, *nev, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
    if r > 0 {
        *np = r;
        *mxiter = iter;
        return Ok(-9999);
    }
    loop {
        // 1000
        iter += 1;
        // Compute NP additional steps of the Arnoldi factorization. Adjust
        // NP since NEV might have been updated by last call to dnapps.
        *np = kplusp - *nev;
        // 20
        let r = dnaitr(par, op, st, arp, bmat_g, n, *nev, *np, resid, &mut rnorm, v, ldv, h, ldh, workd)?;
        if r > 0 {
            *np = r;
            *mxiter = iter;
            return Ok(-9999);
        }
        // Compute the eigenvalues and corresponding error bounds of the
        // current upper Hessenberg matrix.
        if dneigh(rnorm, kplusp, h, ldh, ritzr, ritzi, bounds, q, ldq, workl) != 0 {
            return Ok(-8);
        }
        // Make a copy of eigenvalues and corresponding error bounds
        // obtained from dneigh.
        let k2 = kplusp * kplusp;
        dcopy(kplusp, ritzr, &mut workl[k2..]);
        dcopy(kplusp, ritzi, &mut workl[k2 + kplusp..]);
        dcopy(kplusp, bounds, &mut workl[k2 + 2 * kplusp..]);
        // Select the wanted Ritz values and their bounds to be used in the
        // convergence test. The wanted part of the spectrum and
        // corresponding error bounds are in the last NEV loc. of RITZR,
        // RITZI and BOUNDS respectively.
        *nev = nev0;
        *np = np0;
        numcnv = *nev;
        dngets(ishift, which, nev, np, ritzr, ritzi, bounds)?;
        if *nev == nev0 + 1 {
            numcnv = nev0 + 1;
        }
        // Convergence test.
        {
            let (p, k) = (*np, *nev);
            let b: Vec<f64> = bounds[p..p + k].to_vec();
            workl[2 * p..2 * p + k].copy_from_slice(&b);
            nconv = dnconv(k, &ritzr[p..], &ritzi[p..], &workl[2 * p..], tol);
        }
        // Count the number of unwanted Ritz values that have zero Ritz
        // estimates. If any Ritz estimates are equal to zero, then a
        // leading block of H of order equal to at least the number of Ritz
        // values with zero Ritz estimates has split off.
        let nptemp = *np;
        for &b in &bounds[..nptemp] {
            if b == 0.0 {
                *np -= 1;
                *nev += 1;
            }
        }
        if nconv >= numcnv || iter > *mxiter || *np == 0 {
            // Prepare to exit. Put the converged Ritz values and
            // corresponding bounds in RITZ(1:NCONV) and BOUNDS(1:NCONV)
            // respectively. Then sort. Be careful when NCONV > NP.
            //
            // Use h( 3,1 ) as storage to communicate rnorm to _neupd if
            // needed.
            h[2] = rnorm;
            // To be consistent with dngets, we first do a pre-processing
            // sort in order to keep complex conjugate pairs together.
            let wprime = match which {
                NWhich::LM => NWhich::SR,
                NWhich::SM => NWhich::LR,
                NWhich::LR => NWhich::SM,
                NWhich::SR => NWhich::LM,
                NWhich::LI => NWhich::SM,
                NWhich::SI => NWhich::LM,
            };
            dsortc(wprime, true, kplusp, ritzr, ritzi, bounds);
            // Now sort Ritz values so that converged Ritz values appear
            // within the first NEV locations of ritzr, ritzi and bounds, and
            // the most desired one appears at the front.
            let wprime = match which {
                NWhich::LM => NWhich::SM,
                NWhich::SM => NWhich::LM,
                NWhich::LR => NWhich::SR,
                NWhich::SR => NWhich::LR,
                NWhich::LI => NWhich::SI,
                NWhich::SI => NWhich::LI,
            };
            dsortc(wprime, true, kplusp, ritzr, ritzi, bounds);
            // Scale the Ritz estimate of each Ritz value by
            // 1 / max(eps23, magnitude of the Ritz value).
            // (PARPACK's pdnaup2 scales and sorts numcnv of them.)
            let nscale = if par.is_some() { numcnv } else { nev0 };
            for j in 0..nscale {
                let temp = eps23.max(dlapy2(ritzr[j], ritzi[j]));
                bounds[j] /= temp;
            }
            // Sort the Ritz values according to the scaled Ritz estimates.
            // This will push all the converged ones towards the front of
            // ritzr, ritzi, bounds (in the case when NCONV < NEV.)
            dsortc(NWhich::LR, true, nscale, bounds, ritzr, ritzi);
            // Scale the Ritz estimate back to its original value.
            for j in 0..nscale {
                let temp = eps23.max(dlapy2(ritzr[j], ritzi[j]));
                bounds[j] *= temp;
            }
            // Sort the converged Ritz values again so that the "threshold"
            // value appears at the front of ritzr and ritzi.
            dsortc(which, true, nconv, ritzr, ritzi, bounds);
            if iter > *mxiter && nconv < numcnv {
                // Max iterations have been exceeded.
                info = 1;
            }
            if *np == 0 && nconv < numcnv {
                // No shifts to apply.
                info = 2;
            }
            *np = nconv;
            // 1100
            *mxiter = iter;
            *nev = numcnv;
            return Ok(info);
        } else if nconv < numcnv && ishift == 1 {
            // Do not have all the requested eigenvalues yet. To prevent
            // possible stagnation, adjust the size of NEV.
            let nevbef = *nev;
            *nev += nconv.min(*np / 2);
            if *nev == 1 && kplusp >= 6 {
                *nev = kplusp / 2;
            } else if *nev == 1 && kplusp > 3 {
                *nev = 2;
            }
            *np = kplusp - *nev;
            // If the size of NEV was just increased resort the eigenvalues.
            if nevbef < *nev {
                dngets(ishift, which, nev, np, ritzr, ritzi, bounds)?;
            }
        }
        // 50: Apply the NP implicit shifts by QR bulge chasing. Each shift
        // is applied to the whole upper Hessenberg matrix H. The first 2*N
        // locations of WORKD are used as workspace.
        let (sr, si): (Vec<f64>, Vec<f64>) = (ritzr[..*np].to_vec(), ritzi[..*np].to_vec());
        dnapps(par, arp, n, nev, *np, &sr, &si, v, ldv, h, ldh, resid, q, ldq, workl, workd);
        // Compute the B-norm of the updated residual. Keep B*RESID in
        // WORKD(1:N) to be used in the first step of the next call to
        // dnaitr.
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
}
