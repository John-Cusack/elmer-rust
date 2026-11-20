//! The reference LAPACK 3.12 routines under ARPACK's non-symmetric path
//! (`dnaupd` / `dneupd` and ARPACK's own `dlaqrb`), as the pristine Elmer
//! build links them from the system `liblapack.so.3`: `DLAHQR`, `DLANV2`,
//! `DLANHS('1')`, `DTREVC` (right eigenvectors, `HOWMNY` = 'A' or 'S'),
//! `DLALN2` and `DLADIV`.
//!
//! Matrices are column-major slices with a leading dimension; the routines
//! keep Fortran's 1-based loop indices and address `A(i,j)` as
//! `a[(i-1) + (j-1)*lda]`, so each line can be read against the reference
//! source. `DLABAD` does nothing for IEEE double precision (LAPACK 3.12)
//! and is left out.

use super::blas::{daxpy, dcopy, dscal};
use super::lapack::{dlapy2, dlarfg, dlartg, EPS, PREC, SAFMIN};

#[inline]
fn ix(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

#[inline]
fn sign(a: f64, b: f64) -> f64 {
    // Fortran `SIGN(A, B)`: |A| with the sign bit of B.
    a.abs().copysign(b)
}

/// `DROT(N, DX, INCX, DY, INCY, C, S)` over two strided sequences of one
/// slice: `x` starts at `xs` with stride `incx`, `y` at `ys` with `incy`.
#[allow(clippy::too_many_arguments)]
pub fn drot_within(n: usize, a: &mut [f64], xs: usize, incx: usize, ys: usize, incy: usize, c: f64, s: f64) {
    for k in 0..n {
        let (xi, yi) = (xs + k * incx, ys + k * incy);
        let dtemp = c * a[xi] + s * a[yi];
        a[yi] = c * a[yi] - s * a[xi];
        a[xi] = dtemp;
    }
}

/// `IDAMAX(N, DX, 1)`, 1-based; 0 when `n == 0`.
pub fn idamax(n: usize, x: &[f64]) -> usize {
    if n == 0 {
        return 0;
    }
    let mut imax = 1;
    let mut dmax = x[0].abs();
    for (i, v) in x.iter().enumerate().take(n).skip(1) {
        if v.abs() > dmax {
            imax = i + 1;
            dmax = v.abs();
        }
    }
    imax
}

/// `DTRMM('Right', 'Upper', 'No transpose', 'Non-unit', M, N, ALPHA, A,
/// LDA, B, LDB)`: `B := alpha * B * A`, `A` upper triangular.
#[allow(clippy::too_many_arguments)]
pub fn dtrmm_run(m: usize, n: usize, alpha: f64, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    if m == 0 || n == 0 {
        return;
    }
    if alpha == 0.0 {
        for j in 1..=n {
            for i in 1..=m {
                b[ix(i, j, ldb)] = 0.0;
            }
        }
        return;
    }
    for j in (1..=n).rev() {
        let temp = alpha * a[ix(j, j, lda)];
        for i in 1..=m {
            b[ix(i, j, ldb)] *= temp;
        }
        for k in 1..j {
            let akj = a[ix(k, j, lda)];
            if akj != 0.0 {
                let temp = alpha * akj;
                for i in 1..=m {
                    b[ix(i, j, ldb)] += temp * b[ix(i, k, ldb)];
                }
            }
        }
    }
}

/// `DLANHS('1', N, A, LDA, WORK)`: the largest column sum of the upper
/// Hessenberg part.
pub fn dlanhs_one(n: usize, a: &[f64], lda: usize) -> f64 {
    let mut value = 0.0f64;
    for j in 1..=n {
        let mut sum = 0.0;
        for i in 1..=n.min(j + 1) {
            sum += a[ix(i, j, lda)].abs();
        }
        if value < sum || sum.is_nan() {
            value = sum;
        }
    }
    value
}

/// `DLADIV2` (LAPACK 3.7+).
fn dladiv2(a: f64, b: f64, c: f64, d: f64, r: f64, t: f64) -> f64 {
    if r != 0.0 {
        let br = b * r;
        if br != 0.0 {
            (a + br) * t
        } else {
            a * t + (b * t) * r
        }
    } else {
        (a + d * (b / c)) * t
    }
}

/// `DLADIV1` (LAPACK 3.7+).
fn dladiv1(a: f64, b: f64, c: f64, d: f64) -> (f64, f64) {
    let r = d / c;
    let t = 1.0 / (c + d * r);
    let p = dladiv2(a, b, c, d, r, t);
    let q = dladiv2(b, -a, c, d, r, t);
    (p, q)
}

/// `DLADIV(A, B, C, D, P, Q)`: `p + iq = (a + ib) / (c + id)`, the
/// robust algorithm of Baudin and Smith that LAPACK 3.7 adopted.
pub fn dladiv(a: f64, b: f64, c: f64, d: f64) -> (f64, f64) {
    const BS: f64 = 2.0;
    let (mut aa, mut bb, mut cc, mut dd) = (a, b, c, d);
    let ab = a.abs().max(b.abs());
    let cd = c.abs().max(d.abs());
    let mut s = 1.0;
    let ov = f64::MAX;
    let un = SAFMIN;
    let eps = EPS;
    let be = BS / (eps * eps);
    if ab >= 0.5 * ov {
        aa *= 0.5;
        bb *= 0.5;
        s *= 2.0;
    }
    if cd >= 0.5 * ov {
        cc *= 0.5;
        dd *= 0.5;
        s *= 0.5;
    }
    if ab <= un * BS / eps {
        aa *= be;
        bb *= be;
        s /= be;
    }
    if cd <= un * BS / eps {
        cc *= be;
        dd *= be;
        s *= be;
    }
    let (p, mut q);
    if d.abs() <= c.abs() {
        (p, q) = dladiv1(aa, bb, cc, dd);
    } else {
        (p, q) = dladiv1(bb, aa, dd, cc);
        q = -q;
    }
    (p * s, q * s)
}

/// The Schur factorization of a real 2x2 nonsymmetric matrix in
/// standardized form (`DLANV2`, LAPACK 3.12): returns
/// `(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)`.
#[allow(clippy::type_complexity)]
pub fn dlanv2(a: f64, b: f64, c: f64, d: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    const MULTPL: f64 = 4.0;
    let (mut a, mut b, mut c, mut d) = (a, b, c, d);
    let safmin = SAFMIN;
    let eps = PREC;
    // `DLAMCH('B')**INT(LOG(SAFMIN/EPS)/LOG(DLAMCH('B'))/TWO)`.
    let k = (((safmin / eps).ln() / 2.0f64.ln()) / 2.0) as i32;
    let safmn2 = 2.0f64.powi(k);
    let safmx2 = 1.0 / safmn2;
    let (mut cs, mut sn);
    if c == 0.0 {
        cs = 1.0;
        sn = 0.0;
    } else if b == 0.0 {
        // Swap rows and columns.
        cs = 0.0;
        sn = 1.0;
        std::mem::swap(&mut a, &mut d);
        b = -c;
        c = 0.0;
    } else if (a - d) == 0.0 && sign(1.0, b) != sign(1.0, c) {
        cs = 1.0;
        sn = 0.0;
    } else {
        let mut temp = a - d;
        let mut p = 0.5 * temp;
        let bcmax = b.abs().max(c.abs());
        let bcmis = b.abs().min(c.abs()) * sign(1.0, b) * sign(1.0, c);
        let mut scale = p.abs().max(bcmax);
        let mut z = (p / scale) * p + (bcmax / scale) * bcmis;
        if z >= MULTPL * eps {
            // Real eigenvalues. Compute A and D.
            z = p + sign(scale.sqrt() * z.sqrt(), p);
            a = d + z;
            d -= (bcmax / z) * bcmis;
            // Compute B and the rotation matrix.
            let tau = dlapy2(c, z);
            cs = z / tau;
            sn = c / tau;
            b -= c;
            c = 0.0;
        } else {
            // Complex eigenvalues, or real (almost) equal eigenvalues.
            // Make diagonal elements equal.
            let mut count = 0;
            let mut sigma = b + c;
            loop {
                count += 1;
                scale = temp.abs().max(sigma.abs());
                if scale >= safmx2 {
                    sigma *= safmn2;
                    temp *= safmn2;
                    if count <= 20 {
                        continue;
                    }
                }
                if scale <= safmn2 {
                    sigma *= safmx2;
                    temp *= safmx2;
                    if count <= 20 {
                        continue;
                    }
                }
                break;
            }
            p = 0.5 * temp;
            let mut tau = dlapy2(sigma, temp);
            cs = (0.5 * (1.0 + sigma.abs() / tau)).sqrt();
            sn = -(p / (tau * cs)) * sign(1.0, sigma);
            // [AA BB; CC DD] = [A B; C D] [CS -SN; SN CS]
            let aa = a * cs + b * sn;
            let bb = -a * sn + b * cs;
            let cc = c * cs + d * sn;
            let dd = -c * sn + d * cs;
            // [A B; C D] = [CS SN; -SN CS] [AA BB; CC DD]
            a = aa * cs + cc * sn;
            b = bb * cs + dd * sn;
            c = -aa * sn + cc * cs;
            d = -bb * sn + dd * cs;
            temp = 0.5 * (a + d);
            a = temp;
            d = temp;
            if c != 0.0 {
                if b != 0.0 {
                    if sign(1.0, b) == sign(1.0, c) {
                        // Real eigenvalues: reduce to upper triangular form.
                        let sab = b.abs().sqrt();
                        let sac = c.abs().sqrt();
                        p = sign(sab * sac, c);
                        tau = 1.0 / (b + c).abs().sqrt();
                        a = temp + p;
                        d = temp - p;
                        b -= c;
                        c = 0.0;
                        let cs1 = sab * tau;
                        let sn1 = sac * tau;
                        let temp = cs * cs1 - sn * sn1;
                        sn = cs * sn1 + sn * cs1;
                        cs = temp;
                    }
                } else {
                    b = -c;
                    c = 0.0;
                    let temp = cs;
                    cs = -sn;
                    sn = temp;
                }
            }
        }
    }
    let (rt1r, rt2r) = (a, d);
    let (rt1i, rt2i) = if c == 0.0 {
        (0.0, 0.0)
    } else {
        let r = b.abs().sqrt() * c.abs().sqrt();
        (r, -r)
    };
    (a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
}

/// `DLAHQR(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ,
/// INFO)` of LAPACK 3.12: the eigenvalues of the Hessenberg submatrix in
/// rows and columns `ilo..=ihi` (1-based) by the double-shift QR algorithm,
/// with the Schur form in `h` when `wantt` and the Schur vectors
/// accumulated into rows `iloz..=ihiz` of `z` when `wantz`. Returns `INFO`.
#[allow(clippy::too_many_arguments)]
pub fn dlahqr(
    wantt: bool,
    wantz: bool,
    n: usize,
    ilo: usize,
    ihi: usize,
    h: &mut [f64],
    ldh: usize,
    wr: &mut [f64],
    wi: &mut [f64],
    iloz: usize,
    ihiz: usize,
    z: &mut [f64],
    ldz: usize,
) -> usize {
    const DAT1: f64 = 3.0 / 4.0;
    const DAT2: f64 = -0.4375;
    const KEXSH: usize = 10;
    if n == 0 {
        return 0;
    }
    if ilo == ihi {
        wr[ilo - 1] = h[ix(ilo, ilo, ldh)];
        wi[ilo - 1] = 0.0;
        return 0;
    }
    // Clear out the trash.
    for j in ilo..=ihi.saturating_sub(3) {
        h[ix(j + 2, j, ldh)] = 0.0;
        h[ix(j + 3, j, ldh)] = 0.0;
    }
    if ilo + 2 <= ihi {
        h[ix(ihi, ihi - 2, ldh)] = 0.0;
    }
    let nh = ihi - ilo + 1;
    let nz = ihiz + 1 - iloz;
    let safmin = SAFMIN;
    let ulp = PREC;
    let smlnum = safmin * (nh as f64 / ulp);
    let (mut i1, mut i2) = (0usize, 0usize);
    if wantt {
        i1 = 1;
        i2 = n;
    }
    let itmax = 30 * nh.max(10);
    let mut kdefl = 0usize;
    // `I` indexes the active block's last row; the blocks above deflate
    // as `I` moves up. Signed: the loop ends when it passes `ILO`.
    let mut i = ihi as isize;
    let mut v = [0.0f64; 3];
    loop {
        // 20
        let mut l = ilo;
        if i < ilo as isize {
            return 0;
        }
        let iu = i as usize;
        let mut converged = false;
        for _its in 0..=itmax {
            // Look for a single small subdiagonal element.
            let mut k = iu;
            while k > l {
                let hkk1 = h[ix(k, k - 1, ldh)].abs();
                if hkk1 <= smlnum {
                    break;
                }
                let mut tst = h[ix(k - 1, k - 1, ldh)].abs() + h[ix(k, k, ldh)].abs();
                if tst == 0.0 {
                    if k >= ilo + 2 {
                        tst += h[ix(k - 1, k - 2, ldh)].abs();
                    }
                    if k < ihi {
                        tst += h[ix(k + 1, k, ldh)].abs();
                    }
                }
                // Ahues & Kahan's conservative small subdiagonal criterion.
                if hkk1 <= ulp * tst {
                    let hk1k = h[ix(k - 1, k, ldh)].abs();
                    let ab = hkk1.max(hk1k);
                    let ba = hkk1.min(hk1k);
                    let hkk = h[ix(k, k, ldh)];
                    let dif = (h[ix(k - 1, k - 1, ldh)] - hkk).abs();
                    let aa = hkk.abs().max(dif);
                    let bb = hkk.abs().min(dif);
                    let s = aa + ab;
                    if ba * (ab / s) <= smlnum.max(ulp * (bb * (aa / s))) {
                        break;
                    }
                }
                k -= 1;
            }
            l = k;
            if l > ilo {
                // H(L,L-1) is negligible.
                h[ix(l, l - 1, ldh)] = 0.0;
            }
            // Exit when a submatrix of order 1 or 2 has split off.
            if l + 1 >= iu {
                converged = true;
                break;
            }
            kdefl += 1;
            // The active submatrix is rows/columns L to I.
            if !wantt {
                i1 = l;
                i2 = iu;
            }
            let (mut h11, mut h12, mut h21, mut h22);
            if kdefl.is_multiple_of(2 * KEXSH) {
                // Exceptional shift.
                let s = h[ix(iu, iu - 1, ldh)].abs() + h[ix(iu - 1, iu - 2, ldh)].abs();
                h11 = DAT1 * s + h[ix(iu, iu, ldh)];
                h12 = DAT2 * s;
                h21 = s;
                h22 = h11;
            } else if kdefl.is_multiple_of(KEXSH) {
                // Exceptional shift.
                let s = h[ix(l + 1, l, ldh)].abs() + h[ix(l + 2, l + 1, ldh)].abs();
                h11 = DAT1 * s + h[ix(l, l, ldh)];
                h12 = DAT2 * s;
                h21 = s;
                h22 = h11;
            } else {
                // Prepare to use Francis' double shift (the eigenvalues of
                // the trailing 2x2).
                h11 = h[ix(iu - 1, iu - 1, ldh)];
                h21 = h[ix(iu, iu - 1, ldh)];
                h12 = h[ix(iu - 1, iu, ldh)];
                h22 = h[ix(iu, iu, ldh)];
            }
            let s = h11.abs() + h12.abs() + h21.abs() + h22.abs();
            let (rt1r, rt1i, rt2r, rt2i);
            if s == 0.0 {
                rt1r = 0.0;
                rt1i = 0.0;
                rt2r = 0.0;
                rt2i = 0.0;
            } else {
                h11 /= s;
                h21 /= s;
                h12 /= s;
                h22 /= s;
                let tr = (h11 + h22) / 2.0;
                let det = (h11 - tr) * (h22 - tr) - h12 * h21;
                let rtdisc = det.abs().sqrt();
                if det >= 0.0 {
                    // Complex conjugate shifts.
                    rt1r = tr * s;
                    rt2r = rt1r;
                    rt1i = rtdisc * s;
                    rt2i = -rt1i;
                } else {
                    // Real shifts (use only one of them).
                    let r1 = tr + rtdisc;
                    let r2 = tr - rtdisc;
                    let r = if (r1 - h22).abs() <= (r2 - h22).abs() { r1 * s } else { r2 * s };
                    rt1r = r;
                    rt2r = r;
                    rt1i = 0.0;
                    rt2i = 0.0;
                }
            }
            // Look for two consecutive small subdiagonals.
            let mut m = iu - 2;
            loop {
                // Determine the effect of starting the double-shift QR
                // iteration at row M, and see if this would make H(M,M-1)
                // negligible.
                let mut h21s = h[ix(m + 1, m, ldh)];
                let hmm = h[ix(m, m, ldh)];
                let s = (hmm - rt2r).abs() + rt2i.abs() + h21s.abs();
                h21s = h[ix(m + 1, m, ldh)] / s;
                v[0] = h21s * h[ix(m, m + 1, ldh)] + (hmm - rt1r) * ((hmm - rt2r) / s) - rt1i * (rt2i / s);
                v[1] = h21s * (hmm + h[ix(m + 1, m + 1, ldh)] - rt1r - rt2r);
                v[2] = h21s * h[ix(m + 2, m + 1, ldh)];
                let s = v[0].abs() + v[1].abs() + v[2].abs();
                v[0] /= s;
                v[1] /= s;
                v[2] /= s;
                if m == l {
                    break;
                }
                let h00 = h[ix(m, m - 1, ldh)].abs() * (v[1].abs() + v[2].abs());
                let h01 = ulp * v[0].abs() * (h[ix(m - 1, m - 1, ldh)].abs() + hmm.abs() + h[ix(m + 1, m + 1, ldh)].abs());
                if h00 <= h01 {
                    break;
                }
                m -= 1;
            }
            // Double-shift QR step.
            for k in m..iu {
                // The first iteration of this loop determines a reflection
                // G from the vector V and applies it from left and right to
                // H, thus creating a nonzero bulge below the subdiagonal.
                // Each subsequent iteration determines a reflection G to
                // restore the Hessenberg form in the (K-1)th column, and
                // thus chases the bulge one step toward the bottom of the
                // active submatrix. NR is the order of G.
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
                    // Use the following instead of H(K,K-1) = -H(K,K-1) to
                    // avoid a bug when v(2) and v(3) underflow.
                    h[ix(k, k - 1, ldh)] *= 1.0 - t1;
                }
                let v2 = v[1];
                let t2 = t1 * v2;
                if nr == 3 {
                    let v3 = v[2];
                    let t3 = t1 * v3;
                    // Apply G from the left to transform the rows of the
                    // matrix in columns K to I2.
                    for j in k..=i2 {
                        let sum = h[ix(k, j, ldh)] + v2 * h[ix(k + 1, j, ldh)] + v3 * h[ix(k + 2, j, ldh)];
                        h[ix(k, j, ldh)] -= sum * t1;
                        h[ix(k + 1, j, ldh)] -= sum * t2;
                        h[ix(k + 2, j, ldh)] -= sum * t3;
                    }
                    // Apply G from the right to transform the columns of
                    // the matrix in rows I1 to min(K+3,I).
                    for j in i1..=(k + 3).min(iu) {
                        let sum = h[ix(j, k, ldh)] + v2 * h[ix(j, k + 1, ldh)] + v3 * h[ix(j, k + 2, ldh)];
                        h[ix(j, k, ldh)] -= sum * t1;
                        h[ix(j, k + 1, ldh)] -= sum * t2;
                        h[ix(j, k + 2, ldh)] -= sum * t3;
                    }
                    if wantz {
                        // Accumulate transformations in the matrix Z.
                        for j in iloz..=ihiz {
                            let sum = z[ix(j, k, ldz)] + v2 * z[ix(j, k + 1, ldz)] + v3 * z[ix(j, k + 2, ldz)];
                            z[ix(j, k, ldz)] -= sum * t1;
                            z[ix(j, k + 1, ldz)] -= sum * t2;
                            z[ix(j, k + 2, ldz)] -= sum * t3;
                        }
                    }
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
                    if wantz {
                        for j in iloz..=ihiz {
                            let sum = z[ix(j, k, ldz)] + v2 * z[ix(j, k + 1, ldz)];
                            z[ix(j, k, ldz)] -= sum * t1;
                            z[ix(j, k + 1, ldz)] -= sum * t2;
                        }
                    }
                }
            }
        }
        if !converged {
            // Failure to converge in remaining number of iterations.
            return iu;
        }
        // 150
        if l == iu {
            // H(I,I-1) is negligible: one eigenvalue has converged.
            wr[iu - 1] = h[ix(iu, iu, ldh)];
            wi[iu - 1] = 0.0;
        } else if l + 1 == iu {
            // H(I-1,I-2) is negligible: a pair of eigenvalues have
            // converged. Transform the 2x2 submatrix to standard Schur
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
                // Apply the transformation to the rest of H.
                if i2 > iu {
                    drot_within(i2 - iu, h, ix(iu - 1, iu + 1, ldh), ldh, ix(iu, iu + 1, ldh), ldh, cs, sn);
                }
                drot_within(iu - i1 - 1, h, ix(i1, iu - 1, ldh), 1, ix(i1, iu, ldh), 1, cs, sn);
            }
            if wantz {
                // Apply the transformation to Z.
                drot_within(nz, z, ix(iloz, iu - 1, ldz), 1, ix(iloz, iu, ldz), 1, cs, sn);
            }
        }
        // Reset deflation counter.
        kdefl = 0;
        // Return to start of the main loop with new value of I.
        i = l as isize - 1;
    }
}

/// `DLALN2(LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, LDB, WR, WI, X,
/// LDX, SCALE, XNORM, INFO)`: solves `(ca A - w D) X = s B` for a 1x1 or
/// 2x2 `A` (`na`), real or complex `w` (`nw`), with scaling `s <= 1` to
/// avoid overflow. `a` is `A(1,1)` onward with leading dimension `lda`, `b`
/// the right-hand side as `[b11, b21, b12, b22]` (column-major 2x2, the
/// columns being real and imaginary parts). Returns `(X, scale, xnorm,
/// info)` with `X` column-major 2x2.
#[allow(clippy::too_many_arguments)]
pub fn dlaln2(
    ltrans: bool,
    na: usize,
    nw: usize,
    smin: f64,
    ca: f64,
    a: &[f64],
    lda: usize,
    d1: f64,
    d2: f64,
    b: [f64; 4],
    wr: f64,
    wi: f64,
) -> ([f64; 4], f64, f64, i32) {
    const ZSWAP: [bool; 4] = [false, false, true, true];
    const RSWAP: [bool; 4] = [false, true, false, true];
    // IPIVOT(4,4), column ICMAX = [1..4].
    const IPIVOT: [[usize; 4]; 4] = [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]];
    let (b11, b21, b12, b22) = (b[0], b[1], b[2], b[3]);
    let mut x = [0.0f64; 4]; // X(1,1), X(2,1), X(1,2), X(2,2)
    let smlnum = 2.0 * SAFMIN;
    let bignum = 1.0 / smlnum;
    let smini = smin.max(smlnum);
    let mut info = 0;
    let mut scale = 1.0;
    let xnorm;
    let a11 = a[0];
    if na == 1 {
        if nw == 1 {
            // Real 1x1 system.
            let mut csr = ca * a11 - wr * d1;
            let mut cnorm = csr.abs();
            if cnorm < smini {
                csr = smini;
                cnorm = smini;
                info = 1;
            }
            let bnorm = b11.abs();
            if cnorm < 1.0 && bnorm > 1.0 && bnorm > bignum * cnorm {
                scale = 1.0 / bnorm;
            }
            x[0] = (b11 * scale) / csr;
            xnorm = x[0].abs();
        } else {
            // Complex 1x1 system (w is complex).
            let mut csr = ca * a11 - wr * d1;
            let mut csi = -wi * d1;
            let mut cnorm = csr.abs() + csi.abs();
            if cnorm < smini {
                csr = smini;
                csi = 0.0;
                cnorm = smini;
                info = 1;
            }
            let bnorm = b11.abs() + b12.abs();
            if cnorm < 1.0 && bnorm > 1.0 && bnorm > bignum * cnorm {
                scale = 1.0 / bnorm;
            }
            let (p, q) = dladiv(scale * b11, scale * b12, csr, csi);
            x[0] = p;
            x[2] = q;
            xnorm = x[0].abs() + x[2].abs();
        }
        return (x, scale, xnorm, info);
    }
    // 2x2 system: compute the real part of C = ca A - w D.
    let (a21, a12, a22) = (a[1], a[lda], a[lda + 1]);
    let mut cr = [0.0f64; 4]; // CR(1,1), CR(2,1), CR(1,2), CR(2,2)
    cr[0] = ca * a11 - wr * d1;
    cr[3] = ca * a22 - wr * d2;
    if ltrans {
        cr[2] = ca * a21;
        cr[1] = ca * a12;
    } else {
        cr[1] = ca * a21;
        cr[2] = ca * a12;
    }
    if nw == 1 {
        // Real 2x2 system (w is real): find the largest element in C.
        let mut cmax = 0.0;
        let mut icmax = 0usize;
        for j in 1..=4 {
            if cr[j - 1].abs() > cmax {
                cmax = cr[j - 1].abs();
                icmax = j;
            }
        }
        // If norm(C) < SMINI, use SMINI*identity.
        if cmax < smini {
            let bnorm = b11.abs().max(b21.abs());
            if smini < 1.0 && bnorm > 1.0 && bnorm > bignum * smini {
                scale = 1.0 / bnorm;
            }
            let temp = scale / smini;
            x[0] = temp * b11;
            x[1] = temp * b21;
            return (x, scale, temp * bnorm, 1);
        }
        // Gaussian elimination with complete pivoting.
        let piv = IPIVOT[icmax - 1];
        let ur11 = cr[icmax - 1];
        let cr21 = cr[piv[1] - 1];
        let ur12 = cr[piv[2] - 1];
        let cr22 = cr[piv[3] - 1];
        let ur11r = 1.0 / ur11;
        let lr21 = ur11r * cr21;
        let mut ur22 = cr22 - ur12 * lr21;
        // If smaller pivot < SMINI, use SMINI.
        if ur22.abs() < smini {
            ur22 = smini;
            info = 1;
        }
        let (br1, mut br2) = if RSWAP[icmax - 1] { (b21, b11) } else { (b11, b21) };
        br2 -= lr21 * br1;
        let bbnd = (br1 * (ur22 * ur11r)).abs().max(br2.abs());
        if bbnd > 1.0 && ur22.abs() < 1.0 && bbnd >= bignum * ur22.abs() {
            scale = 1.0 / bbnd;
        }
        let xr2 = (br2 * scale) / ur22;
        let xr1 = (scale * br1) * ur11r - xr2 * (ur11r * ur12);
        if ZSWAP[icmax - 1] {
            x[0] = xr2;
            x[1] = xr1;
        } else {
            x[0] = xr1;
            x[1] = xr2;
        }
        let mut xn = xr1.abs().max(xr2.abs());
        // Further scaling if norm(A) norm(X) > overflow.
        if xn > 1.0 && cmax > 1.0 && xn > bignum / cmax {
            let temp = cmax / bignum;
            x[0] *= temp;
            x[1] *= temp;
            xn *= temp;
            scale *= temp;
        }
        return (x, scale, xn, info);
    }
    // Complex 2x2 system (w is complex): find the largest element in C.
    let mut ci = [0.0f64; 4];
    ci[0] = -wi * d1;
    ci[1] = 0.0;
    ci[2] = 0.0;
    ci[3] = -wi * d2;
    let mut cmax = 0.0;
    let mut icmax = 0usize;
    for j in 1..=4 {
        if cr[j - 1].abs() + ci[j - 1].abs() > cmax {
            cmax = cr[j - 1].abs() + ci[j - 1].abs();
            icmax = j;
        }
    }
    // If norm(C) < SMINI, use SMINI*identity.
    if cmax < smini {
        let bnorm = (b11.abs() + b12.abs()).max(b21.abs() + b22.abs());
        if smini < 1.0 && bnorm > 1.0 && bnorm > bignum * smini {
            scale = 1.0 / bnorm;
        }
        let temp = scale / smini;
        x = [temp * b11, temp * b21, temp * b12, temp * b22];
        return (x, scale, temp * bnorm, 1);
    }
    // Gaussian elimination with complete pivoting.
    let piv = IPIVOT[icmax - 1];
    let ur11 = cr[icmax - 1];
    let ui11 = ci[icmax - 1];
    let cr21 = cr[piv[1] - 1];
    let ci21 = ci[piv[1] - 1];
    let ur12 = cr[piv[2] - 1];
    let ui12 = ci[piv[2] - 1];
    let cr22 = cr[piv[3] - 1];
    let ci22 = ci[piv[3] - 1];
    let (ur11r, ui11r, lr21, li21, ur12s, ui12s, mut ur22, mut ui22);
    if icmax == 1 || icmax == 4 {
        // Code when off-diagonals of pivoted C are real.
        if ur11.abs() > ui11.abs() {
            let temp = ui11 / ur11;
            ur11r = 1.0 / (ur11 * (1.0 + temp.powi(2)));
            ui11r = -temp * ur11r;
        } else {
            let temp = ur11 / ui11;
            ui11r = -1.0 / (ui11 * (1.0 + temp.powi(2)));
            ur11r = -temp * ui11r;
        }
        lr21 = cr21 * ur11r;
        li21 = cr21 * ui11r;
        ur12s = ur12 * ur11r;
        ui12s = ur12 * ui11r;
        ur22 = cr22 - ur12 * lr21;
        ui22 = ci22 - ur12 * li21;
    } else {
        // Code when diagonals of pivoted C are real.
        ur11r = 1.0 / ur11;
        ui11r = 0.0;
        lr21 = cr21 * ur11r;
        li21 = ci21 * ur11r;
        ur12s = ur12 * ur11r;
        ui12s = ui12 * ur11r;
        ur22 = cr22 - ur12 * lr21 + ui12 * li21;
        ui22 = -ur12 * li21 - ui12 * lr21;
    }
    let u22abs = ur22.abs() + ui22.abs();
    // If smaller pivot < SMINI, use SMINI.
    if u22abs < smini {
        ur22 = smini;
        ui22 = 0.0;
        info = 1;
    }
    let (mut br1, mut br2, mut bi1, mut bi2) =
        if RSWAP[icmax - 1] { (b21, b11, b22, b12) } else { (b11, b21, b12, b22) };
    br2 = br2 - lr21 * br1 + li21 * bi1;
    bi2 = bi2 - li21 * br1 - lr21 * bi1;
    let bbnd = ((br1.abs() + bi1.abs()) * (u22abs * (ur11r.abs() + ui11r.abs()))).max(br2.abs() + bi2.abs());
    if bbnd > 1.0 && u22abs < 1.0 && bbnd >= bignum * u22abs {
        scale = 1.0 / bbnd;
        br1 *= scale;
        bi1 *= scale;
        br2 *= scale;
        bi2 *= scale;
    }
    let (xr2, xi2) = dladiv(br2, bi2, ur22, ui22);
    let xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
    let xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
    if ZSWAP[icmax - 1] {
        x = [xr2, xr1, xi2, xi1];
    } else {
        x = [xr1, xr2, xi1, xi2];
    }
    let mut xn = (xr1.abs() + xi1.abs()).max(xr2.abs() + xi2.abs());
    // Further scaling if norm(A) norm(X) > overflow.
    if xn > 1.0 && cmax > 1.0 && xn > bignum / cmax {
        let temp = cmax / bignum;
        for v in &mut x {
            *v *= temp;
        }
        xn *= temp;
        scale *= temp;
    }
    (x, scale, xn, info)
}

/// `DTREVC('R', HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M,
/// WORK, INFO)` for `HOWMNY` = 'A' (`select` = `None`) or 'S': the right
/// eigenvectors of the quasi-triangular `t` into the columns of `vr`
/// (a complex pair as two columns, real then imaginary part), each scaled
/// so its largest component has `|re| + |im| = 1`. With 'S', `select` is
/// updated as the Fortran does (the second of a selected pair cleared).
/// Returns `M`, the number of columns used. `work` holds at least `3n`.
#[allow(clippy::too_many_arguments)]
pub fn dtrevc_right(
    mut select: Option<&mut [bool]>,
    n: usize,
    t: &[f64],
    ldt: usize,
    vr: &mut [f64],
    ldvr: usize,
    work: &mut [f64],
) -> usize {
    // Set M to the number of columns required to store the selected
    // eigenvectors, standardize the array SELECT if necessary.
    let m = if let Some(sel) = select.as_deref_mut() {
        let mut m = 0;
        let mut pair = false;
        for j in 1..=n {
            if pair {
                pair = false;
                sel[j - 1] = false;
            } else if j < n {
                if t[ix(j + 1, j, ldt)] == 0.0 {
                    if sel[j - 1] {
                        m += 1;
                    }
                } else {
                    pair = true;
                    if sel[j - 1] || sel[j] {
                        sel[j - 1] = true;
                        m += 2;
                    }
                }
            } else if sel[n - 1] {
                m += 1;
            }
        }
        m
    } else {
        n
    };
    if n == 0 {
        return m;
    }
    // Set the constants to control overflow.
    let unfl = SAFMIN;
    let ulp = PREC;
    let smlnum = unfl * (n as f64 / ulp);
    let bignum = (1.0 - ulp) / smlnum;
    // Compute 1-norm of each column of strictly upper triangular part of T
    // to control overflow in triangular solver.
    work[0] = 0.0;
    for j in 2..=n {
        work[j - 1] = 0.0;
        for i in 1..j {
            work[j - 1] += t[ix(i, j, ldt)].abs();
        }
    }
    let n2 = 2 * n;
    // WORK(k) as 1-based.
    macro_rules! w {
        ($k:expr) => {
            work[$k - 1]
        };
    }
    let mut ip: i32 = 0;
    let mut is = m;
    for ki in (1..=n).rev() {
        'body: {
            if ip == 1 {
                break 'body;
            }
            if ki != 1 && t[ix(ki, ki - 1, ldt)] != 0.0 {
                ip = -1;
            }
            if let Some(sel) = select.as_deref() {
                if ip == 0 {
                    if !sel[ki - 1] {
                        break 'body;
                    }
                } else if !sel[ki - 2] {
                    break 'body;
                }
            }
            // Compute the KI-th eigenvalue (WR,WI).
            let wr = t[ix(ki, ki, ldt)];
            let mut wi = 0.0;
            if ip != 0 {
                wi = t[ix(ki, ki - 1, ldt)].abs().sqrt() * t[ix(ki - 1, ki, ldt)].abs().sqrt();
            }
            let smin = (ulp * (wr.abs() + wi.abs())).max(smlnum);
            if ip == 0 {
                // Real right eigenvector.
                w!(ki + n) = 1.0;
                // Form right-hand side.
                for k in 1..ki {
                    w!(k + n) = -t[ix(k, ki, ldt)];
                }
                // Solve the upper quasi-triangular system:
                // (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
                let mut jnxt = ki as isize - 1;
                for j in (1..ki).rev() {
                    if j as isize > jnxt {
                        continue;
                    }
                    let mut j1 = j;
                    let j2 = j;
                    jnxt = j as isize - 1;
                    if j > 1 && t[ix(j, j - 1, ldt)] != 0.0 {
                        j1 = j - 1;
                        jnxt = j as isize - 2;
                    }
                    if j1 == j2 {
                        // 1-by-1 diagonal block.
                        let (mut x, mut scale, xnorm, _) =
                            dlaln2(false, 1, 1, smin, 1.0, &t[ix(j, j, ldt)..], ldt, 1.0, 1.0, [w!(j + n), 0.0, 0.0, 0.0], wr, 0.0);
                        // Scale X(1,1) to avoid overflow when updating the
                        // right-hand side.
                        if xnorm > 1.0 && w!(j) > bignum / xnorm {
                            x[0] /= xnorm;
                            scale /= xnorm;
                        }
                        // Scale if necessary.
                        if scale != 1.0 {
                            dscal(ki, scale, &mut work[n..]);
                        }
                        w!(j + n) = x[0];
                        // Update right-hand side.
                        daxpy_col(j - 1, -x[0], t, ix(1, j, ldt), work, n);
                    } else {
                        // 2-by-2 diagonal block.
                        let (mut x, mut scale, xnorm, _) = dlaln2(
                            false,
                            2,
                            1,
                            smin,
                            1.0,
                            &t[ix(j - 1, j - 1, ldt)..],
                            ldt,
                            1.0,
                            1.0,
                            [w!(j - 1 + n), w!(j + n), 0.0, 0.0],
                            wr,
                            0.0,
                        );
                        // Scale X(1,1) and X(2,1) to avoid overflow when
                        // updating the right-hand side.
                        if xnorm > 1.0 {
                            let beta = w!(j - 1).max(w!(j));
                            if beta > bignum / xnorm {
                                x[0] /= xnorm;
                                x[1] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        if scale != 1.0 {
                            dscal(ki, scale, &mut work[n..]);
                        }
                        w!(j - 1 + n) = x[0];
                        w!(j + n) = x[1];
                        daxpy_col(j - 2, -x[0], t, ix(1, j - 1, ldt), work, n);
                        daxpy_col(j - 2, -x[1], t, ix(1, j, ldt), work, n);
                    }
                }
                // Copy the vector x or Q*x to VR and normalize.
                let c = ix(1, is, ldvr);
                dcopy(ki, &work[n..], &mut vr[c..]);
                let ii = idamax(ki, &vr[c..]);
                let remax = 1.0 / vr[c + ii - 1].abs();
                dscal(ki, remax, &mut vr[c..]);
                for k in ki + 1..=n {
                    vr[ix(k, is, ldvr)] = 0.0;
                }
            } else {
                // Complex right eigenvector. Initial solve
                // [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
                //   [ (T(KI,KI-1)   T(KI,KI)   )               ]
                if t[ix(ki - 1, ki, ldt)].abs() >= t[ix(ki, ki - 1, ldt)].abs() {
                    w!(ki - 1 + n) = 1.0;
                    w!(ki + n2) = wi / t[ix(ki - 1, ki, ldt)];
                } else {
                    w!(ki - 1 + n) = -wi / t[ix(ki, ki - 1, ldt)];
                    w!(ki + n2) = 1.0;
                }
                w!(ki + n) = 0.0;
                w!(ki - 1 + n2) = 0.0;
                // Form right-hand side.
                for k in 1..ki - 1 {
                    w!(k + n) = -w!(ki - 1 + n) * t[ix(k, ki - 1, ldt)];
                    w!(k + n2) = -w!(ki + n2) * t[ix(k, ki, ldt)];
                }
                // Solve upper quasi-triangular system:
                // (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
                let mut jnxt = ki as isize - 2;
                for j in (1..ki - 1).rev() {
                    if j as isize > jnxt {
                        continue;
                    }
                    let mut j1 = j;
                    let j2 = j;
                    jnxt = j as isize - 1;
                    if j > 1 && t[ix(j, j - 1, ldt)] != 0.0 {
                        j1 = j - 1;
                        jnxt = j as isize - 2;
                    }
                    if j1 == j2 {
                        // 1-by-1 diagonal block.
                        let (mut x, mut scale, xnorm, _) = dlaln2(
                            false,
                            1,
                            2,
                            smin,
                            1.0,
                            &t[ix(j, j, ldt)..],
                            ldt,
                            1.0,
                            1.0,
                            [w!(j + n), 0.0, w!(j + n2), 0.0],
                            wr,
                            wi,
                        );
                        // Scale X(1,1) and X(1,2) to avoid overflow when
                        // updating the right-hand side.
                        if xnorm > 1.0 && w!(j) > bignum / xnorm {
                            x[0] /= xnorm;
                            x[2] /= xnorm;
                            scale /= xnorm;
                        }
                        // Scale if necessary.
                        if scale != 1.0 {
                            dscal(ki, scale, &mut work[n..]);
                            dscal(ki, scale, &mut work[n2..]);
                        }
                        w!(j + n) = x[0];
                        w!(j + n2) = x[2];
                        // Update the right-hand side.
                        daxpy_col(j - 1, -x[0], t, ix(1, j, ldt), work, n);
                        daxpy_col(j - 1, -x[2], t, ix(1, j, ldt), work, n2);
                    } else {
                        // 2-by-2 diagonal block.
                        let (mut x, mut scale, xnorm, _) = dlaln2(
                            false,
                            2,
                            2,
                            smin,
                            1.0,
                            &t[ix(j - 1, j - 1, ldt)..],
                            ldt,
                            1.0,
                            1.0,
                            [w!(j - 1 + n), w!(j + n), w!(j - 1 + n2), w!(j + n2)],
                            wr,
                            wi,
                        );
                        // Scale X to avoid overflow when updating the
                        // right-hand side.
                        if xnorm > 1.0 {
                            let beta = w!(j - 1).max(w!(j));
                            if beta > bignum / xnorm {
                                let rec = 1.0 / xnorm;
                                for v in &mut x {
                                    *v *= rec;
                                }
                                scale *= rec;
                            }
                        }
                        // Scale if necessary.
                        if scale != 1.0 {
                            dscal(ki, scale, &mut work[n..]);
                            dscal(ki, scale, &mut work[n2..]);
                        }
                        w!(j - 1 + n) = x[0];
                        w!(j + n) = x[1];
                        w!(j - 1 + n2) = x[2];
                        w!(j + n2) = x[3];
                        // Update the right-hand side.
                        daxpy_col(j - 2, -x[0], t, ix(1, j - 1, ldt), work, n);
                        daxpy_col(j - 2, -x[1], t, ix(1, j, ldt), work, n);
                        daxpy_col(j - 2, -x[2], t, ix(1, j - 1, ldt), work, n2);
                        daxpy_col(j - 2, -x[3], t, ix(1, j, ldt), work, n2);
                    }
                }
                // Copy the vector x or Q*x to VR and normalize.
                let (c1, c2) = (ix(1, is - 1, ldvr), ix(1, is, ldvr));
                dcopy(ki, &work[n..], &mut vr[c1..]);
                dcopy(ki, &work[n2..], &mut vr[c2..]);
                let mut emax = 0.0f64;
                for k in 0..ki {
                    emax = emax.max(vr[c1 + k].abs() + vr[c2 + k].abs());
                }
                let remax = 1.0 / emax;
                dscal(ki, remax, &mut vr[c1..]);
                dscal(ki, remax, &mut vr[c2..]);
                for k in ki + 1..=n {
                    vr[ix(k, is - 1, ldvr)] = 0.0;
                    vr[ix(k, is, ldvr)] = 0.0;
                }
            }
            is -= 1;
            if ip != 0 {
                is -= 1;
            }
        }
        // 130
        if ip == 1 {
            ip = 0;
        }
        if ip == -1 {
            ip = 1;
        }
    }
    m
}

/// `DAXPY(N, A, T(1,J), 1, WORK(1+OFF), 1)` with `T(1,J)` at `t[ts]`.
fn daxpy_col(n: usize, a: f64, t: &[f64], ts: usize, work: &mut [f64], off: usize) {
    daxpy(n, a, &t[ts..], &mut work[off..]);
}

/// `DLARFG` on a 3-vector held in a fixed array, for callers that keep
/// the Householder vector in `v`.
pub fn dlarfg_v(nr: usize, v: &mut [f64; 3]) -> f64 {
    let (beta, tau) = dlarfg(nr, v[0], &mut v[1..nr]);
    v[0] = beta;
    tau
}

/// `DLARFX(SIDE, M, N, V, TAU, C, LDC, WORK)` for a reflector of order 3
/// (the unrolled code of the Fortran): `C := H C` (`left`, `m == 3`) or
/// `C := C H` (`n == 3`), `H = I - tau v v'`.
#[allow(clippy::too_many_arguments)]
pub fn dlarfx3(left: bool, m: usize, n: usize, v: &[f64; 3], tau: f64, c: &mut [f64], cs: usize, ldc: usize) {
    if tau == 0.0 {
        return;
    }
    let (v1, v2, v3) = (v[0], v[1], v[2]);
    let (t1, t2, t3) = (tau * v1, tau * v2, tau * v3);
    if left {
        debug_assert_eq!(m, 3);
        for j in 0..n {
            let b = cs + j * ldc;
            let sum = v1 * c[b] + v2 * c[b + 1] + v3 * c[b + 2];
            c[b] -= sum * t1;
            c[b + 1] -= sum * t2;
            c[b + 2] -= sum * t3;
        }
    } else {
        debug_assert_eq!(n, 3);
        for j in 0..m {
            let b = cs + j;
            let sum = v1 * c[b] + v2 * c[b + ldc] + v3 * c[b + 2 * ldc];
            c[b] -= sum * t1;
            c[b + ldc] -= sum * t2;
            c[b + 2 * ldc] -= sum * t3;
        }
    }
}

/// `DLANGE('M', M, N, A, LDA, WORK)`: the largest absolute entry.
pub fn dlange_max(m: usize, n: usize, a: &[f64], lda: usize) -> f64 {
    let mut value = 0.0f64;
    for j in 1..=n {
        for i in 1..=m {
            let temp = a[ix(i, j, lda)].abs();
            if value < temp || temp.is_nan() {
                value = temp;
            }
        }
    }
    value
}

/// `DLASY2(.FALSE., .FALSE., ISGN, N1, N2, TL, LDTL, TR, LDTR, B, LDB,
/// SCALE, X, LDX, XNORM, INFO)`: solves `TL X + isgn X TR = scale B` for
/// `n1`, `n2` in 1..=2. `tl`, `tr` and `b` start at their `(1,1)` entries.
/// Returns `(scale, X, xnorm, info)`, `X` column-major 2x2.
#[allow(clippy::too_many_arguments)]
pub fn dlasy2_nn(
    isgn: i32,
    n1: usize,
    n2: usize,
    tl: &[f64],
    ldtl: usize,
    tr: &[f64],
    ldtr: usize,
    b: &[f64],
    ldb: usize,
) -> (f64, [f64; 4], f64, i32) {
    const LOCU12: [usize; 4] = [3, 4, 1, 2];
    const LOCL21: [usize; 4] = [2, 1, 4, 3];
    const LOCU22: [usize; 4] = [4, 3, 2, 1];
    const XSWPIV: [bool; 4] = [false, false, true, true];
    const BSWPIV: [bool; 4] = [false, true, false, true];
    let tl_ = |i: usize, j: usize| tl[ix(i, j, ldtl)];
    let tr_ = |i: usize, j: usize| tr[ix(i, j, ldtr)];
    let b_ = |i: usize, j: usize| b[ix(i, j, ldb)];
    let mut x = [0.0f64; 4]; // X(1,1), X(2,1), X(1,2), X(2,2)
    let mut info = 0;
    if n1 == 0 || n2 == 0 {
        return (1.0, x, 0.0, 0);
    }
    let eps = PREC;
    let smlnum = SAFMIN / eps;
    let sgn = f64::from(isgn);
    let k = n1 + n1 + n2 - 2;
    if k == 1 {
        // 1 by 1: TL11*X + SGN*X*TR11 = B11
        let mut tau1 = tl_(1, 1) + sgn * tr_(1, 1);
        let mut bet = tau1.abs();
        if bet <= smlnum {
            tau1 = smlnum;
            bet = smlnum;
            info = 1;
        }
        let mut scale = 1.0;
        let gam = b_(1, 1).abs();
        if smlnum * gam > bet {
            scale = 1.0 / gam;
        }
        x[0] = (b_(1, 1) * scale) / tau1;
        return (scale, x, x[0].abs(), info);
    }
    if k == 2 || k == 3 {
        // 1 by 2 (k == 2) or 2 by 1 (k == 3).
        let mut tmp = [0.0f64; 4];
        let mut btmp = [0.0f64; 4];
        let smin;
        if k == 2 {
            smin = (eps * tl_(1, 1).abs().max(tr_(1, 1).abs()).max(tr_(1, 2).abs()).max(tr_(2, 1).abs()).max(tr_(2, 2).abs()))
                .max(smlnum);
            tmp[0] = tl_(1, 1) + sgn * tr_(1, 1);
            tmp[3] = tl_(1, 1) + sgn * tr_(2, 2);
            tmp[1] = sgn * tr_(1, 2);
            tmp[2] = sgn * tr_(2, 1);
            btmp[0] = b_(1, 1);
            btmp[1] = b_(1, 2);
        } else {
            smin = (eps * tr_(1, 1).abs().max(tl_(1, 1).abs()).max(tl_(1, 2).abs()).max(tl_(2, 1).abs()).max(tl_(2, 2).abs()))
                .max(smlnum);
            tmp[0] = tl_(1, 1) + sgn * tr_(1, 1);
            tmp[3] = tl_(2, 2) + sgn * tr_(1, 1);
            tmp[1] = tl_(2, 1);
            tmp[2] = tl_(1, 2);
            btmp[0] = b_(1, 1);
            btmp[1] = b_(2, 1);
        }
        // 40: solve 2 by 2 system using complete pivoting. Set pivots
        // near zero to SMIN.
        let ipiv = idamax(4, &tmp);
        let mut u11 = tmp[ipiv - 1];
        if u11.abs() <= smin {
            info = 1;
            u11 = smin;
        }
        let u12 = tmp[LOCU12[ipiv - 1] - 1];
        let l21 = tmp[LOCL21[ipiv - 1] - 1] / u11;
        let mut u22 = tmp[LOCU22[ipiv - 1] - 1] - u12 * l21;
        let xswap = XSWPIV[ipiv - 1];
        let bswap = BSWPIV[ipiv - 1];
        if u22.abs() <= smin {
            info = 1;
            u22 = smin;
        }
        if bswap {
            let temp = btmp[1];
            btmp[1] = btmp[0] - l21 * temp;
            btmp[0] = temp;
        } else {
            btmp[1] -= l21 * btmp[0];
        }
        let mut scale = 1.0;
        if (2.0 * smlnum) * btmp[1].abs() > u22.abs() || (2.0 * smlnum) * btmp[0].abs() > u11.abs() {
            scale = 0.5 / btmp[0].abs().max(btmp[1].abs());
            btmp[0] *= scale;
            btmp[1] *= scale;
        }
        let mut x2 = [0.0f64; 2];
        x2[1] = btmp[1] / u22;
        x2[0] = btmp[0] / u11 - (u12 / u11) * x2[1];
        if xswap {
            x2.swap(0, 1);
        }
        x[0] = x2[0];
        let xnorm = if n1 == 1 {
            x[2] = x2[1];
            x[0].abs() + x[2].abs()
        } else {
            x[1] = x2[1];
            x[0].abs().max(x[1].abs())
        };
        return (scale, x, xnorm, info);
    }
    // 50: 2 by 2: solve equivalent 4 by 4 system using complete pivoting.
    // Set pivots near zero to SMIN.
    let smin = tr_(1, 1).abs().max(tr_(1, 2).abs()).max(tr_(2, 1).abs()).max(tr_(2, 2).abs());
    let smin = smin.max(tl_(1, 1).abs()).max(tl_(1, 2).abs()).max(tl_(2, 1).abs()).max(tl_(2, 2).abs());
    let smin = (eps * smin).max(smlnum);
    // T16(i,j) column-major 4x4.
    let mut t16 = [0.0f64; 16];
    let t = |i: usize, j: usize| (i - 1) + (j - 1) * 4;
    t16[t(1, 1)] = tl_(1, 1) + sgn * tr_(1, 1);
    t16[t(2, 2)] = tl_(2, 2) + sgn * tr_(1, 1);
    t16[t(3, 3)] = tl_(1, 1) + sgn * tr_(2, 2);
    t16[t(4, 4)] = tl_(2, 2) + sgn * tr_(2, 2);
    t16[t(1, 2)] = tl_(1, 2);
    t16[t(2, 1)] = tl_(2, 1);
    t16[t(3, 4)] = tl_(1, 2);
    t16[t(4, 3)] = tl_(2, 1);
    t16[t(1, 3)] = sgn * tr_(2, 1);
    t16[t(2, 4)] = sgn * tr_(2, 1);
    t16[t(3, 1)] = sgn * tr_(1, 2);
    t16[t(4, 2)] = sgn * tr_(1, 2);
    let mut btmp = [b_(1, 1), b_(2, 1), b_(1, 2), b_(2, 2)];
    let mut jpiv = [0usize; 4];
    for i in 1..=3 {
        let mut xmax = 0.0;
        let (mut ipsv, mut jpsv) = (0usize, 0usize);
        for ip in i..=4 {
            for jp in i..=4 {
                if t16[t(ip, jp)].abs() >= xmax {
                    xmax = t16[t(ip, jp)].abs();
                    ipsv = ip;
                    jpsv = jp;
                }
            }
        }
        if ipsv != i {
            for c in 1..=4 {
                t16.swap(t(ipsv, c), t(i, c));
            }
            btmp.swap(i - 1, ipsv - 1);
        }
        if jpsv != i {
            for r in 1..=4 {
                t16.swap(t(r, jpsv), t(r, i));
            }
        }
        jpiv[i - 1] = jpsv;
        if t16[t(i, i)].abs() < smin {
            info = 1;
            t16[t(i, i)] = smin;
        }
        for j in i + 1..=4 {
            t16[t(j, i)] /= t16[t(i, i)];
            btmp[j - 1] -= t16[t(j, i)] * btmp[i - 1];
            for k in i + 1..=4 {
                t16[t(j, k)] -= t16[t(j, i)] * t16[t(i, k)];
            }
        }
    }
    if t16[t(4, 4)].abs() < smin {
        t16[t(4, 4)] = smin;
    }
    let mut scale = 1.0;
    if (8.0 * smlnum) * btmp[0].abs() > t16[t(1, 1)].abs()
        || (8.0 * smlnum) * btmp[1].abs() > t16[t(2, 2)].abs()
        || (8.0 * smlnum) * btmp[2].abs() > t16[t(3, 3)].abs()
        || (8.0 * smlnum) * btmp[3].abs() > t16[t(4, 4)].abs()
    {
        scale = (1.0 / 8.0) / btmp[0].abs().max(btmp[1].abs()).max(btmp[2].abs()).max(btmp[3].abs());
        for v in &mut btmp {
            *v *= scale;
        }
    }
    let mut tmp = [0.0f64; 4];
    for i in 1..=4 {
        let k = 5 - i;
        let temp = 1.0 / t16[t(k, k)];
        tmp[k - 1] = btmp[k - 1] * temp;
        for j in k + 1..=4 {
            tmp[k - 1] -= (temp * t16[t(k, j)]) * tmp[j - 1];
        }
    }
    for i in 1..=3 {
        let p = jpiv[4 - i - 1];
        if p != 4 - i {
            tmp.swap(4 - i - 1, p - 1);
        }
    }
    x = [tmp[0], tmp[1], tmp[2], tmp[3]];
    let xnorm = (tmp[0].abs() + tmp[2].abs()).max(tmp[1].abs() + tmp[3].abs());
    (scale, x, xnorm, info)
}

/// `DLANV2` applied in place to `T(j,j), T(j,j+1), T(j+1,j), T(j+1,j+1)`;
/// returns `(cs, sn)`.
fn dlanv2_at(t: &mut [f64], ldt: usize, j: usize) -> (f64, f64) {
    let (a, b, c, d, _, _, _, _, cs, sn) = dlanv2(t[ix(j, j, ldt)], t[ix(j, j + 1, ldt)], t[ix(j + 1, j, ldt)], t[ix(j + 1, j + 1, ldt)]);
    t[ix(j, j, ldt)] = a;
    t[ix(j, j + 1, ldt)] = b;
    t[ix(j + 1, j, ldt)] = c;
    t[ix(j + 1, j + 1, ldt)] = d;
    (cs, sn)
}

/// `DLAEXC(WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, INFO)`: swaps the
/// adjacent diagonal blocks of orders `n1` and `n2` starting at row `j1` of
/// the Schur form `t`, updating `q` when `wantq`. Returns `INFO` (1: the
/// swap was rejected, the result would be too far from Schur form).
#[allow(clippy::too_many_arguments)]
pub fn dlaexc(wantq: bool, n: usize, t: &mut [f64], ldt: usize, q: &mut [f64], ldq: usize, j1: usize, n1: usize, n2: usize) -> i32 {
    if n == 0 || n1 == 0 || n2 == 0 {
        return 0;
    }
    if j1 + n1 > n {
        return 0;
    }
    let j2 = j1 + 1;
    let j3 = j1 + 2;
    let j4 = j1 + 3;
    if n1 == 1 && n2 == 1 {
        // Swap two 1-by-1 blocks.
        let t11 = t[ix(j1, j1, ldt)];
        let t22 = t[ix(j2, j2, ldt)];
        // Determine the transformation to perform the interchange.
        let (cs, sn, _) = dlartg(t[ix(j1, j2, ldt)], t22 - t11);
        // Apply transformation to the matrix T.
        if j3 <= n {
            drot_within(n - j1 - 1, t, ix(j1, j3, ldt), ldt, ix(j2, j3, ldt), ldt, cs, sn);
        }
        drot_within(j1 - 1, t, ix(1, j1, ldt), 1, ix(1, j2, ldt), 1, cs, sn);
        t[ix(j1, j1, ldt)] = t22;
        t[ix(j2, j2, ldt)] = t11;
        if wantq {
            drot_within(n, q, ix(1, j1, ldq), 1, ix(1, j2, ldq), 1, cs, sn);
        }
        return 0;
    }
    // Swapping involves at least one 2-by-2 block. Copy the diagonal block
    // of order N1+N2 to the local array D and compute its norm.
    const LDD: usize = 4;
    let nd = n1 + n2;
    let mut d = [0.0f64; 16];
    for j in 1..=nd {
        for i in 1..=nd {
            d[ix(i, j, LDD)] = t[ix(j1 + i - 1, j1 + j - 1, ldt)];
        }
    }
    let dnorm = dlange_max(nd, nd, &d, LDD);
    // Compute machine-dependent threshold for test for accepting swap.
    let eps = PREC;
    let smlnum = SAFMIN / eps;
    let thresh = (10.0 * eps * dnorm).max(smlnum);
    // Solve T11*X - X*T22 = scale*T12 for X.
    let (scale, x, _, _) = dlasy2_nn(-1, n1, n2, &d, LDD, &d[ix(n1 + 1, n1 + 1, LDD)..], LDD, &d[ix(1, n1 + 1, LDD)..], LDD);
    let x_ = |i: usize, j: usize| x[(i - 1) + (j - 1) * 2];
    let k = n1 + n1 + n2 - 3;
    match k {
        1 => {
            // N1 = 1, N2 = 2: generate elementary reflector H so that
            // ( scale, X11, X12 ) H = ( 0, 0, * ).
            let mut u = [scale, x_(1, 1), x_(1, 2)];
            let (beta, tau) = dlarfg(3, u[2], &mut u[0..2]);
            u[2] = beta;
            u[2] = 1.0;
            let t11 = t[ix(j1, j1, ldt)];
            // Perform swap provisionally on diagonal block in D.
            dlarfx3(true, 3, 3, &u, tau, &mut d, 0, LDD);
            dlarfx3(false, 3, 3, &u, tau, &mut d, 0, LDD);
            // Test whether to reject swap.
            if d[ix(3, 1, LDD)].abs().max(d[ix(3, 2, LDD)].abs()).max((d[ix(3, 3, LDD)] - t11).abs()) > thresh {
                return 1;
            }
            // Accept swap: apply transformation to the entire matrix T.
            dlarfx3(true, 3, n - j1 + 1, &u, tau, t, ix(j1, j1, ldt), ldt);
            dlarfx3(false, j2, 3, &u, tau, t, ix(1, j1, ldt), ldt);
            t[ix(j3, j1, ldt)] = 0.0;
            t[ix(j3, j2, ldt)] = 0.0;
            t[ix(j3, j3, ldt)] = t11;
            if wantq {
                // Accumulate transformation in the matrix Q.
                dlarfx3(false, n, 3, &u, tau, q, ix(1, j1, ldq), ldq);
            }
        }
        2 => {
            // N1 = 2, N2 = 1: generate elementary reflector H so that
            // H ( -X11 ) = ( * )
            //   ( -X21 )   ( 0 )
            //   ( scale)   ( 0 )
            let mut u = [-x_(1, 1), -x_(2, 1), scale];
            let (beta, tau) = dlarfg(3, u[0], &mut u[1..3]);
            u[0] = beta;
            u[0] = 1.0;
            let t33 = t[ix(j3, j3, ldt)];
            // Perform swap provisionally on diagonal block in D.
            dlarfx3(true, 3, 3, &u, tau, &mut d, 0, LDD);
            dlarfx3(false, 3, 3, &u, tau, &mut d, 0, LDD);
            // Test whether to reject swap.
            if d[ix(2, 1, LDD)].abs().max(d[ix(3, 1, LDD)].abs()).max((d[ix(1, 1, LDD)] - t33).abs()) > thresh {
                return 1;
            }
            // Accept swap: apply transformation to the entire matrix T.
            dlarfx3(false, j3, 3, &u, tau, t, ix(1, j1, ldt), ldt);
            dlarfx3(true, 3, n - j1, &u, tau, t, ix(j1, j2, ldt), ldt);
            t[ix(j1, j1, ldt)] = t33;
            t[ix(j2, j1, ldt)] = 0.0;
            t[ix(j3, j1, ldt)] = 0.0;
            if wantq {
                dlarfx3(false, n, 3, &u, tau, q, ix(1, j1, ldq), ldq);
            }
        }
        _ => {
            // N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2)
            // so that:
            // H(2) H(1) (  -X11  -X12 ) = (  *  * )
            //           (  -X21  -X22 )   (  0  * )
            //           ( scale    0  )   (  0  0 )
            //           (    0  scale )   (  0  0 )
            let mut u1 = [-x_(1, 1), -x_(2, 1), scale];
            let (beta, tau1) = dlarfg(3, u1[0], &mut u1[1..3]);
            u1[0] = beta;
            u1[0] = 1.0;
            let temp = -tau1 * (x_(1, 2) + u1[1] * x_(2, 2));
            let mut u2 = [-temp * u1[1] - x_(2, 2), -temp * u1[2], scale];
            let (beta, tau2) = dlarfg(3, u2[0], &mut u2[1..3]);
            u2[0] = beta;
            u2[0] = 1.0;
            // Perform swap provisionally on diagonal block in D.
            dlarfx3(true, 3, 4, &u1, tau1, &mut d, 0, LDD);
            dlarfx3(false, 4, 3, &u1, tau1, &mut d, 0, LDD);
            dlarfx3(true, 3, 4, &u2, tau2, &mut d, ix(2, 1, LDD), LDD);
            dlarfx3(false, 4, 3, &u2, tau2, &mut d, ix(1, 2, LDD), LDD);
            // Test whether to reject swap.
            if d[ix(3, 1, LDD)].abs().max(d[ix(3, 2, LDD)].abs()).max(d[ix(4, 1, LDD)].abs()).max(d[ix(4, 2, LDD)].abs()) > thresh {
                return 1;
            }
            // Accept swap: apply transformation to the entire matrix T.
            dlarfx3(true, 3, n - j1 + 1, &u1, tau1, t, ix(j1, j1, ldt), ldt);
            dlarfx3(false, j4, 3, &u1, tau1, t, ix(1, j1, ldt), ldt);
            dlarfx3(true, 3, n - j1 + 1, &u2, tau2, t, ix(j2, j1, ldt), ldt);
            dlarfx3(false, j4, 3, &u2, tau2, t, ix(1, j2, ldt), ldt);
            t[ix(j3, j1, ldt)] = 0.0;
            t[ix(j3, j2, ldt)] = 0.0;
            t[ix(j4, j1, ldt)] = 0.0;
            t[ix(j4, j2, ldt)] = 0.0;
            if wantq {
                dlarfx3(false, n, 3, &u1, tau1, q, ix(1, j1, ldq), ldq);
                dlarfx3(false, n, 3, &u2, tau2, q, ix(1, j2, ldq), ldq);
            }
        }
    }
    // 40
    if n2 == 2 {
        // Standardize new 2-by-2 block T11.
        let (cs, sn) = dlanv2_at(t, ldt, j1);
        drot_within(n - j1 - 1, t, ix(j1, j1 + 2, ldt), ldt, ix(j2, j1 + 2, ldt), ldt, cs, sn);
        drot_within(j1 - 1, t, ix(1, j1, ldt), 1, ix(1, j2, ldt), 1, cs, sn);
        if wantq {
            drot_within(n, q, ix(1, j1, ldq), 1, ix(1, j2, ldq), 1, cs, sn);
        }
    }
    if n1 == 2 {
        // Standardize new 2-by-2 block T22.
        let j3 = j1 + n2;
        let j4 = j3 + 1;
        let (cs, sn) = dlanv2_at(t, ldt, j3);
        if j3 + 2 <= n {
            drot_within(n - j3 - 1, t, ix(j3, j3 + 2, ldt), ldt, ix(j4, j3 + 2, ldt), ldt, cs, sn);
        }
        drot_within(j3 - 1, t, ix(1, j3, ldt), 1, ix(1, j4, ldt), 1, cs, sn);
        if wantq {
            drot_within(n, q, ix(1, j3, ldq), 1, ix(1, j4, ldq), 1, cs, sn);
        }
    }
    0
}

/// `DTREXC(COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, INFO)` with `COMPQ`
/// 'V' (`wantq`) or 'N': moves the diagonal block of `t` at row `ifst` to
/// row `ilst` by adjacent swaps. `ifst` and `ilst` are updated as the
/// Fortran does; returns `INFO` (1: a swap was rejected).
#[allow(clippy::too_many_arguments)]
pub fn dtrexc(wantq: bool, n: usize, t: &mut [f64], ldt: usize, q: &mut [f64], ldq: usize, ifst: &mut usize, ilst: &mut usize) -> i32 {
    if n <= 1 {
        return 0;
    }
    let tt = |t: &[f64], i: usize, j: usize| t[ix(i, j, ldt)];
    // Determine the first row of specified block and find out it is 1 by 1
    // or 2 by 2.
    if *ifst > 1 && tt(t, *ifst, *ifst - 1) != 0.0 {
        *ifst -= 1;
    }
    let mut nbf = 1;
    if *ifst < n && tt(t, *ifst + 1, *ifst) != 0.0 {
        nbf = 2;
    }
    // Determine the first row of the final block and find out it is 1 by
    // 1 or 2 by 2.
    if *ilst > 1 && tt(t, *ilst, *ilst - 1) != 0.0 {
        *ilst -= 1;
    }
    let mut nbl = 1;
    if *ilst < n && tt(t, *ilst + 1, *ilst) != 0.0 {
        nbl = 2;
    }
    if *ifst == *ilst {
        return 0;
    }
    let mut here;
    if *ifst < *ilst {
        // Update ILST.
        if nbf == 2 && nbl == 1 {
            *ilst -= 1;
        }
        if nbf == 1 && nbl == 2 {
            *ilst += 1;
        }
        here = *ifst;
        loop {
            // 10: swap block with next one below.
            if nbf == 1 || nbf == 2 {
                // Current block either 1 by 1 or 2 by 2.
                let mut nbnext = 1;
                if here + nbf < n && tt(t, here + nbf + 1, here + nbf) != 0.0 {
                    nbnext = 2;
                }
                let info = dlaexc(wantq, n, t, ldt, q, ldq, here, nbf, nbnext);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                here += nbnext;
                // Test if 2 by 2 block breaks into two 1 by 1 blocks.
                if nbf == 2 && tt(t, here + 1, here) == 0.0 {
                    nbf = 3;
                }
            } else {
                // Current block consists of two 1 by 1 blocks each of which
                // must be swapped individually.
                let mut nbnext = 1;
                if here + 3 <= n && tt(t, here + 3, here + 2) != 0.0 {
                    nbnext = 2;
                }
                let info = dlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, nbnext);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                if nbnext == 1 {
                    // Swap two 1 by 1 blocks, no problems possible.
                    dlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext);
                    here += 1;
                } else {
                    // Recompute NBNEXT in case 2 by 2 split.
                    if tt(t, here + 2, here + 1) == 0.0 {
                        nbnext = 1;
                    }
                    if nbnext == 2 {
                        // 2 by 2 block did not split.
                        let info = dlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext);
                        if info != 0 {
                            *ilst = here;
                            return info;
                        }
                        here += 2;
                    } else {
                        // 2 by 2 block did split.
                        dlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1);
                        dlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, 1);
                        here += 2;
                    }
                }
            }
            if here >= *ilst {
                break;
            }
        }
    } else {
        here = *ifst;
        loop {
            // 20: swap block with next one above.
            if nbf == 1 || nbf == 2 {
                // Current block either 1 by 1 or 2 by 2.
                let mut nbnext = 1;
                if here >= 3 && tt(t, here - 1, here - 2) != 0.0 {
                    nbnext = 2;
                }
                let info = dlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                here -= nbnext;
                // Test if 2 by 2 block breaks into two 1 by 1 blocks.
                if nbf == 2 && tt(t, here + 1, here) == 0.0 {
                    nbf = 3;
                }
            } else {
                // Current block consists of two 1 by 1 blocks each of which
                // must be swapped individually.
                let mut nbnext = 1;
                if here >= 3 && tt(t, here - 1, here - 2) != 0.0 {
                    nbnext = 2;
                }
                let info = dlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, 1);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                if nbnext == 1 {
                    // Swap two 1 by 1 blocks, no problems possible.
                    dlaexc(wantq, n, t, ldt, q, ldq, here, nbnext, 1);
                    here -= 1;
                } else {
                    // Recompute NBNEXT in case 2 by 2 split.
                    if tt(t, here, here - 1) == 0.0 {
                        nbnext = 1;
                    }
                    if nbnext == 2 {
                        // 2 by 2 block did not split.
                        let info = dlaexc(wantq, n, t, ldt, q, ldq, here - 1, 2, 1);
                        if info != 0 {
                            *ilst = here;
                            return info;
                        }
                        here -= 2;
                    } else {
                        // 2 by 2 block did split.
                        dlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1);
                        dlaexc(wantq, n, t, ldt, q, ldq, here - 1, 1, 1);
                        here -= 2;
                    }
                }
            }
            if here <= *ilst {
                break;
            }
        }
    }
    *ilst = here;
    0
}

/// `DTRSEN('N', 'V', SELECT, N, T, LDT, Q, LDQ, WR, WI, M, S, SEP, WORK,
/// LWORK, IWORK, LIWORK, INFO)`: reorders the Schur form `t` (and the
/// Schur vectors `q`) so the selected eigenvalues lead, a complex pair
/// moving together. Returns `(M, INFO)` with the eigenvalues of the
/// reordered form in `wr`, `wi` (`INFO = 1`: a swap was rejected; the
/// eigenvalues are of the partly reordered form).
#[allow(clippy::too_many_arguments)]
pub fn dtrsen_nv(
    select: &[bool],
    n: usize,
    t: &mut [f64],
    ldt: usize,
    q: &mut [f64],
    ldq: usize,
    wr: &mut [f64],
    wi: &mut [f64],
) -> (usize, i32) {
    // Set M to the dimension of the specified invariant subspace.
    let mut m = 0;
    let mut pair = false;
    for k in 1..=n {
        if pair {
            pair = false;
        } else if k < n {
            if t[ix(k + 1, k, ldt)] == 0.0 {
                if select[k - 1] {
                    m += 1;
                }
            } else {
                pair = true;
                if select[k - 1] || select[k] {
                    m += 2;
                }
            }
        } else if select[n - 1] {
            m += 1;
        }
    }
    let mut info = 0;
    if m != n && m != 0 {
        // Collect the selected blocks at the top-left corner of T.
        let mut ks = 0;
        let mut pair = false;
        for k in 1..=n {
            if pair {
                pair = false;
                continue;
            }
            let mut swap = select[k - 1];
            if k < n && t[ix(k + 1, k, ldt)] != 0.0 {
                pair = true;
                swap = swap || select[k];
            }
            if swap {
                ks += 1;
                // Swap the K-th block to position KS.
                let mut ierr = 0;
                let mut kk = k;
                if k != ks {
                    // KS is DTREXC's ILST, updated in place.
                    ierr = dtrexc(true, n, t, ldt, q, ldq, &mut kk, &mut ks);
                }
                if ierr == 1 || ierr == 2 {
                    // Blocks too close to swap: exit.
                    info = 1;
                    break;
                }
                if pair {
                    ks += 1;
                }
            }
        }
    }
    // 40: store the output eigenvalues in WR and WI.
    for k in 1..=n {
        wr[k - 1] = t[ix(k, k, ldt)];
        wi[k - 1] = 0.0;
    }
    for k in 1..n {
        if t[ix(k + 1, k, ldt)] != 0.0 {
            wi[k - 1] = t[ix(k, k + 1, ldt)].abs().sqrt() * t[ix(k + 1, k, ldt)].abs().sqrt();
            wi[k] = -wi[k - 1];
        }
    }
    (m, info)
}
