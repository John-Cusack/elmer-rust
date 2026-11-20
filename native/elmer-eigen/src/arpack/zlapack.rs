//! The reference LAPACK 3.12 routines under ARPACK's complex driver
//! (`znaupd` / `zneupd`), as the pristine Elmer build links them from the
//! system `liblapack.so.3`: `ZLAHQR`, `ZLARFG`, `ZLARTG`, `ZLARNV`,
//! `ZLACPY`, `ZLASET`, `ZLANHS('1')`, `ZLASCL('G')`, `ZLADIV` and
//! `DLAPY3`, with glibc's `csqrt` for the Fortran `SQRT` of a complex.
//!
//! Matrices are column-major slices with a leading dimension, addressed
//! with Fortran's 1-based indices (`A(i,j)` is `a[(i-1) + (j-1)*lda]`).
//! Complex quotients go through [`cdiv`] (gfortran's Smith division).

// Index loops and argument lists mirror the Fortran's.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use num_complex::Complex64 as C;

use super::lapack::{dlaruv, EPS, PREC, SAFMIN};
use super::lapack_ns::dladiv;
use super::zblas::{cabs, cabs1, cdiv, dznrm2, zcopy, zdscal, zscal, ONE, ZERO};

/// `ZSCAL` on a strided run of `a`: `count` entries from `start`, `step`
/// apart (nothing when `za` is one, as the 3.12 reference does).
fn zscal_strided(a: &mut [C], start: usize, step: usize, count: usize, za: C) {
    if za == ONE {
        return;
    }
    for k in 0..count {
        let p = start + k * step;
        a[p] = za * a[p];
    }
}

#[inline]
fn ix(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

/// glibc's `csqrt` (s_csqrt_template.c) for finite arguments, which
/// gfortran calls for `SQRT` of a `COMPLEX*16`.
pub fn csqrt(x: C) -> C {
    let (mut re, mut im) = (x.re, x.im);
    if im == 0.0 {
        return if re < 0.0 {
            C::new(0.0, (-re).sqrt().copysign(im))
        } else {
            C::new(re.sqrt().abs(), 0f64.copysign(im))
        };
    }
    if re == 0.0 {
        let r = if im.abs() >= 2.0 * f64::MIN_POSITIVE { (0.5 * im.abs()).sqrt() } else { 0.5 * (2.0 * im.abs()).sqrt() };
        return C::new(r, r.copysign(im));
    }
    let scalbn = |v: f64, e: i32| v * 2f64.powi(e);
    let mut scale = 0i32;
    if re.abs() > f64::MAX / 4.0 {
        scale = 1;
        re = scalbn(re, -2);
        im = scalbn(im, -2);
    }
    if im.abs() > f64::MAX / 4.0 {
        scale = 1;
        re = if re.abs() >= 4.0 * f64::MIN_POSITIVE { scalbn(re, -2) } else { 0.0 };
        im = scalbn(im, -2);
    } else if re.abs() < 2.0 * f64::MIN_POSITIVE && im.abs() < 2.0 * f64::MIN_POSITIVE {
        scale = -((f64::MANTISSA_DIGITS as i32 + 1) / 2);
        re = scalbn(re, -2 * scale);
        im = scalbn(im, -2 * scale);
    }
    let d = re.hypot(im);
    let (mut r, mut s);
    if re > 0.0 {
        r = (0.5 * (d + re)).sqrt();
        if scale == 1 && im.abs() < 1.0 {
            s = im / r;
            r = scalbn(r, scale);
            scale = 0;
        } else {
            s = 0.5 * (im / r);
        }
    } else {
        s = (0.5 * (d - re)).sqrt();
        if scale == 1 && im.abs() < 1.0 {
            r = (im / s).abs();
            s = scalbn(s, scale);
            scale = 0;
        } else {
            r = (0.5 * (im / s)).abs();
        }
    }
    if scale != 0 {
        r = scalbn(r, scale);
        s = scalbn(s, scale);
    }
    C::new(r, s.copysign(im))
}

/// `ZLADIV(X, Y)`: `X / Y` with `DLADIV`'s robust real arithmetic.
pub fn zladiv(x: C, y: C) -> C {
    let (zr, zi) = dladiv(x.re, x.im, y.re, y.im);
    C::new(zr, zi)
}

/// `DLAPY3(X, Y, Z)` (3.12): `sqrt(x^2 + y^2 + z^2)` avoiding overflow.
pub fn dlapy3(x: f64, y: f64, z: f64) -> f64 {
    let hugeval = f64::MAX;
    let (xa, ya, za) = (x.abs(), y.abs(), z.abs());
    let w = xa.max(ya).max(za);
    if w == 0.0 || w > hugeval {
        xa + ya + za
    } else {
        w * ((xa / w) * (xa / w) + (ya / w) * (ya / w) + (za / w) * (za / w)).sqrt()
    }
}

/// `ZLARFG(N, ALPHA, X, 1, TAU)`: the elementary reflector with
/// `H^H (alpha; x) = (beta; 0)`, `beta` real. Returns `tau`.
pub fn zlarfg(n: usize, alpha: &mut C, x: &mut [C]) -> C {
    if n == 0 {
        return ZERO;
    }
    let mut xnorm = dznrm2(n - 1, x);
    let mut alphr = alpha.re;
    let mut alphi = alpha.im;
    if xnorm == 0.0 && alphi == 0.0 {
        return ZERO;
    }
    let mut beta = -dlapy3(alphr, alphi, xnorm).copysign(alphr);
    let safmin = SAFMIN / EPS;
    let rsafmn = 1.0 / safmin;
    let mut knt = 0;
    if beta.abs() < safmin {
        loop {
            knt += 1;
            zdscal(n - 1, rsafmn, x);
            beta *= rsafmn;
            alphi *= rsafmn;
            alphr *= rsafmn;
            if !(beta.abs() < safmin && knt < 20) {
                break;
            }
        }
        xnorm = dznrm2(n - 1, x);
        *alpha = C::new(alphr, alphi);
        beta = -dlapy3(alphr, alphi, xnorm).copysign(alphr);
    }
    let tau = C::new((beta - alphr) / beta, -alphi / beta);
    let a = zladiv(ONE, *alpha - beta);
    zscal(n - 1, a, x);
    for _ in 0..knt {
        beta *= safmin;
    }
    *alpha = C::new(beta, 0.0);
    tau
}

/// `ZLARTG(F, G, C, S, R)` (zlartg.f90, LAPACK 3.10+): the plane rotation
/// with `[c s; -conj(s) c] (f; g) = (r; 0)`. Returns `(c, s, r)`.
pub fn zlartg(f: C, g: C) -> (f64, C, C) {
    let safmin = SAFMIN;
    let safmax = 1.0 / safmin;
    let rtmin = safmin.sqrt();
    let abssq = |t: C| t.re * t.re + t.im * t.im;
    if g == ZERO {
        return (1.0, ZERO, f);
    }
    if f == ZERO {
        let (s, r);
        if g.re == 0.0 {
            r = g.im.abs();
            s = g.conj() / r;
        } else if g.im == 0.0 {
            r = g.re.abs();
            s = g.conj() / r;
        } else {
            let g1 = g.re.abs().max(g.im.abs());
            let rtmax = (safmax / 2.0).sqrt();
            if g1 > rtmin && g1 < rtmax {
                let g2 = abssq(g);
                let d = g2.sqrt();
                s = g.conj() / d;
                r = d;
            } else {
                let u = safmax.min(safmin.max(g1));
                let gs = g / u;
                let g2 = abssq(gs);
                let d = g2.sqrt();
                s = gs.conj() / d;
                r = d * u;
            }
        }
        return (0.0, s, C::new(r, 0.0));
    }
    let f1 = f.re.abs().max(f.im.abs());
    let g1 = g.re.abs().max(g.im.abs());
    let mut rtmax = (safmax / 4.0).sqrt();
    let (mut c, s, mut r);
    if f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax {
        let f2 = abssq(f);
        let g2 = abssq(g);
        let h2 = f2 + g2;
        if f2 >= h2 * safmin {
            c = (f2 / h2).sqrt();
            r = f / c;
            rtmax *= 2.0;
            if f2 > rtmin && h2 < rtmax {
                s = g.conj() * (f / (f2 * h2).sqrt());
            } else {
                s = g.conj() * (r / h2);
            }
        } else {
            let d = (f2 * h2).sqrt();
            c = f2 / d;
            if c >= safmin {
                r = f / c;
            } else {
                r = f * (h2 / d);
            }
            s = g.conj() * (f / d);
        }
    } else {
        let u = safmax.min(safmin.max(f1).max(g1));
        let gs = g / u;
        let g2 = abssq(gs);
        let (w, fs, f2, h2);
        if f1 / u < rtmin {
            let v = safmax.min(safmin.max(f1));
            w = v / u;
            fs = f / v;
            f2 = abssq(fs);
            h2 = f2 * (w * w) + g2;
        } else {
            w = 1.0;
            fs = f / u;
            f2 = abssq(fs);
            h2 = f2 + g2;
        }
        let sv;
        if f2 >= h2 * safmin {
            c = (f2 / h2).sqrt();
            r = fs / c;
            rtmax *= 2.0;
            if f2 > rtmin && h2 < rtmax {
                sv = gs.conj() * (fs / (f2 * h2).sqrt());
            } else {
                sv = gs.conj() * (r / h2);
            }
        } else {
            let d = (f2 * h2).sqrt();
            c = f2 / d;
            if c >= safmin {
                r = fs / c;
            } else {
                r = fs * (h2 / d);
            }
            sv = gs.conj() * (fs / d);
        }
        c *= w;
        r *= u;
        s = sv;
    }
    (c, s, r)
}

/// `ZLARNV(2, ISEED, N, X)`: complex numbers with real and imaginary parts
/// uniform on (-1, 1).
pub fn zlarnv2(iseed: &mut [i32; 4], n: usize, x: &mut [C]) {
    const LV: usize = 128;
    let mut u = [0.0f64; LV];
    let mut iv = 0;
    while iv < n {
        let il = (LV / 2).min(n - iv);
        dlaruv(iseed, 2 * il, &mut u);
        for i in 0..il {
            x[iv + i] = C::new(2.0 * u[2 * i] - 1.0, 2.0 * u[2 * i + 1] - 1.0);
        }
        iv += LV / 2;
    }
}

/// `ZLACPY('A', M, N, A, LDA, B, LDB)`.
pub fn zlacpy(m: usize, n: usize, a: &[C], lda: usize, b: &mut [C], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = a[i + j * lda];
        }
    }
}

/// `ZLASET('A', M, N, ALPHA, BETA, A, LDA)`.
pub fn zlaset(m: usize, n: usize, alpha: C, beta: C, a: &mut [C], lda: usize) {
    for j in 0..n {
        for i in 0..m {
            a[i + j * lda] = alpha;
        }
    }
    for i in 0..m.min(n) {
        a[i + i * lda] = beta;
    }
}

/// `ZLANHS('1', N, A, LDA, WORK)`: the one norm of an upper Hessenberg
/// matrix.
pub fn zlanhs_one(n: usize, a: &[C], lda: usize) -> f64 {
    let mut value = 0.0f64;
    for j in 1..=n {
        let mut sum = 0.0;
        for i in 1..=n.min(j + 1) {
            sum += cabs(a[ix(i, j, lda)]);
        }
        if value < sum || sum.is_nan() {
            value = sum;
        }
    }
    value
}

/// `ZLASCL('G', 0, 0, CFROM, CTO, M, N, A, LDA, INFO)`: `A = A * (cto /
/// cfrom)` without over- or underflow.
pub fn zlascl_g(cfrom: f64, cto: f64, m: usize, n: usize, a: &mut [C], lda: usize) {
    let smlnum = SAFMIN;
    let bignum = 1.0 / smlnum;
    let (mut cfromc, mut ctoc) = (cfrom, cto);
    loop {
        let cfrom1 = cfromc * smlnum;
        let (mul, done);
        if cfrom1 == cfromc {
            mul = ctoc / cfromc;
            done = true;
        } else {
            let cto1 = ctoc / bignum;
            if cto1 == ctoc {
                mul = ctoc;
                done = true;
                cfromc = 1.0;
            } else if cfrom1.abs() > ctoc.abs() && ctoc != 0.0 {
                mul = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if cto1.abs() > cfromc.abs() {
                mul = bignum;
                done = false;
                ctoc = cto1;
            } else {
                mul = ctoc / cfromc;
                done = true;
                if mul == 1.0 {
                    return;
                }
            }
        }
        for j in 0..n {
            for i in 0..m {
                let v = a[i + j * lda];
                a[i + j * lda] = v * mul;
            }
        }
        if done {
            return;
        }
    }
}

/// `ZLAHQR(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ,
/// INFO)` (LAPACK 3.12): the single-shift complex QR algorithm on rows and
/// columns `ILO..IHI` of an upper Hessenberg `H`. Returns `INFO`.
pub fn zlahqr(
    wantt: bool,
    wantz: bool,
    n: usize,
    ilo: usize,
    ihi: usize,
    h: &mut [C],
    ldh: usize,
    w: &mut [C],
    iloz: usize,
    ihiz: usize,
    z: &mut [C],
    ldz: usize,
) -> usize {
    const DAT1: f64 = 3.0 / 4.0;
    const KEXSH: usize = 10;
    if n == 0 {
        return 0;
    }
    if ilo == ihi {
        w[ilo - 1] = h[ix(ilo, ilo, ldh)];
        return 0;
    }
    // Clear out the trash.
    for j in ilo..=ihi.saturating_sub(3) {
        h[ix(j + 2, j, ldh)] = ZERO;
        h[ix(j + 3, j, ldh)] = ZERO;
    }
    if ilo + 2 <= ihi {
        h[ix(ihi, ihi - 2, ldh)] = ZERO;
    }
    let (jlo, jhi) = if wantt { (1, n) } else { (ilo, ihi) };
    // Ensure that the subdiagonal entries are real.
    for i in ilo + 1..=ihi {
        let hii1 = h[ix(i, i - 1, ldh)];
        if hii1.im != 0.0 {
            let sc = hii1 / cabs1(hii1);
            let sc = sc.conj() / cabs(sc);
            h[ix(i, i - 1, ldh)] = C::new(cabs(hii1), 0.0);
            // ZSCAL(JHI-I+1, SC, H(I,I), LDH): along row I.
            zscal_strided(h, ix(i, i, ldh), ldh, jhi + 1 - i, sc);
            // ZSCAL(MIN(JHI,I+1)-JLO+1, CONJG(SC), H(JLO,I), 1).
            let cs = sc.conj();
            zscal_strided(h, ix(jlo, i, ldh), 1, jhi.min(i + 1) + 1 - jlo, cs);
            if wantz {
                zscal_strided(z, ix(iloz, i, ldz), 1, ihiz + 1 - iloz, cs);
            }
        }
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
    let mut i = ihi as isize;
    loop {
        if i < ilo as isize {
            return 0;
        }
        let iu = i as usize;
        let mut l = ilo;
        let mut converged = false;
        for _its in 0..=itmax {
            // Look for a single small subdiagonal element.
            let mut k = iu;
            while k > l {
                if cabs1(h[ix(k, k - 1, ldh)]) <= smlnum {
                    break;
                }
                let mut tst = cabs1(h[ix(k - 1, k - 1, ldh)]) + cabs1(h[ix(k, k, ldh)]);
                if tst == 0.0 {
                    if k >= ilo + 2 {
                        tst += h[ix(k - 1, k - 2, ldh)].re.abs();
                    }
                    if k < ihi {
                        tst += h[ix(k + 1, k, ldh)].re.abs();
                    }
                }
                // Ahues & Kressner's conservative deflation criterion.
                if h[ix(k, k - 1, ldh)].re.abs() <= ulp * tst {
                    let (hkk1, hk1k) = (cabs1(h[ix(k, k - 1, ldh)]), cabs1(h[ix(k - 1, k, ldh)]));
                    let ab = hkk1.max(hk1k);
                    let ba = hkk1.min(hk1k);
                    let hkk = h[ix(k, k, ldh)];
                    let dif = cabs1(h[ix(k - 1, k - 1, ldh)] - hkk);
                    let aa = cabs1(hkk).max(dif);
                    let bb = cabs1(hkk).min(dif);
                    let s = aa + ab;
                    if ba * (ab / s) <= smlnum.max(ulp * (bb * (aa / s))) {
                        break;
                    }
                }
                k -= 1;
            }
            l = k;
            if l > ilo {
                h[ix(l, l - 1, ldh)] = ZERO;
            }
            // A submatrix of order 1 has split off.
            if l >= iu {
                converged = true;
                break;
            }
            kdefl += 1;
            if !wantt {
                i1 = l;
                i2 = iu;
            }
            let mut t;
            if kdefl.is_multiple_of(2 * KEXSH) {
                // Exceptional shift.
                let s = DAT1 * h[ix(iu, iu - 1, ldh)].re.abs();
                t = s + h[ix(iu, iu, ldh)];
            } else if kdefl.is_multiple_of(KEXSH) {
                let s = DAT1 * h[ix(l + 1, l, ldh)].re.abs();
                t = s + h[ix(l, l, ldh)];
            } else {
                // Wilkinson's shift.
                t = h[ix(iu, iu, ldh)];
                let u = csqrt(h[ix(iu - 1, iu, ldh)]) * csqrt(h[ix(iu, iu - 1, ldh)]);
                let mut s = cabs1(u);
                if s != 0.0 {
                    let x = 0.5 * (h[ix(iu - 1, iu - 1, ldh)] - t);
                    let sx = cabs1(x);
                    s = s.max(cabs1(x));
                    let xs = x / s;
                    let us = u / s;
                    let mut y = s * csqrt(xs * xs + us * us);
                    if sx > 0.0 {
                        let xx = x / sx;
                        if xx.re * y.re + xx.im * y.im < 0.0 {
                            y = -y;
                        }
                    }
                    t -= u * zladiv(u, x + y);
                }
            }
            // Look for two consecutive small subdiagonal elements.
            let mut v = [ZERO; 2];
            let mut m = iu - 1;
            let mut found = false;
            while m > l {
                let h11 = h[ix(m, m, ldh)];
                let h22 = h[ix(m + 1, m + 1, ldh)];
                let mut h11s = h11 - t;
                let mut h21 = h[ix(m + 1, m, ldh)].re;
                let s = cabs1(h11s) + h21.abs();
                h11s /= s;
                h21 /= s;
                v[0] = h11s;
                v[1] = C::new(h21, 0.0);
                let h10 = h[ix(m, m - 1, ldh)].re;
                if h10.abs() * h21.abs() <= ulp * (cabs1(h11s) * (cabs1(h11) + cabs1(h22))) {
                    found = true;
                    break;
                }
                m -= 1;
            }
            if !found {
                m = l;
                let h11 = h[ix(l, l, ldh)];
                let mut h11s = h11 - t;
                let mut h21 = h[ix(l + 1, l, ldh)].re;
                let s = cabs1(h11s) + h21.abs();
                h11s /= s;
                h21 /= s;
                v[0] = h11s;
                v[1] = C::new(h21, 0.0);
            }
            // Single-shift QR step.
            for k in m..iu {
                if k > m {
                    zcopy(2, &h[ix(k, k - 1, ldh)..ix(k, k - 1, ldh) + 2], &mut v);
                }
                let (mut a0, mut x1) = (v[0], [v[1]]);
                let t1 = zlarfg(2, &mut a0, &mut x1);
                v = [a0, x1[0]];
                if k > m {
                    h[ix(k, k - 1, ldh)] = v[0];
                    h[ix(k + 1, k - 1, ldh)] = ZERO;
                }
                let v2 = v[1];
                let t2 = (t1 * v2).re;
                // G from the left, columns K..I2.
                for j in k..=i2 {
                    let (p, q) = (ix(k, j, ldh), ix(k + 1, j, ldh));
                    let sum = t1.conj() * h[p] + t2 * h[q];
                    h[p] -= sum;
                    h[q] -= sum * v2;
                }
                // G from the right, rows I1..MIN(K+2,I).
                for j in i1..=(k + 2).min(iu) {
                    let (p, q) = (ix(j, k, ldh), ix(j, k + 1, ldh));
                    let sum = t1 * h[p] + t2 * h[q];
                    h[p] -= sum;
                    h[q] -= sum * v2.conj();
                }
                if wantz {
                    for j in iloz..=ihiz {
                        let (p, q) = (ix(j, k, ldz), ix(j, k + 1, ldz));
                        let sum = t1 * z[p] + t2 * z[q];
                        z[p] -= sum;
                        z[q] -= sum * v2.conj();
                    }
                }
                if k == m && m > l {
                    // Keep H(M,M-1) real after a start at row M > L.
                    let mut temp = ONE - t1;
                    temp /= cabs(temp);
                    let p = ix(m + 1, m, ldh);
                    h[p] *= temp.conj();
                    if m + 2 <= iu {
                        let p = ix(m + 2, m + 1, ldh);
                        h[p] *= temp;
                    }
                    for j in m..=iu {
                        if j != m + 1 {
                            if i2 > j {
                                zscal_strided(h, ix(j, j + 1, ldh), ldh, i2 - j, temp);
                            }
                            let ct = temp.conj();
                            zscal_strided(h, ix(i1, j, ldh), 1, j - i1, ct);
                            if wantz {
                                zscal_strided(z, ix(iloz, j, ldz), 1, nz, ct);
                            }
                        }
                    }
                }
            }
            // Ensure that H(I,I-1) is real.
            let temp = h[ix(iu, iu - 1, ldh)];
            if temp.im != 0.0 {
                let rtemp = cabs(temp);
                h[ix(iu, iu - 1, ldh)] = C::new(rtemp, 0.0);
                let temp = temp / rtemp;
                if i2 > iu {
                    zscal_strided(h, ix(iu, iu + 1, ldh), ldh, i2 - iu, temp.conj());
                }
                zscal_strided(h, ix(i1, iu, ldh), 1, iu - i1, temp);
                if wantz {
                    zscal_strided(z, ix(iloz, iu, ldz), 1, nz, temp);
                }
            }
        }
        if !converged {
            return iu;
        }
        // One eigenvalue has converged.
        w[iu - 1] = h[ix(iu, iu, ldh)];
        kdefl = 0;
        i = l as isize - 1;
    }
}

/// `ZTRSV('U', 'N', 'N', N, A, LDA, X, 1)`: `x = A^-1 x`, `A` upper
/// triangular.
pub fn ztrsv_unn(n: usize, a: &[C], lda: usize, x: &mut [C]) {
    for j in (1..=n).rev() {
        if x[j - 1] != ZERO {
            x[j - 1] = cdiv(x[j - 1], a[ix(j, j, lda)]);
            let temp = x[j - 1];
            for i in (1..j).rev() {
                x[i - 1] -= temp * a[ix(i, j, lda)];
            }
        }
    }
}

/// `ZLATRS('Upper', 'No transpose', 'Non-unit', 'Y', N, A, LDA, X, SCALE,
/// CNORM, INFO)`: `A x = s b` with the scale `s` chosen against overflow;
/// `cnorm` holds the column norms of the strictly upper part. Returns
/// `SCALE`.
pub fn zlatrs_unn(n: usize, a: &[C], lda: usize, x: &mut [C], cnorm: &mut [f64]) -> f64 {
    const HALF: f64 = 0.5;
    let cabs2 = |z: C| (z.re / 2.0).abs() + (z.im / 2.0).abs();
    if n == 0 {
        return 1.0;
    }
    let mut smlnum = SAFMIN;
    smlnum /= PREC;
    let bignum = 1.0 / smlnum;
    let mut scale = 1.0;
    // IDAMAX on CNORM.
    let mut imax = 0;
    for j in 1..n {
        if cnorm[j].abs() > cnorm[imax].abs() {
            imax = j;
        }
    }
    let tmax = cnorm[imax];
    let tscal;
    if tmax <= bignum * HALF {
        tscal = 1.0;
    } else if tmax <= f64::MAX {
        tscal = HALF / (smlnum * tmax);
        for c in &mut cnorm[..n] {
            *c *= tscal;
        }
    } else {
        // Entries of CNORM beyond the overflow threshold (DLANGE path).
        panic!("ZLATRS: column norms beyond the overflow threshold are not ported");
    }
    let mut xmax = 0.0f64;
    for j in 0..n {
        xmax = xmax.max(cabs2(x[j]));
    }
    let mut xbnd = xmax;
    let mut grow;
    if tscal != 1.0 {
        grow = 0.0;
    } else {
        grow = HALF / xbnd.max(smlnum);
        xbnd = grow;
        let mut j = n;
        let mut early = false;
        while j >= 1 {
            if grow <= smlnum {
                early = true;
                break;
            }
            let tjj = cabs1(a[ix(j, j, lda)]);
            if tjj >= smlnum {
                xbnd = xbnd.min(1f64.min(tjj) * grow);
            } else {
                xbnd = 0.0;
            }
            if tjj + cnorm[j - 1] >= smlnum {
                grow *= tjj / (tjj + cnorm[j - 1]);
            } else {
                grow = 0.0;
            }
            j -= 1;
        }
        if !early {
            grow = xbnd;
        }
    }
    if grow * tscal > smlnum {
        ztrsv_unn(n, a, lda, x);
    } else {
        if xmax > bignum * HALF {
            scale = (bignum * HALF) / xmax;
            zdscal(n, scale, x);
            xmax = bignum;
        } else {
            xmax *= 2.0;
        }
        for j in (1..=n).rev() {
            let mut xj = cabs1(x[j - 1]);
            let tjjs = a[ix(j, j, lda)] * tscal;
            let tjj = cabs1(tjjs);
            if tjj > smlnum {
                if tjj < 1.0 && xj > tjj * bignum {
                    let rec = 1.0 / xj;
                    zdscal(n, rec, x);
                    scale *= rec;
                    xmax *= rec;
                }
                x[j - 1] = zladiv(x[j - 1], tjjs);
                xj = cabs1(x[j - 1]);
            } else if tjj > 0.0 {
                if xj > tjj * bignum {
                    let mut rec = (tjj * bignum) / xj;
                    if cnorm[j - 1] > 1.0 {
                        rec /= cnorm[j - 1];
                    }
                    zdscal(n, rec, x);
                    scale *= rec;
                    xmax *= rec;
                }
                x[j - 1] = zladiv(x[j - 1], tjjs);
                xj = cabs1(x[j - 1]);
            } else {
                for v in &mut x[..n] {
                    *v = ZERO;
                }
                x[j - 1] = ONE;
                xj = 1.0;
                scale = 0.0;
                xmax = 0.0;
            }
            if xj > 1.0 {
                let mut rec = 1.0 / xj;
                if cnorm[j - 1] > (bignum - xmax) * rec {
                    rec *= HALF;
                    zdscal(n, rec, x);
                    scale *= rec;
                }
            } else if xj * cnorm[j - 1] > bignum - xmax {
                zdscal(n, HALF, x);
                scale *= HALF;
            }
            if j > 1 {
                let za = -x[j - 1] * tscal;
                let col: Vec<C> = (1..j).map(|i| a[ix(i, j, lda)]).collect();
                super::zblas::zaxpy(j - 1, za, &col, x);
                let i = super::zblas::izamax(j - 1, x);
                xmax = cabs1(x[i - 1]);
            }
        }
        scale /= tscal;
    }
    if tscal != 1.0 {
        for c in &mut cnorm[..n] {
            *c *= 1.0 / tscal;
        }
    }
    scale
}

/// `ZTREVC('Right', 'Back', SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M,
/// WORK, RWORK, INFO)`: the right eigenvectors of the upper triangular `T`,
/// back-transformed by the matrix `VR` holds on entry, each normalized to
/// unit largest `|re| + |im|`. `T`'s diagonal is restored on return.
pub fn ztrevc_rb(n: usize, t: &mut [C], ldt: usize, vr: &mut [C], ldvr: usize) {
    if n == 0 {
        return;
    }
    let unfl = SAFMIN;
    let ulp = PREC;
    let smlnum = unfl * (n as f64 / ulp);
    let mut work = vec![ZERO; 2 * n];
    let mut rwork = vec![0.0f64; n];
    for i in 1..=n {
        work[i + n - 1] = t[ix(i, i, ldt)];
    }
    rwork[0] = 0.0;
    for j in 2..=n {
        let col: Vec<C> = (1..j).map(|i| t[ix(i, j, ldt)]).collect();
        rwork[j - 1] = super::zblas::dzasum(j - 1, &col);
    }
    for ki in (1..=n).rev() {
        let smin = (ulp * cabs1(t[ix(ki, ki, ldt)])).max(smlnum);
        work[0] = ONE;
        for k in 1..ki {
            work[k - 1] = -t[ix(k, ki, ldt)];
        }
        let tkk = t[ix(ki, ki, ldt)];
        for k in 1..ki {
            let p = ix(k, k, ldt);
            t[p] -= tkk;
            if cabs1(t[p]) < smin {
                t[p] = C::new(smin, 0.0);
            }
        }
        if ki > 1 {
            let scale = zlatrs_unn(ki - 1, t, ldt, &mut work[..ki - 1], &mut rwork[..ki - 1]);
            work[ki - 1] = C::new(scale, 0.0);
            // ZGEMV('N', N, KI-1, ONE, VR, LDVR, WORK, 1, DCMPLX(SCALE), VR(1,KI), 1).
            let mut y: Vec<C> = (1..=n).map(|i| vr[ix(i, ki, ldvr)]).collect();
            super::zblas::zgemv_n(n, ki - 1, ONE, vr, ldvr, &work[..ki - 1], C::new(scale, 0.0), &mut y);
            for i in 1..=n {
                vr[ix(i, ki, ldvr)] = y[i - 1];
            }
        }
        let col: Vec<C> = (1..=n).map(|i| vr[ix(i, ki, ldvr)]).collect();
        let ii = super::zblas::izamax(n, &col);
        let remax = 1.0 / cabs1(col[ii - 1]);
        let mut col = col;
        zdscal(n, remax, &mut col);
        for i in 1..=n {
            vr[ix(i, ki, ldvr)] = col[i - 1];
        }
        for k in 1..ki {
            t[ix(k, k, ldt)] = work[k + n - 1];
        }
    }
}

/// `ILAZLC(M, N, A, LDA)`: the last column of `A` with a nonzero entry.
pub fn ilazlc(m: usize, n: usize, a: &[C], lda: usize) -> usize {
    if n == 0 {
        return 0;
    }
    if a[ix(1, n, lda)] != ZERO || a[ix(m.max(1), n, lda)] != ZERO {
        return n;
    }
    for j in (1..=n).rev() {
        for i in 1..=m {
            if a[ix(i, j, lda)] != ZERO {
                return j;
            }
        }
    }
    0
}

/// `ILAZLR(M, N, A, LDA)`: the last row of `A` with a nonzero entry.
pub fn ilazlr(m: usize, n: usize, a: &[C], lda: usize) -> usize {
    if m == 0 {
        return 0;
    }
    if a[ix(m, 1, lda)] != ZERO || a[ix(m, n.max(1), lda)] != ZERO {
        return m;
    }
    let mut r = 0;
    for j in 1..=n {
        let mut i = m;
        while i >= 1 && a[ix(i.max(1), j, lda)] == ZERO {
            i -= 1;
        }
        r = r.max(i);
    }
    r
}

/// `ZGERC(M, N, ALPHA, X, 1, Y, 1, A, LDA)`: `A = A + alpha x y^H`.
pub fn zgerc(m: usize, n: usize, alpha: C, x: &[C], y: &[C], a: &mut [C], lda: usize) {
    if m == 0 || n == 0 || alpha == ZERO {
        return;
    }
    for j in 0..n {
        if y[j] != ZERO {
            let temp = alpha * y[j].conj();
            for i in 0..m {
                a[i + j * lda] += x[i] * temp;
            }
        }
    }
}

/// `ZLARF(SIDE, M, N, V, 1, TAU, C, LDC, WORK)` (3.12, with the trailing
/// zeros of `v` and of `C` skipped): `C = H C` (`left`) or `C H`, `H = I -
/// tau v v^H`.
pub fn zlarf(left: bool, m: usize, n: usize, v: &[C], tau: C, c: &mut [C], ldc: usize) {
    let (mut lastv, mut lastc) = (0usize, 0usize);
    if tau != ZERO {
        lastv = if left { m } else { n };
        while lastv > 0 && v[lastv - 1] == ZERO {
            lastv -= 1;
        }
        lastc = if left { ilazlc(lastv, n, c, ldc) } else { ilazlr(m, lastv, c, ldc) };
    }
    if lastv == 0 {
        return;
    }
    if left {
        let mut work = vec![ZERO; lastc];
        super::zblas::zgemv_c(lastv, lastc, ONE, c, ldc, v, ZERO, &mut work);
        zgerc(lastv, lastc, -tau, v, &work, c, ldc);
    } else {
        let mut work = vec![ZERO; lastc];
        super::zblas::zgemv_n(lastc, lastv, ONE, c, ldc, v, ZERO, &mut work);
        zgerc(lastc, lastv, -tau, &work, v, c, ldc);
    }
}

/// `ZGEQR2(M, N, A, LDA, TAU, WORK, INFO)`: the unblocked QR
/// factorization. Returns `TAU`.
pub fn zgeqr2(m: usize, n: usize, a: &mut [C], lda: usize) -> Vec<C> {
    let k = m.min(n);
    let mut tau = vec![ZERO; k];
    for i in 1..=k {
        let p = ix(i, i, lda);
        let mut alpha = a[p];
        let q = ix(m.min(i + 1), i, lda);
        let len = m - i;
        let mut xs: Vec<C> = a[q..q + len].to_vec();
        tau[i - 1] = zlarfg(m - i + 1, &mut alpha, &mut xs);
        a[p] = alpha;
        a[q..q + len].copy_from_slice(&xs);
        if i < n {
            let aii = a[p];
            a[p] = ONE;
            let v: Vec<C> = a[p..p + (m - i + 1)].to_vec();
            let c0 = ix(i, i + 1, lda);
            zlarf(true, m - i + 1, n - i, &v, tau[i - 1].conj(), &mut a[c0..], lda);
            a[p] = aii;
        }
    }
    tau
}

/// `ZUNM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO)`:
/// `C` times the unitary `Q` of [`zgeqr2`]'s `K` reflectors (`left`:
/// `Q C`; `conj`: with `Q^H`).
pub fn zunm2r(left: bool, conj: bool, m: usize, n: usize, k: usize, a: &mut [C], lda: usize, tau: &[C], c: &mut [C], ldc: usize) {
    let notran = !conj;
    let forward = (left && !notran) || (!left && notran);
    let order: Vec<usize> = if forward { (1..=k).collect() } else { (1..=k).rev().collect() };
    for i in order {
        let (mi, ni, ic, jc) = if left { (m - i + 1, n, i, 1) } else { (m, n - i + 1, 1, i) };
        let taui = if notran { tau[i - 1] } else { tau[i - 1].conj() };
        let p = ix(i, i, lda);
        let aii = a[p];
        a[p] = ONE;
        let vlen = if left { mi } else { ni };
        let v: Vec<C> = a[p..p + vlen].to_vec();
        let c0 = ix(ic, jc, ldc);
        zlarf(left, mi, ni, &v, taui, &mut c[c0..], ldc);
        a[p] = aii;
    }
}

/// `ZROT(N, CX, INCX, CY, INCY, C, S)` on strided runs of one array.
fn zrot_within(a: &mut [C], x0: usize, y0: usize, step: usize, n: usize, c: f64, s: C) {
    for k in 0..n {
        let (px, py) = (x0 + k * step, y0 + k * step);
        let (cx, cy) = (a[px], a[py]);
        let stemp = c * cx + s * cy;
        a[py] = c * cy - s.conj() * cx;
        a[px] = stemp;
    }
}

/// `ZTREXC('V', N, T, LDT, Q, LDQ, IFST, ILST, INFO)`: moves the diagonal
/// entry at `IFST` of the upper triangular `T` to `ILST` by unitary
/// similarity, accumulated in `Q`.
pub fn ztrexc_v(n: usize, t: &mut [C], ldt: usize, q: &mut [C], ldq: usize, ifst: usize, ilst: usize) {
    if n <= 1 || ifst == ilst {
        return;
    }
    let ks: Vec<usize> = if ifst < ilst { (ifst..ilst).collect() } else { (ilst..ifst).rev().collect() };
    for k in ks {
        let t11 = t[ix(k, k, ldt)];
        let t22 = t[ix(k + 1, k + 1, ldt)];
        let (cs, sn, _) = zlartg(t[ix(k, k + 1, ldt)], t22 - t11);
        if k + 2 <= n {
            zrot_within(t, ix(k, k + 2, ldt), ix(k + 1, k + 2, ldt), ldt, n - k - 1, cs, sn);
        }
        zrot_within(t, ix(1, k, ldt), ix(1, k + 1, ldt), 1, k - 1, cs, sn.conj());
        t[ix(k, k, ldt)] = t22;
        t[ix(k + 1, k + 1, ldt)] = t11;
        zrot_within(q, ix(1, k, ldq), ix(1, k + 1, ldq), 1, n, cs, sn.conj());
    }
}

/// `ZTRSEN('N', 'V', SELECT, N, T, LDT, Q, LDQ, W, M, S, SEP, WORK,
/// LWORK, INFO)`: the selected eigenvalues moved to the leading block of
/// the Schur form, `Q` updated; `W` gets the reordered diagonal. Returns
/// `M`, the number selected.
pub fn ztrsen_nv(select: &[bool], n: usize, t: &mut [C], ldt: usize, q: &mut [C], ldq: usize, w: &mut [C]) -> usize {
    let m = select[..n].iter().filter(|&&s| s).count();
    if !(m == n || m == 0) {
        let mut ks = 0;
        for k in 1..=n {
            if select[k - 1] {
                ks += 1;
                if k != ks {
                    ztrexc_v(n, t, ldt, q, ldq, k, ks);
                }
            }
        }
    }
    for k in 1..=n {
        w[k - 1] = t[ix(k, k, ldt)];
    }
    m
}

/// `ZTREVC('Right', 'Select', SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM,
/// M, WORK, RWORK, INFO)`: the right eigenvectors of the upper triangular
/// `T` for the selected eigenvalues, not back-transformed, in columns `1..M`
/// of `VR` in increasing order of the selected index, each normalized to
/// unit largest `|re| + |im|` and zero below its index. Returns `M`.
pub fn ztrevc_rs(select: &[bool], n: usize, t: &mut [C], ldt: usize, vr: &mut [C], ldvr: usize) -> usize {
    let m = select[..n].iter().filter(|&&s| s).count();
    if n == 0 {
        return m;
    }
    let unfl = SAFMIN;
    let ulp = PREC;
    let smlnum = unfl * (n as f64 / ulp);
    let mut work = vec![ZERO; 2 * n];
    let mut rwork = vec![0.0f64; n];
    for i in 1..=n {
        work[i + n - 1] = t[ix(i, i, ldt)];
    }
    rwork[0] = 0.0;
    for j in 2..=n {
        let col: Vec<C> = (1..j).map(|i| t[ix(i, j, ldt)]).collect();
        rwork[j - 1] = super::zblas::dzasum(j - 1, &col);
    }
    let mut is = m;
    for ki in (1..=n).rev() {
        if !select[ki - 1] {
            continue;
        }
        let smin = (ulp * cabs1(t[ix(ki, ki, ldt)])).max(smlnum);
        work[0] = ONE;
        for k in 1..ki {
            work[k - 1] = -t[ix(k, ki, ldt)];
        }
        let tkk = t[ix(ki, ki, ldt)];
        for k in 1..ki {
            let p = ix(k, k, ldt);
            t[p] -= tkk;
            if cabs1(t[p]) < smin {
                t[p] = C::new(smin, 0.0);
            }
        }
        if ki > 1 {
            let scale = zlatrs_unn(ki - 1, t, ldt, &mut work[..ki - 1], &mut rwork[..ki - 1]);
            work[ki - 1] = C::new(scale, 0.0);
        }
        // Copy the vector x to VR and normalize.
        let mut col: Vec<C> = work[..ki].to_vec();
        let ii = super::zblas::izamax(ki, &col);
        let remax = 1.0 / cabs1(col[ii - 1]);
        zdscal(ki, remax, &mut col);
        for k in 1..=ki {
            vr[ix(k, is, ldvr)] = col[k - 1];
        }
        for k in ki + 1..=n {
            vr[ix(k, is, ldvr)] = ZERO;
        }
        for k in 1..ki {
            t[ix(k, k, ldt)] = work[k + n - 1];
        }
        is -= 1;
    }
    m
}
