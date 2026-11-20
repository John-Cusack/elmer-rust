//! The reference LAPACK routines ARPACK's symmetric path calls, in the
//! operation order of Netlib LAPACK 3.12 (the `liblapack.so.3` the pristine
//! Elmer build links). Only the argument combinations ARPACK uses are
//! ported. Arrays are column-major; indices in comments are Fortran's.

use super::blas::{dgemv_n, dgemv_t, dger, dnrm2, dscal, dswap};

/// `DLAMCH('E')`: relative machine epsilon (rounding), 2**-53.
pub const EPS: f64 = f64::EPSILON * 0.5;
/// `DLAMCH('P')`: `eps * base`.
pub const PREC: f64 = f64::EPSILON;
/// `DLAMCH('S')`: safe minimum, 2**-1022.
pub const SAFMIN: f64 = f64::MIN_POSITIVE;

/// `DLAPY2(x, y)`: `sqrt(x**2 + y**2)` without needless overflow.
pub fn dlapy2(x: f64, y: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if y.is_nan() {
        return y;
    }
    let (xa, ya) = (x.abs(), y.abs());
    let w = xa.max(ya);
    let z = xa.min(ya);
    if z == 0.0 || w > f64::MAX {
        w
    } else {
        w * (1.0 + (z / w) * (z / w)).sqrt()
    }
}

/// `DLARTG(f, g, c, s, r)` of LAPACK 3.10+ (`la_xlartg`): `[c s; -s c] [f; g] = [r; 0]`.
pub fn dlartg(f: f64, g: f64) -> (f64, f64, f64) {
    let safmin = f64::MIN_POSITIVE; // la_constants dsafmin = 2**-1022
    let safmax = 1.0 / safmin;
    let rtmin = safmin.sqrt();
    let rtmax = (safmax / 2.0).sqrt();
    let f1 = f.abs();
    let g1 = g.abs();
    if g == 0.0 {
        (1.0, 0.0, f)
    } else if f == 0.0 {
        (0.0, 1f64.copysign(g), g1)
    } else if f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax {
        let d = (f * f + g * g).sqrt();
        let c = f1 / d;
        let r = d.copysign(f);
        (c, g / r, r)
    } else {
        let u = safmax.min(safmin.max(f1).max(g1));
        let fs = f / u;
        let gs = g / u;
        let d = (fs * fs + gs * gs).sqrt();
        let c = fs.abs() / d;
        let r = d.copysign(f);
        let s = gs / r;
        (c, s, r * u)
    }
}

/// `DLASCL('G', 0, 0, cfrom, cto, m, n, a, lda, info)`: `A = A * (cto / cfrom)`
/// without over- or underflow, in steps when needed.
#[allow(clippy::too_many_arguments)]
pub fn dlascl_g(cfrom: f64, cto: f64, m: usize, n: usize, a: &mut [f64], lda: usize) {
    if n == 0 || m == 0 || cfrom == 0.0 || cfrom.is_nan() || cto.is_nan() {
        return;
    }
    let smlnum = SAFMIN;
    let bignum = 1.0 / smlnum;
    let mut cfromc = cfrom;
    let mut ctoc = cto;
    loop {
        let cfrom1 = cfromc * smlnum;
        let (mul, done);
        if cfrom1 == cfromc {
            // CFROMC is an inf.
            mul = ctoc / cfromc;
            done = true;
        } else {
            let cto1 = ctoc / bignum;
            if cto1 == ctoc {
                // CTOC is zero or inf.
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
            for v in &mut a[j * lda..j * lda + m] {
                *v *= mul;
            }
        }
        if done {
            return;
        }
    }
}

/// `DLASET('A', m, n, alpha, beta, a, lda)`: off-diagonal `alpha`, diagonal `beta`.
pub fn dlaset(m: usize, n: usize, alpha: f64, beta: f64, a: &mut [f64], lda: usize) {
    for j in 0..n {
        for i in 0..m {
            a[j * lda + i] = if i == j { beta } else { alpha };
        }
    }
}

/// `DLACPY('A', m, n, a, lda, b, ldb)`.
pub fn dlacpy(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        b[j * ldb..j * ldb + m].copy_from_slice(&a[j * lda..j * lda + m]);
    }
}

/// `MM` of DLARUV: the multipliers `a**(2**k)`-style table, 12-bit limbs.
const MM: [[i32; 4]; 128] = [
    [494, 322, 2508, 2549],
    [2637, 789, 3754, 1145],
    [255, 1440, 1766, 2253],
    [2008, 752, 3572, 305],
    [1253, 2859, 2893, 3301],
    [3344, 123, 307, 1065],
    [4084, 1848, 1297, 3133],
    [1739, 643, 3966, 2913],
    [3143, 2405, 758, 3285],
    [3468, 2638, 2598, 1241],
    [688, 2344, 3406, 1197],
    [1657, 46, 2922, 3729],
    [1238, 3814, 1038, 2501],
    [3166, 913, 2934, 1673],
    [1292, 3649, 2091, 541],
    [3422, 339, 2451, 2753],
    [1270, 3808, 1580, 949],
    [2016, 822, 1958, 2361],
    [154, 2832, 2055, 1165],
    [2862, 3078, 1507, 4081],
    [697, 3633, 1078, 2725],
    [1706, 2970, 3273, 3305],
    [491, 637, 17, 3069],
    [931, 2249, 854, 3617],
    [1444, 2081, 2916, 3733],
    [444, 4019, 3971, 409],
    [3577, 1478, 2889, 2157],
    [3944, 242, 3831, 1361],
    [2184, 481, 2621, 3973],
    [1661, 2075, 1541, 1865],
    [3482, 4058, 893, 2525],
    [657, 622, 736, 1409],
    [3023, 3376, 3992, 3445],
    [3618, 812, 787, 3577],
    [1267, 234, 2125, 77],
    [1828, 641, 2364, 3761],
    [164, 4005, 2460, 2149],
    [3798, 1122, 257, 1449],
    [3087, 3135, 1574, 3005],
    [2400, 2640, 3912, 225],
    [2870, 2302, 1216, 85],
    [3876, 40, 3248, 3673],
    [1905, 1832, 3401, 3117],
    [1593, 2247, 2124, 3089],
    [1797, 2034, 2762, 1349],
    [1234, 2637, 149, 2057],
    [3460, 1287, 2245, 413],
    [328, 1691, 166, 65],
    [2861, 496, 466, 1845],
    [1950, 1597, 4018, 697],
    [617, 2394, 1399, 3085],
    [2070, 2584, 190, 3441],
    [3331, 1843, 2879, 1573],
    [769, 336, 153, 3689],
    [1558, 1472, 2320, 2941],
    [2412, 2407, 18, 929],
    [2800, 433, 712, 533],
    [189, 2096, 2159, 2841],
    [287, 1761, 2318, 4077],
    [2045, 2810, 2091, 721],
    [1227, 566, 3443, 2821],
    [2838, 442, 1510, 2249],
    [209, 41, 449, 2397],
    [2770, 1238, 1956, 2817],
    [3654, 1086, 2201, 245],
    [3993, 603, 3137, 1913],
    [192, 840, 3399, 1997],
    [2253, 3168, 1321, 3121],
    [3491, 1499, 2271, 997],
    [2889, 1084, 3667, 1833],
    [2857, 3438, 2703, 2877],
    [2094, 2408, 629, 1633],
    [1818, 1589, 2365, 981],
    [688, 2391, 2431, 2009],
    [1407, 288, 1113, 941],
    [634, 26, 3922, 2449],
    [3231, 512, 2554, 197],
    [815, 1456, 184, 2441],
    [3524, 171, 2099, 285],
    [1914, 1677, 3228, 1473],
    [516, 2657, 4012, 2741],
    [164, 2270, 1921, 3129],
    [303, 2587, 3452, 909],
    [2144, 2961, 3901, 2801],
    [3480, 1970, 572, 421],
    [119, 1817, 3309, 4073],
    [3357, 676, 3171, 2813],
    [837, 1410, 817, 2337],
    [2826, 3723, 3039, 1429],
    [2332, 2803, 1696, 1177],
    [2089, 3185, 1256, 1901],
    [3780, 184, 3715, 81],
    [1700, 663, 2077, 1669],
    [3712, 499, 3019, 2633],
    [150, 3784, 1497, 2269],
    [2000, 1631, 1101, 129],
    [3375, 1925, 717, 1141],
    [1621, 3912, 51, 249],
    [3090, 1398, 981, 3917],
    [3765, 1349, 1978, 2481],
    [1149, 1441, 1813, 3941],
    [3146, 2224, 3881, 2217],
    [33, 2411, 76, 2749],
    [3082, 1907, 3846, 3041],
    [2741, 3192, 3694, 1877],
    [359, 2786, 1682, 345],
    [3316, 382, 124, 2861],
    [1749, 37, 1660, 1809],
    [185, 759, 3997, 3141],
    [2784, 2948, 479, 2825],
    [2202, 1862, 1141, 157],
    [2199, 3802, 886, 2881],
    [1364, 2423, 3514, 3637],
    [1244, 2051, 1301, 1465],
    [2020, 2295, 3604, 2829],
    [3160, 1332, 1888, 2161],
    [2785, 1832, 1836, 3365],
    [2772, 2405, 1990, 361],
    [1217, 3638, 2058, 2685],
    [1822, 3661, 692, 3745],
    [1245, 327, 1194, 2325],
    [2252, 3660, 20, 3609],
    [3904, 716, 3285, 3821],
    [2774, 1842, 2046, 3537],
    [997, 3987, 2107, 517],
    [2573, 1368, 3508, 3017],
    [1148, 1848, 3525, 2141],
    [545, 2366, 3801, 1537],
];

/// `DLARUV(iseed, n, x)`: up to 128 uniform (0, 1) numbers from the
/// multiplicative congruential generator of LAPACK, in 48-bit integer
/// arithmetic on four 12-bit limbs.
pub fn dlaruv(iseed: &mut [i32; 4], n: usize, x: &mut [f64]) {
    const IPW2: i32 = 4096;
    let r = 1.0 / f64::from(IPW2);
    let (mut i1, mut i2, mut i3, mut i4) = (iseed[0], iseed[1], iseed[2], iseed[3]);
    let (mut it1, mut it2, mut it3, mut it4) = (0, 0, 0, 0);
    for i in 0..n.min(128) {
        let m = MM[i];
        loop {
            it4 = i4 * m[3];
            it3 = it4 / IPW2;
            it4 -= IPW2 * it3;
            it3 = it3 + i3 * m[3] + i4 * m[2];
            it2 = it3 / IPW2;
            it3 -= IPW2 * it2;
            it2 = it2 + i2 * m[3] + i3 * m[2] + i4 * m[1];
            it1 = it2 / IPW2;
            it2 -= IPW2 * it1;
            it1 = it1 + i1 * m[3] + i2 * m[2] + i3 * m[1] + i4 * m[0];
            it1 %= IPW2;
            x[i] = r * (f64::from(it1) + r * (f64::from(it2) + r * (f64::from(it3) + r * f64::from(it4))));
            if x[i] != 1.0 {
                break;
            }
            // A 48-bit value whose first 53 bits are all one rounds to
            // exactly 1: draw again from the next seed (LAPACK 3.2+).
            i1 += 2;
            i2 += 2;
            i3 += 2;
            i4 += 2;
        }
    }
    *iseed = [it1, it2, it3, it4];
}

/// `DLARNV(idist, iseed, n, x)` for `idist` 1 (uniform (0, 1)) and 2
/// (uniform (-1, 1)), drawn 64 at a time.
pub fn dlarnv(idist: i32, iseed: &mut [i32; 4], n: usize, x: &mut [f64]) {
    let mut u = [0.0f64; 128];
    let mut iv = 0;
    while iv < n {
        let il = 64.min(n - iv);
        dlaruv(iseed, il, &mut u);
        for i in 0..il {
            x[iv + i] = match idist {
                1 => u[i],
                _ => 2.0 * u[i] - 1.0,
            };
        }
        iv += 64;
    }
}

/// `DLAE2(a, b, c)`: eigenvalues `(rt1, rt2)` of `[a b; b c]`, `|rt1| >= |rt2|`.
pub fn dlae2(a: f64, b: f64, c: f64) -> (f64, f64) {
    let sm = a + c;
    let df = a - c;
    let adf = df.abs();
    let tb = b + b;
    let ab = tb.abs();
    let (acmx, acmn) = if a.abs() > c.abs() { (a, c) } else { (c, a) };
    let rt = if adf > ab {
        adf * (1.0 + (ab / adf) * (ab / adf)).sqrt()
    } else if adf < ab {
        ab * (1.0 + (adf / ab) * (adf / ab)).sqrt()
    } else {
        ab * 2f64.sqrt()
    };
    if sm < 0.0 {
        let rt1 = 0.5 * (sm - rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b)
    } else if sm > 0.0 {
        let rt1 = 0.5 * (sm + rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b)
    } else {
        (0.5 * rt, -0.5 * rt)
    }
}

/// `DLAEV2(a, b, c)`: `(rt1, rt2, cs1, sn1)`, the eigen decomposition of `[a b; b c]`.
pub fn dlaev2(a: f64, b: f64, c: f64) -> (f64, f64, f64, f64) {
    let sm = a + c;
    let df = a - c;
    let adf = df.abs();
    let tb = b + b;
    let ab = tb.abs();
    let (acmx, acmn) = if a.abs() > c.abs() { (a, c) } else { (c, a) };
    let rt = if adf > ab {
        adf * (1.0 + (ab / adf) * (ab / adf)).sqrt()
    } else if adf < ab {
        ab * (1.0 + (adf / ab) * (adf / ab)).sqrt()
    } else {
        ab * 2f64.sqrt()
    };
    let (rt1, rt2, sgn1);
    if sm < 0.0 {
        rt1 = 0.5 * (sm - rt);
        sgn1 = -1;
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else if sm > 0.0 {
        rt1 = 0.5 * (sm + rt);
        sgn1 = 1;
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else {
        rt1 = 0.5 * rt;
        rt2 = -0.5 * rt;
        sgn1 = 1;
    }
    let (cs, sgn2) = if df >= 0.0 { (df + rt, 1) } else { (df - rt, -1) };
    let acs = cs.abs();
    let (mut cs1, mut sn1);
    if acs > ab {
        let ct = -tb / cs;
        sn1 = 1.0 / (1.0 + ct * ct).sqrt();
        cs1 = ct * sn1;
    } else if ab == 0.0 {
        cs1 = 1.0;
        sn1 = 0.0;
    } else {
        let tn = -cs / tb;
        cs1 = 1.0 / (1.0 + tn * tn).sqrt();
        sn1 = tn * cs1;
    }
    if sgn1 == sgn2 {
        let tn = cs1;
        cs1 = -sn1;
        sn1 = tn;
    }
    (rt1, rt2, cs1, sn1)
}

/// `DLANST('M' | 'I', n, d, e)`: max abs entry, or the infinity (= one) norm
/// of the symmetric tridiagonal matrix.
pub fn dlanst(norm: u8, n: usize, d: &[f64], e: &[f64]) -> f64 {
    if n == 0 {
        return 0.0;
    }
    // MAX that keeps a NaN (`ANORM.LT.SUM .OR. DISNAN(SUM)`).
    let upd = |anorm: f64, sum: f64| if anorm < sum || sum.is_nan() { sum } else { anorm };
    match norm {
        b'M' => {
            let mut anorm = d[n - 1].abs();
            for i in 0..n - 1 {
                anorm = upd(anorm, d[i].abs());
                anorm = upd(anorm, e[i].abs());
            }
            anorm
        }
        _ => {
            if n == 1 {
                return d[0].abs();
            }
            let mut anorm = d[0].abs() + e[0].abs();
            anorm = upd(anorm, e[n - 2].abs() + d[n - 1].abs());
            for i in 1..n - 1 {
                anorm = upd(anorm, d[i].abs() + e[i].abs() + e[i - 1].abs());
            }
            anorm
        }
    }
}

/// `DLASR('R', 'V', direct, m, n, c, s, a, lda)`: `A = A P^T` with the plane
/// rotations `(c(j), s(j))` acting on columns `j, j+1`, forward (`b'F'`) or
/// backward (`b'B'`).
#[allow(clippy::too_many_arguments)]
pub fn dlasr_rv(direct: u8, m: usize, n: usize, c: &[f64], s: &[f64], a: &mut [f64], lda: usize) {
    if m == 0 || n == 0 {
        return;
    }
    let mut rot = |j: usize| {
        let (ct, st) = (c[j], s[j]);
        if ct != 1.0 || st != 0.0 {
            for i in 0..m {
                let temp = a[(j + 1) * lda + i];
                let aij = a[j * lda + i];
                a[(j + 1) * lda + i] = ct * temp - st * aij;
                a[j * lda + i] = st * temp + ct * aij;
            }
        }
    };
    if direct == b'F' {
        for j in 0..n - 1 {
            rot(j);
        }
    } else {
        for j in (0..n - 1).rev() {
            rot(j);
        }
    }
}

/// `DLASRT('I', n, d)`: ascending order (values only, so any correct sort
/// gives the same result).
pub fn dlasrt_inc(d: &mut [f64]) {
    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}

/// `DSTEQR('I', n, d, e, z, ldz, work, info)`: all eigenvalues and
/// eigenvectors of the symmetric tridiagonal matrix (diagonal `d`,
/// off-diagonal `e`) by implicit QL/QR; `z` starts as the identity. `work`
/// holds `2*(n-1)`. Returns `info` (0, or the number of unconverged `e`).
#[allow(clippy::too_many_arguments)]
pub fn dsteqr_i(n: usize, d: &mut [f64], e: &mut [f64], z: &mut [f64], ldz: usize, work: &mut [f64]) -> usize {
    if n == 0 {
        return 0;
    }
    if n == 1 {
        z[0] = 1.0;
        return 0;
    }
    dlaset(n, n, 0.0, 1.0, z, ldz);
    steqr_core(n, d, e, z, n, ldz, work)
}

/// ARPACK's `DSTQRB(n, d, e, z, work, info)`: `DSTEQR` that keeps only the
/// last row of the eigenvector matrix, `z` (length `n`, starting as `e_n`).
pub fn dstqrb(n: usize, d: &mut [f64], e: &mut [f64], z: &mut [f64], work: &mut [f64]) -> usize {
    if n == 0 {
        return 0;
    }
    if n == 1 {
        z[0] = 1.0;
        return 0;
    }
    z[..n - 1].iter_mut().for_each(|v| *v = 0.0);
    z[n - 1] = 1.0;
    steqr_core(n, d, e, z, 1, 1, work)
}

/// The implicit QL/QR iteration of `DSTEQR` / `DSTQRB` with the rotations
/// applied to the `nrz` rows of `z` (leading dimension `ldz`).
#[allow(clippy::too_many_arguments)]
fn steqr_core(n: usize, d: &mut [f64], e: &mut [f64], z: &mut [f64], nrz: usize, ldz: usize, work: &mut [f64]) -> usize {
    let eps = EPS;
    let eps2 = eps * eps;
    let safmin = SAFMIN;
    let safmax = 1.0 / safmin;
    let ssfmax = safmax.sqrt() / 3.0;
    let ssfmin = safmin.sqrt() / eps2;
    let nmaxit = n * 30;
    let mut jtot = 0usize;
    // Fortran 1-based indices throughout; d(i) = d[i-1].
    let mut l1 = 1usize;
    let nm1 = n - 1;
    macro_rules! dd { ($i:expr) => { d[$i - 1] } }
    macro_rules! ee { ($i:expr) => { e[$i - 1] } }
    let mut info = 0usize;
    'outer: loop {
        // 10
        if l1 > n {
            break 'outer; // 160
        }
        if l1 > 1 {
            ee!(l1 - 1) = 0.0;
        }
        let mut m = n;
        if l1 <= nm1 {
            let mut mm = l1;
            loop {
                let tst = ee!(mm).abs();
                if tst == 0.0 {
                    m = mm;
                    break;
                }
                if tst <= (dd!(mm).abs().sqrt() * dd!(mm + 1).abs().sqrt()) * eps {
                    ee!(mm) = 0.0;
                    m = mm;
                    break;
                }
                if mm == nm1 {
                    m = n;
                    break;
                }
                mm += 1;
            }
        }
        // 30
        let mut l = l1;
        let lsv = l;
        let mut lend = m;
        let lendsv = lend;
        l1 = m + 1;
        if lend == l {
            continue 'outer;
        }
        let anorm = dlanst(b'I', lend - l + 1, &d[l - 1..], &e[l - 1..]);
        let mut iscale = 0;
        if anorm == 0.0 {
            continue 'outer;
        }
        if anorm > ssfmax {
            iscale = 1;
            dlascl_g(anorm, ssfmax, lend - l + 1, 1, &mut d[l - 1..], n);
            dlascl_g(anorm, ssfmax, lend - l, 1, &mut e[l - 1..], n);
        } else if anorm < ssfmin {
            iscale = 2;
            dlascl_g(anorm, ssfmin, lend - l + 1, 1, &mut d[l - 1..], n);
            dlascl_g(anorm, ssfmin, lend - l, 1, &mut e[l - 1..], n);
        }
        if dd!(lend).abs() < dd!(l).abs() {
            lend = lsv;
            l = lendsv;
        }
        if lend > l {
            // QL iteration.
            loop {
                // 40
                let mut m = lend;
                if l != lend {
                    let lendm1 = lend - 1;
                    for mm in l..=lendm1 {
                        let tst = ee!(mm).abs() * ee!(mm).abs();
                        if tst <= (eps2 * dd!(mm).abs()) * dd!(mm + 1).abs() + safmin {
                            m = mm;
                            break;
                        }
                    }
                }
                // 60
                if m < lend {
                    ee!(m) = 0.0;
                }
                let mut p = dd!(l);
                if m == l {
                    // 80
                    dd!(l) = p;
                    l += 1;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if m == l + 1 {
                    let (rt1, rt2, c, s) = dlaev2(dd!(l), ee!(l), dd!(l + 1));
                    work[l - 1] = c;
                    work[n - 1 + l - 1] = s;
                    let (wc, ws) = work.split_at_mut(n - 1);
                    dlasr_rv(b'B', nrz, 2, &wc[l - 1..], &ws[l - 1..], &mut z[(l - 1) * ldz..], ldz);
                    dd!(l) = rt1;
                    dd!(l + 1) = rt2;
                    ee!(l) = 0.0;
                    l += 2;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let mut g = (dd!(l + 1) - p) / (2.0 * ee!(l));
                let mut r = dlapy2(g, 1.0);
                g = dd!(m) - p + (ee!(l) / (g + r.copysign(g)));
                let mut s = 1.0;
                let mut c = 1.0;
                p = 0.0;
                let mm1 = m - 1;
                let mut i = mm1;
                loop {
                    let f = s * ee!(i);
                    let b = c * ee!(i);
                    let (cc, ss, rr) = dlartg(g, f);
                    c = cc;
                    s = ss;
                    r = rr;
                    if i != m - 1 {
                        ee!(i + 1) = r;
                    }
                    g = dd!(i + 1) - p;
                    r = (dd!(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    dd!(i + 1) = g + p;
                    g = c * r - b;
                    work[i - 1] = c;
                    work[n - 1 + i - 1] = -s;
                    if i == l {
                        break;
                    }
                    i -= 1;
                }
                let mm = m - l + 1;
                let (wc, ws) = work.split_at_mut(n - 1);
                dlasr_rv(b'B', nrz, mm, &wc[l - 1..], &ws[l - 1..], &mut z[(l - 1) * ldz..], ldz);
                dd!(l) -= p;
                ee!(l) = g;
            }
        } else {
            // QR iteration.
            loop {
                // 90
                let mut m = lend;
                if l != lend {
                    let lendp1 = lend + 1;
                    let mut mm = l;
                    loop {
                        let tst = ee!(mm - 1).abs() * ee!(mm - 1).abs();
                        if tst <= (eps2 * dd!(mm).abs()) * dd!(mm - 1).abs() + safmin {
                            m = mm;
                            break;
                        }
                        if mm == lendp1 {
                            break;
                        }
                        mm -= 1;
                    }
                }
                // 110
                if m > lend {
                    ee!(m - 1) = 0.0;
                }
                let mut p = dd!(l);
                if m == l {
                    // 130
                    dd!(l) = p;
                    if l == lend {
                        break;
                    }
                    l -= 1;
                    if l >= lend {
                        continue;
                    }
                    break;
                }
                if m + 1 == l {
                    let (rt1, rt2, c, s) = dlaev2(dd!(l - 1), ee!(l - 1), dd!(l));
                    work[m - 1] = c;
                    work[n - 1 + m - 1] = s;
                    let (wc, ws) = work.split_at_mut(n - 1);
                    dlasr_rv(b'F', nrz, 2, &wc[m - 1..], &ws[m - 1..], &mut z[(l - 2) * ldz..], ldz);
                    dd!(l - 1) = rt1;
                    dd!(l) = rt2;
                    ee!(l - 1) = 0.0;
                    if l < lend + 2 {
                        break;
                    }
                    l -= 2;
                    continue;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let mut g = (dd!(l - 1) - p) / (2.0 * ee!(l - 1));
                let mut r = dlapy2(g, 1.0);
                g = dd!(m) - p + (ee!(l - 1) / (g + r.copysign(g)));
                let mut s = 1.0;
                let mut c = 1.0;
                p = 0.0;
                let lm1 = l - 1;
                for i in m..=lm1 {
                    let f = s * ee!(i);
                    let b = c * ee!(i);
                    let (cc, ss, rr) = dlartg(g, f);
                    c = cc;
                    s = ss;
                    r = rr;
                    if i != m {
                        ee!(i - 1) = r;
                    }
                    g = dd!(i) - p;
                    r = (dd!(i + 1) - g) * s + 2.0 * c * b;
                    p = s * r;
                    dd!(i) = g + p;
                    g = c * r - b;
                    work[i - 1] = c;
                    work[n - 1 + i - 1] = s;
                }
                let mm = l - m + 1;
                let (wc, ws) = work.split_at_mut(n - 1);
                dlasr_rv(b'F', nrz, mm, &wc[m - 1..], &ws[m - 1..], &mut z[(m - 1) * ldz..], ldz);
                dd!(l) -= p;
                ee!(lm1) = g;
            }
        }
        // 140
        if iscale == 1 {
            dlascl_g(ssfmax, anorm, lendsv - lsv + 1, 1, &mut d[lsv - 1..], n);
            dlascl_g(ssfmax, anorm, lendsv - lsv, 1, &mut e[lsv - 1..], n);
        } else if iscale == 2 {
            dlascl_g(ssfmin, anorm, lendsv - lsv + 1, 1, &mut d[lsv - 1..], n);
            dlascl_g(ssfmin, anorm, lendsv - lsv, 1, &mut e[lsv - 1..], n);
        }
        if jtot < nmaxit {
            continue 'outer;
        }
        for i in 1..n {
            if ee!(i) != 0.0 {
                info += 1;
            }
        }
        return info;
    }
    // 160: selection sort with the vectors.
    for ii in 2..=n {
        let i = ii - 1;
        let mut k = i;
        let mut p = dd!(i);
        for j in ii..=n {
            if dd!(j) < p {
                k = j;
                p = dd!(j);
            }
        }
        if k != i {
            dd!(k) = dd!(i);
            dd!(i) = p;
            let (lo, hi) = z.split_at_mut((k - 1) * ldz);
            dswap(nrz, &mut lo[(i - 1) * ldz..], hi);
        }
    }
    info
}

/// `DLARFG(n, alpha, x, 1, tau)`: the elementary reflector `H` with
/// `H [alpha; x] = [beta; 0]`. Returns `(beta -> new alpha, tau)`; `x` is
/// overwritten with `v(2:n)`.
pub fn dlarfg(n: usize, alpha: f64, x: &mut [f64]) -> (f64, f64) {
    if n <= 1 {
        return (alpha, 0.0);
    }
    let mut xnorm = dnrm2(n - 1, x);
    if xnorm == 0.0 {
        return (alpha, 0.0);
    }
    let mut alpha = alpha;
    let mut beta = -dlapy2(alpha, xnorm).copysign(alpha);
    let safmin = SAFMIN / EPS;
    let mut knt = 0;
    if beta.abs() < safmin {
        let rsafmn = 1.0 / safmin;
        loop {
            knt += 1;
            dscal(n - 1, rsafmn, x);
            beta *= rsafmn;
            alpha *= rsafmn;
            if !(beta.abs() < safmin && knt < 20) {
                break;
            }
        }
        xnorm = dnrm2(n - 1, x);
        beta = -dlapy2(alpha, xnorm).copysign(alpha);
    }
    let tau = (beta - alpha) / beta;
    dscal(n - 1, 1.0 / (alpha - beta), x);
    for _ in 0..knt {
        beta *= safmin;
    }
    (beta, tau)
}

/// `DLARF('L' | 'R', m, n, v, 1, tau, c, ldc, work)`: `C = H C` or `C H`
/// with `H = I - tau v v^T`.
#[allow(clippy::too_many_arguments)]
pub fn dlarf(left: bool, m: usize, n: usize, v: &[f64], tau: f64, c: &mut [f64], ldc: usize, work: &mut [f64]) {
    if tau == 0.0 {
        return;
    }
    if left {
        dgemv_t(m, n, 1.0, c, ldc, v, 0.0, work);
        dger(m, n, -tau, v, work, c, ldc);
    } else {
        dgemv_n(m, n, 1.0, c, ldc, v, 0.0, work);
        dger(m, n, -tau, work, v, c, ldc);
    }
}

/// `DGEQR2(m, n, a, lda, tau, work, info)`: unblocked QR factorisation.
pub fn dgeqr2(m: usize, n: usize, a: &mut [f64], lda: usize, tau: &mut [f64], work: &mut [f64]) {
    let k = m.min(n);
    for i in 0..k {
        let aii = a[i * lda + i];
        let xs = (i + 1).min(m - 1);
        let (alpha_new, t) = {
            let (head, tail) = a.split_at_mut(i * lda + xs);
            let _ = head;
            dlarfg(m - i, aii, &mut tail[..m - xs])
        };
        a[i * lda + i] = alpha_new;
        tau[i] = t;
        if i + 1 < n {
            let aii = a[i * lda + i];
            a[i * lda + i] = 1.0;
            let v: Vec<f64> = a[i * lda + i..i * lda + m].to_vec();
            dlarf(true, m - i, n - i - 1, &v, t, &mut a[(i + 1) * lda + i..], lda, work);
            a[i * lda + i] = aii;
        }
    }
}

/// `DORM2R(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)`: `C = Q C`,
/// `Q^T C`, `C Q` or `C Q^T` (`left`, `trans`), with `Q = H(1) .. H(k)`.
#[allow(clippy::too_many_arguments)]
pub fn dorm2r(left: bool, trans: bool, m: usize, n: usize, k: usize, a: &mut [f64], lda: usize, tau: &[f64], c: &mut [f64], ldc: usize, work: &mut [f64]) {
    if m == 0 || n == 0 || k == 0 {
        return;
    }
    let forward = (left && trans) || (!left && !trans);
    let order: Vec<usize> = if forward { (0..k).collect() } else { (0..k).rev().collect() };
    for i in order {
        let (mi, ni, ic, jc) = if left { (m - i, n, i, 0) } else { (m, n - i, 0, i) };
        let aii = a[i * lda + i];
        a[i * lda + i] = 1.0;
        let len = if left { mi } else { ni };
        let v: Vec<f64> = a[i * lda + i..i * lda + i + len].to_vec();
        dlarf(left, mi, ni, &v, tau[i], &mut c[jc * ldc + ic..], ldc, work);
        a[i * lda + i] = aii;
    }
}

/// [`dorm2r`] without transposition.
#[allow(clippy::too_many_arguments)]
pub fn dorm2r_n(left: bool, m: usize, n: usize, k: usize, a: &mut [f64], lda: usize, tau: &[f64], c: &mut [f64], ldc: usize, work: &mut [f64]) {
    dorm2r(left, false, m, n, k, a, lda, tau, c, ldc, work)
}
