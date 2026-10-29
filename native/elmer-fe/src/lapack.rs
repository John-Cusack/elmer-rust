//! The `DSYGV( 1, 'N', 'U', ... )` call of `StabParam`, transliterated from
//! reference LAPACK/BLAS 3.12.0 (the `liblapack3`/`libblas3` the pristine
//! binary links) for the orders `StabParam` uses (`n <= 31`, so every
//! blocked routine takes its unblocked path):
//! `DPOTRF` -> `DPOTRF2`, `DSYGST` -> `DSYGS2`, `DSYEV` -> `DSYTRD` ->
//! `DSYTD2` (with `DLARFG`, `DNRM2`, `DLAPY2`, `DSYMV`, `DDOT`, `DAXPY`,
//! `DSYR2`) and `DSTERF` (with `DLAE2`). Matrices are column-major with
//! leading dimension `n`, 0-based: `A(i,j)` is `a[i + j*n]`.
//!
//! The scaling branches (`DSYEV` when `max|A|` leaves
//! `[2^-485, 2^485]`, `DSTERF` when a block norm leaves
//! `[2^-405, 2^511/3]`) are not transliterated; those inputs return `None`.

const HALF: f64 = 0.5;

/// `DLAMCH('S')`, `DLAMCH('E')` and `DLAMCH('P')`.
const SAFMIN: f64 = f64::MIN_POSITIVE;
const EPS_E: f64 = f64::EPSILON * 0.5;
const EPS_P: f64 = f64::EPSILON;

/// Largest eigenvalue of `A x = lambda B x` (`A` symmetric, `B` symmetric
/// positive definite, upper triangles read), i.e. `W(N)` of
/// `DSYGV( 1, 'N', 'U', N, A, N, B, N, W, ... )`. `None` when `B` is not
/// positive definite, `DSTERF` fails to converge, or a scaling branch is hit.
pub(crate) fn dsygv_max(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<f64> {
    debug_assert!(n < 32, "unblocked LAPACK paths only");
    if n == 0 {
        return None;
    }
    dpotrf2(&mut b, 0, n, n)?;
    dsygs2(&mut a, &b, n);
    let w = dsyev_values(&mut a, n)?;
    w.into_iter().reduce(f64::max)
}

/// `DPOTRF2( 'U', n, A(base), lda )`, recursive.
fn dpotrf2(a: &mut [f64], base: usize, n: usize, lda: usize) -> Option<()> {
    if n == 1 {
        if a[base] <= 0.0 || a[base].is_nan() {
            return None;
        }
        a[base] = a[base].sqrt();
        return Some(());
    }
    let n1 = n / 2;
    let n2 = n - n1;
    dpotrf2(a, base, n1, lda)?;
    // DTRSM( 'L', 'U', 'T', 'N', n1, n2, ONE, A11, A12 )
    let b12 = base + n1 * lda;
    for j in 0..n2 {
        for i in 0..n1 {
            let mut temp = 1.0 * a[b12 + i + j * lda];
            for k in 0..i {
                temp -= a[base + k + i * lda] * a[b12 + k + j * lda];
            }
            temp /= a[base + i + i * lda];
            a[b12 + i + j * lda] = temp;
        }
    }
    // DSYRK( 'U', 'T', n2, n1, -ONE, A12, ONE, A22 )
    let b22 = base + n1 + n1 * lda;
    for j in 0..n2 {
        for i in 0..=j {
            let mut temp = 0.0;
            for l in 0..n1 {
                temp += a[b12 + l + i * lda] * a[b12 + l + j * lda];
            }
            let c = &mut a[b22 + i + j * lda];
            // ALPHA*TEMP + BETA*C with ALPHA = -1, BETA = 1: (-temp) + c rounds as c - temp.
            *c -= temp;
        }
    }
    dpotrf2(a, b22, n2, lda)
}

/// `DSYGS2( 1, 'U', n, A, n, B, n )`: `A := inv(U^T) A inv(U)`.
fn dsygs2(a: &mut [f64], b: &[f64], n: usize) {
    let ix = |i: usize, j: usize| i + j * n;
    for k in 0..n {
        let bkk = b[ix(k, k)];
        let akk = a[ix(k, k)] / (bkk * bkk);
        a[ix(k, k)] = akk;
        if k + 1 < n {
            // DSCAL( N-K, ONE/BKK, A(K,K+1), LDA )
            let r = 1.0 / bkk;
            for j in k + 1..n {
                a[ix(k, j)] *= r;
            }
            let ct = -HALF * akk;
            daxpy_row(a, b, k, n, ct);
            // DSYR2( 'U', N-K, -ONE, A(K,K+1), LDA, B(K,K+1), LDB, A(K+1,K+1), LDA )
            let alpha = -1.0;
            for j in k + 1..n {
                let (xj, yj) = (a[ix(k, j)], b[ix(k, j)]);
                if xj != 0.0 || yj != 0.0 {
                    let temp1 = alpha * yj;
                    let temp2 = alpha * xj;
                    for i in k + 1..=j {
                        a[ix(i, j)] = a[ix(i, j)] + a[ix(k, i)] * temp1 + b[ix(k, i)] * temp2;
                    }
                }
            }
            daxpy_row(a, b, k, n, ct);
            // DTRSV( 'U', 'T', 'N', N-K, B(K+1,K+1), LDB, A(K,K+1), LDA )
            for j in k + 1..n {
                let mut temp = a[ix(k, j)];
                for i in k + 1..j {
                    temp -= b[ix(i, j)] * a[ix(k, i)];
                }
                temp /= b[ix(j, j)];
                a[ix(k, j)] = temp;
            }
        }
    }
}

/// `DAXPY( N-K, CT, B(K,K+1), LDB, A(K,K+1), LDA )`.
fn daxpy_row(a: &mut [f64], b: &[f64], k: usize, n: usize, ct: f64) {
    if ct == 0.0 {
        return;
    }
    for j in k + 1..n {
        a[k + j * n] += ct * b[k + j * n];
    }
}

/// `DSYEV( 'N', 'U', n, A, ... )`: the eigenvalues, ascending order not
/// needed by the caller.
fn dsyev_values(a: &mut [f64], n: usize) -> Option<Vec<f64>> {
    if n == 1 {
        return Some(vec![a[0]]);
    }
    let smlnum = SAFMIN / EPS_P;
    let bignum = 1.0 / smlnum;
    let (rmin, rmax) = (smlnum.sqrt(), bignum.sqrt());
    // DLANSY( 'M', 'U', ... )
    let mut anrm = 0.0f64;
    for j in 0..n {
        for i in 0..=j {
            let v = a[i + j * n].abs();
            if anrm < v || v.is_nan() {
                anrm = v;
            }
        }
    }
    if (anrm > 0.0 && anrm < rmin) || anrm > rmax {
        return None;
    }
    let (d, e) = dsytd2(a, n);
    dsterf(d, e)
}

/// `DSYTD2( 'U', n, A, n, D, E, TAU )`: returns `(D, E)`.
fn dsytd2(a: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let ix = |i: usize, j: usize| i + j * n;
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n - 1];
    let mut tau = vec![0.0; n];
    // Fortran I = N-1 .. 1; here i = I-1, the Householder vector is
    // column I+1 (0-based i+1), rows 1..I (0-based 0..=i).
    for i in (0..n - 1).rev() {
        let len = i + 1;
        let col = i + 1;
        let taui = dlarfg(a, len, ix(i, col), ix(0, col));
        e[i] = a[ix(i, col)];
        if taui != 0.0 {
            a[ix(i, col)] = 1.0;
            // DSYMV( 'U', I, TAUI, A, LDA, A(1,I+1), 1, ZERO, TAU, 1 )
            for t in tau.iter_mut().take(len) {
                *t = 0.0;
            }
            for j in 0..len {
                let temp1 = taui * a[ix(j, col)];
                let mut temp2 = 0.0;
                for r in 0..j {
                    tau[r] += temp1 * a[ix(r, j)];
                    temp2 += a[ix(r, j)] * a[ix(r, col)];
                }
                tau[j] = tau[j] + temp1 * a[ix(j, j)] + taui * temp2;
            }
            // ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A(1,I+1), 1 )
            let mut dot = 0.0;
            for r in 0..len {
                dot += tau[r] * a[ix(r, col)];
            }
            let alpha = -HALF * taui * dot;
            // DAXPY( I, ALPHA, A(1,I+1), 1, TAU, 1 )
            if alpha != 0.0 {
                for r in 0..len {
                    tau[r] += alpha * a[ix(r, col)];
                }
            }
            // DSYR2( 'U', I, -ONE, A(1,I+1), 1, TAU, 1, A, LDA )
            for j in 0..len {
                let (xj, yj) = (a[ix(j, col)], tau[j]);
                if xj != 0.0 || yj != 0.0 {
                    let temp1 = -yj;
                    let temp2 = -xj;
                    for r in 0..=j {
                        a[ix(r, j)] = a[ix(r, j)] + a[ix(r, col)] * temp1 + tau[r] * temp2;
                    }
                }
            }
            a[ix(i, col)] = e[i];
        }
        d[i + 1] = a[ix(i + 1, i + 1)];
        tau[i] = taui;
    }
    d[0] = a[0];
    (d, e)
}

/// `DLARFG( n, ALPHA = a[alpha], X = a[x..x+n-1], 1, TAU )`; returns `TAU`.
fn dlarfg(a: &mut [f64], n: usize, alpha_at: usize, x: usize) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    let m = n - 1;
    let mut xnorm = dnrm2(&a[x..x + m]);
    if xnorm == 0.0 {
        return 0.0;
    }
    let mut alpha = a[alpha_at];
    let mut beta = -dlapy2(alpha, xnorm).copysign(alpha);
    let safmin = SAFMIN / EPS_E;
    let mut knt = 0;
    if beta.abs() < safmin {
        let rsafmn = 1.0 / safmin;
        loop {
            knt += 1;
            for v in &mut a[x..x + m] {
                *v *= rsafmn;
            }
            beta *= rsafmn;
            alpha *= rsafmn;
            if !(beta.abs() < safmin && knt < 20) {
                break;
            }
        }
        xnorm = dnrm2(&a[x..x + m]);
        beta = -dlapy2(alpha, xnorm).copysign(alpha);
    }
    let tau = (beta - alpha) / beta;
    let r = 1.0 / (alpha - beta);
    for v in &mut a[x..x + m] {
        *v *= r;
    }
    for _ in 0..knt {
        beta *= SAFMIN;
    }
    a[alpha_at] = beta;
    tau
}

/// `DNRM2` (BLAS 3.10+ Blue's-algorithm version, dnrm2.f90).
fn dnrm2(x: &[f64]) -> f64 {
    let tsml = 2f64.powi(-511);
    let tbig = 2f64.powi(486);
    let ssml = 2f64.powi(537);
    let sbig = 2f64.powi(-538);
    if x.is_empty() {
        return 0.0;
    }
    let (mut asml, mut amed, mut abig) = (0.0f64, 0.0f64, 0.0f64);
    let mut notbig = true;
    for &v in x {
        let ax = v.abs();
        if ax > tbig {
            let t = ax * sbig;
            abig += t * t;
            notbig = false;
        } else if ax < tsml {
            if notbig {
                let t = ax * ssml;
                asml += t * t;
            }
        } else {
            amed += ax * ax;
        }
    }
    let (scl, sumsq);
    if abig > 0.0 {
        if amed > 0.0 || amed > f64::MAX || amed.is_nan() {
            abig += (amed * sbig) * sbig;
        }
        scl = 1.0 / sbig;
        sumsq = abig;
    } else if asml > 0.0 {
        if amed > 0.0 || amed > f64::MAX || amed.is_nan() {
            let amed = amed.sqrt();
            let asml = asml.sqrt() / ssml;
            let (ymin, ymax) = if asml > amed {
                (amed, asml)
            } else {
                (asml, amed)
            };
            scl = 1.0;
            let q = ymin / ymax;
            sumsq = ymax * ymax * (1.0 + q * q);
        } else {
            scl = 1.0 / ssml;
            sumsq = asml;
        }
    } else {
        scl = 1.0;
        sumsq = amed;
    }
    scl * sumsq.sqrt()
}

/// `DLAPY2( x, y )` = `sqrt(x^2 + y^2)` avoiding overflow.
fn dlapy2(x: f64, y: f64) -> f64 {
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
        let q = z / w;
        w * (1.0 + q * q).sqrt()
    }
}

/// `DLAE2( A, B, C, RT1, RT2 )`: eigenvalues of `[[a, b], [b, c]]`.
fn dlae2(a: f64, b: f64, c: f64) -> (f64, f64) {
    let sm = a + c;
    let df = a - c;
    let adf = df.abs();
    let tb = b + b;
    let ab = tb.abs();
    let (acmx, acmn) = if a.abs() > c.abs() { (a, c) } else { (c, a) };
    let rt = if adf > ab {
        let q = ab / adf;
        adf * (1.0 + q * q).sqrt()
    } else if adf < ab {
        let q = adf / ab;
        ab * (1.0 + q * q).sqrt()
    } else {
        ab * 2f64.sqrt()
    };
    if sm < 0.0 {
        let rt1 = HALF * (sm - rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b)
    } else if sm > 0.0 {
        let rt1 = HALF * (sm + rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b)
    } else {
        (HALF * rt, -HALF * rt)
    }
}

/// `DSTERF( N, D, E, INFO )`: eigenvalues of the symmetric tridiagonal
/// matrix by the Pal-Walker-Kahan QL/QR iteration (sorting omitted).
#[allow(clippy::many_single_char_names)]
fn dsterf(d0: Vec<f64>, e0: Vec<f64>) -> Option<Vec<f64>> {
    const MAXIT: usize = 30;
    let n = d0.len();
    if n <= 1 {
        return Some(d0);
    }
    // 1-based copies: d[1..=n], e[1..=n-1].
    let mut d = vec![0.0; n + 1];
    let mut e = vec![0.0; n + 1];
    d[1..].copy_from_slice(&d0);
    e[1..n].copy_from_slice(&e0);

    let eps = EPS_E;
    let eps2 = eps * eps;
    let safmax = 1.0 / SAFMIN;
    let ssfmax = safmax.sqrt() / 3.0;
    let ssfmin = SAFMIN.sqrt() / eps2;
    let nmaxit = n * MAXIT;
    let mut jtot = 0;
    let mut l1 = 1;

    'outer: loop {
        // 10
        if l1 > n {
            break;
        }
        if l1 > 1 {
            e[l1 - 1] = 0.0;
        }
        let mut m = n;
        for mm in l1..n {
            if e[mm].abs() <= (d[mm].abs().sqrt() * d[mm + 1].abs().sqrt()) * eps {
                e[mm] = 0.0;
                m = mm;
                break;
            }
        }
        let mut l = l1;
        let lsv = l;
        let mut lend = m;
        let lendsv = lend;
        l1 = m + 1;
        if lend == l {
            continue 'outer;
        }
        // DLANST( 'M', ... )
        let mut anorm = d[lend].abs();
        for i in l..lend {
            let s = d[i].abs();
            if anorm < s || s.is_nan() {
                anorm = s;
            }
            let s = e[i].abs();
            if anorm < s || s.is_nan() {
                anorm = s;
            }
        }
        if anorm == 0.0 {
            continue 'outer;
        }
        if anorm > ssfmax || anorm < ssfmin {
            return None;
        }
        for v in &mut e[l..lend] {
            *v *= *v;
        }
        if d[lend].abs() < d[l].abs() {
            lend = lsv;
            l = lendsv;
        }
        if lend >= l {
            // QL iteration.
            'ql: loop {
                // 50
                let mut m = lend;
                if l != lend {
                    for mm in l..lend {
                        if e[mm].abs() <= eps2 * (d[mm] * d[mm + 1]).abs() {
                            m = mm;
                            break;
                        }
                    }
                }
                if m < lend {
                    e[m] = 0.0;
                }
                let mut p = d[l];
                if m == l {
                    // 90
                    d[l] = p;
                    l += 1;
                    if l <= lend {
                        continue 'ql;
                    }
                    break 'ql;
                }
                if m == l + 1 {
                    let rte = e[l].sqrt();
                    let (rt1, rt2) = dlae2(d[l], rte, d[l + 1]);
                    d[l] = rt1;
                    d[l + 1] = rt2;
                    e[l] = 0.0;
                    l += 2;
                    if l <= lend {
                        continue 'ql;
                    }
                    break 'ql;
                }
                if jtot == nmaxit {
                    break 'ql;
                }
                jtot += 1;
                let rte = e[l].sqrt();
                let mut sigma = (d[l + 1] - p) / (2.0 * rte);
                let r = dlapy2(sigma, 1.0);
                sigma = p - (rte / (sigma + r.copysign(sigma)));
                let mut c = 1.0;
                let mut s = 0.0;
                let mut gamma = d[m] - sigma;
                p = gamma * gamma;
                for i in (l..m).rev() {
                    let bb = e[i];
                    let r = p + bb;
                    if i != m - 1 {
                        e[i + 1] = s * r;
                    }
                    let oldc = c;
                    c = p / r;
                    s = bb / r;
                    let oldgam = gamma;
                    let alpha = d[i];
                    gamma = c * (alpha - sigma) - s * oldgam;
                    d[i + 1] = oldgam + (alpha - gamma);
                    if c != 0.0 {
                        p = (gamma * gamma) / c;
                    } else {
                        p = oldc * bb;
                    }
                }
                e[l] = s * p;
                d[l] = sigma + gamma;
            }
        } else {
            // QR iteration.
            'qr: loop {
                // 100
                let mut m = lend;
                let mut mm = l;
                while mm > lend {
                    if e[mm - 1].abs() <= eps2 * (d[mm] * d[mm - 1]).abs() {
                        m = mm;
                        break;
                    }
                    mm -= 1;
                }
                if m > lend {
                    e[m - 1] = 0.0;
                }
                let mut p = d[l];
                if m == l {
                    // 140
                    d[l] = p;
                    l -= 1;
                    if l >= lend {
                        continue 'qr;
                    }
                    break 'qr;
                }
                if m + 1 == l {
                    let rte = e[l - 1].sqrt();
                    let (rt1, rt2) = dlae2(d[l], rte, d[l - 1]);
                    d[l] = rt1;
                    d[l - 1] = rt2;
                    e[l - 1] = 0.0;
                    if l < lend + 2 {
                        break 'qr;
                    }
                    l -= 2;
                    continue 'qr;
                }
                if jtot == nmaxit {
                    break 'qr;
                }
                jtot += 1;
                let rte = e[l - 1].sqrt();
                let mut sigma = (d[l - 1] - p) / (2.0 * rte);
                let r = dlapy2(sigma, 1.0);
                sigma = p - (rte / (sigma + r.copysign(sigma)));
                let mut c = 1.0;
                let mut s = 0.0;
                let mut gamma = d[m] - sigma;
                p = gamma * gamma;
                for i in m..l {
                    let bb = e[i];
                    let r = p + bb;
                    if i != m {
                        e[i - 1] = s * r;
                    }
                    let oldc = c;
                    c = p / r;
                    s = bb / r;
                    let oldgam = gamma;
                    let alpha = d[i + 1];
                    gamma = c * (alpha - sigma) - s * oldgam;
                    d[i] = oldgam + (alpha - gamma);
                    if c != 0.0 {
                        p = (gamma * gamma) / c;
                    } else {
                        p = oldc * bb;
                    }
                }
                e[l - 1] = s * p;
                d[l] = sigma + gamma;
            }
        }
        // 150
        if jtot >= nmaxit {
            return None;
        }
    }
    Some(d[1..].to_vec())
}
