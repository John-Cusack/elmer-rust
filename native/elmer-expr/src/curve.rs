//! Table interpolation, cubic splines and curve integration.
//!
//! Direct transliterations of the curve routines in `fem/src/GeneralUtils.F90`
//! (`CubicSpline` :1576, `SolveTriDiag` :1524, `SearchInterval`, `InterpolateCurve`
//! :1756, `DerivateCurve`, `CumulativeIntegral`, `IntegrateCurve` :1928). The
//! floating-point operation order follows the Fortran source so results match
//! bit for bit. Fortran indices are 1-based; here every slice is 0-based and the
//! comments name the Fortran index where the shift matters.

/// `AEPS` from `Types.F90`: `10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// `CheckMonotone`: abscissae strictly increasing.
pub fn is_strictly_increasing(x: &[f64]) -> bool {
    x.windows(2).all(|w| w[1] > w[0])
}

/// `SolveTriDiag` (GeneralUtils.F90:1524). `r[0]` and `r[n-1]` hold the end
/// slopes on entry; `h[i] = x[i+1]-x[i]`.
fn solve_tridiag(n: usize, y: &[f64], h: &[f64], r: &mut [f64]) {
    let mut b = vec![0.0_f64; n];
    // DO i=2,n-1
    for i in 1..n - 1 {
        b[i] = 2.0 * (h[i - 1] + h[i]);
        r[i] = 3.0 * (h[i] * (y[i] - y[i - 1]) / h[i - 1] + h[i - 1] * (y[i + 1] - y[i]) / h[i]);
    }
    // r(2) = r(2) - h(2) * r(1)
    r[1] -= h[1] * r[0];
    // DO i=2,n-2
    for i in 1..n.saturating_sub(2) {
        let s = -h[i + 1] / b[i];
        r[i + 1] += s * r[i];
        // b(i+1) = b(i+1) + s*h(i-1)
        b[i + 1] += s * h[i - 1];
    }
    // DO i=n-1,2,-1
    for i in (1..n - 1).rev() {
        r[i] = (r[i] - h[i - 1] * r[i + 1]) / b[i];
    }
}

/// `CubicSpline`: slope coefficients `r` of the interpolating cubic through
/// `(x, y)`. With `monotone` the slopes are limited so the curve is monotone
/// between data points (the `Monotone Cubic` keyword). Needs `n >= 2`; Elmer
/// only builds splines for `n > 3`.
pub fn cubic_spline(x: &[f64], y: &[f64], monotone: bool) -> Vec<f64> {
    let n = x.len();
    assert!(n >= 2 && y.len() == n, "cubic_spline needs n >= 2 matching x, y");
    let mut h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let mut r = vec![0.0_f64; n];
    r[0] = (y[1] - y[0]) / h[0];
    r[n - 1] = (y[n - 1] - y[n - 2]) / h[n - 2];

    if monotone {
        for i in 0..n - 1 {
            h[i] = (y[i + 1] - y[i]) / h[i];
        }
        for i in 1..n - 1 {
            r[i] = (h[i - 1] + h[i]) / 2.0;
        }
        for i in 0..n - 1 {
            if h[i].abs() < 10.0 * AEPS {
                r[i] = 0.0;
                r[i + 1] = 0.0;
                continue;
            }
            let alpha = r[i] / h[i];
            let beta = r[i + 1] / h[i];
            if alpha < 0.0 || beta < 0.0 {
                r[i] = 0.0;
                continue;
            }
            let mut tau = (alpha * alpha + beta * beta).sqrt();
            if tau > 3.0 {
                tau = 3.0 / tau;
                r[i] = alpha * tau * h[i];
                r[i + 1] = beta * tau * h[i];
            }
        }
    } else if n > 2 {
        // Fortran reads h(2) unconditionally; with two points the end slopes stand.
        solve_tridiag(n, y, &h, &mut r);
    }
    r
}

/// Cubic Hermite coefficients of one interval, `CubicSplineVal` (:1634).
#[inline]
fn hermite(x: [f64; 2], y: [f64; 2], r: [f64; 2]) -> (f64, f64, f64, f64, f64) {
    let h = x[1] - x[0];
    let a = -2.0 * (y[1] - y[0]) + (r[0] + r[1]) * h;
    let b = 3.0 * (y[1] - y[0]) - (2.0 * r[0] + r[1]) * h;
    let c = r[0] * h;
    (a, b, c, y[0], h)
}

/// `CubicSplineVal`: evaluate the spline of one interval at `t`.
pub fn cubic_spline_val(x: [f64; 2], y: [f64; 2], r: [f64; 2], t: f64) -> f64 {
    let (a, b, c, d, h) = hermite(x, y, r);
    let lt = (t - x[0]) / h;
    ((a * lt + b) * lt + c) * lt + d
}

/// `CubicSplinedVal`: derivative of the spline of one interval at `t`.
pub fn cubic_spline_dval(x: [f64; 2], y: [f64; 2], r: [f64; 2], t: f64) -> f64 {
    let (a, b, c, _d, h) = hermite(x, y, r);
    let lt = (t - x[0]) / h;
    ((3.0 * a * lt + 2.0 * b) * lt + c) / h
}

/// `SearchInterval`: 0-based `i` with `tval[i] <= t < tval[i+1]`, clamped to the
/// first and last interval. Needs `tval.len() >= 2`.
pub fn search_interval(tval: &[f64], t: f64) -> usize {
    let n = tval.len() as i64;
    debug_assert!(n >= 2);
    // Fortran 1-based `i`; tval(k) is tval[k-1].
    let tv = |k: i64| tval[(k - 1) as usize];
    let mut i: i64;
    if t < tv(2) {
        i = 1;
    } else if t >= tv(n - 1) {
        i = n - 1;
    } else {
        let mut n0 = 1_i64;
        let mut n1 = n;
        i = (n0 + n1) / 2;
        loop {
            if tv(i) <= t && tv(i + 1) > t {
                break;
            }
            if tv(i) > t {
                n1 = i - 1;
            } else {
                n0 = i + 1;
            }
            // NaN comparisons are all false and would walk off the table.
            if n0 > n1 {
                break;
            }
            i = (n0 + n1) / 2;
        }
    }
    if i > n - 1 {
        i = n - 1;
    }
    (i - 1) as usize
}

/// `InterpolateCurve`: linear table lookup, or cubic when `cubic` holds the
/// spline coefficients from [`cubic_spline`]. A one-point table means `y = f0 * t`.
pub fn interpolate_curve(tv: &[f64], fv: &[f64], t: f64, cubic: Option<&[f64]>) -> f64 {
    let n = tv.len();
    if n == 1 {
        return fv[0] * t;
    }
    let i = search_interval(tv, t);
    let use_cubic = cubic.is_some() && t >= tv[0] && t <= tv[n - 1];
    if use_cubic {
        let r = cubic.unwrap();
        cubic_spline_val(
            [tv[i], tv[i + 1]],
            [fv[i], fv[i + 1]],
            [r[i], r[i + 1]],
            t,
        )
    } else {
        let f = (t - tv[i]) / (tv[i + 1] - tv[i]);
        (1.0 - f) * fv[i] + f * fv[i + 1]
    }
}

/// `DerivateCurve`: slope of the table or spline at `t`.
pub fn derivate_curve(tv: &[f64], fv: &[f64], t: f64, cubic: Option<&[f64]>) -> f64 {
    let n = tv.len();
    let i = search_interval(tv, t);
    let use_cubic = cubic.is_some() && t >= tv[0] && t <= tv[n - 1];
    if use_cubic {
        let r = cubic.unwrap();
        cubic_spline_dval(
            [tv[i], tv[i + 1]],
            [fv[i], fv[i + 1]],
            [r[i], r[i + 1]],
            t,
        )
    } else {
        (fv[i + 1] - fv[i]) / (tv[i + 1] - tv[i])
    }
}

/// `CumulativeIntegral`: integral from the first abscissa to each table point
/// (whole intervals only).
pub fn cumulative_integral(tv: &[f64], fv: &[f64], cubic: Option<&[f64]>) -> Vec<f64> {
    let n = tv.len();
    let mut cum = vec![0.0_f64; n];
    for i in 0..n.saturating_sub(1) {
        let h = tv[i + 1] - tv[i];
        let (y1, y2) = (fv[i], fv[i + 1]);
        cum[i + 1] = match cubic {
            Some(r) => {
                let (r1, r2) = (r[i], r[i + 1]);
                let a = (-2.0 * (y2 - y1) + (r1 + r2) * h) / 4.0;
                let b = (3.0 * (y2 - y1) - (2.0 * r1 + r2) * h) / 3.0;
                let c = (r1 * h) / 2.0;
                cum[i] + h * (a + b + c + y1)
            }
            None => {
                let c = (y2 - y1) / 2.0;
                cum[i] + h * (c + y1)
            }
        };
    }
    cum
}

/// Antiderivative of one interval's polynomial between local coordinates
/// `s0..s1` (`s` in units of the interval width).
#[inline]
fn interval_integral(
    h: f64,
    y: [f64; 2],
    r: Option<[f64; 2]>,
    s0: f64,
    s1: f64,
) -> f64 {
    match r {
        Some(r) => {
            let a = (-2.0 * (y[1] - y[0]) + (r[0] + r[1]) * h) / 4.0;
            let b = (3.0 * (y[1] - y[0]) - (2.0 * r[0] + r[1]) * h) / 3.0;
            let c = (r[0] * h) / 2.0;
            let d = y[0];
            h * ((((a * s1 + b) * s1 + c) * s1 + d) * s1 - (((a * s0 + b) * s0 + c) * s0 + d) * s0)
        }
        None => {
            let c = (y[1] - y[0]) / 2.0;
            let d = y[0];
            h * ((c * s1 + d) * s1 - (c * s0 + d) * s0)
        }
    }
}

/// `IntegrateCurve`: integral of the table or spline over `[t0, t1]` (defaulting
/// to the whole table). Outside the table the end intervals are extrapolated
/// linearly. Returns `None` where Fortran reports `Found = .FALSE.`: fewer than
/// two points, mismatched lengths, or `t0 >= t1`. `cumulative` is the optional
/// precomputed [`cumulative_integral`].
pub fn integrate_curve(
    tv: &[f64],
    fv: &[f64],
    cubic: Option<&[f64]>,
    t0: Option<f64>,
    t1: Option<f64>,
    cumulative: Option<&[f64]>,
) -> Option<f64> {
    let n = tv.len();
    if n < 2 || fv.len() != n {
        return None;
    }
    let mut tt0 = t0.unwrap_or(tv[0]);
    let mut tt1 = t1.unwrap_or(tv[n - 1]);
    if tt0 >= tt1 {
        return None;
    }
    let mut sumf = 0.0_f64;

    // Whole range left of the table (t1 <= first).
    if tt1 <= tv[0] {
        let h = tv[1] - tv[0];
        let s0 = (tt0 - tv[0]) / h;
        let s1 = (tt1 - tv[0]) / h;
        let c = (fv[1] - fv[0]) / 2.0;
        let d = fv[0];
        return Some(sumf + h * ((c * s1 + d) * s1 - (c * s0 + d) * s0));
    }
    // Whole range right of the table (t0 >= last).
    if tt0 >= tv[n - 1] {
        let h = tv[n - 1] - tv[n - 2];
        let s0 = (tt0 - tv[n - 2]) / h;
        let s1 = (tt1 - tv[n - 2]) / h;
        let c = (fv[n - 1] - fv[n - 2]) / 2.0;
        let d = fv[n - 2];
        return Some(sumf + h * ((c * s1 + d) * s1 - (c * s0 + d) * s0));
    }
    // Part left of the table.
    if tt0 < tv[0] {
        let h = tv[1] - tv[0];
        let s0 = (tt0 - tv[0]) / h;
        let c = (fv[1] - fv[0]) / 2.0;
        let d = fv[0];
        sumf -= h * (c * s0 + d) * s0;
        tt0 = tv[0];
    }
    // Part right of the table.
    if tt1 > tv[n - 1] {
        let h = tv[n - 1] - tv[n - 2];
        let s1 = (tt1 - tv[n - 2]) / h;
        let c = (fv[n - 1] - fv[n - 2]) / 2.0;
        let d = fv[n - 2];
        sumf += h * ((c * s1 + d) * s1 - (c + d));
        tt1 = tv[n - 1];
    }
    if tt0 >= tt1 {
        return Some(sumf);
    }

    let coeffs = |i: usize| cubic.map(|r| [r[i], r[i + 1]]);

    let mut i0 = search_interval(tv, tt0);
    // First (possibly partial) interval.
    let h = tv[i0 + 1] - tv[i0];
    let s0 = (tt0 - tv[i0]) / h;
    let s1 = ((tt1 - tv[i0]) / h).min(1.0);
    if s0 > 0.0 || s1 < 1.0 {
        sumf += interval_integral(h, [fv[i0], fv[i0 + 1]], coeffs(i0), s0, s1);
        i0 += 1;
        tt0 = tv[i0];
        if tt0 >= tt1 {
            return Some(sumf);
        }
    }

    let mut i1 = search_interval(tv, tt1);
    // Last (possibly partial) interval.
    let h = tv[i1 + 1] - tv[i1];
    let s0 = ((tt0 - tv[i1]) / h).max(0.0);
    let s1 = (tt1 - tv[i1]) / h;
    if s0 > 0.0 || s1 < 1.0 {
        sumf += interval_integral(h, [fv[i1], fv[i1 + 1]], coeffs(i1), s0, s1);
        // Fortran: i1 = i1 - 1; tt1 = Tvalues(i1+1)  (1-based)
        if i1 == 0 {
            return Some(sumf);
        }
        i1 -= 1;
        tt1 = tv[i1 + 1];
        if tt0 >= tt1 {
            return Some(sumf);
        }
    }

    // Whole intervals i0..=i1 (0-based, inclusive).
    if let Some(cum) = cumulative {
        return Some(sumf + cum[i1 + 1] - cum[i0]);
    }
    for i in i0..=i1 {
        let h = tv[i + 1] - tv[i];
        let (y1, y2) = (fv[i], fv[i + 1]);
        match cubic {
            Some(r) => {
                let (r1, r2) = (r[i], r[i + 1]);
                let a = (-2.0 * (y2 - y1) + (r1 + r2) * h) / 4.0;
                let b = (3.0 * (y2 - y1) - (2.0 * r1 + r2) * h) / 3.0;
                let c = (r1 * h) / 2.0;
                sumf += h * (a + b + c + y1);
            }
            None => {
                let c = (y2 - y1) / 2.0;
                sumf += h * (c + y1);
            }
        }
    }
    Some(sumf)
}
