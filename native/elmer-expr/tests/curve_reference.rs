//! Bitwise comparison of `elmer_expr::curve` against a literal, 1-based
//! transliteration of the Fortran in `fem/src/GeneralUtils.F90`.
//! Arrays below carry a dummy element at index 0 so Fortran indices read as-is.

#![allow(clippy::assign_op_pattern, clippy::needless_range_loop)] // keep the Fortran statement shapes

use elmer_expr::curve;

const AEPS: f64 = 10.0 * f64::EPSILON;

fn pad(v: &[f64]) -> Vec<f64> {
    let mut o = vec![f64::NAN];
    o.extend_from_slice(v);
    o
}

fn solve_tri_diag(n: usize, y: &[f64], h: &[f64], r: &mut [f64]) {
    let mut b = vec![0.0; n + 1];
    for i in 2..=n - 1 {
        b[i] = 2.0 * (h[i - 1] + h[i]);
        r[i] = 3.0 * (h[i] * (y[i] - y[i - 1]) / h[i - 1] + h[i - 1] * (y[i + 1] - y[i]) / h[i]);
    }
    r[2] = r[2] - h[2] * r[1];
    for i in 2..=n.saturating_sub(2) {
        let s = -h[i + 1] / b[i];
        r[i + 1] = r[i + 1] + s * r[i];
        b[i + 1] = b[i + 1] + s * h[i - 1];
    }
    for i in (2..=n - 1).rev() {
        r[i] = (r[i] - h[i - 1] * r[i + 1]) / b[i];
    }
}

fn cubic_spline(n: usize, x: &[f64], y: &[f64], mono: bool) -> Vec<f64> {
    let mut h = vec![0.0; n + 1];
    let mut r = vec![0.0; n + 1];
    for i in 1..=n - 1 {
        h[i] = x[i + 1] - x[i];
    }
    r[1] = (y[2] - y[1]) / h[1];
    r[n] = (y[n] - y[n - 1]) / h[n - 1];
    if mono {
        for i in 1..=n - 1 {
            h[i] = (y[i + 1] - y[i]) / h[i];
        }
        for i in 2..=n - 1 {
            r[i] = (h[i - 1] + h[i]) / 2.0;
        }
        for i in 1..=n - 1 {
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
            let mut tau = (alpha.powi(2) + beta.powi(2)).sqrt();
            if tau > 3.0 {
                tau = 3.0 / tau;
                r[i] = alpha * tau * h[i];
                r[i + 1] = beta * tau * h[i];
            }
        }
    } else {
        solve_tri_diag(n, y, &h, &mut r);
    }
    r
}

fn cubic_spline_val(x: [f64; 2], y: [f64; 2], r: [f64; 2], t: f64) -> f64 {
    let h = x[1] - x[0];
    let a = -(2.0 * (y[1] - y[0])) + (r[0] + r[1]) * h;
    let b = 3.0 * (y[1] - y[0]) - (2.0 * r[0] + r[1]) * h;
    let c = r[0] * h;
    let d = y[0];
    let lt = (t - x[0]) / h;
    ((a * lt + b) * lt + c) * lt + d
}

fn search_interval(tval: &[f64], n: usize, t: f64) -> usize {
    let mut i;
    if t < tval[2] {
        i = 1;
    } else if t >= tval[n - 1] {
        i = n - 1;
    } else {
        let mut n0 = 1;
        let mut n1 = n;
        i = (n0 + n1) / 2;
        loop {
            if tval[i] <= t && tval[i + 1] > t {
                break;
            }
            if tval[i] > t {
                n1 = i - 1;
            } else {
                n0 = i + 1;
            }
            i = (n0 + n1) / 2;
        }
    }
    if i > n - 1 {
        i = n - 1;
    }
    i
}

fn interpolate_curve(tv: &[f64], fv: &[f64], n: usize, t: f64, cubic: Option<&[f64]>) -> f64 {
    if n == 1 {
        return fv[1] * t;
    }
    let i = search_interval(tv, n, t);
    let use_cubic = cubic.is_some() && t >= tv[1] && t <= tv[n];
    if use_cubic {
        let r = cubic.unwrap();
        cubic_spline_val([tv[i], tv[i + 1]], [fv[i], fv[i + 1]], [r[i], r[i + 1]], t)
    } else {
        let mut f = (t - tv[i]) / (tv[i + 1] - tv[i]);
        f = (1.0 - f) * fv[i] + f * fv[i + 1];
        f
    }
}

/// Deterministic generator, no dependencies.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn table(rng: &mut Lcg, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = Vec::with_capacity(n);
    let mut t = -3.0 + 10.0 * rng.next();
    for _ in 0..n {
        x.push(t);
        t += 0.05 + 2.0 * rng.next();
    }
    let y = (0..n).map(|_| 100.0 * rng.next() - 30.0).collect();
    (x, y)
}

#[test]
fn spline_coefficients_and_interpolation_are_bitwise_equal() {
    let mut rng = Lcg(12345);
    for &mono in &[false, true] {
        for n in 3..=40 {
            let (x, y) = table(&mut rng, n);
            let (px, py) = (pad(&x), pad(&y));
            let r_ref = cubic_spline(n, &px, &py, mono);
            let r = curve::cubic_spline(&x, &y, mono);
            for i in 0..n {
                assert_eq!(r[i].to_bits(), r_ref[i + 1].to_bits(), "n={n} mono={mono} r[{i}]");
            }
            let lo = x[0] - 1.5;
            let hi = x[n - 1] + 1.5;
            for k in 0..200 {
                let t = lo + (hi - lo) * (k as f64) / 199.0;
                let want_lin = interpolate_curve(&px, &py, n, t, None);
                let got_lin = curve::interpolate_curve(&x, &y, t, None);
                assert_eq!(got_lin.to_bits(), want_lin.to_bits(), "n={n} t={t} linear");
                let want_cub = interpolate_curve(&px, &py, n, t, Some(&r_ref));
                let got_cub = curve::interpolate_curve(&x, &y, t, Some(&r));
                assert_eq!(got_cub.to_bits(), want_cub.to_bits(), "n={n} t={t} cubic");
            }
        }
    }
}

#[test]
fn table_nodes_and_knots_reproduce_data() {
    let x = [0.0, 1.0, 2.5, 4.0, 5.0];
    let y = [1.0, 3.0, 2.0, 8.0, 7.0];
    let r = curve::cubic_spline(&x, &y, false);
    for i in 0..5 {
        let v = curve::interpolate_curve(&x, &y, x[i], Some(&r));
        assert!((v - y[i]).abs() < 1e-13, "knot {i}: {v}");
    }
}

#[test]
fn one_point_table_is_proportional() {
    assert_eq!(curve::interpolate_curve(&[1.0], &[4.0], 2.5, None), 10.0);
}

#[test]
fn integral_matches_trapezoid_and_cumulative() {
    let x = [0.0, 1.0, 3.0, 4.0];
    let y = [0.0, 2.0, 2.0, 6.0];
    let whole = curve::integrate_curve(&x, &y, None, None, None, None).unwrap();
    assert_eq!(whole, 1.0 + 4.0 + 4.0);
    let cum = curve::cumulative_integral(&x, &y, None);
    assert_eq!(cum, vec![0.0, 1.0, 5.0, 9.0]);
    let part = curve::integrate_curve(&x, &y, None, Some(0.5), Some(3.5), Some(&cum)).unwrap();
    let direct = curve::integrate_curve(&x, &y, None, Some(0.5), Some(3.5), None).unwrap();
    assert!((part - direct).abs() < 1e-14);
    // Linear extrapolation on both ends: region left of table.
    let left = curve::integrate_curve(&x, &y, None, Some(-1.0), Some(0.0), None).unwrap();
    assert!((left - (-1.0)).abs() < 1e-14, "left {left}");
    assert!(curve::integrate_curve(&x, &y, None, Some(2.0), Some(1.0), None).is_none());
}

#[test]
fn cubic_integral_is_antiderivative_of_cubic_interpolant() {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 0.5, 2.0, 3.0, 2.0];
    let r = curve::cubic_spline(&x, &y, false);
    let (a, b) = (0.3, 4.4);
    let exact = curve::integrate_curve(&x, &y, Some(&r), Some(a), Some(b), None).unwrap();
    // composite Simpson on the interpolant
    let n = 20000;
    let h = (b - a) / n as f64;
    let mut s = curve::interpolate_curve(&x, &y, a, Some(&r)) + curve::interpolate_curve(&x, &y, b, Some(&r));
    for k in 1..n {
        let w = if k % 2 == 1 { 4.0 } else { 2.0 };
        s += w * curve::interpolate_curve(&x, &y, a + k as f64 * h, Some(&r));
    }
    s *= h / 3.0;
    assert!((exact - s).abs() < 1e-9, "{exact} vs {s}");
}
