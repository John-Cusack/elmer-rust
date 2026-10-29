//! gfortran list-directed-free edit descriptors that Elmer's messages use.
//!
//! Only the descriptors that appear in the stdout the GUI parses are
//! implemented: `Gw.d` (SolverBasics.F90:4146, `g15.8`) and `ESw.d`
//! (ElmerSolver.F90:3127, `ES12.3`).

/// Right-justifies `s` in `w` columns, or fills the field with `*` when it
/// does not fit (Fortran output overflow).
fn justify(s: &str, w: usize) -> String {
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{s:>w$}")
    }
}

/// Splits `{:.*e}` output into (digits without the point, decimal exponent).
fn sci_digits(x: f64, d: usize) -> (String, i32) {
    let s = format!("{:.*e}", d.saturating_sub(1), x);
    let (mant, exp) = s.split_once('e').expect("exponent form");
    (mant.replace('.', ""), exp.parse().expect("exponent"))
}

/// Exponent field of an `E`/`ES` edit: `E+02`, or `+100` once three digits are
/// needed (the `E` is dropped, as gfortran does).
fn exp_field(e: i32) -> String {
    if e.abs() <= 99 {
        format!("E{e:+03}")
    } else {
        format!("{e:+04}")
    }
}

fn nonfinite(x: f64, w: usize) -> Option<String> {
    if x.is_nan() {
        Some(justify("NaN", w))
    } else if x.is_infinite() {
        let neg = x < 0.0;
        let text = match (w >= 8, neg) {
            (true, false) => "Infinity",
            (true, true) if w >= 9 => "-Infinity",
            (false, false) => "Inf",
            _ => "-Inf",
        };
        Some(justify(text, w))
    } else {
        None
    }
}

/// `Gw.d`: `F(w-4).(d-n)` followed by four blanks when `0.1 <= |x| < 10**d`
/// (after rounding to `d` significant digits), else `Ew.d`. Zero is written as
/// `F(w-4).(d-1)`.
pub fn fmt_g(x: f64, w: usize, d: usize) -> String {
    if let Some(s) = nonfinite(x, w) {
        return s;
    }
    let neg = x.is_sign_negative();
    let sign = if neg { "-" } else { "" };
    if x == 0.0 {
        let body = format!("{sign}{:.*}", d - 1, 0.0);
        return format!("{}    ", justify(&body, w - 4));
    }
    let (digits, e) = sci_digits(x.abs(), d);
    let n = e + 1; // digits before the decimal point
    if n >= 0 && n as usize <= d {
        let decimals = d - n as usize;
        let mut body = format!("{sign}{:.*}", decimals, x.abs());
        if decimals == 0 {
            body.push('.');
        }
        format!("{}    ", justify(&body, w - 4))
    } else {
        justify(&format!("{sign}0.{digits}{}", exp_field(n)), w)
    }
}

/// `ESw.d`: scientific form with one digit before the point.
pub fn fmt_es(x: f64, w: usize, d: usize) -> String {
    if let Some(s) = nonfinite(x, w) {
        return s;
    }
    let sign = if x.is_sign_negative() { "-" } else { "" };
    let (digits, e) = if x == 0.0 {
        ("0".repeat(d + 1), 0)
    } else {
        sci_digits(x.abs(), d + 1)
    };
    let (int, frac) = digits.split_at(1);
    justify(&format!("{sign}{int}.{frac}{}", exp_field(e)), w)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g15_8_matches_gfortran_captures() {
        assert_eq!(fmt_g(0.75696824e-2, 15, 8), " 0.75696824E-02");
        assert_eq!(fmt_g(2.0, 15, 8), "  2.0000000    ");
        assert_eq!(fmt_g(0.93065875, 15, 8), " 0.93065875    ");
        assert_eq!(fmt_g(0.0, 15, 8), "  0.0000000    ");
        assert_eq!(fmt_g(0.62474118, 15, 8), " 0.62474118    ");
    }

    #[test]
    fn g_range_edges() {
        // 0.1 is the first F-branch value, just under it switches to E.
        assert_eq!(fmt_g(0.1, 15, 8), " 0.10000000    ");
        assert_eq!(fmt_g(0.09999, 15, 8), " 0.99990000E-01");
        // Rounds up across 10**d: E branch.
        assert_eq!(fmt_g(99999999.7, 15, 8), " 0.10000000E+09");
        assert_eq!(fmt_g(12345678.0, 15, 8), "  12345678.    ");
        assert_eq!(fmt_g(-3.5, 15, 8), " -3.5000000    ");
        assert_eq!(fmt_g(1.0e100, 15, 8), " 0.10000000+101");
        assert_eq!(fmt_g(f64::NAN, 15, 8), "            NaN");
        assert_eq!(fmt_g(f64::INFINITY, 15, 8), "       Infinity");
    }

    #[test]
    fn es12_3() {
        assert_eq!(fmt_es(0.1, 12, 3), "   1.000E-01");
        assert_eq!(fmt_es(12.0, 12, 3), "   1.200E+01");
        assert_eq!(fmt_es(0.0, 12, 3), "   0.000E+00");
        assert_eq!(fmt_es(-2.5e-7, 12, 3), "  -2.500E-07");
        assert_eq!(fmt_es(9.9996, 12, 3), "   1.000E+01");
    }
}
