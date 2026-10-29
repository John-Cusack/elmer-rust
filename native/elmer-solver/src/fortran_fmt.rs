//! The few gfortran edit descriptors Elmer's solver messages use, so the
//! messages are byte-identical (`ESw.d`, `Gw.d`, `I0`, list-directed integers).
//!
//! Digits come from Rust's `{:e}` formatting, which is correctly rounded from
//! the exact binary value, as libgfortran's default (`RN`-like) rounding is.

/// `ESw.d`: one nonzero digit before the point, `E+dd` (or `+ddd` past 99).
pub fn es(x: f64, w: usize, d: usize) -> String {
    let body = if !x.is_finite() {
        non_finite(x)
    } else {
        let (digits, e) = round_sig(x, d + 1);
        let mut s = String::new();
        if x.is_sign_negative() {
            s.push('-');
        }
        s.push_str(&digits[..1]);
        s.push('.');
        s.push_str(&digits[1..]);
        // round_sig's exponent is for 0.d1d2..., ES puts one digit before the point.
        s.push_str(&exponent(if x == 0.0 { 0 } else { e - 1 }));
        s
    };
    fit(body, w)
}

/// `Gw.d` without an exponent width: `F(w-4).(d-k)` plus four blanks when
/// `0.1 <= |x| < 10**d` after rounding to `d` digits, else `Ew.d`.
pub fn g(x: f64, w: usize, d: usize) -> String {
    if !x.is_finite() {
        return fit(non_finite(x), w);
    }
    let neg = x.is_sign_negative();
    let (digits, k) = if x == 0.0 {
        // Zero is edited as F(w-4).(d-1).
        ("0".repeat(d), 1)
    } else {
        round_sig(x, d)
    };
    if (0..=d as i32).contains(&k) {
        let k = k as usize;
        let mut s = String::new();
        if neg {
            s.push('-');
        }
        if k == 0 {
            s.push('0');
        } else {
            s.push_str(&digits[..k]);
        }
        s.push('.');
        s.push_str(&digits[k..]);
        let mut out = fit(s, w - 4);
        out.push_str("    ");
        out
    } else {
        let mut s = String::new();
        if neg {
            s.push('-');
        }
        s.push_str("0.");
        s.push_str(&digits);
        s.push_str(&exponent(k));
        if s.len() > w {
            // The optional leading zero is dropped before giving up.
            s = s.replacen("0.", ".", 1);
        }
        fit(s, w)
    }
}

/// Elmer `I2S`: an integer with no padding (`I0`).
pub fn i2s(i: i64) -> String {
    i.to_string()
}

/// A default-kind integer as gfortran's list-directed output writes it after
/// a character item: right-justified in 12 columns.
pub fn list_int(i: i32) -> String {
    format!("{i:>12}")
}

/// `|x|` rounded to `n` significant digits: (`n` digit characters, `k`) with
/// `|x| = 0.d1d2...dn * 10**k`.
fn round_sig(x: f64, n: usize) -> (String, i32) {
    if x == 0.0 {
        return ("0".repeat(n), 0);
    }
    let s = format!("{:.*e}", n - 1, x.abs());
    let (mant, exp) = s.split_once('e').expect("exponent in {:e}");
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let e: i32 = exp.parse().expect("integer exponent");
    (digits, e + 1)
}

fn exponent(e: i32) -> String {
    let sign = if e < 0 { '-' } else { '+' };
    if e.abs() <= 99 {
        format!("E{sign}{:02}", e.abs())
    } else {
        format!("{sign}{:03}", e.abs())
    }
}

fn non_finite(x: f64) -> String {
    if x.is_nan() {
        "NaN".into()
    } else if x > 0.0 {
        "Infinity".into()
    } else {
        "-Infinity".into()
    }
}

/// Right-justify in `w` columns, or `w` asterisks when it does not fit.
fn fit(s: String, w: usize) -> String {
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{s:>w$}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es_matches_gfortran() {
        assert_eq!(es(0.622811933, 15, 8), " 6.22811933E-01");
        assert_eq!(es(1.546336307357751, 15, 8), " 1.54633631E+00");
        assert_eq!(es(1.069559e-15, 13, 6), " 1.069559E-15");
        assert_eq!(es(0.0, 15, 8), " 0.00000000E+00");
        assert_eq!(es(-2.5e-120, 15, 8), "-2.50000000-120");
    }

    #[test]
    fn g_matches_gfortran() {
        assert_eq!(g(1.1666851, 15, 8), "  1.1666851    ");
        assert_eq!(g(2.0, 15, 8), "  2.0000000    ");
        assert_eq!(g(0.70704951, 15, 8), " 0.70704951    ");
        assert_eq!(g(0.089679829, 15, 8), " 0.89679829E-01");
        assert_eq!(g(0.0, 15, 8), "  0.0000000    ");
        // Rounding decides the range: 9.999999999 rounds to 10.000000.
        assert_eq!(g(9.999999999, 15, 8), "  10.000000    ");
        assert_eq!(g(0.09999999999, 15, 8), " 0.10000000    ");
        assert_eq!(g(123456789.0, 15, 8), " 0.12345679E+09");
    }

    #[test]
    fn list_int_width() {
        assert_eq!(format!(" Coupled system iteration: {}", list_int(1)),
            " Coupled system iteration:            1");
    }
}
