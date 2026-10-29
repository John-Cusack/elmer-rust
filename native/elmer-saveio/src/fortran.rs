//! Fortran formatted-I/O primitives used by every Elmer text file.
//!
//! Only the edit descriptors Elmer actually writes are implemented: `ESw.d`,
//! `ESw.dEe`, `Ew.d`, `Iw`, `I0` and the list-directed (`WRITE(u,*)`) forms of
//! `REAL(8)` and `INTEGER(4)` as gfortran prints them. A value that does not
//! fit its field becomes a row of `*`, as in Fortran.
//!
//! The list-directed layouts were read off files written by the pristine
//! ElmerSolver (`.result` bodies) and are checked against them in
//! `tests/format_golden.rs`.

/// Right-justify `s` in `w` columns; `*` fill when it does not fit.
/// `w == 0` means "minimum width" (no padding).
fn fit(s: String, w: usize) -> String {
    if w == 0 || s.len() == w {
        s
    } else if s.len() < w {
        format!("{s:>w$}")
    } else {
        "*".repeat(w)
    }
}

/// `NaN` / `Infinity` as gfortran writes them for a field of width `w`.
fn non_finite(x: f64, w: usize) -> String {
    let s = if x.is_nan() {
        "NaN".to_string()
    } else if w == 0 || w >= 8 + usize::from(x < 0.0) {
        if x < 0.0 { "-Infinity" } else { "Infinity" }.to_string()
    } else if x < 0.0 {
        "-Inf".to_string()
    } else {
        "Inf".to_string()
    };
    fit(s, w)
}

/// Correctly rounded `ndig` significant digits of `|x|` (`x` finite) and the
/// decimal exponent of the first digit (`x = d.ddd * 10^exp`). Zero gives
/// all-zero digits and exponent 0.
fn sci_digits(x: f64, ndig: usize) -> (String, i32) {
    debug_assert!(x.is_finite() && ndig >= 1);
    let s = format!("{:.*e}", ndig - 1, x.abs());
    let (mant, exp) = s.split_once('e').expect("exponent marker");
    (mant.replace('.', ""), exp.parse().expect("exponent"))
}

/// Exponent field. `e = Some(n)`: `E` and exactly `n` digits. `None`: `E` and
/// two digits, or sign and three digits (no `E`) once `|exp| > 99`.
fn exp_field(exp: i32, e: Option<usize>) -> Option<String> {
    let sign = if exp < 0 { '-' } else { '+' };
    let a = exp.unsigned_abs();
    match e {
        Some(n) => {
            let digits = a.to_string();
            (digits.len() <= n).then(|| format!("E{sign}{digits:0>n$}"))
        }
        None if a <= 99 => Some(format!("E{sign}{a:02}")),
        None if a <= 999 => Some(format!("{sign}{a:03}")),
        None => None,
    }
}

fn sign_prefix(x: f64) -> &'static str {
    if x.is_sign_negative() {
        "-"
    } else {
        ""
    }
}

/// `ESw.dEe` (scientific, one digit before the point). `e = None` is `ESw.d`.
pub fn fortran_es_e(x: f64, w: usize, d: usize, e: Option<usize>) -> String {
    if !x.is_finite() {
        return non_finite(x, w);
    }
    let (digits, exp) = sci_digits(x, d + 1);
    let Some(ef) = exp_field(exp, e) else {
        return fit("*".repeat(w.max(1) + 1), w.max(1));
    };
    let body = format!("{}{}.{}{}", sign_prefix(x), &digits[..1], &digits[1..], ef);
    fit(body, w)
}

/// `ESw.d`: scientific notation, e.g. `fortran_es(0.5, 12, 3) == "   5.000E-01"`.
pub fn fortran_es(x: f64, w: usize, d: usize) -> String {
    fortran_es_e(x, w, d, None)
}

/// `Ew.dEe` (normalized `0.ddd` mantissa). `e = None` is `Ew.d`.
pub fn fortran_e_e(x: f64, w: usize, d: usize, e: Option<usize>) -> String {
    if !x.is_finite() {
        return non_finite(x, w);
    }
    let (digits, exp) = sci_digits(x, d.max(1));
    let (digits, exp10) = if x == 0.0 {
        ("0".repeat(d), 0)
    } else {
        (digits, exp + 1)
    };
    let Some(ef) = exp_field(exp10, e) else {
        return fit("*".repeat(w.max(1) + 1), w.max(1));
    };
    // `E` keeps a leading zero only when the field is wide enough; Elmer never
    // narrows it, so the zero is always written.
    let body = format!("{}0.{}{}", sign_prefix(x), digits, ef);
    fit(body, w)
}

/// `Ew.d`: e.g. `fortran_e(0.5, 12, 3) == "   0.500E+00"`.
pub fn fortran_e(x: f64, w: usize, d: usize) -> String {
    fortran_e_e(x, w, d, None)
}

/// `Iw` (`w == 0` is `I0`).
pub fn fortran_i(n: i64, w: usize) -> String {
    fit(n.to_string(), w)
}

/// `WRITE(u,*) x` for `REAL(8)`: the item with its leading blank, exactly as
/// gfortran lays it out (17 significant digits).
///
/// * `0.1 <= |x| < 1e17` or `x == 0`: fixed notation, number right-justified
///   in 21 columns, followed by 5 blanks;
/// * otherwise `d.dddddddddddddddd` + `E±ddd`, right-justified in 26 columns.
pub fn list_real(x: f64) -> String {
    if !x.is_finite() {
        return fit(non_finite(x, 0), 26);
    }
    let (digits, exp) = sci_digits(x, 17);
    let sign = sign_prefix(x);
    if x == 0.0 {
        return format!("{:>21}     ", format!("{sign}0.0000000000000000"));
    }
    if (-1..=16).contains(&exp) {
        // 17 significant digits laid out in fixed notation.
        let body = if exp < 0 {
            format!("0.{digits}")
        } else {
            let int_len = exp as usize + 1;
            if int_len >= 17 {
                format!("{}{}.", digits, "0".repeat(int_len - 17))
            } else {
                format!("{}.{}", &digits[..int_len], &digits[int_len..])
            }
        };
        return format!("{:>21}     ", format!("{sign}{body}"));
    }
    let ef = format!(
        "E{}{:03}",
        if exp < 0 { '-' } else { '+' },
        exp.unsigned_abs()
    );
    format!(
        "{:>26}",
        format!("{sign}{}.{}{ef}", &digits[..1], &digits[1..])
    )
}

/// `WRITE(u,*) n` for `INTEGER(4)`: a blank separator and `I11`.
pub fn list_i32(n: i32) -> String {
    format!("{n:>12}")
}

/// Fortran list-directed read of the first real on `s` (blank or comma
/// separated; `D` exponents; `NaN`, `Inf`, `Infinity`).
pub fn parse_list_real(s: &str) -> Option<f64> {
    let tok = s
        .split(|c: char| c.is_whitespace() || c == ',')
        .find(|t| !t.is_empty())?;
    parse_real_token(tok)
}

/// One real token (no separators), see [`parse_list_real`].
pub fn parse_real_token(tok: &str) -> Option<f64> {
    let lower = tok.to_ascii_lowercase();
    match lower.trim_start_matches(['+', '-']) {
        "nan" => return Some(f64::NAN),
        "inf" | "infinity" => {
            return Some(if lower.starts_with('-') {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            })
        }
        _ => {}
    }
    let t = lower.replace('d', "e");
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    // `1.0-100`: three-digit exponent printed without `E`.
    let split = t
        .rfind(['+', '-'])
        .filter(|&i| i > 0 && !t[..i].ends_with('e'))?;
    format!("{}e{}", &t[..split], &t[split..]).parse().ok()
}

/// Elmer's `FormatDate()`: `yyyy/mm/dd hh:mm:ss`, here from UTC seconds since
/// the Unix epoch (the Fortran code uses local time; callers that need a
/// specific stamp pass their own string).
pub fn format_elmer_date(unix_secs: u64) -> String {
    let days = (unix_secs / 86_400) as i64;
    let rem = unix_secs % 86_400;
    // Civil-from-days (Howard Hinnant).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{y:04}/{m:02}/{d:02} {:02}:{:02}:{:02}",
        rem / 3_600,
        rem % 3_600 / 60,
        rem % 60
    )
}

/// `format_elmer_date` for the current wall clock.
pub fn now_elmer_date() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    format_elmer_date(secs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es_basic() {
        assert_eq!(fortran_es(0.5, 12, 3), "   5.000E-01");
        assert_eq!(fortran_es(-12345.678, 12, 3), "  -1.235E+04");
        assert_eq!(fortran_es(0.0, 12, 3), "   0.000E+00");
        assert_eq!(fortran_es(1e100, 12, 3), "   1.000+100");
        assert_eq!(fortran_es(1.0, 3, 3), "***");
    }

    #[test]
    fn es_e3() {
        assert_eq!(fortran_es_e(1.0, 17, 8, Some(3)), "  1.00000000E+000");
        assert_eq!(fortran_es_e(0.5, 17, 8, Some(3)), "  5.00000000E-001");
        assert_eq!(
            fortran_es_e(-2.5e-7, 22, 12, Some(3)),
            "  -2.500000000000E-007"
        );
        assert_eq!(fortran_es_e(f64::MAX, 17, 8, Some(3)), "  1.79769313E+308");
    }

    #[test]
    fn e_basic() {
        assert_eq!(fortran_e(0.5, 12, 3), "   0.500E+00");
        assert_eq!(fortran_e(123.456, 12, 4), "  0.1235E+03");
        assert_eq!(fortran_e(0.0, 12, 3), "   0.000E+00");
    }

    #[test]
    fn list_real_layout() {
        assert_eq!(list_real(0.0), "   0.0000000000000000     ");
        assert_eq!(list_real(0.5), "  0.50000000000000000     ");
        assert_eq!(list_real(1.0).len(), 26);
        assert_eq!(list_real(1.0e-6), "   9.9999999999999995E-007".to_string());
        assert_eq!(list_real(-1.0e-6), "  -9.9999999999999995E-007".to_string());
        assert_eq!(list_real(1.0), "   1.0000000000000000     ");
    }

    #[test]
    fn parse_variants() {
        assert_eq!(parse_list_real("  0.5000E+00  "), Some(0.5));
        assert_eq!(parse_list_real("1.5D+02"), Some(150.0));
        assert_eq!(parse_list_real("1.0-100"), Some(1.0e-100));
        assert!(parse_list_real("NaN").unwrap().is_nan());
        assert_eq!(parse_list_real("-Infinity"), Some(f64::NEG_INFINITY));
        assert_eq!(parse_list_real("   "), None);
    }

    #[test]
    fn date() {
        assert_eq!(format_elmer_date(0), "1970/01/01 00:00:00");
        assert_eq!(format_elmer_date(1_793_000_000), "2026/10/26 07:33:20");
    }
}
