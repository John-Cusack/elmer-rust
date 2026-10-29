//! The few gfortran edit descriptors the Elmer VTU writers rely on.
//!
//! The behaviours here were checked against output of the pristine
//! `ElmerSolver` (see `tests/fixtures/g0_times.tsv` and the ASCII `.vtu`
//! fixtures), not assumed.

/// `ES<width>.<digits>E3` of `val`, as `AscBinRealWrite` (GeneralUtils.F90)
/// uses it: `ES16.7E3` in double precision, `ES12.3E3` in single.
///
/// The caller handles the `|val| <= TINY` case, which `AscBinRealWrite` writes
/// as the literal `" 0.0"`.
pub(crate) fn es_e3(val: f64, digits: usize, width: usize) -> String {
    let body = if val.is_nan() {
        "NaN".to_string()
    } else if val.is_infinite() {
        if val < 0.0 { "-Infinity" } else { "Infinity" }.to_string()
    } else {
        let s = format!("{:.*e}", digits, val.abs());
        let (mant, exp) = s.split_once('e').expect("exponent form");
        let exp: i32 = exp.parse().expect("integer exponent");
        let sign = if val.is_sign_negative() { "-" } else { "" };
        format!(
            "{sign}{mant}E{}{:03}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    };
    format!("{body:>width$}")
}

/// One value as `AscBinRealWrite` writes it in ASCII mode.
pub(crate) fn ascii_real(val: f64, single_precision: bool) -> String {
    if val.abs() <= f64::MIN_POSITIVE {
        " 0.0".to_string()
    } else if single_precision {
        es_e3(val, 3, 12)
    } else {
        es_e3(val, 7, 16)
    }
}

/// gfortran `G0` for a `REAL(8)`, used for the `timestep` attribute of the
/// PVD file (`WritePvdFile`, VtuOutputSolver.F90:2227).
///
/// Seventeen significant digits. Values in `[0.1, 1e17)` print in fixed
/// notation (`0.25000000000000000`, `1234.5678000000000`, `10000000000000000.`);
/// everything else prints as `0.ddddddddddddddddd` followed by `E`, a sign and
/// the exponent without padding (`0.10000000000000000E-1`).
///
/// Zero has no sample in the Fortran runs (a transient time is never zero when
/// written), so it prints as seventeen fractional zeros.
pub(crate) fn g0(val: f64) -> String {
    if val.is_nan() {
        return "NaN".to_string();
    }
    if val.is_infinite() {
        return if val < 0.0 { "-Infinity" } else { "Infinity" }.to_string();
    }
    let sign = if val.is_sign_negative() { "-" } else { "" };
    if val == 0.0 {
        return format!("{sign}0.{}", "0".repeat(17));
    }
    // d.dddddddddddddddde<X>: correctly rounded to 17 significant digits.
    let s = format!("{:.16e}", val.abs());
    let (mant, exp) = s.split_once('e').expect("exponent form");
    let x: i32 = exp.parse().expect("integer exponent");
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    debug_assert_eq!(digits.len(), 17);
    if (-1..=16).contains(&x) {
        if x == -1 {
            format!("{sign}0.{digits}")
        } else {
            let ip = (x + 1) as usize;
            format!("{sign}{}.{}", &digits[..ip], &digits[ip..])
        }
    } else {
        let e = x + 1;
        format!("{sign}0.{digits}E{}{}", if e < 0 { '-' } else { '+' }, e.abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es_matches_fortran_samples() {
        assert_eq!(es_e3(0.5, 7, 16), "  5.0000000E-001");
        assert_eq!(es_e3(-9.0006453e-5, 7, 16), " -9.0006453E-005");
        assert_eq!(es_e3(1.0, 3, 12), "  1.000E+000");
        assert_eq!(ascii_real(0.0, false), " 0.0");
        assert_eq!(ascii_real(1.0e-310, false), " 0.0");
    }

    #[test]
    fn g0_basic() {
        assert_eq!(g0(0.25), "0.25000000000000000");
        assert_eq!(g0(1.5), "1.5000000000000000");
        assert_eq!(g0(0.01), "0.10000000000000000E-1");
        assert_eq!(g0(-2.0), "-2.0000000000000000");
    }
}
