//! HUTI per-iteration stdout lines and the Fortran edit descriptors behind them.

/// The standard HUTI iteration line, `WRITE (*, '(I8, E11.4)') iter, residual`
/// (huti_cg.F90:472 and the CG/CGS/BiCGStab/BiCGStab(2)/QMR solvers).
pub fn iteration_line(it: usize, res: f64) -> String {
    format!("{}{}", format_i(it, 8), format_e(res, 11, 4))
}

/// GMRES prefixes the line: `'(A, I8, E11.4)') '   gmres:'` (`'   gmresz:'`
/// for complex), huti_gmres.F90:779 and :1399.
pub(crate) fn gmres_line(complex: bool, it: usize, res: f64) -> String {
    let tag = if complex { "   gmresz:" } else { "   gmres:" };
    format!("{tag}{}", iteration_line(it, res))
}

/// TFQMR prints `'(I8, E17.7)'`, or `'(I8, 2E17.7)'` with the upper bound
/// under the upper-bound stopping rule (huti_tfqmr.F90:739-744).
pub(crate) fn tfqmr_line(it: usize, res: f64, upperb: Option<f64>) -> String {
    let mut s = format!("{}{}", format_i(it, 8), format_e(res, 17, 7));
    if let Some(u) = upperb {
        s.push_str(&format_e(u, 17, 7));
    }
    s
}

/// Fortran `Iw`: right-justified, all asterisks when it does not fit.
pub fn format_i(v: usize, w: usize) -> String {
    let s = v.to_string();
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{s:>w$}")
    }
}

/// Fortran `Ew.d` as gfortran writes it: `0.ddddE+xx`, a three-digit exponent
/// drops the `E` (`0.dddd+xxx`), NaN and Infinity are right-justified words,
/// and the optional leading zero goes first when the field is too narrow.
pub fn format_e(x: f64, w: usize, d: usize) -> String {
    let body = if x.is_nan() {
        "NaN".to_string()
    } else if x.is_infinite() {
        let word = if w >= 8 + usize::from(x < 0.0) {
            "Infinity"
        } else {
            "Inf"
        };
        format!("{}{word}", if x < 0.0 { "-" } else { "" })
    } else {
        let (digits, exp) = if x == 0.0 {
            ("0".repeat(d), 0)
        } else {
            // Correctly rounded d significant digits, as gfortran gets them
            // from printf under the default (processor-defined) rounding.
            let s = format!("{:.*e}", d.saturating_sub(1), x.abs());
            let (mant, e) = s.split_once('e').expect("exponent form");
            (mant.replace('.', ""), e.parse::<i32>().expect("exponent") + 1)
        };
        let sign = if x.is_sign_negative() { "-" } else { "" };
        let e = if exp.abs() <= 99 {
            format!("E{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
        } else {
            format!("{}{:03}", if exp < 0 { '-' } else { '+' }, exp.abs())
        };
        let full = format!("{sign}0.{digits}{e}");
        if full.len() > w {
            format!("{sign}.{digits}{e}")
        } else {
            full
        }
    };
    if body.len() > w {
        "*".repeat(w)
    } else {
        format!("{body:>w$}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_match_captured_elmer_output() {
        // Lines copied from pristine ElmerSolver runs (HeatAniso, CG/GMRES/TFQMR).
        assert_eq!(iteration_line(191, 9.5304e-9), "     191 0.9530E-08");
        assert_eq!(iteration_line(1, 1.50797), "       1 0.1508E+01");
        assert_eq!(iteration_line(297, f64::INFINITY), "     297   Infinity");
        assert_eq!(gmres_line(false, 768, 1.0e-8), "   gmres:     768 0.1000E-07");
        assert_eq!(
            tfqmr_line(2001, 0.73548563, None),
            "    2001    0.7354856E+00"
        );
    }

    #[test]
    fn edge_values() {
        assert_eq!(format_e(0.0, 11, 4), " 0.0000E+00");
        assert_eq!(format_e(-0.0, 11, 4), "-0.0000E+00");
        assert_eq!(format_e(f64::NAN, 11, 4), "        NaN");
        assert_eq!(format_e(1.0e100, 11, 4), " 0.1000+101");
        assert_eq!(format_e(-1.5e-120, 11, 4), "-0.1500-119");
        assert_eq!(format_e(0.99996, 11, 4), " 0.1000E+01");
        assert_eq!(format_i(123_456_789, 8), "********");
    }
}
