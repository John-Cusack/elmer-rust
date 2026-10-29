//! Arithmetic that must match gfortran's, beyond what `elmer_core::Scalar` offers.

use elmer_core::{CsrMatrix, Scalar};
use num_complex::Complex64;

/// `AEPS = 10 * EPSILON(1.0_dp)` (Types.F90:71): pivots at or below it are not divided by.
pub const AEPS: f64 = 10.0 * f64::EPSILON;

/// `TINY(1.0_dp)`.
pub const TINY: f64 = f64::MIN_POSITIVE;

/// Row-`i` diagonal position in `a.values`, `None` when the row stores no diagonal.
pub fn diag_index<T: Scalar>(a: &CsrMatrix<T>, i: usize) -> Option<usize> {
    if let Some(d) = &a.diag {
        let k = d[i];
        return (k != crate::ops::NO_DIAG).then_some(k);
    }
    let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
    a.col_idx[lo..hi]
        .binary_search(&(i as u32))
        .ok()
        .map(|p| lo + p)
}

/// A [`Scalar`] whose operations reproduce the Fortran build bit for bit.
///
/// The Fortran complex code works on an interleaved real matrix; the native
/// complex CSR is the same arithmetic on `Complex64`. Complex division is
/// libgcc's `__divdc3`, which differs from `num_complex`'s `Div` in the last bit.
pub trait ElmerScalar: Scalar {
    /// Real and imaginary parts are separate reals in Elmer's storage (1 or 2).
    const COMPONENTS: usize;
    /// ILUT keeps `|s| >= tol*norm` for real (`CRS_ILUT`) but `|s| > tol*norm` for complex.
    const DROP_INCLUSIVE: bool;
    /// `CRS_DiagPrecondition` skips tiny pivots; the complex version does not.
    const JACOBI_GUARD: bool;

    /// gfortran's `self / d`.
    fn elmer_div(self, d: Self) -> Self;
    /// Visit the real components in Elmer's storage order (`re`, `im`).
    fn for_each_component(self, f: &mut dyn FnMut(f64));
    /// Multiply every real component by `s`.
    fn scale_real(self, s: f64) -> Self;
    /// Divide every real component by `s`.
    fn unscale_real(self, s: f64) -> Self;
    /// `Diag(i)` of `ScaleLinearSystemDiagonal`: one positive factor per row.
    fn diagonal_scaling(a: &CsrMatrix<Self>) -> Vec<f64>;
    /// The terms `CRS`-style diagonal sums add for one diagonal entry.
    fn diag_abs_terms(self, f: &mut dyn FnMut(f64));
}

impl ElmerScalar for f64 {
    const COMPONENTS: usize = 1;
    const DROP_INCLUSIVE: bool = true;
    const JACOBI_GUARD: bool = true;

    fn elmer_div(self, d: f64) -> f64 {
        self / d
    }
    fn for_each_component(self, f: &mut dyn FnMut(f64)) {
        f(self)
    }
    fn scale_real(self, s: f64) -> f64 {
        self * s
    }
    fn unscale_real(self, s: f64) -> f64 {
        self / s
    }
    fn diag_abs_terms(self, f: &mut dyn FnMut(f64)) {
        f(self.abs())
    }

    fn diagonal_scaling(a: &CsrMatrix<f64>) -> Vec<f64> {
        let n = a.nrows;
        let mut d = vec![0.0; n];
        for (i, di) in d.iter_mut().enumerate() {
            if let Some(j) = diag_index(a, i) {
                *di = a.values[j];
            }
        }
        if d.iter().any(|v| v.abs() <= TINY) {
            for (i, di) in d.iter_mut().enumerate() {
                if di.abs() <= TINY {
                    let mut s = 0.0;
                    for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                        s += a.values[k].abs();
                    }
                    *di = s;
                }
            }
        }
        for v in d.iter_mut() {
            *v = if v.abs() > TINY {
                1.0 / v.abs().sqrt()
            } else {
                1.0
            };
        }
        d
    }
}

/// libgcc `__divdc3` (the pre-GCC-12 formulas plus the subnormal-ratio branch).
fn divdc3(a: f64, b: f64, c: f64, d: f64) -> (f64, f64) {
    let (x, y);
    if c.abs() < d.abs() {
        let ratio = c / d;
        let denom = (c * ratio) + d;
        if ratio.abs() > f64::MIN_POSITIVE {
            x = ((a * ratio) + b) / denom;
            y = ((b * ratio) - a) / denom;
        } else {
            x = ((c * (a / d)) + b) / denom;
            y = ((c * (b / d)) - a) / denom;
        }
    } else {
        let ratio = d / c;
        let denom = (d * ratio) + c;
        if ratio.abs() > f64::MIN_POSITIVE {
            x = ((b * ratio) + a) / denom;
            y = (b - (a * ratio)) / denom;
        } else {
            x = ((d * (b / c)) + a) / denom;
            y = (b - (d * (a / c))) / denom;
        }
    }
    if x.is_nan() && y.is_nan() && c == 0.0 && d == 0.0 && (!a.is_nan() || !b.is_nan()) {
        let inf = f64::INFINITY.copysign(c);
        return (inf * a, inf * b);
    }
    (x, y)
}

impl ElmerScalar for Complex64 {
    const COMPONENTS: usize = 2;
    const DROP_INCLUSIVE: bool = false;
    const JACOBI_GUARD: bool = false;

    fn elmer_div(self, d: Complex64) -> Complex64 {
        let (re, im) = divdc3(self.re, self.im, d.re, d.im);
        Complex64::new(re, im)
    }
    fn for_each_component(self, f: &mut dyn FnMut(f64)) {
        f(self.re);
        f(self.im);
    }
    fn scale_real(self, s: f64) -> Complex64 {
        Complex64::new(self.re * s, self.im * s)
    }
    fn unscale_real(self, s: f64) -> Complex64 {
        Complex64::new(self.re / s, self.im / s)
    }
    fn diag_abs_terms(self, f: &mut dyn FnMut(f64)) {
        f(self.re.abs());
        f(self.re.abs());
    }

    fn diagonal_scaling(a: &CsrMatrix<Complex64>) -> Vec<f64> {
        (0..a.nrows)
            .map(|i| {
                let c = diag_index(a, i).map_or(Complex64::new(0.0, 0.0), |j| a.values[j]);
                let s = c.norm().sqrt();
                if s > TINY {
                    1.0 / s
                } else {
                    1.0
                }
            })
            .collect()
    }
}

/// gfortran's `NORM2` over a sequence of reals (scaled sum of squares).
pub fn norm2(vals: impl Iterator<Item = f64>) -> f64 {
    let mut scale = 1.0_f64;
    let mut result = 0.0_f64;
    for x in vals {
        if x != 0.0 {
            let abs_x = x.abs();
            if scale < abs_x {
                let val = scale / abs_x;
                result = 1.0 + result * val * val;
                scale = abs_x;
            } else {
                let val = abs_x / scale;
                result += val * val;
            }
        }
    }
    scale * result.sqrt()
}
