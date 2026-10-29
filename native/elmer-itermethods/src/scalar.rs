//! Real/complex differences the Fortran methods rely on, on top of `elmer_core::Scalar`.

use elmer_core::Scalar;
use num_complex::Complex64;

/// `Scalar` plus what the real (`itermethod_*`) and complex (`itermethod_z_*`)
/// code paths need to tell apart. Implemented for `f64` and `Complex64` only.
pub trait IterScalar: Scalar {
    /// Selects the complex variant's quirks (no NaN fatals, `|.|` for the
    /// small-system quadratic forms, no robust mode, ...).
    const IS_COMPLEX: bool;
    fn re(self) -> f64;
    fn im(self) -> f64;
    /// Real scalars drop `im`.
    fn from_parts(re: f64, im: f64) -> Self;
    /// `self * r` with a real factor, component-wise like Fortran's mixed-mode product.
    fn scale(self, r: f64) -> Self;
    /// `self / r` with a real divisor, component-wise.
    fn unscale(self, r: f64) -> Self;
    fn is_nan(self) -> bool;
}

impl IterScalar for f64 {
    const IS_COMPLEX: bool = false;
    fn re(self) -> f64 {
        self
    }
    fn im(self) -> f64 {
        0.0
    }
    fn from_parts(re: f64, _im: f64) -> Self {
        re
    }
    fn scale(self, r: f64) -> Self {
        self * r
    }
    fn unscale(self, r: f64) -> Self {
        self / r
    }
    fn is_nan(self) -> bool {
        f64::is_nan(self)
    }
}

impl IterScalar for Complex64 {
    const IS_COMPLEX: bool = true;
    fn re(self) -> f64 {
        self.re
    }
    fn im(self) -> f64 {
        self.im
    }
    fn from_parts(re: f64, im: f64) -> Self {
        Complex64::new(re, im)
    }
    fn scale(self, r: f64) -> Self {
        Complex64::new(self.re * r, self.im * r)
    }
    fn unscale(self, r: f64) -> Self {
        Complex64::new(self.re / r, self.im / r)
    }
    fn is_nan(self) -> bool {
        self.re.is_nan() || self.im.is_nan()
    }
}
