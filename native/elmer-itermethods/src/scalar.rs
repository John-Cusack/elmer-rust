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
    /// `self / d` as gfortran evaluates it: the plain quotient for reals;
    /// for complex operands the range-reduced (Smith) division GCC expands
    /// inline under the Fortran rules (`-fcx-fortran-rules`, the gfortran
    /// default), not `num_complex`'s `(self * conj(d)) / |d|^2`.
    fn fdiv(self, d: Self) -> Self;
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
    fn fdiv(self, d: Self) -> Self {
        self / d
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
    fn fdiv(self, d: Self) -> Self {
        // tree-complex.cc expand_complex_div_wide.
        let (ar, ai, br, bi) = (self.re, self.im, d.re, d.im);
        if br.abs() < bi.abs() {
            let ratio = br / bi;
            let div = (br * ratio) + bi;
            let tr = (ar * ratio) + ai;
            let ti = (ai * ratio) - ar;
            Complex64::new(tr / div, ti / div)
        } else {
            let ratio = bi / br;
            let div = (bi * ratio) + br;
            let tr = (ai * ratio) + ar;
            let ti = ai - (ar * ratio);
            Complex64::new(tr / div, ti / div)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Quotients gfortran 13 (-O3) prints for the same operands (`Z17`).
    #[test]
    fn complex_division_rounds_like_gfortran() {
        let q = Complex64::new(0.3, 0.7).fdiv(Complex64::new(1.1, 1e-3));
        assert_eq!((q.re.to_bits(), q.im.to_bits()), (0x3FD17DD699B33D82, 0x3FE45B0E371624DB));
        let q = Complex64::new(0.123456789, -9.87654321).fdiv(Complex64::new(3.3, 7.7));
        assert_eq!((q.re.to_bits(), q.im.to_bits()), (0xBFF13EC8AAAA27B6, 0xBFDE96E42B591715));
        // The naive quotient differs in the last bit here.
        let n = Complex64::new(0.123456789, -9.87654321) / Complex64::new(3.3, 7.7);
        assert_ne!(n.re.to_bits(), 0xBFF13EC8AAAA27B6);
        assert_eq!(2.0f64.fdiv(8.0), 0.25);
    }
}
