//! The scalar field a matrix or vector is built over (real or complex).

use num_complex::Complex64;

/// Real (`f64`) or complex (`Complex64`) entry type, as Elmer's real and
/// complex matrix paths.
pub trait Scalar:
    Copy
    + Send
    + Sync
    + core::fmt::Debug
    + PartialEq
    + 'static
    + core::ops::Add<Output = Self>
    + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self>
    + core::ops::Div<Output = Self>
    + core::ops::Neg<Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn conj(self) -> Self;
    /// Modulus (`|x|`), always real.
    fn abs(self) -> f64;
    fn from_f64(x: f64) -> Self;
}

impl Scalar for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn conj(self) -> Self {
        self
    }
    fn abs(self) -> f64 {
        f64::abs(self)
    }
    fn from_f64(x: f64) -> Self {
        x
    }
}

impl Scalar for Complex64 {
    fn zero() -> Self {
        Complex64::new(0.0, 0.0)
    }
    fn one() -> Self {
        Complex64::new(1.0, 0.0)
    }
    fn conj(self) -> Self {
        Complex64::conj(&self)
    }
    fn abs(self) -> f64 {
        self.norm()
    }
    fn from_f64(x: f64) -> Self {
        Complex64::new(x, 0.0)
    }
}
