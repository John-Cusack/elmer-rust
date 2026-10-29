//! Scalar types the direct solvers run over.
//!
//! `elmer_core::Scalar` is deliberately minimal; the factorisations here also
//! need faer's `ComplexField` (for the sparse and dense paths) and the two
//! LAPACK-specific measures the band solvers use. The trait is sealed: it is
//! implemented for exactly `f64` and `Complex64`.

use elmer_core::Scalar;
use faer::traits::ComplexField;
use num_complex::Complex64;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex64 {}
}

/// `f64` or `Complex64`, usable by every solver in this crate.
pub trait DirectScalar: Scalar + ComplexField + sealed::Sealed {
    /// LAPACK `CABS1`: `|re| + |im|`, the measure `IZAMAX` pivots on
    /// (`IDAMAX` for real data).
    fn abs1(self) -> f64;
    /// Real part (the diagonal of a Hermitian factor is real).
    fn re(self) -> f64;
}

impl DirectScalar for f64 {
    fn abs1(self) -> f64 {
        self.abs()
    }
    fn re(self) -> f64 {
        self
    }
}

impl DirectScalar for Complex64 {
    fn abs1(self) -> f64 {
        self.re.abs() + self.im.abs()
    }
    fn re(self) -> f64 {
        self.re
    }
}
