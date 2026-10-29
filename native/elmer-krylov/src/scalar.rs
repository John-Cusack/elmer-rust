//! Scalar arithmetic as the reference HUTI build evaluates it.
//!
//! The reference libfhuti is gfortran at `-O0` linked to reference BLAS 3.12,
//! so every primitive here mirrors one of those: plain left-to-right sums for
//! `Otmp_ddot`/`Otmp_zdotc`/`Otmp_zdotu` (IterSolve.F90:370-545), Blue's
//! algorithm for `dnrm2`/`dznrm2`, Smith's range-reduced complex division
//! (`-fcx-fortran-rules`) and glibc `csqrt`.

use elmer_core::Scalar;
use num_complex::Complex64;

/// What HUTI needs from a scalar beyond [`Scalar`].
pub trait HutiScalar: Scalar {
    /// `true` for the `z` (double complex) HUTI variants.
    const COMPLEX: bool;

    /// Fortran `a / b` as gfortran compiles it.
    fn fdiv(self, d: Self) -> Self;

    /// Fortran `SQRT`.
    fn fsqrt(self) -> Self;

    /// Fortran `REAL(a)` without `KIND`: a real argument becomes default
    /// (single precision) real, a complex one its real part at full precision.
    fn freal(self) -> f64;

    /// `(re, im)`; `im` is 0 for reals.
    fn parts(self) -> [f64; 2];

    fn from_parts(p: [f64; 2]) -> Self;

    /// Sequential `sum conj(x_i) * y_i` (`Otmp_ddot` / `Otmp_zdotc`).
    fn dotc(x: &[Self], y: &[Self]) -> Self;

    /// Sequential `sum x_i * y_i` (`Otmp_zdotu`, used by complex CG).
    fn dotu(x: &[Self], y: &[Self]) -> Self;

    /// Reference BLAS `dnrm2` / `dznrm2`.
    fn nrm2(x: &[Self]) -> f64;

    /// Local `sum |x_i|^2`, for the multi-rank norm.
    fn sumsq(x: &[Self]) -> f64;

    /// `acc + self * y` as one step of libgfortran's `MATMUL` inner loop. The
    /// library dispatches to a variant built for AVX2/AVX-512 with FMA
    /// (`matmul_r8_avx2`), where GCC contracts the real update into a fused
    /// multiply-add.
    fn mul_acc(self, y: Self, acc: Self) -> Self;

    fn is_zero(self) -> bool {
        self == Self::zero()
    }
}

impl HutiScalar for f64 {
    const COMPLEX: bool = false;

    fn fdiv(self, d: Self) -> Self {
        self / d
    }

    fn fsqrt(self) -> Self {
        self.sqrt()
    }

    fn freal(self) -> f64 {
        self as f32 as f64
    }

    fn parts(self) -> [f64; 2] {
        [self, 0.0]
    }

    fn from_parts(p: [f64; 2]) -> Self {
        p[0]
    }

    fn dotc(x: &[Self], y: &[Self]) -> Self {
        let mut s = 0.0;
        for i in 0..x.len() {
            s += x[i] * y[i];
        }
        s
    }

    fn dotu(x: &[Self], y: &[Self]) -> Self {
        Self::dotc(x, y)
    }

    fn nrm2(x: &[Self]) -> f64 {
        blue_nrm2(x.iter().copied())
    }

    fn sumsq(x: &[Self]) -> f64 {
        x.iter().fold(0.0, |s, v| s + v * v)
    }

    fn mul_acc(self, y: Self, acc: Self) -> Self {
        self.mul_add(y, acc)
    }
}

impl HutiScalar for Complex64 {
    const COMPLEX: bool = true;

    fn fdiv(self, d: Self) -> Self {
        smith_div(self, d)
    }

    fn fsqrt(self) -> Self {
        glibc_csqrt(self)
    }

    fn freal(self) -> f64 {
        self.re
    }

    fn parts(self) -> [f64; 2] {
        [self.re, self.im]
    }

    fn from_parts(p: [f64; 2]) -> Self {
        Complex64::new(p[0], p[1])
    }

    fn dotc(x: &[Self], y: &[Self]) -> Self {
        let mut s = Complex64::new(0.0, 0.0);
        for i in 0..x.len() {
            s += x[i].conj() * y[i];
        }
        s
    }

    fn dotu(x: &[Self], y: &[Self]) -> Self {
        let mut s = Complex64::new(0.0, 0.0);
        for i in 0..x.len() {
            s += x[i] * y[i];
        }
        s
    }

    fn nrm2(x: &[Self]) -> f64 {
        blue_nrm2(x.iter().flat_map(|z| [z.re, z.im]))
    }

    fn sumsq(x: &[Self]) -> f64 {
        x.iter().fold(0.0, |s, z| s + z.re * z.re + z.im * z.im)
    }

    /// How GCC contracts `matmul_c8_avx2` has not been checked against the
    /// reference; this is the unfused complex product.
    fn mul_acc(self, y: Self, acc: Self) -> Self {
        acc + self * y
    }
}

/// Blue's scaled sum of squares, as reference BLAS 3.10+ `dnrm2.f90` (and
/// `dznrm2.f90`, which feeds real then imaginary part of each entry).
fn blue_nrm2(xs: impl Iterator<Item = f64>) -> f64 {
    // Constants for IEEE double, as the .f90 derives them from radix/digits.
    let tsml = 2f64.powi(-511);
    let tbig = 2f64.powi(486);
    let ssml = 2f64.powi(537);
    let sbig = 2f64.powi(-538);
    let maxn = f64::MAX;

    let mut notbig = true;
    let (mut asml, mut amed, mut abig) = (0.0, 0.0, 0.0);
    for x in xs {
        let ax = x.abs();
        if ax > tbig {
            abig += (ax * sbig) * (ax * sbig);
            notbig = false;
        } else if ax < tsml {
            if notbig {
                asml += (ax * ssml) * (ax * ssml);
            }
        } else {
            amed += ax * ax;
        }
    }
    let (scl, sumsq);
    if abig > 0.0 {
        if amed > 0.0 || amed > maxn || amed.is_nan() {
            abig += (amed * sbig) * sbig;
        }
        scl = 1.0 / sbig;
        sumsq = abig;
    } else if asml > 0.0 {
        if amed > 0.0 || amed > maxn || amed.is_nan() {
            let amed = amed.sqrt();
            let asml = asml.sqrt() / ssml;
            let (ymin, ymax) = if asml > amed {
                (amed, asml)
            } else {
                (asml, amed)
            };
            scl = 1.0;
            sumsq = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax));
        } else {
            scl = 1.0 / ssml;
            sumsq = asml;
        }
    } else {
        scl = 1.0;
        sumsq = amed;
    }
    scl * sumsq.sqrt()
}

/// GCC's inline complex division under `-fcx-fortran-rules`
/// (`expand_complex_div_wide`, Smith's method without NaN recovery).
fn smith_div(a: Complex64, b: Complex64) -> Complex64 {
    let (ar, ai, br, bi) = (a.re, a.im, b.re, b.im);
    if br.abs() < bi.abs() {
        let ratio = br / bi;
        let div = br * ratio + bi;
        Complex64::new((ar * ratio + ai) / div, (ai * ratio - ar) / div)
    } else {
        let ratio = bi / br;
        let div = bi * ratio + br;
        Complex64::new((ai * ratio + ar) / div, (ai - ar * ratio) / div)
    }
}

/// glibc `csqrt` (`s_csqrt_template.c`) for finite arguments; non-finite
/// arguments only propagate NaN/Inf, which HUTI then reports as divergence.
fn glibc_csqrt(x: Complex64) -> Complex64 {
    let (mut re, mut im) = (x.re, x.im);
    if !re.is_finite() || !im.is_finite() {
        return x.sqrt();
    }
    if im == 0.0 {
        return if re < 0.0 {
            Complex64::new(0.0, (-re).sqrt().copysign(im))
        } else {
            Complex64::new(re.sqrt().abs(), 0f64.copysign(im))
        };
    }
    if re == 0.0 {
        let r = if im.abs() >= 2.0 * f64::MIN_POSITIVE {
            (0.5 * im.abs()).sqrt()
        } else {
            0.5 * (2.0 * im.abs()).sqrt()
        };
        return Complex64::new(r, r.copysign(im));
    }
    let mut scale: i32 = 0;
    if re.abs() > f64::MAX / 4.0 {
        scale = 1;
        re *= 0.25;
        im *= 0.25;
    } else if im.abs() > f64::MAX / 4.0 {
        scale = 1;
        re = if re.abs() >= 4.0 * f64::MIN_POSITIVE {
            re * 0.25
        } else {
            0.0
        };
        im *= 0.25;
    } else if re.abs() < 2.0 * f64::MIN_POSITIVE && im.abs() < 2.0 * f64::MIN_POSITIVE {
        scale = -((f64::MANTISSA_DIGITS as i32 + 1) / 2);
        re *= 2f64.powi(-2 * scale);
        im *= 2f64.powi(-2 * scale);
    }
    let d = re.hypot(im);
    let (mut r, mut s);
    if re > 0.0 {
        r = (0.5 * (d + re)).sqrt();
        if scale == 1 && im.abs() < 1.0 {
            s = im / r;
            r *= 2f64.powi(scale);
            scale = 0;
        } else {
            s = 0.5 * (im / r);
        }
    } else {
        s = (0.5 * (d - re)).sqrt();
        if scale == 1 && im.abs() < 1.0 {
            r = (im / s).abs();
            s *= 2f64.powi(scale);
            scale = 0;
        } else {
            r = (0.5 * (im / s)).abs();
        }
    }
    if scale != 0 {
        r *= 2f64.powi(scale);
        s *= 2f64.powi(scale);
    }
    Complex64::new(r, s.copysign(im))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nrm2_matches_plain_sum_in_mid_range() {
        let x = [3.0, -4.0, 12.0];
        assert_eq!(f64::nrm2(&x), 13.0);
        let z = [Complex64::new(3.0, 4.0), Complex64::new(0.0, 12.0)];
        assert_eq!(Complex64::nrm2(&z), 13.0);
    }

    #[test]
    fn nrm2_survives_tiny_and_huge_entries() {
        // A naive sum of squares underflows to 0 / overflows to inf here.
        let tiny = f64::nrm2(&[3e-300, 4e-300]);
        assert!((tiny / 5e-300 - 1.0).abs() < 1e-15);
        let big = f64::nrm2(&[3e300, 4e300]);
        assert!((big / 5e300 - 1.0).abs() < 1e-15);
    }

    #[test]
    fn smith_division_is_exact_on_simple_quotients() {
        let q = smith_div(Complex64::new(1.0, 1.0), Complex64::new(0.0, 2.0));
        assert_eq!(q, Complex64::new(0.5, -0.5));
        // Dividing by a real promotes it, and is then componentwise.
        let q = smith_div(Complex64::new(3.0, -6.0), Complex64::new(3.0, 0.0));
        assert_eq!(q, Complex64::new(1.0, -2.0));
    }

    #[test]
    fn csqrt_principal_branch() {
        assert_eq!(glibc_csqrt(Complex64::new(-4.0, 0.0)), Complex64::new(0.0, 2.0));
        assert_eq!(glibc_csqrt(Complex64::new(0.0, 2.0)), Complex64::new(1.0, 1.0));
        let z = glibc_csqrt(Complex64::new(3.0, 4.0));
        assert_eq!(z, Complex64::new(2.0, 1.0));
        let z = glibc_csqrt(Complex64::new(-3.0, -4.0));
        assert_eq!(z, Complex64::new(1.0, -2.0));
    }

    #[test]
    fn freal_rounds_reals_to_single() {
        assert_eq!(0.1f64.freal(), 0.1f32 as f64);
        assert_eq!(Complex64::new(0.1, 7.0).freal(), 0.1);
    }
}
