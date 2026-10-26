//! lane-fhutiter port: Elmer `huti_aux` (`fhutiter/src/huti_aux.F90` @ `9f6af2f85`)
//! through the gfortran ABI shim (ADR-027, T014/T017 admitted shape).
//!
//! Twelve auxiliary routines in four precisions (single, double, single
//! complex, double complex): dummy preconditioner copy (`*_dummy_pcondfun`),
//! pseudo-random fill (`*randvec`), and LU factorisation + solve without
//! pivoting (`*lusolve`, Saad Algorithm 10.4). All dummies are scalars by
//! reference or assumed-size arrays (raw base address, T017 shape B); the
//! module carries no descriptors, no `COMMON`, no allocatables.
//!
//! Calling convention pinned to gfortran 13.3.0 / x86-64 (`-O0`, the
//! campaign Debug default):
//! - complex division reproduces gfortran's inline Smith algorithm
//!   (decoded from `cdiv.s`; naive `(ac+bd)/d` diverges by ulps — measured
//!   `0x3f1be41c` vs naive `0x3f1be3xx` on the corpus matrix). NaN branch:
//!   `comiss`+`seta` takes the `.L2` path exactly when `!(|bi| > |br|)`,
//!   which Rust `>` reproduces (NaN compares false).
//! - complex multiply/subtract evaluation order is V1 (`t = a*b` with
//!   `re = ar*br - ai*bi`, then `c - t`), confirmed against `-O0` asm;
//!   no FMA contraction in the measured codegen (`mulss`/`subss`/`addss`).
//! - the small-pivot diagnostic reproduces gfortran list-directed bytes:
//!   scalar `E15.8`-class field (`{:>15.8E}`) for `REAL`, `E23.16`-class
//!   (`{:>23}` of a 3-digit-exponent expansion) for `DOUBLE PRECISION`,
//!   minimal `(re,im)` for complex with the same digit/exponent rules.
//! - `*randvec` fills from a documented xorshift64* stream (fixed seed):
//!   the orig `random_number` intrinsic is process-seeded and
//!   non-reproducible run to run, so stream identity is unachievable by
//!   construction; the corpus checks shape/range (and the observed
//!   pure-imaginary quirk of the complex variants, where the first harvest
//!   is overwritten: `u(i) = (0, h2)`), never the stream.
//!
//! Upstream license: LGPL-2.1 (see crate NOTICE / repo LICENSE).

use std::os::raw::c_int;
use std::sync::atomic::{AtomicU64, Ordering};

/// Fortran `COMPLEX` (2 x `REAL`): C interoperable, matches gfortran layout.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct C32 {
    pub re: f32,
    pub im: f32,
}

/// Fortran `DOUBLE COMPLEX` (2 x `DOUBLE PRECISION`).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

// ---------------------------------------------------------------------------
// Complex arithmetic: gfortran 13.3.0 -O0 emission, bit-exact.
// ---------------------------------------------------------------------------

/// Complex multiply, V1 order: `re = ar*br - ai*bi`, `im = ar*bi + ai*br`.
#[inline]
fn c32_mul(a: C32, b: C32) -> C32 {
    C32 { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }
}

/// gfortran inline Smith division (decoded `cdiv.s`).
#[inline]
fn c32_div(a: C32, b: C32) -> C32 {
    if b.im.abs() > b.re.abs() {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        C32 { re: (a.re * r + a.im) / den, im: (a.im * r - a.re) / den }
    } else {
        let r = b.im / b.re;
        let den = b.im * r + b.re;
        C32 { re: (a.im * r + a.re) / den, im: (a.im - a.re * r) / den }
    }
}

/// `c - a*b`, V1 order (matches `-O0` array-update codegen).
#[inline]
fn c32_msub(c: C32, a: C32, b: C32) -> C32 {
    let t = c32_mul(a, b);
    C32 { re: c.re - t.re, im: c.im - t.im }
}

/// Complex multiply, V1 order.
#[inline]
fn c64_mul(a: C64, b: C64) -> C64 {
    C64 { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }
}

/// gfortran inline Smith division, double precision.
#[inline]
fn c64_div(a: C64, b: C64) -> C64 {
    if b.im.abs() > b.re.abs() {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        C64 { re: (a.re * r + a.im) / den, im: (a.im * r - a.re) / den }
    } else {
        let r = b.im / b.re;
        let den = b.im * r + b.re;
        C64 { re: (a.im * r + a.re) / den, im: (a.im - a.re * r) / den }
    }
}

/// `c - a*b`, V1 order.
#[inline]
fn c64_msub(c: C64, a: C64, b: C64) -> C64 {
    let t = c64_mul(a, b);
    C64 { re: c.re - t.re, im: c.im - t.im }
}

// ---------------------------------------------------------------------------
// List-directed pivot diagnostic: gfortran byte emulation (pinned by corpus
// + direct format probes, `proof/complex-probe/{pivot,cfmt,dfmt,efmt}.F90`).
// Reachable pivot magnitudes are all < 1e-16 (the branch condition), so
// only the E-range + zero rules below are live; the F-range switch gfortran
// uses for moderate magnitudes is documented, not emulated.
// Measured rules (gfortran 13.3.0, list-directed):
// - scalar: field width 17 (`REAL`) / 26 (`DOUBLE PRECISION`) after `pivot`.
// - complex: field width 36 / 54, parts unpadded inside `(re,im)`.
// - nonzero parts: E with 8/16 decimals (scalar, 9/17 sig) but 9/17 decimals
//   (complex, 10/18 sig) — one extra digit for complex parts, measured.
// - zero parts: `0.00000000` / `0.0000000000000000` (sign kept for `-0.0`).
// - exponents padded to 2 digits (`REAL`: `E-09`) / 3 digits (`DOUBLE`).
// ---------------------------------------------------------------------------

/// Pad a Rust E-exponent (`E-8`) to gfortran width (`E-08` / `E-008`).
fn exp_pad(s: &str, width: usize) -> String {
    if let Some(e) = s.find('E') {
        let (mant, exp) = s.split_at(e + 1);
        if exp.len() >= 2 {
            let (sign, digits) = exp.split_at(1);
            if digits.len() < width {
                return format!("{mant}{sign}{:0>w$}", digits, w = width);
            }
        }
    }
    s.to_string()
}

/// `REAL(4)` scalar value: E with 8 decimals (E-range + zero only).
fn f32_scal(x: f32) -> String {
    if x == 0.0 {
        return (if x.is_sign_negative() { "-0.00000000" } else { "0.00000000" }).into();
    }
    exp_pad(&format!("{x:.8E}"), 2)
}

/// `DOUBLE PRECISION` scalar value: E with 16 decimals (E-range + zero only).
fn f64_scal(x: f64) -> String {
    if x == 0.0 {
        return (if x.is_sign_negative() { "-0.0000000000000000" } else { "0.0000000000000000" }).into();
    }
    exp_pad(&format!("{x:.16E}"), 3)
}

/// `COMPLEX` part: E with 9 decimals (E-range + zero only).
fn f32_part(x: f32) -> String {
    if x == 0.0 {
        return (if x.is_sign_negative() { "-0.00000000" } else { "0.00000000" }).into();
    }
    exp_pad(&format!("{x:.9E}"), 2)
}

/// `DOUBLE COMPLEX` part: E with 17 decimals (E-range + zero only).
fn f64_part(x: f64) -> String {
    if x == 0.0 {
        return (if x.is_sign_negative() { "-0.0000000000000000" } else { "0.0000000000000000" }).into();
    }
    exp_pad(&format!("{x:.17E}"), 3)
}

/// `REAL(4)` pivot line.
fn fmt_pivot_r(x: f32) -> String {
    format!(" (libhuti.a) GMRES: small pivot{:>17}", f32_scal(x))
}

/// `DOUBLE PRECISION` pivot line.
fn fmt_pivot_d(x: f64) -> String {
    format!(" (libhuti.a) GMRES: small pivot{:>26}", f64_scal(x))
}

/// `COMPLEX` pivot line.
fn fmt_pivot_c(x: C32) -> String {
    format!(" (libhuti.a) GMRES: small pivot{:>36}", format!("({},{})", f32_part(x.re), f32_part(x.im)))
}

/// `DOUBLE COMPLEX` pivot line.
fn fmt_pivot_z(x: C64) -> String {
    format!(" (libhuti.a) GMRES: small pivot{:>54}", format!("({},{})", f64_part(x.re), f64_part(x.im)))
}

// ---------------------------------------------------------------------------
// PRNG: xorshift64*, fixed seed. The orig `random_number` intrinsic is
// process-seeded (non-reproducible across runs), so no port can match its
// stream; this provides the same shape (count + `[0,1)` range) with a
// documented, process-deterministic filler.
// ---------------------------------------------------------------------------

/// Fixed seed (SplitMix-style nothing-up-my-sleeve constant).
static RNG_STATE: AtomicU64 = AtomicU64::new(0x9E3779B97F4A7C15);

fn next_u64() -> u64 {
    let mut x = RNG_STATE.load(Ordering::Relaxed);
    loop {
        let mut y = x;
        y ^= y >> 12;
        y ^= y << 25;
        y ^= y >> 27;
        match RNG_STATE.compare_exchange_weak(x, y, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => return y.wrapping_mul(0x2545_F491_4F6C_DD1D),
            Err(actual) => x = actual,
        }
    }
}

/// One `[0,1)` draw at single precision (high 24 bits / 2^24).
#[inline]
fn draw_f32() -> f32 {
    ((next_u64() >> 40) as f32) / 16777216.0
}

// ---------------------------------------------------------------------------
// HUTI_NDIM contract: `HUTI_NDIM` is `ipar(3)` (1-based) = `ipar[2]`.
// Assumed-size arrays arrive as raw base addresses; the valid prefix length
// is `ndim` elements (`workdim`-style contracts do not apply to huti_aux).
// ---------------------------------------------------------------------------

/// Read `HUTI_NDIM` (`ipar(3)`); non-positive means the Fortran loops skip.
#[inline]
unsafe fn ndim(ipar: *const c_int) -> usize {
    let n = unsafe { *ipar.add(2) };
    if n <= 0 {
        0
    } else {
        n as usize
    }
}

macro_rules! pcondfun {
    ($name:ident, $export:literal, $t:ty) => {
        /// Dummy preconditioner: pure copy `u(i) = v(i)`, `i = 1..ndim`.
        ///
        /// # Safety
        /// `u`/`v` must be valid for `ndim` element writes/reads and `ipar`
        /// valid for 3 reads.
        #[export_name = $export]
        pub unsafe extern "C" fn $name(u: *mut $t, v: *const $t, ipar: *const c_int) {
            let n = unsafe { ndim(ipar) };
            for i in 0..n {
                unsafe {
                    *u.add(i) = *v.add(i);
                }
            }
        }
    };
}

pcondfun!(port_sdummy_pcondfun, "__huti_aux_MOD_huti_sdummy_pcondfun", f32);
pcondfun!(port_ddummy_pcondfun, "__huti_aux_MOD_huti_ddummy_pcondfun", f64);
pcondfun!(port_cdummy_pcondfun, "__huti_aux_MOD_huti_cdummy_pcondfun", C32);
pcondfun!(port_zdummy_pcondfun, "__huti_aux_MOD_huti_zdummy_pcondfun", C64);

/// Single-precision random fill: `u(i) = harvest`, `harvest` in `[0,1)`.
///
/// # Safety
/// `u` must be valid for `ndim` writes and `ipar` valid for 3 reads.
#[export_name = "__huti_aux_MOD_huti_srandvec"]
pub unsafe extern "C" fn port_srandvec(u: *mut f32, ipar: *const c_int) {
    let n = unsafe { ndim(ipar) };
    for i in 0..n {
        unsafe {
            *u.add(i) = draw_f32();
        }
    }
}

/// Double-precision random fill (orig draws single-precision `harvest`).
///
/// # Safety
/// `u` must be valid for `ndim` writes and `ipar` valid for 3 reads.
#[export_name = "__huti_aux_MOD_huti_drandvec"]
pub unsafe extern "C" fn port_drandvec(u: *mut f64, ipar: *const c_int) {
    let n = unsafe { ndim(ipar) };
    for i in 0..n {
        unsafe {
            *u.add(i) = draw_f32() as f64;
        }
    }
}

/// Complex random fill. Observed upstream shape: the first harvest is
/// overwritten (`u(i) = (0, h2)`), so the real part is always `0.0` —
/// reproduced exactly (shape parity, not stream parity).
///
/// # Safety
/// `u` must be valid for `ndim` writes and `ipar` valid for 3 reads.
#[export_name = "__huti_aux_MOD_huti_crandvec"]
pub unsafe extern "C" fn port_crandvec(u: *mut C32, ipar: *const c_int) {
    let n = unsafe { ndim(ipar) };
    for i in 0..n {
        let h2 = draw_f32();
        unsafe {
            *u.add(i) = C32 { re: 0.0, im: h2 };
        }
    }
}

/// Double-complex random fill (same observed pure-imaginary shape).
///
/// # Safety
/// `u` must be valid for `ndim` writes and `ipar` valid for 3 reads.
#[export_name = "__huti_aux_MOD_huti_zrandvec"]
pub unsafe extern "C" fn port_zrandvec(u: *mut C64, ipar: *const c_int) {
    let n = unsafe { ndim(ipar) };
    for i in 0..n {
        let h2 = draw_f32();
        unsafe {
            *u.add(i) = C64 { re: 0.0, im: h2 as f64 };
        }
    }
}

macro_rules! lusolve_real {
    ($name:ident, $export:literal, $t:ty, $fmt:ident, $pivot:literal) => {
        /// LU factorisation + forward/backward solve, no pivoting
        /// (Saad Algorithm 10.4). Column-major `n`x`n`; prints the
        /// small-pivot diagnostic when `|lumat(k,k)| < 1e-16`.
        ///
        /// # Safety
        /// `lumat` must be valid for `n*n` read/writes, `u` for `n`
        /// writes, `v` for `n` reads, `n` points to one readable `c_int`.
        #[export_name = $export]
        pub unsafe extern "C" fn $name(
            n_: *const c_int,
            lumat: *mut $t,
            u: *mut $t,
            v: *const $t,
        ) {
            let n = unsafe { *n_ };
            if n <= 0 {
                return;
            }
            let n = n as usize;
            let at = |i: usize, j: usize| j * n + i;
            for i in 1..n {
                for k in 0..i {
                    let kk = unsafe { *lumat.add(at(k, k)) };
                    if kk.abs() < $pivot {
                        println!("{}", $fmt(kk));
                    }
                    unsafe {
                        *lumat.add(at(i, k)) = *lumat.add(at(i, k)) / kk;
                    }
                    for j in (k + 1)..n {
                        unsafe {
                            let aik = *lumat.add(at(i, k));
                            *lumat.add(at(i, j)) = *lumat.add(at(i, j)) - aik * *lumat.add(at(k, j));
                        }
                    }
                }
            }
            for i in 0..n {
                unsafe {
                    *u.add(i) = *v.add(i);
                    for k in 0..i {
                        *u.add(i) = *u.add(i) - *lumat.add(at(i, k)) * *u.add(k);
                    }
                }
            }
            for i in (0..n).rev() {
                unsafe {
                    for k in (i + 1)..n {
                        *u.add(i) = *u.add(i) - *lumat.add(at(i, k)) * *u.add(k);
                    }
                    *u.add(i) = *u.add(i) / *lumat.add(at(i, i));
                }
            }
        }
    };
}

lusolve_real!(
    port_slusolve,
    "__huti_aux_MOD_huti_slusolve",
    f32,
    fmt_pivot_r,
    1.0e-16f32
);
lusolve_real!(
    port_dlusolve,
    "__huti_aux_MOD_huti_dlusolve",
    f64,
    fmt_pivot_d,
    1.0e-16f64
);

/// Complex LU solve (Smith divisions, V1 updates).
///
/// # Safety
/// Same contract as the real variants with `C32` elements.
#[export_name = "__huti_aux_MOD_huti_clusolve"]
pub unsafe extern "C" fn port_clusolve(
    n_: *const c_int,
    lumat: *mut C32,
    u: *mut C32,
    v: *const C32,
) {
    let n = unsafe { *n_ };
    if n <= 0 {
        return;
    }
    let n = n as usize;
    let at = |i: usize, j: usize| j * n + i;
    unsafe {
        for i in 1..n {
            for k in 0..i {
                let kk = *lumat.add(at(k, k));
                if kk.re.hypot(kk.im) < 1.0e-16f32 {
                    println!("{}", fmt_pivot_c(kk));
                }
                *lumat.add(at(i, k)) = c32_div(*lumat.add(at(i, k)), kk);
                for j in (k + 1)..n {
                    let aik = *lumat.add(at(i, k));
                    *lumat.add(at(i, j)) = c32_msub(*lumat.add(at(i, j)), aik, *lumat.add(at(k, j)));
                }
            }
        }
        for i in 0..n {
            *u.add(i) = *v.add(i);
            for k in 0..i {
                *u.add(i) = c32_msub(*u.add(i), *lumat.add(at(i, k)), *u.add(k));
            }
        }
        for i in (0..n).rev() {
            for k in (i + 1)..n {
                *u.add(i) = c32_msub(*u.add(i), *lumat.add(at(i, k)), *u.add(k));
            }
            *u.add(i) = c32_div(*u.add(i), *lumat.add(at(i, i)));
        }
    }
}

/// Double-complex LU solve (Smith divisions, V1 updates).
///
/// # Safety
/// Same contract as the real variants with `C64` elements.
#[export_name = "__huti_aux_MOD_huti_zlusolve"]
pub unsafe extern "C" fn port_zlusolve(
    n_: *const c_int,
    lumat: *mut C64,
    u: *mut C64,
    v: *const C64,
) {
    let n = unsafe { *n_ };
    if n <= 0 {
        return;
    }
    let n = n as usize;
    let at = |i: usize, j: usize| j * n + i;
    unsafe {
        for i in 1..n {
            for k in 0..i {
                let kk = *lumat.add(at(k, k));
                if kk.re.hypot(kk.im) < 1.0e-16f64 {
                    println!("{}", fmt_pivot_z(kk));
                }
                *lumat.add(at(i, k)) = c64_div(*lumat.add(at(i, k)), kk);
                for j in (k + 1)..n {
                    let aik = *lumat.add(at(i, k));
                    *lumat.add(at(i, j)) = c64_msub(*lumat.add(at(i, j)), aik, *lumat.add(at(k, j)));
                }
            }
        }
        for i in 0..n {
            *u.add(i) = *v.add(i);
            for k in 0..i {
                *u.add(i) = c64_msub(*u.add(i), *lumat.add(at(i, k)), *u.add(k));
            }
        }
        for i in (0..n).rev() {
            for k in (i + 1)..n {
                *u.add(i) = c64_msub(*u.add(i), *lumat.add(at(i, k)), *u.add(k));
            }
            *u.add(i) = c64_div(*u.add(i), *lumat.add(at(i, i)));
        }
    }
}
