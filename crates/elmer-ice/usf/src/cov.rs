//! `tenpowera_`, `tenpowera_d_`, `derivative_tenpowera_`, `log10a_`,
//! `asquare_`, `asquare_d_`, `derivative_asquare_`, `sqrta_`,
//! `tenpowera_masked_`, `tenpowera_d_masked_`, `asquare_scaled_`,
//! `asquare_d_scaled_` — port of USF_CoV.F90 (variable changes for inverse
//! methods; beta=10^a or beta=a^2).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - The file USEs DefUtils only for types/dp; no Elmer call crosses.
//!   Entries ignore `Model`/`nodenumber` exactly where the source does.
//! - `**` evaluates in the Fortran accessor (`usf_pow_dp`): gfortran lowers
//!   real-exponent `**` to a libm `pow` call, and the reference was built
//!   at the same -O0; routing it through gfortran keeps the rounding
//!   identical by construction. `log`/`sqrt` are same-libm single calls.
//! - Masked variants branch on `mask <= 0` (grounded-mask convention:
//!   -1 floating, 0 grounding line, 1 grounded), verbatim.

use crate::support::*;

/// Shared entry prelude: the source reads no Elmer state; the handles
/// only fix the `(Model, nodenumber, slot) -> REAL(8)` ABI.
#[inline]
unsafe fn slot(v: *const c_double, i: usize) -> c_double {
    // SAFETY: caller (Elmer `ElmerRealFn`, `val(*)` assumed-size) passes
    // the address of at least i+1 doubles — the same contract the
    // reference relies on.
    unsafe { *v.add(i) }
}

/// `LOG(10.0)` with a default-REAL literal evaluates in SINGLE precision
/// in gfortran, then promotes: `(10.0f32).ln() as f64` is bit-identical,
/// `10.0f64.ln()` is not (verified against the reference entry).
#[inline]
fn ln10_sp() -> c_double {
    (10.0f32).ln() as c_double
}

pub unsafe fn tenpowera(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    unsafe { usf_pow_dp(10.0, slot(v, 0)) }
}

pub unsafe fn tenpowera_d(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    unsafe { usf_pow_dp(10.0, slot(v, 0)) * ln10_sp() }
}

pub unsafe fn derivative_tenpowera(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { slot(v, 0) * usf_pow_dp(10.0, slot(v, 1)) * ln10_sp() }
}

pub unsafe fn log10a(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    unsafe { slot(v, 0).log10() }
}

pub unsafe fn asquare(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let x = unsafe { slot(v, 0) };
    x * x
}

pub unsafe fn asquare_d(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    2.0 * unsafe { slot(v, 0) }
}

pub unsafe fn derivative_asquare(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    2.0 * unsafe { slot(v, 0) } * unsafe { slot(v, 1) }
}

pub unsafe fn sqrta(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    unsafe { slot(v, 0).sqrt() }
}

pub unsafe fn tenpowera_masked(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let (x, mask) = unsafe { (slot(v, 0), slot(v, 1)) };
    if mask <= 0.0 {
        0.0
    } else {
        unsafe { usf_pow_dp(10.0, x) }
    }
}

pub unsafe fn tenpowera_d_masked(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let (x, mask) = unsafe { (slot(v, 0), slot(v, 1)) };
    if mask <= 0.0 {
        0.0
    } else {
        unsafe { usf_pow_dp(10.0, x) * ln10_sp() }
    }
}

pub unsafe fn asquare_scaled(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let (x, s) = unsafe { (slot(v, 0), slot(v, 1)) };
    x * x * s
}

pub unsafe fn asquare_d_scaled(
    _m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let (x, s) = unsafe { (slot(v, 0), slot(v, 1)) };
    2.0 * x * s
}
