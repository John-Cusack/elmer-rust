//! lane-fem-b port: Elmer `HarmUtils` (`SinSum`/`CosSum`,
//! `fem/src/modules/HarmonicSignals.F90`) through the gfortran ABI shim
//! (ADR-027, T014 admitted shape). Clean-room reimplementation of two tiny
//! numerical kernels; calling convention pinned to gfortran 13.3.0 / x86-64.
//! Only the `BT_REAL` (`REAL(KIND=dp)`, 8-byte) dtype is admitted by the
//! gate; the descriptor layout below is the T003-measured one.
//!
//! Upstream semantics (GPL-2.0-or-later, see crate NOTICE):
//! `sumA = sum over rows of amplitudes(i,2) * sqrt(2) * sin/cos(
//! amplitudes(i,1) * 2 * PI * fundamental_f * t + phase)`, iterating
//! `SIZE(amplitudes(:,1))` rows.

/// One dimension of a gfortran array descriptor: element stride (signed),
/// lower bound, upper bound.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct GfDim {
    pub stride: i64,
    pub lower: i64,
    pub upper: i64,
}

/// gfortran dtype word: version 0, rank, BT_* type code, attribute.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct GfDtype {
    pub version: i32,
    pub rank: u8,
    pub type_code: u8,
    pub attribute: u16,
}

/// gfortran array descriptor for `REAL(8)` rank 2: 40-byte header +
/// 24 bytes per dimension (ADR-027 §2). The callee never retains it.
#[repr(C)]
pub struct GfDesc2 {
    pub base_addr: *const f64,
    pub offset: i64,
    pub elem_size: i64,
    pub dtype: GfDtype,
    pub spare: i64,
    pub dims: [GfDim; 2],
}

const SQRT2: f64 = std::f64::consts::SQRT_2;

fn elem(d: &GfDesc2, i: i64, j: i64) -> f64 {
    let idx = d.offset + i * d.dims[0].stride + j * d.dims[1].stride;
    unsafe { *d.base_addr.offset(idx as isize) }
}

fn rows(d: &GfDesc2) -> std::ops::RangeInclusive<i64> {
    d.dims[0].lower..=d.dims[0].upper
}

/// `SinSum`: sum over rows of `amp(i,2)*sqrt(2)*sin(amp(i,1)*w*t+phase)`.
/// Scalar dummies travel by reference; `amps` is an assumed-shape
/// descriptor (whole arrays and sections alike: base + offset + strides).
#[export_name = "__harmutils_MOD_sinsum"]
pub extern "C" fn port_sinsum(
    f: *const f64,
    amps: *const GfDesc2,
    t: *const f64,
    phase: *const f64,
) -> f64 {
    let d = unsafe { &*amps };
    debug_assert_eq!(d.elem_size, 8);
    debug_assert_eq!(d.dtype.version, 0);
    debug_assert_eq!(d.dtype.rank, 2);
    debug_assert_eq!(d.dtype.type_code, 3);
    let (f, t, phase) = unsafe { (*f, *t, *phase) };
    let w = 2.0 * std::f64::consts::PI * f;
    let mut s = 0.0;
    for i in rows(d) {
        s += elem(d, i, 2) * SQRT2 * (elem(d, i, 1) * w * t + phase).sin();
    }
    s
}

/// `CosSum`: same with `cos`.
#[export_name = "__harmutils_MOD_cossum"]
pub extern "C" fn port_cossum(
    f: *const f64,
    amps: *const GfDesc2,
    t: *const f64,
    phase: *const f64,
) -> f64 {
    let d = unsafe { &*amps };
    debug_assert_eq!(d.elem_size, 8);
    debug_assert_eq!(d.dtype.version, 0);
    debug_assert_eq!(d.dtype.rank, 2);
    debug_assert_eq!(d.dtype.type_code, 3);
    let (f, t, phase) = unsafe { (*f, *t, *phase) };
    let w = 2.0 * std::f64::consts::PI * f;
    let mut s = 0.0;
    for i in rows(d) {
        s += elem(d, i, 2) * SQRT2 * (elem(d, i, 1) * w * t + phase).cos();
    }
    s
}
