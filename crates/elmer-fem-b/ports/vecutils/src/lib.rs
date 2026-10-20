//! Port of `fem/src/modules/CoordinateTransform.F90`, module `VecUtils`
//! (Elmer @ 9f6af2f85) for lane-fem-b, through the gfortran ABI shim
//! (ADR-027, T017 shape B).
//!
//! Unit: `fortran:fem/src/modules/CoordinateTransform.F90#vecutils` — four
//! module procedures whose array dummies are **explicit-shape**
//! (`a(3)`, `OldCoord(3,3)`, `A(3,3)`, `jac(3,3)`).
//!
//! Wire (measured on the call site and the callee, ADR-050):
//! - a module procedure is reached through `__<mod>_MOD_<proc>`;
//! - an explicit-shape array dummy travels as a **raw base address**, in
//!   Fortran's **column-major** order;
//! - these four functions return **arrays**, so the result travels back
//!   through a hidden **first** argument: a gfortran **array
//!   descriptor** (base address at offset 0, the per-dimension stride at
//!   offset 40), exactly the descriptor layout of ADR-027 §2. The callee
//!   writes element `i` at `base + (i-1) * stride`.
//!
//! Faithfulness notes:
//! - Summation and accumulation order is kept verbatim (`SUM(a**2)`,
//!   `DOT_PRODUCT`, and the `k`/`l` loops of `transform2rank`); the
//!   equivalence corpus grades bit-exact output.
//! - `x**2` is gfortran's `__powidf2` (repeated multiplication), not
//!   `pow(x, 2)`; the two agree here because the square of a double is
//!   correctly rounded, but the port multiplies so the order matches.

use std::os::raw::c_double;

/// One dimension of a gfortran array descriptor: element stride (signed),
/// lower bound, upper bound.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct GfDim {
    pub stride: i64,
    pub lower: i64,
    pub upper: i64,
}

/// gfortran dtype word: version, rank, BT_* type code, attribute.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct GfDtype {
    pub version: i32,
    pub rank: u8,
    pub type_code: u8,
    pub attribute: u16,
}

/// gfortran array descriptor, rank 1 (ADR-027 §2): 40-byte header plus
/// 24 bytes per dimension.
#[repr(C)]
pub struct GfDesc1 {
    pub base_addr: *const f64,
    pub offset: i64,
    pub elem_size: i64,
    pub dtype: GfDtype,
    pub spare: i64,
    pub dims: [GfDim; 1],
}

/// gfortran array descriptor, rank 2.
#[repr(C)]
pub struct GfDesc2 {
    pub base_addr: *const f64,
    pub offset: i64,
    pub elem_size: i64,
    pub dtype: GfDtype,
    pub spare: i64,
    pub dims: [GfDim; 2],
}

/// Column-major element of a 3x3 array (Fortran storage), 0-based.
#[inline]
fn at3x3(a: *const c_double, i: usize, j: usize) -> f64 {
    unsafe { *a.add(i + 3 * j) }
}

/// `NORM2`: `norm2 = sqrt(SUM(a**2))`, assigned to all three elements.
#[export_name = "__vecutils_MOD_norm2"]
pub extern "C" fn port_norm2(result: *const GfDesc1, a: *const c_double) {
    let d = unsafe { &*result };
    let base = d.base_addr as *mut c_double;
    let stride = d.dims[0].stride;
    let mut s = 0.0f64;
    for k in 0..3 {
        let v = unsafe { *a.add(k) };
        s += v * v;
    }
    let n = s.sqrt();
    for k in 0..3 {
        unsafe {
            *base.offset((k as isize) * (stride as isize)) = n;
        }
    }
}

/// `NORMALIZED`: `normalized = a/norm2(a)`.
#[export_name = "__vecutils_MOD_normalized"]
pub extern "C" fn port_normalized(result: *const GfDesc1, a: *const c_double) {
    let d = unsafe { &*result };
    let base = d.base_addr as *mut c_double;
    let stride = d.dims[0].stride;
    let mut s = 0.0f64;
    for k in 0..3 {
        let v = unsafe { *a.add(k) };
        s += v * v;
    }
    let n = s.sqrt();
    for k in 0..3 {
        let v = unsafe { *a.add(k) };
        unsafe {
            *base.offset((k as isize) * (stride as isize)) = v / n;
        }
    }
}

/// `JAC`: the Jacobian of the coordinate transform between two systems,
/// `jac(i,j) = DOT_PRODUCT(NewCoord(i,:), OldCoord(j,:))
///             / sqrt(DOT_PRODUCT(OldCoord(j,:), OldCoord(j,:)))`.
#[export_name = "__vecutils_MOD_jac"]
pub extern "C" fn port_jac(
    result: *const GfDesc2,
    oldcoord: *const c_double,
    newcoord: *const c_double,
) {
    let d = unsafe { &*result };
    let base = d.base_addr as *mut c_double;
    let stride0 = d.dims[0].stride;
    let stride1 = d.dims[1].stride;
    for i in 0..3 {
        for j in 0..3 {
            let mut num = 0.0f64;
            for k in 0..3 {
                num += at3x3(newcoord, i, k) * at3x3(oldcoord, j, k);
            }
            let mut den = 0.0f64;
            for k in 0..3 {
                let v = at3x3(oldcoord, j, k);
                den += v * v;
            }
            let v = num / den.sqrt();
            unsafe {
                *base.offset((i as isize) * (stride0 as isize) + (j as isize) * (stride1 as isize)) = v;
            }
        }
    }
}

/// `TRANSFORM2RANK`: `B(i,j) = sum_k sum_l jac(i,l)*jac(j,k)*A(k,l)`.
#[export_name = "__vecutils_MOD_transform2rank"]
pub extern "C" fn port_transform2rank(
    result: *const GfDesc2,
    a: *const c_double,
    jac: *const c_double,
) {
    let d = unsafe { &*result };
    let base = d.base_addr as *mut c_double;
    let stride0 = d.dims[0].stride;
    let stride1 = d.dims[1].stride;
    for i in 0..3 {
        for j in 0..3 {
            let mut acc = 0.0f64;
            for k in 0..3 {
                for l in 0..3 {
                    acc += at3x3(jac, i, l) * at3x3(jac, j, k) * at3x3(a, k, l);
                }
            }
            unsafe {
                *base.offset((i as isize) * (stride0 as isize) + (j as isize) * (stride1 as isize)) = acc;
            }
        }
    }
}
