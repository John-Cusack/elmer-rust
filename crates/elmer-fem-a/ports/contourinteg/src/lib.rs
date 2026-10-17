#![allow(static_mut_refs)]
//! `static mut` mirrors C file-static mutable state.

// Port of `fem/src/view3d/ContourInteg.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/ContourInteg.c.
//
// True object symbols: `ContourCountInit`, `ContourCountSum`,
// `ContourPrepare`, `ContourPoly`, `ContourEvaluate` (plus file-local
// `ContourCount`, `ClipToFrontHalf`).
//
// Faithfulness notes:
// - `_OPENMP` is NOT defined in the grade build: `ClosedFormNThreads = 1`,
//   `tid = 0` (the `#else` branches). Recorded in EVIDENCE.md.
// - `ClosedFormCount`, `ClosedFormNThreads`, `ClosedFormInteg` live in the
//   ViewFactors unit; declared extern here.
// - `TriangleValue`/`BiLinearValue` are header macros, expanded verbatim.
// - `CI_MAXV = CL_MAXV + 8` with `CL_MAXV = 24`.
// - malloc/calloc/free go through libc externs.

use std::os::raw::{c_double, c_int, c_long, c_void};

const CI_MAXV: usize = 32;
const CI_PLANAR_TOL: c_double = 1.0e-6;
const CI_NORMAL_TOL: c_double = 1.0e-6;
const CI_EPS: c_double = 1.0e-30;
const CF_STRIDE: usize = 8;

const GEOMETRY_TRIANGLE: c_int = 1;
const GEOMETRY_BILINEAR: c_int = 2;

#[repr(C)]
pub struct ContourTarget {
    pub nv: c_int,
    pub v: [[c_double; 3]; 4],
    pub n: [c_double; 3],
}

#[repr(C)]
pub struct Geometry {
    pub geometrytype: c_int,
    pub _bbox: [c_double; 6],
    pub entry: *mut c_void,
    pub flags: c_int,
    pub n: c_int,
    pub _rest: [c_double; 4],
    pub _link: *mut c_void,
    pub _left: *mut Geometry,
    pub _right: *mut Geometry,
}

#[link(name = "c")]
extern "C" {
    static mut ClosedFormCount: *mut c_long;
    static mut ClosedFormNThreads: c_int;
    static mut ClosedFormInteg: c_int;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn sqrt(x: c_double) -> c_double;
    fn atan2(y: c_double, x: c_double) -> c_double;
}

#[inline]
fn abs_c(x: c_double) -> c_double {
    if x > 0.0 { x } else { -x }
}

#[export_name = "ContourCountInit"]
pub unsafe extern "C" fn port_contour_count_init() {
    // _OPENMP off in the grade build.
    ClosedFormNThreads = 1;
    if ClosedFormNThreads < 1 {
        ClosedFormNThreads = 1;
    }
    if !ClosedFormCount.is_null() {
        free(ClosedFormCount as *mut c_void);
    }
    ClosedFormCount = calloc(
        2 * ClosedFormNThreads as usize * CF_STRIDE,
        std::mem::size_of::<c_long>(),
    ) as *mut c_long;
}

unsafe fn contour_count(used: c_int) {
    let tid: c_int = 0;
    if ClosedFormInteg == 0 || ClosedFormCount.is_null() {
        return;
    }
    if tid < 0 || tid >= ClosedFormNThreads {
        return;
    }
    let slot = if used != 0 { 0 } else { 1 };
    *ClosedFormCount.add((2 * tid as usize + slot) * CF_STRIDE) += 1;
}

#[export_name = "ContourCountSum"]
pub unsafe extern "C" fn port_contour_count_sum(hits: *mut c_long, miss: *mut c_long) {
    *hits = 0;
    *miss = 0;
    if ClosedFormCount.is_null() {
        return;
    }
    for i in 0..ClosedFormNThreads as usize {
        *hits += *ClosedFormCount.add((2 * i) * CF_STRIDE);
        *miss += *ClosedFormCount.add((2 * i + 1) * CF_STRIDE);
    }
}

unsafe fn clip_to_front_half(
    v: *const [c_double; 3], n: c_int,
    fx: c_double, fy: c_double, fz: c_double,
    nx: c_double, ny: c_double, nz: c_double,
    c: *mut [c_double; 3],
) -> c_int {
    let mut dc: c_double;
    let mut m: c_int = 0;
    let mut j = n - 1;
    dc = ((*v.add(j as usize))[0] - fx) * nx
        + ((*v.add(j as usize))[1] - fy) * ny
        + ((*v.add(j as usize))[2] - fz) * nz;
    for i in 0..n {
        let dn = ((*v.add(i as usize))[0] - fx) * nx
            + ((*v.add(i as usize))[1] - fy) * ny
            + ((*v.add(i as usize))[2] - fz) * nz;
        if dn >= 0.0 {
            if dc < 0.0 {
                if m as usize >= CI_MAXV {
                    return -1;
                }
                let t = dc / (dc - dn);
                for k in 0..3 {
                    (*c.add(m as usize))[k] =
                        (*v.add(j as usize))[k] + t * ((*v.add(i as usize))[k] - (*v.add(j as usize))[k]);
                }
                m += 1;
            }
            if m as usize >= CI_MAXV {
                return -1;
            }
            for k in 0..3 {
                (*c.add(m as usize))[k] = (*v.add(i as usize))[k];
            }
            m += 1;
        } else if dc >= 0.0 {
            if m as usize >= CI_MAXV {
                return -1;
            }
            let t = dc / (dc - dn);
            for k in 0..3 {
                (*c.add(m as usize))[k] =
                    (*v.add(j as usize))[k] + t * ((*v.add(i as usize))[k] - (*v.add(j as usize))[k]);
            }
            m += 1;
        }
        dc = dn;
        j = i;
    }
    m
}

#[inline]
unsafe fn triangle_value(u: c_double, v: c_double, nn: *const c_double) -> c_double {
    *nn.add(2) * v + *nn.add(1) * u + *nn.add(0)
}

#[inline]
unsafe fn bilinear_value(u: c_double, v: c_double, nn: *const c_double) -> c_double {
    (*nn.add(3) * u + *nn.add(2)) * v + *nn.add(1) * u + *nn.add(0)
}

#[export_name = "ContourPrepare"]
pub unsafe extern "C" fn port_contour_prepare(gb: *const Geometry, t: *mut ContourTarget) -> c_int {
    const QU: [c_double; 4] = [0.0, 1.0, 1.0, 0.0];
    const QV: [c_double; 4] = [0.0, 0.0, 1.0, 1.0];
    const TU: [c_double; 3] = [0.0, 1.0, 0.0];
    const TV: [c_double; 3] = [0.0, 0.0, 1.0];
    // PolyFactors[6][K]: entry points at row 0; row r strides K.
    let (bx, by, bz, nbx, nby, nbz, k): (*const c_double, *const c_double, *const c_double, *const c_double, *const c_double, *const c_double, usize);
    match (*gb).geometrytype {
        GEOMETRY_TRIANGLE => {
            let e = (*gb).entry as *const c_double;
            bx = e.add(0 * 3); by = e.add(1 * 3); bz = e.add(2 * 3);
            nbx = e.add(3 * 3); nby = e.add(4 * 3); nbz = e.add(5 * 3);
            k = 3;
            (*t).nv = 3;
            for i in 0..3 {
                (*t).v[i][0] = triangle_value(TU[i], TV[i], bx);
                (*t).v[i][1] = triangle_value(TU[i], TV[i], by);
                (*t).v[i][2] = triangle_value(TU[i], TV[i], bz);
            }
        }
        GEOMETRY_BILINEAR => {
            let e = (*gb).entry as *const c_double;
            bx = e.add(0 * 4); by = e.add(1 * 4); bz = e.add(2 * 4);
            nbx = e.add(3 * 4); nby = e.add(4 * 4); nbz = e.add(5 * 4);
            k = 4;
            (*t).nv = 4;
            for i in 0..4 {
                (*t).v[i][0] = bilinear_value(QU[i], QV[i], bx);
                (*t).v[i][1] = bilinear_value(QU[i], QV[i], by);
                (*t).v[i][2] = bilinear_value(QU[i], QV[i], bz);
            }
        }
        _ => return 0,
    }
    let _ = k;
    let (mut ntx, mut nty, mut ntz): (c_double, c_double, c_double);
    if (*gb).geometrytype == GEOMETRY_TRIANGLE {
        ntx = triangle_value(1.0 / 3.0, 1.0 / 3.0, nbx);
        nty = triangle_value(1.0 / 3.0, 1.0 / 3.0, nby);
        ntz = triangle_value(1.0 / 3.0, 1.0 / 3.0, nbz);
    } else {
        ntx = bilinear_value(0.5, 0.5, nbx);
        nty = bilinear_value(0.5, 0.5, nby);
        ntz = bilinear_value(0.5, 0.5, nbz);
    }
    let mut d = ntx * ntx + nty * nty + ntz * ntz;
    if d <= 0.0 {
        contour_count(0);
        return 0;
    }
    d = 1.0 / sqrt(d);
    ntx *= d;
    nty *= d;
    ntz *= d;
    let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
    for i in 0..(*t).nv as usize {
        let j = (i + 1) % (*t).nv as usize;
        gx += (*t).v[i][1] * (*t).v[j][2] - (*t).v[i][2] * (*t).v[j][1];
        gy += (*t).v[i][2] * (*t).v[j][0] - (*t).v[i][0] * (*t).v[j][2];
        gz += (*t).v[i][0] * (*t).v[j][1] - (*t).v[i][1] * (*t).v[j][0];
    }
    let mut a2 = sqrt(gx * gx + gy * gy + gz * gz);
    if a2 < CI_EPS {
        contour_count(0);
        return 0;
    }
    gx /= a2;
    gy /= a2;
    gz /= a2;
    let size = sqrt(a2);
    let _ = &mut a2;
    for i in 3..(*t).nv as usize {
        d = ((*t).v[i][0] - (*t).v[0][0]) * gx
            + ((*t).v[i][1] - (*t).v[0][1]) * gy
            + ((*t).v[i][2] - (*t).v[0][2]) * gz;
        if abs_c(d) > CI_PLANAR_TOL * size {
            contour_count(0);
            return 0;
        }
    }
    d = gx * ntx + gy * nty + gz * ntz;
    if abs_c(d) < 1.0 - CI_NORMAL_TOL {
        contour_count(0);
        return 0;
    }
    if d < 0.0 {
        gx = -gx;
        gy = -gy;
        gz = -gz;
    }
    (*t).n[0] = gx;
    (*t).n[1] = gy;
    (*t).n[2] = gz;
    contour_count(1);
    1
}

#[export_name = "ContourPoly"]
pub unsafe extern "C" fn port_contour_poly(
    v: *const [c_double; 3], nv: c_int, nn: *const c_double,
    fx: c_double, fy: c_double, fz: c_double,
    nfx: c_double, nfy: c_double, nfz: c_double,
    result: *mut c_double,
) -> c_int {
    let mut c = [[0.0; 3]; CI_MAXV];
    *result = 0.0;
    if nv < 3 {
        return 1;
    }
    // Back-facing target: every Gauss point is dropped, answer stays 0.
    if (fx - (*v.add(0))[0]) * *nn.add(0)
        + (fy - (*v.add(0))[1]) * *nn.add(1)
        + (fz - (*v.add(0))[2]) * *nn.add(2)
        <= 0.0
    {
        return 1;
    }
    // Front-side test over all vertices.
    let (p, n): (*const [c_double; 3], c_int);
    let mut s = 0.0;
    for i in 0..nv as usize {
        let d = ((*v.add(i))[0] - fx) * nfx + ((*v.add(i))[1] - fy) * nfy + ((*v.add(i))[2] - fz) * nfz;
        if d < s {
            s = d;
        }
    }
    if s < 0.0 {
        n = clip_to_front_half(v, nv, fx, fy, fz, nfx, nfy, nfz, c.as_mut_ptr());
        if n < 0 {
            return 0;
        }
        if n < 3 {
            return 1;
        }
        p = c.as_ptr();
    } else {
        p = v;
        n = nv;
    }
    let mut f = 0.0;
    for i in 0..n as usize {
        let j = (i + 1) % n as usize;
        let ax = (*p.add(i))[0] - fx;
        let ay = (*p.add(i))[1] - fy;
        let az = (*p.add(i))[2] - fz;
        let bx = (*p.add(j))[0] - fx;
        let by = (*p.add(j))[1] - fy;
        let bz = (*p.add(j))[2] - fz;
        let rx = ay * bz - az * by;
        let ry = az * bx - ax * bz;
        let rz = ax * by - ay * bx;
        let s = sqrt(rx * rx + ry * ry + rz * rz);
        if s < CI_EPS {
            continue;
        }
        let beta = atan2(s, ax * bx + ay * by + az * bz);
        f += beta * (nfx * rx + nfy * ry + nfz * rz) / s;
    }
    *result = 0.5 * abs_c(f);
    1
}

#[export_name = "ContourEvaluate"]
pub unsafe extern "C" fn port_contour_evaluate(
    t: *mut ContourTarget,
    fx: c_double, fy: c_double, fz: c_double,
    nfx: c_double, nfy: c_double, nfz: c_double,
    result: *mut c_double,
) -> c_int {
    port_contour_poly(
        (*t).v.as_ptr(), (*t).nv, (*t).n.as_ptr(),
        fx, fy, fz, nfx, nfy, nfz, result,
    )
}
