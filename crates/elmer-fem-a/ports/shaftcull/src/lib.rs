#![allow(static_mut_refs)]
//! `static mut` mirrors C file-static mutable state.

// Port of `fem/src/view3d/ShaftCull.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/ShaftCull.c.
//
// True object symbols: `ShaftCountInit`, `ShaftCountAdd`, `ShaftCountSum`,
// `ShaftInitBoxes`, `ShaftInit`, `ShaftCandidates` (plus file-local
// `BoxBehind`, `BoxDisjoint`, `BoxInShaft`, `TryHullPlane`, `CannotBlock`,
// `Walk`, and a header-static copy of `FunctionValue` — all inlined).
//
// Faithfulness notes:
// - `_OPENMP` is NOT defined in the grade build: `ShaftNThreads = 1`,
//   `tid = 0` (the `#else` branches). Recorded in EVIDENCE.md.
// - `ShaftCount`, `RTElementBBox`, `RTElementNof`, `VolumeBounds` and
//   `ContourPrepare` live in other units; declared extern here.
// - `MIN`/`MAX`/`ABS` and the `BiLinearValue`/`TriangleValue`/…
//   polynomials are expanded exactly as in ViewFactors.h.
// - malloc/calloc/free go through libc externs (pointer-identity and
//   cross-call lifetimes match the original).

use std::os::raw::{c_double, c_int, c_long, c_void};

const SC_EPS: c_double = 1.0e-12;
const SH_MAXPLANE: usize = 24;
const SC_NBUCKET: usize = 9;
const CF_STRIDE: usize = 8;

const GEOMETRY_LINE: c_int = 0;
const GEOMETRY_TRIANGLE: c_int = 1;
const GEOMETRY_BILINEAR: c_int = 2;
const GEOMETRY_BICUBIC: c_int = 3;
const GEOMETRY_BIQUADRATIC: c_int = 4;

#[repr(C)]
pub struct BBox {
    pub xmin: c_double, pub ymin: c_double, pub zmin: c_double,
    pub xmax: c_double, pub ymax: c_double, pub zmax: c_double,
}

#[repr(C)]
pub struct ContourTarget {
    pub nv: c_int,
    pub v: [[c_double; 3]; 4],
    pub n: [c_double; 3],
}

#[repr(C)]
pub struct Shaft {
    pub bbox: BBox,
    pub n: [[c_double; 3]; SH_MAXPLANE],
    pub d: [c_double; SH_MAXPLANE],
    pub np: c_int,
    pub tol: c_double,
    pub a: ContourTarget,
    pub b: ContourTarget,
}

#[repr(C)]
pub struct VolumeBounds {
    pub left: *mut VolumeBounds,
    pub right: *mut VolumeBounds,
    pub bbox: BBox,
    pub n: c_int,
    pub elements: *mut c_int,
}

// Geometry_t mirror (only the fields this unit touches through pointers;
// full layout matches ViewFactors.h by construction — see EVIDENCE.md).
#[repr(C)]
pub struct Geometry {
    pub geometrytype: c_int,
    pub bbox: BBox,
    pub entry: *mut c_void, // union of factor pointers; indexed per type
    pub flags: c_int,
    pub n: c_int,
    pub area: c_double, pub b_: c_double, pub e: c_double, pub m: c_double,
    pub link: *mut c_void,
    pub left: *mut Geometry, pub right: *mut Geometry,
}

#[link(name = "c")]
extern "C" {
    static mut ShaftCount: *mut c_long;
    static mut ShaftNThreads: c_int;
    static mut RTElementBBox: *mut BBox;
    static mut RTElementNof: c_int;
    static mut VolumeBounds: *mut VolumeBounds;
    static mut RTElements: *mut Geometry;
    fn ContourPrepare(gb: *const Geometry, t: *mut ContourTarget) -> c_int;
    fn malloc(n: usize) -> *mut c_void;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn sqrt(x: c_double) -> c_double;
}

#[inline]
fn min_c(a: c_double, b: c_double) -> c_double {
    if a > b { b } else { a }
}
#[inline]
fn max_c(a: c_double, b: c_double) -> c_double {
    if a > b { a } else { b }
}

unsafe fn function_value(geom: *const Geometry, u: c_double, v: c_double, n: c_int) -> c_double {
    let entry = (*geom).entry as *const c_double;
    // PolyFactors[6][K]: row N holds the factor vector.
    match (*geom).geometrytype {
        GEOMETRY_LINE => {
            let p = entry.add(n as usize * 2);
            *p.add(1) * u + *p.add(0)
        }
        GEOMETRY_TRIANGLE => {
            let p = entry.add(n as usize * 3);
            *p.add(2) * v + *p.add(1) * u + *p.add(0)
        }
        GEOMETRY_BILINEAR => {
            let p = entry.add(n as usize * 4);
            (*p.add(3) * u + *p.add(2)) * v + *p.add(1) * u + *p.add(0)
        }
        GEOMETRY_BICUBIC => {
            let p = entry.add(n as usize * 16);
            let get = |i: usize| *p.add(i);
            (((((get(15) * u + get(14)) * u + get(13)) * u + get(12)) * v
                + ((get(11) * u + get(10)) * u + get(9)) * u
                + get(8))
                * v
                + ((get(7) * u + get(6)) * u + get(5)) * u
                + get(4))
                * v
                + ((get(3) * u + get(2)) * u + get(1)) * u
                + get(0)
        }
        GEOMETRY_BIQUADRATIC => {
            let p = entry.add(n as usize * 9);
            let get = |i: usize| *p.add(i);
            let a = ((get(8) * u + get(7)) * u + get(6)) * v
                + (get(5) * u + get(4)) * u
                + get(3);
            a * v + (get(2) * u + get(1)) * u + get(0)
        }
        _ => 0.0,
    }
}

#[export_name = "ShaftCountInit"]
pub unsafe extern "C" fn port_shaft_count_init() {
    // _OPENMP off in the grade build: single thread.
    let n: c_int = 1;
    ShaftNThreads = n;
    if !ShaftCount.is_null() {
        free(ShaftCount as *mut c_void);
    }
    ShaftCount = calloc(
        (n as usize) * SC_NBUCKET * CF_STRIDE,
        std::mem::size_of::<c_long>(),
    ) as *mut c_long;
}

#[export_name = "ShaftCountAdd"]
pub unsafe extern "C" fn port_shaft_count_add(n: c_int) {
    let tid: c_int = 0;
    if ShaftCount.is_null() {
        return;
    }
    if tid < 0 || tid >= ShaftNThreads {
        return;
    }
    let b: usize = if n == -3 {
        8
    } else if n == -2 {
        7
    } else if n < 0 {
        6
    } else if n == 0 {
        0
    } else if n <= 2 {
        1
    } else if n <= 4 {
        2
    } else if n <= 8 {
        3
    } else if n <= 16 {
        4
    } else {
        5
    };
    *ShaftCount.add((tid as usize * SC_NBUCKET + b) * CF_STRIDE) += 1;
}

#[export_name = "ShaftCountSum"]
pub unsafe extern "C" fn port_shaft_count_sum(buckets: *mut c_long) {
    for b in 0..SC_NBUCKET {
        *buckets.add(b) = 0;
    }
    if ShaftCount.is_null() {
        return;
    }
    for i in 0..ShaftNThreads as usize {
        for b in 0..SC_NBUCKET {
            *buckets.add(b) += *ShaftCount.add((i * SC_NBUCKET + b) * CF_STRIDE);
        }
    }
}

#[export_name = "ShaftInitBoxes"]
pub unsafe extern "C" fn port_shaft_init_boxes(n: c_int, el: *mut Geometry) {
    const U: [c_double; 4] = [0.0, 1.0, 0.0, 1.0];
    const V: [c_double; 4] = [0.0, 0.0, 1.0, 1.0];
    if !RTElementBBox.is_null() {
        free(RTElementBBox as *mut c_void);
    }
    RTElementBBox =
        malloc(n as usize * std::mem::size_of::<BBox>()) as *mut BBox;
    RTElementNof = n;
    for i in 0..n as usize {
        let b = &mut *RTElementBBox.add(i);
        let e = &*el.add(i);
        b.xmin = 1.0e20; b.ymin = 1.0e20; b.zmin = 1.0e20;
        b.xmax = -1.0e20; b.ymax = -1.0e20; b.zmax = -1.0e20;
        // Circle layout: RMin,RMax,Axis(3),CenterPoint(3),IdentMatrix,RotationMatrix
        if e.geometrytype == 7 {
            // GEOMETRY_CIRCLE
            let rp = e.entry as *const c_double;
            let rmax = *rp.add(1);
            let cx = *rp.add(2 + 3);
            let cy = *rp.add(2 + 4);
            let cz = *rp.add(2 + 5);
            b.xmin = cx - rmax; b.xmax = cx + rmax;
            b.ymin = cy - rmax; b.ymax = cy + rmax;
            b.zmin = cz; b.zmax = cz;
            continue;
        }
        let nc: usize = match e.geometrytype {
            0 => 2,
            1 => 3,
            _ => 4,
        };
        for j in 0..nc {
            let x = function_value(e, U[j], V[j], 0);
            let y = function_value(e, U[j], V[j], 1);
            let z = function_value(e, U[j], V[j], 2);
            b.xmin = min_c(x, b.xmin); b.xmax = max_c(x, b.xmax);
            b.ymin = min_c(y, b.ymin); b.ymax = max_c(y, b.ymax);
            b.zmin = min_c(z, b.zmin); b.zmax = max_c(z, b.zmax);
        }
    }
}

unsafe fn box_behind(b: *const BBox, nv: *const c_double, d: c_double) -> bool {
    let x = if *nv.add(0) > 0.0 { (*b).xmax } else { (*b).xmin };
    let y = if *nv.add(1) > 0.0 { (*b).ymax } else { (*b).ymin };
    let z = if *nv.add(2) > 0.0 { (*b).zmax } else { (*b).zmin };
    *nv.add(0) * x + *nv.add(1) * y + *nv.add(2) * z < d
}

unsafe fn box_disjoint(a: *const BBox, b: *const BBox) -> bool {
    (*a).xmax < (*b).xmin
        || (*a).xmin > (*b).xmax
        || (*a).ymax < (*b).ymin
        || (*a).ymin > (*b).ymax
        || (*a).zmax < (*b).zmin
        || (*a).zmin > (*b).zmax
}

unsafe fn box_in_shaft(s: *const Shaft, b: *const BBox) -> bool {
    if box_disjoint(b, &(*s).bbox) {
        return false;
    }
    for i in 0..(*s).np as usize {
        if box_behind(b, (*s).n[i].as_ptr(), (*s).d[i]) {
            return false;
        }
    }
    true
}

unsafe fn try_hull_plane(
    s: *mut Shaft, a: *const ContourTarget, b: *const ContourTarget,
    p: *const c_double, q: *const c_double, r: *const c_double,
) {
    let mut nv = [0.0; 3];
    let mut u = [0.0; 3];
    let mut vv = [0.0; 3];
    if (*s).np >= SH_MAXPLANE as c_int {
        return;
    }
    for k in 0..3 {
        u[k] = *q.add(k) - *p.add(k);
        vv[k] = *r.add(k) - *p.add(k);
    }
    nv[0] = u[1] * vv[2] - u[2] * vv[1];
    nv[1] = u[2] * vv[0] - u[0] * vv[2];
    nv[2] = u[0] * vv[1] - u[1] * vv[0];
    let l = sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
    if l < SC_EPS {
        return;
    }
    for k in 0..3 {
        nv[k] /= l;
    }
    let mut d = nv[0] * *p.add(0) + nv[1] * *p.add(1) + nv[2] * *p.add(2);
    let mut pos = 0;
    let mut neg = 0;
    let na = (*a).nv as usize;
    let nb = (*b).nv as usize;
    for i in 0..na + nb {
        let x: *const c_double = if i < na {
            (*a).v[i].as_ptr()
        } else {
            (*b).v[i - na].as_ptr()
        };
        let t = nv[0] * *x.add(0) + nv[1] * *x.add(1) + nv[2] * *x.add(2) - d;
        if t > (*s).tol {
            pos += 1;
        }
        if t < -(*s).tol {
            neg += 1;
        }
    }
    if pos != 0 && neg != 0 {
        return;
    }
    if pos == 0 && neg == 0 {
        return;
    }
    if neg != 0 {
        for k in 0..3 {
            nv[k] = -nv[k];
        }
        d = -d;
    }
    let np = (*s).np as usize;
    for k in 0..3 {
        (*s).n[np][k] = nv[k];
    }
    (*s).d[np] = d - (*s).tol;
    (*s).np += 1;
}

#[export_name = "ShaftInit"]
pub unsafe extern "C" fn port_shaft_init(
    s: *mut Shaft, a: *const ContourTarget, b: *const ContourTarget,
) -> c_int {
    let mut ca = [0.0; 3];
    let mut cb = [0.0; 3];
    (*s).bbox.xmin = 1.0e20; (*s).bbox.ymin = 1.0e20; (*s).bbox.zmin = 1.0e20;
    (*s).bbox.xmax = -1.0e20; (*s).bbox.ymax = -1.0e20; (*s).bbox.zmax = -1.0e20;
    for i in 0..(*a).nv as usize {
        (*s).bbox.xmin = min_c((*a).v[i][0], (*s).bbox.xmin);
        (*s).bbox.xmax = max_c((*a).v[i][0], (*s).bbox.xmax);
        (*s).bbox.ymin = min_c((*a).v[i][1], (*s).bbox.ymin);
        (*s).bbox.ymax = max_c((*a).v[i][1], (*s).bbox.ymax);
        (*s).bbox.zmin = min_c((*a).v[i][2], (*s).bbox.zmin);
        (*s).bbox.zmax = max_c((*a).v[i][2], (*s).bbox.zmax);
        for k in 0..3 {
            ca[k] += (*a).v[i][k] / (*a).nv as c_double;
        }
    }
    for i in 0..(*b).nv as usize {
        (*s).bbox.xmin = min_c((*b).v[i][0], (*s).bbox.xmin);
        (*s).bbox.xmax = max_c((*b).v[i][0], (*s).bbox.xmax);
        (*s).bbox.ymin = min_c((*b).v[i][1], (*s).bbox.ymin);
        (*s).bbox.ymax = max_c((*b).v[i][1], (*s).bbox.ymax);
        (*s).bbox.zmin = min_c((*b).v[i][2], (*s).bbox.zmin);
        (*s).bbox.zmax = max_c((*b).v[i][2], (*s).bbox.zmax);
        for k in 0..3 {
            cb[k] += (*b).v[i][k] / (*b).nv as c_double;
        }
    }
    let tol = max_c(
        (*s).bbox.xmax - (*s).bbox.xmin,
        max_c(
            (*s).bbox.ymax - (*s).bbox.ymin,
            (*s).bbox.zmax - (*s).bbox.zmin,
        ),
    );
    (*s).tol = 1.0e-8 * tol + SC_EPS;
    std::ptr::copy_nonoverlapping(a, &mut (*s).a, 1);
    std::ptr::copy_nonoverlapping(b, &mut (*s).b, 1);
    (*s).np = 0;
    try_hull_plane(s, a, b, (*a).v[0].as_ptr(), (*a).v[1].as_ptr(), (*a).v[2].as_ptr());
    try_hull_plane(s, a, b, (*b).v[0].as_ptr(), (*b).v[1].as_ptr(), (*b).v[2].as_ptr());
    for i in 0..(*a).nv {
        for k in 0..(*b).nv {
            try_hull_plane(
                s, a, b,
                (*a).v[i as usize].as_ptr(),
                (*a).v[((i + 1) % (*a).nv) as usize].as_ptr(),
                (*b).v[k as usize].as_ptr(),
            );
        }
    }
    for i in 0..(*b).nv {
        for k in 0..(*a).nv {
            try_hull_plane(
                s, a, b,
                (*b).v[i as usize].as_ptr(),
                (*b).v[((i + 1) % (*b).nv) as usize].as_ptr(),
                (*a).v[k as usize].as_ptr(),
            );
        }
    }
    1
}

unsafe fn cannot_block(s: *const Shaft, j: c_int) -> bool {
    let mut t = std::mem::zeroed::<ContourTarget>();
    if ContourPrepare(rt_elements_ptr(j), &mut t) == 0 {
        return false;
    }
    let mut pos = 0;
    let mut neg = 0;
    let na = (*s).a.nv as usize;
    let nb = (*s).b.nv as usize;
    for i in 0..na + nb {
        let x: *const c_double = if i < na {
            (*s).a.v[i].as_ptr()
        } else {
            (*s).b.v[i - na].as_ptr()
        };
        let mut tt = 0.0;
        for k in 0..3 {
            tt += t.n[k] * (*x.add(k) - t.v[0][k]);
        }
        if tt > (*s).tol {
            pos += 1;
        }
        if tt < -(*s).tol {
            neg += 1;
        }
    }
    !(pos != 0 && neg != 0)
}

unsafe fn rt_elements_ptr(j: c_int) -> *const Geometry {
    RTElements.add(j as usize)
}

unsafe fn walk(
    s: *const Shaft, volume: *const VolumeBounds, list: *mut c_int,
    n: c_int, maxn: c_int,
) -> c_int {
    if volume.is_null() {
        return n;
    }
    if !box_in_shaft(s, &(*volume).bbox) {
        return n;
    }
    if !(*volume).left.is_null() {
        let n = walk(s, (*volume).left, list, n, maxn);
        if n < 0 {
            return n;
        }
        return walk(s, (*volume).right, list, n, maxn);
    }
    let mut n = n;
    for i in 0..(*volume).n {
        // NOTE: the original reuses `j` for the dedup scan (clobbering
        // the element index) — kept verbatim via shadowing discipline.
        let jj = *(*volume).elements.add(i as usize);
        if jj < 0 || jj >= RTElementNof {
            continue;
        }
        if !box_in_shaft(s, &*RTElementBBox.add(jj as usize)) {
            continue;
        }
        if cannot_block(s, jj) {
            continue;
        }
        let mut j: c_int = 0;
        while j < n {
            if *list.add(j as usize) == *(*volume).elements.add(i as usize) {
                break;
            }
            j += 1;
        }
        if j < n {
            continue;
        }
        if n >= maxn {
            return -1;
        }
        *list.add(n as usize) = *(*volume).elements.add(i as usize);
        n += 1;
    }
    n
}

#[export_name = "ShaftCandidates"]
pub unsafe extern "C" fn port_shaft_candidates(
    s: *const Shaft, list: *mut c_int, maxn: c_int,
) -> c_int {
    if VolumeBounds.is_null() || RTElementBBox.is_null() {
        return 0;
    }
    walk(s, VolumeBounds, list, 0, maxn)
}
