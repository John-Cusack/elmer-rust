// Port of `fem/src/view3d/ClipShadow.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/ClipShadow.c.
//
// True object symbols: `ClipVisible` (plus file-local `ClipHalf`,
// `PolyArea`, `ShadowWedge`, `ConvexDifference`, `PolyIntegNumeric` and
// the file-local `Wedge_t`).
//
// Faithfulness notes:
// - `ContourPrepare`, `ContourPoly` live in the contourinteg unit;
//   `RTElements`, `ClosedFormInteg`, and the `N/U/V/S_Integ3` triangle
//   rule tables live in other view3d units (ViewFactors); all extern.
// - Struct copies (`Out[nout++] = rest`, `A[i] = B[i]`, `rest = *P`)
//   are plain bitwise copies — `Clone, Copy` types.
// - Evaluation order and early-out structure kept verbatim.

use std::os::raw::{c_double, c_int, c_void};

const CL_EPS: c_double = 1.0e-30;
const CL_AREPS: c_double = 1.0e-14;
const CL_MAXV: usize = 24;
const CL_MAXW: usize = CL_MAXV + 1;
const CL_MAXPIECE: usize = 64;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct ClipPoly {
    pub n: c_int,
    pub v: [[c_double; 3]; CL_MAXV],
}

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

#[repr(C)]
#[derive(Clone, Copy)]
struct Wedge {
    n: c_int,
    nn: [[c_double; 3]; CL_MAXW],
    d: [c_double; CL_MAXW],
}

#[link(name = "c")]
extern "C" {
    static mut RTElements: *mut Geometry;
    static mut ClosedFormInteg: c_int;
    static N_Integ3: c_int;
    static U_Integ3: *const c_double;
    static V_Integ3: *const c_double;
    static S_Integ3: *const c_double;
    fn ContourPrepare(gb: *const Geometry, t: *mut ContourTarget) -> c_int;
    fn ContourPoly(
        v: *const [c_double; 3], nv: c_int, n: *const c_double,
        fx: c_double, fy: c_double, fz: c_double,
        nfx: c_double, nfy: c_double, nfz: c_double,
        result: *mut c_double,
    ) -> c_int;
    fn sqrt(x: c_double) -> c_double;
}

#[inline]
fn abs_c(x: c_double) -> c_double {
    if x > 0.0 { x } else { -x }
}

unsafe fn clip_half(
    inn: *const ClipPoly, nrm: *const c_double, d: c_double, side: c_double,
    out: *mut ClipPoly,
) -> c_int {
    let mut dc: c_double;
    let mut m: c_int = 0;
    if (*inn).n < 3 {
        (*out).n = 0;
        return 0;
    }
    let mut j = (*inn).n - 1;
    dc = side
        * (*nrm.add(0) * (*inn).v[j as usize][0] + *nrm.add(1) * (*inn).v[j as usize][1]
            + *nrm.add(2) * (*inn).v[j as usize][2]
            - d);
    for i in 0..(*inn).n {
        let dn = side
            * (*nrm.add(0) * (*inn).v[i as usize][0]
                + *nrm.add(1) * (*inn).v[i as usize][1]
                + *nrm.add(2) * (*inn).v[i as usize][2]
                - d);
        if dn >= 0.0 {
            if dc < 0.0 {
                if m as usize >= CL_MAXV {
                    (*out).n = 0;
                    return -1;
                }
                let t = dc / (dc - dn);
                for k in 0..3 {
                    (*out).v[m as usize][k] = (*inn).v[j as usize][k]
                        + t * ((*inn).v[i as usize][k] - (*inn).v[j as usize][k]);
                }
                m += 1;
            }
            if m as usize >= CL_MAXV {
                (*out).n = 0;
                return -1;
            }
            for k in 0..3 {
                (*out).v[m as usize][k] = (*inn).v[i as usize][k];
            }
            m += 1;
        } else if dc >= 0.0 {
            if m as usize >= CL_MAXV {
                (*out).n = 0;
                return -1;
            }
            let t = dc / (dc - dn);
            for k in 0..3 {
                (*out).v[m as usize][k] = (*inn).v[j as usize][k]
                    + t * ((*inn).v[i as usize][k] - (*inn).v[j as usize][k]);
            }
            m += 1;
        }
        dc = dn;
        j = i;
    }
    (*out).n = m;
    m
}

unsafe fn poly_area(p: *const ClipPoly) -> c_double {
    let mut g = [0.0; 3];
    if (*p).n < 3 {
        return 0.0;
    }
    for i in 0..(*p).n as usize {
        let j = (i + 1) % (*p).n as usize;
        g[0] += (*p).v[i][1] * (*p).v[j][2] - (*p).v[i][2] * (*p).v[j][1];
        g[1] += (*p).v[i][2] * (*p).v[j][0] - (*p).v[i][0] * (*p).v[j][2];
        g[2] += (*p).v[i][0] * (*p).v[j][1] - (*p).v[i][1] * (*p).v[j][0];
    }
    0.5 * sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2])
}

unsafe fn shadow_wedge(occ: *const ContourTarget, x: *const c_double, w: *mut Wedge) -> c_int {
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    let mut nn = [0.0; 3];
    let mut c = [0.0; 3];
    (*w).n = 0;
    let fx = (*occ).n[0] * (*x.add(0) - (*occ).v[0][0])
        + (*occ).n[1] * (*x.add(1) - (*occ).v[0][1])
        + (*occ).n[2] * (*x.add(2) - (*occ).v[0][2]);
    if abs_c(fx) < CL_EPS {
        return 0;
    }
    let mut s = if fx > 0.0 { -1.0 } else { 1.0 };
    for k in 0..3 {
        (*w).nn[0][k] = s * (*occ).n[k];
    }
    (*w).d[0] = s
        * ((*occ).n[0] * (*occ).v[0][0] + (*occ).n[1] * (*occ).v[0][1] + (*occ).n[2] * (*occ).v[0][2]);
    (*w).n = 1;
    for k in 0..3 {
        c[k] = 0.0;
    }
    for i in 0..(*occ).nv as usize {
        for k in 0..3 {
            c[k] += (*occ).v[i][k] / (*occ).nv as c_double;
        }
    }
    for i in 0..(*occ).nv as usize {
        let j = (i + 1) % (*occ).nv as usize;
        for k in 0..3 {
            a[k] = (*occ).v[i][k] - *x.add(k);
            b[k] = (*occ).v[j][k] - *x.add(k);
        }
        nn[0] = a[1] * b[2] - a[2] * b[1];
        nn[1] = a[2] * b[0] - a[0] * b[2];
        nn[2] = a[0] * b[1] - a[1] * b[0];
        let l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
        if l < CL_EPS {
            return 0;
        }
        for k in 0..3 {
            nn[k] /= l;
        }
        s = nn[0] * (c[0] - *x.add(0)) + nn[1] * (c[1] - *x.add(1)) + nn[2] * (c[2] - *x.add(2));
        if abs_c(s) < CL_EPS {
            return 0;
        }
        if s < 0.0 {
            for k in 0..3 {
                nn[k] = -nn[k];
            }
        }
        if (*w).n as usize >= CL_MAXW {
            return 0;
        }
        let wn = (*w).n as usize;
        for k in 0..3 {
            (*w).nn[wn][k] = nn[k];
        }
        (*w).d[wn] = nn[0] * *x.add(0) + nn[1] * *x.add(1) + nn[2] * *x.add(2);
        (*w).n += 1;
    }
    (*w).n
}

unsafe fn convex_difference(
    p: *const ClipPoly, w: *const Wedge, out: *mut ClipPoly, mut nout: c_int, maxout: c_int,
) -> c_int {
    let mut rest = *p;
    let mut keep = std::mem::zeroed::<ClipPoly>();
    let mut drop = std::mem::zeroed::<ClipPoly>();
    for i in 0..(*w).n as usize {
        let mut dmin = 1.0e30;
        let mut dmax = -1.0e30;
        let mut scale = 0.0;
        for j in 0..rest.n as usize {
            let mut d = (*w).nn[i][0] * rest.v[j][0]
                + (*w).nn[i][1] * rest.v[j][1]
                + (*w).nn[i][2] * rest.v[j][2];
            if abs_c(d) > scale {
                scale = abs_c(d);
            }
            d -= (*w).d[i];
            if d < dmin {
                dmin = d;
            }
            if d > dmax {
                dmax = d;
            }
        }
        let tol = 1.0e-12 * (scale + abs_c((*w).d[i])) + CL_EPS;
        if dmax - dmin <= tol {
            continue;
        }
        if dmin >= -tol {
            continue;
        }
        if dmax <= tol {
            if nout >= maxout {
                return -1;
            }
            *out.add(nout as usize) = rest;
            nout += 1;
            return nout;
        }
        if clip_half(&rest, (*w).nn[i].as_ptr(), (*w).d[i], -1.0, &mut drop) < 0 {
            return -1;
        }
        if drop.n >= 3 && poly_area(&drop) > CL_AREPS {
            if nout >= maxout {
                return -1;
            }
            *out.add(nout as usize) = drop;
            nout += 1;
        }
        if clip_half(&rest, (*w).nn[i].as_ptr(), (*w).d[i], 1.0, &mut keep) < 0 {
            return -1;
        }
        rest = keep;
        if rest.n < 3 {
            break;
        }
    }
    nout
}

unsafe fn poly_integ_numeric(
    v: *const [c_double; 3], nv: c_int, nn: *const c_double,
    fx: c_double, fy: c_double, fz: c_double,
    nfx: c_double, nfy: c_double, nfz: c_double,
    result: *mut c_double,
) -> c_int {
    let mut f = 0.0;
    *result = 0.0;
    if nv < 3 {
        return 1;
    }
    if (fx - (*v.add(0))[0]) * *nn.add(0)
        + (fy - (*v.add(0))[1]) * *nn.add(1)
        + (fz - (*v.add(0))[2]) * *nn.add(2)
        <= 0.0
    {
        return 1;
    }
    for i in 1..nv as usize - 1 {
        let mut e1 = [0.0; 3];
        let mut e2 = [0.0; 3];
        let mut cr = [0.0; 3];
        for k in 0..3 {
            e1[k] = (*v.add(i))[k] - (*v.add(0))[k];
            e2[k] = (*v.add(i + 1))[k] - (*v.add(0))[k];
        }
        cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
        cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
        cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
        let a2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
        if a2 <= 0.0 {
            continue;
        }
        for j in 0..N_Integ3 as usize {
            let mut d = [0.0; 3];
            let u = *U_Integ3.add(j);
            let w = *V_Integ3.add(j);
            for k in 0..3 {
                d[k] = (*v.add(0))[k] + u * e1[k] + w * e2[k];
            }
            d[0] -= fx;
            d[1] -= fy;
            d[2] -= fz;
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if r2 <= 0.0 {
                continue;
            }
            let r = sqrt(r2);
            let cosa = (d[0] * nfx + d[1] * nfy + d[2] * nfz) / r;
            if cosa < 1.0e-8 {
                continue;
            }
            let cosb = -(d[0] * *nn.add(0) + d[1] * *nn.add(1) + d[2] * *nn.add(2)) / r;
            if cosb < 1.0e-8 {
                continue;
            }
            f += a2 * cosa * cosb * *S_Integ3.add(j) / r2;
        }
    }
    *result = f;
    1
}

#[export_name = "ClipVisible"]
pub unsafe extern "C" fn port_clip_visible(
    tgt: *const ContourTarget, cand: *const c_int, nc: c_int,
    x: *const c_double, nf: *const c_double, result: *mut c_double,
) -> c_int {
    let mut aa = [ClipPoly { n: 0, v: [[0.0; 3]; CL_MAXV] }; CL_MAXPIECE];
    let mut bb = [ClipPoly { n: 0, v: [[0.0; 3]; CL_MAXV] }; CL_MAXPIECE];
    let mut occ = std::mem::zeroed::<ContourTarget>();
    let mut w = std::mem::zeroed::<Wedge>();
    let mut f: c_double;
    let mut ff: c_double = 0.0;
    *result = 0.0;
    aa[0].n = (*tgt).nv;
    for i in 0..(*tgt).nv as usize {
        for k in 0..3 {
            aa[0].v[i][k] = (*tgt).v[i][k];
        }
    }
    let mut na: c_int = 1;
    for c in 0..nc {
        let mut dmin = 1.0e30;
        let mut dmax = -1.0e30;
        let mut scale = 0.0;
        if ContourPrepare(RTElements.add(*cand.add(c as usize) as usize), &mut occ) == 0 {
            return 0;
        }
        for i in 0..occ.nv as usize {
            let mut d = 0.0;
            for k in 0..3 {
                d += (*tgt).n[k] * (occ.v[i][k] - (*tgt).v[0][k]);
            }
            let m = abs_c(occ.v[i][0]) + abs_c(occ.v[i][1]) + abs_c(occ.v[i][2]);
            if m > scale {
                scale = m;
            }
            if d < dmin {
                dmin = d;
            }
            if d > dmax {
                dmax = d;
            }
        }
        let tol = 1.0e-10 * scale + CL_EPS;
        if dmax - dmin <= tol && abs_c(dmin) <= tol {
            continue;
        }
        if shadow_wedge(&occ, x, &mut w) == 0 {
            continue;
        }
        let mut nb: c_int = 0;
        for i in 0..na as usize {
            nb = convex_difference(&aa[i], &w, bb.as_mut_ptr(), nb, CL_MAXPIECE as c_int);
            if nb < 0 {
                return 0;
            }
        }
        for i in 0..nb as usize {
            aa[i] = bb[i];
        }
        na = nb;
        if na == 0 {
            return 1;
        }
    }
    f = 0.0;
    for i in 0..na as usize {
        if ClosedFormInteg != 0 {
            if ContourPoly(
                aa[i].v.as_ptr(), aa[i].n, (*tgt).n.as_ptr(),
                *x.add(0), *x.add(1), *x.add(2),
                *nf.add(0), *nf.add(1), *nf.add(2),
                &mut ff,
            ) == 0
            {
                return 0;
            }
        } else if poly_integ_numeric(
            aa[i].v.as_ptr(), aa[i].n, (*tgt).n.as_ptr(),
            *x.add(0), *x.add(1), *x.add(2),
            *nf.add(0), *nf.add(1), *nf.add(2),
            &mut ff,
        ) == 0
        {
            return 0;
        }
        f += ff;
    }
    *result = f;
    1
}
