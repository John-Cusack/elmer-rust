// Port of `fem/src/view3d/LineUtil.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/LineUtil.c.
//
// True object symbols: `LinearSubdivideHalfU`, `LinearSubdivide`,
// `LinearEofL`, `LinearEofA`, `LinearLength`, `LinearArea`,
// `LinearIntegrateDiffToArea`, `LinearViewFactor`,
// `LinearComputeViewFactors`, `LinearComputeRadiatorFactors`.
//
// Faithfulness notes:
// - `LinearValue`/`LinearPartialU` header macros expanded verbatim.
// - `LinearViewFactor` divides the not-yet-assigned `NTX`/`NTY` when the
//   normal needs renormalizing; the quotient is immediately overwritten,
//   so the port performs the same dead division from 0.0 (x86 FP division
//   never traps; the value is unobservable). Recorded in EVIDENCE.md.
// - `F` in `LinearComputeViewFactors` is uninitialized when `GA == GB`
//   reaches the `Hit == Nrays` path with `#if 1`; the port uses 0.0 —
//   the upstream value is indeterminate there (documented, unprobed).
// - `goto subdivide` is structured with a flag, same statement order.
// - `ShapeFunctionMatrix2`, the `N/U/S_Integ1d` tables, `AreaCompute`/
//   `IntegrateDiffToArea`/`Subdivide` tables, `AreaEPS`/`FactorEPS`,
//   `Nrays`, `RayHitGeometry`, `vrand`, `FunctionValue` live in other
//   units; declared extern. malloc/calloc/free via libc.

use std::os::raw::{c_double, c_int, c_void};

const GEOMETRY_TRIANGLE: c_int = 1;
const GEOMETRY_LINE: c_int = 0;
const GEOMETRY_FLAG_LEAF: c_int = 2;

#[repr(C)]
pub struct BBox {
    pub xmin: c_double, pub ymin: c_double, pub zmin: c_double,
    pub xmax: c_double, pub ymax: c_double, pub zmax: c_double,
}

#[repr(C)]
pub struct Geometry {
    pub geometrytype: c_int,
    pub bbox: BBox,
    pub entry: *mut c_void,
    pub flags: c_int,
    pub n: c_int,
    pub area: c_double, pub b_: c_double, pub e: c_double, pub m: c_double,
    pub link: *mut GeometryList,
    pub left: *mut Geometry, pub right: *mut Geometry,
}

#[repr(C)]
pub struct GeometryList {
    pub next: *mut GeometryList,
    pub viewfactor: c_double,
    pub entry: *mut Geometry,
}

#[repr(C)]
pub struct Cylinder {
    pub radius: c_double,
    pub length: c_double,
    pub axis: [c_double; 3],
    pub centerpoint: [c_double; 3],
    pub identmatrix: c_int,
    pub rotationmatrix: [[c_double; 3]; 3],
}

#[link(name = "c")]
extern "C" {
    static mut ShapeFunctionMatrix2: [[c_double; 2]; 2];
    static N_Integ1d: c_int;
    static U_Integ1d: *const c_double;
    static S_Integ1d: *const c_double;
    static mut AreaEPS: c_double;
    static mut FactorEPS: c_double;
    static mut Nrays: c_int;
    static mut AreaCompute: [*mut c_void; 9];
    static mut IntegrateDiffToArea: [*mut c_void; 9];
    static mut Subdivide: [*mut c_void; 9];
    fn RayHitGeometry(
        fx: c_double, fy: c_double, fz: c_double,
        dx: c_double, dy: c_double, dz: c_double,
    ) -> c_int;
    fn vrand() -> c_double;
    fn FunctionValue(geom: *const Geometry, u: c_double, v: c_double, n: c_int) -> c_double;
    fn malloc(n: usize) -> *mut c_void;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn sqrt(x: c_double) -> c_double;
    fn acos(x: c_double) -> c_double;
}

#[inline]
fn abs_c(x: c_double) -> c_double {
    if x > 0.0 { x } else { -x }
}

#[inline]
unsafe fn linear_value(u: c_double, nn: *const c_double) -> c_double {
    *nn.add(1) * u + *nn.add(0)
}

#[inline]
unsafe fn linear_partial_u(_u: c_double, nn: *const c_double) -> c_double {
    *nn.add(1)
}

type IntegrateFn = unsafe extern "C" fn(
    *const Geometry, *const Cylinder,
    c_double, c_double, c_double, c_double, c_double, c_double,
) -> c_double;
type AreaFn = unsafe extern "C" fn(*const Geometry) -> c_double;
type SubdivideFn = unsafe extern "C" fn(*mut Geometry, c_int, c_int);

#[export_name = "LinearSubdivideHalfU"]
pub unsafe extern "C" fn port_linear_subdivide_half_u(
    i: *const c_double, l: *mut c_double, r: *mut c_double,
) {
    let mut x = [0.0; 2];
    x[0] = linear_value(0.0, i);
    x[1] = linear_value(0.5, i);
    for ii in 0..2 {
        for jj in 0..2 {
            *l.add(ii) += ShapeFunctionMatrix2[jj][ii] * x[jj];
        }
    }
    x[0] = linear_value(0.5, i);
    x[1] = linear_value(1.0, i);
    for ii in 0..2 {
        for jj in 0..2 {
            *r.add(ii) += ShapeFunctionMatrix2[jj][ii] * x[jj];
        }
    }
}

#[export_name = "LinearSubdivide"]
pub unsafe extern "C" fn port_linear_subdivide(
    geometry: *mut Geometry, _sublev: c_int, _where: c_int,
) {
    (*geometry).left = calloc(1, std::mem::size_of::<Geometry>()) as *mut Geometry;
    (*(*geometry).left).entry =
        calloc(1, 6 * 2 * std::mem::size_of::<c_double>()) as *mut c_void;
    (*geometry).right = calloc(1, std::mem::size_of::<Geometry>()) as *mut Geometry;
    (*(*geometry).right).entry =
        calloc(1, 6 * 2 * std::mem::size_of::<c_double>()) as *mut c_void;
    (*(*geometry).left).geometrytype = GEOMETRY_LINE;
    (*(*geometry).right).geometrytype = GEOMETRY_LINE;
    // Linear_t PolyFactors[6][2]: row j, entries 0..1.
    let src = (*geometry).entry as *const c_double;
    for j in 0..6 {
        port_linear_subdivide_half_u(
            src.add(j * 2),
            ((*(*geometry).left).entry as *mut c_double).add(j * 2),
            ((*(*geometry).right).entry as *mut c_double).add(j * 2),
        );
    }
}

#[export_name = "LinearEofL"]
pub unsafe extern "C" fn port_linear_eofl(
    u: c_double, _v: c_double,
    x: *const c_double, y: *const c_double, _z: *const c_double, _unotv: c_int,
) -> c_double {
    let dxdu = linear_partial_u(u, x);
    let dydu = linear_partial_u(u, y);
    sqrt(dxdu * dxdu + dydu * dydu)
}

#[export_name = "LinearEofA"]
pub unsafe extern "C" fn port_linear_eofa(
    u: c_double, _v: c_double, x: *const c_double, y: *const c_double, _z: *const c_double,
) -> c_double {
    let dxdu = linear_partial_u(u, x);
    let dydu = linear_partial_u(u, y);
    let deta = dxdu * dxdu + dydu * dydu;
    sqrt(deta)
}

#[export_name = "LinearLength"]
pub unsafe extern "C" fn port_linear_length(geometry: *const Geometry, unotv: c_int) -> c_double {
    let base = (*geometry).entry as *const c_double;
    // PolyFactors[6][2]: X=row0, Y=row1, Z=row2.
    let x = base.add(0 * 2);
    let y = base.add(1 * 2);
    let z = base.add(2 * 2);
    let mut length = 0.0;
    for j in 0..2 {
        for i in 0..N_Integ1d as usize {
            let eofl = if unotv != 0 {
                port_linear_eofl(*U_Integ1d.add(i), j as c_double, x, y, z, 1)
            } else {
                port_linear_eofl(j as c_double, *U_Integ1d.add(i), x, y, z, 0)
            };
            length += *S_Integ1d.add(i) * eofl;
        }
    }
    length
}

#[export_name = "LinearArea"]
pub unsafe extern "C" fn port_linear_area(geometry: *const Geometry) -> c_double {
    let base = (*geometry).entry as *const c_double;
    let x = base.add(0 * 2);
    let y = base.add(1 * 2);
    let z = base.add(2 * 2);
    let mut area = 0.0;
    for i in 0..N_Integ1d as usize {
        let eofa = port_linear_eofa(*U_Integ1d.add(i), 0.0, x, y, z);
        area += *S_Integ1d.add(i) * eofa;
    }
    area
}

#[export_name = "LinearIntegrateDiffToArea"]
pub unsafe extern "C" fn port_linear_integrate_diff_to_area(
    gb: *const Geometry, _cyl: *const Cylinder,
    fx: c_double, fy: c_double, _fz: c_double,
    nfx: c_double, nfy: c_double, _nfz: c_double,
) -> c_double {
    let base = (*gb).entry as *const c_double;
    let bx = base.add(0 * 2);
    let by = base.add(1 * 2);
    let nbx = base.add(3 * 2);
    let nby = base.add(4 * 2);
    let mut nfx = nfx;
    let mut nfy = nfy;
    let rs = nfx * nfx + nfy * nfy;
    if rs != 0.0 && abs_c(1.0 - rs) > 1.0e-9 {
        let r = 1.0 / sqrt(rs);
        let _ = r;
        nfx *= rs;
        nfy *= rs;
    }
    let mut f = 0.0;
    let mut cosa = 1.0;
    for i in 0..N_Integ1d as usize {
        let u = *U_Integ1d.add(i);
        let v = 0.0;
        let dx = linear_value(u, bx) - fx;
        let dy = linear_value(u, by) - fy;
        let r = sqrt(dx * dx + dy * dy);
        if rs > 0.0 {
            cosa = (dx * nfx + dy * nfy) / r;
            if cosa < 1.0e-9 {
                continue;
            }
        }
        let mut ntx = linear_value(u, nbx);
        let mut nty = linear_value(u, nby);
        let mut d = ntx * ntx + nty * nty;
        if abs_c(1.0 - r) > 1.0e-9 {
            d = 1.0 / sqrt(d);
            ntx *= d;
            nty *= d;
        }
        let cosb = (-dx * ntx - dy * nty) / r;
        if cosb < 1.0e-9 {
            continue;
        }
        let ea = port_linear_eofa(u, v, bx, by, base.add(2 * 2));
        f += ea * *S_Integ1d.add(i) * cosa * cosb / (2.0 * r);
    }
    f
}

#[allow(unused_assignments)]
#[export_name = "LinearViewFactor"]
pub unsafe extern "C" fn port_linear_view_factor(ga: *const Geometry, gb: *const Geometry) -> c_double {
    let abase = (*ga).entry as *const c_double;
    let bbase = (*gb).entry as *const c_double;
    let ax = abase.add(0 * 2);
    let ay = abase.add(1 * 2);
    let nax = abase.add(3 * 2);
    let nay = abase.add(4 * 2);
    let bx = bbase.add(0 * 2);
    let by = bbase.add(1 * 2);
    let nbx = bbase.add(3 * 2);
    let nby = bbase.add(4 * 2);
    let eps = 1e-8;
    let fx = linear_value(0.5, ax);
    let fy = linear_value(0.5, ay);
    let nfx = linear_value(0.5, nax);
    let nfy = linear_value(0.5, nay);
    let mut d = nfx * nfx + nfy * nfy;
    // Dead division kept verbatim: NTX/NTY are overwritten below.
    // (x86 FP division never traps; the quotient is unobservable.)
    #[allow(unused_assignments)]
    let (mut ntx, mut nty): (c_double, c_double) = (0.0, 0.0);
    if abs_c(1.0 - d) > eps {
        d = sqrt(d);
        ntx /= d;
        nty /= d;
    }
    ntx = linear_value(0.5, nbx);
    nty = linear_value(0.5, nby);
    d = ntx * ntx + nty * nty;
    if abs_c(1.0 - d) > eps {
        d = sqrt(d);
        ntx /= d;
        nty /= d;
    }
    let mut dx = linear_value(0.5, bx) - fx;
    let mut dy = linear_value(0.5, by) - fy;
    d = dx * dx + dy * dy;
    if d != 0.0 {
        d = 1.0 / sqrt(d);
        dx *= d;
        dy *= d;
    }
    let cosa = nfx * dx + nfy * dy;
    let cosb = -ntx * dx - nty * dy;
    if cosa < eps || cosb < eps {
        return 0.0;
    }
    dx = linear_value(0.0, bx) - linear_value(0.0, ax);
    dy = linear_value(0.0, by) - linear_value(0.0, ay);
    let mut fact = sqrt(dx * dx + dy * dy);
    dx = linear_value(1.0, bx) - linear_value(1.0, ax);
    dy = linear_value(1.0, by) - linear_value(1.0, ay);
    fact += sqrt(dx * dx + dy * dy);
    dx = linear_value(0.0, bx) - linear_value(1.0, ax);
    dy = linear_value(0.0, by) - linear_value(1.0, ay);
    fact -= sqrt(dx * dx + dy * dy);
    dx = linear_value(1.0, bx) - linear_value(0.0, ax);
    dy = linear_value(1.0, by) - linear_value(0.0, ay);
    fact -= sqrt(dx * dx + dy * dy);
    fact /= 2.0;
    abs_c(fact)
}

unsafe fn do_subdivide_aa(gb: *mut Geometry, levelb: c_int) {
    if (*gb).left.is_null() {
        port_linear_subdivide(gb, levelb, 1);
    }
    if (*gb).flags & GEOMETRY_FLAG_LEAF != 0 {
        (*gb).flags &= !GEOMETRY_FLAG_LEAF;
        (*(*gb).left).flags |= GEOMETRY_FLAG_LEAF;
        (*(*gb).right).flags |= GEOMETRY_FLAG_LEAF;
        if (*(*gb).left).area == 0.0 {
            (*(*gb).left).area = port_linear_area((*gb).left);
            (*(*gb).right).area = port_linear_area((*gb).right);
        }
    }
}

#[export_name = "LinearComputeViewFactors"]
pub unsafe extern "C" fn port_linear_compute_view_factors(
    ga: *mut Geometry, gb: *mut Geometry, levela: c_int, levelb: c_int,
) {
    let pi = 2.0 * acos(0.0);
    let _ = pi;
    let x = ((*ga).entry as *const c_double).add(0 * 2);
    let y = ((*ga).entry as *const c_double).add(1 * 2);
    let nx = ((*ga).entry as *const c_double).add(3 * 2);
    let ny = ((*ga).entry as *const c_double).add(4 * 2);
    // Initialized: every read follows a real assignment on reachable paths.
    let (mut fa, mut fb): (c_double, c_double) = (0.0, 0.0);
    let mut subdivide = false;
    if levela & 1 != 0 {
        fa = 0.0;
        fb = 1.0;
        subdivide = true;
    } else if (levelb & 1 != 0) && (*gb).geometrytype != GEOMETRY_TRIANGLE {
        fa = 1.0;
        fb = 0.0;
        subdivide = true;
    }
    if !subdivide {
        // (sole path into the block below is fallthrough)
        let u = 0.5;
        let fx = linear_value(u, x);
        let fy = linear_value(u, y);
        let fz = 0.0;
        let dx = linear_value(u, nx);
        let dy = linear_value(u, ny);
        let dz = 0.0;
        let integrate: IntegrateFn =
            std::mem::transmute(IntegrateDiffToArea[(*gb).geometrytype as usize]);
        fa = integrate(gb, std::ptr::null(), fx, fy, fz, dx, dy, dz);
        fb = fa;
        if ga != gb {
            let (mut u, mut v) = (0.5, 0.0);
            if (*gb).geometrytype == GEOMETRY_TRIANGLE {
                u = 1.0 / 3.0;
                v = 1.0 / 3.0;
            }
            let fx = FunctionValue(gb, u, v, 0);
            let fy = FunctionValue(gb, u, v, 1);
            let fz = 0.0;
            let dx = FunctionValue(gb, u, v, 3);
            let dy = FunctionValue(gb, u, v, 4);
            let dz = 0.0;
            fb = port_linear_integrate_diff_to_area(ga, std::ptr::null(), fx, fy, fz, dx, dy, dz);
        }
        if fa < 1.0e-10 && fb < 1.0e-10 {
            return;
        }
        if (fa < FactorEPS || (*gb).area < AreaEPS) && (fb < FactorEPS || (*ga).area < AreaEPS) {
            let mut hit = Nrays as c_double;
            for _ in 0..Nrays {
                let mut u = vrand();
                let mut v = vrand();
                let fx = linear_value(u, x);
                let fy = linear_value(u, y);
                let fz = 0.0;
                let w = u;
                u = 1.0 - v;
                v = 1.0 - w;
                if (*gb).geometrytype == GEOMETRY_TRIANGLE {
                    while u + v > 1.0 {
                        u = 1.0 - u;
                        v = 1.0 - v;
                    }
                }
                let dx = FunctionValue(gb, u, v, 0) - fx;
                let dy = FunctionValue(gb, u, v, 1) - fy;
                let dz = 0.0;
                if RayHitGeometry(fx, fy, fz, dx, dy, dz) != 0 {
                    hit -= 1.0;
                }
            }
            if hit == 0.0 {
                return;
            }
            if hit == Nrays as c_double
                || ((fa < FactorEPS / 2.0 || (*gb).area < AreaEPS / 2.0)
                    && (fb < FactorEPS / 2.0 || (*ga).area < AreaEPS / 2.0))
            {
                // NOTE: F is indeterminate when GA == GB (upstream); 0.0 here.
                let mut f: c_double = 0.0;
                if ga != gb {
                    f = port_linear_view_factor(ga, gb);
                }
                f = hit * f / (1.0 * Nrays as c_double);
                fb = f / (*gb).area;
                fa = f / (*ga).area;
                let link = calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
                (*link).next = (*ga).link;
                (*ga).link = link;
                (*link).entry = gb;
                (*link).viewfactor = fa;
                if ga != gb {
                    let link =
                        calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
                    (*link).next = (*gb).link;
                    (*gb).link = link;
                    (*link).entry = ga;
                    (*link).viewfactor = fb;
                }
                return;
            }
        }
        // Falls through to subdivide when the criteria above do not return.
    }
    if ga == gb {
        do_subdivide_aa(gb, levelb);
        port_linear_compute_view_factors((*gb).left, (*gb).left, levela + 1, levelb + 1);
        port_linear_compute_view_factors((*gb).right, (*gb).right, levela + 1, levelb + 1);
        port_linear_compute_view_factors((*gb).left, (*gb).right, levela + 1, levelb + 1);
    } else if fa > fb {
        let gb_type = (*gb).geometrytype as usize;
        if (*gb).left.is_null() {
            let sub: SubdivideFn = std::mem::transmute(Subdivide[gb_type]);
            sub(gb, levelb, 1);
        }
        if (*gb).flags & GEOMETRY_FLAG_LEAF != 0 {
            (*gb).flags &= !GEOMETRY_FLAG_LEAF;
            (*(*gb).left).flags |= GEOMETRY_FLAG_LEAF;
            (*(*gb).right).flags |= GEOMETRY_FLAG_LEAF;
            if (*(*gb).left).area == 0.0 {
                let area: AreaFn = std::mem::transmute(AreaCompute[gb_type]);
                (*(*gb).left).area = area((*gb).left);
                (*(*gb).right).area = area((*gb).right);
            }
        }
        port_linear_compute_view_factors(ga, (*gb).left, levela, levelb + 1);
        port_linear_compute_view_factors(ga, (*gb).right, levela, levelb + 1);
    } else {
        do_subdivide_aa(ga, levela);
        port_linear_compute_view_factors((*ga).left, gb, levela + 1, levelb);
        port_linear_compute_view_factors((*ga).right, gb, levela + 1, levelb);
    }
}

#[allow(unused_assignments)]
#[export_name = "LinearComputeRadiatorFactors"]
pub unsafe extern "C" fn port_linear_compute_radiator_factors(
    ga: *mut Geometry, _lineflag: c_int,
    dx: c_double, dy: c_double, dz: c_double,
    nx: c_double, ny: c_double, nz: c_double,
    levela: c_int,
) {
    let pi = 2.0 * acos(0.0);
    let base = (*ga).entry as *const c_double;
    let x = base.add(0 * 2);
    let y = base.add(1 * 2);
    // Initialized: every read follows a real assignment on reachable paths.
    let (mut fa, mut fb): (c_double, c_double) = (0.0, 0.0);
    let mut subdivide = false;
    if levela & 1 != 0 {
        fa = 0.0;
        fb = 0.0;
        subdivide = true;
    }
    if !subdivide {
        // (sole path into the block below is fallthrough)
        fa = port_linear_integrate_diff_to_area(ga, std::ptr::null(), dx, dy, dz, nx, ny, nz);
        fb = fa;
        if fa < 1.0e-10 && fb < 1.0e-10 {
            return;
        }
        if fa < FactorEPS || (*ga).area < AreaEPS {
            let mut hit = Nrays as c_double;
            for _ in 0..Nrays {
                let u = vrand();
                let fx = linear_value(u, x);
                let fy = linear_value(u, y);
                let fz = 0.0;
                let gx = dx - fx;
                let gy = dy - fy;
                let gz = 0.0;
                if RayHitGeometry(fx, fy, fz, gx, gy, gz) != 0 {
                    hit -= 1.0;
                }
            }
            if hit == 0.0 {
                return;
            }
            if hit == Nrays as c_double || fa < FactorEPS / 2.0 || (*ga).area < AreaEPS / 2.0 {
                let f = hit * fa / (1.0 * Nrays as c_double);
                fa = f / (*ga).area / pi;
                fb = fa;
                let link = calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
                (*link).next = (*ga).link;
                (*ga).link = link;
                (*link).entry = ga;
                (*link).viewfactor = fa;
                return;
            }
        }
        subdivide = true;
    }
    do_subdivide_aa(ga, levela);
    port_linear_compute_radiator_factors(
        (*ga).left, _lineflag, dx, dy, dz, nx, ny, nz, levela + 1,
    );
    port_linear_compute_radiator_factors(
        (*ga).right, _lineflag, dx, dy, dz, nx, ny, nz, levela + 1,
    );
    let _ = fb;
}
