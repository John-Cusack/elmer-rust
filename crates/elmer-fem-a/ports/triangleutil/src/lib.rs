// Port of `fem/src/view3d/TriangleUtil.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/TriangleUtil.c.
//
// True object symbols: `TriangleSubdivide`, `TriangleEofL`,
// `TriangleEofA`, `TriangleLength`, `TriangleArea`,
// `TriangleIntegrateDiffToArea`, `TriangleComputeViewFactors`,
// `TriangleComputeRadiatorFactors`.
//
// Faithfulness notes (see also lineutil, same family):
// - `TriangleValue`/`TrianglePartialU/V` header macros expanded verbatim.
// - `goto subdivide` / `goto store` structured with an explicit flow
//   enum; statement order inside each region is unchanged.
// - `Cylinder_t CylS` zero-initialized (upstream leaves Axis/Length/
//   IdentMatrix indeterminate; `CylinderNormal` reads only Radius,
//   CenterPoint, RotationMatrix — all assigned).
// - Cross-unit externs: `CylinderNormal`,
//   `GetMatrixToRotateVectorToZAxis` (vectorutil unit), `ContourPrepare`,
//   `ContourEvaluate` (contourinteg), `ShaftInit`, `ShaftCandidates`,
//   `ShaftCountAdd` (shaftcull), `ClipVisible` (clipshadow),
//   `RayHitGeometry`, `RayHitCandidates`, `FunctionValue`, `vrand`,
//   the `N/U/V/S_Integ1d` + `N/U/V/S_Integ3` tables, `AreaCompute`/
//   `IntegrateDiffToArea`/`Subdivide` tables, `AreaEPS`, `FactorEPS`,
//   `Nrays`, `ShaftRayCull`, `ShaftStats`, `ClipShadows`,
//   `ShapeFunctionMatrix3`. malloc/calloc/free via libc.

use std::os::raw::{c_double, c_int, c_void};

const GEOMETRY_TRIANGLE: c_int = 1;
const GEOMETRY_FLAG_LEAF: c_int = 2;
const SC_MAXCAND: usize = 256;

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
pub struct Point {
    pub x: c_double, pub y: c_double, pub z: c_double,
}

#[repr(C)]
pub struct Cylinder {
    pub radius: c_double,
    pub length: c_double,
    pub axis: Point,
    pub centerpoint: Point,
    pub identmatrix: c_int,
    pub rotationmatrix: [[c_double; 3]; 3],
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
    pub n: [[[c_double; 3]; 24]; 1],
    pub d: [c_double; 24],
    pub np: c_int,
    pub tol: c_double,
    pub a: ContourTarget,
    pub b: ContourTarget,
}

#[link(name = "c")]
extern "C" {
    static N_Integ1d: c_int;
    static U_Integ1d: *const c_double;
    static S_Integ1d: *const c_double;
    static N_Integ3: c_int;
    static U_Integ3: *const c_double;
    static V_Integ3: *const c_double;
    static S_Integ3: *const c_double;
    static mut AreaEPS: c_double;
    static mut FactorEPS: c_double;
    static mut Nrays: c_int;
    static mut ShaftRayCull: c_int;
    static mut ShaftStats: c_int;
    static mut ClipShadows: c_int;
    static mut ClosedFormInteg: c_int;
    static mut AreaCompute: [*mut c_void; 9];
    static mut IntegrateDiffToArea: [*mut c_void; 9];
    static mut Subdivide: [*mut c_void; 9];
    fn RayHitGeometry(
        fx: c_double, fy: c_double, fz: c_double,
        dx: c_double, dy: c_double, dz: c_double,
    ) -> c_int;
    fn RayHitCandidates(
        cand: *const c_int, nc: c_int,
        fx: c_double, fy: c_double, fz: c_double,
        dx: c_double, dy: c_double, dz: c_double,
    ) -> c_int;
    fn vrand() -> c_double;
    fn CylinderNormal(
        fx: c_double, fy: c_double, fz: c_double,
        dx: c_double, dy: c_double, dz: c_double,
        cyl: *const Cylinder, nx: *mut c_double, ny: *mut c_double, nz: *mut c_double,
    );
    fn GetMatrixToRotateVectorToZAxis(
        x: c_double, y: c_double, z: c_double,
        matrix: *mut c_double, ident: *mut c_int,
    );
    fn ContourPrepare(gb: *const Geometry, t: *mut ContourTarget) -> c_int;
    fn ContourEvaluate(
        t: *const ContourTarget,
        fx: c_double, fy: c_double, fz: c_double,
        nfx: c_double, nfy: c_double, nfz: c_double,
        result: *mut c_double,
    ) -> c_int;
    fn ShaftInit(s: *mut Shaft, a: *const ContourTarget, b: *const ContourTarget) -> c_int;
    fn ShaftCandidates(s: *const Shaft, list: *mut c_int, maxn: c_int) -> c_int;
    fn ShaftCountAdd(n: c_int);
    fn ClipVisible(
        tgt: *const ContourTarget, cand: *const c_int, nc: c_int,
        x: *const c_double, nf: *const c_double, result: *mut c_double,
    ) -> c_int;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn sqrt(x: c_double) -> c_double;
    fn acos(x: c_double) -> c_double;
}

#[inline]
fn triangle_value(u: c_double, v: c_double, nn: *const c_double) -> c_double {
    unsafe { *nn.add(2) * v + *nn.add(1) * u + *nn.add(0) }
}

#[inline]
fn triangle_partial_u(_u: c_double, _v: c_double, nn: *const c_double) -> c_double {
    unsafe { *nn.add(1) }
}

#[inline]
fn triangle_partial_v(_u: c_double, _v: c_double, nn: *const c_double) -> c_double {
    unsafe { *nn.add(2) }
}

type IntegrateFn = unsafe extern "C" fn(
    *const Geometry, *const Cylinder,
    c_double, c_double, c_double, c_double, c_double, c_double,
) -> c_double;
type AreaFn = unsafe extern "C" fn(*const Geometry) -> c_double;
type SubdivideFn = unsafe extern "C" fn(*mut Geometry, c_int, c_int);

const GEOMETRY_LINE: c_int = 0;
const GEOMETRY_BILINEAR: c_int = 2;
const GEOMETRY_BICUBIC: c_int = 3;
const GEOMETRY_BIQUADRATIC: c_int = 4;

/// Header-static `FunctionValue` (ViewFactors.h) inlined per translation
/// unit, exactly like the C object this port replaces (local symbol, no
/// cross-unit linkage).
unsafe fn function_value(geom: *const Geometry, u: c_double, v: c_double, n: c_int) -> c_double {
    let entry = (*geom).entry as *const c_double;
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

#[export_name = "TriangleSubdivide"]
pub unsafe extern "C" fn port_triangle_subdivide(
    geometry: *mut Geometry, _sublev: c_int, _where: c_int,
) {
    (*geometry).left = calloc(1, std::mem::size_of::<Geometry>()) as *mut Geometry;
    let left_tri = calloc(1, 6 * 3 * std::mem::size_of::<c_double>()) as *mut c_double;
    (*(*geometry).left).entry = left_tri as *mut c_void;
    (*geometry).right = calloc(1, std::mem::size_of::<Geometry>()) as *mut Geometry;
    let right_tri = calloc(1, 6 * 3 * std::mem::size_of::<c_double>()) as *mut c_double;
    (*(*geometry).right).entry = right_tri as *mut c_void;
    (*(*geometry).left).geometrytype = GEOMETRY_TRIANGLE;
    (*(*geometry).right).geometrytype = GEOMETRY_TRIANGLE;
    let src = (*geometry).entry as *const c_double;
    let pf = |r: usize| src.add(r * 3);
    let x0 = triangle_value(0.0, 0.0, pf(0));
    let y0 = triangle_value(0.0, 0.0, pf(1));
    let z0 = triangle_value(0.0, 0.0, pf(2));
    let x1 = triangle_value(1.0, 0.0, pf(0));
    let y1 = triangle_value(1.0, 0.0, pf(1));
    let z1 = triangle_value(1.0, 0.0, pf(2));
    let x2 = triangle_value(0.0, 1.0, pf(0));
    let y2 = triangle_value(0.0, 1.0, pf(1));
    let z2 = triangle_value(0.0, 1.0, pf(2));
    let l1 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0);
    let l2 = (x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0) + (z2 - z0) * (z2 - z0);
    let l3 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
    // (xn,yn,zn): midpoint of the longest edge; (lx,ly,lz) layout per branch.
    let (xn, yn, zn, left_cols, right_cols): (
        c_double, c_double, c_double,
        [[c_double; 3]; 3], [[c_double; 3]; 3],
    );
    if l1 > l2 && l1 > l3 {
        xn = triangle_value(0.5, 0.0, pf(0));
        yn = triangle_value(0.5, 0.0, pf(1));
        zn = triangle_value(0.5, 0.0, pf(2));
        left_cols = [[x0, y0, z0], [xn - x0, yn - y0, zn - z0], [x2 - x0, y2 - y0, z2 - z0]];
        right_cols = [[xn, yn, zn], [x1 - xn, y1 - yn, z1 - zn], [x2 - xn, y2 - yn, z2 - zn]];
    } else if l2 > l3 {
        xn = triangle_value(0.0, 0.5, pf(0));
        yn = triangle_value(0.0, 0.5, pf(1));
        zn = triangle_value(0.0, 0.5, pf(2));
        left_cols = [[x0, y0, z0], [x1 - x0, y1 - y0, z1 - z0], [xn - x0, yn - y0, zn - z0]];
        right_cols = [[xn, yn, zn], [x1 - xn, y1 - yn, z1 - zn], [x2 - xn, y2 - yn, z2 - zn]];
    } else {
        xn = triangle_value(0.5, 0.5, pf(0));
        yn = triangle_value(0.5, 0.5, pf(1));
        zn = triangle_value(0.5, 0.5, pf(2));
        left_cols = [[x0, y0, z0], [x1 - x0, y1 - y0, z1 - z0], [xn - x0, yn - y0, zn - z0]];
        right_cols = [[x0, y0, z0], [xn - x0, yn - y0, zn - z0], [x2 - x0, y2 - y0, z2 - z0]];
    }
    for r in 0..3 {
        for cc in 0..3 {
            *left_tri.add(r * 3 + cc) = left_cols[cc][r];
            *right_tri.add(r * 3 + cc) = right_cols[cc][r];
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let v = *src.add((i + 3) * 3 + j);
            *left_tri.add((i + 3) * 3 + j) = v;
            *right_tri.add((i + 3) * 3 + j) = v;
        }
    }
    let _ = (xn, yn, zn);
}

#[export_name = "TriangleEofL"]
pub unsafe extern "C" fn port_triangle_eofl(
    u: c_double, v: c_double,
    x: *const c_double, y: *const c_double, z: *const c_double, unotv: c_int,
) -> c_double {
    let dxdu = triangle_partial_u(u, v, x);
    let dxdv = triangle_partial_v(u, v, x);
    let dydu = triangle_partial_u(u, v, y);
    let dydv = triangle_partial_v(u, v, y);
    let dzdu = triangle_partial_u(u, v, z);
    let dzdv = triangle_partial_v(u, v, z);
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    if unotv != 0 {
        sqrt(auu)
    } else {
        sqrt(avv)
    }
}

#[export_name = "TriangleEofA"]
pub unsafe extern "C" fn port_triangle_eofa(
    u: c_double, v: c_double,
    x: *const c_double, y: *const c_double, z: *const c_double,
) -> c_double {
    let dxdu = triangle_partial_u(u, v, x);
    let dxdv = triangle_partial_v(u, v, x);
    let dydu = triangle_partial_u(u, v, y);
    let dydv = triangle_partial_v(u, v, y);
    let dzdu = triangle_partial_u(u, v, z);
    let dzdv = triangle_partial_v(u, v, z);
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let deta = auu * avv - auv * auv;
    sqrt(deta)
}

#[export_name = "TriangleLength"]
pub unsafe extern "C" fn port_triangle_length(geometry: *const Geometry, unotv: c_int) -> c_double {
    let base = (*geometry).entry as *const c_double;
    let x = base.add(0 * 3);
    let y = base.add(1 * 3);
    let z = base.add(2 * 3);
    let mut length = 0.0;
    for j in 0..2 {
        for i in 0..N_Integ1d as usize {
            let eofl = if unotv != 0 {
                port_triangle_eofl(*U_Integ1d.add(i), j as c_double, x, y, z, 1)
            } else {
                port_triangle_eofl(j as c_double, *U_Integ1d.add(i), x, y, z, 0)
            };
            length += *S_Integ1d.add(i) * eofl;
        }
    }
    length
}

#[export_name = "TriangleArea"]
pub unsafe extern "C" fn port_triangle_area(geometry: *const Geometry) -> c_double {
    let base = (*geometry).entry as *const c_double;
    let x = base.add(0 * 3);
    let y = base.add(1 * 3);
    let z = base.add(2 * 3);
    let mut area = 0.0;
    for i in 0..N_Integ3 as usize {
        let eofa = port_triangle_eofa(*U_Integ3.add(i), *V_Integ3.add(i), x, y, z);
        area += *S_Integ3.add(i) * eofa;
    }
    area
}

#[export_name = "TriangleIntegrateDiffToArea"]
pub unsafe extern "C" fn port_triangle_integrate_diff_to_area(
    gb: *const Geometry, cyl: *const Cylinder,
    fx: c_double, fy: c_double, fz: c_double,
    mut nfx: c_double, mut nfy: c_double, mut nfz: c_double,
) -> c_double {
    let base = (*gb).entry as *const c_double;
    let bx = base.add(0 * 3);
    let by = base.add(1 * 3);
    let bz = base.add(2 * 3);
    let nbx = base.add(3 * 3);
    let nby = base.add(4 * 3);
    let nbz = base.add(5 * 3);
    let mut r = nfx * nfx + nfy * nfy + nfz * nfz;
    if cyl.is_null() && r != 0.0 {
        r = 1.0 / sqrt(r);
        nfx *= r;
        nfy *= r;
        nfz *= r;
    }
    let rs = r;
    let (u0, v0) = (1.0 / 3.0, 1.0 / 3.0);
    let mut ntx = triangle_value(u0, v0, nbx);
    let mut nty = triangle_value(u0, v0, nby);
    let mut ntz = triangle_value(u0, v0, nbz);
    r = ntx * ntx + nty * nty + ntz * ntz;
    r = 1.0 / sqrt(r);
    ntx *= r;
    nty *= r;
    ntz *= r;
    if ClosedFormInteg != 0 && cyl.is_null() && rs != 0.0 {
        let mut ct = std::mem::zeroed::<ContourTarget>();
        let mut cf: c_double = 0.0;
        if ContourPrepare(gb, &mut ct) != 0
            && ContourEvaluate(&ct, fx, fy, fz, nfx, nfy, nfz, &mut cf) != 0
        {
            return cf;
        }
    }
    let mut f = 0.0;
    let mut cosa = 1.0;
    for i in 0..N_Integ3 as usize {
        let u = *U_Integ3.add(i);
        let v = *V_Integ3.add(i);
        let dx = triangle_value(u, v, bx) - fx;
        let dy = triangle_value(u, v, by) - fy;
        let dz = triangle_value(u, v, bz) - fz;
        let r = sqrt(dx * dx + dy * dy + dz * dz);
        if !cyl.is_null() {
            CylinderNormal(fx, fy, fz, dx, dy, dz, cyl, &mut nfx, &mut nfy, &mut nfz);
        }
        if rs != 0.0 {
            cosa = (dx * nfx + dy * nfy + dz * nfz) / r;
            if cosa < 1.0e-8 {
                continue;
            }
        }
        let cosb = (-dx * ntx - dy * nty - dz * ntz) / r;
        if cosb < 1.0e-8 {
            continue;
        }
        f += 2.0 * (*gb).area * cosa * cosb * *S_Integ3.add(i) / (r * r);
    }
    f
}

#[allow(clippy::too_many_arguments)]
unsafe fn subdivide_aa(gb: *mut Geometry, levelb: c_int, is_triangle: bool, gb_type: usize) {
    if (*gb).left.is_null() {
        if is_triangle {
            port_triangle_subdivide(gb, levelb, 1);
        } else {
            let sub: SubdivideFn = std::mem::transmute(Subdivide[gb_type]);
            sub(gb, levelb, 1);
        }
    }
    if (*gb).flags & GEOMETRY_FLAG_LEAF != 0 {
        (*gb).flags &= !GEOMETRY_FLAG_LEAF;
        (*(*gb).left).flags |= GEOMETRY_FLAG_LEAF;
        (*(*gb).right).flags |= GEOMETRY_FLAG_LEAF;
        if (*(*gb).left).area == 0.0 {
            if is_triangle {
                (*(*gb).left).area = port_triangle_area((*gb).left);
                (*(*gb).right).area = port_triangle_area((*gb).right);
            } else {
                let area: AreaFn = std::mem::transmute(AreaCompute[gb_type]);
                (*(*gb).left).area = area((*gb).left);
                (*(*gb).right).area = area((*gb).right);
            }
        }
    }
}

// C initializes locals unconditionally; the allow covers initializers that
// reads on every path dominate (same values C would hold there).
#[allow(unused_assignments)]
#[export_name = "TriangleComputeViewFactors"]
pub unsafe extern "C" fn port_triangle_compute_view_factors(
    ga: *mut Geometry, gb: *mut Geometry, levela: c_int, levelb: c_int,
) {
    let pi = 2.0 * acos(0.0);
    let abase = (*ga).entry as *const c_double;
    let ax = abase.add(0 * 3);
    let ay = abase.add(1 * 3);
    let az = abase.add(2 * 3);
    let nx = abase.add(3 * 3);
    let ny = abase.add(4 * 3);
    let nz = abase.add(5 * 3);
    if ga == gb {
        return;
    }
    let (u0, v0) = (1.0 / 3.0, 1.0 / 3.0);
    let fx = triangle_value(u0, v0, ax);
    let fy = triangle_value(u0, v0, ay);
    let fz = triangle_value(u0, v0, az);
    let dx = triangle_value(u0, v0, nx);
    let dy = triangle_value(u0, v0, ny);
    let dz = triangle_value(u0, v0, nz);
    let integrate: IntegrateFn =
        std::mem::transmute(IntegrateDiffToArea[(*gb).geometrytype as usize]);
    let fa = integrate(gb, std::ptr::null(), fx, fy, fz, dx, dy, dz);
    let mut fb = fa;
    if ga != gb {
        let (mut u, mut v) = (0.5, 0.5);
        if (*gb).geometrytype == GEOMETRY_TRIANGLE {
            u = 1.0 / 3.0;
            v = 1.0 / 3.0;
        }
        let fx = function_value(gb, u, v, 0);
        let fy = function_value(gb, u, v, 1);
        let fz = function_value(gb, u, v, 2);
        let dx = function_value(gb, u, v, 3);
        let dy = function_value(gb, u, v, 4);
        let dz = function_value(gb, u, v, 5);
        fb = port_triangle_integrate_diff_to_area(ga, std::ptr::null(), fx, fy, fz, dx, dy, dz);
    }
    if fa < 1.0e-10 && fb < 1.0e-10 {
        return;
    }
    if !((fa < FactorEPS || (*gb).area < AreaEPS) && (fb < FactorEPS || (*ga).area < AreaEPS)) {
        // Criteria not met: falls to `subdivide:`.
        return subdivide_branch(ga, gb, levela, levelb, fa, fb);
    }
    // Shaft block (computed once, reused by statistics/clipping/rays).
    let mut clipped = false;
    let mut have_shaft = false;
    let mut nc: c_int = -1;
    let mut cand = [0 as c_int; SC_MAXCAND];
    let mut ta = std::mem::zeroed::<ContourTarget>();
    let mut tb = std::mem::zeroed::<ContourTarget>();
    let mut sh = std::mem::zeroed::<Shaft>();
    if (ShaftRayCull != 0 || ShaftStats != 0 || ClipShadows != 0) && ga != gb {
        if ContourPrepare(ga, &mut ta) != 0
            && ContourPrepare(gb, &mut tb) != 0
            && ShaftInit(&mut sh, &ta, &tb) != 0
        {
            nc = ShaftCandidates(&sh, cand.as_mut_ptr(), SC_MAXCAND as c_int);
            have_shaft = nc >= 0;
        }
    }
    if ShaftStats != 0 {
        ShaftCountAdd(if have_shaft { nc } else { -3 });
    }
    if ClipShadows != 0 && have_shaft {
        // `goto subdivide` when partly blocked but factors not small.
        if !(nc > 0
            && !((fa < FactorEPS / 2.0 || (*gb).area < AreaEPS / 2.0)
                && (fb < FactorEPS / 2.0 || (*ga).area < AreaEPS / 2.0)))
        {
            let mut f: c_double = 0.0;
            clipped = true;
            for i in 0..N_Integ3 as usize {
                let mut xp = [0.0; 3];
                let mut nf = [0.0; 3];
                let u = *U_Integ3.add(i);
                let v = *V_Integ3.add(i);
                xp[0] = triangle_value(u, v, ax);
                xp[1] = triangle_value(u, v, ay);
                xp[2] = triangle_value(u, v, az);
                nf[0] = triangle_value(u, v, nx);
                nf[1] = triangle_value(u, v, ny);
                nf[2] = triangle_value(u, v, nz);
                let mut dr = nf[0] * nf[0] + nf[1] * nf[1] + nf[2] * nf[2];
                if dr <= 0.0 {
                    clipped = false;
                    break;
                }
                dr = 1.0 / sqrt(dr);
                nf[0] *= dr;
                nf[1] *= dr;
                nf[2] *= dr;
                let mut df: c_double = 0.0;
                if ClipVisible(&tb, cand.as_ptr(), nc, xp.as_ptr(), nf.as_ptr(), &mut df) == 0 {
                    clipped = false;
                    break;
                }
                f += *S_Integ3.add(i) * 2.0 * (*ga).area * df;
            }
            if clipped {
                f /= pi;
                if f <= 0.0 {
                    return;
                }
                // `goto store`
                do_store(ga, gb, f);
                return;
            }
            // Clipped false: fall through to the ray path below.
        } else {
            // `goto subdivide`
            return subdivide_branch(ga, gb, levela, levelb, fa, fb);
        }
        let _ = clipped;
    }
    // Ray path.
    let mut hit = Nrays as c_double;
    for _ in 0..Nrays {
        let (mut u, mut v) = (vrand(), vrand());
        while u + v > 1.0 {
            u = 1.0 - u;
            v = 1.0 - v;
        }
        let fx = triangle_value(u, v, ax);
        let fy = triangle_value(u, v, ay);
        let fz = triangle_value(u, v, az);
        let w = u;
        u = 1.0 - v;
        v = 1.0 - w;
        if (*gb).geometrytype == GEOMETRY_TRIANGLE {
            while u + v > 1.0 {
                u = 1.0 - u;
                v = 1.0 - v;
            }
        }
        let dx = function_value(gb, u, v, 0) - fx;
        let dy = function_value(gb, u, v, 1) - fy;
        let dz = function_value(gb, u, v, 2) - fz;
        if ShaftRayCull == 0 || !have_shaft {
            hit -= RayHitGeometry(fx, fy, fz, dx, dy, dz) as c_double;
        } else if nc > 0 {
            hit -= RayHitCandidates(cand.as_ptr(), nc, fx, fy, fz, dx, dy, dz) as c_double;
        }
    }
    if ShaftStats != 0 && have_shaft && nc == 0 && hit < Nrays as c_double {
        ShaftCountAdd(-2);
    }
    if hit == 0.0 {
        return;
    }
    if hit == Nrays as c_double
        || ((fa < FactorEPS / 2.0 || (*gb).area < AreaEPS / 2.0)
            && (fb < FactorEPS / 2.0 || (*ga).area < AreaEPS / 2.0))
    {
        let hitn = hit / Nrays as c_double;
        let mut ct = std::mem::zeroed::<ContourTarget>();
        let use_cf = ClosedFormInteg != 0 && ContourPrepare(gb, &mut ct) != 0;
        let mut f: c_double = 0.0;
        for i in 0..N_Integ3 as usize {
            let u = *U_Integ3.add(i);
            let v = *V_Integ3.add(i);
            let fx = triangle_value(u, v, ax);
            let fy = triangle_value(u, v, ay);
            let fz = triangle_value(u, v, az);
            let dx = triangle_value(u, v, nx);
            let dy = triangle_value(u, v, ny);
            let dz = triangle_value(u, v, nz);
            let df: c_double;
            if use_cf {
                let dr = dx * dx + dy * dy + dz * dz;
                if dr > 0.0 {
                    let dr = 1.0 / sqrt(dr);
                    let mut dfr: c_double = 0.0;
                    ContourEvaluate(&ct, fx, fy, fz, dx * dr, dy * dr, dz * dr, &mut dfr);
                    df = dfr;
                } else {
                    df = 0.0;
                }
            } else {
                let integrate: IntegrateFn =
                    std::mem::transmute(IntegrateDiffToArea[(*gb).geometrytype as usize]);
                df = integrate(gb, std::ptr::null(), fx, fy, fz, dx, dy, dz);
            }
            let ea = 2.0 * (*ga).area;
            f += *S_Integ3.add(i) * ea * df;
        }
        let f = hitn * f / pi;
        // `store:` label
        do_store(ga, gb, f);
        return;
    }
    // Falls to `subdivide:`.
    subdivide_branch(ga, gb, levela, levelb, fa, fb);
    let _ = fa;
}

/// Shared store block (`store:` label upstream).
unsafe fn do_store(ga: *mut Geometry, gb: *mut Geometry, f: c_double) {
    let fb = f / (*gb).area;
    let fa = f / (*ga).area;
    let link = calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
    (*link).next = (*ga).link;
    (*ga).link = link;
    (*link).entry = gb;
    (*link).viewfactor = fa;
    if ga != gb {
        let link = calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
        (*link).next = (*gb).link;
        (*gb).link = link;
        (*link).entry = ga;
        (*link).viewfactor = fb;
    }
}

#[allow(clippy::too_many_arguments)]
unsafe fn subdivide_branch(
    ga: *mut Geometry, gb: *mut Geometry, levela: c_int, levelb: c_int, fa: c_double, fb: c_double,
) {
    if ga == gb {
        subdivide_aa(gb, levelb, true, 0);
        port_triangle_compute_view_factors((*gb).left, (*gb).left, levela + 1, levelb + 1);
        port_triangle_compute_view_factors((*gb).right, (*gb).right, levela + 1, levelb + 1);
        port_triangle_compute_view_factors((*gb).left, (*gb).right, levela + 1, levelb + 1);
    } else if fa > fb {
        let gb_type = (*gb).geometrytype as usize;
        subdivide_aa(gb, levelb, (*gb).geometrytype == GEOMETRY_TRIANGLE, gb_type);
        port_triangle_compute_view_factors(ga, (*gb).left, levela, levelb + 1);
        port_triangle_compute_view_factors(ga, (*gb).right, levela, levelb + 1);
    } else {
        subdivide_aa(ga, levela, true, 0);
        port_triangle_compute_view_factors((*ga).left, gb, levela + 1, levelb);
        port_triangle_compute_view_factors((*ga).right, gb, levela + 1, levelb);
    }
}

#[allow(unused_assignments)]
#[export_name = "TriangleComputeRadiatorFactors"]
pub unsafe extern "C" fn port_triangle_compute_radiator_factors(
    ga: *mut Geometry, lineflag: c_int,
    dx: c_double, dy: c_double, dz: c_double,
    nx: c_double, ny: c_double, nz: c_double,
    levela: c_int,
) {
    let pi = 2.0 * acos(0.0);
    let base = (*ga).entry as *const c_double;
    let x = base.add(0 * 3);
    let y = base.add(1 * 3);
    let z = base.add(2 * 3);
    let (mut fa, mut fb): (c_double, c_double) = (0.0, 0.0);
    let mut cyl = std::mem::zeroed::<Cylinder>();
    let mut cylp: *const Cylinder = std::ptr::null();
    let mut ident: c_int = 0;
    let mut subdivide = false;
    if levela & 1 != 0 {
        subdivide = true;
    }
    if !subdivide {
        let r = nx * nx + ny * ny + nz * nz;
        if lineflag != 0 && r != 0.0 {
            let r = sqrt(r);
            cylp = &cyl;
            cyl.radius = r / 25.0;
            cyl.centerpoint.x = (2.0 * dx + nx) / 2.0;
            cyl.centerpoint.y = (2.0 * dy + ny) / 2.0;
            cyl.centerpoint.z = (2.0 * dz + nz) / 2.0;
            GetMatrixToRotateVectorToZAxis(nx / r, ny / r, nz / r, cyl.rotationmatrix.as_mut_ptr() as *mut c_double, &mut ident);
            fa = 0.0;
            for i in 0..N_Integ1d as usize {
                let u = *U_Integ1d.add(i);
                fa += *S_Integ1d.add(i)
                    * port_triangle_integrate_diff_to_area(
                        ga, cylp, dx + u * nx, dy + u * ny, dz + u * nz, nx, ny, nz,
                    );
            }
            fb = fa;
        } else {
            fa = port_triangle_integrate_diff_to_area(ga, cylp, dx, dy, dz, nx, ny, nz);
            fb = fa;
        }
        if fa < 1.0e-10 {
            return;
        }
        if fa < FactorEPS || (*ga).area < AreaEPS {
            let mut hit = Nrays as c_double;
            for _ in 0..Nrays {
                let (mut u, mut v) = (vrand(), vrand());
                while u + v > 1.0 {
                    u = 1.0 - u;
                    v = 1.0 - v;
                }
                let fx = triangle_value(u, v, x);
                let fy = triangle_value(u, v, y);
                let fz = triangle_value(u, v, z);
                let mut gx = dx - fx;
                let mut gy = dy - fy;
                let mut gz = dz - fz;
                if !cylp.is_null() {
                    u = vrand();
                    gx += u * nx;
                    gy += u * ny;
                    gz += u * nz;
                }
                if RayHitGeometry(fx, fy, fz, gx, gy, gz) != 0 {
                    hit -= 1.0;
                }
            }
            if hit == 0.0 {
                return;
            }
            if hit == Nrays as c_double || fa < FactorEPS / 2.0 || (*ga).area < AreaEPS / 2.0 {
                let f = hit * fa / (1.0 * Nrays as c_double);
                fa = f / (*ga).area / (4.0 * pi);
                fb = fa;
                let link = calloc(1, std::mem::size_of::<GeometryList>()) as *mut GeometryList;
                (*link).next = (*ga).link;
                (*ga).link = link;
                (*link).entry = ga;
                (*link).viewfactor = fa;
                return;
            }
        }
    }
    subdivide_aa(ga, levela, true, 0);
    port_triangle_compute_radiator_factors(
        (*ga).left, lineflag, dx, dy, dz, nx, ny, nz, levela + 1,
    );
    port_triangle_compute_radiator_factors(
        (*ga).right, lineflag, dx, dy, dz, nx, ny, nz, levela + 1,
    );
    let _ = (fb, ident);
}
