#![allow(static_mut_refs)]
//! `static mut` mirrors C file-static mutable state, accessed the same
//! single-threaded way the original is.

// Port of `fem/src/view3d/16node_quad.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/16node_quad.c. Crate name `_16node_quad` is the
// exact `scaffold_crate_name("16node_quad")` (leading-digit rule).
//
// True object symbols: `elm_16node_quad_shape_functions`,
// `elm_16node_quad_fit`, `elm_16node_quad_fvalue`,
// `elm_16node_quad_nvalue`, `elm_16node_quad_dndu_fvalue`,
// `elm_16node_quad_dndu_nvalue`, `elm_16node_quad_dndv_fvalue`,
// `elm_16node_quad_dndv_nvalue`, `ElementOfArea_16node`.
// (`elm_16node_quad_add` is only declared, never defined, in this unit.)
//
// Faithfulness notes: same as `_2d_4node` — static N/A state, extern
// `mtrinv`, libc `fprintf` for `%g` paths, verbatim evaluation order.

use std::os::raw::{c_char, c_double, c_int, c_void};

static mut N: [[c_double; 16]; 16] = [[0.0; 16]; 16];
static mut A: [[c_double; 16]; 16] = [[0.0; 16]; 16];
static NODE_U: [c_double; 16] = [
    0.0, 1.0, 1.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0,
    2.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0,
];
static NODE_V: [c_double; 16] = [
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0,
    1.0, 1.0, 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0,
];

#[link(name = "c")]
extern "C" {
    static stderr: *mut c_void;
    fn fprintf(stream: *mut c_void, format: *const c_char, ...) -> c_int;
    fn mtrinv(a: *mut c_double, n: c_int);
}

#[export_name = "elm_16node_quad_shape_functions"]
pub unsafe extern "C" fn port_shape_functions(b: *mut c_double) {
    for i in 0..16 {
        let u = NODE_U[i];
        let v = NODE_V[i];
        A[i][0] = 1.0;
        A[i][1] = u;
        A[i][2] = u * u;
        A[i][3] = u * u * u;
        A[i][4] = v;
        A[i][5] = u * v;
        A[i][6] = u * u * v;
        A[i][7] = u * u * u * v;
        A[i][8] = v * v;
        A[i][9] = u * v * v;
        A[i][10] = u * u * v * v;
        A[i][11] = u * u * u * v * v;
        A[i][12] = v * v * v;
        A[i][13] = u * v * v * v;
        A[i][14] = u * u * v * v * v;
        A[i][15] = u * u * u * v * v * v;
    }
    mtrinv(A.as_mut_ptr() as *mut c_double, 16);
    let mut p = b;
    for i in 0..16 {
        for j in 0..16 {
            N[i][j] = A[j][i];
            *p = N[i][j];
            p = p.add(1);
        }
    }
}

#[export_name = "elm_16node_quad_fit"]
pub unsafe extern "C" fn port_fit(
    x: *const c_double, y: *const c_double, z: *const c_double, r: *mut c_double,
) {
    for i in 0..16 {
        let u = *x.add(i);
        let v = *y.add(i);
        A[i][0] = 1.0;
        A[i][1] = u;
        A[i][2] = v;
        A[i][3] = u * u;
        A[i][4] = u * v;
        A[i][5] = v * v;
        A[i][6] = u * u * u;
        A[i][7] = u * u * v;
        A[i][8] = u * v * v;
        A[i][9] = v * v * v;
        A[i][10] = u * u * u * v;
        A[i][11] = u * u * v * v;
        A[i][12] = u * v * v * v;
        A[i][13] = u * u * u * v * v;
        A[i][14] = u * u * v * v * v;
        A[i][15] = u * u * u * v * v * v;
    }
    mtrinv(A.as_mut_ptr() as *mut c_double, 16);
    for i in 0..16 {
        let mut s = 0.0;
        for j in 0..16 {
            s += A[i][j] * *z.add(j);
        }
        *r.add(i) = s;
    }
}

#[export_name = "elm_16node_quad_fvalue"]
pub unsafe extern "C" fn port_fvalue(f: *const c_double, u: c_double, v: c_double) -> c_double {
    let uv = u * v;
    let uu = u * u;
    let vv = v * v;
    let uuu = u * u * u;
    let vvv = v * v * v;
    let uuv = u * u * v;
    let uvv = u * v * v;
    let uuuv = uuu * v;
    let uvvv = u * vvv;
    let uuvv = uu * vv;
    let uuuvv = uuu * vv;
    let uuvvv = uu * vvv;
    let uuuvvv = uuu * vvv;
    let mut r = 0.0;
    for i in 0..16 {
        r += *f.add(i)
            * (N[i][0] + N[i][1] * u + N[i][2] * uu + N[i][3] * uuu + N[i][4] * v
                + N[i][5] * uv + N[i][6] * uuv + N[i][7] * uuuv + N[i][8] * vv
                + N[i][9] * uvv + N[i][10] * uuvv + N[i][11] * uuuvv + N[i][12] * vvv
                + N[i][13] * uvvv + N[i][14] * uuvvv + N[i][15] * uuuvvv);
    }
    r
}

#[export_name = "elm_16node_quad_nvalue"]
pub unsafe extern "C" fn port_nvalue(f: *mut c_double, u: c_double, v: c_double) {
    let uv = u * v;
    let uu = u * u;
    let vv = v * v;
    let uuu = u * u * u;
    let vvv = v * v * v;
    let uuv = u * u * v;
    let uvv = u * v * v;
    let uuuv = uuu * v;
    let uvvv = u * vvv;
    let uuvv = uu * vv;
    let uuuvv = uuu * vv;
    let uuvvv = uu * vvv;
    let uuuvvv = uuu * vvv;
    // NOTE: the original computes into a local R and assigns F[i] from the
    // parenthesized sum (different parenthesization than fvalue above —
    // kept per-function verbatim).
    for i in 0..16 {
        *f.add(i) = N[i][0] + N[i][1] * u + N[i][2] * v + N[i][3] * uu + N[i][4] * uv
            + N[i][5] * vv + N[i][6] * uuu + N[i][7] * uuv + N[i][8] * uvv
            + N[i][9] * vvv + N[i][10] * uuuv + N[i][11] * uuvv + N[i][12] * uvvv
            + N[i][13] * uuuvv + N[i][14] * uuvvv + N[i][15] * uuuvvv;
    }
}

#[export_name = "elm_16node_quad_dndu_fvalue"]
pub unsafe extern "C" fn port_dndu_fvalue(
    f: *const c_double, u: c_double, v: c_double,
) -> c_double {
    let u2 = 2.0 * u;
    let uu3 = 3.0 * u * u;
    let vv = v * v;
    let u2v = 2.0 * u * v;
    let uu3v = uu3 * v;
    let u2vv = u2 * vv;
    let vvv = vv * v;
    let uu3vv = uu3 * vv;
    let u2vvv = u2 * vvv;
    let uu3vvv = uu3 * vvv;
    let mut r = 0.0;
    for i in 0..16 {
        r += *f.add(i)
            * (N[i][1] + N[i][3] * u2 + N[i][4] * v + N[i][6] * uu3 + N[i][7] * u2v
                + N[i][8] * vv + N[i][10] * uu3v + N[i][11] * u2vv + N[i][12] * vvv
                + N[i][13] * uu3vv + N[i][14] * u2vvv + N[i][15] * uu3vvv);
    }
    r
}

#[export_name = "elm_16node_quad_dndu_nvalue"]
pub unsafe extern "C" fn port_dndu_nvalue(f: *mut c_double, u: c_double, v: c_double) {
    let u2 = 2.0 * u;
    let uu3 = 3.0 * u * u;
    let vv = v * v;
    let u2v = 2.0 * u * v;
    let uu3v = uu3 * v;
    let u2vv = u2 * vv;
    let vvv = vv * v;
    let uu3vv = uu3 * vv;
    let u2vvv = u2 * vvv;
    let uu3vvv = uu3 * vvv;
    for i in 0..16 {
        *f.add(i) = N[i][1] + N[i][3] * u2 + N[i][4] * v + N[i][6] * uu3 + N[i][7] * u2v
            + N[i][8] * vv + N[i][10] * uu3v + N[i][11] * u2vv + N[i][12] * vvv
            + N[i][13] * uu3vv + N[i][14] * u2vvv + N[i][15] * uu3vvv;
    }
}

#[export_name = "elm_16node_quad_dndv_fvalue"]
pub unsafe extern "C" fn port_dndv_fvalue(
    f: *const c_double, u: c_double, v: c_double,
) -> c_double {
    let v2 = 2.0 * v;
    let uu = u * u;
    let uv2 = u * v2;
    let vv3 = 3.0 * v * v;
    let uuu = uu * u;
    let uuv2 = uu * v2;
    let uvv3 = u * vv3;
    let uuuv2 = uuu * v2;
    let uuvv3 = uu * vv3;
    let uuuvv3 = uuu * vv3;
    let mut r = 0.0;
    for i in 0..16 {
        r += *f.add(i)
            * (N[i][2] + N[i][4] * u + N[i][5] * v2 + N[i][7] * uu + N[i][8] * uv2
                + N[i][9] * vv3 + N[i][10] * uuu + N[i][11] * uuv2 + N[i][12] * uvv3
                + N[i][13] * uuuv2 + N[i][14] * uuvv3 + N[i][15] * uuuvv3);
    }
    r
}

#[export_name = "elm_16node_quad_dndv_nvalue"]
pub unsafe extern "C" fn port_dndv_nvalue(f: *mut c_double, u: c_double, v: c_double) {
    let v2 = 2.0 * v;
    let uu = u * u;
    let uv2 = u * v2;
    let vv3 = 3.0 * v * v;
    let uuu = uu * u;
    let uuv2 = uu * v2;
    let uvv3 = u * vv3;
    let uuuv2 = uuu * v2;
    let uuvv3 = uu * vv3;
    let uuuvv3 = uuu * vv3;
    for i in 0..16 {
        *f.add(i) = N[i][2] + N[i][4] * u + N[i][5] * v2 + N[i][7] * uu + N[i][8] * uv2
            + N[i][9] * vv3 + N[i][10] * uuu + N[i][11] * uuv2 + N[i][12] * uvv3
            + N[i][13] * uuuv2 + N[i][14] * uuvv3 + N[i][15] * uuuvv3;
    }
}

#[export_name = "ElementOfArea_16node"]
pub unsafe extern "C" fn port_element_of_area(
    x: *const c_double, y: *const c_double, z: *const c_double,
    u: c_double, v: c_double,
) -> c_double {
    let dxdu = port_dndu_fvalue(x, u, v);
    let dydu = port_dndu_fvalue(y, u, v);
    let dzdu = port_dndu_fvalue(z, u, v);
    let dxdv = port_dndv_fvalue(x, u, v);
    let dydv = port_dndv_fvalue(y, u, v);
    let dzdv = port_dndv_fvalue(z, u, v);
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let deta = auu * avv - auv * auv;
    if deta < 1.0e-12 {
        for i in 0..16 {
            fprintf(stderr, c" %8g %8g\n".as_ptr(), *x.add(i), *y.add(i));
        }
        fprintf(stderr, c"area: Element Jacobian singular at: %g %g\n".as_ptr(), u, v);
        return 0.0;
    }
    deta.sqrt()
}
