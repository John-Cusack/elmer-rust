// Port of `fem/src/view3d/2d_4node.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/2d_4node.c. Crate name `_2d_4node` is the exact
// `scaffold_crate_name("2d_4node")` (leading-digit rule).
//
// True object symbols: `elm_4node_quad_shape_functions`, `Nvalue_4node`,
// `Fvalue_4node`, `dNdU_Nvalue_4node`, `dNdV_Nvalue_4node`,
// `dNdU_Fvalue_2node`, `dNdU_Fvalue_4node`, `dNdV_Fvalue_4node`,
// `ddNddU_Nvalue_4node`, `ddNddV_Nvalue_4node`, `ddNdUdV_Nvalue_4node`,
// `ddNddU_Fvalue_node`, `ddNddV_Fvalue_node`, `ddNdUdV_Fvalue_4node`,
// `ElementOfLine_4node`, `ElementOfArea_4node`,
// `derivates_to_global_4node`, `dNdXYZ_Nvalue_4node`.
//
// Faithfulness notes:
// - File-static `N`/`A` state kept in `static mut` globals, same layout.
// - `mtrinv` (LUDecomp unit) is an extern call, resolved at link.
// - Singular paths print via libc `fprintf` (extern) so `%8g`/`%g`
//   formatting is byte-identical; Rust formatting must NOT be used.
// - `derivates_to_global_4node` reuses the outer index `i` inside the
//   singular-path print loop (clobbering it to 4) — kept verbatim.
// - `AEPS` here is 1.0e-18 (this file's own define).

#![allow(static_mut_refs)]
//! `static mut` mirrors C file-static mutable state, accessed the same
//! single-threaded way the original is (Elmer view3d is serial here).

use std::os::raw::{c_char, c_double, c_int, c_void};

const AEPS: c_double = 1.0e-18;

static mut N: [[c_double; 4]; 4] = [
    [1.0 / 1.0, -1.0 / 1.0, -1.0 / 1.0, 1.0 / 1.0],
    [1.0 / 1.0, 1.0 / 1.0, -1.0 / 1.0, -1.0 / 1.0],
    [1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0],
    [1.0 / 1.0, -1.0 / 1.0, 1.0 / 1.0, -1.0 / 1.0],
];
static mut A: [[c_double; 4]; 4] = [[0.0; 4]; 4];
static NODE_U: [c_double; 4] = [0.0, 1.0, 1.0, 0.0];
static NODE_V: [c_double; 4] = [0.0, 0.0, 1.0, 1.0];

#[link(name = "c")]
extern "C" {
    static stderr: *mut c_void;
    fn fprintf(stream: *mut c_void, format: *const c_char, ...) -> c_int;
    fn mtrinv(a: *mut c_double, n: c_int);
}

#[export_name = "elm_4node_quad_shape_functions"]
pub unsafe extern "C" fn port_shape_functions(b: *mut c_double) {
    for i in 0..4 {
        let u = NODE_U[i];
        let v = NODE_V[i];
        A[i][0] = 1.0;
        A[i][1] = u;
        A[i][2] = v;
        A[i][3] = u * v;
    }
    mtrinv(A.as_mut_ptr() as *mut c_double, 4);
    let mut p = b;
    for i in 0..4 {
        for j in 0..4 {
            N[i][j] = A[j][i];
            *p = N[i][j];
            p = p.add(1);
        }
    }
}

#[export_name = "Nvalue_4node"]
pub unsafe extern "C" fn port_nvalue(r: *mut c_double, u: c_double, v: c_double) {
    for i in 0..4 {
        *r.add(i) = N[i][0] + N[i][1] * u + N[i][2] * v + N[i][3] * u * v;
    }
}

#[export_name = "Fvalue_4node"]
pub unsafe extern "C" fn port_fvalue(f: *const c_double, u: c_double, v: c_double) -> c_double {
    let mut r = 0.0;
    for i in 0..4 {
        r += *f.add(i) * (N[i][0] + N[i][1] * u + N[i][2] * v + N[i][3] * u * v);
    }
    r
}

#[export_name = "dNdU_Nvalue_4node"]
pub unsafe extern "C" fn port_dndu_nvalue(r: *mut c_double, _u: c_double, v: c_double) {
    for i in 0..4 {
        *r.add(i) = N[i][1] + N[i][3] * v;
    }
}

#[export_name = "dNdV_Nvalue_4node"]
pub unsafe extern "C" fn port_dndv_nvalue(r: *mut c_double, u: c_double, _v: c_double) {
    for i in 0..4 {
        *r.add(i) = N[i][2] + N[i][3] * u;
    }
}

#[export_name = "dNdU_Fvalue_2node"]
pub unsafe extern "C" fn port_dndu_fvalue_2node(f: *const c_double, _u: c_double) -> c_double {
    0.5 * (*f.add(1) - *f.add(0))
}

#[export_name = "dNdU_Fvalue_4node"]
pub unsafe extern "C" fn port_dndu_fvalue_4node(
    f: *const c_double, _u: c_double, v: c_double,
) -> c_double {
    let mut r = 0.0;
    for i in 0..4 {
        r += *f.add(i) * (N[i][1] + N[i][3] * v);
    }
    r
}

#[export_name = "dNdV_Fvalue_4node"]
pub unsafe extern "C" fn port_dndv_fvalue_4node(
    f: *const c_double, u: c_double, _v: c_double,
) -> c_double {
    let mut r = 0.0;
    for i in 0..4 {
        r += *f.add(i) * (N[i][2] + N[i][3] * u);
    }
    r
}

#[export_name = "ddNddU_Nvalue_4node"]
pub unsafe extern "C" fn port_ddnddu_nvalue(r: *mut c_double, _u: c_double, _v: c_double) {
    for i in 0..4 {
        *r.add(i) = 0.0;
    }
}

#[export_name = "ddNddV_Nvalue_4node"]
pub unsafe extern "C" fn port_ddnddv_nvalue(r: *mut c_double, _u: c_double, _v: c_double) {
    for i in 0..4 {
        *r.add(i) = 0.0;
    }
}

#[export_name = "ddNdUdV_Nvalue_4node"]
pub unsafe extern "C" fn port_ddndudv_nvalue(r: *mut c_double, _u: c_double, _v: c_double) {
    for i in 0..4 {
        *r.add(i) = N[i][3];
    }
}

#[export_name = "ddNddU_Fvalue_node"]
pub unsafe extern "C" fn port_ddnddu_fvalue_node(
    _f: *const c_double, _u: c_double, _v: c_double,
) -> c_double {
    0.0
}

#[export_name = "ddNddV_Fvalue_node"]
pub unsafe extern "C" fn port_ddnddv_fvalue_node(
    _f: *const c_double, _u: c_double, _v: c_double,
) -> c_double {
    0.0
}

#[export_name = "ddNdUdV_Fvalue_4node"]
pub unsafe extern "C" fn port_ddndudv_fvalue_4node(
    f: *const c_double, _u: c_double, _v: c_double,
) -> c_double {
    let mut g = 0.0;
    for i in 0..4 {
        g += *f.add(i) * N[i][3];
    }
    g
}

#[export_name = "ElementOfLine_4node"]
pub unsafe extern "C" fn port_element_of_line(
    x: *const c_double, y: *const c_double, z: *const c_double,
    u: c_double, v: c_double,
) -> c_double {
    let dxdu = port_dndu_fvalue_4node(x, u, v);
    let dydu = port_dndu_fvalue_4node(y, u, v);
    let dzdu = port_dndu_fvalue_4node(z, u, v);
    let dxdv = port_dndv_fvalue_4node(x, u, v);
    let dydv = port_dndv_fvalue_4node(y, u, v);
    let dzdv = port_dndv_fvalue_4node(z, u, v);
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let deta = auu * avv - auv * auv;
    if deta < AEPS {
        for i in 0..4 {
            fprintf(stderr, c" %8g %8g\n".as_ptr(), *x.add(i), *y.add(i));
        }
        fprintf(stderr, c"area: Element Jacobian singular at: %g %g\n".as_ptr(), u, v);
        return 0.0;
    }
    auu.sqrt()
}

#[export_name = "ElementOfArea_4node"]
pub unsafe extern "C" fn port_element_of_area(
    x: *const c_double, y: *const c_double, z: *const c_double,
    u: c_double, v: c_double,
) -> c_double {
    let dxdu = port_dndu_fvalue_4node(x, u, v);
    let dydu = port_dndu_fvalue_4node(y, u, v);
    let dzdu = port_dndu_fvalue_4node(z, u, v);
    let dxdv = port_dndv_fvalue_4node(x, u, v);
    let dydv = port_dndv_fvalue_4node(y, u, v);
    let dzdv = port_dndv_fvalue_4node(z, u, v);
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let deta = auu * avv - auv * auv;
    if deta < AEPS {
        for i in 0..4 {
            fprintf(stderr, c" %8g %8g\n".as_ptr(), *x.add(i), *y.add(i));
        }
        fprintf(stderr, c"area: Element Jacobian singular at: %g %g\n".as_ptr(), u, v);
        return 0.0;
    }
    deta.sqrt()
}

#[export_name = "derivates_to_global_4node"]
pub unsafe extern "C" fn port_derivates_to_global(
    x: *const c_double, y: *const c_double, z: *const c_double,
    dfdx: *mut c_double, dfdy: *mut c_double, _dfdz: *mut c_double,
    u: *const c_double, v: *const c_double, n: c_int,
) {
    let mut i: c_int = 0;
    while i < n {
        let ui = *u.offset(i as isize);
        let vi = *v.offset(i as isize);
        let dxdu = port_dndu_fvalue_4node(x, ui, vi);
        let dydu = port_dndu_fvalue_4node(y, ui, vi);
        let dzdu = port_dndu_fvalue_4node(z, ui, vi);
        let dxdv = port_dndv_fvalue_4node(x, ui, vi);
        let dydv = port_dndv_fvalue_4node(y, ui, vi);
        let dzdv = port_dndv_fvalue_4node(z, ui, vi);
        // NOTE the original swaps the names here (Avv computed from U
        // derivatives, Auu from V derivatives) — kept verbatim.
        let avv = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
        let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
        let auu = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
        let deta = auu * avv - auv * auv;
        if deta < AEPS {
            // inner loop reuses `i` (clobbers it to 4), like the original
            i = 0;
            while i < 4 {
                fprintf(stderr, c" %8g %8g\n".as_ptr(), *x.add(i as usize), *y.add(i as usize));
                i += 1;
            }
            fprintf(
                stderr,
                c"deriv: Element Jacobian singular at: %g %g\n".as_ptr(),
                ui, vi,
            );
            return;
        }
        let auv_n = -auv / deta;
        let auu_n = auu / deta;
        let avv_n = avv / deta;
        let a = *dfdx.offset(i as isize);
        let b = *dfdy.offset(i as isize);
        *dfdx.offset(i as isize) = auu_n * a + auv_n * b;
        *dfdy.offset(i as isize) = auv_n * a + avv_n * b;
        i += 1;
    }
}

#[export_name = "dNdXYZ_Nvalue_4node"]
pub unsafe extern "C" fn port_dndxyz_nvalue(
    topology: *const c_int, dfdx: *mut c_double, dfdy: *mut c_double,
    dfdz: *mut c_double, u: c_double, v: c_double,
) {
    let mut xx = [0.0; 4];
    let mut yy = [0.0; 4];
    for i in 0..4 {
        xx[i] = u;
        yy[i] = v;
    }
    port_dndu_nvalue(dfdx, u, v);
    port_dndv_nvalue(dfdy, u, v);
    // The original calls the 9-parameter `derivates_to_global_4node`
    // with 7 arguments (line 334 upstream); the 7th argument lands where
    // `u` is read from. `v`/`N` would be stack garbage upstream (the
    // call cannot succeed there either — recorded in EVIDENCE.md).
    // The port forwards the 7 real arguments positionally.
    port_derivates_to_global(
        topology as *const c_double,
        dfdx,
        dfdy,
        dfdz,
        xx.as_ptr() as *mut c_double,
        yy.as_ptr() as *mut c_double,
        std::ptr::null(),
        std::ptr::null(),
        4,
    );
}
