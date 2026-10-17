// Port of `fem/src/view3d/VectorUtil.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/VectorUtil.c. True object symbols:
// `RotateVector`, `RotateVectorTranspose`,
// `GetMatrixToRotateVectorToZAxis`, `CylinderNormal`.
//
// Faithfulness notes:
// - `ABS(s)` is `((s)>0?(s):(-(s)))` (so ABS(0.0) is -0.0); kept verbatim.
// - `Matrix_t` (double[3][3]) decays to a row pointer at the ABI; the port
//   takes the same single pointer and indexes `m[i*3+j]`.
// - libm calls (atan2/sin/cos/sqrt) lower to the same libm.so routines;
//   the differential driver checks bitwise identity (see EVIDENCE.md).
// - `exit(1)` on a negative cylinder discriminant is kept via
//   `std::process::exit(1)`.

use std::os::raw::{c_double, c_int};

#[repr(C)]
pub struct Point {
    pub x: c_double,
    pub y: c_double,
    pub z: c_double,
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

#[inline]
fn abs_c(x: c_double) -> c_double {
    if x > 0.0 { x } else { -x }
}

#[inline]
unsafe fn m_at(m: *mut c_double, i: isize, j: isize) -> c_double {
    *m.offset(i * 3 + j)
}

#[inline]
unsafe fn m_set(m: *mut c_double, i: isize, j: isize, v: c_double) {
    *m.offset(i * 3 + j) = v;
}

#[export_name = "RotateVector"]
pub unsafe extern "C" fn port_rotate_vector(
    x: *mut c_double,
    y: *mut c_double,
    z: *mut c_double,
    m: *mut c_double,
) {
    let rx = *x;
    let ry = *y;
    let rz = *z;
    *x = m_at(m, 0, 0) * rx + m_at(m, 0, 1) * ry + m_at(m, 0, 2) * rz;
    *y = m_at(m, 1, 0) * rx + m_at(m, 1, 1) * ry + m_at(m, 1, 2) * rz;
    *z = m_at(m, 2, 0) * rx + m_at(m, 2, 1) * ry + m_at(m, 2, 2) * rz;
}

#[export_name = "RotateVectorTranspose"]
pub unsafe extern "C" fn port_rotate_vector_transpose(
    x: *mut c_double,
    y: *mut c_double,
    z: *mut c_double,
    m: *mut c_double,
) {
    let rx = *x;
    let ry = *y;
    let rz = *z;
    *x = m_at(m, 0, 0) * rx + m_at(m, 1, 0) * ry + m_at(m, 2, 0) * rz;
    *y = m_at(m, 0, 1) * rx + m_at(m, 1, 1) * ry + m_at(m, 2, 1) * rz;
    *z = m_at(m, 0, 2) * rx + m_at(m, 1, 2) * ry + m_at(m, 2, 2) * rz;
}

#[export_name = "GetMatrixToRotateVectorToZAxis"]
pub unsafe extern "C" fn port_get_matrix_to_z(
    x: c_double,
    y: c_double,
    z: c_double,
    matrix: *mut c_double,
    ident: *mut c_int,
) {
    let a = y.atan2(z);
    let r = y * a.sin() + z * a.cos();
    let b = x.atan2(r);

    if abs_c(a) < 1.0e-12 && abs_c(b) < 1.0e-12 {
        *ident = 1;
        for k in 0..9 {
            *matrix.offset(k) = 0.0;
        }
        m_set(matrix, 0, 0, 1.0);
        m_set(matrix, 1, 1, 1.0);
        m_set(matrix, 2, 2, 1.0);
        return;
    }

    m_set(matrix, 0, 0, b.cos());
    m_set(matrix, 0, 1, -(a.sin()) * b.sin());
    m_set(matrix, 0, 2, -(a.cos()) * b.sin());

    m_set(matrix, 1, 0, 0.0);
    m_set(matrix, 1, 1, a.cos());
    m_set(matrix, 1, 2, -(a.sin()));

    m_set(matrix, 2, 0, b.sin());
    m_set(matrix, 2, 1, a.sin() * b.cos());
    m_set(matrix, 2, 2, a.cos() * b.cos());

    *ident = 0;
}

#[export_name = "CylinderNormal"]
pub unsafe extern "C" fn port_cylinder_normal(
    fx: c_double,
    fy: c_double,
    fz: c_double,
    dx: c_double,
    dy: c_double,
    dz: c_double,
    cyl: *const Cylinder,
    nx: *mut c_double,
    ny: *mut c_double,
    nz: *mut c_double,
) {
    let eps = 1e-12;
    let mut fx = fx - (*cyl).centerpoint.x;
    let mut fy = fy - (*cyl).centerpoint.y;
    let mut fz = fz - (*cyl).centerpoint.z;
    let mut dx = dx;
    let mut dy = dy;
    let mut dz = dz;
    let r0 = (*cyl).radius;

    port_rotate_vector(
        &mut fx,
        &mut fy,
        &mut fz,
        (*cyl).rotationmatrix.as_ptr() as *mut c_double,
    );
    port_rotate_vector(
        &mut dx,
        &mut dy,
        &mut dz,
        (*cyl).rotationmatrix.as_ptr() as *mut c_double,
    );

    let a = dx * dx + dy * dy;
    let b = 2.0 * (fx * dx + fy * dy);
    let c = fx * fx + fy * fy - r0 * r0;

    let t: c_double;
    if a < eps {
        t = -c / b;
    } else if abs_c(b) < eps {
        let mut tt = (-c / a).sqrt();
        if tt < 0.0 || tt > 1.0 {
            tt = -tt;
        }
        t = tt;
    } else {
        let d = b * b - 4.0 * a * c;
        if d >= 0.0 {
            let ds = d.sqrt();
            let mut tt = (-b + ds) / (2.0 * a);
            if tt < 0.0 || tt > 1.0 {
                tt = (-b - ds) / (2.0 * a);
            }
            t = tt;
        } else {
            std::process::exit(1);
        }
    }

    if t >= 0.0 && t <= 1.0 {
        *nx = fx + t * dx;
        *ny = fy + t * dy;
        *nz = 0.0;
        port_rotate_vector_transpose(
            nx,
            ny,
            nz,
            (*cyl).rotationmatrix.as_ptr() as *mut c_double,
        );
        let r = (*nx * (*nx) + *ny * (*ny) + *nz * (*nz)).sqrt();
        *nx /= r;
        *ny /= r;
        *nz /= r;
    } else {
        *nx = 0.0;
        *ny = 0.0;
        *nz = 0.0;
    }
}
