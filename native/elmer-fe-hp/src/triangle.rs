//! Triangle (303) p basis on the equilateral reference triangle: the `H1Basis_*`
//! routines `TriangleNodalP`, `TriangleL`, `TriangleEdgeP`, `TriangleBubbleP`
//! and their derivatives, plus the linear-only `TriangleNodal`. Conventions:
//! see [`crate::types`].

use crate::legendre::{d_legendre_p, d_var_phi, legendre_p, var_phi};
use crate::types::pow_int;

/// `H1Basis_TriangleNodalP` (barycentric on the equilateral triangle).
pub fn triangle_nodal_p(u: f64, v: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    val[*nb] = c * (1.0 - u - d * v);
    val[*nb + 1] = c * (1.0 + u - d * v);
    val[*nb + 2] = d * v;
    *nb += 3;
}

/// `H1Basis_TriangleL(node, u, v)`: barycentric coordinate of local node
/// `node` (1-based).
pub fn triangle_l(node: usize, u: f64, v: f64) -> f64 {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    match node {
        1 => c * (1.0 - u - d * v),
        2 => c * (1.0 + u - d * v),
        3 => d * v,
        _ => panic!("triangle_l: node {node} outside 1..=3"),
    }
}

/// `H1Basis_dTriangleNodalP`.
pub fn d_triangle_nodal_p(_u: f64, _v: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let half = 1.0 / 2.0;
    let a = 1.0 / 3f64.sqrt();
    let b = -a / 2.0;
    grad[*nb][0] = -half;
    grad[*nb + 1][0] = half;
    grad[*nb + 2][0] = 0.0;
    grad[*nb][1] = b;
    grad[*nb + 1][1] = b;
    grad[*nb + 2][1] = a;
    *nb += 3;
}

/// `H1Basis_dTriangleL(node)`.
pub fn d_triangle_l(node: usize) -> [f64; 2] {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    match node {
        1 => [-c, c * -d],
        2 => [c, c * -d],
        3 => [0.0, d],
        _ => panic!("d_triangle_l: node {node} outside 1..=3"),
    }
}

/// `H1Basis_TriangleEdgeP`: edge functions `La*Lb*varphi_{j+1}(Lb-La)` for the
/// three edges, `pmax[i]-1` of each.
pub fn triangle_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 0..3 {
        let la = triangle_l(edgedir[i][0], u, v);
        let lb = triangle_l(edgedir[i][1], u, v);
        for j in 1..pmax[i] {
            val[*nb + j as usize - 1] = la * lb * var_phi(j as usize + 1, lb - la);
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_dTriangleEdgeP`.
pub fn d_triangle_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 0..3 {
        let d_la = d_triangle_l(edgedir[i][0]);
        let d_lb = d_triangle_l(edgedir[i][1]);
        let la = triangle_l(edgedir[i][0], u, v);
        let lb = triangle_l(edgedir[i][1], u, v);
        for j in 1..pmax[i] {
            let v_phi = var_phi(j as usize + 1, lb - la);
            let d_v_phi = d_var_phi(j as usize + 1, lb - la);
            let idx = *nb + j as usize - 1;
            grad[idx][0] = d_la[0] * lb * v_phi
                + la * d_lb[0] * v_phi
                + la * lb * d_v_phi * (d_lb[0] - d_la[0]);
            grad[idx][1] = d_la[1] * lb * v_phi
                + la * d_lb[1] * v_phi
                + la * lb * d_v_phi * (d_lb[1] - d_la[1]);
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_TriangleBubbleP`. With `localnumbers` the function is the
/// directed triangular-face bubble (Legendre form), without it the element
/// bubble (power form).
pub fn triangle_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    val: &mut [f64],
    nb: &mut usize,
    localnumbers: Option<[usize; 3]>,
) {
    match localnumbers {
        None => {
            let la = triangle_l(1, u, v);
            let lb = triangle_l(2, u, v);
            let lc = triangle_l(3, u, v);
            for i in 0..=pmax - 3 {
                for j in 1..=pmax - i - 2 {
                    val[*nb + j as usize - 1] =
                        la * lb * lc * pow_int(lb - la, i) * pow_int(2.0 * lc - 1.0, j - 1);
                }
                *nb += (pmax - i - 2).max(0) as usize;
            }
        }
        Some(ln) => {
            let la = triangle_l(ln[0], u, v);
            let lb = triangle_l(ln[1], u, v);
            let lc = triangle_l(ln[2], u, v);
            for i in 0..=pmax - 3 {
                for j in 1..=pmax - i - 2 {
                    val[*nb + j as usize - 1] = la
                        * lb
                        * lc
                        * legendre_p(i as usize, lb - la)
                        * legendre_p((j - 1) as usize, 2.0 * lc - 1.0);
                }
                *nb += (pmax - i - 2).max(0) as usize;
            }
        }
    }
}

/// `H1Basis_dTriangleBubbleP`.
pub fn d_triangle_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    localnumbers: Option<[usize; 3]>,
) {
    match localnumbers {
        None => {
            let c = 1.0 / 2.0;
            let d = -(3f64.sqrt()) / 6.0;
            let e = 3f64.sqrt() / 3.0;
            let la = triangle_l(1, u, v);
            let lb = triangle_l(2, u, v);
            let lc = triangle_l(3, u, v);
            for i in 0..=pmax - 3 {
                for j in 1..=pmax - i - 2 {
                    let lb_lai = pow_int(lb - la, i);
                    let lc_1n = pow_int(2.0 * lc - 1.0, j - 1);
                    let idx = *nb + j as usize - 1;
                    grad[idx][0] = -c * lb * lc * lb_lai * lc_1n
                        + la * c * lc * lb_lai * lc_1n
                        + la * lb * lc * f64::from(i) * pow_int(lb - la, i - 1) * lc_1n;
                    grad[idx][1] = d * lb * lc * lb_lai * lc_1n
                        + la * d * lc * lb_lai * lc_1n
                        + la * lb * e * lb_lai * lc_1n
                        + la * lb
                            * lc
                            * lb_lai
                            * f64::from(j - 1)
                            * pow_int(2.0 * lc - 1.0, j - 2)
                            * 2.0
                            * e;
                }
                *nb += (pmax - i - 2).max(0) as usize;
            }
        }
        Some(ln) => {
            let d_la = d_triangle_l(ln[0]);
            let d_lb = d_triangle_l(ln[1]);
            let d_lc = d_triangle_l(ln[2]);
            let la = triangle_l(ln[0], u, v);
            let lb = triangle_l(ln[1], u, v);
            let lc = triangle_l(ln[2], u, v);
            for i in 0..=pmax - 3 {
                for j in 1..=pmax - i - 2 {
                    let lb_la = lb - la;
                    let lc_1 = 2.0 * lc - 1.0;
                    let a = legendre_p(i as usize, lb_la);
                    let b = legendre_p((j - 1) as usize, lc_1);
                    let idx = *nb + j as usize - 1;
                    for c in 0..2 {
                        grad[idx][c] = d_la[c] * lb * lc * a * b
                            + la * d_lb[c] * lc * a * b
                            + la * lb * d_lc[c] * a * b
                            + la * lb * lc * d_legendre_p(i as usize, lb_la) * (d_lb[c] - d_la[c]) * b
                            + la * lb * lc * a * d_legendre_p((j - 1) as usize, lc_1) * 2.0 * d_lc[c];
                    }
                }
                *nb += (pmax - i - 2).max(0) as usize;
            }
        }
    }
}

/// `H1Basis_TriangleNodal`: linear-only nodal functions on the *standard*
/// (not equilateral) triangle `1-u-v, u, v`. Fortran quirk: written to slots
/// `0..3` and does not advance a counter.
pub fn triangle_nodal(u: f64, v: f64, val: &mut [f64]) {
    val[0] = 1.0 - u - v;
    val[1] = u;
    val[2] = v;
}

/// `H1Basis_dTriangleNodal`.
pub fn d_triangle_nodal(grad: &mut [[f64; 3]]) {
    grad[0][0] = -1.0;
    grad[1][0] = 1.0;
    grad[2][0] = 0.0;
    grad[0][1] = -1.0;
    grad[1][1] = 0.0;
    grad[2][1] = 1.0;
}
