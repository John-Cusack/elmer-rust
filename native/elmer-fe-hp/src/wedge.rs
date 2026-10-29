//! Wedge (706) p-element basis: nodal, edge, face and bubble functions with
//! their gradients, in the full-polynomial and serendipity (`sd_*`) variants.
//!
//! Scalar transliteration of the wedge section of `fem/src/H1Basis.F90`
//! (`H1Basis_WedgeNodalP` ... `H1Basis_dWedgeBubbleP`, `H1Basis_SD_WedgeEdgeP`
//! ... `H1Basis_SD_dWedgeBubbleP`, `H1Basis_WedgeNodal`, `H1Basis_dWedgeNodal`).
//! See `types.rs` for the argument conventions. Local nodes are 1-based; the
//! reference wedge is the biunit prism with a regular triangle,
//! `u,w in [-1,1]`, `v in [0, sqrt(3)]`.
#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::many_single_char_names
)]

use crate::legendre::{d_legendre_p, d_phi, d_var_phi, legendre_p, phi, var_phi};

/// `nbasis = nbasis + delta` (the Fortran counter may in principle decrease).
fn advance(nb: &mut usize, delta: i32) {
    *nb = (*nb as isize + delta as isize) as usize;
}

/// `1/sqrt(3)` as the Fortran parameter `d`/`a`.
fn inv_sqrt3() -> f64 {
    1.0 / 3.0f64.sqrt()
}

/// `H1Basis_WedgeL(node, u, v)`: triangle barycentric coordinate of a wedge node.
pub fn wedge_l(node: usize, u: f64, v: f64) -> f64 {
    let a = inv_sqrt3();
    match node {
        1 | 4 => (1.0 - u - a * v) / 2.0,
        2 | 5 => (1.0 + u - a * v) / 2.0,
        3 | 6 => a * v,
        _ => panic!("wedge_l: bad node {node}"),
    }
}

/// `H1Basis_dWedgeL(node)`: gradient of [`wedge_l`].
pub fn d_wedge_l(node: usize) -> [f64; 3] {
    let c = 1.0 / 2.0;
    let d = inv_sqrt3();
    match node {
        1 | 4 => [-c, c * -d, c * 0.0],
        2 | 5 => [c * 1.0, c * -d, c * 0.0],
        3 | 6 => [0.0, d, 0.0],
        _ => panic!("d_wedge_l: bad node {node}"),
    }
}

/// `H1Basis_WedgeH(node, w)`: linear coordinate of a wedge node along `w`.
pub fn wedge_h(node: usize, w: f64) -> f64 {
    match node {
        1..=3 => -w / 2.0,
        4..=6 => w / 2.0,
        _ => panic!("wedge_h: bad node {node}"),
    }
}

/// `H1Basis_dWedgeH(node)`: gradient of [`wedge_h`].
pub fn d_wedge_h(node: usize) -> [f64; 3] {
    let c = 1.0 / 2.0;
    match node {
        1..=3 => [0.0, 0.0, -c],
        4..=6 => [0.0, 0.0, c],
        _ => panic!("d_wedge_h: bad node {node}"),
    }
}

/// `H1Basis_WedgeNodalP`: six nodal functions (regular-triangle wedge).
pub fn wedge_nodal_p(u: f64, v: f64, w: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 4.0;
    let d = inv_sqrt3();
    let e = 3.0f64.sqrt() / 6.0;
    let n = *nb;
    val[n] = c * (1.0 - u - d * v) * (1.0 - w);
    val[n + 1] = c * (1.0 + u - d * v) * (1.0 - w);
    val[n + 2] = e * v * (1.0 - w);
    val[n + 3] = c * (1.0 - u - d * v) * (1.0 + w);
    val[n + 4] = c * (1.0 + u - d * v) * (1.0 + w);
    val[n + 5] = e * v * (1.0 + w);
    *nb += 6;
}

/// `H1Basis_dWedgeNodalP`: gradients of [`wedge_nodal_p`].
pub fn d_wedge_nodal_p(u: f64, v: f64, w: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let c = 1.0 / 4.0;
    let d = inv_sqrt3();
    let e = 3.0f64.sqrt() / 6.0;
    let f = 3.0f64.sqrt() / 12.0;
    let n = *nb;
    grad[n][0] = -c * (1.0 - w);
    grad[n + 1][0] = c * (1.0 - w);
    grad[n + 2][0] = 0.0;
    grad[n + 3][0] = -c * (1.0 + w);
    grad[n + 4][0] = c * (1.0 + w);
    grad[n + 5][0] = 0.0;
    grad[n][1] = -f * (1.0 - w);
    grad[n + 1][1] = -f * (1.0 - w);
    grad[n + 2][1] = e * (1.0 - w);
    grad[n + 3][1] = -f * (1.0 + w);
    grad[n + 4][1] = -f * (1.0 + w);
    grad[n + 5][1] = e * (1.0 + w);
    grad[n][2] = -c * (1.0 - u - d * v);
    grad[n + 1][2] = -c * (1.0 + u - d * v);
    grad[n + 2][2] = -e * v;
    grad[n + 3][2] = c * (1.0 - u - d * v);
    grad[n + 4][2] = c * (1.0 + u - d * v);
    grad[n + 5][2] = e * v;
    *nb += 6;
}

/// `H1Basis_WedgeNodal`: linear wedge nodal functions on the unit-triangle
/// wedge (`H1Basis_WedgeNodal` writes from slot 0; "not a barycentric wedge").
pub fn wedge_nodal(u: f64, v: f64, w: f64, val: &mut [f64]) {
    let c = 1.0 / 2.0;
    val[0] = c * (1.0 - u - v) * (1.0 - w);
    val[1] = c * u * (1.0 - w);
    val[2] = c * v * (1.0 - w);
    val[3] = c * (1.0 - u - v) * (1.0 + w);
    val[4] = c * u * (1.0 + w);
    val[5] = c * v * (1.0 + w);
}

/// `H1Basis_dWedgeNodal`: gradients of [`wedge_nodal`].
pub fn d_wedge_nodal(u: f64, v: f64, w: f64, grad: &mut [[f64; 3]]) {
    let c = 1.0 / 2.0;
    grad[0][0] = -c * (1.0 - w);
    grad[1][0] = c * (1.0 - w);
    grad[2][0] = 0.0;
    grad[3][0] = -c * (1.0 + w);
    grad[4][0] = c * (1.0 + w);
    grad[5][0] = 0.0;
    grad[0][1] = -c * (1.0 - w);
    grad[1][1] = 0.0;
    grad[2][1] = c * (1.0 - w);
    grad[3][1] = -c * (1.0 + w);
    grad[4][1] = 0.0;
    grad[5][1] = c * (1.0 + w);
    grad[0][2] = -c * (1.0 - u - v);
    grad[1][2] = -c * u;
    grad[2][2] = -c * v;
    grad[3][2] = c * (1.0 - u - v);
    grad[4][2] = c * u;
    grad[5][2] = c * v;
}

/// `H1Basis_WedgeEdgeP`: edge functions; the nine edges are the six triangle
/// edges (1-6) then the three vertical edges (7-9). Reads the nodal values
/// the caller stored in `val[0..6]` (Fortran `fval(k,node)`).
pub fn wedge_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 1..=6 {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        let lav = wedge_l(node1, u, v);
        let lbv = wedge_l(node2, u, v);
        for j in 1..pmax[i - 1] {
            let la = lav;
            let lb = lbv;
            let na = val[node1 - 1];
            let nbv = val[node2 - 1];
            val[*nb + j as usize - 1] = na * nbv * var_phi(j as usize + 1, lb - la);
        }
        advance(nb, pmax[i - 1] - 1);
    }
    for i in 7..=9 {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        let lav = wedge_h(node1, w);
        let lbv = wedge_h(node2, w);
        for j in 1..pmax[i - 1] {
            let la = lav;
            let lb = lbv;
            let na = val[node1 - 1];
            let nbv = val[node2 - 1];
            val[*nb + j as usize - 1] = na * nbv * var_phi(j as usize + 1, lb - la);
        }
        advance(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_dWedgeEdgeP`: gradients of [`wedge_edge_p`]. Reads the nodal
/// gradients the caller stored in `grad[0..6]`.
pub fn d_wedge_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let mut n = [0.0f64; 6];
    let mut nnb = 0;
    wedge_nodal_p(u, v, w, &mut n, &mut nnb);
    for i in 1..=6 {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        let d_la = d_wedge_l(node1);
        let d_lb = d_wedge_l(node2);
        let lav = wedge_l(node1, u, v);
        let lbv = wedge_l(node2, u, v);
        for j in 1..pmax[i - 1] {
            let la = lav;
            let lb = lbv;
            let na = n[node1 - 1];
            let nbv = n[node2 - 1];
            let d_na = grad[node1 - 1];
            let d_nb = grad[node2 - 1];
            let v_phi = var_phi(j as usize + 1, lb - la);
            let d_v_phi = d_var_phi(j as usize + 1, lb - la);
            let idx = *nb + j as usize - 1;
            for c in 0..3 {
                grad[idx][c] = d_na[c] * nbv * v_phi
                    + na * d_nb[c] * v_phi
                    + na * nbv * d_v_phi * (d_lb[c] - d_la[c]);
            }
        }
        advance(nb, pmax[i - 1] - 1);
    }
    for i in 7..=9 {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        let d_la = d_wedge_h(node1);
        let d_lb = d_wedge_h(node2);
        let lav = wedge_h(node1, w);
        let lbv = wedge_h(node2, w);
        for j in 1..pmax[i - 1] {
            let la = lav;
            let lb = lbv;
            let v_phi = var_phi(j as usize + 1, lb - la);
            let d_v_phi = d_var_phi(j as usize + 1, lb - la);
            let na = n[node1 - 1];
            let nbv = n[node2 - 1];
            let d_na = grad[node1 - 1];
            let d_nb = grad[node2 - 1];
            let idx = *nb + j as usize - 1;
            for c in 0..3 {
                grad[idx][c] = d_na[c] * nbv * v_phi
                    + na * d_nb[c] * v_phi
                    + na * nbv * d_v_phi * (d_lb[c] - d_la[c]);
            }
        }
        advance(nb, pmax[i - 1] - 1);
    }
}

/// Whether the first two nodes of a square face form an edge of the lower or
/// upper triangle (`nonpermuted` in the Fortran).
fn square_face_nonpermuted(node1: usize, node2: usize) -> bool {
    ((1..=3).contains(&node1) && (1..=3).contains(&node2))
        || ((4..=6).contains(&node1) && (4..=6).contains(&node2))
}

/// `H1Basis_WedgeFaceP`: two triangle faces (1-2) then three square faces
/// (3-5). Reads the nodal values the caller stored in `val[0..6]`.
pub fn wedge_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let c = 1.0 / 2.0;
    for i in 1..=2 {
        let fd = facedir[i - 1];
        let lav = wedge_l(fd[0], u, v);
        let lbv = wedge_l(fd[1], u, v);
        let lcv = wedge_l(fd[2], u, v);
        let nav = wedge_h(fd[0], w);
        for j in 0..=(pmax[i - 1] - 3) {
            for k in 1..=(pmax[i - 1] - j - 2) {
                let (la, lb, lc, na) = (lav, lbv, lcv, nav);
                val[*nb + k as usize - 1] = c
                    * (1.0 + 2.0 * na)
                    * la
                    * lb
                    * lc
                    * legendre_p(j as usize, lb - la)
                    * legendre_p((k - 1) as usize, 2.0 * lc - 1.0);
            }
            *nb += (pmax[i - 1] - j - 2).max(0) as usize;
        }
    }
    for i in 3..=5 {
        let node1 = facedir[i - 1][0];
        let node2 = facedir[i - 1][1];
        let node3 = facedir[i - 1][2];
        let node4 = facedir[i - 1][3];
        let nonpermuted = square_face_nonpermuted(node1, node2);
        for j in 0..=(pmax[i - 1] - 2) {
            for k in 1..pmax[i - 1] {
                if nonpermuted {
                    let la = wedge_l(node1, u, v);
                    let lb = wedge_l(node2, u, v);
                    let lha = wedge_h(node1, w);
                    let lhc = wedge_h(node4, w);
                    let na = val[node1 - 1];
                    let nbv = val[node3 - 1];
                    val[*nb + k as usize - 1] = na
                        * nbv
                        * legendre_p(j as usize, lb - la)
                        * legendre_p((k - 1) as usize, lhc - lha);
                } else {
                    // Numbering permuted, switch indices j,k and nodes 2 and 4
                    let la = wedge_l(node1, u, v);
                    let lb = wedge_l(node4, u, v);
                    let lha = wedge_h(node1, w);
                    let lhc = wedge_h(node2, w);
                    let na = val[node1 - 1];
                    let nbv = val[node3 - 1];
                    val[*nb + k as usize - 1] = na
                        * nbv
                        * legendre_p((k - 1) as usize, lb - la)
                        * legendre_p(j as usize, lhc - lha);
                }
            }
            *nb += (pmax[i - 1] - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_dWedgeFaceP`: gradients of [`wedge_face_p`]. Reads the nodal
/// gradients the caller stored in `grad[0..6]`.
pub fn d_wedge_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let c = 1.0 / 2.0;
    let mut n = [0.0f64; 6];
    let mut nnb = 0;
    wedge_nodal_p(u, v, w, &mut n, &mut nnb);
    for i in 1..=2 {
        let fd = facedir[i - 1];
        let d_la = d_wedge_l(fd[0]);
        let d_lb = d_wedge_l(fd[1]);
        let d_lc = d_wedge_l(fd[2]);
        let d_na = d_wedge_h(fd[0]);
        let lav = wedge_l(fd[0], u, v);
        let lbv = wedge_l(fd[1], u, v);
        let lcv = wedge_l(fd[2], u, v);
        let nav = wedge_h(fd[0], w);
        for j in 0..=(pmax[i - 1] - 3) {
            for k in 1..=(pmax[i - 1] - j - 2) {
                let (la, lb, lc, na) = (lav, lbv, lcv, nav);
                let leg_p_lb_la = legendre_p(j as usize, lb - la);
                let leg_p2_lc1 = legendre_p((k - 1) as usize, 2.0 * lc - 1.0);
                let c_na = c * (1.0 + 2.0 * na);
                let idx = *nb + k as usize - 1;
                for q in 0..2 {
                    grad[idx][q] = c_na * d_la[q] * lb * lc * leg_p_lb_la * leg_p2_lc1
                        + c_na * la * d_lb[q] * lc * leg_p_lb_la * leg_p2_lc1
                        + c_na * la * lb * d_lc[q] * leg_p_lb_la * leg_p2_lc1
                        + c_na
                            * la
                            * lb
                            * lc
                            * d_legendre_p(j as usize, lb - la)
                            * (d_lb[q] - d_la[q])
                            * leg_p2_lc1
                        + c_na
                            * la
                            * lb
                            * lc
                            * leg_p_lb_la
                            * d_legendre_p((k - 1) as usize, 2.0 * lc - 1.0)
                            * (2.0 * d_lc[q]);
                }
                grad[idx][2] =
                    2.0 * c * d_na[2] * la * lb * lc * leg_p_lb_la * leg_p2_lc1;
            }
            *nb += (pmax[i - 1] - j - 2).max(0) as usize;
        }
    }
    for i in 3..=5 {
        let node1 = facedir[i - 1][0];
        let node2 = facedir[i - 1][1];
        let node3 = facedir[i - 1][2];
        let node4 = facedir[i - 1][3];
        let nonpermuted = square_face_nonpermuted(node1, node2);
        let (d_la, d_lb, d_lha, d_lhc);
        if nonpermuted {
            d_la = d_wedge_l(node1);
            d_lb = d_wedge_l(node2);
            d_lha = d_wedge_h(node1);
            d_lhc = d_wedge_h(node4);
        } else {
            d_la = d_wedge_l(node1);
            d_lb = d_wedge_l(node4);
            d_lha = d_wedge_h(node1);
            d_lhc = d_wedge_h(node2);
        }
        for j in 0..=(pmax[i - 1] - 2) {
            for k in 1..pmax[i - 1] {
                let (la, lb, lha, lhc, leg_j, leg_k, d_leg_j, d_leg_k);
                let na = n[node1 - 1];
                let nbv = n[node3 - 1];
                let d_na = grad[node1 - 1];
                let d_nb = grad[node3 - 1];
                if nonpermuted {
                    la = wedge_l(node1, u, v);
                    lb = wedge_l(node2, u, v);
                    lha = wedge_h(node1, w);
                    lhc = wedge_h(node4, w);
                    leg_j = legendre_p(j as usize, lb - la);
                    leg_k = legendre_p((k - 1) as usize, lhc - lha);
                    let sj = d_legendre_p(j as usize, lb - la);
                    let sk = d_legendre_p((k - 1) as usize, lhc - lha);
                    d_leg_j = [0, 1, 2].map(|q| sj * (d_lb[q] - d_la[q]));
                    d_leg_k = [0, 1, 2].map(|q| sk * (d_lhc[q] - d_lha[q]));
                } else {
                    la = wedge_l(node1, u, v);
                    lb = wedge_l(node4, u, v);
                    lha = wedge_h(node1, w);
                    lhc = wedge_h(node2, w);
                    leg_j = legendre_p((k - 1) as usize, lb - la);
                    leg_k = legendre_p(j as usize, lhc - lha);
                    let sj = d_legendre_p((k - 1) as usize, lb - la);
                    let sk = d_legendre_p(j as usize, lhc - lha);
                    d_leg_j = [0, 1, 2].map(|q| sj * (d_lb[q] - d_la[q]));
                    d_leg_k = [0, 1, 2].map(|q| sk * (d_lhc[q] - d_lha[q]));
                }
                let idx = *nb + k as usize - 1;
                for q in 0..3 {
                    grad[idx][q] = d_na[q] * nbv * leg_j * leg_k
                        + na * d_nb[q] * leg_j * leg_k
                        + na * nbv * d_leg_j[q] * leg_k
                        + na * nbv * leg_j * d_leg_k[q];
                }
            }
            *nb += (pmax[i - 1] - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_WedgeBubbleP`: element bubbles of degree `pmax`.
pub fn wedge_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    let l1v = wedge_l(1, u, v);
    let l2v = wedge_l(2, u, v);
    let l3v = wedge_l(3, u, v);
    for i in 0..=(pmax - 3) {
        for j in 0..=(pmax - 3 - i) {
            for k in 1..pmax {
                let (l1, l2, l3) = (l1v, l2v, l3v);
                let s = l2 - l1;
                let t = 2.0 * l3 - 1.0;
                val[*nb + k as usize - 1] = l1
                    * l2
                    * l3
                    * legendre_p(i as usize, s)
                    * legendre_p(j as usize, t)
                    * phi(k as usize + 1, w);
            }
            *nb += (pmax - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_dWedgeBubbleP`: gradients of [`wedge_bubble_p`].
pub fn d_wedge_bubble_p(u: f64, v: f64, w: f64, pmax: i32, grad: &mut [[f64; 3]], nb: &mut usize) {
    let d_l1 = d_wedge_l(1);
    let d_l2 = d_wedge_l(2);
    let d_l3 = d_wedge_l(3);
    let l1v = wedge_l(1, u, v);
    let l2v = wedge_l(2, u, v);
    let l3v = wedge_l(3, u, v);
    for i in 0..=(pmax - 3) {
        for j in 0..=(pmax - 3 - i) {
            for k in 1..pmax {
                let (l1, l2, l3) = (l1v, l2v, l3v);
                let s = l2 - l1;
                let t = 2.0 * l3 - 1.0;
                let ds = [0, 1, 2].map(|q| d_l2[q] - d_l1[q]);
                let dt = [0, 1, 2].map(|q| 2.0 * d_l3[q]);
                let legi = legendre_p(i as usize, s);
                let legj = legendre_p(j as usize, t);
                let legk = phi(k as usize + 1, w);
                let sd = d_legendre_p(i as usize, s);
                let td = d_legendre_p(j as usize, t);
                let wd = d_phi(k as usize + 1, w);
                let d_legi = [0, 1, 2].map(|q| sd * ds[q]);
                let d_legj = [0, 1, 2].map(|q| td * dt[q]);
                let d_legk = [wd * 0.0, wd * 0.0, wd * 1.0];
                let idx = *nb + k as usize - 1;
                for q in 0..3 {
                    grad[idx][q] = d_l1[q] * l2 * l3 * legi * legj * legk
                        + l1 * d_l2[q] * l3 * legi * legj * legk
                        + l1 * l2 * d_l3[q] * legi * legj * legk
                        + l1 * l2 * l3 * d_legi[q] * legj * legk
                        + l1 * l2 * l3 * legi * d_legj[q] * legk
                        + l1 * l2 * l3 * legi * legj * d_legk[q];
                }
            }
            *nb += (pmax - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_SD_WedgeEdgeP`: serendipity edge functions (nine edges).
pub fn sd_wedge_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0;
    for i in 1..=6 {
        let lav = wedge_l(edgedir[i - 1][0], u, v);
        let lbv = wedge_l(edgedir[i - 1][1], u, v);
        let nav = wedge_h(edgedir[i - 1][0], w);
        let nbv = wedge_h(edgedir[i - 1][1], w);
        for j in 1..pmax[i - 1] {
            let (la, lb, na, nbb) = (lav, lbv, nav, nbv);
            val[*nb + j as usize - 1] =
                c * la * lb * var_phi(j as usize + 1, lb - la) * (1.0 + na + nbb);
        }
        advance(nb, pmax[i - 1] - 1);
    }
    for i in 7..=9 {
        let lav = wedge_l(edgedir[i - 1][0], u, v);
        let nav = wedge_h(edgedir[i - 1][0], w);
        let nbv = wedge_h(edgedir[i - 1][1], w);
        for j in 1..pmax[i - 1] {
            let (la, na, nbb) = (lav, nav, nbv);
            val[*nb + j as usize - 1] = la * phi(j as usize + 1, nbb - na);
        }
        advance(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_SD_dWedgeEdgeP`: gradients of [`sd_wedge_edge_p`].
pub fn sd_d_wedge_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0;
    for i in 1..=6 {
        let d_la = d_wedge_l(edgedir[i - 1][0]);
        let d_lb = d_wedge_l(edgedir[i - 1][1]);
        let d_na = d_wedge_h(edgedir[i - 1][0]);
        let d_nb = d_wedge_h(edgedir[i - 1][1]);
        let lav = wedge_l(edgedir[i - 1][0], u, v);
        let lbv = wedge_l(edgedir[i - 1][1], u, v);
        let nav = wedge_h(edgedir[i - 1][0], w);
        let nbv = wedge_h(edgedir[i - 1][1], w);
        for j in 1..pmax[i - 1] {
            let (la, lb, na, nbb) = (lav, lbv, nav, nbv);
            let v_phi = var_phi(j as usize + 1, lb - la);
            let d_v_phi = d_var_phi(j as usize + 1, lb - la);
            let na_nb = 1.0 + na + nbb;
            let idx = *nb + j as usize - 1;
            for q in 0..2 {
                grad[idx][q] = c * d_la[q] * lb * v_phi * na_nb
                    + c * la * d_lb[q] * v_phi * na_nb
                    + c * la * lb * d_v_phi * (d_lb[q] - d_la[q]) * na_nb;
            }
            grad[idx][2] = c * la * lb * v_phi * (d_nb[2] + d_na[2]);
        }
        advance(nb, pmax[i - 1] - 1);
    }
    for i in 7..=9 {
        let d_la = d_wedge_l(edgedir[i - 1][0]);
        let d_na = d_wedge_h(edgedir[i - 1][0]);
        let d_nb = d_wedge_h(edgedir[i - 1][1]);
        let lav = wedge_l(edgedir[i - 1][0], u, v);
        let nav = wedge_h(edgedir[i - 1][0], w);
        let nbv = wedge_h(edgedir[i - 1][1], w);
        for j in 1..pmax[i - 1] {
            let (la, na, nbb) = (lav, nav, nbv);
            let ph = phi(j as usize + 1, nbb - na);
            let d_ph = d_phi(j as usize + 1, nbb - na);
            let idx = *nb + j as usize - 1;
            for q in 0..3 {
                grad[idx][q] = d_la[q] * ph + la * d_ph * (d_nb[q] - d_na[q]);
            }
        }
        advance(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_SD_WedgeFaceP`: serendipity face functions.
pub fn sd_wedge_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let c = 1.0 / 2.0;
    for i in 1..=2 {
        let fd = facedir[i - 1];
        let lav = wedge_l(fd[0], u, v);
        let lbv = wedge_l(fd[1], u, v);
        let lcv = wedge_l(fd[2], u, v);
        let nav = wedge_h(fd[0], w);
        for j in 0..=(pmax[i - 1] - 3) {
            for k in 1..=(pmax[i - 1] - j - 2) {
                let (la, lb, lc, na) = (lav, lbv, lcv, nav);
                val[*nb + k as usize - 1] = c
                    * (1.0 + 2.0 * na)
                    * la
                    * lb
                    * lc
                    * legendre_p(j as usize, lb - la)
                    * legendre_p((k - 1) as usize, 2.0 * lc - 1.0);
            }
            *nb += (pmax[i - 1] - j - 2).max(0) as usize;
        }
    }
    for i in 3..=5 {
        let fd = facedir[i - 1];
        let nonpermuted = square_face_nonpermuted(fd[0], fd[1]);
        for j in 2..=(pmax[i - 1] - 2) {
            for k in 1..=(pmax[i - 1] - j - 1) {
                if nonpermuted {
                    let la = wedge_l(fd[0], u, v);
                    let lb = wedge_l(fd[1], u, v);
                    let na = wedge_h(fd[0], w);
                    let nc = wedge_h(fd[3], w);
                    val[*nb + k as usize - 1] = la
                        * lb
                        * var_phi(j as usize, lb - la)
                        * phi(k as usize + 1, nc - na);
                } else {
                    // Numbering permuted, switch indices j,k and nodes 2 and 4
                    let la = wedge_l(fd[0], u, v);
                    let lb = wedge_l(fd[3], u, v);
                    let na = wedge_h(fd[0], w);
                    let nc = wedge_h(fd[1], w);
                    val[*nb + k as usize - 1] = la
                        * lb
                        * var_phi(k as usize + 1, lb - la)
                        * phi(j as usize, nc - na);
                }
            }
            *nb += (pmax[i - 1] - j - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_SD_dWedgeFaceP`: gradients of [`sd_wedge_face_p`].
pub fn sd_d_wedge_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let c = 1.0 / 2.0;
    for i in 1..=2 {
        let fd = facedir[i - 1];
        let d_la = d_wedge_l(fd[0]);
        let d_lb = d_wedge_l(fd[1]);
        let d_lc = d_wedge_l(fd[2]);
        let d_na = d_wedge_h(fd[0]);
        let lav = wedge_l(fd[0], u, v);
        let lbv = wedge_l(fd[1], u, v);
        let lcv = wedge_l(fd[2], u, v);
        let nav = wedge_h(fd[0], w);
        for j in 0..=(pmax[i - 1] - 3) {
            for k in 1..=(pmax[i - 1] - j - 2) {
                let (la, lb, lc, na) = (lav, lbv, lcv, nav);
                let leg_p_lb_la = legendre_p(j as usize, lb - la);
                let leg_p2_lc1 = legendre_p((k - 1) as usize, 2.0 * lc - 1.0);
                let c_na = c * (1.0 + 2.0 * na);
                let idx = *nb + k as usize - 1;
                for q in 0..2 {
                    grad[idx][q] = c_na * d_la[q] * lb * lc * leg_p_lb_la * leg_p2_lc1
                        + c_na * la * d_lb[q] * lc * leg_p_lb_la * leg_p2_lc1
                        + c_na * la * lb * d_lc[q] * leg_p_lb_la * leg_p2_lc1
                        + c_na
                            * la
                            * lb
                            * lc
                            * d_legendre_p(j as usize, lb - la)
                            * (d_lb[q] - d_la[q])
                            * leg_p2_lc1
                        + c_na
                            * la
                            * lb
                            * lc
                            * leg_p_lb_la
                            * d_legendre_p((k - 1) as usize, 2.0 * lc - 1.0)
                            * (2.0 * d_lc[q]);
                }
                grad[idx][2] =
                    2.0 * c * d_na[2] * la * lb * lc * leg_p_lb_la * leg_p2_lc1;
            }
            *nb += (pmax[i - 1] - j - 2).max(0) as usize;
        }
    }
    for i in 3..=5 {
        let fd = facedir[i - 1];
        let nonpermuted = square_face_nonpermuted(fd[0], fd[1]);
        let (d_la, d_lb, d_na, d_nc);
        if nonpermuted {
            d_la = d_wedge_l(fd[0]);
            d_lb = d_wedge_l(fd[1]);
            d_na = d_wedge_h(fd[0]);
            d_nc = d_wedge_h(fd[3]);
        } else {
            d_la = d_wedge_l(fd[0]);
            d_lb = d_wedge_l(fd[3]);
            d_na = d_wedge_h(fd[0]);
            d_nc = d_wedge_h(fd[1]);
        }
        for j in 2..=(pmax[i - 1] - 2) {
            for k in 1..=(pmax[i - 1] - j - 1) {
                let idx = *nb + k as usize - 1;
                if nonpermuted {
                    let la = wedge_l(fd[0], u, v);
                    let lb = wedge_l(fd[1], u, v);
                    let na = wedge_h(fd[0], w);
                    let nc = wedge_h(fd[3], w);
                    let v_ph = var_phi(j as usize, lb - la);
                    let ph = phi(k as usize + 1, nc - na);
                    for q in 0..2 {
                        grad[idx][q] = d_la[q] * lb * v_ph * ph
                            + la * d_lb[q] * v_ph * ph
                            + la * lb * d_var_phi(j as usize, lb - la) * (d_lb[q] - d_la[q]) * ph;
                    }
                    grad[idx][2] =
                        la * lb * v_ph * d_phi(k as usize + 1, nc - na) * (d_nc[2] - d_na[2]);
                } else {
                    let la = wedge_l(fd[0], u, v);
                    let lb = wedge_l(fd[3], u, v);
                    let na = wedge_h(fd[0], w);
                    let nc = wedge_h(fd[1], w);
                    let v_ph = var_phi(k as usize + 1, lb - la);
                    let ph = phi(j as usize, nc - na);
                    for q in 0..2 {
                        grad[idx][q] = d_la[q] * lb * v_ph * ph
                            + la * d_lb[q] * v_ph * ph
                            + la * lb
                                * d_var_phi(k as usize + 1, lb - la)
                                * (d_lb[q] - d_la[q])
                                * ph;
                    }
                    grad[idx][2] = la * lb * v_ph * d_phi(j as usize, nc - na) * (d_nc[2] - d_na[2]);
                }
            }
            *nb += (pmax[i - 1] - j - 1).max(0) as usize;
        }
    }
}

/// `H1Basis_SD_WedgeBubbleP`: serendipity element bubbles.
pub fn sd_wedge_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    let l1v = wedge_l(1, u, v);
    let l2v = wedge_l(2, u, v);
    let l3v = wedge_l(3, u, v);
    for i in 0..=(pmax - 5) {
        for j in 0..=(pmax - 5 - i) {
            for k in 1..=(pmax - 4 - i - j) {
                let (l1, l2, l3) = (l1v, l2v, l3v);
                let l2_l1 = l2 - l1;
                let l3_1 = 2.0 * l3 - 1.0;
                val[*nb + k as usize - 1] = l1
                    * l2
                    * l3
                    * legendre_p(i as usize, l2_l1)
                    * legendre_p(j as usize, l3_1)
                    * phi(k as usize + 1, w);
            }
            *nb += (pmax - 4 - i - j).max(0) as usize;
        }
    }
}

/// `H1Basis_SD_dWedgeBubbleP`: gradients of [`sd_wedge_bubble_p`].
pub fn sd_d_wedge_bubble_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
) {
    let sq3_6 = 3.0f64.sqrt() / 6.0;
    let sq3_3 = 3.0f64.sqrt() / 3.0;
    let l1v = wedge_l(1, u, v);
    let l2v = wedge_l(2, u, v);
    let l3v = wedge_l(3, u, v);
    for i in 0..=(pmax - 5) {
        for j in 0..=(pmax - 5 - i) {
            for k in 1..=(pmax - 4 - i - j) {
                let (l1, l2, l3) = (l1v, l2v, l3v);
                let l2_l1 = l2 - l1;
                let l3_1 = 2.0 * l3 - 1.0;
                let legi = legendre_p(i as usize, l2_l1);
                let legj = legendre_p(j as usize, l3_1);
                let phi_w = phi(k as usize + 1, w);
                let idx = *nb + k as usize - 1;
                grad[idx][0] = ((-1.0 / 2.0) * l2 * l3 * legi * legj
                    + l1 * (1.0 / 2.0) * l3 * legi * legj
                    + l1 * l2 * l3 * d_legendre_p(i as usize, l2_l1) * legj)
                    * phi_w;
                grad[idx][1] = ((-sq3_6) * l2 * l3 * legi * legj
                    + l1 * (-sq3_6) * l3 * legi * legj
                    + l1 * l2 * sq3_3 * legi * legj
                    + l1 * l2 * l3 * legi * d_legendre_p(j as usize, l3_1) * (2.0 * sq3_3))
                    * phi_w;
                grad[idx][2] = l1 * l2 * l3 * legi * legj * d_phi(k as usize + 1, w);
            }
            *nb += (pmax - 4 - i - j).max(0) as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EDGE_MAP: [[usize; 2]; 9] = [
        [1, 2],
        [2, 3],
        [3, 1],
        [4, 5],
        [5, 6],
        [6, 4],
        [1, 4],
        [2, 5],
        [3, 6],
    ];
    const FACE_MAP: [[usize; 4]; 5] = [
        [1, 2, 3, 0],
        [4, 5, 6, 0],
        [1, 2, 5, 4],
        [2, 3, 6, 5],
        [3, 1, 4, 6],
    ];

    /// Interior points of the reference wedge (u,w in [-1,1], v in [0,sqrt3]).
    const POINTS: [[f64; 3]; 4] = [
        [0.1, 0.6, -0.3],
        [-0.35, 0.45, 0.55],
        [0.0, 0.9, 0.1],
        [0.4, 0.3, -0.7],
    ];

    fn edge_dirs(swap: bool) -> Vec<[usize; 2]> {
        let mut d = EDGE_MAP.to_vec();
        if swap {
            for i in [0usize, 4, 7] {
                d[i].swap(0, 1);
            }
        }
        d
    }

    fn face_dirs(variant: u8) -> Vec<[usize; 4]> {
        let mut d = FACE_MAP.to_vec();
        match variant {
            0 => {}
            1 => {
                // triangle faces with reversed first two nodes, permuted square faces
                d[0] = [2, 1, 3, 0];
                d[1] = [6, 4, 5, 0];
                d[2] = [1, 4, 5, 2];
                d[3] = [3, 2, 5, 6];
                d[4] = [3, 6, 4, 1];
            }
            _ => unreachable!(),
        }
        d
    }

    #[derive(Clone, Copy, PartialEq, Debug)]
    enum Kind {
        Edge,
        Face,
        Bubble,
    }

    /// Values of nodal functions followed by the `kind` functions.
    fn values(kind: Kind, sd: bool, p: i32, e: &[[usize; 2]], f: &[[usize; 4]], x: [f64; 3]) -> Vec<f64> {
        let mut val = vec![0.0; 400];
        let mut nb = 0;
        wedge_nodal_p(x[0], x[1], x[2], &mut val, &mut nb);
        let pm = [p; 9];
        let n0 = nb;
        match (kind, sd) {
            (Kind::Edge, false) => wedge_edge_p(x[0], x[1], x[2], &pm, &mut val, &mut nb, e),
            (Kind::Edge, true) => sd_wedge_edge_p(x[0], x[1], x[2], &pm, &mut val, &mut nb, e),
            (Kind::Face, false) => wedge_face_p(x[0], x[1], x[2], &pm, &mut val, &mut nb, f),
            (Kind::Face, true) => sd_wedge_face_p(x[0], x[1], x[2], &pm, &mut val, &mut nb, f),
            (Kind::Bubble, false) => wedge_bubble_p(x[0], x[1], x[2], p, &mut val, &mut nb),
            (Kind::Bubble, true) => sd_wedge_bubble_p(x[0], x[1], x[2], p, &mut val, &mut nb),
        }
        val.truncate(nb);
        val.drain(..n0);
        val
    }

    fn grads(kind: Kind, sd: bool, p: i32, e: &[[usize; 2]], f: &[[usize; 4]], x: [f64; 3]) -> Vec<[f64; 3]> {
        let mut g = vec![[0.0; 3]; 400];
        let mut nb = 0;
        d_wedge_nodal_p(x[0], x[1], x[2], &mut g, &mut nb);
        let pm = [p; 9];
        let n0 = nb;
        match (kind, sd) {
            (Kind::Edge, false) => d_wedge_edge_p(x[0], x[1], x[2], &pm, &mut g, &mut nb, e),
            (Kind::Edge, true) => sd_d_wedge_edge_p(x[0], x[1], x[2], &pm, &mut g, &mut nb, e),
            (Kind::Face, false) => d_wedge_face_p(x[0], x[1], x[2], &pm, &mut g, &mut nb, f),
            (Kind::Face, true) => sd_d_wedge_face_p(x[0], x[1], x[2], &pm, &mut g, &mut nb, f),
            (Kind::Bubble, false) => d_wedge_bubble_p(x[0], x[1], x[2], p, &mut g, &mut nb),
            (Kind::Bubble, true) => sd_d_wedge_bubble_p(x[0], x[1], x[2], p, &mut g, &mut nb),
        }
        g.truncate(nb);
        g.drain(..n0);
        g
    }

    fn fd_check(kind: Kind, sd: bool, p: i32, e: &[[usize; 2]], f: &[[usize; 4]]) -> usize {
        let h = 1e-6;
        let mut total = 0;
        for x in POINTS {
            let g = grads(kind, sd, p, e, f, x);
            let v0 = values(kind, sd, p, e, f, x);
            assert_eq!(g.len(), v0.len(), "{kind:?} sd={sd} p={p}: value/gradient counts");
            total = v0.len();
            for c in 0..3 {
                let mut xp = x;
                let mut xm = x;
                xp[c] += h;
                xm[c] -= h;
                let vp = values(kind, sd, p, e, f, xp);
                let vm = values(kind, sd, p, e, f, xm);
                for i in 0..vp.len() {
                    let fd = (vp[i] - vm[i]) / (2.0 * h);
                    let tol = 1e-6 * (1.0 + g[i][c].abs());
                    assert!(
                        (fd - g[i][c]).abs() < tol,
                        "{kind:?} sd={sd} p={p} fn {i} comp {c} at {x:?}: fd {fd} vs {}",
                        g[i][c]
                    );
                }
            }
        }
        total
    }

    #[test]
    fn nodal_partition_of_unity() {
        for x in POINTS {
            let mut v = [0.0; 6];
            let mut nb = 0;
            wedge_nodal_p(x[0], x[1], x[2], &mut v, &mut nb);
            assert_eq!(nb, 6);
            assert!((v.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let mut g = [[0.0; 3]; 6];
            let mut nb = 0;
            d_wedge_nodal_p(x[0], x[1], x[2], &mut g, &mut nb);
            assert_eq!(nb, 6);
            for c in 0..3 {
                assert!(g.iter().map(|r| r[c]).sum::<f64>().abs() < 1e-14);
            }
            // Linear unit-triangle wedge.
            let (u, w) = (0.2, x[2]);
            let v1 = 0.3;
            wedge_nodal(u, v1, w, &mut v);
            assert!((v.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let mut g = [[0.0; 3]; 6];
            d_wedge_nodal(u, v1, w, &mut g);
            for c in 0..3 {
                assert!(g.iter().map(|r| r[c]).sum::<f64>().abs() < 1e-14);
            }
        }
    }

    #[test]
    fn linear_nodal_is_p_nodal_in_other_coordinates() {
        // wedge_nodal and wedge_nodal_p describe different triangles but must
        // each be the shape functions of their own coordinates: check FD.
        let (u, v, w) = (0.2, 0.3, -0.4);
        let mut g = [[0.0; 3]; 6];
        d_wedge_nodal(u, v, w, &mut g);
        let h = 1e-6;
        for c in 0..3 {
            let mut xp = [u, v, w];
            let mut xm = [u, v, w];
            xp[c] += h;
            xm[c] -= h;
            let (mut a, mut b) = ([0.0; 6], [0.0; 6]);
            wedge_nodal(xp[0], xp[1], xp[2], &mut a);
            wedge_nodal(xm[0], xm[1], xm[2], &mut b);
            for i in 0..6 {
                assert!(((a[i] - b[i]) / (2.0 * h) - g[i][c]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn nodal_p_gradient_matches_fd_and_edge_functions_vanish_on_other_nodes() {
        let h = 1e-6;
        let x = POINTS[0];
        let mut g = [[0.0; 3]; 6];
        let mut nb = 0;
        d_wedge_nodal_p(x[0], x[1], x[2], &mut g, &mut nb);
        for c in 0..3 {
            let (mut xp, mut xm) = (x, x);
            xp[c] += h;
            xm[c] -= h;
            let (mut a, mut b) = ([0.0; 6], [0.0; 6]);
            let mut n = 0;
            wedge_nodal_p(xp[0], xp[1], xp[2], &mut a, &mut n);
            n = 0;
            wedge_nodal_p(xm[0], xm[1], xm[2], &mut b, &mut n);
            for i in 0..6 {
                assert!(((a[i] - b[i]) / (2.0 * h) - g[i][c]).abs() < 1e-8);
            }
        }
        // Nodal functions are 1 at their own node (node 1 at (-1,0,-1)).
        let mut v = [0.0; 6];
        let mut n = 0;
        wedge_nodal_p(-1.0, 0.0, -1.0, &mut v, &mut n);
        assert!((v[0] - 1.0).abs() < 1e-14 && v[1..].iter().all(|a| a.abs() < 1e-14));
    }

    #[test]
    fn full_and_serendipity_gradients_match_finite_differences() {
        for sd in [false, true] {
            for p in 2..=6 {
                if !sd && p > 5 {
                    continue;
                }
                for swap in [false, true] {
                    let e = edge_dirs(swap);
                    let f = face_dirs(swap as u8);
                    for kind in [Kind::Edge, Kind::Face, Kind::Bubble] {
                        fd_check(kind, sd, p, &e, &f);
                    }
                }
            }
        }
    }

    #[test]
    fn counts_match_element_info_formulas() {
        for sd in [false, true] {
            for p in 2..=7i32 {
                let e = edge_dirs(false);
                let f = face_dirs(0);
                // edges: SUM(EdgeDegree-1)
                let ne = values(Kind::Edge, sd, p, &e, &f, POINTS[0]).len() as i32;
                assert_eq!(ne, 9 * (p - 1));
                // faces: ElementInfoVec_ComputePElementBasis (706)
                let mut nbq = 0;
                for _ in 1..=2 {
                    for j in 0..=(p - 3) {
                        nbq += (p - j - 2).max(0);
                    }
                }
                for _ in 3..=5 {
                    if sd {
                        for j in 2..=(p - 2) {
                            nbq += (p - j - 1).max(0);
                        }
                    } else {
                        for _ in 0..=(p - 2) {
                            nbq += (p - 1).max(0);
                        }
                    }
                }
                let nf = values(Kind::Face, sd, p, &e, &f, POINTS[0]).len() as i32;
                assert_eq!(nf, nbq, "face count sd={sd} p={p}");
                // bubbles: getBubbleDOFs
                let bd = if sd {
                    ((p - 2) * (p - 3) * (p - 4) / 6).max(0)
                } else {
                    (p - 1) * (p - 1) * (p - 2) / 2
                };
                let nbub = values(Kind::Bubble, sd, p, &e, &f, POINTS[0]).len() as i32;
                assert_eq!(nbub, bd, "bubble count sd={sd} p={p}");
                let ng = grads(Kind::Bubble, sd, p, &e, &f, POINTS[0]).len() as i32;
                assert_eq!(ng, bd);
            }
        }
    }

    #[test]
    fn edge_functions_vanish_off_their_edge() {
        // The bottom-triangle edge-1 function (nodes 1-2) vanishes on the
        // opposite vertical edge (node 3 column, v = sqrt3).
        let e = edge_dirs(false);
        let f = face_dirs(0);
        for p in 2..=4 {
            let vals = values(Kind::Edge, false, p, &e, &f, [-1.0, 0.0, 0.0]);
            // At the left vertical edge at (u,v)=(-1,0), every function that
            // is not attached to nodes 1 or 4 must vanish: edges 2 (2-3), 5, 8, 9.
            let per = (p - 1) as usize;
            for edge in [1usize, 4, 7, 8] {
                for j in 0..per {
                    assert!(vals[edge * per + j].abs() < 1e-13, "edge {edge} p={p}");
                }
            }
        }
    }
}
