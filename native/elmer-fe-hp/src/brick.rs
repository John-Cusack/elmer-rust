//! p-hierarchic H1 basis on the brick (element 808).
//!
//! Scalar transliteration of the brick routines of `fem/src/H1Basis.F90`
//! (`H1Basis_BrickNodal` .. `H1Basis_dBrickBubbleP`, including the serendipity
//! `H1Basis_SD_*` variants). See [`crate::types`] for the shared conventions:
//! 1-based local node numbers in `edgedir` / `facedir`, a running `nb` count
//! of stored functions, and caller-provided nodal values / gradients in the
//! first eight slots of `val` / `grad` for the edge and face routines that
//! read them (`fval(k,node)` / `grad(k,node,:)` in the Fortran).

#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::many_single_char_names
)]

use crate::legendre::{d_legendre_p, d_phi, d_var_phi, legendre_p, phi, var_phi};

/// Fortran `nbasis = nbasis + pmax - 1` where `pmax` may be zero or negative.
#[inline]
fn adv(nb: &mut usize, count: i32) {
    *nb = (*nb as isize + count as isize) as usize;
}

/// Fortran `nbasis + k` (`k` 1-based) as a 0-based slot index.
#[inline]
fn slot(nb: usize, k: i32) -> usize {
    (nb as isize + k as isize - 1) as usize
}

/// `H1Basis_BrickNodal`: the eight trilinear nodal functions.
pub fn brick_nodal(u: f64, v: f64, w: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 8.0;
    let b = *nb;
    val[b] = c * (1.0 - u) * (1.0 - v) * (1.0 - w);
    val[b + 1] = c * (1.0 + u) * (1.0 - v) * (1.0 - w);
    val[b + 2] = c * (1.0 + u) * (1.0 + v) * (1.0 - w);
    val[b + 3] = c * (1.0 - u) * (1.0 + v) * (1.0 - w);
    val[b + 4] = c * (1.0 - u) * (1.0 - v) * (1.0 + w);
    val[b + 5] = c * (1.0 + u) * (1.0 - v) * (1.0 + w);
    val[b + 6] = c * (1.0 + u) * (1.0 + v) * (1.0 + w);
    val[b + 7] = c * (1.0 - u) * (1.0 + v) * (1.0 + w);
    *nb += 8;
}

/// `H1Basis_dBrickNodal`: local gradients of the eight nodal functions.
pub fn d_brick_nodal(u: f64, v: f64, w: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let c = 1.0 / 8.0;
    let b = *nb;
    // First coordinate (xi)
    grad[b][0] = -c * (1.0 - v) * (1.0 - w);
    grad[b + 1][0] = c * (1.0 - v) * (1.0 - w);
    grad[b + 2][0] = c * (1.0 + v) * (1.0 - w);
    grad[b + 3][0] = -c * (1.0 + v) * (1.0 - w);
    grad[b + 4][0] = -c * (1.0 - v) * (1.0 + w);
    grad[b + 5][0] = c * (1.0 - v) * (1.0 + w);
    grad[b + 6][0] = c * (1.0 + v) * (1.0 + w);
    grad[b + 7][0] = -c * (1.0 + v) * (1.0 + w);
    // Second coordinate (eta)
    grad[b][1] = -c * (1.0 - u) * (1.0 - w);
    grad[b + 1][1] = -c * (1.0 + u) * (1.0 - w);
    grad[b + 2][1] = c * (1.0 + u) * (1.0 - w);
    grad[b + 3][1] = c * (1.0 - u) * (1.0 - w);
    grad[b + 4][1] = -c * (1.0 - u) * (1.0 + w);
    grad[b + 5][1] = -c * (1.0 + u) * (1.0 + w);
    grad[b + 6][1] = c * (1.0 + u) * (1.0 + w);
    grad[b + 7][1] = c * (1.0 - u) * (1.0 + w);
    // Third coordinate (zeta)
    grad[b][2] = -c * (1.0 - u) * (1.0 - v);
    grad[b + 1][2] = -c * (1.0 + u) * (1.0 - v);
    grad[b + 2][2] = -c * (1.0 + u) * (1.0 + v);
    grad[b + 3][2] = -c * (1.0 - u) * (1.0 + v);
    grad[b + 4][2] = c * (1.0 - u) * (1.0 - v);
    grad[b + 5][2] = c * (1.0 + u) * (1.0 - v);
    grad[b + 6][2] = c * (1.0 + u) * (1.0 + v);
    grad[b + 7][2] = c * (1.0 - u) * (1.0 + v);
    *nb += 8;
}

/// `H1Basis_BrickL(node, u, v, w)`: affine "lambda" of local node `node` (1..=8).
pub fn brick_l(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let c = 1.0 / 2.0;
    match node {
        1 => c * (3.0 - u - v - w),
        2 => c * (3.0 + u - v - w),
        3 => c * (3.0 + u + v - w),
        4 => c * (3.0 - u + v - w),
        5 => c * (3.0 - u - v + w),
        6 => c * (3.0 + u - v + w),
        7 => c * (3.0 + u + v + w),
        8 => c * (3.0 - u + v + w),
        _ => panic!("brick_l: node {node} outside 1..=8"),
    }
}

/// `H1Basis_dBrickL(node)`: gradient of [`brick_l`].
pub fn d_brick_l(node: usize) -> [f64; 3] {
    // c = 1/2 (exact), so c*(+-1) is +-0.5.
    match node {
        1 => [-0.5, -0.5, -0.5],
        2 => [0.5, -0.5, -0.5],
        3 => [0.5, 0.5, -0.5],
        4 => [-0.5, 0.5, -0.5],
        5 => [-0.5, -0.5, 0.5],
        6 => [0.5, -0.5, 0.5],
        7 => [0.5, 0.5, 0.5],
        8 => [-0.5, 0.5, 0.5],
        _ => panic!("d_brick_l: node {node} outside 1..=8"),
    }
}

/// `H1Basis_BrickEdgeL(edge, u, v, w, La, Lb)`: the two affine factors `(La, Lb)`
/// that vanish on the faces adjacent to edge `edge` (1..=12).
pub fn brick_edge_l(edge: usize, u: f64, v: f64, w: f64) -> (f64, f64) {
    match edge {
        1 => (1.0 - v, 1.0 - w),
        2 => (1.0 + u, 1.0 - w),
        3 => (1.0 + v, 1.0 - w),
        4 => (1.0 - u, 1.0 - w),
        5 => (1.0 - v, 1.0 + w),
        6 => (1.0 + u, 1.0 + w),
        7 => (1.0 + v, 1.0 + w),
        8 => (1.0 - u, 1.0 + w),
        9 => (1.0 - u, 1.0 - v),
        10 => (1.0 + u, 1.0 - v),
        11 => (1.0 + u, 1.0 + v),
        12 => (1.0 - u, 1.0 + v),
        _ => panic!("brick_edge_l: edge {edge} outside 1..=12"),
    }
}

/// `H1Basis_dBrickEdgeL(edge, dLa, dLb)`: gradients `(dLa, dLb)` of [`brick_edge_l`].
pub fn d_brick_edge_l(edge: usize) -> ([f64; 3], [f64; 3]) {
    match edge {
        1 => ([0.0, -1.0, 0.0], [0.0, 0.0, -1.0]),
        2 => ([1.0, 0.0, 0.0], [0.0, 0.0, -1.0]),
        3 => ([0.0, 1.0, 0.0], [0.0, 0.0, -1.0]),
        4 => ([-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]),
        5 => ([0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
        6 => ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        7 => ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        8 => ([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        9 => ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]),
        10 => ([1.0, 0.0, 0.0], [0.0, -1.0, 0.0]),
        11 => ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        12 => ([-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        _ => panic!("d_brick_edge_l: edge {edge} outside 1..=12"),
    }
}

// --- start serendipity brick

/// `H1Basis_SD_BrickEdgeP`: serendipity edge functions, `pmax[i]-1` per edge.
pub fn sd_brick_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 4.0;
    // For each edge
    for i in 1..=12usize {
        let la = brick_l(edgedir[i - 1][0], u, v, w);
        let lb = brick_l(edgedir[i - 1][1], u, v, w);
        let (aa, ba) = brick_edge_l(i, u, v, w);
        for j in 1..=pmax[i - 1] - 1 {
            val[slot(*nb, j)] = c * phi((j + 1) as usize, lb - la) * aa * ba;
        }
        adv(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_SD_dBrickEdgeP`: gradients of [`sd_brick_edge_p`].
pub fn d_sd_brick_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 4.0;
    // For each edge
    for i in 1..=12usize {
        let dla = d_brick_l(edgedir[i - 1][0]);
        let dlb = d_brick_l(edgedir[i - 1][1]);
        let (daa, dba) = d_brick_edge_l(i);
        let la = brick_l(edgedir[i - 1][0], u, v, w);
        let lb = brick_l(edgedir[i - 1][1], u, v, w);
        let (aa, ba) = brick_edge_l(i, u, v, w);
        for j in 1..=pmax[i - 1] - 1 {
            let p = phi((j + 1) as usize, lb - la);
            let dp = d_phi((j + 1) as usize, lb - la);
            let g = &mut grad[slot(*nb, j)];
            g[0] = c * dp * (dlb[0] - dla[0]) * aa * ba + c * p * daa[0] * ba + c * p * aa * dba[0];
            g[1] = c * dp * (dlb[1] - dla[1]) * aa * ba + c * p * daa[1] * ba + c * p * aa * dba[1];
            g[2] = c * dp * (dlb[2] - dla[2]) * aa * ba + c * p * daa[2] * ba + c * p * aa * dba[2];
        }
        adv(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_SD_BrickFaceP`: serendipity face functions of face degrees `pmax[i]`.
pub fn sd_brick_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let a = 1.0 / 4.0;
    // For each face
    for i in 1..=6usize {
        let la = brick_l(facedir[i - 1][0], u, v, w);
        let lb = brick_l(facedir[i - 1][1], u, v, w);
        let lc = brick_l(facedir[i - 1][2], u, v, w);
        let ld = brick_l(facedir[i - 1][3], u, v, w);
        let pm = pmax[i - 1];
        for j in 2..=pm {
            for k in 1..=pm - j - 1 {
                val[slot(*nb, k)] = (a * (la + lb + lc + ld) - 1.0)
                    * phi(j as usize, lb - la)
                    * phi((k + 1) as usize, ld - la);
            }
            adv(nb, (pm - j - 1).max(0));
        }
    }
}

/// `H1Basis_SD_dBrickFaceP`: gradients of [`sd_brick_face_p`].
pub fn d_sd_brick_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let a = 1.0 / 4.0;
    // For each face
    for i in 1..=6usize {
        let dla = d_brick_l(facedir[i - 1][0]);
        let dlb = d_brick_l(facedir[i - 1][1]);
        let dlc = d_brick_l(facedir[i - 1][2]);
        let dld = d_brick_l(facedir[i - 1][3]);
        let la = brick_l(facedir[i - 1][0], u, v, w);
        let lb = brick_l(facedir[i - 1][1], u, v, w);
        let lc = brick_l(facedir[i - 1][2], u, v, w);
        let ld = brick_l(facedir[i - 1][3], u, v, w);
        let pm = pmax[i - 1];
        for j in 2..=pm {
            for k in 1..=pm - j - 1 {
                let phi_lb_la = phi(j as usize, lb - la);
                let phi_ld_la = phi((k + 1) as usize, ld - la);
                let la_lb_lc_ld = a * (la + lb + lc + ld) - 1.0;
                let g = &mut grad[slot(*nb, k)];
                for c in 0..3 {
                    g[c] = a * (dla[c] + dlb[c] + dlc[c] + dld[c]) * phi_lb_la * phi_ld_la
                        + la_lb_lc_ld * d_phi(j as usize, lb - la) * (dlb[c] - dla[c]) * phi_ld_la
                        + la_lb_lc_ld
                            * phi_lb_la
                            * d_phi((k + 1) as usize, ld - la)
                            * (dld[c] - dla[c]);
                }
            }
            adv(nb, (pm - j - 1).max(0));
        }
    }
}

/// `H1Basis_SD_BrickBubbleP`: serendipity bubble functions of degree `pmax`.
pub fn sd_brick_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    // NOTE (Fortran): getBubbleDOFs is the volume of these loops.
    for i in 2..=pmax - 4 {
        for j in 2..=pmax - i - 2 {
            for k in 1..=pmax - i - j - 1 {
                val[slot(*nb, k)] =
                    phi(i as usize, u) * phi(j as usize, v) * phi((k + 1) as usize, w);
            }
            adv(nb, (pmax - i - j - 1).max(0));
        }
    }
}

/// `H1Basis_SD_dBrickBubbleP`: gradients of [`sd_brick_bubble_p`].
pub fn d_sd_brick_bubble_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
) {
    for i in 2..=pmax - 4 {
        for j in 2..=pmax - i - 2 {
            for k in 1..=pmax - i - j - 1 {
                let phi_u = phi(i as usize, u);
                let phi_v = phi(j as usize, v);
                let phi_w = phi((k + 1) as usize, w);
                let g = &mut grad[slot(*nb, k)];
                g[0] = d_phi(i as usize, u) * phi_v * phi_w;
                g[1] = phi_u * d_phi(j as usize, v) * phi_w;
                g[2] = phi_u * phi_v * d_phi((k + 1) as usize, w);
            }
            adv(nb, (pmax - i - j - 1).max(0));
        }
    }
}

// --- end serendipity brick

/// `H1Basis_BrickEdgeP`: full-polynomial edge functions, `pmax[i]-1` per edge.
///
/// Reads the nodal values stored in `val[0..8]` by [`brick_nodal`] (`fval(k,node)`).
pub fn brick_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let mut lp = [0.0; 8];
    for i in 1..=8usize {
        lp[i - 1] = brick_l(i, u, v, w);
    }
    // For each edge
    for i in 1..=12usize {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        for j in 1..=pmax[i - 1] - 1 {
            let na = val[node1 - 1];
            let nbv = val[node2 - 1];
            let la = lp[node1 - 1];
            let lb = lp[node2 - 1];
            val[slot(*nb, j)] = na * nbv * var_phi((j + 1) as usize, lb - la);
        }
        adv(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_dBrickEdgeP`: gradients of [`brick_edge_p`].
///
/// Recomputes the nodal values internally (the Fortran local `N`) and reads the
/// nodal gradients stored in `grad[0..8]` by [`d_brick_nodal`].
pub fn d_brick_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let mut n = [0.0; 8];
    let mut nnb = 0usize;
    brick_nodal(u, v, w, &mut n, &mut nnb);
    let mut lp = [0.0; 8];
    for i in 1..=8usize {
        lp[i - 1] = brick_l(i, u, v, w);
    }
    // For each edge
    for i in 1..=12usize {
        let node1 = edgedir[i - 1][0];
        let node2 = edgedir[i - 1][1];
        let dla = d_brick_l(node1);
        let dlb = d_brick_l(node2);
        for j in 1..=pmax[i - 1] - 1 {
            let la = lp[node1 - 1];
            let lb = lp[node2 - 1];
            let na = n[node1 - 1];
            let nbv = n[node2 - 1];
            let dna = grad[node1 - 1];
            let dnb = grad[node2 - 1];
            let p = var_phi((j + 1) as usize, lb - la);
            let dp = d_var_phi((j + 1) as usize, lb - la);
            let g = &mut grad[slot(*nb, j)];
            for c in 0..3 {
                g[c] = dp * (dlb[c] - dla[c]) * na * nbv + p * dna[c] * nbv + p * na * dnb[c];
            }
        }
        adv(nb, pmax[i - 1] - 1);
    }
}

/// `H1Basis_BrickFaceP`: full-polynomial face functions, `(pmax[i]-1)^2` per face.
///
/// Reads the nodal values stored in `val[0..8]` by [`brick_nodal`].
pub fn brick_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let mut lp = [0.0; 8];
    for i in 1..=8usize {
        lp[i - 1] = brick_l(i, u, v, w);
    }
    // For each face
    for i in 1..=6usize {
        let node1 = facedir[i - 1][0];
        let node2 = facedir[i - 1][1];
        let node3 = facedir[i - 1][2];
        let node4 = facedir[i - 1][3];
        let pm = pmax[i - 1];
        for j in 0..=pm - 2 {
            for k in 1..=pm - 1 {
                let na = val[node1 - 1];
                let nbv = val[node3 - 1];
                let la = lp[node1 - 1];
                let lb = lp[node2 - 1];
                let ld = lp[node4 - 1];
                val[slot(*nb, k)] = na
                    * nbv
                    * legendre_p(j as usize, lb - la)
                    * legendre_p((k - 1) as usize, ld - la);
            }
            adv(nb, pm - 1);
        }
    }
}

/// `H1Basis_dBrickFaceP`: gradients of [`brick_face_p`].
///
/// Recomputes the nodal values internally and reads the nodal gradients stored
/// in `grad[0..8]` by [`d_brick_nodal`].
pub fn d_brick_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let mut n = [0.0; 8];
    let mut nnb = 0usize;
    brick_nodal(u, v, w, &mut n, &mut nnb);
    let mut lp = [0.0; 8];
    for i in 1..=8usize {
        lp[i - 1] = brick_l(i, u, v, w);
    }
    // For each face
    for i in 1..=6usize {
        let node1 = facedir[i - 1][0];
        let node2 = facedir[i - 1][1];
        let node3 = facedir[i - 1][2];
        let node4 = facedir[i - 1][3];
        let dla = d_brick_l(node1);
        let dlb = d_brick_l(node2);
        // Fortran also evaluates dLc = dBrickL(node3), which is not used.
        let dld = d_brick_l(node4);
        let pm = pmax[i - 1];
        for j in 0..=pm - 2 {
            for k in 1..=pm - 1 {
                let la = lp[node1 - 1];
                let lb = lp[node2 - 1];
                let ld = lp[node4 - 1];
                let phi_u = legendre_p(j as usize, lb - la);
                let phi_v = legendre_p((k - 1) as usize, ld - la);
                let dphi_u = d_legendre_p(j as usize, lb - la);
                let dphi_v = d_legendre_p((k - 1) as usize, ld - la);
                let na = n[node1 - 1];
                let nbv = n[node3 - 1];
                let dna = grad[node1 - 1];
                let dnb = grad[node3 - 1];
                let g = &mut grad[slot(*nb, k)];
                for c in 0..3 {
                    g[c] = dna[c] * nbv * phi_u * phi_v
                        + na * dnb[c] * phi_u * phi_v
                        + na * nbv * dphi_u * (dlb[c] - dla[c]) * phi_v
                        + na * nbv * phi_u * dphi_v * (dld[c] - dla[c]);
                }
            }
            adv(nb, pm - 1);
        }
    }
}

/// `H1Basis_BrickBubbleP`: full-polynomial bubble functions, `(pmax-1)^3` of them.
pub fn brick_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    for i in 2..=pmax {
        for j in 2..=pmax {
            for k in 1..=pmax - 1 {
                val[slot(*nb, k)] =
                    phi(i as usize, u) * phi(j as usize, v) * phi((k + 1) as usize, w);
            }
            adv(nb, pmax - 1);
        }
    }
}

/// `H1Basis_dBrickBubbleP`: gradients of [`brick_bubble_p`].
pub fn d_brick_bubble_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
) {
    for i in 2..=pmax {
        for j in 2..=pmax {
            for k in 1..=pmax - 1 {
                let phi_u = phi(i as usize, u);
                let phi_v = phi(j as usize, v);
                let phi_w = phi((k + 1) as usize, w);
                let g = &mut grad[slot(*nb, k)];
                g[0] = d_phi(i as usize, u) * phi_v * phi_w;
                g[1] = phi_u * d_phi(j as usize, v) * phi_w;
                g[2] = phi_u * phi_v * d_phi((k + 1) as usize, w);
            }
            adv(nb, pmax - 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `H1Basis_GetEdgeMap(808)`.
    const EDGE_MAP: [[usize; 2]; 12] = [
        [1, 2],
        [2, 3],
        [4, 3],
        [1, 4],
        [5, 6],
        [6, 7],
        [8, 7],
        [5, 8],
        [1, 5],
        [2, 6],
        [3, 7],
        [4, 8],
    ];
    /// `H1Basis_GetFaceMap(808)`.
    const FACE_MAP: [[usize; 4]; 6] = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [4, 3, 7, 8],
        [1, 4, 8, 5],
    ];

    /// `H1Basis_GetEdgeDirection` for global node indices `g` (1 per local node).
    fn edge_direction(g: &[i32; 8]) -> [[usize; 2]; 12] {
        let mut d = EDGE_MAP;
        for e in d.iter_mut() {
            if g[e[1] - 1] < g[e[0] - 1] {
                e.swap(0, 1);
            }
        }
        d
    }

    /// `H1Basis_GetFaceDirection` (square faces) for global node indices `g`.
    fn face_direction(g: &[i32; 8]) -> [[usize; 4]; 6] {
        let mut d = FACE_MAP;
        for f in d.iter_mut() {
            let mut min_i = 0;
            for i in 1..4 {
                if g[f[i] - 1] < g[f[min_i] - 1] {
                    min_i = i;
                }
            }
            match min_i {
                0 => {}
                1 => {
                    let t = f[0];
                    f[0] = f[1];
                    f[1] = f[2];
                    f[2] = f[3];
                    f[3] = t;
                }
                2 => {
                    f.swap(0, 2);
                    f.swap(1, 3);
                }
                _ => {
                    let t = f[3];
                    f[3] = f[2];
                    f[2] = f[1];
                    f[1] = f[0];
                    f[0] = t;
                }
            }
            if g[f[1] - 1] > g[f[3] - 1] {
                f.swap(1, 3);
            }
        }
        d
    }

    #[derive(Clone, Copy)]
    struct Cfg {
        sd: bool,
        pe: [i32; 12],
        pf: [i32; 6],
        pb: i32,
        edir: [[usize; 2]; 12],
        fdir: [[usize; 4]; 6],
    }

    fn total(cfg: &Cfg) -> usize {
        // Generous upper bound for the buffers.
        8 + 12 * 6 + 6 * 40 + 400 + (cfg.pb.max(0) as usize).pow(3) + 8
    }

    fn values(cfg: &Cfg, p: [f64; 3]) -> (Vec<f64>, usize) {
        let mut val = vec![0.0; total(cfg)];
        let mut nb = 0;
        brick_nodal(p[0], p[1], p[2], &mut val, &mut nb);
        if cfg.sd {
            sd_brick_edge_p(p[0], p[1], p[2], &cfg.pe, &mut val, &mut nb, &cfg.edir);
            sd_brick_face_p(p[0], p[1], p[2], &cfg.pf, &mut val, &mut nb, &cfg.fdir);
            sd_brick_bubble_p(p[0], p[1], p[2], cfg.pb, &mut val, &mut nb);
        } else {
            brick_edge_p(p[0], p[1], p[2], &cfg.pe, &mut val, &mut nb, &cfg.edir);
            brick_face_p(p[0], p[1], p[2], &cfg.pf, &mut val, &mut nb, &cfg.fdir);
            brick_bubble_p(p[0], p[1], p[2], cfg.pb, &mut val, &mut nb);
        }
        val.truncate(nb);
        (val, nb)
    }

    fn gradients(cfg: &Cfg, p: [f64; 3]) -> (Vec<[f64; 3]>, usize) {
        let mut grad = vec![[0.0; 3]; total(cfg)];
        let mut nb = 0;
        d_brick_nodal(p[0], p[1], p[2], &mut grad, &mut nb);
        if cfg.sd {
            d_sd_brick_edge_p(p[0], p[1], p[2], &cfg.pe, &mut grad, &mut nb, &cfg.edir);
            d_sd_brick_face_p(p[0], p[1], p[2], &cfg.pf, &mut grad, &mut nb, &cfg.fdir);
            d_sd_brick_bubble_p(p[0], p[1], p[2], cfg.pb, &mut grad, &mut nb);
        } else {
            d_brick_edge_p(p[0], p[1], p[2], &cfg.pe, &mut grad, &mut nb, &cfg.edir);
            d_brick_face_p(p[0], p[1], p[2], &cfg.pf, &mut grad, &mut nb, &cfg.fdir);
            d_brick_bubble_p(p[0], p[1], p[2], cfg.pb, &mut grad, &mut nb);
        }
        grad.truncate(nb);
        (grad, nb)
    }

    fn cfg(p: i32, sd: bool, g: &[i32; 8]) -> Cfg {
        Cfg {
            sd,
            pe: [p; 12],
            pf: [p; 6],
            pb: p,
            edir: edge_direction(g),
            fdir: face_direction(g),
        }
    }

    const POINTS: [[f64; 3]; 4] = [
        [0.3, -0.45, 0.62],
        [-0.71, 0.18, -0.33],
        [0.05, 0.9, 0.27],
        [-0.5, -0.6, -0.8],
    ];

    fn fd_check(cfg: &Cfg) {
        let h = 1e-6;
        for p in POINTS {
            let (g, n) = gradients(cfg, p);
            let (_, n0) = values(cfg, p);
            assert_eq!(n, n0);
            for c in 0..3 {
                let mut pp = p;
                let mut pm = p;
                pp[c] += h;
                pm[c] -= h;
                let (vp, _) = values(cfg, pp);
                let (vm, _) = values(cfg, pm);
                for i in 0..n {
                    let fd = (vp[i] - vm[i]) / (2.0 * h);
                    let err = (fd - g[i][c]).abs();
                    assert!(
                        err <= 1e-6 * (1.0 + g[i][c].abs()),
                        "sd={} p={} fn {} comp {}: fd {} vs analytic {}",
                        cfg.sd,
                        cfg.pb,
                        i,
                        c,
                        fd,
                        g[i][c]
                    );
                }
            }
        }
    }

    #[test]
    fn nodal_partition_of_unity() {
        for p in POINTS {
            let mut val = [0.0; 8];
            let mut g = [[0.0; 3]; 8];
            let mut nb = 0;
            brick_nodal(p[0], p[1], p[2], &mut val, &mut nb);
            assert_eq!(nb, 8);
            let mut nb = 0;
            d_brick_nodal(p[0], p[1], p[2], &mut g, &mut nb);
            assert!((val.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            for c in 0..3 {
                assert!(g.iter().map(|x| x[c]).sum::<f64>().abs() < 1e-14);
            }
            // nodal reproduces the coordinates
            let xs = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
            let ys = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
            let zs = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
            let x: f64 = (0..8).map(|i| val[i] * xs[i]).sum();
            let y: f64 = (0..8).map(|i| val[i] * ys[i]).sum();
            let z: f64 = (0..8).map(|i| val[i] * zs[i]).sum();
            assert!((x - p[0]).abs() < 1e-14 && (y - p[1]).abs() < 1e-14 && (z - p[2]).abs() < 1e-14);
        }
    }

    #[test]
    fn lambdas_and_gradients_are_consistent() {
        let p = [0.21, -0.37, 0.58];
        let h = 1e-6;
        for node in 1..=8 {
            let g = d_brick_l(node);
            for c in 0..3 {
                let mut a = p;
                let mut b = p;
                a[c] += h;
                b[c] -= h;
                let fd = (brick_l(node, a[0], a[1], a[2]) - brick_l(node, b[0], b[1], b[2])) / (2.0 * h);
                assert!((fd - g[c]).abs() < 1e-8);
            }
        }
        for e in 1..=12 {
            let (da, db) = d_brick_edge_l(e);
            for c in 0..3 {
                let mut a = p;
                let mut b = p;
                a[c] += h;
                b[c] -= h;
                let (la1, lb1) = brick_edge_l(e, a[0], a[1], a[2]);
                let (la0, lb0) = brick_edge_l(e, b[0], b[1], b[2]);
                assert!(((la1 - la0) / (2.0 * h) - da[c]).abs() < 1e-8);
                assert!(((lb1 - lb0) / (2.0 * h) - db[c]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn finite_difference_natural_directions() {
        let nat = [1, 2, 3, 4, 5, 6, 7, 8];
        for p in 2..=5 {
            fd_check(&cfg(p, false, &nat));
            fd_check(&cfg(p, true, &nat));
        }
    }

    #[test]
    fn finite_difference_permuted_global_numbering() {
        // Swapped edge directions and rotated / reflected face orderings.
        let numberings: [[i32; 8]; 4] = [
            [8, 7, 6, 5, 4, 3, 2, 1],
            [3, 9, 1, 12, 7, 2, 11, 5],
            [20, 4, 15, 8, 2, 30, 6, 11],
            [5, 1, 8, 3, 9, 4, 7, 2],
        ];
        for g in &numberings {
            for p in 2..=4 {
                fd_check(&cfg(p, false, g));
                fd_check(&cfg(p, true, g));
            }
        }
    }

    #[test]
    fn finite_difference_mixed_degrees() {
        let g = [3, 9, 1, 12, 7, 2, 11, 5];
        for sd in [false, true] {
            let mut c = cfg(4, sd, &g);
            c.pe = [2, 3, 4, 5, 1, 2, 3, 4, 5, 2, 3, 4];
            c.pf = [2, 3, 4, 5, 3, 4];
            c.pb = 5;
            fd_check(&c);
        }
    }

    #[test]
    fn counts_match_element_info_formulas() {
        let nat = [1, 2, 3, 4, 5, 6, 7, 8];
        for p in 1..=6i32 {
            for sd in [false, true] {
                let c = cfg(p, sd, &nat);
                let (_, n) = values(&c, POINTS[0]);
                let (_, ng) = gradients(&c, POINTS[0]);
                assert_eq!(n, ng);
                let edges = 12 * (p - 1).max(0) as usize;
                let faces: usize = if sd {
                    // nbq in ElementInfoVec: sum_j=2..p max(p-j-1,0)
                    6 * (2..=p).map(|j| (p - j - 1).max(0) as usize).sum::<usize>()
                } else {
                    6 * ((p - 1).max(0) as usize).pow(2)
                };
                let bubble: usize = if sd {
                    let mut s = 0usize;
                    for i in 2..=p - 4 {
                        for j in 2..=p - i - 2 {
                            s += (p - i - j - 1).max(0) as usize;
                        }
                    }
                    s
                } else {
                    ((p - 1).max(0) as usize).pow(3)
                };
                assert_eq!(n, 8 + edges + faces + bubble, "p={p} sd={sd}");
            }
        }
        // Bubble counts quoted in ElemInfo.F90: i^3 for the full scheme ...
        for (p, expect) in [(1, 0), (2, 1), (3, 8), (4, 27), (5, 64)] {
            let mut val = vec![0.0; 100];
            let mut nb = 0;
            brick_bubble_p(0.1, 0.2, 0.3, p, &mut val, &mut nb);
            assert_eq!(nb, expect);
        }
        // ... and the serendipity bubbles start at p=6 with 1 function.
        let mut val = vec![0.0; 100];
        let mut nb = 0;
        sd_brick_bubble_p(0.1, 0.2, 0.3, 6, &mut val, &mut nb);
        assert_eq!(nb, 1);
    }

    #[test]
    fn edge_functions_vanish_on_other_edges_and_nodes() {
        // The hierarchic edge function of a full-polynomial scheme vanishes at all nodes.
        let nat = [1, 2, 3, 4, 5, 6, 7, 8];
        let c = cfg(3, false, &nat);
        let corners = [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        for p in corners {
            let (v, _) = values(&c, p);
            for x in &v[8..] {
                assert!(x.abs() < 1e-13, "{x}");
            }
        }
    }
}
