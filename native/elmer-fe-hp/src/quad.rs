//! Quadrilateral (404) p basis: the `H1Basis_*` routines `QuadNodal`,
//! `QuadL`, `QuadEdgeP`, `QuadBubbleP`, their serendipity `SD_*` variants and
//! derivatives. Conventions: see [`crate::types`].

use crate::legendre::{d_legendre_p, d_phi, d_var_phi, legendre_p, phi, var_phi};

/// `H1Basis_QuadNodal`. Fortran quirk: the four values are stored at slots
/// `0..4`, not `nb..nb+4`, while `nb` is still advanced by 4.
pub fn quad_nodal(u: f64, v: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 4.0;
    val[0] = c * (1.0 - u) * (1.0 - v);
    val[1] = c * (1.0 + u) * (1.0 - v);
    val[2] = c * (1.0 + u) * (1.0 + v);
    val[3] = c * (1.0 - u) * (1.0 + v);
    *nb += 4;
}

/// `H1Basis_dQuadNodal`.
pub fn d_quad_nodal(u: f64, v: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let c = 1.0 / 4.0;
    let n = *nb;
    grad[n][0] = -c * (1.0 - v);
    grad[n + 1][0] = c * (1.0 - v);
    grad[n + 2][0] = c * (1.0 + v);
    grad[n + 3][0] = -c * (1.0 + v);
    grad[n][1] = -c * (1.0 - u);
    grad[n + 1][1] = -c * (1.0 + u);
    grad[n + 2][1] = c * (1.0 + u);
    grad[n + 3][1] = c * (1.0 - u);
    *nb += 4;
}

/// `H1Basis_QuadL(node, u, v)`: the affine "L" coordinate of local node
/// `node` (1-based), `(2 -/+ u -/+ v)/2`.
pub fn quad_l(node: usize, u: f64, v: f64) -> f64 {
    let c = 1.0 / 2.0;
    match node {
        1 => c * (2.0 - u - v),
        2 => c * (2.0 + u - v),
        3 => c * (2.0 + u + v),
        4 => c * (2.0 - u + v),
        _ => panic!("quad_l: node {node} outside 1..=4"),
    }
}

/// `H1Basis_dQuadL(node)`.
pub fn d_quad_l(node: usize) -> [f64; 2] {
    let c = 1.0 / 2.0;
    match node {
        1 => [-c, -c],
        2 => [c, -c],
        3 => [c, c],
        4 => [-c, c],
        _ => panic!("d_quad_l: node {node} outside 1..=4"),
    }
}

/// `H1Basis_SD_QuadEdgeP`: serendipity edge functions.
pub fn sd_quad_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0;
    for i in 0..4 {
        let la = quad_l(edgedir[i][0], u, v);
        let lb = quad_l(edgedir[i][1], u, v);
        for j in 1..pmax[i] {
            val[*nb + j as usize - 1] = c * (la + lb - 1.0) * phi(j as usize + 1, lb - la);
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_SD_dQuadEdgeP`.
pub fn sd_d_quad_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0;
    for i in 0..4 {
        let d_la = d_quad_l(edgedir[i][0]);
        let d_lb = d_quad_l(edgedir[i][1]);
        let la = quad_l(edgedir[i][0], u, v);
        let lb = quad_l(edgedir[i][1], u, v);
        for j in 1..pmax[i] {
            let p = phi(j as usize + 1, lb - la);
            let dp = d_phi(j as usize + 1, lb - la);
            let idx = *nb + j as usize - 1;
            grad[idx][0] = c * ((d_la[0] + d_lb[0]) * p + (la + lb - 1.0) * dp * (d_lb[0] - d_la[0]));
            grad[idx][1] = c * ((d_la[1] + d_lb[1]) * p + (la + lb - 1.0) * dp * (d_lb[1] - d_la[1]));
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_SD_QuadBubbleP`: serendipity bubbles (undirected or, with
/// `local_numbers`, directed by the face's local node order).
pub fn sd_quad_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    val: &mut [f64],
    nb: &mut usize,
    local_numbers: Option<[usize; 4]>,
) {
    match local_numbers {
        None => {
            for i in 2..=pmax - 2 {
                for j in 1..pmax - i {
                    val[*nb + j as usize - 1] = phi(i as usize, u) * phi(j as usize + 1, v);
                }
                *nb += (pmax - i - 1).max(0) as usize;
            }
        }
        Some(ln) => {
            let la = quad_l(ln[0], u, v);
            let lb = quad_l(ln[1], u, v);
            let lc = quad_l(ln[3], u, v);
            for i in 2..=pmax - 2 {
                for j in 1..pmax - i {
                    val[*nb + j as usize - 1] =
                        phi(i as usize, lb - la) * phi(j as usize + 1, lc - la);
                }
                *nb += (pmax - i - 1).max(0) as usize;
            }
        }
    }
}

/// `H1Basis_SD_dQuadBubbleP`.
pub fn sd_d_quad_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    local_numbers: Option<[usize; 4]>,
) {
    match local_numbers {
        None => {
            for i in 2..=pmax - 2 {
                for j in 1..pmax - i {
                    let idx = *nb + j as usize - 1;
                    grad[idx][0] = d_phi(i as usize, u) * phi(j as usize + 1, v);
                    grad[idx][1] = phi(i as usize, u) * d_phi(j as usize + 1, v);
                }
                *nb += (pmax - i - 1).max(0) as usize;
            }
        }
        Some(ln) => {
            let d_la = d_quad_l(ln[0]);
            let d_lb = d_quad_l(ln[1]);
            let d_lc = d_quad_l(ln[3]);
            let d_lb_d_la = [d_lb[0] - d_la[0], d_lb[1] - d_la[1]];
            let d_lc_d_la = [d_lc[0] - d_la[0], d_lc[1] - d_la[1]];
            let la = quad_l(ln[0], u, v);
            let lb = quad_l(ln[1], u, v);
            let lc = quad_l(ln[3], u, v);
            for i in 2..=pmax - 2 {
                for j in 1..pmax - i {
                    let idx = *nb + j as usize - 1;
                    for c in 0..2 {
                        grad[idx][c] = d_phi(i as usize, lb - la)
                            * d_lb_d_la[c]
                            * phi(j as usize + 1, lc - la)
                            + phi(i as usize, lb - la)
                                * d_phi(j as usize + 1, lc - la)
                                * d_lc_d_la[c];
                    }
                }
                *nb += (pmax - i - 1).max(0) as usize;
            }
        }
    }
}

/// `H1Basis_QuadEdgeP`: edge functions `N_a*N_b*varphi_{j+1}(Lb-La)`.
/// Fortran quirk: `N_a`, `N_b` are read from `val[node-1]`, so the nodal
/// values must have been stored there first.
pub fn quad_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0 * 2.0;
    for i in 0..4 {
        let node1 = edgedir[i][0];
        let node2 = edgedir[i][1];
        let la = quad_l(node1, u, v);
        let lb = quad_l(node2, u, v);
        for j in 1..pmax[i] {
            let na = val[node1 - 1];
            let nbv = val[node2 - 1];
            val[*nb + j as usize - 1] = c * na * nbv * var_phi(j as usize + 1, lb - la);
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_dQuadEdgeP`. The nodal values come from a fresh `QuadNodal`
/// evaluation, the nodal gradients (`dNa`, `dNb`) from `grad[node-1]`, which
/// the caller must have filled with [`d_quad_nodal`].
pub fn d_quad_edge_p(
    u: f64,
    v: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let c = 1.0 / 2.0 * 2.0;
    let mut n = [0.0; 4];
    let mut nnb = 0;
    quad_nodal(u, v, &mut n, &mut nnb);
    for i in 0..4 {
        let node1 = edgedir[i][0];
        let node2 = edgedir[i][1];
        let d_la = d_quad_l(node1);
        let d_lb = d_quad_l(node2);
        let la = quad_l(node1, u, v);
        let lb = quad_l(node2, u, v);
        for j in 1..pmax[i] {
            let na = n[node1 - 1];
            let nbv = n[node2 - 1];
            let d_na = [grad[node1 - 1][0], grad[node1 - 1][1]];
            let d_nb = [grad[node2 - 1][0], grad[node2 - 1][1]];
            let p = var_phi(j as usize + 1, lb - la);
            let dp = d_var_phi(j as usize + 1, lb - la);
            let idx = *nb + j as usize - 1;
            for k in 0..2 {
                grad[idx][k] =
                    c * (d_na[k] * nbv * p + na * d_nb[k] * p + na * nbv * dp * (d_lb[k] - d_la[k]));
            }
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_QuadBubbleP`. With `local_numbers` the directed face bubble
/// reads `Pa`, `Pb` from the nodal values stored at `val[ln[0]-1]` and
/// `val[ln[2]-1]` (Fortran quirk, see [`quad_edge_p`]).
pub fn quad_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    val: &mut [f64],
    nb: &mut usize,
    local_numbers: Option<[usize; 4]>,
) {
    match local_numbers {
        None => {
            for i in 2..=pmax {
                for j in 1..pmax {
                    val[*nb + j as usize - 1] = phi(i as usize, u) * phi(j as usize + 1, v);
                }
                *nb += (pmax - 1).max(0) as usize;
            }
        }
        Some(ln) => {
            let la = quad_l(ln[0], u, v);
            let lb = quad_l(ln[1], u, v);
            let lc = quad_l(ln[3], u, v);
            for i in 0..=pmax - 2 {
                for j in 1..pmax {
                    let pa = val[ln[0] - 1];
                    let pb = val[ln[2] - 1];
                    val[*nb + j as usize - 1] = pa
                        * pb
                        * legendre_p(i as usize, lb - la)
                        * legendre_p((j - 1) as usize, lc - la);
                }
                *nb += (pmax - 1).max(0) as usize;
            }
        }
    }
}

/// `H1Basis_dQuadBubbleP`; with `local_numbers`, nodal values come from a
/// fresh `QuadNodal` evaluation and nodal gradients from `grad[ln[k]-1]`.
pub fn d_quad_bubble_p(
    u: f64,
    v: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    local_numbers: Option<[usize; 4]>,
) {
    match local_numbers {
        None => {
            for i in 2..=pmax {
                for j in 1..pmax {
                    let idx = *nb + j as usize - 1;
                    grad[idx][0] = d_phi(i as usize, u) * phi(j as usize + 1, v);
                    grad[idx][1] = phi(i as usize, u) * d_phi(j as usize + 1, v);
                }
                *nb += (pmax - 1).max(0) as usize;
            }
        }
        Some(ln) => {
            let mut n = [0.0; 4];
            let mut nnb = 0;
            quad_nodal(u, v, &mut n, &mut nnb);
            let d_la = d_quad_l(ln[0]);
            let d_lb = d_quad_l(ln[1]);
            let d_lc = d_quad_l(ln[3]);
            let la = quad_l(ln[0], u, v);
            let lb = quad_l(ln[1], u, v);
            let lc = quad_l(ln[3], u, v);
            for i in 0..=pmax - 2 {
                for j in 1..pmax {
                    let pa = n[ln[0] - 1];
                    let pb = n[ln[2] - 1];
                    let d_pa = [grad[ln[0] - 1][0], grad[ln[0] - 1][1]];
                    let d_pb = [grad[ln[2] - 1][0], grad[ln[2] - 1][1]];
                    let legi = legendre_p(i as usize, lb - la);
                    let legj = legendre_p((j - 1) as usize, lc - la);
                    let dl_i = d_legendre_p(i as usize, lb - la);
                    let dl_j = d_legendre_p((j - 1) as usize, lc - la);
                    let idx = *nb + j as usize - 1;
                    for k in 0..2 {
                        let d_legi = dl_i * (d_lb[k] - d_la[k]);
                        let d_legj = dl_j * (d_lc[k] - d_la[k]);
                        grad[idx][k] = d_pa[k] * pb * legi * legj
                            + pa * d_pb[k] * legi * legj
                            + pa * pb * d_legi * legj
                            + pa * pb * legi * d_legj;
                    }
                }
                *nb += (pmax - 1).max(0) as usize;
            }
        }
    }
}
