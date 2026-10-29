//! Pyramid p-element basis (`H1Basis_Pyramid*` in `H1Basis.F90`).
//!
//! Scalar transliterations; see [`crate::types`] for the calling conventions.
//! Several routines here read the nodal values (`val[0..5]`) or nodal
//! gradients (`grad[0..5]`) that the caller stored first, as the Fortran does.

#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::many_single_char_names,
    clippy::excessive_precision
)]

use crate::legendre::{d_legendre_p, d_var_phi, legendre_p, var_phi};

/// `H1Basis_PyramidNodalP`: the five nodal functions, appended at `*nb`.
pub fn pyramid_nodal_p(u: f64, v: f64, w: f64, val: &mut [f64], nb: &mut usize) {
    let sq2 = 2.0f64.sqrt();
    let s = w / sq2;
    val[*nb] = (1.0 - u - v - s + u * v / (1.0 - s)) / 4.0;
    val[*nb + 1] = (1.0 + u - v - s - u * v / (1.0 - s)) / 4.0;
    val[*nb + 2] = (1.0 + u + v - s + u * v / (1.0 - s)) / 4.0;
    val[*nb + 3] = (1.0 - u + v - s - u * v / (1.0 - s)) / 4.0;
    val[*nb + 4] = s;
    *nb += 5;
}

/// `H1Basis_dPyramidNodalP`: gradients of the five nodal functions.
pub fn d_pyramid_nodal_p(u: f64, v: f64, w: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let sq2 = 2.0f64.sqrt();
    let sq42 = 4.0 * 2.0f64.sqrt();
    let s = v / (1.0 - w / sq2);
    grad[*nb][0] = (-1.0 + s) / 4.0;
    grad[*nb + 1][0] = (1.0 - s) / 4.0;
    grad[*nb + 2][0] = (1.0 + s) / 4.0;
    grad[*nb + 3][0] = (-1.0 - s) / 4.0;
    grad[*nb + 4][0] = 0.0;
    let s = u / (1.0 - w / sq2);
    grad[*nb][1] = (-1.0 + s) / 4.0;
    grad[*nb + 1][1] = (-1.0 - s) / 4.0;
    grad[*nb + 2][1] = (1.0 + s) / 4.0;
    grad[*nb + 3][1] = (1.0 - s) / 4.0;
    grad[*nb + 4][1] = 0.0;
    let t = 1.0 - w / sq2;
    let s = u * v / (t * t);
    grad[*nb][2] = (-1.0 + s) / sq42;
    grad[*nb + 1][2] = (-1.0 - s) / sq42;
    grad[*nb + 2][2] = (-1.0 + s) / sq42;
    grad[*nb + 3][2] = (-1.0 - s) / sq42;
    grad[*nb + 4][2] = 1.0 / sq2;
    *nb += 5;
}

/// `H1Basis_PyramidL`: affine coordinate of square-face node `which` (1..=4).
pub fn pyramid_l(which: usize, u: f64, v: f64) -> f64 {
    match which {
        1 => ((1.0 - u) + (1.0 - v)) / 2.0,
        2 => ((1.0 + u) + (1.0 - v)) / 2.0,
        3 => ((1.0 + u) + (1.0 + v)) / 2.0,
        4 => ((1.0 - u) + (1.0 + v)) / 2.0,
        _ => panic!("pyramid_l: unknown affine coordinate {which} for square face"),
    }
}

/// `H1Basis_dPyramidL`: gradient of [`pyramid_l`].
pub fn d_pyramid_l(which: usize) -> [f64; 3] {
    match which {
        1 => [-1.0 / 2.0, -1.0 / 2.0, 0.0],
        2 => [1.0 / 2.0, -1.0 / 2.0, 0.0],
        3 => [1.0 / 2.0, 1.0 / 2.0, 0.0],
        4 => [-1.0 / 2.0, 1.0 / 2.0, 0.0],
        _ => panic!("d_pyramid_l: unknown affine coordinate {which} for square face"),
    }
}

/// `H1Basis_PyramidTL`: affine coordinate `which` (1..=5) for triangular faces.
pub fn pyramid_tl(which: usize, u: f64, v: f64, w: f64) -> f64 {
    let s = w / 2.0f64.sqrt();
    match which {
        1 => (2.0 - u - v - s) / 2.0,
        2 => (2.0 + u - v - s) / 2.0,
        3 => (2.0 + u + v - s) / 2.0,
        4 => (2.0 - u + v - s) / 2.0,
        5 => s,
        _ => panic!("pyramid_tl: unknown function {which}"),
    }
}

/// `H1Basis_dPyramidTL`: gradient of [`pyramid_tl`].
pub fn d_pyramid_tl(which: usize) -> [f64; 3] {
    let mut g = match which {
        1 => [-1.0, -1.0, -1.0],
        2 => [1.0, -1.0, -1.0],
        3 => [1.0, 1.0, -1.0],
        4 => [-1.0, 1.0, -1.0],
        5 => [0.0, 0.0, 2.0],
        _ => panic!("d_pyramid_tl: unknown function {which}"),
    };
    for x in g.iter_mut() {
        *x /= 2.0;
    }
    g[2] /= 2.0f64.sqrt();
    g
}

/// `H1Basis_PyramidEdgeP`: the eight edge functions. Reads the nodal values
/// from `val[0..5]`, which the caller must have stored with
/// [`pyramid_nodal_p`].
pub fn pyramid_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 0..8 {
        let node1 = edgedir[i][0];
        let node2 = edgedir[i][1];
        let (la, lb) = if i < 4 {
            // square face edges
            (pyramid_l(node1, u, v), pyramid_l(node2, u, v))
        } else {
            // triangular face edges
            (pyramid_tl(node1, u, v, w), pyramid_tl(node2, u, v, w))
        };
        for j in 1..=pmax[i] - 1 {
            let na = val[node1 - 1];
            let nb_ = val[node2 - 1];
            val[*nb + j as usize - 1] = na * nb_ * var_phi(j as usize + 1, lb - la);
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_dPyramidEdgeP`: gradients matching [`pyramid_edge_p`]. Reads the
/// nodal gradients from `grad[0..5]` (stored by [`d_pyramid_nodal_p`]); the
/// nodal values are recomputed internally.
pub fn d_pyramid_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    let mut n = [0.0; 5];
    let mut nnb = 0;
    pyramid_nodal_p(u, v, w, &mut n, &mut nnb);
    for i in 0..8 {
        let node1 = edgedir[i][0];
        let node2 = edgedir[i][1];
        let (dla, dlb, la, lb) = if i < 4 {
            (
                d_pyramid_l(node1),
                d_pyramid_l(node2),
                pyramid_l(node1, u, v),
                pyramid_l(node2, u, v),
            )
        } else {
            (
                d_pyramid_tl(node1),
                d_pyramid_tl(node2),
                pyramid_tl(node1, u, v, w),
                pyramid_tl(node2, u, v, w),
            )
        };
        for j in 1..=pmax[i] - 1 {
            let na = n[node1 - 1];
            let nb_ = n[node2 - 1];
            let dna = grad[node1 - 1];
            let dnb = grad[node2 - 1];
            let phi = var_phi(j as usize + 1, lb - la);
            let dv = d_var_phi(j as usize + 1, lb - la);
            let g = &mut grad[*nb + j as usize - 1];
            for c in 0..3 {
                let dphi = dv * (dlb[c] - dla[c]);
                g[c] = dna[c] * nb_ * phi + na * dnb[c] * phi + na * nb_ * dphi;
            }
        }
        *nb += (pmax[i] - 1).max(0) as usize;
    }
}

/// `H1Basis_PyramidFaceP`: face functions (face 1 square, faces 2..=5
/// triangular). Reads nodal values from `val[0..5]`.
pub fn pyramid_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    // square face
    {
        let i = 0;
        let node1 = facedir[i][0];
        let node2 = facedir[i][1];
        let node3 = facedir[i][2];
        let node4 = facedir[i][3];
        let la = pyramid_l(node1, u, v);
        let lb = pyramid_l(node2, u, v);
        let lc = pyramid_l(node4, u, v);
        for j in 0..=pmax[i] - 2 {
            for k in 1..=pmax[i] - 1 {
                let pa = val[node1 - 1];
                let pb = val[node3 - 1];
                val[*nb + k as usize - 1] =
                    pa * pb * legendre_p(j as usize, lb - la) * legendre_p(k as usize - 1, lc - la);
            }
            *nb = (*nb as i64 + pmax[i] as i64 - 1) as usize;
        }
    }
    // triangular faces
    for i in 1..5 {
        let node1 = facedir[i][0];
        let node2 = facedir[i][1];
        let node3 = facedir[i][2];
        let la = pyramid_tl(node1, u, v, w);
        let lb = pyramid_tl(node2, u, v, w);
        let lc = pyramid_tl(node3, u, v, w);
        for j in 0..=pmax[i] - 3 {
            for k in 1..=pmax[i] - 2 - j {
                let pa = val[node1 - 1];
                let pb = val[node2 - 1];
                let pc = val[node3 - 1];
                val[*nb + k as usize - 1] = pa
                    * pb
                    * pc
                    * legendre_p(j as usize, lb - la)
                    * legendre_p(k as usize - 1, 2.0 * lc - 1.0);
            }
            *nb += (pmax[i] - j - 2).max(0) as usize;
        }
    }
}

/// `H1Basis_dPyramidFaceP`: gradients matching [`pyramid_face_p`]. Reads
/// nodal gradients from `grad[0..5]`; nodal values are recomputed.
pub fn d_pyramid_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    let mut n = [0.0; 5];
    let mut nnb = 0;
    pyramid_nodal_p(u, v, w, &mut n, &mut nnb);
    // square face
    {
        let i = 0;
        let node1 = facedir[i][0];
        let node2 = facedir[i][1];
        let node3 = facedir[i][2];
        let node4 = facedir[i][3];
        let dla = d_pyramid_l(node1);
        let dlb = d_pyramid_l(node2);
        let dlc = d_pyramid_l(node4);
        let la = pyramid_l(node1, u, v);
        let lb = pyramid_l(node2, u, v);
        let lc = pyramid_l(node4, u, v);
        for j in 0..=pmax[i] - 2 {
            for k in 1..=pmax[i] - 1 {
                let pa = n[node1 - 1];
                let pb = n[node3 - 1];
                let dpa = grad[node1 - 1];
                let dpb = grad[node3 - 1];
                let legi = legendre_p(j as usize, lb - la);
                let legj = legendre_p(k as usize - 1, lc - la);
                let g = &mut grad[*nb + k as usize - 1];
                for c in 0..3 {
                    let dlegi = d_legendre_p(j as usize, lb - la) * (dlb[c] - dla[c]);
                    let dlegj = d_legendre_p(k as usize - 1, lc - la) * (dlc[c] - dla[c]);
                    g[c] = dpa[c] * pb * legi * legj
                        + pa * dpb[c] * legi * legj
                        + pa * pb * dlegi * legj
                        + pa * pb * legi * dlegj;
                }
            }
            *nb = (*nb as i64 + pmax[i] as i64 - 1) as usize;
        }
    }
    // triangular faces
    for i in 1..5 {
        let node1 = facedir[i][0];
        let node2 = facedir[i][1];
        let node3 = facedir[i][2];
        let dla = d_pyramid_tl(node1);
        let dlb = d_pyramid_tl(node2);
        let dlc = d_pyramid_tl(node3);
        let la = pyramid_tl(node1, u, v, w);
        let lb = pyramid_tl(node2, u, v, w);
        let lc = pyramid_tl(node3, u, v, w);
        for j in 0..=pmax[i] - 3 {
            for k in 1..=pmax[i] - 2 - j {
                let pa = n[node1 - 1];
                let pb = n[node2 - 1];
                let pc = n[node3 - 1];
                let dpa = grad[node1 - 1];
                let dpb = grad[node2 - 1];
                let dpc = grad[node3 - 1];
                let legi = legendre_p(j as usize, lb - la);
                let legj = legendre_p(k as usize - 1, 2.0 * lc - 1.0);
                let g = &mut grad[*nb + k as usize - 1];
                for c in 0..3 {
                    let dlegi = d_legendre_p(j as usize, lb - la) * (dlb[c] - dla[c]);
                    let dlegj = d_legendre_p(k as usize - 1, 2.0 * lc - 1.0) * 2.0 * dlc[c];
                    g[c] = dpa[c] * pb * pc * legi * legj
                        + pa * dpb[c] * pc * legi * legj
                        + pa * pb * dpc[c] * legi * legj
                        + pa * pb * pc * dlegi * legj
                        + pa * pb * pc * legi * dlegj;
                }
            }
            *nb += (pmax[i] - j - 2).max(0) as usize;
        }
    }
}

/// `H1Basis_PyramidBubbleP`: interior bubbles of degree `pmax`. Reads the
/// nodal values 1, 3 and 5 from `val[0]`, `val[2]`, `val[4]`.
pub fn pyramid_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    for i in 0..=pmax - 3 {
        for j in 0..=pmax - i - 3 {
            for k in 1..=pmax - i - j - 2 {
                *nb += 1;
                let s = w / 2.0f64.sqrt();
                let pa = val[0];
                let pb = val[2];
                let pc = val[4];
                let legi = legendre_p(i as usize, u);
                let legj = legendre_p(j as usize, v);
                let legk = legendre_p(k as usize - 1, 2.0 * s - 1.0);
                val[*nb - 1] = pa * pb * pc * legi * legj * legk;
            }
        }
    }
}

/// `H1Basis_dPyramidBubbleP`: gradients matching [`pyramid_bubble_p`]. Reads
/// nodal gradients 1, 3, 5 from `grad[0]`, `grad[2]`, `grad[4]`.
pub fn d_pyramid_bubble_p(u: f64, v: f64, w: f64, pmax: i32, grad: &mut [[f64; 3]], nb: &mut usize) {
    let a = 1.0 / 2.0f64.sqrt();
    let mut n = [0.0; 5];
    let mut nnb = 0;
    pyramid_nodal_p(u, v, w, &mut n, &mut nnb);
    for i in 0..=pmax - 3 {
        for j in 0..=pmax - i - 3 {
            for k in 1..=pmax - i - j - 2 {
                *nb += 1;
                let s = a * w;
                let pa = n[0];
                let pb = n[2];
                let pc = n[4];
                let dpa = grad[0];
                let dpb = grad[2];
                let dpc = grad[4];
                let legi = legendre_p(i as usize, u);
                let legj = legendre_p(j as usize, v);
                let legk = legendre_p(k as usize - 1, 2.0 * s - 1.0);
                let dlegi = [d_legendre_p(i as usize, u), 0.0, 0.0];
                let dlegj = [0.0, d_legendre_p(j as usize, v), 0.0];
                let dlegk = [0.0, 0.0, d_legendre_p(k as usize - 1, 2.0 * s - 1.0) * 2.0 * a];
                let g = &mut grad[*nb - 1];
                for c in 0..3 {
                    g[c] = dpa[c] * pb * pc * legi * legj * legk
                        + pa * dpb[c] * pc * legi * legj * legk
                        + pa * pb * dpc[c] * legi * legj * legk
                        + pa * pb * pc * dlegi[c] * legj * legk
                        + pa * pb * pc * legi * dlegj[c] * legk
                        + pa * pb * pc * legi * legj * dlegk[c];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `H1Basis_GetEdgeMap(605)`.
    const EDGES: [[usize; 2]; 8] = [
        [1, 2],
        [2, 3],
        [4, 3],
        [1, 4],
        [1, 5],
        [2, 5],
        [3, 5],
        [4, 5],
    ];
    /// `H1Basis_GetFaceMap(605)`.
    const FACES: [[usize; 4]; 5] = [
        [1, 2, 3, 4],
        [1, 2, 5, 0],
        [2, 3, 5, 0],
        [3, 4, 5, 0],
        [4, 1, 5, 0],
    ];

    /// Interior points of the reference pyramid |u|,|v| < 1 - w/sqrt2.
    fn points() -> Vec<[f64; 3]> {
        let sq2 = 2f64.sqrt();
        let mut pts = Vec::new();
        for (a, b, t) in [
            (0.1, 0.2, 0.2),
            (-0.3, 0.25, 0.4),
            (0.35, -0.15, 0.1),
            (0.05, 0.05, 0.6),
        ] {
            let w = t * sq2;
            let s = 1.0 - t;
            pts.push([a * s, b * s, w]);
        }
        pts
    }

    /// Fill `o[0..5]` with nodal values then call `f`.
    type ValFn<'a> = &'a dyn Fn(f64, f64, f64, &mut Vec<f64>);
    type GradFn<'a> = &'a dyn Fn(f64, f64, f64, &mut Vec<[f64; 3]>);

    fn check_fd(
        f: ValFn,
        df: GradFn,
        first: usize,
        n: usize,
    ) {
        let h = 1e-6;
        for p in points() {
            let mut g = vec![[0.0; 3]; 200];
            df(p[0], p[1], p[2], &mut g);
            for c in 0..3 {
                let mut pp = p;
                let mut pm = p;
                pp[c] += h;
                pm[c] -= h;
                let mut vp = vec![0.0; 200];
                let mut vm = vec![0.0; 200];
                f(pp[0], pp[1], pp[2], &mut vp);
                f(pm[0], pm[1], pm[2], &mut vm);
                for m in first..first + n {
                    let fd = (vp[m] - vm[m]) / (2.0 * h);
                    let tol = 1e-6 * (1.0 + g[m][c].abs());
                    assert!(
                        (fd - g[m][c]).abs() < tol,
                        "fn {m} comp {c}: fd {fd} vs {}",
                        g[m][c]
                    );
                }
            }
        }
    }

    #[test]
    fn nodal_partition_of_unity_and_gradients() {
        for p in points() {
            let mut v = [0.0; 5];
            let mut nb = 0;
            pyramid_nodal_p(p[0], p[1], p[2], &mut v, &mut nb);
            assert_eq!(nb, 5);
            assert!((v.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let mut g = [[0.0; 3]; 5];
            let mut nb = 0;
            d_pyramid_nodal_p(p[0], p[1], p[2], &mut g, &mut nb);
            assert_eq!(nb, 5);
            for c in 0..3 {
                assert!(g.iter().map(|x| x[c]).sum::<f64>().abs() < 1e-13);
            }
        }
        check_fd(
            &|u, v, w, o| {
                let mut nb = 0;
                pyramid_nodal_p(u, v, w, o, &mut nb);
            },
            &|u, v, w, g| {
                let mut nb = 0;
                d_pyramid_nodal_p(u, v, w, g, &mut nb);
            },
            0,
            5,
        );
    }

    #[test]
    fn affine_coordinate_gradients() {
        let h = 1e-6;
        let (u, v, w) = (0.1, 0.2, 0.3);
        for which in 1..=4 {
            let g = d_pyramid_l(which);
            let fd = [
                (pyramid_l(which, u + h, v) - pyramid_l(which, u - h, v)) / (2.0 * h),
                (pyramid_l(which, u, v + h) - pyramid_l(which, u, v - h)) / (2.0 * h),
                0.0,
            ];
            for c in 0..3 {
                assert!((fd[c] - g[c]).abs() < 1e-8);
            }
        }
        for which in 1..=5 {
            let g = d_pyramid_tl(which);
            let fd = [
                (pyramid_tl(which, u + h, v, w) - pyramid_tl(which, u - h, v, w)) / (2.0 * h),
                (pyramid_tl(which, u, v + h, w) - pyramid_tl(which, u, v - h, w)) / (2.0 * h),
                (pyramid_tl(which, u, v, w + h) - pyramid_tl(which, u, v, w - h)) / (2.0 * h),
            ];
            for c in 0..3 {
                assert!((fd[c] - g[c]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn edge_face_bubble_gradients_and_counts() {
        for p in 2..=5i32 {
            let pe = [p; 8];
            let pf = [p; 5];
            let ne = 8 * (p as usize - 1);
            // ElemInfo nbq formulae for the pyramid
            let mut nf = 0usize;
            for _ in 0..1 {
                for _ in 0..=p - 2 {
                    nf += (p - 1).max(0) as usize;
                }
            }
            for _ in 2..=5 {
                for j in 0..=p - 3 {
                    nf += (p - j - 2).max(0) as usize;
                }
            }
            let nbub: usize = (0..=p - 3)
                .map(|i| (0..=p - i - 3).map(|j| (p - i - j - 2).max(0) as usize).sum::<usize>())
                .sum();

            let nodal = |u: f64, v: f64, w: f64, o: &mut Vec<f64>| {
                let mut nb = 0;
                pyramid_nodal_p(u, v, w, o, &mut nb);
                nb
            };
            let dnodal = |u: f64, v: f64, w: f64, g: &mut Vec<[f64; 3]>| {
                let mut nb = 0;
                d_pyramid_nodal_p(u, v, w, g, &mut nb);
                nb
            };

            // counts
            let mut v = vec![0.0; 200];
            let mut nb = nodal(0.1, 0.1, 0.3, &mut v);
            pyramid_edge_p(0.1, 0.1, 0.3, &pe, &mut v, &mut nb, &EDGES);
            assert_eq!(nb, 5 + ne);
            let mut nb = nodal(0.1, 0.1, 0.3, &mut v);
            pyramid_face_p(0.1, 0.1, 0.3, &pf, &mut v, &mut nb, &FACES);
            assert_eq!(nb, 5 + nf);
            let mut nb = nodal(0.1, 0.1, 0.3, &mut v);
            pyramid_bubble_p(0.1, 0.1, 0.3, p, &mut v, &mut nb);
            assert_eq!(nb, 5 + nbub);

            check_fd(
                &|u, v, w, o| {
                    let mut nb = nodal(u, v, w, o);
                    pyramid_edge_p(u, v, w, &pe, o, &mut nb, &EDGES);
                },
                &|u, v, w, g| {
                    let mut nb = dnodal(u, v, w, g);
                    d_pyramid_edge_p(u, v, w, &pe, g, &mut nb, &EDGES);
                },
                5,
                ne,
            );
            if nf > 0 {
                check_fd(
                    &|u, v, w, o| {
                        let mut nb = nodal(u, v, w, o);
                        pyramid_face_p(u, v, w, &pf, o, &mut nb, &FACES);
                    },
                    &|u, v, w, g| {
                        let mut nb = dnodal(u, v, w, g);
                        d_pyramid_face_p(u, v, w, &pf, g, &mut nb, &FACES);
                    },
                    5,
                    nf,
                );
            }
            if nbub > 0 {
                check_fd(
                    &|u, v, w, o| {
                        let mut nb = nodal(u, v, w, o);
                        pyramid_bubble_p(u, v, w, p, o, &mut nb);
                    },
                    &|u, v, w, g| {
                        let mut nb = dnodal(u, v, w, g);
                        d_pyramid_bubble_p(u, v, w, p, g, &mut nb);
                    },
                    5,
                    nbub,
                );
            }
        }
    }
}
