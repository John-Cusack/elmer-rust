//! Tetrahedron p-element basis (`H1Basis_Tetra*` in `H1Basis.F90`).
//!
//! Scalar transliterations; see [`crate::types`] for the calling conventions.
//! The affine coordinates `L1..L4` live on Elmer's tetrahedron with the base
//! triangle of side 2 in the `(u, v)` plane (the "p" reference tetrahedron),
//! not on the unit tetrahedron of the linear nodal basis.

#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::many_single_char_names,
    clippy::excessive_precision
)]

use crate::legendre::{d_legendre_p, d_var_phi, legendre_p, var_phi};

/// `H1Basis_TetraNodalP`: the four affine nodal functions, appended at `*nb`.
pub fn tetra_nodal_p(u: f64, v: f64, w: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    let e = 1.0 / 6f64.sqrt();
    let f = 1.0 / 8f64.sqrt();
    let g = 3f64.sqrt() * f;
    val[*nb] = c * (1.0 - u - d * v - e * w);
    val[*nb + 1] = c * (1.0 + u - d * v - e * w);
    val[*nb + 2] = d * (v - f * w);
    val[*nb + 3] = g * w;
    *nb += 4;
}

/// `H1Basis_TetraL`: affine coordinate of local node `node` (1-based).
pub fn tetra_l(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    let e = 1.0 / 6f64.sqrt();
    let f = 1.0 / 8f64.sqrt();
    let g = 3f64.sqrt() * f;
    match node {
        1 => c * (1.0 - u - d * v - e * w),
        2 => c * (1.0 + u - d * v - e * w),
        3 => d * (v - f * w),
        4 => g * w,
        _ => panic!("tetra_l: node {node} outside 1..=4"),
    }
}

/// `H1Basis_dTetraNodalP`: gradients of the four affine nodal functions.
pub fn d_tetra_nodal_p(grad: &mut [[f64; 3]], nb: &mut usize) {
    let half = 1.0 / 2.0;
    let a = 1.0 / 3f64.sqrt();
    let b = -a / 2.0;
    let c = -6f64.sqrt() / 12.0;
    let d = -3.0 * c;
    grad[*nb] = [-half, b, c];
    grad[*nb + 1] = [half, b, c];
    grad[*nb + 2] = [0.0, a, c];
    grad[*nb + 3] = [0.0, 0.0, d];
    *nb += 4;
}

/// `H1Basis_dTetraL`: gradient of the affine coordinate of node `node`.
pub fn d_tetra_l(node: usize) -> [f64; 3] {
    let c = 1.0 / 2.0;
    let d = 1.0 / 3f64.sqrt();
    let e = 1.0 / 6f64.sqrt();
    let f = 1.0 / 8f64.sqrt();
    let g = 3f64.sqrt() * f;
    match node {
        1 => [-c, c * -d, c * -e],
        2 => [c * 1.0, c * -d, c * -e],
        3 => [d * 0.0, d * 1.0, d * -f],
        4 => [0.0, 0.0, g],
        _ => panic!("d_tetra_l: node {node} outside 1..=4"),
    }
}

/// `H1Basis_TetraEdgeP`: edge functions of the six edges, `pmax[i]` = degree
/// of edge `i` (`pmax-1` functions each).
pub fn tetra_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 0..6 {
        let la = tetra_l(edgedir[i][0], u, v, w);
        let lb = tetra_l(edgedir[i][1], u, v, w);
        for j in 1..=pmax[i] - 1 {
            val[*nb + j as usize - 1] = la * lb * var_phi(j as usize + 1, lb - la);
        }
        *nb = (*nb as i64 + pmax[i] as i64 - 1) as usize;
    }
}

/// `H1Basis_dTetraEdgeP`: gradients matching [`tetra_edge_p`].
pub fn d_tetra_edge_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    edgedir: &[[usize; 2]],
) {
    for i in 0..6 {
        let dla = d_tetra_l(edgedir[i][0]);
        let dlb = d_tetra_l(edgedir[i][1]);
        let la = tetra_l(edgedir[i][0], u, v, w);
        let lb = tetra_l(edgedir[i][1], u, v, w);
        for j in 1..=pmax[i] - 1 {
            let vphi = var_phi(j as usize + 1, lb - la);
            let dvphi = d_var_phi(j as usize + 1, lb - la);
            let g = &mut grad[*nb + j as usize - 1];
            for c in 0..3 {
                g[c] = dla[c] * lb * vphi + la * dlb[c] * vphi + la * lb * dvphi * (dlb[c] - dla[c]);
            }
        }
        *nb = (*nb as i64 + pmax[i] as i64 - 1) as usize;
    }
}

/// `H1Basis_TetraFaceP`: face functions of the four faces, `pmax[i]` = degree
/// of face `i`; `facedir[i][0..3]` are 1-based local nodes.
pub fn tetra_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    val: &mut [f64],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    for i in 0..4 {
        let la = tetra_l(facedir[i][0], u, v, w);
        let lb = tetra_l(facedir[i][1], u, v, w);
        let lc = tetra_l(facedir[i][2], u, v, w);
        for j in 0..=pmax[i] - 3 {
            for k in 1..=pmax[i] - j - 2 {
                val[*nb + k as usize - 1] = la
                    * lb
                    * lc
                    * legendre_p(j as usize, lb - la)
                    * legendre_p(k as usize - 1, 2.0 * lc - 1.0);
            }
            *nb += (pmax[i] - j - 2).max(0) as usize;
        }
    }
}

/// `H1Basis_dTetraFaceP`: gradients matching [`tetra_face_p`].
pub fn d_tetra_face_p(
    u: f64,
    v: f64,
    w: f64,
    pmax: &[i32],
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    facedir: &[[usize; 4]],
) {
    for i in 0..4 {
        let dla = d_tetra_l(facedir[i][0]);
        let dlb = d_tetra_l(facedir[i][1]);
        let dlc = d_tetra_l(facedir[i][2]);
        let la = tetra_l(facedir[i][0], u, v, w);
        let lb = tetra_l(facedir[i][1], u, v, w);
        let lc = tetra_l(facedir[i][2], u, v, w);
        for j in 0..=pmax[i] - 3 {
            for k in 1..=pmax[i] - j - 2 {
                let lb_la = lb - la;
                let lc_1 = 2.0 * lc - 1.0;
                let a = legendre_p(j as usize, lb_la);
                let b = legendre_p(k as usize - 1, lc_1);
                let g = &mut grad[*nb + k as usize - 1];
                for c in 0..3 {
                    g[c] = dla[c] * lb * lc * a * b
                        + la * dlb[c] * lc * a * b
                        + la * lb * dlc[c] * a * b
                        + la * lb * lc * d_legendre_p(j as usize, lb_la) * (dlb[c] - dla[c]) * b
                        + la * lb * lc * a * d_legendre_p(k as usize - 1, lc_1) * 2.0 * dlc[c];
                }
            }
            *nb += (pmax[i] - j - 2).max(0) as usize;
        }
    }
}

/// `H1Basis_TetraBubbleP`: interior bubble functions of degree `pmax`.
pub fn tetra_bubble_p(u: f64, v: f64, w: f64, pmax: i32, val: &mut [f64], nb: &mut usize) {
    let l1 = tetra_l(1, u, v, w);
    let l2 = tetra_l(2, u, v, w);
    let l3 = tetra_l(3, u, v, w);
    let l4 = tetra_l(4, u, v, w);
    for i in 0..=pmax - 4 {
        for j in 0..=pmax - i - 4 {
            for k in 1..=pmax - i - j - 3 {
                let l2_l1 = l2 - l1;
                let l3_1 = 2.0 * l3 - 1.0;
                let l4_1 = 2.0 * l4 - 1.0;
                val[*nb + k as usize - 1] = l1
                    * l2
                    * l3
                    * l4
                    * legendre_p(i as usize, l2_l1)
                    * legendre_p(j as usize, l3_1)
                    * legendre_p(k as usize - 1, l4_1);
            }
            *nb += (pmax - i - j - 3).max(0) as usize;
        }
    }
}

/// `H1Basis_dTetraBubbleP`: gradients matching [`tetra_bubble_p`].
pub fn d_tetra_bubble_p(u: f64, v: f64, w: f64, pmax: i32, grad: &mut [[f64; 3]], nb: &mut usize) {
    let sq3_6 = 3f64.sqrt() / 6.0;
    let sq3_3 = 3f64.sqrt() / 3.0;
    let sq6_12 = 6f64.sqrt() / 12.0;
    let sq6_6 = 6f64.sqrt() / 6.0;
    let sq6_4 = 6f64.sqrt() / 4.0;
    let sq6_2 = 6f64.sqrt() / 2.0;
    let l1 = tetra_l(1, u, v, w);
    let l2 = tetra_l(2, u, v, w);
    let l3 = tetra_l(3, u, v, w);
    let l4 = tetra_l(4, u, v, w);
    for i in 0..=pmax - 4 {
        for j in 0..=pmax - i - 4 {
            for k in 1..=pmax - i - j - 3 {
                let l2_l1 = l2 - l1;
                let l3_1 = 2.0 * l3 - 1.0;
                let l4_1 = 2.0 * l4 - 1.0;
                let a = legendre_p(i as usize, l2_l1);
                let b = legendre_p(j as usize, l3_1);
                let c = legendre_p(k as usize - 1, l4_1);
                let g = &mut grad[*nb + k as usize - 1];
                g[0] = -1.0 / 2.0 * l2 * l3 * l4 * a * b * c
                    + 1.0 / 2.0 * l1 * l3 * l4 * a * b * c
                    + l1 * l2 * l3 * l4 * d_legendre_p(i as usize, l2_l1) * b * c;
                g[1] = -sq3_6 * l2 * l3 * l4 * a * b * c - sq3_6 * l1 * l3 * l4 * a * b * c
                    + sq3_3 * l1 * l2 * l4 * a * b * c
                    + 2.0 * sq3_3 * l1 * l2 * l3 * l4 * a * d_legendre_p(j as usize, l3_1) * c;
                g[2] = -sq6_12 * l2 * l3 * l4 * a * b * c
                    - sq6_12 * l1 * l3 * l4 * a * b * c
                    - sq6_12 * l1 * l2 * l4 * a * b * c
                    + sq6_4 * l1 * l2 * l3 * a * b * c
                    - sq6_6 * l1 * l2 * l3 * l4 * a * d_legendre_p(j as usize, l3_1) * c
                    + sq6_2 * l1 * l2 * l3 * l4 * a * b * d_legendre_p(k as usize - 1, l4_1);
            }
            *nb += (pmax - i - j - 3).max(0) as usize;
        }
    }
}

/// `H1Basis_TetraNodal` (linear-only, "not a barycentric tetra"): values at
/// `val[0..4]` of the unit-tetrahedron barycentric functions.
pub fn tetra_nodal(u: f64, v: f64, w: f64, val: &mut [f64]) {
    val[0] = 1.0 - u - v - w;
    val[1] = u;
    val[2] = v;
    val[3] = w;
}

/// `H1Basis_dTetraNodal`: constant gradients of [`tetra_nodal`] into
/// `grad[0..4]`.
pub fn d_tetra_nodal(grad: &mut [[f64; 3]]) {
    grad[0] = [-1.0, -1.0, -1.0];
    grad[1] = [1.0, 0.0, 0.0];
    grad[2] = [0.0, 1.0, 0.0];
    grad[3] = [0.0, 0.0, 1.0];
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `H1Basis_GetTetraEdgeDirection` type 1.
    const EDGES1: [[usize; 2]; 6] = [[1, 2], [2, 3], [1, 3], [1, 4], [2, 4], [3, 4]];
    /// `H1Basis_GetTetraFaceDirection` type 1.
    const FACES1: [[usize; 4]; 4] = [[1, 2, 3, 0], [1, 2, 4, 0], [2, 3, 4, 0], [1, 3, 4, 0]];
    /// `H1Basis_GetTetraEdgeDirection` type 2.
    const EDGES2: [[usize; 2]; 6] = [[1, 2], [3, 2], [1, 3], [1, 4], [2, 4], [3, 4]];
    /// `H1Basis_GetTetraFaceDirection` type 2.
    const FACES2: [[usize; 4]; 4] = [[1, 3, 2, 0], [1, 2, 4, 0], [3, 2, 4, 0], [1, 3, 4, 0]];

    /// Points inside the reference tetrahedron (L1..L4 > 0).
    fn points() -> Vec<[f64; 3]> {
        let mut pts = Vec::new();
        // pick barycentric combos and map: x = sum L_i X_i with the
        // reference vertices (-1,0,0),(1,0,0),(0,sqrt3,0),(0,sqrt3/3,2sqrt6/3)
        let verts = [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 3f64.sqrt(), 0.0],
            [0.0, 3f64.sqrt() / 3.0, 2.0 * 6f64.sqrt() / 3.0],
        ];
        for bary in [
            [0.25, 0.25, 0.25, 0.25],
            [0.4, 0.3, 0.2, 0.1],
            [0.1, 0.2, 0.3, 0.4],
            [0.35, 0.15, 0.4, 0.1],
        ] {
            let mut p = [0.0; 3];
            for (b, vtx) in bary.iter().zip(verts.iter()) {
                for c in 0..3 {
                    p[c] += b * vtx[c];
                }
            }
            pts.push(p);
        }
        pts
    }

    type ValFn<'a> = &'a dyn Fn(f64, f64, f64, &mut Vec<f64>);
    type GradFn<'a> = &'a dyn Fn(f64, f64, f64, &mut Vec<[f64; 3]>);

    fn check_fd(
        f: ValFn,
        df: GradFn,
        n: usize,
    ) {
        let h = 1e-6;
        for p in points() {
            let mut g = vec![[0.0; 3]; n + 4];
            df(p[0], p[1], p[2], &mut g);
            for c in 0..3 {
                let mut pp = p;
                let mut pm = p;
                pp[c] += h;
                pm[c] -= h;
                let mut vp = vec![0.0; n + 4];
                let mut vm = vec![0.0; n + 4];
                f(pp[0], pp[1], pp[2], &mut vp);
                f(pm[0], pm[1], pm[2], &mut vm);
                for m in 0..n {
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
    fn nodal_partition_of_unity() {
        for p in points() {
            let mut v = [0.0; 4];
            let mut nb = 0;
            tetra_nodal_p(p[0], p[1], p[2], &mut v, &mut nb);
            assert_eq!(nb, 4);
            assert!((v.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let mut g = [[0.0; 3]; 4];
            let mut nb = 0;
            d_tetra_nodal_p(&mut g, &mut nb);
            for c in 0..3 {
                assert!(g.iter().map(|x| x[c]).sum::<f64>().abs() < 1e-14);
            }
            for n in 1..=4 {
                assert_eq!(tetra_l(n, p[0], p[1], p[2]), v[n - 1]);
            }
        }
        let mut v = [0.0; 4];
        tetra_nodal(0.1, 0.2, 0.3, &mut v);
        assert!((v.iter().sum::<f64>() - 1.0).abs() < 1e-15);
        let mut g = [[0.0; 3]; 4];
        d_tetra_nodal(&mut g);
        for c in 0..3 {
            assert_eq!(g.iter().map(|x| x[c]).sum::<f64>(), 0.0);
        }
    }

    #[test]
    fn nodal_gradients_match_finite_differences() {
        check_fd(
            &|u, v, w, o| {
                let mut nb = 0;
                tetra_nodal_p(u, v, w, o, &mut nb);
            },
            &|_, _, _, g| {
                let mut nb = 0;
                d_tetra_nodal_p(g, &mut nb);
            },
            4,
        );
        for node in 1..=4 {
            let h = 1e-6;
            let (u, v, w) = (0.1, 0.2, 0.3);
            let g = d_tetra_l(node);
            let fd = [
                (tetra_l(node, u + h, v, w) - tetra_l(node, u - h, v, w)) / (2.0 * h),
                (tetra_l(node, u, v + h, w) - tetra_l(node, u, v - h, w)) / (2.0 * h),
                (tetra_l(node, u, v, w + h) - tetra_l(node, u, v, w - h)) / (2.0 * h),
            ];
            for c in 0..3 {
                assert!((fd[c] - g[c]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn edge_face_bubble_gradients_and_counts() {
        for p in 2..=5i32 {
            for (edges, faces) in [(&EDGES1, &FACES1), (&EDGES2, &FACES2)] {
                let pe = [p, p, p, p, p, p];
                let ne = 6 * (p as usize - 1);
                let nf = 4 * (p as usize - 2) * (p as usize - 1) / 2;
                let nbub = (0..=p - 4)
                    .map(|i| (0..=p - i - 4).map(|j| (p - i - j - 3).max(0)).sum::<i32>())
                    .sum::<i32>() as usize;
                // count checks (ElemInfo nbq formulae)
                let mut nb = 0;
                let mut v = vec![0.0; 100];
                tetra_edge_p(0.1, 0.2, 0.3, &pe, &mut v, &mut nb, edges);
                assert_eq!(nb, ne);
                let mut nb = 0;
                tetra_face_p(0.1, 0.2, 0.3, &[p; 4], &mut v, &mut nb, faces);
                let mut nbq = 0;
                for _ in 0..4 {
                    for j in 0..=p {
                        nbq += (p - j - 2).max(0) as usize;
                    }
                }
                assert_eq!(nb, nbq);
                assert_eq!(nb, nf);
                let mut nb = 0;
                tetra_bubble_p(0.1, 0.2, 0.3, p, &mut v, &mut nb);
                assert_eq!(nb, nbub);

                check_fd(
                    &|u, v, w, o| {
                        let mut nb = 0;
                        tetra_edge_p(u, v, w, &pe, o, &mut nb, edges);
                    },
                    &|u, v, w, g| {
                        let mut nb = 0;
                        d_tetra_edge_p(u, v, w, &pe, g, &mut nb, edges);
                    },
                    ne,
                );
                if nf > 0 {
                    check_fd(
                        &|u, v, w, o| {
                            let mut nb = 0;
                            tetra_face_p(u, v, w, &[p; 4], o, &mut nb, faces);
                        },
                        &|u, v, w, g| {
                            let mut nb = 0;
                            d_tetra_face_p(u, v, w, &[p; 4], g, &mut nb, faces);
                        },
                        nf,
                    );
                }
                if nbub > 0 {
                    check_fd(
                        &|u, v, w, o| {
                            let mut nb = 0;
                            tetra_bubble_p(u, v, w, p, o, &mut nb);
                        },
                        &|u, v, w, g| {
                            let mut nb = 0;
                            d_tetra_bubble_p(u, v, w, p, g, &mut nb);
                        },
                        nbub,
                    );
                }
            }
        }
    }
}
