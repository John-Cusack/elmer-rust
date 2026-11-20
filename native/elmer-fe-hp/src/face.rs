//! Face (div-conforming) element bases: `FaceElementInfo` (ElemInfo.F90:2747),
//! `PiolaTransformationData` and `FaceElementOrientation` for the cases
//! `ModelMixedPoisson` uses on triangles: the Raviart-Thomas space RT(0)
//! (`BasisDegree` absent) and the BDM(1) space (`BDM = .TRUE.`), on the p
//! reference triangle, optionally Piola-transformed to the physical element.
//!
//! Every expression keeps the Fortran's operand order and grouping, so the
//! values are those the pristine build computes.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

/// `TriangleNodalPBasis(node, u, v)` (PElementBase.F90): the scalar nodal
/// basis of the p reference triangle (`FaceElementInfo` uses this, not the
/// `H1Basis_*` variant of [`crate::triangle`]).
pub fn triangle_nodal_p_basis(node: usize, u: f64, v: f64) -> f64 {
    let s3 = 3f64.sqrt();
    match node {
        1 => (1.0 - u - v / s3) / 2.0,
        2 => (1.0 + u - v / s3) / 2.0,
        3 => v / s3,
        _ => 0.0,
    }
}

/// `dTriangleNodalPBasis(node, u, v)`.
pub fn d_triangle_nodal_p_basis(node: usize) -> [f64; 2] {
    let s3 = 3f64.sqrt();
    match node {
        1 => [-1.0 / 2.0, -s3 / 6.0],
        2 => [1.0 / 2.0, -s3 / 6.0],
        3 => [0.0, s3 / 3.0],
        _ => [0.0, 0.0],
    }
}

/// `GetEdgeMap(3)`, 1-based local nodes of each triangle edge.
pub const TRIANGLE_EDGE_MAP: [[usize; 2]; 3] = [[1, 2], [2, 3], [3, 1]];

/// `FaceElementOrientation(Element, ReverseSign)` for a triangle: the
/// face (edge) dof of edge `q` is reversed when its second node has the
/// smaller global index. `g` are the element's global node indices.
pub fn triangle_face_orientation(g: &[u32; 3]) -> [bool; 3] {
    let mut r = [false; 3];
    for q in 0..3 {
        let [a, b] = TRIANGLE_EDGE_MAP[q];
        r[q] = g[b - 1] < g[a - 1];
    }
    r
}

/// What `FaceElementInfo` returns.
#[derive(Clone, Debug, Default)]
pub struct FaceInfo {
    /// The nodal basis `Basis(1:n)`.
    pub basis: Vec<f64>,
    /// `FBasis(1:DOFs, 1:3)`.
    pub fbasis: Vec<[f64; 3]>,
    /// `DivFBasis(1:DOFs)`.
    pub div: Vec<f64>,
    /// `detF` (its absolute value after the Piola transform).
    pub det_f: f64,
    /// `F`, the gradient of the element map.
    pub f: [[f64; 3]; 3],
}

/// `FaceElementInfo(Element, Nodes, u, v, w, detF=, Basis=, FBasis=,
/// DivFBasis=, BDM=bdm, ApplyPiolaTransform=piola)` for a triangle (303)
/// with nodes `x`, `y`, `z` and global node indices `g`.
#[allow(clippy::too_many_arguments)]
pub fn face_element_info_triangle(x: &[f64; 3], y: &[f64; 3], z: &[f64; 3], g: &[u32; 3], u: f64, v: f64, bdm: bool, piola: bool) -> FaceInfo {
    let n = 3;
    let dim = 2;
    let mut basis = vec![0.0; n];
    let mut dl = [[0.0f64; 3]; 3];
    for q in 0..n {
        basis[q] = triangle_nodal_p_basis(q + 1, u, v);
        let d = d_triangle_nodal_p_basis(q + 1);
        dl[q][0] = d[0];
        dl[q][1] = d[1];
    }
    // PiolaTransformationData: F(k,i) = SUM( x(1:n) * dLBasisdx(1:n,i) ).
    let mut lf = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for (k, c) in [x, y, z].into_iter().enumerate() {
            let mut s = 0.0;
            for q in 0..n {
                s += c[q] * dl[q][i];
            }
            lf[k][i] = s;
        }
    }
    let mut det_f = lf[0][0] * lf[1][1] - lf[0][1] * lf[1][0];

    let reverse = triangle_face_orientation(g);
    let s3 = 3f64.sqrt();
    let (mut fb, mut divb): (Vec<[f64; 3]>, Vec<f64>);
    if bdm {
        fb = vec![[0.0; 3]; 6];
        divb = vec![0.0; 6];
        // Two basis functions on face 12.
        fb[0][0] = s3 / 6.0 * (-s3 + u + v);
        fb[0][1] = s3 / 6.0 * (-s3 + 3.0 * u + v);
        divb[0] = s3 / 3.0;
        fb[1][0] = s3 / 6.0 * (s3 + u - v);
        fb[1][1] = s3 / 6.0 * (-s3 - 3.0 * u + v);
        divb[1] = s3 / 3.0;
        // Two on face 23.
        fb[2][0] = 1.0 / (3.0 + s3) * (2.0 + s3 + (2.0 + s3) * u - (1.0 + s3) * v);
        fb[2][1] = 1.0 / 6.0 * (-3.0 + s3) * v;
        divb[2] = s3 / 3.0;
        fb[3][0] = 1.0 / 6.0 * (-3.0 + s3 + (-3.0 + s3) * u + 2.0 * s3 * v);
        fb[3][1] = 1.0 / 6.0 * (3.0 + s3) * v;
        divb[3] = s3 / 3.0;
        // Two on face 31.
        fb[4][0] = 1.0 / (3.0 + s3) * (1.0 - u - v - s3 * v);
        fb[4][1] = (3.0 + 2.0 * s3) * v / (3.0 * (1.0 + s3));
        divb[4] = s3 / 3.0;
        fb[5][0] = 1.0 / 6.0 * (-3.0 - s3 + (3.0 + s3) * u + 2.0 * s3 * v);
        fb[5][1] = 1.0 / 6.0 * (-3.0 + s3) * v;
        divb[5] = s3 / 3.0;
        // Reordering and sign reversion.
        for q in 0..3 {
            if reverse[q] {
                let (w1, w2) = (fb[2 * q], fb[2 * q + 1]);
                let (d1, d2) = (divb[2 * q], divb[2 * q + 1]);
                fb[2 * q] = [-w2[0], -w2[1], 0.0];
                divb[2 * q] = -d2;
                fb[2 * q + 1] = [-w1[0], -w1[1], 0.0];
                divb[2 * q + 1] = -d1;
            }
        }
    } else {
        fb = vec![[0.0; 3]; 3];
        divb = vec![0.0; 3];
        fb[0][0] = s3 / 6.0 * u;
        fb[0][1] = -0.5 + s3 / 6.0 * v;
        divb[0] = s3 / 3.0;
        fb[1][0] = s3 / 6.0 * (1.0 + u);
        fb[1][1] = s3 / 6.0 * v;
        divb[1] = s3 / 3.0;
        fb[2][0] = s3 / 6.0 * (-1.0 + u);
        fb[2][1] = s3 / 6.0 * v;
        divb[2] = s3 / 3.0;
        for q in 0..3 {
            if reverse[q] {
                fb[q] = [-fb[q][0], -fb[q][1], -fb[q][2]];
                divb[q] = -divb[q];
            }
        }
    }

    if piola {
        for j in 0..fb.len() {
            let mut w = [0.0f64; 3];
            for k in 0..dim {
                let mut s = 0.0;
                for i in 0..dim {
                    s += lf[k][i] * fb[j][i];
                }
                w[k] = s;
            }
            for k in 0..dim {
                fb[j][k] = 1.0 / det_f * w[k];
            }
            divb[j] *= 1.0 / det_f;
        }
        det_f = det_f.abs();
    }
    FaceInfo { basis, fbasis: fb, div: divb, det_f, f: lf }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The RT(0) functions have unit outward flux through their own edge
    /// and none through the others (the defining degrees of freedom), on a
    /// physical triangle after the Piola transform.
    #[test]
    fn rt0_fluxes_are_the_dofs() {
        let (x, y, z) = ([0.0, 2.0, 0.5], [0.0, 0.3, 1.5], [0.0; 3]);
        let g = [4, 9, 1]; // edge 3 (3 -> 1) and 1 (1 -> 2) keep, edge 2 reverses
        let r = triangle_face_orientation(&g);
        assert_eq!(r, [false, true, false]);
        // Integrate the normal flux along each edge with the 3-point
        // Gauss rule on the reference edge.
        let pts = [(-0.774_596_669_241_483_4, 5.0 / 9.0), (0.0, 8.0 / 9.0), (0.774_596_669_241_483_4, 5.0 / 9.0)];
        let (pu, pv) = ([-1.0, 1.0, 0.0], [0.0, 0.0, 3f64.sqrt()]);
        for (q, [a, b]) in TRIANGLE_EDGE_MAP.iter().enumerate() {
            let (a, b) = (a - 1, b - 1);
            let (tx, ty) = (x[b] - x[a], y[b] - y[a]);
            let normal = [ty, -tx]; // outward for counter-clockwise nodes, length = edge length
            for k in 0..3 {
                let mut flux = 0.0;
                for (t, w) in pts {
                    let s = (t + 1.0) / 2.0;
                    let (u, v) = (pu[a] + s * (pu[b] - pu[a]), pv[a] + s * (pv[b] - pv[a]));
                    let fi = face_element_info_triangle(&x, &y, &z, &g, u, v, false, true);
                    flux += w / 2.0 * (fi.fbasis[k][0] * normal[0] + fi.fbasis[k][1] * normal[1]);
                }
                let want = if k == q { if r[q] { -1.0 } else { 1.0 } } else { 0.0 };
                assert!((flux - want).abs() < 1e-14, "edge {q} function {k}: {flux}");
            }
        }
    }

    /// The divergence is constant and integrates to the edge flux sum
    /// (`div` of each RT(0) function times the area equals its total flux).
    #[test]
    fn divergence_matches_the_flux() {
        let (x, y, z) = ([0.0, 2.0, 0.5], [0.0, 0.3, 1.5], [0.0; 3]);
        let fi = face_element_info_triangle(&x, &y, &z, &[1, 2, 3], 0.1, 0.2, false, true);
        let area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
        let r = triangle_face_orientation(&[1, 2, 3]);
        assert_eq!(r, [false, false, true]);
        for k in 0..3 {
            let want = if r[k] { -1.0 } else { 1.0 };
            assert!((fi.div[k] * area - want).abs() < 1e-14, "{k}: {}", fi.div[k] * area);
        }
        let bdm = face_element_info_triangle(&x, &y, &z, &[1, 2, 3], 0.1, 0.2, true, true);
        assert_eq!(bdm.fbasis.len(), 6);
    }
}
