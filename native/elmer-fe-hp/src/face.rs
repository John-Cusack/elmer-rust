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

/// `TetraNodalPBasis(node, u, v, w)` (PElementBase.F90).
pub fn tetra_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let (s3, s6, s8) = (3f64.sqrt(), 6f64.sqrt(), 8f64.sqrt());
    match node {
        1 => 1.0 / 2.0 * (1.0 - u - v / s3 - w / s6),
        2 => 1.0 / 2.0 * (1.0 + u - v / s3 - w / s6),
        3 => s3 / 3.0 * (v - w / s8),
        4 => (3.0f64 / 8.0).sqrt() * w,
        _ => 0.0,
    }
}

/// `dTetraNodalPBasis(node, u, v, w)`.
pub fn d_tetra_nodal_p_basis(node: usize) -> [f64; 3] {
    let (s3, s6) = (3f64.sqrt(), 6f64.sqrt());
    match node {
        1 => [-1.0 / 2.0, -s3 / 6.0, -s6 / 12.0],
        2 => [1.0 / 2.0, -s3 / 6.0, -s6 / 12.0],
        3 => [0.0, s3 / 3.0, -s6 / 12.0],
        4 => [0.0, 0.0, s6 / 4.0],
        _ => [0.0; 3],
    }
}

/// The `TetraFaceMap` of `FaceElementOrientation` and
/// `FaceElementBasisOrdering`, 1-based local nodes of each face.
pub const TETRA_FACE_MAP: [[usize; 3]; 4] = [[2, 1, 3], [1, 2, 4], [2, 3, 4], [3, 1, 4]];

/// `FaceElementOrientation(Element, ReverseSign)` for a tetrahedron.
pub fn tetra_face_orientation(g: &[u32; 4]) -> [bool; 4] {
    let mut r = [false; 4];
    for q in 0..4 {
        let f = TETRA_FACE_MAP[q].map(|k| g[k - 1]);
        r[q] = if f[0] < f[1] && f[0] < f[2] {
            f[2] < f[1]
        } else if f[1] < f[0] && f[1] < f[2] {
            f[0] < f[2]
        } else {
            f[1] < f[0]
        };
    }
    r
}

/// `FaceElementBasisOrdering(Element, FDofMap)` for a tetrahedron: the
/// 1-based position, in the default table, of each of the three dofs of
/// each face.
pub fn tetra_face_basis_ordering(g: &[u32; 4]) -> [[usize; 3]; 4] {
    let mut m = [[0usize; 3]; 4];
    for q in 0..4 {
        let f = TETRA_FACE_MAP[q].map(|k| g[k - 1]);
        m[q] = if f[0] < f[1] && f[0] < f[2] {
            if f[1] < f[2] { [1, 2, 3] } else { [1, 3, 2] }
        } else if f[1] < f[0] && f[1] < f[2] {
            if f[0] < f[2] { [2, 1, 3] } else { [2, 3, 1] }
        } else if f[0] < f[1] {
            [3, 1, 2]
        } else {
            [3, 2, 1]
        };
    }
    m
}

/// `FaceElementInfo` for a tetrahedron (504): RT(0) or, with `bdm`, the
/// BDM(1) space of the first kind's twelve functions, as
/// [`face_element_info_triangle`].
#[allow(clippy::too_many_arguments)]
pub fn face_element_info_tetra(x: &[f64; 4], y: &[f64; 4], z: &[f64; 4], g: &[u32; 4], u: f64, v: f64, w: f64, bdm: bool, piola: bool) -> FaceInfo {
    let n = 4;
    let dim = 3;
    let mut basis = vec![0.0; n];
    let mut dl = [[0.0f64; 3]; 4];
    for q in 0..n {
        basis[q] = tetra_nodal_p_basis(q + 1, u, v, w);
        dl[q] = d_tetra_nodal_p_basis(q + 1);
    }
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
    let f = &lf;
    let mut det_f = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
        + f[0][1] * (f[1][2] * f[2][0] - f[1][0] * f[2][2])
        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);

    let reverse = tetra_face_orientation(g);
    let (s2, s3, s6) = (2f64.sqrt(), 3f64.sqrt(), 6f64.sqrt());
    let (mut fb, mut divb): (Vec<[f64; 3]>, Vec<f64>);
    if bdm {
        // The default table; `5*Sqrt(6.0)` in WorkBasis(8,2) is single
        // precision in the Fortran.
        let c8 = (5.0f32 * 6.0f32.sqrt()) as f64;
        let wb: [[f64; 3]; 12] = [
            [
                (3.0 * s6 + 2.0 * s6 * u - 3.0 * s2 * v - 3.0 * w) / 12.0,
                (-2.0 * s2 - 3.0 * s2 * u + s3 * w) / 12.0,
                (-8.0 - 12.0 * u + 4.0 * s3 * v + s6 * w) / 12.0,
            ],
            [
                (2.0 * s6 * u + 3.0 * (-s6 + s2 * v + w)) / 12.0,
                (-2.0 * s2 + 3.0 * s2 * u + s3 * w) / 12.0,
                u + (-8.0 + 4.0 * s3 * v + s6 * w) / 12.0,
            ],
            [-u / (2.0 * s6), (s2 + 3.0 * s6 * v - 2.0 * s3 * w) / 12.0, (4.0 - 8.0 * s3 * v + s6 * w) / 12.0],
            [
                (2.0 * s6 * u + 3.0 * (-s6 + s2 * v + w)) / 12.0,
                (-6.0 * s2 + 9.0 * s2 * u + 2.0 * s6 * v + 3.0 * s3 * w) / 12.0,
                -w / (2.0 * s6),
            ],
            [
                (3.0 * s6 + 2.0 * s6 * u - 3.0 * s2 * v - 3.0 * w) / 12.0,
                (-6.0 * s2 - 9.0 * s2 * u + 2.0 * s6 * v + 3.0 * s3 * w) / 12.0,
                -w / (2.0 * s6),
            ],
            [-u / (2.0 * s6), (3.0 * s2 - s6 * v - 6.0 * s3 * w) / 12.0, (5.0 * w) / (2.0 * s6)],
            [(5.0 * s6 + 5.0 * s6 * u - 6.0 * s2 * v - 6.0 * w) / 12.0, -v / (2.0 * s6), -w / (2.0 * s6)],
            [(-s6 - s6 * u + 6.0 * s2 * v - 3.0 * w) / 12.0, (c8 * v - 3.0 * s3 * w) / 12.0, -w / (2.0 * s6)],
            [(-s6 - s6 * u + 9.0 * w) / 12.0, (-(s6 * v) + 3.0 * s3 * w) / 12.0, (5.0 * w) / (2.0 * s6)],
            [(s6 - s6 * u - 6.0 * s2 * v + 3.0 * w) / 12.0, (5.0 * s6 * v - 3.0 * s3 * w) / 12.0, -w / (2.0 * s6)],
            [(-5.0 * s6 + 5.0 * s6 * u + 6.0 * s2 * v + 6.0 * w) / 12.0, -v / (2.0 * s6), -w / (2.0 * s6)],
            [(s6 - s6 * u - 9.0 * w) / 12.0, (-(s6 * v) + 3.0 * s3 * w) / 12.0, (5.0 * w) / (2.0 * s6)],
        ];
        let map = tetra_face_basis_ordering(g);
        fb = vec![[0.0; 3]; 12];
        divb = vec![0.0; 12];
        for q in 0..4 {
            let sgn = if reverse[q] { -1.0 } else { 1.0 };
            for j in 0..3 {
                let k = map[q][j];
                let i = q * 3 + j;
                let src = wb[q * 3 + k - 1];
                fb[i] = [sgn * src[0], sgn * src[1], sgn * src[2]];
                divb[i] = sgn * s3 / (2.0 * s2);
            }
        }
    } else {
        let a = s2 / 4.0;
        fb = vec![
            [a * u, -s6 / 12.0 + a * v, -1.0 / s3 + a * w],
            [a * u, -s6 / 4.0 + a * v, a * w],
            [a + a * u, a * v, a * w],
            [-s2 / 4.0 + a * u, a * v, a * w],
        ];
        divb = vec![3.0 * s2 / 4.0; 4];
        for q in 0..4 {
            if reverse[q] {
                fb[q] = [-fb[q][0], -fb[q][1], -fb[q][2]];
                divb[q] = -divb[q];
            }
        }
    }
    if piola {
        for j in 0..fb.len() {
            let mut wv = [0.0f64; 3];
            for k in 0..dim {
                let mut s = 0.0;
                for i in 0..dim {
                    s += lf[k][i] * fb[j][i];
                }
                wv[k] = s;
            }
            for k in 0..dim {
                fb[j][k] = 1.0 / det_f * wv[k];
            }
            divb[j] *= 1.0 / det_f;
        }
        det_f = det_f.abs();
    }
    FaceInfo { basis, fbasis: fb, div: divb, det_f, f: lf }
}

/// `QuadNodalPBasis(node, u, v)` (PElementBase.F90).
pub fn quad_nodal_p_basis(node: usize, u: f64, v: f64) -> f64 {
    match node {
        1 => (1.0 - u) * (1.0 - v) / 4.0,
        2 => (1.0 + u) * (1.0 - v) / 4.0,
        3 => (1.0 + u) * (1.0 + v) / 4.0,
        4 => (1.0 - u) * (1.0 + v) / 4.0,
        _ => 0.0,
    }
}

/// `LineNodalPBasis(node, u)`.
pub fn line_nodal_p_basis(node: usize, u: f64) -> f64 {
    match node {
        1 => (1.0 - u) / 2.0,
        2 => (1.0 + u) / 2.0,
        _ => 0.0,
    }
}

/// `dLineNodalPBasis(node, u)`.
pub fn d_line_nodal_p_basis(node: usize) -> f64 {
    match node {
        1 => -1.0 / 2.0,
        2 => 1.0 / 2.0,
        _ => 0.0,
    }
}

/// `BrickNodalPBasis(node, u, v, w)`.
pub fn brick_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let value = match node {
        1 => (1.0 - u) * (1.0 - v) * (1.0 - w),
        2 => (1.0 + u) * (1.0 - v) * (1.0 - w),
        3 => (1.0 + u) * (1.0 + v) * (1.0 - w),
        4 => (1.0 - u) * (1.0 + v) * (1.0 - w),
        5 => (1.0 - u) * (1.0 - v) * (1.0 + w),
        6 => (1.0 + u) * (1.0 - v) * (1.0 + w),
        7 => (1.0 + u) * (1.0 + v) * (1.0 + w),
        8 => (1.0 - u) * (1.0 + v) * (1.0 + w),
        _ => 0.0,
    };
    value / 8.0
}

/// `dBrickNodalPBasis(node, u, v, w)`.
pub fn d_brick_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> [f64; 3] {
    let g = match node {
        1 => [-((1.0 - v) * (1.0 - w)), -((1.0 - u) * (1.0 - w)), -((1.0 - u) * (1.0 - v))],
        2 => [(1.0 - v) * (1.0 - w), -((1.0 + u) * (1.0 - w)), -((1.0 + u) * (1.0 - v))],
        3 => [(1.0 + v) * (1.0 - w), (1.0 + u) * (1.0 - w), -((1.0 + u) * (1.0 + v))],
        4 => [-((1.0 + v) * (1.0 - w)), (1.0 - u) * (1.0 - w), -((1.0 - u) * (1.0 + v))],
        5 => [-((1.0 - v) * (1.0 + w)), -((1.0 - u) * (1.0 + w)), (1.0 - u) * (1.0 - v)],
        6 => [(1.0 - v) * (1.0 + w), -((1.0 + u) * (1.0 + w)), (1.0 + u) * (1.0 - v)],
        7 => [(1.0 + v) * (1.0 + w), (1.0 + u) * (1.0 + w), (1.0 + u) * (1.0 + v)],
        8 => [-((1.0 + v) * (1.0 + w)), (1.0 - u) * (1.0 + w), (1.0 - u) * (1.0 + v)],
        _ => [0.0; 3],
    };
    [g[0] / 8.0, g[1] / 8.0, g[2] / 8.0]
}

/// The `BrickFaceMap` of `FaceElementOrientation` and
/// `FaceElementBasisOrdering`.
pub const BRICK_FACE_MAP: [[usize; 4]; 6] = [[2, 1, 4, 3], [5, 6, 7, 8], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]];

/// The 1-based corner of a square face with the smallest global index, as
/// `SquareFaceDofsOrdering` and `FaceElementBasisOrdering` find it.
fn square_first(ind: &[u32; 4]) -> usize {
    let mut k = if ind[0] < ind[1] { 1 } else { 2 };
    let l = if ind[3] < ind[2] { 4 } else { 3 };
    if ind[k - 1] > ind[l - 1] {
        k = l;
    }
    k
}

/// `FaceElementOrientation(Element, ReverseSign)` for a brick: the
/// `ReverseSign` output of `SquareFaceDofsOrdering` for each face.
pub fn brick_face_orientation(g: &[u32; 8]) -> [bool; 6] {
    let mut r = [false; 6];
    for q in 0..6 {
        let ind = BRICK_FACE_MAP[q].map(|k| g[k - 1]);
        r[q] = match square_first(&ind) {
            1 => ind[1] >= ind[3],
            2 => ind[2] >= ind[0],
            3 => ind[3] >= ind[1],
            _ => ind[0] >= ind[2],
        };
    }
    r
}

/// `FaceElementBasisOrdering(Element, FDofMap)` for a brick.
pub fn brick_face_basis_ordering(g: &[u32; 8]) -> [[usize; 4]; 6] {
    let mut m = [[0usize; 4]; 6];
    for q in 0..6 {
        let ind = BRICK_FACE_MAP[q].map(|k| g[k - 1]);
        let f = &mut m[q];
        match square_first(&ind) {
            1 => {
                (f[0], f[2]) = (1, 3);
                (f[1], f[3]) = if ind[1] < ind[3] { (2, 4) } else { (4, 2) };
            }
            2 => {
                (f[1], f[3]) = (1, 3);
                (f[0], f[2]) = if ind[2] < ind[0] { (4, 2) } else { (2, 4) };
            }
            3 => {
                (f[2], f[0]) = (1, 3);
                (f[1], f[3]) = if ind[3] < ind[1] { (4, 2) } else { (2, 4) };
            }
            _ => {
                (f[3], f[1]) = (1, 3);
                (f[0], f[2]) = if ind[0] < ind[2] { (2, 4) } else { (4, 2) };
            }
        }
    }
    m
}

/// `FaceElementInfo` for a brick (808): the 48-function space of four
/// dofs per face and 24 interior functions.
#[allow(clippy::too_many_arguments)]
pub fn face_element_info_brick(x: &[f64; 8], y: &[f64; 8], z: &[f64; 8], g: &[u32; 8], u: f64, v: f64, w: f64, piola: bool) -> FaceInfo {
    let n = 8;
    let dim = 3;
    let mut basis = vec![0.0; n];
    let mut dl = [[0.0f64; 3]; 8];
    for q in 0..n {
        basis[q] = brick_nodal_p_basis(q + 1, u, v, w);
        dl[q] = d_brick_nodal_p_basis(q + 1, u, v, w);
    }
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
    let f = &lf;
    let mut det_f = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
        + f[0][1] * (f[1][2] * f[2][0] - f[1][0] * f[2][2])
        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);

    let reverse = brick_face_orientation(g);
    let qb = quad_nodal_p_basis;
    let (l1, l2, d1, d2) = (line_nodal_p_basis(1, u), line_nodal_p_basis(2, u), d_line_nodal_p_basis(1), d_line_nodal_p_basis(2));
    let lv = (line_nodal_p_basis(1, v), line_nodal_p_basis(2, v));
    let lw = (line_nodal_p_basis(1, w), line_nodal_p_basis(2, w));
    let sq = [2usize, 1, 4, 3];
    let mut wb = [[0.0f64; 3]; 24];
    let mut wd = [0.0f64; 24];
    for q in 0..4 {
        // Face 2143.
        wb[q][2] = -qb(sq[q], u, v) * lw.0;
        wd[q] = -qb(sq[q], u, v) * d1;
        // Face 5678.
        wb[4 + q][2] = qb(q + 1, u, v) * lw.1;
        wd[4 + q] = qb(q + 1, u, v) * d2;
        // Face 1265.
        wb[8 + q][1] = -qb(q + 1, u, w) * lv.0;
        wd[8 + q] = -qb(q + 1, u, w) * d1;
        // Face 2376.
        wb[12 + q][0] = qb(q + 1, v, w) * l2;
        wd[12 + q] = qb(q + 1, v, w) * d2;
        // Face 3487.
        wb[16 + q][1] = qb(sq[q], u, w) * lv.1;
        wd[16 + q] = qb(sq[q], u, w) * d2;
        // Face 4152.
        wb[20 + q][0] = -qb(sq[q], v, w) * l1;
        wd[20 + q] = -qb(sq[q], v, w) * d1;
    }
    let map = brick_face_basis_ordering(g);
    let mut fb = vec![[0.0f64; 3]; 48];
    let mut divb = vec![0.0f64; 48];
    for q in 0..6 {
        let s = if reverse[q] { -1.0 } else { 1.0 };
        for j in 0..4 {
            let k = map[q][j];
            let i = q * 4 + j;
            let src = wb[q * 4 + k - 1];
            fb[i] = [s * src[0], s * src[1], s * src[2]];
            divb[i] = s * wd[q * 4 + k - 1];
        }
    }
    // 24 interior functions, 8 per coordinate direction.
    let mut k = 24;
    for (c, t, (a1, a2), (p, q2)) in [(0usize, u, (l1, l2), (1usize, 2usize)), (1, v, lv, (0, 2)), (2, w, lw, (0, 1))] {
        let uvw = [u, v, w];
        for j in 0..2 {
            let (fun, dfun) = if j == 0 { (1.0, 0.0) } else { (2.0 * t, 2.0) };
            for q in 0..4 {
                let qn = qb(q + 1, uvw[p], uvw[q2]);
                fb[k][c] = qn * a1 * a2 * fun;
                divb[k] = qn * (d1 * a2 * fun + a1 * d2 * fun + a1 * a2 * dfun);
                k += 1;
            }
        }
    }
    if piola {
        for j in 0..fb.len() {
            let mut wv = [0.0f64; 3];
            for kk in 0..dim {
                let mut s = 0.0;
                for i in 0..dim {
                    s += lf[kk][i] * fb[j][i];
                }
                wv[kk] = s;
            }
            for kk in 0..dim {
                fb[j][kk] = 1.0 / det_f * wv[kk];
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

    /// Tetrahedra: each RT(0) function's divergence times the volume is its
    /// unit flux (with the orientation sign); a BDM(1) function's is a third
    /// of its face's area on the p reference element (`sqrt(3)/3`), its dofs
    /// being face-Lagrange weights.
    #[test]
    fn tetra_divergence_matches_the_flux() {
        let (x, y, z) = ([0.0, 1.0, 0.2, 0.1], [0.0, 0.1, 1.2, 0.3], [0.0, 0.0, 0.1, 0.9]);
        let g = [7, 3, 9, 1];
        let a = [x[1] - x[0], y[1] - y[0], z[1] - z[0]];
        let b = [x[2] - x[0], y[2] - y[0], z[2] - z[0]];
        let c = [x[3] - x[0], y[3] - y[0], z[3] - z[0]];
        let vol = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0;
        let r = tetra_face_orientation(&g);
        let rt = face_element_info_tetra(&x, &y, &z, &g, 0.1, 0.2, 0.3, false, true);
        let bdm = face_element_info_tetra(&x, &y, &z, &g, 0.1, 0.2, 0.3, true, true);
        for q in 0..4 {
            let want = if r[q] { -1.0 } else { 1.0 };
            assert!((rt.div[q] * vol.abs() - want).abs() < 1e-13, "{q}: {}", rt.div[q] * vol.abs());
            for j in 0..3 {
                assert!((bdm.div[3 * q + j] * vol.abs() - want * 3f64.sqrt() / 3.0).abs() < 1e-13);
            }
        }
        let mut seen = tetra_face_basis_ordering(&g).map(|m| {
            let mut m = m;
            m.sort();
            m
        });
        seen.sort();
        assert!(seen.iter().all(|m| *m == [1, 2, 3]));
    }

    /// Bricks: the divergence of each face function integrates to its flux
    /// (a quarter of the face's reference area 4, with the orientation sign)
    /// and that of each interior function to zero.
    #[test]
    fn brick_divergence_matches_the_flux() {
        let c = |a: f64, b: f64, d: f64| [a + 0.1 * b, b + 0.05 * d, d + 0.02 * a];
        let corners = [[-1.0, -1.0, -1.0], [1.0, -1.0, -1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 1.0], [-1.0, 1.0, 1.0]];
        let (mut x, mut y, mut z) = ([0.0; 8], [0.0; 8], [0.0; 8]);
        for i in 0..8 {
            let p = c(corners[i][0], corners[i][1], corners[i][2]);
            (x[i], y[i], z[i]) = (p[0], p[1], p[2]);
        }
        let g = [11, 4, 9, 2, 7, 15, 3, 8];
        let r = brick_face_orientation(&g);
        let pts = [-0.861_136_311_594_052_6, -0.339_981_043_584_856_3, 0.339_981_043_584_856_3, 0.861_136_311_594_052_6];
        let wts = [0.347_854_845_137_453_9, 0.652_145_154_862_546_1, 0.652_145_154_862_546_1, 0.347_854_845_137_453_9];
        let mut total = vec![0.0f64; 48];
        for (a, wa) in pts.iter().zip(wts) {
            for (b, wb) in pts.iter().zip(wts) {
                for (d, wd) in pts.iter().zip(wts) {
                    let fi = face_element_info_brick(&x, &y, &z, &g, *a, *b, *d, true);
                    for k in 0..48 {
                        total[k] += fi.div[k] * fi.det_f * wa * wb * wd;
                    }
                }
            }
        }
        for q in 0..6 {
            let want = if r[q] { -1.0 } else { 1.0 };
            for j in 0..4 {
                assert!((total[4 * q + j] - want).abs() < 1e-12, "face {q} dof {j}: {}", total[4 * q + j]);
            }
        }
        for k in 24..48 {
            assert!(total[k].abs() < 1e-12, "interior {k}: {}", total[k]);
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
