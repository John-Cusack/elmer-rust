//! H(curl)-conforming edge element bases: Whitney/Nédélec forms and the
//! Piola-transformed first/second-degree Nédélec bases of Elmer.
//!
//! Transliterated from `fem/src/ElemInfo.F90`:
//!
//! * [`edge_element_info`] is `EdgeElementInfo` (Piola version, degrees 1 and
//!   2, Nédélec first kind) for the element types 303, 404, 504, 605, 706 and
//!   808;
//! * [`get_edge_basis`] is `GetEdgeBasis` (the older lowest-order design that
//!   works with the *global* nodal basis derivatives), together with the part
//!   of `ElementInfo` it relies on (`NodalBasisFunctions`,
//!   `NodalFirstDerivatives`, `ElementMetric`) for those linear element types;
//! * [`edge_basis`] is the reference-element view of both.
//!
//! Orientation of the edge functions follows the global node numbers handed
//! in as `gindex` (an edge points towards the larger global number), exactly
//! as `Element % NodeIndexes` does in Elmer.
//!
//! The family bodies (`*_family`) were produced mechanically from the Fortran
//! statements (array indices shifted to 0-based, operation order preserved);
//! they are deliberately left in the shape of the original.
#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::many_single_char_names,
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::collapsible_else_if,
    clippy::manual_memcpy,
    clippy::neg_multiply,
    clippy::identity_op,
    clippy::erasing_op,
    clippy::type_complexity,
    clippy::if_same_then_else,
    clippy::precedence,
    clippy::assign_op_pattern,
    clippy::manual_swap,
    clippy::let_and_return,
    clippy::unnecessary_cast,
    clippy::int_plus_one,
    clippy::overly_complex_bool_expr,
    unused_parens,
    unused_assignments,
    unused_variables,
    unused_mut,
    unused_comparisons,
    non_snake_case
)]

use elmer_core::mesh::ElementType;

/// Edge basis storage as in Elmer's `EdgeBasis(:,3)` / `CurlBasis(54,3)`.
type Mat = [[f64; 3]; 64];

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

#[inline]
fn cube(x: f64) -> f64 {
    x * x * x
}

// ---------------------------------------------------------------------------
// p-reference element nodal bases (PElementBase.F90)
// ---------------------------------------------------------------------------

fn line_nodal_p_basis(node: usize, u: f64) -> f64 {
    match node {
        1 => (1.0 - u) / 2.0,
        2 => (1.0 + u) / 2.0,
        _ => panic!("line_nodal_p_basis: unknown node {node}"),
    }
}

fn d_line_nodal_p_basis(node: usize, _u: f64) -> f64 {
    match node {
        1 => -1.0 / 2.0,
        2 => 1.0 / 2.0,
        _ => panic!("d_line_nodal_p_basis: unknown node {node}"),
    }
}

/// `TriangleNodalPBasis`.
fn triangle_nodal_p_basis(node: usize, u: f64, v: f64) -> f64 {
    match node {
        1 => (1.0 - u - v / f64::sqrt(3.0)) / 2.0,
        2 => (1.0 + u - v / f64::sqrt(3.0)) / 2.0,
        3 => v / f64::sqrt(3.0),
        _ => panic!("triangle_nodal_p_basis: unknown node {node}"),
    }
}

/// `dTriangleNodalPBasis` (third component zero).
fn d_triangle_nodal_p_basis(node: usize, _u: f64, _v: f64) -> [f64; 3] {
    match node {
        1 => [-1.0 / 2.0, -f64::sqrt(3.0) / 6.0, 0.0],
        2 => [1.0 / 2.0, -f64::sqrt(3.0) / 6.0, 0.0],
        3 => [0.0, f64::sqrt(3.0) / 3.0, 0.0],
        _ => panic!("d_triangle_nodal_p_basis: unknown node {node}"),
    }
}

/// `QuadNodalPBasis`.
fn quad_nodal_p_basis(node: usize, u: f64, v: f64) -> f64 {
    match node {
        1 => (1.0 - u) * (1.0 - v) / 4.0,
        2 => (1.0 + u) * (1.0 - v) / 4.0,
        3 => (1.0 + u) * (1.0 + v) / 4.0,
        4 => (1.0 - u) * (1.0 + v) / 4.0,
        _ => panic!("quad_nodal_p_basis: unknown node {node}"),
    }
}

/// `dQuadNodalPBasis`.
fn d_quad_nodal_p_basis(node: usize, u: f64, v: f64) -> [f64; 3] {
    match node {
        1 => [-(1.0 - v) / 4.0, -(1.0 - u) / 4.0, 0.0],
        2 => [(1.0 - v) / 4.0, -(1.0 + u) / 4.0, 0.0],
        3 => [(1.0 + v) / 4.0, (1.0 + u) / 4.0, 0.0],
        4 => [-(1.0 + v) / 4.0, (1.0 - u) / 4.0, 0.0],
        _ => panic!("d_quad_nodal_p_basis: unknown node {node}"),
    }
}

/// `TetraNodalPBasis`.
fn tetra_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    match node {
        1 => 1.0 / 2.0 * (1.0 - u - v / f64::sqrt(3.0) - w / f64::sqrt(6.0)),
        2 => 1.0 / 2.0 * (1.0 + u - v / f64::sqrt(3.0) - w / f64::sqrt(6.0)),
        3 => f64::sqrt(3.0) / 3.0 * (v - w / f64::sqrt(8.0)),
        4 => f64::sqrt(3.0 / 8.0) * w,
        _ => panic!("tetra_nodal_p_basis: unknown node {node}"),
    }
}

/// `dTetraNodalPBasis`.
fn d_tetra_nodal_p_basis(node: usize, _u: f64, _v: f64, _w: f64) -> [f64; 3] {
    match node {
        1 => [-1.0 / 2.0, -f64::sqrt(3.0) / 6.0, -f64::sqrt(6.0) / 12.0],
        2 => [1.0 / 2.0, -f64::sqrt(3.0) / 6.0, -f64::sqrt(6.0) / 12.0],
        3 => [0.0, f64::sqrt(3.0) / 3.0, -f64::sqrt(6.0) / 12.0],
        4 => [0.0, 0.0, f64::sqrt(6.0) / 4.0],
        _ => panic!("d_tetra_nodal_p_basis: unknown node {node}"),
    }
}

/// `PyramidNodalPBasis`.
fn pyramid_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let sq2 = f64::sqrt(2.0);
    let s = w / sq2;
    match node {
        1 => (1.0 - u - v - s + u * v / (1.0 - s)) / 4.0,
        2 => (1.0 + u - v - s - u * v / (1.0 - s)) / 4.0,
        3 => (1.0 + u + v - s + u * v / (1.0 - s)) / 4.0,
        4 => (1.0 - u + v - s - u * v / (1.0 - s)) / 4.0,
        5 => s,
        _ => panic!("pyramid_nodal_p_basis: unknown node {node}"),
    }
}

/// `dPyramidNodalPBasis`.
fn d_pyramid_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> [f64; 3] {
    let sq2 = f64::sqrt(2.0);
    let s = w / sq2;
    let mut grad = [0.0f64; 3];
    match node {
        1 => {
            grad[0] = -1.0 + v / (1.0 - s);
            grad[1] = -1.0 + u / (1.0 - s);
            grad[2] = -1.0 + u * v / sq(1.0 - s);
        }
        2 => {
            grad[0] = 1.0 - v / (1.0 - s);
            grad[1] = -1.0 - u / (1.0 - s);
            grad[2] = -1.0 - u * v / sq(1.0 - s);
        }
        3 => {
            grad[0] = 1.0 + v / (1.0 - s);
            grad[1] = 1.0 + u / (1.0 - s);
            grad[2] = -1.0 + u * v / sq(1.0 - s);
        }
        4 => {
            grad[0] = -1.0 - v / (1.0 - s);
            grad[1] = 1.0 - u / (1.0 - s);
            grad[2] = -1.0 - u * v / sq(1.0 - s);
        }
        5 => grad[2] = 4.0,
        _ => panic!("d_pyramid_nodal_p_basis: unknown node {node}"),
    }
    for g in grad.iter_mut() {
        *g /= 4.0;
    }
    grad[2] /= sq2;
    grad
}

/// `WedgeL`.
fn wedge_l(which: usize, u: f64, v: f64) -> f64 {
    match which {
        1 | 4 => 1.0 / 2.0 * (1.0 - u - v / f64::sqrt(3.0)),
        2 | 5 => 1.0 / 2.0 * (1.0 + u - v / f64::sqrt(3.0)),
        3 | 6 => f64::sqrt(3.0) / 3.0 * v,
        _ => 0.0,
    }
}

/// `dWedgeL`.
fn d_wedge_l(which: usize, _u: f64, _v: f64) -> [f64; 3] {
    match which {
        1 | 4 => [-1.0 / 2.0, -f64::sqrt(3.0) / 6.0, 0.0],
        2 | 5 => [1.0 / 2.0, -f64::sqrt(3.0) / 6.0, 0.0],
        3 | 6 => [0.0, f64::sqrt(3.0) / 3.0, 0.0],
        _ => [0.0; 3],
    }
}

/// `WedgeNodalPBasis`.
fn wedge_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let value = match node {
        1..=3 => wedge_l(node, u, v) * (1.0 - w),
        4..=6 => wedge_l(node, u, v) * (1.0 + w),
        _ => panic!("wedge_nodal_p_basis: unknown node {node}"),
    };
    value / 2.0
}

/// `dWedgeNodalPBasis`.
fn d_wedge_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> [f64; 3] {
    let sign_w = if node <= 3 { -1.0 } else { 1.0 };
    let dl = d_wedge_l(node, u, v);
    let l = wedge_l(node, u, v);
    [
        1.0 / 2.0 * dl[0] * (1.0 + sign_w * w),
        1.0 / 2.0 * dl[1] * (1.0 + sign_w * w),
        sign_w * 1.0 / 2.0 * l,
    ]
}

/// `BrickNodalPBasis`.
fn brick_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> f64 {
    let value = match node {
        1 => (1.0 - u) * (1.0 - v) * (1.0 - w),
        2 => (1.0 + u) * (1.0 - v) * (1.0 - w),
        3 => (1.0 + u) * (1.0 + v) * (1.0 - w),
        4 => (1.0 - u) * (1.0 + v) * (1.0 - w),
        5 => (1.0 - u) * (1.0 - v) * (1.0 + w),
        6 => (1.0 + u) * (1.0 - v) * (1.0 + w),
        7 => (1.0 + u) * (1.0 + v) * (1.0 + w),
        8 => (1.0 - u) * (1.0 + v) * (1.0 + w),
        _ => panic!("brick_nodal_p_basis: unknown node {node}"),
    };
    value / 8.0
}

/// `dBrickNodalPBasis`.
fn d_brick_nodal_p_basis(node: usize, u: f64, v: f64, w: f64) -> [f64; 3] {
    let g = match node {
        1 => [-(1.0 - v) * (1.0 - w), -(1.0 - u) * (1.0 - w), -(1.0 - u) * (1.0 - v)],
        2 => [(1.0 - v) * (1.0 - w), -(1.0 + u) * (1.0 - w), -(1.0 + u) * (1.0 - v)],
        3 => [(1.0 + v) * (1.0 - w), (1.0 + u) * (1.0 - w), -(1.0 + u) * (1.0 + v)],
        4 => [-(1.0 + v) * (1.0 - w), (1.0 - u) * (1.0 - w), -(1.0 - u) * (1.0 + v)],
        5 => [-(1.0 - v) * (1.0 + w), -(1.0 - u) * (1.0 + w), (1.0 - u) * (1.0 - v)],
        6 => [(1.0 - v) * (1.0 + w), -(1.0 + u) * (1.0 + w), (1.0 + u) * (1.0 - v)],
        7 => [(1.0 + v) * (1.0 + w), (1.0 + u) * (1.0 + w), (1.0 + u) * (1.0 + v)],
        8 => [-(1.0 + v) * (1.0 + w), (1.0 - u) * (1.0 + w), (1.0 - u) * (1.0 + v)],
        _ => panic!("d_brick_nodal_p_basis: unknown node {node}"),
    };
    [g[0] / 8.0, g[1] / 8.0, g[2] / 8.0]
}

// ---------------------------------------------------------------------------
// Face orientation helpers (ElemInfo.F90)
// ---------------------------------------------------------------------------

/// `SquareFaceDofsOrdering`: returns `(I1, I2, D1, D2)`.
fn square_face_dofs_ordering(ind: &[i64]) -> (usize, usize, f64, f64) {
    let k = if ind[0] < ind[1] { 1 } else { 2 };
    let l = if ind[3] < ind[2] { 4 } else { 3 };
    let a = if ind[k - 1] > ind[l - 1] { l } else { k };
    match a {
        1 => {
            if ind[1] < ind[3] {
                (1, 2, 1.0, 1.0)
            } else {
                (2, 1, 1.0, 1.0)
            }
        }
        2 => {
            if ind[2] < ind[0] {
                (2, 1, 1.0, -1.0)
            } else {
                (1, 2, -1.0, 1.0)
            }
        }
        3 => {
            if ind[3] < ind[1] {
                (1, 2, -1.0, -1.0)
            } else {
                (2, 1, -1.0, -1.0)
            }
        }
        _ => {
            if ind[0] < ind[2] {
                (2, 1, -1.0, 1.0)
            } else {
                (1, 2, 1.0, -1.0)
            }
        }
    }
}

/// `TriangleFaceDofsOrdering`: returns `(I1, I2, D1, D2, A, B, C)`.
fn triangle_face_dofs_ordering(ind: &[i64]) -> (usize, usize, f64, f64, usize, usize, usize) {
    let mut d1 = 1.0;
    let mut d2 = 1.0;
    let mut i = if ind[0] < ind[1] { 1 } else { 2 };
    if ind[i - 1] > ind[2] {
        i = 3;
    }
    let (i1, i2, j, k);
    match i {
        1 => {
            if ind[2] > ind[1] {
                j = 2;
                k = 3;
                i1 = 1;
                i2 = 2;
            } else {
                j = 3;
                k = 2;
                i1 = 2;
                i2 = 1;
            }
        }
        2 => {
            if ind[2] > ind[0] {
                j = 1;
                k = 3;
                i1 = 1;
                i2 = 3;
                d1 = -1.0;
            } else {
                j = 3;
                k = 1;
                i1 = 3;
                i2 = 1;
                d2 = -1.0;
            }
        }
        _ => {
            if ind[1] > ind[0] {
                j = 1;
                k = 2;
                i1 = 2;
                i2 = 3;
            } else {
                j = 2;
                k = 1;
                i1 = 3;
                i2 = 2;
            }
            d1 = -1.0;
            d2 = -1.0;
        }
    }
    (i1, i2, d1, d2, i, j, k)
}

/// `TriangleFaceDofsOrdering2nd`: returns `(I1, I2, I3)`.
#[allow(dead_code)]
fn triangle_face_dofs_ordering_2nd(ind: &[i64]) -> (usize, usize, usize) {
    let mut i1 = if ind[0] < ind[1] { 1 } else { 2 };
    if ind[i1 - 1] > ind[2] {
        i1 = 3;
    }
    let (i2, i3);
    match i1 {
        1 => {
            if ind[2] > ind[1] {
                i2 = 2;
                i3 = 3;
            } else {
                i2 = 3;
                i3 = 2;
            }
        }
        2 => {
            if ind[2] > ind[0] {
                i2 = 1;
                i3 = 3;
            } else {
                i2 = 3;
                i3 = 1;
            }
        }
        _ => {
            if ind[1] > ind[0] {
                i2 = 1;
                i3 = 2;
            } else {
                i2 = 2;
                i3 = 1;
            }
        }
    }
    (i1, i2, i3)
}


// ---------------------------------------------------------------------------
// Family bodies (generated from EdgeElementInfo)
// ---------------------------------------------------------------------------

fn tri_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 3] = [[1, 2], [2, 3], [3, 1]];

    i = edgemap[0][0];
    j = edgemap[0][1];
    edgebasis[0][0] = (3.0e0 - f64::sqrt(3.0e0) * v) / 6.0e0;
    edgebasis[0][1] = u / (2.0e0 * f64::sqrt(3.0e0));
    curlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0);
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[0][s0] = -edgebasis[0][s0];
        }
        curlbasis[0][2] = -curlbasis[0][2];
    }
    if secondorder {
        edgebasis[1][0] = -(u * (-3.0e0 + f64::sqrt(3.0e0) * v) / 2.0e0);
        edgebasis[1][1] = f64::sqrt(3.0e0) * sq(u) / 2.0e0;
        curlbasis[1][2] = 3.0e0 * f64::sqrt(3.0e0) * u / 2.0e0;
    }
    if secondorder {
        k = 3;
        edgebasis[3][0] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * v / 4.0e0;
        edgebasis[3][1] = f64::sqrt(3.0e0) * (1.0e0 + u) * (-1.0e0 - u + f64::sqrt(3.0e0) * v) / 4.0e0;
        curlbasis[3][2] = -(3.0e0 * (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v)) / 4.0e0;
    } else {
        k = 2;
    }
    i = edgemap[1][0];
    j = edgemap[1][1];
    edgebasis[k - 1][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)));
    edgebasis[k - 1][1] = (1.0 + u) / (2.0e0 * f64::sqrt(3.0e0));
    curlbasis[k - 1][2] = 1.0e0 / f64::sqrt(3.0e0);
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        curlbasis[k - 1][2] = -curlbasis[k - 1][2];
    }
    if secondorder {
        k = 5;
        edgebasis[5][0] = v * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) / 4.0e0;
        edgebasis[5][1] = -(f64::sqrt(3.0e0) * (-1.0e0 + u) * (-1.0e0 + u + f64::sqrt(3.0e0) * v) / 4.0e0);
        curlbasis[5][2] = -(3.0e0 * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v)) / 4.0e0;
    } else {
        k = 3;
    }
    i = edgemap[2][0];
    j = edgemap[2][1];
    edgebasis[k - 1][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)));
    edgebasis[k - 1][1] = ((-1 as f64) + u) / (2.0e0 * f64::sqrt(3.0e0));
    curlbasis[k - 1][2] = 1.0e0 / f64::sqrt(3.0e0);
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        curlbasis[k - 1][2] = -curlbasis[k - 1][2];
    }
    if secondorder {
        trianglefacemap = [1,2,3];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        workbasis[0][0] = (f64::sqrt(3.0e0) - v) * v / 6.0e0;
        workbasis[0][1] = u * v / 6.0e0;
        workcurlbasis[0][2] = (-f64::sqrt(3.0e0) + 3.0e0 * v) / 6.0e0;
        workbasis[1][0] = v * (1.0e0 + u - v / f64::sqrt(3.0e0)) / (4.0e0 * f64::sqrt(3.0e0));
        workbasis[1][1] = (-1.0e0 + u) * (-3.0e0 - 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0));
        workcurlbasis[1][2] = (-f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * v) / 12.0e0;
        workbasis[2][0] = v * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0));
        workbasis[2][1] = -((1.0e0 + u) * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0)));
        workcurlbasis[2][2] = (f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * v) / 12.0e0;
        if redefinefacebasis {
            for s0 in 0..3 {
                edgebasis[6][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
            }
            curlbasis[6][2] = 0.5e0 * d1 * workcurlbasis[i1 - 1][2] + 0.5e0 * d2 * workcurlbasis[i2 - 1][2];
            for s0 in 0..3 {
                edgebasis[7][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
            }
            curlbasis[7][2] = 0.5e0 * d2 * workcurlbasis[i2 - 1][2] - 0.5e0 * d1 * workcurlbasis[i1 - 1][2];
        } else {
            for s0 in 0..3 {
                edgebasis[6][s0] = d1 * workbasis[i1 - 1][s0];
            }
            curlbasis[6][2] = d1 * workcurlbasis[i1 - 1][2];
            for s0 in 0..3 {
                edgebasis[7][s0] = d2 * workbasis[i2 - 1][s0];
            }
            curlbasis[7][2] = d2 * workcurlbasis[i2 - 1][2];
        }
        if scalefacebasis {
            for s0 in 0..3 {
                edgebasis[6][s0] = f64::sqrt(fs1) * edgebasis[6][s0];
            }
            for s0 in 0..3 {
                edgebasis[7][s0] = f64::sqrt(fs2) * edgebasis[7][s0];
            }
            curlbasis[6][2] = f64::sqrt(fs1) * curlbasis[6][2];
            curlbasis[7][2] = f64::sqrt(fs2) * curlbasis[7][2];
        }
    }
    Ok(())
}

fn quad_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 4] = [[1, 2], [2, 3], [3, 4], [4, 1]];

    if secondorder {
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = 0.1e1 / 0.4e1 - v / 0.4e1;
        curlbasis[0][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            curlbasis[0][2] = -curlbasis[0][2];
        }
        edgebasis[1][0] = 0.3e1 * u * (0.1e1 / 0.4e1 - v / 0.4e1);
        curlbasis[1][2] = 0.3e1 / 0.4e1 * u;
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[2][1] = 0.1e1 / 0.4e1 + u / 0.4e1;
        curlbasis[2][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            curlbasis[2][2] = -curlbasis[2][2];
        }
        edgebasis[3][1] = 0.3e1 * v * (0.1e1 / 0.4e1 + u / 0.4e1);
        curlbasis[3][2] = 0.3e1 / 0.4e1 * v;
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[4][0] = -(0.1e1 / 0.4e1) - v / 0.4e1;
        curlbasis[4][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[4][s0] = -edgebasis[4][s0];
            }
            curlbasis[4][2] = -curlbasis[4][2];
        }
        edgebasis[5][0] = -(0.3e1 * u * (-(0.1e1 / 0.4e1) - v / 0.4e1));
        curlbasis[5][2] = -(0.3e1 / 0.4e1 * u);
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[6][1] = -(0.1e1 / 0.4e1) + u / 0.4e1;
        curlbasis[6][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[6][s0] = -edgebasis[6][s0];
            }
            curlbasis[6][2] = -curlbasis[6][2];
        }
        edgebasis[7][1] = -(0.3e1 * v * (-(0.1e1 / 0.4e1) + u / 0.4e1));
        curlbasis[7][2] = -(0.3e1 / 0.4e1 * v);
        squarefacemap = [1,2,3,4];
        for s0 in 0..8 {
            for s1 in 0..3 {
                workbasis[s0][s1] = 0.0e0;
            }
        }
        for s0 in 0..8 {
            for s1 in 0..3 {
                workcurlbasis[s0][s1] = 0.0e0;
            }
        }
        workbasis[0][0] = 0.2e1 * (0.1e1 / 0.2e1 - v / 0.2e1) * (0.1e1 / 0.2e1 + v / 0.2e1);
        workcurlbasis[0][2] = v;
        workbasis[1][0] = 0.12e2 * u * (0.1e1 / 0.2e1 - v / 0.2e1) * (0.1e1 / 0.2e1 + v / 0.2e1);
        workcurlbasis[1][2] = 0.6e1 * u * (0.1e1 / 0.2e1 + v / 0.2e1) - 0.6e1 * u * (0.1e1 / 0.2e1 - v / 0.2e1);
        workbasis[2][1] = 0.2e1 * (0.1e1 / 0.2e1 - u / 0.2e1) * (0.1e1 / 0.2e1 + u / 0.2e1);
        workcurlbasis[2][2] = -u;
        workbasis[3][1] = 0.12e2 * v * (0.1e1 / 0.2e1 - u / 0.2e1) * (0.1e1 / 0.2e1 + u / 0.2e1);
        workcurlbasis[3][2] = -(0.6e1 * v * (0.1e1 / 0.2e1 + u / 0.2e1)) + 0.6e1 * v * (0.1e1 / 0.2e1 - u / 0.2e1);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[8][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[8][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[9][s0] = 0.5e0 * workbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[9][s0] = 0.5e0 * workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[10][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[10][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[11][s0] = 0.5e0 * workbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[11][s0] = 0.5e0 * workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
    } else {
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = (-1.0e0 + v) * v / 4.0e0;
        edgebasis[0][1] = 0.0e0;
        curlbasis[0][2] = (1.0e0 - 2.0 * v) / 4.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            curlbasis[0][2] = -curlbasis[0][2];
        }
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[1][0] = 0.0e0;
        edgebasis[1][1] = u * (1.0e0 + u) / 4.0e0;
        curlbasis[1][2] = (1.0e0 + 2.0 * u) / 4.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[1][s0] = -edgebasis[1][s0];
            }
            curlbasis[1][2] = -curlbasis[1][2];
        }
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[2][0] = -(v * (1.0e0 + v) / 4.0e0);
        edgebasis[2][1] = 0.0e0;
        curlbasis[2][2] = (1.0e0 + 2.0 * v) / 4.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            curlbasis[2][2] = -curlbasis[2][2];
        }
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[3][0] = 0.0e0;
        edgebasis[3][1] = -(((-1 as f64) + u) * u / 4.0e0);
        curlbasis[3][2] = (1.0e0 - 2.0 * u) / 4.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[3][s0] = -edgebasis[3][s0];
            }
            curlbasis[3][2] = -curlbasis[3][2];
        }
        squarefacemap = [1,2,3,4];
        for s0 in 0..3 {
            workbasis[0][s0] = 0.0e0;
        }
        for s0 in 0..3 {
            workbasis[1][s0] = 0.0e0;
        }
        for s0 in 0..3 {
            workcurlbasis[0][s0] = 0.0e0;
        }
        for s0 in 0..3 {
            workcurlbasis[1][s0] = 0.0e0;
        }
        workbasis[0][0] = (1.0e0 - sq(v)) / 2.0e0;
        workbasis[0][1] = 0.0e0;
        workcurlbasis[0][2] = v;
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = (1.0e0 - sq(u)) / 2.0e0;
        workcurlbasis[1][2] = -u;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[4][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[4][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[5][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[5][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
    }
    Ok(())
}

fn tetra_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 6] = [[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]];

    i = edgemap[0][0];
    j = edgemap[0][1];
    edgebasis[0][0] = (6.0e0 - 2.0e0 * f64::sqrt(3.0e0) * v - f64::sqrt(6.0e0) * w) / 24.0e0;
    edgebasis[0][1] = u / (4.0e0 * f64::sqrt(3.0e0));
    edgebasis[0][2] = u / (4.0e0 * f64::sqrt(6.0e0));
    curlbasis[0][0] = 0.0e0;
    curlbasis[0][1] = -(1.0e0 / (2.0e0 * f64::sqrt(6.0e0)));
    curlbasis[0][2] = 1.0e0 / (2.0e0 * f64::sqrt(3.0e0));
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[0][s0] = -edgebasis[0][s0];
        }
        for s0 in 0..3 {
            curlbasis[0][s0] = -curlbasis[0][s0];
        }
    }
    if secondorder {
        edgebasis[1][0] = -(u * (-6.0e0 + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / 4.0e0);
        edgebasis[1][1] = f64::sqrt(3.0e0) * sq(u) / 2.0e0;
        edgebasis[1][2] = f64::sqrt(1.5e0) * sq(u) / 2.0e0;
        curlbasis[1][0] = 0.0e0;
        curlbasis[1][1] = -(3.0e0 * f64::sqrt(1.5e0) * u) / 2.0e0;
        curlbasis[1][2] = 3.0e0 * f64::sqrt(3.0e0) * u / 2.0e0;
    }
    if secondorder {
        k = 3;
        edgebasis[3][0] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * (4.0e0 * v - f64::sqrt(2.0e0) * w) / 16.0e0;
        edgebasis[3][1] = -((1.0e0 + u - f64::sqrt(3.0e0) * v) * (4.0e0 * f64::sqrt(3.0e0) + 4.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * f64::sqrt(2.0e0) * w) / 16.0e0);
        edgebasis[3][2] = -((f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * (-1.0e0 - u + f64::sqrt(3.0e0) * v) / (8.0e0 * f64::sqrt(2.0e0)));
        curlbasis[3][0] = -(9.0e0 * (1.0e0 + u - f64::sqrt(3.0e0) * v)) / (8.0e0 * f64::sqrt(2.0e0));
        curlbasis[3][1] = -(3.0e0 * (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v)) / (8.0e0 * f64::sqrt(2.0e0));
        curlbasis[3][2] = -(3.0e0 * (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v)) / 4.0e0;
    } else {
        k = 2;
    }
    i = edgemap[1][0];
    j = edgemap[1][1];
    edgebasis[k - 1][0] = (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(3.0e0));
    edgebasis[k - 1][1] = (4.0e0 * f64::sqrt(3.0e0) + 4.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * f64::sqrt(2.0e0) * w) / 48.0e0;
    edgebasis[k - 1][2] = -((f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) / (24.0e0 * f64::sqrt(2.0e0)));
    curlbasis[k - 1][0] = 1.0e0 / (4.0e0 * f64::sqrt(2.0e0));
    curlbasis[k - 1][1] = 1.0e0 / (4.0e0 * f64::sqrt(6.0e0));
    curlbasis[k - 1][2] = 1.0e0 / (2.0e0 * f64::sqrt(3.0e0));
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
        }
    }
    if secondorder {
        k = 5;
        edgebasis[5][0] = (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) * (4.0e0 * v - f64::sqrt(2.0e0) * w) / 16.0e0;
        edgebasis[5][1] = -((-1.0e0 + u + f64::sqrt(3.0e0) * v) * (-(4.0e0 * f64::sqrt(3.0e0)) + 4.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * f64::sqrt(2.0e0) * w) / 16.0e0);
        edgebasis[5][2] = (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) * (-1.0e0 + u + f64::sqrt(3.0e0) * v) / (8.0e0 * f64::sqrt(2.0e0));
        curlbasis[5][0] = 9.0e0 * (-1.0e0 + u + f64::sqrt(3.0e0) * v) / (8.0e0 * f64::sqrt(2.0e0));
        curlbasis[5][1] = -(3.0e0 * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v)) / (8.0e0 * f64::sqrt(2.0e0));
        curlbasis[5][2] = -(3.0e0 * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v)) / 4.0e0;
    } else {
        k = 3;
    }
    i = edgemap[2][0];
    j = edgemap[2][1];
    edgebasis[k - 1][0] = (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(3.0e0));
    edgebasis[k - 1][1] = (-(4.0e0 * f64::sqrt(3.0e0)) + 4.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * f64::sqrt(2.0e0) * w) / 48.0e0;
    edgebasis[k - 1][2] = (f64::sqrt(6.0e0) - f64::sqrt(6.0e0) * u - 3.0e0 * f64::sqrt(2.0e0) * v) / 48.0e0;
    curlbasis[k - 1][0] = -(1.0e0 / (4.0e0 * f64::sqrt(2.0e0)));
    curlbasis[k - 1][1] = 1.0e0 / (4.0e0 * f64::sqrt(6.0e0));
    curlbasis[k - 1][2] = 1.0e0 / (2.0e0 * f64::sqrt(3.0e0));
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
        }
    }
    if secondorder {
        k = 7;
        edgebasis[7][0] = 3.0e0 * w * (-f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u + f64::sqrt(2.0e0) * v + 4.0e0 * w) / 16.0e0;
        edgebasis[7][1] = w * (-(3.0e0 * f64::sqrt(2.0e0)) + 3.0e0 * f64::sqrt(2.0e0) * u + f64::sqrt(6.0e0) * v + 4.0e0 * f64::sqrt(3.0e0) * w) / 16.0e0;
        edgebasis[7][2] = -((-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + v) * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v + 2.0e0 * f64::sqrt(6.0e0) * w) / (8.0e0 * f64::sqrt(2.0e0)));
        curlbasis[7][0] = -(3.0e0 * (-(3.0e0 * f64::sqrt(2.0e0)) + 3.0e0 * f64::sqrt(2.0e0) * u + f64::sqrt(6.0e0) * v + 4.0e0 * f64::sqrt(3.0e0) * w)) / 16.0e0;
        curlbasis[7][1] = 9.0e0 * (-f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u + f64::sqrt(2.0e0) * v + 4.0e0 * w) / 16.0e0;
        curlbasis[7][2] = 0.0e0;
    } else {
        k = 4;
    }
    i = edgemap[3][0];
    j = edgemap[3][1];
    edgebasis[k - 1][0] = f64::sqrt(1.5e0) * w / 8.0e0;
    edgebasis[k - 1][1] = w / (8.0e0 * f64::sqrt(2.0e0));
    edgebasis[k - 1][2] = (f64::sqrt(6.0e0) - f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v) / 16.0e0;
    curlbasis[k - 1][0] = -(1.0e0 / (4.0e0 * f64::sqrt(2.0e0)));
    curlbasis[k - 1][1] = f64::sqrt(1.5e0) / 4.0e0;
    curlbasis[k - 1][2] = 0.0e0;
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
        }
    }
    if secondorder {
        k = 9;
        edgebasis[9][0] = 3.0e0 * (f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v - 4.0e0 * w) * w / 16.0e0;
        edgebasis[9][1] = w * (-(3.0e0 * f64::sqrt(2.0e0)) - 3.0e0 * f64::sqrt(2.0e0) * u + f64::sqrt(6.0e0) * v + 4.0e0 * f64::sqrt(3.0e0) * w) / 16.0e0;
        edgebasis[9][2] = (f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v) * (-3.0e0 - 3.0e0 * u + f64::sqrt(3.0e0) * v + 2.0e0 * f64::sqrt(6.0e0) * w) / 16.0e0;
        curlbasis[9][0] = 3.0e0 * (3.0e0 * f64::sqrt(2.0e0) + 3.0e0 * f64::sqrt(2.0e0) * u - f64::sqrt(6.0e0) * v - 4.0e0 * f64::sqrt(3.0e0) * w) / 16.0e0;
        curlbasis[9][1] = 9.0e0 * (f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v - 4.0e0 * w) / 16.0e0;
        curlbasis[9][2] = 0.0e0;
    } else {
        k = 5;
    }
    i = edgemap[4][0];
    j = edgemap[4][1];
    edgebasis[k - 1][0] = -(f64::sqrt(1.5e0) * w / 8.0e0);
    edgebasis[k - 1][1] = w / (8.0e0 * f64::sqrt(2.0e0));
    edgebasis[k - 1][2] = (f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v) / 16.0e0;
    curlbasis[k - 1][0] = -(1.0e0 / (4.0e0 * f64::sqrt(2.0e0)));
    curlbasis[k - 1][1] = -(f64::sqrt(1.5e0) / 4.0e0);
    curlbasis[k - 1][2] = 0.0e0;
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
        }
    }
    if secondorder {
        k = 11;
        edgebasis[11][0] = 0.0e0;
        edgebasis[11][1] = f64::sqrt(3.0e0) * (f64::sqrt(2.0e0) * v - 2.0e0 * w) * w / 4.0e0;
        edgebasis[11][2] = f64::sqrt(1.5e0) * v * (-v + f64::sqrt(2.0e0) * w) / 2.0e0;
        curlbasis[11][0] = -(3.0e0 * (f64::sqrt(6.0e0) * v - 2.0e0 * f64::sqrt(3.0e0) * w)) / 4.0e0;
        curlbasis[11][1] = 0.0e0;
        curlbasis[11][2] = 0.0e0;
    } else {
        k = 6;
    }
    i = edgemap[5][0];
    j = edgemap[5][1];
    edgebasis[k - 1][0] = 0.0e0;
    edgebasis[k - 1][1] = -(w / (4.0e0 * f64::sqrt(2.0e0)));
    edgebasis[k - 1][2] = v / (4.0e0 * f64::sqrt(2.0e0));
    curlbasis[k - 1][0] = 1.0e0 / (2.0e0 * f64::sqrt(2.0e0));
    curlbasis[k - 1][1] = 0.0e0;
    curlbasis[k - 1][2] = 0.0e0;
    if gindexes[j - 1] < gindexes[i - 1] {
        for s0 in 0..3 {
            edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
        }
    }
    if secondorder {
        for k in 1..=6 {
            for s0 in 0..3 {
                edgebasis[2 * (k - 1) + 1 - 1][s0] = 2.0e0 * edgebasis[2 * (k - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[2 * (k - 1) + 1 - 1][s0] = 2.0e0 * curlbasis[2 * (k - 1) + 1 - 1][s0];
            }
        }
    } else {
        for k in 1..=6 {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = 2.0e0 * edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = 2.0e0 * curlbasis[k - 1][s0];
            }
        }
    }
    if secondorder {
        for k in 1..=4 {
            match k {
                1 => {
                    trianglefacemap = [2,1,3];
                }
                2 => {
                    trianglefacemap = [1,2,4];
                }
                3 => {
                    trianglefacemap = [2,3,4];
                }
                4 => {
                    trianglefacemap = [3,1,4];
                }
                _ => {}
            }
            for s0 in 0..3 {
                faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
            }
            (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
            weighted_whitney_forms(&mut workbasis, &mut workcurlbasis, k, u, v, w);
            if redefinefacebasis {
                for s0 in 0..3 {
                    edgebasis[12 + 2 * (k - 1) + 1 - 1][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[12 + 2 * (k - 1) + 1 - 1][s0] = 0.5e0 * d1 * workcurlbasis[i1 - 1][s0] + 0.5e0 * d2 * workcurlbasis[i2 - 1][s0];
                }
                for s0 in 0..3 {
                    edgebasis[12 + 2 * k - 1][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[12 + 2 * k - 1][s0] = 0.5e0 * d2 * workcurlbasis[i2 - 1][s0] - 0.5e0 * d1 * workcurlbasis[i1 - 1][s0];
                }
            } else {
                for s0 in 0..3 {
                    edgebasis[12 + 2 * (k - 1) + 1 - 1][s0] = d1 * workbasis[i1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[12 + 2 * (k - 1) + 1 - 1][s0] = d1 * workcurlbasis[i1 - 1][s0];
                }
                for s0 in 0..3 {
                    edgebasis[12 + 2 * k - 1][s0] = d2 * workbasis[i2 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[12 + 2 * k - 1][s0] = d2 * workcurlbasis[i2 - 1][s0];
                }
            }
        }
        if scalefacebasis {
            for s0 in 0..4 {
                for s1 in 0..3 {
                    edgebasis[13 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * edgebasis[13 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    curlbasis[13 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * curlbasis[13 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    edgebasis[14 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * edgebasis[14 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    curlbasis[14 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * curlbasis[14 - 1 + 2 * s0][s1];
                }
            }
        }
    }
    Ok(())
}

fn pyramid_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 8] = [[1, 2], [2, 3], [4, 3], [1, 4], [1, 5], [2, 5], [3, 5], [4, 5]];

    if secondorder {
        for s0 in 0..16 {
            edgesign[s0] = 1.0e0;
        }
        lbasis[0] = 0.1e1 / 0.4e1 - u / 0.4e1 - v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1 + u * v / ((0.1e1 - w * f64::sqrt(0.2e1) / 0.2e1) * 0.4e1);
        lbasis[1] = 0.1e1 / 0.4e1 + u / 0.4e1 - v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1 - u * v / ((0.1e1 - w * f64::sqrt(0.2e1) / 0.2e1) * 0.4e1);
        lbasis[2] = 0.1e1 / 0.4e1 + u / 0.4e1 + v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1 + u * v / ((0.1e1 - w * f64::sqrt(0.2e1) / 0.2e1) * 0.4e1);
        lbasis[3] = 0.1e1 / 0.4e1 - u / 0.4e1 + v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1 - u * v / ((0.1e1 - w * f64::sqrt(0.2e1) / 0.2e1) * 0.4e1);
        lbasis[4] = w * f64::sqrt(0.2e1) / 0.2e1;
        beta[0] = 0.1e1 / 0.2e1 - u / 0.2e1 - w * f64::sqrt(0.2e1) / 0.4e1;
        beta[1] = 0.1e1 / 0.2e1 - v / 0.2e1 - w * f64::sqrt(0.2e1) / 0.4e1;
        beta[2] = 0.1e1 / 0.2e1 + u / 0.2e1 - w * f64::sqrt(0.2e1) / 0.4e1;
        beta[3] = 0.1e1 / 0.2e1 + v / 0.2e1 - w * f64::sqrt(0.2e1) / 0.4e1;
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = 0.1e1 / 0.4e1 - v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1;
        edgebasis[0][1] = 0.0e0;
        edgebasis[0][2] = f64::sqrt(0.2e1) * u * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[0][0] = f64::sqrt(0.2e1) * u / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.4e1);
        curlbasis[0][1] = -(f64::sqrt(0.2e1) / 0.8e1) - f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[0][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            for s0 in 0..3 {
                curlbasis[0][s0] = -curlbasis[0][s0];
            }
            edgesign[0] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[1][s0] = 3.0e0 * u * edgebasis[0][s0];
        }
        curlbasis[1][0] = 0.3e1 / 0.4e1 * sq(u) * f64::sqrt(0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[1][1] = -(0.3e1 / 0.8e1 * u * f64::sqrt(0.2e1) * (0.3e1 * w * f64::sqrt(0.2e1) + 4.0e0 * v - 0.6e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[1][2] = 0.3e1 / 0.4e1 * u;
        k = 3;
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[k - 1][0] = 0.0e0;
        edgebasis[k - 1][1] = 0.1e1 / 0.4e1 + u / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1;
        edgebasis[k - 1][2] = f64::sqrt(0.2e1) * v * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[k - 1][0] = f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1) + f64::sqrt(0.2e1) / 0.8e1;
        curlbasis[k - 1][1] = f64::sqrt(0.2e1) * v / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.4e1);
        curlbasis[k - 1][2] = 0.1e1 / 0.4e1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * v * edgebasis[k - 1][s0];
        }
        curlbasis[k + 1 - 1][0] = 0.3e1 / 0.8e1 * v * f64::sqrt(0.2e1) * (0.3e1 * w * f64::sqrt(0.2e1) - 4.0e0 * u - 0.6e1) / (w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][1] = 0.3e1 / 0.4e1 * sq(v) * f64::sqrt(0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][2] = 0.3e1 / 0.4e1 * v;
        k = 5;
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[k - 1][0] = 0.1e1 / 0.4e1 + v / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1;
        edgebasis[k - 1][1] = 0.0e0;
        edgebasis[k - 1][2] = f64::sqrt(0.2e1) * u * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[k - 1][0] = -(f64::sqrt(0.2e1) * u / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.4e1));
        curlbasis[k - 1][1] = -(f64::sqrt(0.2e1) / 0.8e1) - f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[k - 1][2] = -(0.1e1 / 0.4e1);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * u * edgebasis[k - 1][s0];
        }
        curlbasis[k + 1 - 1][0] = -(0.3e1 / 0.4e1 * sq(u) * f64::sqrt(0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][1] = -(0.3e1 / 0.8e1 * u * f64::sqrt(0.2e1) * (0.3e1 * w * f64::sqrt(0.2e1) - 4.0e0 * v - 0.6e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][2] = -(0.3e1 / 0.4e1 * u);
        k = 7;
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[k - 1][0] = 0.0e0;
        edgebasis[k - 1][1] = 0.1e1 / 0.4e1 - u / 0.4e1 - w * f64::sqrt(0.2e1) / 0.8e1;
        edgebasis[k - 1][2] = f64::sqrt(0.2e1) * v * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        curlbasis[k - 1][0] = f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1) + f64::sqrt(0.2e1) / 0.8e1;
        curlbasis[k - 1][1] = -(f64::sqrt(0.2e1) * v / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.4e1));
        curlbasis[k - 1][2] = -(0.1e1 / 0.4e1);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * v * edgebasis[k - 1][s0];
        }
        curlbasis[k + 1 - 1][0] = 0.3e1 / 0.8e1 * v * f64::sqrt(0.2e1) * (0.3e1 * w * f64::sqrt(0.2e1) + 4.0e0 * u - 0.6e1) / (w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][1] = -(0.3e1 / 0.4e1 * sq(v) * f64::sqrt(0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][2] = -(0.3e1 / 0.4e1 * v);
        k = 9;
        i = edgemap[4][0];
        j = edgemap[4][1];
        edgebasis[k - 1][0] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        edgebasis[k - 1][1] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        edgebasis[k - 1][2] = -(f64::sqrt(0.2e1) / 0.4e1 * (0.2e1 * f64::sqrt(0.2e1) * u * v * w - 0.2e1 * f64::sqrt(0.2e1) * u * w - 0.2e1 * f64::sqrt(0.2e1) * v * w + u * sq(w) + v * sq(w) + 0.2e1 * w * f64::sqrt(0.2e1) - 0.2e1 * u * v - sq(w) + 0.2e1 * u + 0.2e1 * v - 0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k - 1][0] = (-(f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.2e1 * u * w - 0.2e1 * f64::sqrt(0.2e1) + 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1);
        curlbasis[k - 1][1] = -((-(f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.2e1 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1));
        curlbasis[k - 1][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * edgesign[k - 1] * edgebasis[k - 1][s0] * (lbasis[4] - lbasis[0] + lbasis[2]);
        }
        curlbasis[k + 1 - 1][0] = 0.3e1 / 0.8e1 * (-(0.9e1 * f64::sqrt(0.2e1) * u * sq(w)) - 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) + 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w - 0.8e1 * u * v * w - 0.6e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) + 0.24e2 * u * w + 0.12e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][1] = -(0.3e1 / 0.8e1 * (-(0.3e1 * f64::sqrt(0.2e1) * u * sq(w)) - 0.9e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.4e1 * sq(v) * w - 0.6e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) + 0.12e2 * u * w + 0.24e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][2] = 0.3e1 / 0.8e1 * w * f64::sqrt(0.2e1) * (u - v) / (w * f64::sqrt(0.2e1) - 0.2e1);
        k = 11;
        i = edgemap[5][0];
        j = edgemap[5][1];
        edgebasis[k - 1][0] = -(w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1));
        edgebasis[k - 1][1] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        edgebasis[k - 1][2] = f64::sqrt(0.2e1) / 0.4e1 * (0.2e1 * f64::sqrt(0.2e1) * u * v * w - 0.2e1 * f64::sqrt(0.2e1) * u * w + 0.2e1 * f64::sqrt(0.2e1) * v * w + u * sq(w) - v * sq(w) - 0.2e1 * w * f64::sqrt(0.2e1) - 0.2e1 * u * v + sq(w) + 0.2e1 * u - 0.2e1 * v + 0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k - 1][0] = -((f64::sqrt(0.2e1) * sq(w) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.2e1 * u * w + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1));
        curlbasis[k - 1][1] = (-(f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.2e1 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1);
        curlbasis[k - 1][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * edgesign[k - 1] * edgebasis[k - 1][s0] * (lbasis[4] - lbasis[1] + lbasis[3]);
        }
        curlbasis[k + 1 - 1][0] = 0.3e1 / 0.8e1 * (0.9e1 * f64::sqrt(0.2e1) * u * sq(w) - 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.8e1 * u * v * w - 0.6e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) - 0.24e2 * u * w + 0.12e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][1] = -(0.3e1 / 0.8e1 * (-(0.3e1 * f64::sqrt(0.2e1) * u * sq(w)) + 0.9e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) - 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.4e1 * sq(v) * w + 0.6e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) + 0.12e2 * u * w - 0.24e2 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][2] = 0.3e1 / 0.8e1 * w * f64::sqrt(0.2e1) * (u + v) / (w * f64::sqrt(0.2e1) - 0.2e1);
        k = 13;
        i = edgemap[6][0];
        j = edgemap[6][1];
        edgebasis[k - 1][0] = -(w * f64::sqrt(0.2e1) / 0.8e1 * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        edgebasis[k - 1][1] = -(w * f64::sqrt(0.2e1) / 0.8e1 * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        edgebasis[k - 1][2] = -(f64::sqrt(0.2e1) / 0.4e1 * (0.2e1 * f64::sqrt(0.2e1) * u * v * w + 0.2e1 * f64::sqrt(0.2e1) * u * w + 0.2e1 * f64::sqrt(0.2e1) * v * w - u * sq(w) - v * sq(w) + 0.2e1 * w * f64::sqrt(0.2e1) - 0.2e1 * u * v - sq(w) - 0.2e1 * u - 0.2e1 * v - 0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k - 1][0] = (f64::sqrt(0.2e1) * sq(w) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.2e1 * u * w + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1);
        curlbasis[k - 1][1] = -((f64::sqrt(0.2e1) * sq(w) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.2e1 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1));
        curlbasis[k - 1][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * edgesign[k - 1] * edgebasis[k - 1][s0] * (lbasis[4] - lbasis[2] + lbasis[0]);
        }
        curlbasis[k + 1 - 1][0] = -(0.3e1 / 0.8e1 * (0.9e1 * f64::sqrt(0.2e1) * u * sq(w) + 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) + 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w - 0.8e1 * u * v * w - 0.6e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) - 0.24e2 * u * w - 0.12e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][1] = 0.3e1 / 0.8e1 * (0.3e1 * f64::sqrt(0.2e1) * u * sq(w) + 0.9e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.4e1 * sq(v) * w - 0.6e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) - 0.12e2 * u * w - 0.24e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][2] = -(0.3e1 / 0.8e1 * w * f64::sqrt(0.2e1) * (u - v) / (w * f64::sqrt(0.2e1) - 0.2e1));
        k = 15;
        i = edgemap[7][0];
        j = edgemap[7][1];
        edgebasis[k - 1][0] = w * f64::sqrt(0.2e1) / 0.8e1 * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1);
        edgebasis[k - 1][1] = -(w * f64::sqrt(0.2e1) / 0.8e1 * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / (w * f64::sqrt(0.2e1) - 0.2e1));
        edgebasis[k - 1][2] = f64::sqrt(0.2e1) / 0.4e1 * (0.2e1 * f64::sqrt(0.2e1) * u * v * w + 0.2e1 * f64::sqrt(0.2e1) * u * w - 0.2e1 * f64::sqrt(0.2e1) * v * w - u * sq(w) + v * sq(w) - 0.2e1 * w * f64::sqrt(0.2e1) - 0.2e1 * u * v + sq(w) - 0.2e1 * u + 0.2e1 * v + 0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k - 1][0] = -((-(f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.2e1 * u * w - 0.2e1 * f64::sqrt(0.2e1) + 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1));
        curlbasis[k - 1][1] = (f64::sqrt(0.2e1) * sq(w) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.2e1 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / (sq(w * f64::sqrt(0.2e1) - 0.2e1) * 0.2e1);
        curlbasis[k - 1][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[k - 1][s0] = -edgebasis[k - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k - 1][s0] = -curlbasis[k - 1][s0];
            }
            edgesign[k - 1] = -1.0e0;
        }
        for s0 in 0..3 {
            edgebasis[k + 1 - 1][s0] = 3.0e0 * edgesign[k - 1] * edgebasis[k - 1][s0] * (lbasis[4] - lbasis[3] + lbasis[1]);
        }
        curlbasis[k + 1 - 1][0] = -(0.3e1 / 0.8e1 * (-(0.9e1 * f64::sqrt(0.2e1) * u * sq(w)) + 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.8e1 * u * v * w - 0.6e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) + 0.24e2 * u * w - 0.12e2 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[k + 1 - 1][1] = 0.3e1 / 0.8e1 * (0.3e1 * f64::sqrt(0.2e1) * u * sq(w) - 0.9e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) - 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.13e2 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.4e1 * sq(v) * w + 0.6e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) - 0.12e2 * u * w + 0.24e2 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.16e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        curlbasis[k + 1 - 1][2] = -(0.3e1 / 0.8e1 * w * f64::sqrt(0.2e1) * (u + v) / (w * f64::sqrt(0.2e1) - 0.2e1));
        squarefacemap = [1,2,3,4];
        for s0 in 0..3 {
            workbasis[0][s0] = 2.0e0 * (edgesign[0] * edgebasis[0][s0] * beta[3] + edgesign[4] * edgebasis[4][s0] * beta[1]) / (1.0e0 - lbasis[4]);
        }
        workcurlbasis[0][0] = -(0.2e1 * u * v * f64::sqrt(0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[0][1] = -((f64::sqrt(0.2e1) * sq(w) + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[0][2] = -(0.2e1 * v / (w * f64::sqrt(0.2e1) - 0.2e1));
        for s0 in 0..3 {
            workbasis[1][s0] = 3.0e0 * workbasis[0][s0] * u;
        }
        workcurlbasis[1][0] = -(0.6e1 * sq(u) * f64::sqrt(0.2e1) * v / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[1][1] = 0.3e1 / 0.2e1 * u * (0.2e1 * f64::sqrt(0.2e1) * sq(v) - 0.3e1 * f64::sqrt(0.2e1) * sq(w) - 0.6e1 * f64::sqrt(0.2e1) + 0.12e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        workcurlbasis[1][2] = -(0.6e1 * u * v / (w * f64::sqrt(0.2e1) - 0.2e1));
        for s0 in 0..3 {
            workbasis[2][s0] = 2.0e0 * (edgesign[2] * edgebasis[2][s0] * beta[0] + edgesign[6] * edgebasis[6][s0] * beta[2]) / (1.0e0 - lbasis[4]);
        }
        workcurlbasis[2][0] = (f64::sqrt(0.2e1) * sq(w) + 0.2e1 * f64::sqrt(0.2e1) - 0.4e1 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        workcurlbasis[2][1] = 0.2e1 * u * v * f64::sqrt(0.2e1) / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        workcurlbasis[2][2] = 0.2e1 * u / (w * f64::sqrt(0.2e1) - 0.2e1);
        for s0 in 0..3 {
            workbasis[3][s0] = 3.0e0 * workbasis[2][s0] * v;
        }
        workcurlbasis[3][0] = -(0.3e1 / 0.2e1 * v * (0.2e1 * f64::sqrt(0.2e1) * sq(u) - 0.3e1 * f64::sqrt(0.2e1) * sq(w) - 0.6e1 * f64::sqrt(0.2e1) + 0.12e2 * w) / sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[3][1] = 0.6e1 * f64::sqrt(0.2e1) * sq(v) * u / sq(w * f64::sqrt(0.2e1) - 0.2e1);
        workcurlbasis[3][2] = 0.6e1 * u * v / (w * f64::sqrt(0.2e1) - 0.2e1);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[16][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[16][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[17][s0] = workbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[17][s0] = workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[18][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[18][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[19][s0] = workbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[19][s0] = workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        trianglefacemap = [1,2,5];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            workbasis[0][s0] = lbasis[4] * edgesign[0] * edgebasis[0][s0];
        }
        workcurlbasis[0][0] = w * u / (w * f64::sqrt(0.2e1) - 0.2e1) / 0.4e1;
        workcurlbasis[0][1] = (-(0.3e1 * f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.4e1 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.8e1 * w) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        workcurlbasis[0][2] = w * f64::sqrt(0.2e1) / 0.8e1;
        for s0 in 0..3 {
            workbasis[1][s0] = beta[2] * edgesign[8] * edgebasis[8][s0];
        }
        workcurlbasis[1][0] = (f64::sqrt(0.2e1) * u * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.3e1 * cube(w) - 0.2e1 * u * w - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[1][1] = -((-(0.3e1 * f64::sqrt(0.2e1) * u * sq(w)) + 0.2e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) - 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.3e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) + 0.2e1 * v * f64::sqrt(0.2e1) + 0.12e2 * u * w - 0.6e1 * v * w - 0.2e1 * f64::sqrt(0.2e1) + 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[1][2] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.16e2);
        for s0 in 0..3 {
            workbasis[2][s0] = beta[0] * edgesign[10] * edgebasis[10][s0];
        }
        workcurlbasis[2][0] = (-(f64::sqrt(0.2e1) * u * sq(w)) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.3e1 * cube(w) + 0.2e1 * u * w - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[2][1] = -((-(0.3e1 * f64::sqrt(0.2e1) * u * sq(w)) - 0.2e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * u * v * f64::sqrt(0.2e1) + 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.3e1 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) - 0.2e1 * v * f64::sqrt(0.2e1) + 0.12e2 * u * w + 0.6e1 * v * w + 0.2e1 * f64::sqrt(0.2e1) - 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[2][2] = -(w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1)));
        if redefinefacebasis {
            for s0 in 0..3 {
                edgebasis[20][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[20][s0] = 0.5e0 * d1 * workcurlbasis[i1 - 1][s0] + 0.5e0 * d2 * workcurlbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[21][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[21][s0] = 0.5e0 * d2 * workcurlbasis[i2 - 1][s0] - 0.5e0 * d1 * workcurlbasis[i1 - 1][s0];
            }
        } else {
            for s0 in 0..3 {
                edgebasis[20][s0] = d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[20][s0] = d1 * workcurlbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[21][s0] = d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[21][s0] = d2 * workcurlbasis[i2 - 1][s0];
            }
        }
        trianglefacemap = [2,3,5];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            workbasis[0][s0] = lbasis[4] * edgesign[2] * edgebasis[2][s0];
        }
        workcurlbasis[0][0] = (0.3e1 * f64::sqrt(0.2e1) * sq(w) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.4e1 * u * w + 0.2e1 * f64::sqrt(0.2e1) - 0.8e1 * w) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.8e1);
        workcurlbasis[0][1] = w * v / (w * f64::sqrt(0.2e1) - 0.2e1) / 0.4e1;
        workcurlbasis[0][2] = w * f64::sqrt(0.2e1) / 0.8e1;
        for s0 in 0..3 {
            workbasis[1][s0] = beta[3] * edgesign[10] * edgebasis[10][s0];
        }
        workcurlbasis[1][0] = -((0.2e1 * f64::sqrt(0.2e1) * u * sq(w) + 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v + 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.3e1 * cube(w) + 0.2e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) - 0.6e1 * u * w - 0.12e2 * w * v + 0.2e1 * f64::sqrt(0.2e1) - 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[1][1] = (f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(v) * w + 0.3e1 * cube(w) - 0.2e1 * w * v - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[1][2] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1));
        for s0 in 0..3 {
            workbasis[2][s0] = beta[1] * edgesign[12] * edgebasis[12][s0];
        }
        workcurlbasis[2][0] = -((-(0.2e1 * f64::sqrt(0.2e1) * u * sq(w)) + 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v - 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.3e1 * cube(w) - 0.2e1 * u * f64::sqrt(0.2e1) + 0.6e1 * v * f64::sqrt(0.2e1) + 0.6e1 * u * w - 0.12e2 * w * v - 0.2e1 * f64::sqrt(0.2e1) + 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[2][1] = (-(f64::sqrt(0.2e1) * v * sq(w)) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(v) * w + 0.3e1 * cube(w) + 0.2e1 * w * v - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[2][2] = -(w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / ((w * f64::sqrt(0.2e1) - 0.2e1) * 0.16e2));
        if redefinefacebasis {
            for s0 in 0..3 {
                edgebasis[22][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[22][s0] = 0.5e0 * d1 * workcurlbasis[i1 - 1][s0] + 0.5e0 * d2 * workcurlbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[23][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[23][s0] = 0.5e0 * d2 * workcurlbasis[i2 - 1][s0] - 0.5e0 * d1 * workcurlbasis[i1 - 1][s0];
            }
        } else {
            for s0 in 0..3 {
                edgebasis[22][s0] = d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[22][s0] = d1 * workcurlbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[23][s0] = d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[23][s0] = d2 * workcurlbasis[i2 - 1][s0];
            }
        }
        trianglefacemap = [3,4,5];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            workbasis[0][s0] = -(lbasis[4] * edgesign[4] * edgebasis[4][s0]);
        }
        workcurlbasis[0][0] = w * u / (w * f64::sqrt(0.2e1) - 0.2e1) / 0.4e1;
        workcurlbasis[0][1] = (0.3e1 * f64::sqrt(0.2e1) * sq(w) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.4e1 * w * v + 0.2e1 * f64::sqrt(0.2e1) - 0.8e1 * w) / (0.8e1 * (w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[0][2] = w * f64::sqrt(0.2e1) / 0.8e1;
        for s0 in 0..3 {
            workbasis[1][s0] = beta[0] * edgesign[12] * edgebasis[12][s0];
        }
        workcurlbasis[1][0] = -((-(f64::sqrt(0.2e1) * u * sq(w)) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.3e1 * cube(w) + 0.2e1 * u * w - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[1][1] = (0.3e1 * f64::sqrt(0.2e1) * u * sq(w) - 0.2e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v - 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.3e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) - 0.2e1 * v * f64::sqrt(0.2e1) - 0.12e2 * u * w + 0.6e1 * w * v - 0.2e1 * f64::sqrt(0.2e1) + 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[1][2] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * u - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1));
        for s0 in 0..3 {
            workbasis[2][s0] = beta[2] * edgesign[14] * edgebasis[14][s0];
        }
        workcurlbasis[2][0] = -((f64::sqrt(0.2e1) * u * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(u) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.3e1 * cube(w) - 0.2e1 * u * w - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[2][1] = (0.3e1 * f64::sqrt(0.2e1) * u * sq(w) + 0.2e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v + 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.3e1 * cube(w) + 0.6e1 * u * f64::sqrt(0.2e1) + 0.2e1 * v * f64::sqrt(0.2e1) - 0.12e2 * u * w - 0.6e1 * w * v + 0.2e1 * f64::sqrt(0.2e1) - 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[2][2] = -(w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * u - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1)));
        if redefinefacebasis {
            for s0 in 0..3 {
                edgebasis[24][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[24][s0] = 0.5e0 * d1 * workcurlbasis[i1 - 1][s0] + 0.5e0 * d2 * workcurlbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[25][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[25][s0] = 0.5e0 * d2 * workcurlbasis[i2 - 1][s0] - 0.5e0 * d1 * workcurlbasis[i1 - 1][s0];
            }
        } else {
            for s0 in 0..3 {
                edgebasis[24][s0] = d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[24][s0] = d1 * workcurlbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[25][s0] = d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[25][s0] = d2 * workcurlbasis[i2 - 1][s0];
            }
        }
        trianglefacemap = [4,1,5];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            workbasis[0][s0] = -(lbasis[4] * edgesign[6] * edgebasis[6][s0]);
        }
        workcurlbasis[0][0] = (-(0.3e1 * f64::sqrt(0.2e1) * sq(w)) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.4e1 * u * w - 0.2e1 * f64::sqrt(0.2e1) + 0.8e1 * w) / (0.8e1 * (w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[0][1] = w * v / (w * f64::sqrt(0.2e1) - 0.2e1) / 0.4e1;
        workcurlbasis[0][2] = w * f64::sqrt(0.2e1) / 0.8e1;
        for s0 in 0..3 {
            workbasis[1][s0] = beta[1] * edgesign[14] * edgebasis[14][s0];
        }
        workcurlbasis[1][0] = (-(0.2e1 * f64::sqrt(0.2e1) * u * sq(w)) - 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v + 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w - 0.3e1 * cube(w) - 0.2e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) + 0.6e1 * u * w + 0.12e2 * w * v + 0.2e1 * f64::sqrt(0.2e1) - 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[1][1] = -((-(f64::sqrt(0.2e1) * v * sq(w)) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(v) * w + 0.3e1 * cube(w) + 0.2e1 * w * v - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[1][2] = w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) - 2.0e0 * v - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1));
        for s0 in 0..3 {
            workbasis[2][s0] = beta[3] * edgesign[8] * edgebasis[8][s0];
        }
        workcurlbasis[2][0] = (0.2e1 * f64::sqrt(0.2e1) * u * sq(w) - 0.3e1 * f64::sqrt(0.2e1) * v * sq(w) + 0.6e1 * f64::sqrt(0.2e1) * u * v - 0.7e1 * f64::sqrt(0.2e1) * sq(w) - 0.8e1 * u * v * w + 0.3e1 * cube(w) + 0.2e1 * u * f64::sqrt(0.2e1) - 0.6e1 * v * f64::sqrt(0.2e1) - 0.6e1 * u * w + 0.12e2 * w * v - 0.2e1 * f64::sqrt(0.2e1) + 0.10e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        workcurlbasis[2][1] = -((f64::sqrt(0.2e1) * v * sq(w) + 0.4e1 * f64::sqrt(0.2e1) * sq(v) - 0.8e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(v) * w + 0.3e1 * cube(w) - 0.2e1 * w * v - 0.4e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.8e1 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        workcurlbasis[2][2] = -(w * f64::sqrt(0.2e1) * (w * f64::sqrt(0.2e1) + 2.0e0 * v - 0.2e1) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1)));
        if redefinefacebasis {
            for s0 in 0..3 {
                edgebasis[26][s0] = 0.5e0 * d1 * workbasis[i1 - 1][s0] + 0.5e0 * d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[26][s0] = 0.5e0 * d1 * workcurlbasis[i1 - 1][s0] + 0.5e0 * d2 * workcurlbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[27][s0] = 0.5e0 * d2 * workbasis[i2 - 1][s0] - 0.5e0 * d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[27][s0] = 0.5e0 * d2 * workcurlbasis[i2 - 1][s0] - 0.5e0 * d1 * workcurlbasis[i1 - 1][s0];
            }
        } else {
            for s0 in 0..3 {
                edgebasis[26][s0] = d1 * workbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[26][s0] = d1 * workcurlbasis[i1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[27][s0] = d2 * workbasis[i2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[27][s0] = d2 * workcurlbasis[i2 - 1][s0];
            }
        }
        for s0 in 0..3 {
            edgebasis[28][s0] = lbasis[4] * beta[3] * edgesign[0] * edgebasis[0][s0];
        }
        curlbasis[28][0] = u * v * w / (0.4e1 * (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[28][1] = (0.2e1 * f64::sqrt(0.2e1) * sq(v) - 0.9e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(v) * w + 0.4e1 * cube(w) - 0.2e1 * f64::sqrt(0.2e1) + 0.12e2 * w) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[28][2] = f64::sqrt(0.2e1) * v * w / 0.8e1;
        for s0 in 0..3 {
            edgebasis[29][s0] = lbasis[4] * beta[2] * edgesign[6] * edgebasis[6][s0];
        }
        curlbasis[29][0] = -((0.2e1 * f64::sqrt(0.2e1) * sq(u) - 0.9e1 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * sq(u) * w + 0.4e1 * cube(w) - 0.2e1 * f64::sqrt(0.2e1) + 0.12e2 * w) / (0.16e2 * (w * f64::sqrt(0.2e1) - 0.2e1)));
        curlbasis[29][1] = -(u * v * w / (0.4e1 * (w * f64::sqrt(0.2e1) - 0.2e1)));
        curlbasis[29][2] = -(f64::sqrt(0.2e1) * u * w / 0.8e1);
        for s0 in 0..3 {
            edgebasis[30][s0] = beta[2] * beta[3] * edgesign[8] * edgebasis[8][s0];
        }
        curlbasis[30][0] = (0.2e1 * f64::sqrt(0.2e1) * sq(u) * sq(w) + 0.2e1 * f64::sqrt(0.2e1) * u * v * sq(w) - 0.2e1 * f64::sqrt(0.2e1) * sq(sq(w)) + 0.6e1 * f64::sqrt(0.2e1) * sq(u) * v - 0.11e2 * f64::sqrt(0.2e1) * v * sq(w) - 0.8e1 * sq(u) * v * w + 0.4e1 * v * cube(w) + 0.2e1 * f64::sqrt(0.2e1) * sq(u) - 0.15e2 * f64::sqrt(0.2e1) * sq(w) - 0.6e1 * sq(u) * w - 0.4e1 * u * v * w + 0.13e2 * cube(w) - 0.6e1 * v * f64::sqrt(0.2e1) + 0.20e2 * w * v - 0.2e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.16e2 * sq(w * f64::sqrt(0.2e1) - 0.2e1));
        curlbasis[30][1] = -((0.2e1 * f64::sqrt(0.2e1) * u * v * sq(w) + 0.2e1 * f64::sqrt(0.2e1) * sq(v) * sq(w) - 0.2e1 * f64::sqrt(0.2e1) * sq(sq(w)) + 0.6e1 * f64::sqrt(0.2e1) * u * sq(v) - 0.11e2 * f64::sqrt(0.2e1) * u * sq(w) - 0.8e1 * u * sq(v) * w + 0.4e1 * u * cube(w) + 0.2e1 * f64::sqrt(0.2e1) * sq(v) - 0.15e2 * f64::sqrt(0.2e1) * sq(w) - 0.4e1 * u * v * w - 0.6e1 * sq(v) * w + 0.13e2 * cube(w) - 0.6e1 * u * f64::sqrt(0.2e1) + 0.20e2 * u * w - 0.2e1 * f64::sqrt(0.2e1) + 0.14e2 * w) / (0.16e2 * sq(w * f64::sqrt(0.2e1) - 0.2e1)));
        curlbasis[30][2] = -((u - v) * w * f64::sqrt(0.2e1) / 0.16e2);
        if scalefacebasis {
            for s0 in 0..4 {
                for s1 in 0..3 {
                    edgebasis[21 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * edgebasis[21 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    curlbasis[21 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * curlbasis[21 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    edgebasis[22 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * edgebasis[22 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..4 {
                for s1 in 0..3 {
                    curlbasis[22 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * curlbasis[22 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    edgebasis[28 + s0][s1] = f64::sqrt(506.9e0) * edgebasis[28 + s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    curlbasis[28 + s0][s1] = f64::sqrt(506.9e0) * curlbasis[28 + s0][s1];
                }
            }
            for s0 in 0..3 {
                edgebasis[30][s0] = f64::sqrt(167.8e0) * edgebasis[30][s0];
            }
            for s0 in 0..3 {
                curlbasis[30][s0] = f64::sqrt(167.8e0) * curlbasis[30][s0];
            }
        }
    } else {
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = v * ((-1 as f64) + 2.0 * v / (2.0 - f64::sqrt(2.0e0) * w)) / 4.0e0;
        edgebasis[0][1] = 0.0e0;
        edgebasis[0][2] = u * v * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * v + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[0][0] = u * (-f64::sqrt(2.0e0) + 2.0 * f64::sqrt(2.0e0) * v + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[0][1] = v * (f64::sqrt(2.0e0) - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[0][2] = ((-2 as f64) + 4.0 * v + f64::sqrt(2.0e0) * w) / ((-8 as f64) + 4.0 * f64::sqrt(2.0e0) * w);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            for s0 in 0..3 {
                curlbasis[0][s0] = -curlbasis[0][s0];
            }
        }
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[1][0] = 0.0e0;
        edgebasis[1][1] = u * (1.0 + 2.0 * u / (2.0 - f64::sqrt(2.0e0) * w)) / 4.0e0;
        edgebasis[1][2] = u * v * (f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * u - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[1][0] = u * (f64::sqrt(2.0e0) - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[1][1] = -(v * (f64::sqrt(2.0e0) + 2.0 * f64::sqrt(2.0e0) * u - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w)));
        curlbasis[1][2] = (2.0 + 4.0 * u - f64::sqrt(2.0e0) * w) / (8.0 - 4.0 * f64::sqrt(2.0e0) * w);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[1][s0] = -edgebasis[1][s0];
            }
            for s0 in 0..3 {
                curlbasis[1][s0] = -curlbasis[1][s0];
            }
        }
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[2][0] = v * (1.0 + 2.0 * v / (2.0 - f64::sqrt(2.0e0) * w)) / 4.0e0;
        edgebasis[2][1] = 0.0e0;
        edgebasis[2][2] = u * v * (f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * v - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[2][0] = u * (f64::sqrt(2.0e0) + 2.0 * f64::sqrt(2.0e0) * v - w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[2][1] = v * (-f64::sqrt(2.0e0) + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[2][2] = (2.0 + 4.0 * v - f64::sqrt(2.0e0) * w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            for s0 in 0..3 {
                curlbasis[2][s0] = -curlbasis[2][s0];
            }
        }
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[3][0] = 0.0e0;
        edgebasis[3][1] = u * ((-1 as f64) + 2.0 * u / (2.0 - f64::sqrt(2.0e0) * w)) / 4.0e0;
        edgebasis[3][2] = u * v * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * u + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[3][0] = u * (-f64::sqrt(2.0e0) + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[3][1] = -(v * (-f64::sqrt(2.0e0) + 2.0 * f64::sqrt(2.0e0) * u + w) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w)));
        curlbasis[3][2] = (2.0 - 4.0 * u - f64::sqrt(2.0e0) * w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[3][s0] = -edgebasis[3][s0];
            }
            for s0 in 0..3 {
                curlbasis[3][s0] = -curlbasis[3][s0];
            }
        }
        i = edgemap[4][0];
        j = edgemap[4][1];
        edgebasis[4][0] = w * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * v + w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[4][1] = w * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * u + w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[4][2] = (u * (-(2.0 * f64::sqrt(2.0e0)) + 2.0 * v * (f64::sqrt(2.0e0) - 2.0 * w) + 4.0 * w - f64::sqrt(2.0e0) * sq(w)) - ((-1 as f64) + v) * (2.0 * f64::sqrt(2.0e0) - 4.0 * w + f64::sqrt(2.0e0) * sq(w))) / (4.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[4][0] = (-(2.0 * f64::sqrt(2.0e0)) + 2.0 * u * (f64::sqrt(2.0e0) - w) + 4.0 * w - f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[4][1] = (2.0 * f64::sqrt(2.0e0) - 2.0 * f64::sqrt(2.0e0) * v - 4.0 * w + 2.0 * v * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[4][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[4][s0] = -edgebasis[4][s0];
            }
            for s0 in 0..3 {
                curlbasis[4][s0] = -curlbasis[4][s0];
            }
        }
        i = edgemap[5][0];
        j = edgemap[5][1];
        edgebasis[5][0] = w * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * v + w) / (8.0e0 - 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[5][1] = w * (-f64::sqrt(2.0e0) - f64::sqrt(2.0e0) * u + w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[5][2] = (-(((-1 as f64) + v) * (2.0 * f64::sqrt(2.0e0) - 4.0 * w + f64::sqrt(2.0e0) * sq(w))) + u * (2.0 * f64::sqrt(2.0e0) - 2.0 * f64::sqrt(2.0e0) * v - 4.0 * w + 4.0 * v * w + f64::sqrt(2.0e0) * sq(w))) / (4.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[5][0] = -((2.0 * f64::sqrt(2.0e0) + 2.0 * u * (f64::sqrt(2.0e0) - w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w)));
        curlbasis[5][1] = (-(2.0 * f64::sqrt(2.0e0)) + 2.0 * v * (f64::sqrt(2.0e0) - w) + 4.0 * w - f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[5][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[5][s0] = -edgebasis[5][s0];
            }
            for s0 in 0..3 {
                curlbasis[5][s0] = -curlbasis[5][s0];
            }
        }
        i = edgemap[6][0];
        j = edgemap[6][1];
        edgebasis[6][0] = (f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * v - w) * w / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[6][1] = (f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * u - w) * w / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[6][2] = ((1.0 + v) * (2.0 * f64::sqrt(2.0e0) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) + u * (2.0 * f64::sqrt(2.0e0) + 2.0 * v * (f64::sqrt(2.0e0) - 2.0 * w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w))) / (4.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[6][0] = (2.0 * f64::sqrt(2.0e0) + 2.0 * u * (f64::sqrt(2.0e0) - w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w));
        curlbasis[6][1] = -((2.0 * f64::sqrt(2.0e0) + 2.0 * v * (f64::sqrt(2.0e0) - w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq((-2 as f64) + f64::sqrt(2.0e0) * w)));
        curlbasis[6][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[6][s0] = -edgebasis[6][s0];
            }
            for s0 in 0..3 {
                curlbasis[6][s0] = -curlbasis[6][s0];
            }
        }
        i = edgemap[7][0];
        j = edgemap[7][1];
        edgebasis[7][0] = w * (-f64::sqrt(2.0e0) - f64::sqrt(2.0e0) * v + w) / (-8.0e0 + 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[7][1] = w * (-f64::sqrt(2.0e0) + f64::sqrt(2.0e0) * u + w) / (8.0e0 - 4.0 * f64::sqrt(2.0e0) * w);
        edgebasis[7][2] = ((1.0 + v) * (2.0 * f64::sqrt(2.0e0) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) - u * (2.0 * f64::sqrt(2.0e0) + 2.0 * v * (f64::sqrt(2.0e0) - 2.0 * w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w))) / (4.0e0 * sq(-2.0e0 + f64::sqrt(2.0e0) * w));
        curlbasis[7][0] = (2.0 * f64::sqrt(2.0e0) - 2.0 * f64::sqrt(2.0e0) * u - 4.0 * w + 2.0 * u * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq(-2.0e0 + f64::sqrt(2.0e0) * w));
        curlbasis[7][1] = (2.0 * f64::sqrt(2.0e0) + 2.0 * v * (f64::sqrt(2.0e0) - w) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) / (2.0e0 * sq(-2.0e0 + f64::sqrt(2.0e0) * w));
        curlbasis[7][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[7][s0] = -edgebasis[7][s0];
            }
            for s0 in 0..3 {
                curlbasis[7][s0] = -curlbasis[7][s0];
            }
        }
        squarefacemap = [1,2,3,4];
        workbasis[0][0] = (2.0e0 - 2.0 * sq(v) - 2.0 * f64::sqrt(2.0e0) * w + sq(w)) / (4.0e0 - 2.0 * f64::sqrt(2.0e0) * w);
        workbasis[0][1] = 0.0e0;
        workbasis[0][2] = u * (1.0e0 - 4.0 * sq(v) / sq(-2.0e0 + f64::sqrt(2.0e0) * w)) / (2.0e0 * f64::sqrt(2.0e0));
        workcurlbasis[0][0] = -(2.0 * f64::sqrt(2.0e0) * u * v) / sq(-2.0e0 + f64::sqrt(2.0e0) * w);
        workcurlbasis[0][1] = (-(2.0 * f64::sqrt(2.0e0)) + 4.0 * w - f64::sqrt(2.0e0) * sq(w)) / sq(-2.0e0 + f64::sqrt(2.0e0) * w);
        workcurlbasis[0][2] = 2.0e0 * v / (2.0e0 - f64::sqrt(2.0e0) * w);
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = (2.0e0 - 2.0 * sq(u) - 2.0 * f64::sqrt(2.0e0) * w + sq(w)) / (4.0e0 - 2.0 * f64::sqrt(2.0e0) * w);
        workbasis[1][2] = v * (1.0e0 - 4.0 * sq(u) / sq(-2.0e0 + f64::sqrt(2.0e0) * w)) / (2.0e0 * f64::sqrt(2.0e0));
        workcurlbasis[1][0] = (2.0 * f64::sqrt(2.0e0) - 4.0 * w + f64::sqrt(2.0e0) * sq(w)) / sq(-2.0e0 + f64::sqrt(2.0e0) * w);
        workcurlbasis[1][1] = 2.0 * f64::sqrt(2.0e0) * u * v / sq(-2.0e0 + f64::sqrt(2.0e0) * w);
        workcurlbasis[1][2] = 2.0 * u / (-2.0e0 + f64::sqrt(2.0e0) * w);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[8][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[8][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[9][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[9][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
    }
    Ok(())
}

fn prism_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 9] = [[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4], [1, 4], [2, 5], [3, 6]];

    if secondorder {
        h1 = 0.5e0 * (1.0 - w);
        dh1 = -0.5e0;
        h2 = 0.5e0 * (1.0 + w);
        dh2 = 0.5e0;
        h3 = h1 * h2;
        dh3 = -(0.5e0 * w);
        workbasis[0][0] = (3.0e0 - f64::sqrt(3.0e0) * v) / 6.0e0;
        workbasis[0][1] = u / (2.0e0 * f64::sqrt(3.0e0));
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0);
        workbasis[1][0] = -(u * (-3.0e0 + f64::sqrt(3.0e0) * v) / 2.0e0);
        workbasis[1][1] = f64::sqrt(3.0e0) * sq(u) / 2.0e0;
        workcurlbasis[1][2] = 3.0e0 * f64::sqrt(3.0e0) * u / 2.0e0;
        i = edgemap[0][0];
        j = edgemap[0][1];
        for s0 in 0..2 {
            edgebasis[0][s0] = workbasis[0][s0] * h1;
        }
        curlbasis[0][0] = -(workbasis[0][1] * dh1);
        curlbasis[0][1] = workbasis[0][0] * dh1;
        curlbasis[0][2] = workcurlbasis[0][2] * h1;
        for s0 in 0..2 {
            edgebasis[1][s0] = workbasis[1][s0] * h1;
        }
        curlbasis[1][0] = -(workbasis[1][1] * dh1);
        curlbasis[1][1] = workbasis[1][0] * dh1;
        curlbasis[1][2] = workcurlbasis[1][2] * h1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            for s0 in 0..3 {
                curlbasis[0][s0] = -curlbasis[0][s0];
            }
        }
        i = edgemap[3][0];
        j = edgemap[3][1];
        for s0 in 0..2 {
            edgebasis[6][s0] = workbasis[0][s0] * h2;
        }
        curlbasis[6][0] = -(workbasis[0][1] * dh2);
        curlbasis[6][1] = workbasis[0][0] * dh2;
        curlbasis[6][2] = workcurlbasis[0][2] * h2;
        for s0 in 0..2 {
            edgebasis[7][s0] = workbasis[1][s0] * h2;
        }
        curlbasis[7][0] = -(workbasis[1][1] * dh2);
        curlbasis[7][1] = workbasis[1][0] * dh2;
        curlbasis[7][2] = workcurlbasis[1][2] * h2;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[6][s0] = -edgebasis[6][s0];
            }
            for s0 in 0..3 {
                curlbasis[6][s0] = -curlbasis[6][s0];
            }
        }
        workbasis[0][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)));
        workbasis[0][1] = (1.0 + u) / (2.0e0 * f64::sqrt(3.0e0));
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0);
        workbasis[1][0] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * v / 4.0e0;
        workbasis[1][1] = f64::sqrt(3.0e0) * (1.0e0 + u) * (-1.0e0 - u + f64::sqrt(3.0e0) * v) / 4.0e0;
        workcurlbasis[1][2] = -(3.0e0 * (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v)) / 4.0e0;
        i = edgemap[1][0];
        j = edgemap[1][1];
        for s0 in 0..2 {
            edgebasis[2][s0] = workbasis[0][s0] * h1;
        }
        curlbasis[2][0] = -(workbasis[0][1] * dh1);
        curlbasis[2][1] = workbasis[0][0] * dh1;
        curlbasis[2][2] = workcurlbasis[0][2] * h1;
        for s0 in 0..2 {
            edgebasis[3][s0] = workbasis[1][s0] * h1;
        }
        curlbasis[3][0] = -(workbasis[1][1] * dh1);
        curlbasis[3][1] = workbasis[1][0] * dh1;
        curlbasis[3][2] = workcurlbasis[1][2] * h1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            for s0 in 0..3 {
                curlbasis[2][s0] = -curlbasis[2][s0];
            }
        }
        i = edgemap[4][0];
        j = edgemap[4][1];
        for s0 in 0..2 {
            edgebasis[8][s0] = workbasis[0][s0] * h2;
        }
        curlbasis[8][0] = -(workbasis[0][1] * dh2);
        curlbasis[8][1] = workbasis[0][0] * dh2;
        curlbasis[8][2] = workcurlbasis[0][2] * h2;
        for s0 in 0..2 {
            edgebasis[9][s0] = workbasis[1][s0] * h2;
        }
        curlbasis[9][0] = -(workbasis[1][1] * dh2);
        curlbasis[9][1] = workbasis[1][0] * dh2;
        curlbasis[9][2] = workcurlbasis[1][2] * h2;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[8][s0] = -edgebasis[8][s0];
            }
            for s0 in 0..3 {
                curlbasis[8][s0] = -curlbasis[8][s0];
            }
        }
        workbasis[0][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)));
        workbasis[0][1] = ((-1 as f64) + u) / (2.0e0 * f64::sqrt(3.0e0));
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0);
        workbasis[1][0] = v * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) / 4.0e0;
        workbasis[1][1] = -(f64::sqrt(3.0e0) * (-1.0e0 + u) * (-1.0e0 + u + f64::sqrt(3.0e0) * v) / 4.0e0);
        workcurlbasis[1][2] = -(3.0e0 * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v)) / 4.0e0;
        i = edgemap[2][0];
        j = edgemap[2][1];
        for s0 in 0..2 {
            edgebasis[4][s0] = workbasis[0][s0] * h1;
        }
        curlbasis[4][0] = -(workbasis[0][1] * dh1);
        curlbasis[4][1] = workbasis[0][0] * dh1;
        curlbasis[4][2] = workcurlbasis[0][2] * h1;
        for s0 in 0..2 {
            edgebasis[5][s0] = workbasis[1][s0] * h1;
        }
        curlbasis[5][0] = -(workbasis[1][1] * dh1);
        curlbasis[5][1] = workbasis[1][0] * dh1;
        curlbasis[5][2] = workcurlbasis[1][2] * h1;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[4][s0] = -edgebasis[4][s0];
            }
            for s0 in 0..3 {
                curlbasis[4][s0] = -curlbasis[4][s0];
            }
        }
        i = edgemap[5][0];
        j = edgemap[5][1];
        for s0 in 0..2 {
            edgebasis[10][s0] = workbasis[0][s0] * h2;
        }
        curlbasis[10][0] = -(workbasis[0][1] * dh2);
        curlbasis[10][1] = workbasis[0][0] * dh2;
        curlbasis[10][2] = workcurlbasis[0][2] * h2;
        for s0 in 0..2 {
            edgebasis[11][s0] = workbasis[1][s0] * h2;
        }
        curlbasis[11][0] = -(workbasis[1][1] * dh2);
        curlbasis[11][1] = workbasis[1][0] * dh2;
        curlbasis[11][2] = workcurlbasis[1][2] * h2;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..2 {
                edgebasis[10][s0] = -edgebasis[10][s0];
            }
            for s0 in 0..3 {
                curlbasis[10][s0] = -curlbasis[10][s0];
            }
        }
        for q in 1..=3 {
            i = edgemap[6 + q - 1][0];
            j = edgemap[6 + q - 1][1];
            for s0 in 0..2 {
                grad[s0] = d_triangle_nodal_p_basis(q, u, v)[s0];
            }
            edgebasis[12 + (q - 1) * 2 + 1 - 1][2] = 0.5e0 * triangle_nodal_p_basis(q, u, v);
            curlbasis[12 + (q - 1) * 2 + 1 - 1][0] = 0.5e0 * grad[1];
            curlbasis[12 + (q - 1) * 2 + 1 - 1][1] = -(0.5e0 * grad[0]);
            edgebasis[12 + (q - 1) * 2 + 2 - 1][2] = 3.0e0 * edgebasis[12 + (q - 1) * 2 + 1 - 1][2] * w;
            curlbasis[12 + (q - 1) * 2 + 2 - 1][0] = 1.5e0 * grad[1] * w;
            curlbasis[12 + (q - 1) * 2 + 2 - 1][1] = -(1.5e0 * grad[0] * w);
            if gindexes[j - 1] < gindexes[i - 1] {
                edgebasis[12 + (q - 1) * 2 + 1 - 1][2] = -edgebasis[12 + (q - 1) * 2 + 1 - 1][2];
                for s0 in 0..2 {
                    curlbasis[12 + (q - 1) * 2 + 1 - 1][s0] = -curlbasis[12 + (q - 1) * 2 + 1 - 1][s0];
                }
            }
        }
        trianglefacemap = [1,2,3];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        workbasis[0][0] = (f64::sqrt(3.0e0) - v) * v / 6.0e0;
        workbasis[0][1] = u * v / 6.0e0;
        workcurlbasis[0][2] = (-f64::sqrt(3.0e0) + 3.0e0 * v) / 6.0e0;
        workbasis[1][0] = v * (1.0e0 + u - v / f64::sqrt(3.0e0)) / (4.0e0 * f64::sqrt(3.0e0));
        workbasis[1][1] = (-1.0e0 + u) * (-3.0e0 - 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0));
        workcurlbasis[1][2] = (-f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * v) / 12.0e0;
        workbasis[2][0] = v * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0));
        workbasis[2][1] = -((1.0e0 + u) * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0)));
        workcurlbasis[2][2] = (f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * v) / 12.0e0;
        if redefinefacebasis {
            for s0 in 0..2 {
                edgebasis[18][s0] = (d1 * workbasis[i1 - 1][s0] + d2 * workbasis[i2 - 1][s0]) * 0.5e0 * h1;
            }
            for s0 in 0..2 {
                edgebasis[19][s0] = (-(d1 * workbasis[i1 - 1][s0]) + d2 * workbasis[i2 - 1][s0]) * 0.5e0 * h1;
            }
            curlbasis[18][0] = -((d1 * workbasis[i1 - 1][1] + d2 * workbasis[i2 - 1][1]) * 0.5e0 * dh1);
            curlbasis[18][1] = (d1 * workbasis[i1 - 1][0] + d2 * workbasis[i2 - 1][0]) * 0.5e0 * dh1;
            curlbasis[18][2] = (d1 * workcurlbasis[i1 - 1][2] + d2 * workcurlbasis[i2 - 1][2]) * 0.5e0 * h1;
            curlbasis[19][0] = -((-(d1 * workbasis[i1 - 1][1]) + d2 * workbasis[i2 - 1][1]) * 0.5e0 * dh1);
            curlbasis[19][1] = (-(d1 * workbasis[i1 - 1][0]) + d2 * workbasis[i2 - 1][0]) * 0.5e0 * dh1;
            curlbasis[19][2] = (-(d1 * workcurlbasis[i1 - 1][2]) + d2 * workcurlbasis[i2 - 1][2]) * 0.5e0 * h1;
        } else {
            for s0 in 0..2 {
                edgebasis[18][s0] = d1 * workbasis[i1 - 1][s0] * h1;
            }
            curlbasis[18][0] = -(d1 * workbasis[i1 - 1][1] * dh1);
            curlbasis[18][1] = d1 * workbasis[i1 - 1][0] * dh1;
            curlbasis[18][2] = d1 * workcurlbasis[i1 - 1][2] * h1;
            for s0 in 0..2 {
                edgebasis[19][s0] = d2 * workbasis[i2 - 1][s0] * h1;
            }
            curlbasis[19][0] = -(d2 * workbasis[i2 - 1][1] * dh1);
            curlbasis[19][1] = d2 * workbasis[i2 - 1][0] * dh1;
            curlbasis[19][2] = d2 * workcurlbasis[i2 - 1][2] * h1;
        }
        trianglefacemap = [4,5,6];
        for s0 in 0..3 {
            faceindices[s0] = gindexes[trianglefacemap[s0] - 1];
        }
        (i1, i2, d1, d2, _, _, _) = triangle_face_dofs_ordering(&faceindices);
        if redefinefacebasis {
            for s0 in 0..2 {
                edgebasis[20][s0] = (d1 * workbasis[i1 - 1][s0] + d2 * workbasis[i2 - 1][s0]) * 0.5e0 * h2;
            }
            for s0 in 0..2 {
                edgebasis[21][s0] = (-(d1 * workbasis[i1 - 1][s0]) + d2 * workbasis[i2 - 1][s0]) * 0.5e0 * h2;
            }
            curlbasis[20][0] = -((d1 * workbasis[i1 - 1][1] + d2 * workbasis[i2 - 1][1]) * 0.5e0 * dh2);
            curlbasis[20][1] = (d1 * workbasis[i1 - 1][0] + d2 * workbasis[i2 - 1][0]) * 0.5e0 * dh2;
            curlbasis[20][2] = (d1 * workcurlbasis[i1 - 1][2] + d2 * workcurlbasis[i2 - 1][2]) * 0.5e0 * h2;
            curlbasis[21][0] = -((-(d1 * workbasis[i1 - 1][1]) + d2 * workbasis[i2 - 1][1]) * 0.5e0 * dh2);
            curlbasis[21][1] = (-(d1 * workbasis[i1 - 1][0]) + d2 * workbasis[i2 - 1][0]) * 0.5e0 * dh2;
            curlbasis[21][2] = (-(d1 * workcurlbasis[i1 - 1][2]) + d2 * workcurlbasis[i2 - 1][2]) * 0.5e0 * h2;
        } else {
            for s0 in 0..2 {
                edgebasis[20][s0] = d1 * workbasis[i1 - 1][s0] * h2;
            }
            curlbasis[20][0] = -(d1 * workbasis[i1 - 1][1] * dh2);
            curlbasis[20][1] = d1 * workbasis[i1 - 1][0] * dh2;
            curlbasis[20][2] = d1 * workcurlbasis[i1 - 1][2] * h2;
            for s0 in 0..2 {
                edgebasis[21][s0] = d2 * workbasis[i2 - 1][s0] * h2;
            }
            curlbasis[21][0] = -(d2 * workbasis[i2 - 1][1] * dh2);
            curlbasis[21][1] = d2 * workbasis[i2 - 1][0] * dh2;
            curlbasis[21][2] = d2 * workcurlbasis[i2 - 1][2] * h2;
        }
        if scalefacebasis {
            for s0 in 0..2 {
                for s1 in 0..3 {
                    edgebasis[19 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * edgebasis[19 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    curlbasis[19 - 1 + 2 * s0][s1] = f64::sqrt(fs1) * curlbasis[19 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    edgebasis[20 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * edgebasis[20 - 1 + 2 * s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    curlbasis[20 - 1 + 2 * s0][s1] = f64::sqrt(fs2) * curlbasis[20 - 1 + 2 * s0][s1];
                }
            }
        }
        squarefacemap = [1,2,5,4];
        for s0 in 0..8 {
            for s1 in 0..3 {
                workbasis[s0][s1] = 0.0e0;
            }
        }
        for s0 in 0..8 {
            for s1 in 0..3 {
                workcurlbasis[s0][s1] = 0.0e0;
            }
        }
        workbasis[0][0] = (3.0e0 - f64::sqrt(3.0e0) * v) / 6.0e0 * 4.0e0 * h3;
        workbasis[0][1] = u / (2.0e0 * f64::sqrt(3.0e0)) * 4.0e0 * h3;
        workcurlbasis[0][0] = -(workbasis[0][1] / h3 * dh3);
        workcurlbasis[0][1] = workbasis[0][0] / h3 * dh3;
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0) * 4.0e0 * h3;
        workbasis[1][0] = -(u * (-3.0e0 + f64::sqrt(3.0e0) * v) / 2.0e0 * 4.0e0 * h3);
        workbasis[1][1] = f64::sqrt(3.0e0) * sq(u) / 2.0e0 * 4.0e0 * h3;
        workcurlbasis[1][0] = -(workbasis[1][1] / h3 * dh3);
        workcurlbasis[1][1] = workbasis[1][0] / h3 * dh3;
        workcurlbasis[1][2] = 3.0e0 * f64::sqrt(3.0e0) * u / 2.0e0 * 4.0e0 * h3;
        workbasis[2][2] = 2.0e0 * triangle_nodal_p_basis(1, u, v) * triangle_nodal_p_basis(2, u, v);
        for s0 in 0..2 {
            grad[s0] = d_triangle_nodal_p_basis(1, u, v)[s0] * triangle_nodal_p_basis(2, u, v) + triangle_nodal_p_basis(1, u, v) * d_triangle_nodal_p_basis(2, u, v)[s0];
        }
        workcurlbasis[2][0] = 2.0e0 * grad[1];
        workcurlbasis[2][1] = -(2.0e0 * grad[0]);
        workbasis[3][2] = 3.0e0 * workbasis[2][2] * w;
        workcurlbasis[3][0] = 6.0e0 * grad[1] * w;
        workcurlbasis[3][1] = -(6.0e0 * grad[0] * w);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[22][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[22][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[23][s0] = workbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[23][s0] = workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[24][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[24][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[25][s0] = workbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[25][s0] = workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        squarefacemap = [2,3,6,5];
        for s0 in 0..8 {
            for s1 in 0..3 {
                workbasis[s0][s1] = 0.0e0;
            }
        }
        for s0 in 0..8 {
            for s1 in 0..3 {
                workcurlbasis[s0][s1] = 0.0e0;
            }
        }
        workbasis[0][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)) * 4.0e0 * h3);
        workbasis[0][1] = (1.0 + u) / (2.0e0 * f64::sqrt(3.0e0)) * 4.0e0 * h3;
        workcurlbasis[0][0] = -(workbasis[0][1] / h3 * dh3);
        workcurlbasis[0][1] = workbasis[0][0] / h3 * dh3;
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0) * 4.0e0 * h3;
        workbasis[1][0] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * v / 4.0e0 * 4.0e0 * h3;
        workbasis[1][1] = f64::sqrt(3.0e0) * (1.0e0 + u) * (-1.0e0 - u + f64::sqrt(3.0e0) * v) / 4.0e0 * 4.0e0 * h3;
        workcurlbasis[1][0] = -(workbasis[1][1] / h3 * dh3);
        workcurlbasis[1][1] = workbasis[1][0] / h3 * dh3;
        workcurlbasis[1][2] = -(3.0e0 * (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v)) / 4.0e0 * 4.0e0 * h3;
        workbasis[2][2] = 2.0e0 * triangle_nodal_p_basis(2, u, v) * triangle_nodal_p_basis(3, u, v);
        for s0 in 0..2 {
            grad[s0] = d_triangle_nodal_p_basis(2, u, v)[s0] * triangle_nodal_p_basis(3, u, v) + triangle_nodal_p_basis(2, u, v) * d_triangle_nodal_p_basis(3, u, v)[s0];
        }
        workcurlbasis[2][0] = 2.0e0 * grad[1];
        workcurlbasis[2][1] = -(2.0e0 * grad[0]);
        workbasis[3][2] = 3.0e0 * workbasis[2][2] * w;
        workcurlbasis[3][0] = 6.0e0 * grad[1] * w;
        workcurlbasis[3][1] = -(6.0e0 * grad[0] * w);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[26][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[26][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[27][s0] = workbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[27][s0] = workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[28][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[28][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[29][s0] = workbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[29][s0] = workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        squarefacemap = [3,1,4,6];
        for s0 in 0..8 {
            for s1 in 0..3 {
                workbasis[s0][s1] = 0.0e0;
            }
        }
        for s0 in 0..8 {
            for s1 in 0..3 {
                workcurlbasis[s0][s1] = 0.0e0;
            }
        }
        workbasis[0][0] = -(v / (2.0e0 * f64::sqrt(3.0e0)) * 4.0e0 * h3);
        workbasis[0][1] = ((-1 as f64) + u) / (2.0e0 * f64::sqrt(3.0e0)) * 4.0e0 * h3;
        workcurlbasis[0][0] = -(workbasis[0][1] / h3 * dh3);
        workcurlbasis[0][1] = workbasis[0][0] / h3 * dh3;
        workcurlbasis[0][2] = 1.0e0 / f64::sqrt(3.0e0) * 4.0e0 * h3;
        workbasis[1][0] = v * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) / 4.0e0 * 4.0e0 * h3;
        workbasis[1][1] = -(f64::sqrt(3.0e0) * (-1.0e0 + u) * (-1.0e0 + u + f64::sqrt(3.0e0) * v) / 4.0e0 * 4.0e0 * h3);
        workcurlbasis[1][0] = -(workbasis[1][1] / h3 * dh3);
        workcurlbasis[1][1] = workbasis[1][0] / h3 * dh3;
        workcurlbasis[1][2] = -(3.0e0 * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v)) / 4.0e0 * 4.0e0 * h3;
        workbasis[2][2] = 2.0e0 * triangle_nodal_p_basis(3, u, v) * triangle_nodal_p_basis(1, u, v);
        for s0 in 0..2 {
            grad[s0] = d_triangle_nodal_p_basis(3, u, v)[s0] * triangle_nodal_p_basis(1, u, v) + triangle_nodal_p_basis(3, u, v) * d_triangle_nodal_p_basis(1, u, v)[s0];
        }
        workcurlbasis[2][0] = 2.0e0 * grad[1];
        workcurlbasis[2][1] = -(2.0e0 * grad[0]);
        workbasis[3][2] = 3.0e0 * workbasis[2][2] * w;
        workcurlbasis[3][0] = 6.0e0 * grad[1] * w;
        workcurlbasis[3][1] = -(6.0e0 * grad[0] * w);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[squarefacemap[s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[30][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[30][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[31][s0] = workbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[31][s0] = workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[32][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[32][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[33][s0] = workbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[33][s0] = workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
        }
        edgebasis[34][0] = v * (1.0e0 + u - v / f64::sqrt(3.0e0)) / (4.0e0 * f64::sqrt(3.0e0)) * h3;
        edgebasis[34][1] = (-1.0e0 + u) * (-3.0e0 - 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0)) * h3;
        curlbasis[34][0] = -(edgebasis[34][1] / h3 * dh3);
        curlbasis[34][1] = edgebasis[34][0] / h3 * dh3;
        curlbasis[34][2] = (-f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * v) / 12.0e0 * h3;
        edgebasis[35][0] = v * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0)) * h3;
        edgebasis[35][1] = -((1.0e0 + u) * (-3.0e0 + 3.0e0 * u + f64::sqrt(3.0e0) * v) / (12.0e0 * f64::sqrt(3.0e0)) * h3);
        curlbasis[35][0] = -(edgebasis[35][1] / h3 * dh3);
        curlbasis[35][1] = edgebasis[35][0] / h3 * dh3;
        curlbasis[35][2] = (f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * v) / 12.0e0 * h3;
        if scalefacebasis {
            for s0 in 0..2 {
                for s1 in 0..2 {
                    edgebasis[34 + s0][s1] = f64::sqrt(150.0e0) * edgebasis[34 + s0][s1];
                }
            }
            for s0 in 0..2 {
                for s1 in 0..3 {
                    curlbasis[34 + s0][s1] = f64::sqrt(150.0e0) * curlbasis[34 + s0][s1];
                }
            }
        }
    } else {
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = -((-3.0e0 + f64::sqrt(3.0e0) * v) * (-1.0e0 + w) * w / 12.0e0);
        edgebasis[0][1] = u * (-1.0e0 + w) * w / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[0][2] = 0.0e0;
        curlbasis[0][0] = u * (1.0e0 - 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0));
        curlbasis[0][1] = -((-3.0e0 + f64::sqrt(3.0e0) * v) * (-1.0e0 + 2.0 * w) / 12.0e0);
        curlbasis[0][2] = (-1.0e0 + w) * w / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            for s0 in 0..3 {
                curlbasis[0][s0] = -curlbasis[0][s0];
            }
        }
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[1][0] = -(v * (-1.0e0 + w) * w / (4.0e0 * f64::sqrt(3.0e0)));
        edgebasis[1][1] = (1.0e0 + u) * (-1.0e0 + w) * w / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[1][2] = 0.0e0;
        curlbasis[1][0] = (1.0e0 + u) * (1.0e0 - 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0));
        curlbasis[1][1] = v * (1.0e0 - 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0));
        curlbasis[1][2] = (-1.0e0 + w) * w / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[1][s0] = -edgebasis[1][s0];
            }
            for s0 in 0..3 {
                curlbasis[1][s0] = -curlbasis[1][s0];
            }
        }
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[2][0] = -(v * (-1.0e0 + w) * w / (4.0e0 * f64::sqrt(3.0e0)));
        edgebasis[2][1] = (-1.0e0 + u) * (-1.0e0 + w) * w / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[2][2] = 0.0e0;
        curlbasis[2][0] = (-1.0e0 + u) * (1.0e0 - 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0));
        curlbasis[2][1] = v * (1.0e0 - 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0));
        curlbasis[2][2] = (-1.0e0 + w) * w / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            for s0 in 0..3 {
                curlbasis[2][s0] = -curlbasis[2][s0];
            }
        }
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[3][0] = -((-3.0e0 + f64::sqrt(3.0e0) * v) * w * (1.0e0 + w) / 12.0e0);
        edgebasis[3][1] = u * w * (1.0e0 + w) / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[3][2] = 0.0e0;
        curlbasis[3][0] = -(u * (1.0e0 + 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0)));
        curlbasis[3][1] = -((-3.0e0 + f64::sqrt(3.0e0) * v) * (1.0e0 + 2.0e0 * w) / 12.0e0);
        curlbasis[3][2] = w * (1.0e0 + w) / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[3][s0] = -edgebasis[3][s0];
            }
            for s0 in 0..3 {
                curlbasis[3][s0] = -curlbasis[3][s0];
            }
        }
        i = edgemap[4][0];
        j = edgemap[4][1];
        edgebasis[4][0] = -(v * w * (1.0e0 + w) / (4.0e0 * f64::sqrt(3.0e0)));
        edgebasis[4][1] = (1.0e0 + u) * w * (1.0e0 + w) / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[4][2] = 0.0e0;
        curlbasis[4][0] = -((1.0e0 + u) * (1.0e0 + 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0)));
        curlbasis[4][1] = -(v * (1.0e0 + 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0)));
        curlbasis[4][2] = w * (1.0e0 + w) / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[4][s0] = -edgebasis[4][s0];
            }
            for s0 in 0..3 {
                curlbasis[4][s0] = -curlbasis[4][s0];
            }
        }
        i = edgemap[5][0];
        j = edgemap[5][1];
        edgebasis[5][0] = -(v * w * (1.0e0 + w) / (4.0e0 * f64::sqrt(3.0e0)));
        edgebasis[5][1] = (-1.0e0 + u) * w * (1.0e0 + w) / (4.0e0 * f64::sqrt(3.0e0));
        edgebasis[5][2] = 0.0e0;
        curlbasis[5][0] = -((-1.0e0 + u) * (1.0e0 + 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0)));
        curlbasis[5][1] = -(v * (1.0e0 + 2.0e0 * w) / (4.0e0 * f64::sqrt(3.0e0)));
        curlbasis[5][2] = w * (1.0e0 + w) / (2.0e0 * f64::sqrt(3.0e0));
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[5][s0] = -edgebasis[5][s0];
            }
            for s0 in 0..3 {
                curlbasis[5][s0] = -curlbasis[5][s0];
            }
        }
        i = edgemap[6][0];
        j = edgemap[6][1];
        edgebasis[6][0] = 0.0e0;
        edgebasis[6][1] = 0.0e0;
        edgebasis[6][2] = (3.0 * sq(u) + v * (-f64::sqrt(3.0e0) + v) + u * (-3.0e0 + 2.0 * f64::sqrt(3.0e0) * v)) / 12.0e0;
        curlbasis[6][0] = (-f64::sqrt(3.0e0) + 2.0 * f64::sqrt(3.0e0) * u + 2.0 * v) / 12.0e0;
        curlbasis[6][1] = (3.0e0 - 6.0 * u - 2.0 * f64::sqrt(3.0e0) * v) / 12.0e0;
        curlbasis[6][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[6][s0] = -edgebasis[6][s0];
            }
            for s0 in 0..3 {
                curlbasis[6][s0] = -curlbasis[6][s0];
            }
        }
        i = edgemap[7][0];
        j = edgemap[7][1];
        edgebasis[7][0] = 0.0e0;
        edgebasis[7][1] = 0.0e0;
        edgebasis[7][2] = (3.0 * sq(u) + v * (-f64::sqrt(3.0e0) + v) + u * (3.0e0 - 2.0 * f64::sqrt(3.0e0) * v)) / 12.0e0;
        curlbasis[7][0] = (-f64::sqrt(3.0e0) - 2.0 * f64::sqrt(3.0e0) * u + 2.0 * v) / 12.0e0;
        curlbasis[7][1] = (-3.0e0 - 6.0 * u + 2.0 * f64::sqrt(3.0e0) * v) / 12.0e0;
        curlbasis[7][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[7][s0] = -edgebasis[7][s0];
            }
            for s0 in 0..3 {
                curlbasis[7][s0] = -curlbasis[7][s0];
            }
        }
        i = edgemap[8][0];
        j = edgemap[8][1];
        edgebasis[8][0] = 0.0e0;
        edgebasis[8][1] = 0.0e0;
        edgebasis[8][2] = v * (-f64::sqrt(3.0e0) + 2.0 * v) / 6.0e0;
        curlbasis[8][0] = (-f64::sqrt(3.0e0) + 4.0 * v) / 6.0e0;
        curlbasis[8][1] = 0.0e0;
        curlbasis[8][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[8][s0] = -edgebasis[8][s0];
            }
            for s0 in 0..3 {
                curlbasis[8][s0] = -curlbasis[8][s0];
            }
        }
        prismsquarefacemap[0] = [1,2,5,4];
        prismsquarefacemap[1] = [2,3,6,5];
        prismsquarefacemap[2] = [3,1,4,6];
        workbasis[0][0] = (-3.0e0 + f64::sqrt(3.0e0) * v) * (-1.0e0 + sq(w)) / 6.0e0;
        workbasis[0][1] = -(u * (-1.0e0 + sq(w)) / (2.0e0 * f64::sqrt(3.0e0)));
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = u * w / f64::sqrt(3.0e0);
        workcurlbasis[0][1] = (-1.0e0 + v / f64::sqrt(3.0e0)) * w;
        workcurlbasis[0][2] = -((-1.0e0 + sq(w)) / f64::sqrt(3.0e0));
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = (3.0e0 - 3.0 * sq(u) - 2.0 * f64::sqrt(3.0e0) * v + sq(v)) / 6.0e0;
        workcurlbasis[1][0] = (-f64::sqrt(3.0e0) + v) / 3.0e0;
        workcurlbasis[1][1] = u;
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[prismsquarefacemap[0][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[9][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[9][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[10][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[10][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = v * (-1.0e0 + sq(w)) / (2.0e0 * f64::sqrt(3.0e0));
        workbasis[0][1] = -((1.0e0 + u) * (-1.0e0 + sq(w)) / (2.0e0 * f64::sqrt(3.0e0)));
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = (1.0e0 + u) * w / f64::sqrt(3.0e0);
        workcurlbasis[0][1] = v * w / f64::sqrt(3.0e0);
        workcurlbasis[0][2] = -((-1.0e0 + sq(w)) / f64::sqrt(3.0e0));
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - v) * v / 3.0e0;
        workcurlbasis[1][0] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 2.0 * v) / 3.0e0;
        workcurlbasis[1][1] = -(v / f64::sqrt(3.0e0));
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[prismsquarefacemap[1][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[11][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[11][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[12][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[12][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = v * (-1.0e0 + sq(w)) / (2.0e0 * f64::sqrt(3.0e0));
        workbasis[0][1] = -((-1.0e0 + u) * (-1.0e0 + sq(w)) / (2.0e0 * f64::sqrt(3.0e0)));
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = (-1.0e0 + u) * w / f64::sqrt(3.0e0);
        workcurlbasis[0][1] = v * w / f64::sqrt(3.0e0);
        workcurlbasis[0][2] = -((-1.0e0 + sq(w)) / f64::sqrt(3.0e0));
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = -(v * (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + v) / 3.0e0);
        workcurlbasis[1][0] = (f64::sqrt(3.0e0) - f64::sqrt(3.0e0) * u - 2.0 * v) / 3.0e0;
        workcurlbasis[1][1] = v / f64::sqrt(3.0e0);
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[prismsquarefacemap[2][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[13][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[13][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[14][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[14][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
    }
    Ok(())
}

fn brick_family(second_order: bool, u: f64, v: f64, w: f64, gindexes: &[i64],
    edgebasis: &mut Mat, curlbasis: &mut Mat) -> Result<(), String> {
    let (mut i, mut j, mut k, mut l, mut q) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut i1, mut i2, mut i3) = (0usize, 0usize, 0usize);
    let (mut d1, mut d2) = (0.0f64, 0.0f64);
    let (mut sfun, mut tfun, mut hfun) = (0.0f64, 0.0f64, 0.0f64);
    let (mut h1, mut h2, mut h3, mut dh1, mut dh2, mut dh3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut grad_sfun = [0.0f64; 3];
    let mut grad_tfun = [0.0f64; 3];
    let mut grad_hfun = [0.0f64; 3];
    let mut workbasis = [[0.0f64; 3]; 8];
    let mut workcurlbasis = [[0.0f64; 3]; 8];
    let mut lbasis = [0.0f64; 8];
    let mut beta = [0.0f64; 4];
    let mut edgesign = [0.0f64; 16];
    let mut grad = [0.0f64; 2];
    let mut faceindices = [0i64; 4];
    let mut trianglefacemap = [1usize, 2, 3];
    let mut squarefacemap = [1usize, 2, 3, 4];
    let mut brickfacemap = [[0usize; 4]; 6];
    let mut prismsquarefacemap = [[0usize; 4]; 3];
    let fs1 = 28.0f64;
    let fs2 = 84.0f64;
    let redefinefacebasis = true;
    let scalefacebasis = true;

    let secondorder = second_order;
    let edgemap: [[usize; 2]; 12] = [[1, 2], [2, 3], [4, 3], [1, 4], [5, 6], [6, 7], [8, 7], [5, 8], [1, 5], [2, 6], [3, 7], [4, 8]];

    if secondorder {
        for q in 1..=2 {
            k = 2 * q - 1;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][0] = 0.5e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][1] = 0.5e0 * -0.5e0 * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][2] = -(0.5e0 * line_nodal_p_basis(1, w) * d_line_nodal_p_basis(q, v));
            edgebasis[2 * (k - 1) + 2 - 1][0] = 1.5e0 * line_nodal_p_basis(1, w) * u * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][1] = 1.5e0 * -0.5e0 * u * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][2] = -(1.5e0 * line_nodal_p_basis(1, w) * u * d_line_nodal_p_basis(q, v));
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            k = 4 + 2 * q - 1;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][0] = 0.5e0 * line_nodal_p_basis(2, w) * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][1] = 0.5e0 * 0.5e0 * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][2] = -(0.5e0 * line_nodal_p_basis(2, w) * d_line_nodal_p_basis(q, v));
            edgebasis[2 * (k - 1) + 2 - 1][0] = 1.5e0 * line_nodal_p_basis(2, w) * u * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][1] = 1.5e0 * 0.5e0 * u * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][2] = -(1.5e0 * line_nodal_p_basis(2, w) * u * d_line_nodal_p_basis(q, v));
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            k = 2 * q;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][1] = 0.5e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(3 - q, u);
            curlbasis[2 * (k - 1) + 1 - 1][0] = -(0.5e0 * -0.5e0 * line_nodal_p_basis(3 - q, u));
            curlbasis[2 * (k - 1) + 1 - 1][2] = 0.5e0 * line_nodal_p_basis(1, w) * d_line_nodal_p_basis(3 - q, u);
            edgebasis[2 * (k - 1) + 2 - 1][1] = 1.5e0 * line_nodal_p_basis(1, w) * v * line_nodal_p_basis(3 - q, u);
            curlbasis[2 * (k - 1) + 2 - 1][0] = -(1.5e0 * -0.5e0 * v * line_nodal_p_basis(3 - q, u));
            curlbasis[2 * (k - 1) + 2 - 1][2] = 1.5e0 * line_nodal_p_basis(1, w) * v * d_line_nodal_p_basis(3 - q, u);
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            k = 4 + 2 * q;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][1] = 0.5e0 * line_nodal_p_basis(2, w) * line_nodal_p_basis(3 - q, u);
            curlbasis[2 * (k - 1) + 1 - 1][0] = -(0.5e0 * 0.5e0 * line_nodal_p_basis(3 - q, u));
            curlbasis[2 * (k - 1) + 1 - 1][2] = 0.5e0 * line_nodal_p_basis(2, w) * d_line_nodal_p_basis(3 - q, u);
            edgebasis[2 * (k - 1) + 2 - 1][1] = 1.5e0 * line_nodal_p_basis(2, w) * v * line_nodal_p_basis(3 - q, u);
            curlbasis[2 * (k - 1) + 2 - 1][0] = -(1.5e0 * 0.5e0 * v * line_nodal_p_basis(3 - q, u));
            curlbasis[2 * (k - 1) + 2 - 1][2] = 1.5e0 * line_nodal_p_basis(2, w) * v * d_line_nodal_p_basis(3 - q, u);
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            k = 8 + 3 * (q - 1) + 1;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][2] = 0.5e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][0] = 0.5e0 * line_nodal_p_basis(1, u) * d_line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][1] = -(0.5e0 * d_line_nodal_p_basis(1, u) * line_nodal_p_basis(q, v));
            edgebasis[2 * (k - 1) + 2 - 1][2] = 1.5e0 * line_nodal_p_basis(1, u) * w * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][0] = 1.5e0 * line_nodal_p_basis(1, u) * w * d_line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][1] = -(1.5e0 * d_line_nodal_p_basis(1, u) * w * line_nodal_p_basis(q, v));
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            k = 9 + q;
            i = edgemap[k - 1][0];
            j = edgemap[k - 1][1];
            edgebasis[2 * (k - 1) + 1 - 1][2] = 0.5e0 * line_nodal_p_basis(2, u) * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][0] = 0.5e0 * line_nodal_p_basis(2, u) * d_line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 1 - 1][1] = -(0.5e0 * d_line_nodal_p_basis(2, u) * line_nodal_p_basis(q, v));
            edgebasis[2 * (k - 1) + 2 - 1][2] = 1.5e0 * line_nodal_p_basis(2, u) * w * line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][0] = 1.5e0 * line_nodal_p_basis(2, u) * w * d_line_nodal_p_basis(q, v);
            curlbasis[2 * (k - 1) + 2 - 1][1] = -(1.5e0 * d_line_nodal_p_basis(2, u) * w * line_nodal_p_basis(q, v));
            if gindexes[j - 1] < gindexes[i - 1] {
                for s0 in 0..3 {
                    edgebasis[2 * (k - 1) + 1 - 1][s0] = -edgebasis[2 * (k - 1) + 1 - 1][s0];
                }
                for s0 in 0..3 {
                    curlbasis[2 * (k - 1) + 1 - 1][s0] = -curlbasis[2 * (k - 1) + 1 - 1][s0];
                }
            }
        }
        for q in 1..=2 {
            match q {
                1 => {
                    squarefacemap = [1,2,3,4];
                }
                2 => {
                    squarefacemap = [5,6,7,8];
                }
                _ => {}
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workbasis[s0][s1] = 0.0e0;
                }
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workcurlbasis[s0][s1] = 0.0e0;
                }
            }
            workbasis[0][0] = 2.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * line_nodal_p_basis(q, w);
            workcurlbasis[0][1] = 2.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * d_line_nodal_p_basis(q, w);
            workcurlbasis[0][2] = v * line_nodal_p_basis(q, w);
            workbasis[1][0] = 12.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * u * line_nodal_p_basis(q, w);
            workcurlbasis[1][1] = 12.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * u * d_line_nodal_p_basis(q, w);
            workcurlbasis[1][2] = -(12.0e0 * -(0.5e0 * v) * u * d_line_nodal_p_basis(q, w));
            workbasis[2][1] = 2.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * line_nodal_p_basis(q, w);
            workcurlbasis[2][0] = -(2.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * d_line_nodal_p_basis(q, w));
            workcurlbasis[2][2] = -(u * line_nodal_p_basis(q, w));
            workbasis[3][1] = 12.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * v * line_nodal_p_basis(q, w);
            workcurlbasis[3][0] = -(12.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * v * d_line_nodal_p_basis(q, w));
            workcurlbasis[3][2] = 12.0e0 * -(0.5e0 * u) * v * line_nodal_p_basis(q, w);
            for s0 in 0..4 {
                faceindices[s0] = gindexes[squarefacemap[s0] - 1];
            }
            (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
            k = 24;
            for s0 in 0..3 {
                edgebasis[k + 4 * (q - 1) + 1 - 1][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 * (q - 1) + 1 - 1][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 4 * (q - 1) + 2 - 1][s0] = 0.5e0 * workbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 * (q - 1) + 2 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 4 * (q - 1) + 3 - 1][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 * (q - 1) + 3 - 1][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 4 * (q - 1) + 4 - 1][s0] = 0.5e0 * workbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 * (q - 1) + 4 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
        }
        for q in 1..=2 {
            match q {
                1 => {
                    squarefacemap = [1,2,6,5];
                    k = 32;
                }
                2 => {
                    squarefacemap = [4,3,7,8];
                    k = 40;
                }
                _ => {}
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workbasis[s0][s1] = 0.0e0;
                }
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workcurlbasis[s0][s1] = 0.0e0;
                }
            }
            workbasis[0][0] = 2.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * line_nodal_p_basis(q, v);
            workcurlbasis[0][1] = 2.0e0 * -(0.5e0 * w) * line_nodal_p_basis(q, v);
            workcurlbasis[0][2] = -(2.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * d_line_nodal_p_basis(q, v));
            workbasis[1][0] = 12.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * u * line_nodal_p_basis(q, v);
            workcurlbasis[1][1] = 12.0e0 * -(0.5e0 * w) * u * line_nodal_p_basis(q, v);
            workcurlbasis[1][2] = -(12.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * u * d_line_nodal_p_basis(q, v));
            workbasis[2][2] = 2.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * line_nodal_p_basis(q, v);
            workcurlbasis[2][0] = 2.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * d_line_nodal_p_basis(q, v);
            workcurlbasis[2][1] = u * line_nodal_p_basis(q, v);
            workbasis[3][2] = 12.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * w * line_nodal_p_basis(q, v);
            workcurlbasis[3][0] = 12.0e0 * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u) * w * d_line_nodal_p_basis(q, v);
            workcurlbasis[3][1] = -(12.0e0 * -(0.5e0 * u) * w * line_nodal_p_basis(q, v));
            for s0 in 0..4 {
                faceindices[s0] = gindexes[squarefacemap[s0] - 1];
            }
            (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
            for s0 in 0..3 {
                edgebasis[k + 1 - 1][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 1 - 1][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 2 - 1][s0] = 0.5e0 * workbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 2 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 3 - 1][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 3 - 1][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 4 - 1][s0] = 0.5e0 * workbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
        }
        for q in 1..=2 {
            match q {
                1 => {
                    squarefacemap = [1,4,8,5];
                    k = 44;
                }
                2 => {
                    squarefacemap = [2,3,7,6];
                    k = 36;
                }
                _ => {}
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workbasis[s0][s1] = 0.0e0;
                }
            }
            for s0 in 0..8 {
                for s1 in 0..3 {
                    workcurlbasis[s0][s1] = 0.0e0;
                }
            }
            workbasis[0][1] = 2.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * line_nodal_p_basis(q, u);
            workcurlbasis[0][0] = -(2.0e0 * -(0.5e0 * w) * line_nodal_p_basis(q, u));
            workcurlbasis[0][2] = 2.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * d_line_nodal_p_basis(q, u);
            workbasis[1][1] = 12.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * v * line_nodal_p_basis(q, u);
            workcurlbasis[1][0] = -(12.0e0 * -(0.5e0 * w) * v * line_nodal_p_basis(q, u));
            workcurlbasis[1][2] = 12.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * v * d_line_nodal_p_basis(q, u);
            workbasis[2][2] = 2.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * line_nodal_p_basis(q, u);
            workcurlbasis[2][0] = 2.0e0 * -(0.5e0 * v) * line_nodal_p_basis(q, u);
            workcurlbasis[2][1] = -(2.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * d_line_nodal_p_basis(q, u));
            workbasis[3][2] = 12.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * w * line_nodal_p_basis(q, u);
            workcurlbasis[3][0] = 12.0e0 * -(0.5e0 * v) * w * line_nodal_p_basis(q, u);
            workcurlbasis[3][1] = -(12.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * w * d_line_nodal_p_basis(q, u));
            for s0 in 0..4 {
                faceindices[s0] = gindexes[squarefacemap[s0] - 1];
            }
            (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
            for s0 in 0..3 {
                edgebasis[k + 1 - 1][s0] = d1 * workbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 1 - 1][s0] = d1 * workcurlbasis[2 * (i1 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 2 - 1][s0] = 0.5e0 * workbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 2 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i1 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 3 - 1][s0] = d2 * workbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 3 - 1][s0] = d2 * workcurlbasis[2 * (i2 - 1) + 1 - 1][s0];
            }
            for s0 in 0..3 {
                edgebasis[k + 4 - 1][s0] = 0.5e0 * workbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
            for s0 in 0..3 {
                curlbasis[k + 4 - 1][s0] = 0.5e0 * workcurlbasis[2 * (i2 - 1) + 2 - 1][s0];
            }
        }
        edgebasis[48][0] = 8.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v);
        curlbasis[48][1] = 8.0e0 * -(0.5e0 * w) * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v);
        curlbasis[48][2] = -(8.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * -(0.5e0 * v));
        edgebasis[49][0] = 24.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * u * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v);
        curlbasis[49][1] = 24.0e0 * -(0.5e0 * w) * u * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v);
        curlbasis[49][2] = -(24.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * u * -(0.5e0 * v));
        edgebasis[50][1] = 8.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[50][0] = -(8.0e0 * -(0.5e0 * w) * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u));
        curlbasis[50][2] = 8.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * -(0.5e0 * u);
        edgebasis[51][1] = 24.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * v * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[51][0] = -(24.0e0 * -(0.5e0 * w) * v * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u));
        curlbasis[51][2] = 24.0e0 * line_nodal_p_basis(1, w) * line_nodal_p_basis(2, w) * v * -(0.5e0 * u);
        edgebasis[52][2] = 8.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[52][0] = 8.0e0 * -(0.5e0 * v) * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[52][1] = -(8.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * -(0.5e0 * u));
        edgebasis[53][2] = 24.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * w * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[53][0] = 24.0e0 * -(0.5e0 * v) * w * line_nodal_p_basis(1, u) * line_nodal_p_basis(2, u);
        curlbasis[53][1] = -(24.0e0 * line_nodal_p_basis(1, v) * line_nodal_p_basis(2, v) * w * -(0.5e0 * u));
    } else {
        i = edgemap[0][0];
        j = edgemap[0][1];
        edgebasis[0][0] = (-1.0e0 + v) * v * (-1.0e0 + w) * w / 8.0e0;
        edgebasis[0][1] = 0.0e0;
        edgebasis[0][2] = 0.0e0;
        curlbasis[0][0] = 0.0e0;
        curlbasis[0][1] = (-1.0e0 + v) * v * (-1.0e0 + 2.0 * w) / 8.0e0;
        curlbasis[0][2] = -((-1.0e0 + 2.0 * v) * (-1.0e0 + w) * w / 8.0e0);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[0][s0] = -edgebasis[0][s0];
            }
            for s0 in 0..3 {
                curlbasis[0][s0] = -curlbasis[0][s0];
            }
        }
        i = edgemap[1][0];
        j = edgemap[1][1];
        edgebasis[1][0] = 0.0e0;
        edgebasis[1][1] = u * (1.0e0 + u) * (-1.0e0 + w) * w / 8.0e0;
        edgebasis[1][2] = 0.0e0;
        curlbasis[1][0] = -(u * (1.0e0 + u) * (-1.0e0 + 2.0 * w) / 8.0e0);
        curlbasis[1][1] = 0.0e0;
        curlbasis[1][2] = (1.0e0 + 2.0 * u) * (-1.0e0 + w) * w / 8.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[1][s0] = -edgebasis[1][s0];
            }
            for s0 in 0..3 {
                curlbasis[1][s0] = -curlbasis[1][s0];
            }
        }
        i = edgemap[2][0];
        j = edgemap[2][1];
        edgebasis[2][0] = v * (1.0e0 + v) * (-1.0e0 + w) * w / 8.0e0;
        edgebasis[2][1] = 0.0e0;
        edgebasis[2][2] = 0.0e0;
        curlbasis[2][0] = 0.0e0;
        curlbasis[2][1] = v * (1.0e0 + v) * (-1.0e0 + 2.0 * w) / 8.0e0;
        curlbasis[2][2] = -((1.0e0 + 2.0 * v) * (-1.0e0 + w) * w / 8.0e0);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[2][s0] = -edgebasis[2][s0];
            }
            for s0 in 0..3 {
                curlbasis[2][s0] = -curlbasis[2][s0];
            }
        }
        i = edgemap[3][0];
        j = edgemap[3][1];
        edgebasis[3][0] = 0.0e0;
        edgebasis[3][1] = (-1.0e0 + u) * u * (-1.0e0 + w) * w / 8.0e0;
        edgebasis[3][2] = 0.0e0;
        curlbasis[3][0] = -((-1.0e0 + u) * u * (-1.0e0 + 2.0 * w) / 8.0e0);
        curlbasis[3][1] = 0.0e0;
        curlbasis[3][2] = (-1.0e0 + 2.0 * u) * (-1.0e0 + w) * w / 8.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[3][s0] = -edgebasis[3][s0];
            }
            for s0 in 0..3 {
                curlbasis[3][s0] = -curlbasis[3][s0];
            }
        }
        i = edgemap[4][0];
        j = edgemap[4][1];
        edgebasis[4][0] = (-1.0e0 + v) * v * w * (1.0e0 + w) / 8.0e0;
        edgebasis[4][1] = 0.0e0;
        edgebasis[4][2] = 0.0e0;
        curlbasis[4][0] = 0.0e0;
        curlbasis[4][1] = (-1.0e0 + v) * v * (1.0e0 + 2.0 * w) / 8.0e0;
        curlbasis[4][2] = -((-1.0e0 + 2.0 * v) * w * (1.0e0 + w) / 8.0e0);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[4][s0] = -edgebasis[4][s0];
            }
            for s0 in 0..3 {
                curlbasis[4][s0] = -curlbasis[4][s0];
            }
        }
        i = edgemap[5][0];
        j = edgemap[5][1];
        edgebasis[5][0] = 0.0e0;
        edgebasis[5][1] = u * (1.0e0 + u) * w * (1.0e0 + w) / 8.0e0;
        edgebasis[5][2] = 0.0e0;
        curlbasis[5][0] = -(u * (1.0e0 + u) * (1.0e0 + 2.0 * w) / 8.0e0);
        curlbasis[5][1] = 0.0e0;
        curlbasis[5][2] = (1.0e0 + 2.0 * u) * w * (1.0e0 + w) / 8.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[5][s0] = -edgebasis[5][s0];
            }
            for s0 in 0..3 {
                curlbasis[5][s0] = -curlbasis[5][s0];
            }
        }
        i = edgemap[6][0];
        j = edgemap[6][1];
        edgebasis[6][0] = v * (1.0e0 + v) * w * (1.0e0 + w) / 8.0e0;
        edgebasis[6][1] = 0.0e0;
        edgebasis[6][2] = 0.0e0;
        curlbasis[6][0] = 0.0e0;
        curlbasis[6][1] = v * (1.0e0 + v) * (1.0e0 + 2.0 * w) / 8.0e0;
        curlbasis[6][2] = -((1.0e0 + 2.0 * v) * w * (1.0e0 + w) / 8.0e0);
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[6][s0] = -edgebasis[6][s0];
            }
            for s0 in 0..3 {
                curlbasis[6][s0] = -curlbasis[6][s0];
            }
        }
        i = edgemap[7][0];
        j = edgemap[7][1];
        edgebasis[7][0] = 0.0e0;
        edgebasis[7][1] = (-1.0e0 + u) * u * w * (1.0e0 + w) / 8.0e0;
        edgebasis[7][2] = 0.0e0;
        curlbasis[7][0] = -((-1.0e0 + u) * u * (1.0e0 + 2.0 * w) / 8.0e0);
        curlbasis[7][1] = 0.0e0;
        curlbasis[7][2] = (-1.0e0 + 2.0 * u) * w * (1.0e0 + w) / 8.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[7][s0] = -edgebasis[7][s0];
            }
            for s0 in 0..3 {
                curlbasis[7][s0] = -curlbasis[7][s0];
            }
        }
        i = edgemap[8][0];
        j = edgemap[8][1];
        edgebasis[8][0] = 0.0e0;
        edgebasis[8][1] = 0.0e0;
        edgebasis[8][2] = (-1.0e0 + u) * u * (-1.0e0 + v) * v / 8.0e0;
        curlbasis[8][0] = (-1.0e0 + u) * u * (-1.0e0 + 2.0 * v) / 8.0e0;
        curlbasis[8][1] = -((-1.0e0 + 2.0 * u) * (-1.0e0 + v) * v / 8.0e0);
        curlbasis[8][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[8][s0] = -edgebasis[8][s0];
            }
            for s0 in 0..3 {
                curlbasis[8][s0] = -curlbasis[8][s0];
            }
        }
        i = edgemap[9][0];
        j = edgemap[9][1];
        edgebasis[9][0] = 0.0e0;
        edgebasis[9][1] = 0.0e0;
        edgebasis[9][2] = u * (1.0e0 + u) * (-1.0e0 + v) * v / 8.0e0;
        curlbasis[9][0] = u * (1.0e0 + u) * (-1.0e0 + 2.0 * v) / 8.0e0;
        curlbasis[9][1] = -((1.0e0 + 2.0 * u) * (-1.0e0 + v) * v / 8.0e0);
        curlbasis[9][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[9][s0] = -edgebasis[9][s0];
            }
            for s0 in 0..3 {
                curlbasis[9][s0] = -curlbasis[9][s0];
            }
        }
        i = edgemap[10][0];
        j = edgemap[10][1];
        edgebasis[10][0] = 0.0e0;
        edgebasis[10][1] = 0.0e0;
        edgebasis[10][2] = u * (1.0e0 + u) * v * (1.0e0 + v) / 8.0e0;
        curlbasis[10][0] = u * (1.0e0 + u) * (1.0e0 + 2.0 * v) / 8.0e0;
        curlbasis[10][1] = -((1.0e0 + 2.0 * u) * v * (1.0e0 + v) / 8.0e0);
        curlbasis[10][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[10][s0] = -edgebasis[10][s0];
            }
            for s0 in 0..3 {
                curlbasis[10][s0] = -curlbasis[10][s0];
            }
        }
        i = edgemap[11][0];
        j = edgemap[11][1];
        edgebasis[11][0] = 0.0e0;
        edgebasis[11][1] = 0.0e0;
        edgebasis[11][2] = (-1.0e0 + u) * u * v * (1.0e0 + v) / 8.0e0;
        curlbasis[11][0] = (-1.0e0 + u) * u * (1.0e0 + 2.0 * v) / 8.0e0;
        curlbasis[11][1] = -((-1.0e0 + 2.0 * u) * v * (1.0e0 + v) / 8.0e0);
        curlbasis[11][2] = 0.0e0;
        if gindexes[j - 1] < gindexes[i - 1] {
            for s0 in 0..3 {
                edgebasis[11][s0] = -edgebasis[11][s0];
            }
            for s0 in 0..3 {
                curlbasis[11][s0] = -curlbasis[11][s0];
            }
        }
        brickfacemap[0] = [1,2,3,4];
        brickfacemap[1] = [5,6,7,8];
        brickfacemap[2] = [1,2,6,5];
        brickfacemap[3] = [2,3,7,6];
        brickfacemap[4] = [4,3,7,8];
        brickfacemap[5] = [1,4,8,5];
        workbasis[0][0] = -((-1.0e0 + sq(v)) * (-1.0e0 + w) * w / 4.0e0);
        workbasis[0][1] = 0.0e0;
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = 0.0e0;
        workcurlbasis[0][1] = -((-1.0e0 + sq(v)) * (-1.0e0 + 2.0 * w) / 4.0e0);
        workcurlbasis[0][2] = v * (-1.0e0 + w) * w / 2.0e0;
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = -((-1.0e0 + sq(u)) * (-1.0e0 + w) * w / 4.0e0);
        workbasis[1][2] = 0.0e0;
        workcurlbasis[1][0] = (-1.0e0 + sq(u)) * (-1.0e0 + 2.0 * w) / 4.0e0;
        workcurlbasis[1][1] = 0.0e0;
        workcurlbasis[1][2] = -(u * (-1.0e0 + w) * w / 2.0e0);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[0][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[12][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[12][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[13][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[13][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = -((-1.0e0 + sq(v)) * w * (1.0e0 + w) / 4.0e0);
        workbasis[0][1] = 0.0e0;
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = 0.0e0;
        workcurlbasis[0][1] = -((-1.0e0 + sq(v)) * (1.0e0 + 2.0 * w) / 4.0e0);
        workcurlbasis[0][2] = v * w * (1.0e0 + w) / 2.0e0;
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = -((-1.0e0 + sq(u)) * w * (1.0e0 + w) / 4.0e0);
        workbasis[1][2] = 0.0e0;
        workcurlbasis[1][0] = (-1.0e0 + sq(u)) * (1.0e0 + 2.0 * w) / 4.0e0;
        workcurlbasis[1][1] = 0.0e0;
        workcurlbasis[1][2] = -(u * w * (1.0e0 + w) / 2.0e0);
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[1][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[14][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[14][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[15][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[15][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = -((-1.0e0 + v) * v * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[0][1] = 0.0e0;
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = 0.0e0;
        workcurlbasis[0][1] = -((-1.0e0 + v) * v * w / 2.0e0);
        workcurlbasis[0][2] = (-1.0e0 + 2.0 * v) * (-1.0e0 + sq(w)) / 4.0e0;
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = -((-1.0e0 + sq(u)) * (-1.0e0 + v) * v / 4.0e0);
        workcurlbasis[1][0] = -((-1.0e0 + sq(u)) * (-1.0e0 + 2.0 * v) / 4.0e0);
        workcurlbasis[1][1] = u * (-1.0e0 + v) * v / 2.0e0;
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[2][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[16][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[16][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[17][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[17][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = 0.0e0;
        workbasis[0][1] = -(u * (1.0e0 + u) * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = u * (1.0e0 + u) * w / 2.0e0;
        workcurlbasis[0][1] = 0.0e0;
        workcurlbasis[0][2] = -((1.0e0 + 2.0 * u) * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = -(u * (1.0e0 + u) * ((-1 as f64) + sq(v)) / 4.0e0);
        workcurlbasis[1][0] = -(u * (1.0e0 + u) * v / 2.0e0);
        workcurlbasis[1][1] = (1.0e0 + 2.0 * u) * (-1.0e0 + sq(v)) / 4.0e0;
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[3][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[18][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[18][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[19][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[19][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = -(v * (1.0e0 + v) * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[0][1] = 0.0e0;
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = 0.0e0;
        workcurlbasis[0][1] = -(v * (1.0e0 + v) * w / 2.0e0);
        workcurlbasis[0][2] = (1.0e0 + 2.0 * v) * (-1.0e0 + sq(w)) / 4.0e0;
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = -((-1.0e0 + sq(u)) * v * (1.0e0 + v) / 4.0e0);
        workcurlbasis[1][0] = -((-1.0e0 + sq(u)) * (1.0e0 + 2.0 * v) / 4.0e0);
        workcurlbasis[1][1] = u * v * (1.0e0 + v) / 2.0e0;
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[4][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[20][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[20][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[21][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[21][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        workbasis[0][0] = 0.0e0;
        workbasis[0][1] = -((-1.0e0 + u) * u * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[0][2] = 0.0e0;
        workcurlbasis[0][0] = (-1.0e0 + u) * u * w / 2.0e0;
        workcurlbasis[0][1] = 0.0e0;
        workcurlbasis[0][2] = -((-1.0e0 + 2.0 * u) * (-1.0e0 + sq(w)) / 4.0e0);
        workbasis[1][0] = 0.0e0;
        workbasis[1][1] = 0.0e0;
        workbasis[1][2] = -((-1.0e0 + u) * u * (-1.0e0 + sq(v)) / 4.0e0);
        workcurlbasis[1][0] = -((-1.0e0 + u) * u * v / 2.0e0);
        workcurlbasis[1][1] = (-1.0e0 + 2.0 * u) * (-1.0e0 + sq(v)) / 4.0e0;
        workcurlbasis[1][2] = 0.0e0;
        for s0 in 0..4 {
            faceindices[s0] = gindexes[brickfacemap[5][s0] - 1];
        }
        (i1, i2, d1, d2) = square_face_dofs_ordering(&faceindices);
        for s0 in 0..3 {
            edgebasis[22][s0] = d1 * workbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[22][s0] = d1 * workcurlbasis[i1 - 1][s0];
        }
        for s0 in 0..3 {
            edgebasis[23][s0] = d2 * workbasis[i2 - 1][s0];
        }
        for s0 in 0..3 {
            curlbasis[23][s0] = d2 * workcurlbasis[i2 - 1][s0];
        }
        edgebasis[24][0] = (-1.0e0 + sq(v)) * (-1.0e0 + sq(w)) / 2.0e0;
        edgebasis[24][1] = 0.0e0;
        edgebasis[24][2] = 0.0e0;
        curlbasis[24][0] = 0.0e0;
        curlbasis[24][1] = (-1.0e0 + sq(v)) * w;
        curlbasis[24][2] = v - v * sq(w);
        edgebasis[25][0] = 0.0e0;
        edgebasis[25][1] = (-1.0e0 + sq(u)) * (-1.0e0 + sq(w)) / 2.0e0;
        edgebasis[25][2] = 0.0e0;
        curlbasis[25][0] = w - sq(u) * w;
        curlbasis[25][1] = 0.0e0;
        curlbasis[25][2] = u * ((-1 as f64) + sq(w));
        edgebasis[26][0] = 0.0e0;
        edgebasis[26][1] = 0.0e0;
        edgebasis[26][2] = (-1.0e0 + sq(u)) * (-1.0e0 + sq(v)) / 2.0e0;
        curlbasis[26][0] = (-1.0e0 + sq(u)) * v;
        curlbasis[26][1] = u - u * sq(v);
        curlbasis[26][2] = 0.0e0;
    }
    Ok(())
}

fn weighted_whitney_forms(workbasis: &mut [[f64; 3]; 8], workcurlbasis: &mut [[f64; 3]; 8], k: usize, u: f64, v: f64, w: f64) {
    match k {
        1 => {
            workbasis[0][0] = (4.0e0 * v - f64::sqrt(2.0e0) * w) * (-6.0e0 + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (48.0e0 * f64::sqrt(3.0e0));
            workbasis[0][1] = -(u * (4.0e0 * v - f64::sqrt(2.0e0) * w) / 24.0e0);
            workbasis[0][2] = u * (-(2.0e0 * f64::sqrt(2.0e0) * v) + w) / 24.0e0;
            workcurlbasis[0][0] = -(u / (4.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[0][1] = (f64::sqrt(6.0e0) + 3.0e0 * f64::sqrt(2.0e0) * v - 3.0e0 * w) / 24.0e0;
            workcurlbasis[0][2] = (f64::sqrt(3.0e0) - 3.0e0 * v) / 6.0e0;
            workbasis[1][0] = (4.0e0 * v - f64::sqrt(2.0e0) * w) * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (96.0e0 * f64::sqrt(3.0e0));
            workbasis[1][1] = -((4.0e0 * f64::sqrt(3.0e0) + 4.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * f64::sqrt(2.0e0) * w) * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / 288.0e0);
            workbasis[1][2] = (f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u - 3.0e0 * v) * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (144.0e0 * f64::sqrt(2.0e0));
            workcurlbasis[1][0] = -((-6.0e0 + 2.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (16.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[1][1] = (2.0e0 * f64::sqrt(3.0e0) - 6.0e0 * f64::sqrt(3.0e0) * u + 6.0e0 * v - 3.0e0 * f64::sqrt(2.0e0) * w) / (48.0e0 * f64::sqrt(2.0e0));
            workcurlbasis[1][2] = (f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * v) / 12.0e0;
            workbasis[2][0] = -((4.0e0 * v - f64::sqrt(2.0e0) * w) * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (96.0e0 * f64::sqrt(3.0e0)));
            workbasis[2][1] = (-(4.0e0 * f64::sqrt(3.0e0)) + 4.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * f64::sqrt(2.0e0) * w) * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / 288.0e0;
            workbasis[2][2] = -((-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + 3.0e0 * v) * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (144.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[2][0] = -((-6.0e0 - 2.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (16.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[2][1] = (-(2.0e0 * f64::sqrt(3.0e0)) - 6.0e0 * f64::sqrt(3.0e0) * u - 6.0e0 * v + 3.0e0 * f64::sqrt(2.0e0) * w) / (48.0e0 * f64::sqrt(2.0e0));
            workcurlbasis[2][2] = (-f64::sqrt(3.0e0) - 3.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * v) / 12.0e0;
        }
        2 => {
            workbasis[0][0] = -(w * (-6.0e0 + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (8.0e0 * f64::sqrt(6.0e0)));
            workbasis[0][1] = u * w / (4.0e0 * f64::sqrt(2.0e0));
            workbasis[0][2] = u * w / 8.0e0;
            workcurlbasis[0][0] = -(u / (4.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[0][1] = (f64::sqrt(6.0e0) - f64::sqrt(2.0e0) * v - 3.0e0 * w) / 8.0e0;
            workcurlbasis[0][2] = w / (2.0e0 * f64::sqrt(2.0e0));
            workbasis[1][0] = -(w * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0)));
            workbasis[1][1] = w * (1.0e0 + u - v / f64::sqrt(3.0e0) - w / f64::sqrt(6.0e0)) / (8.0e0 * f64::sqrt(2.0e0));
            workbasis[1][2] = (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + v) * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (48.0e0 * f64::sqrt(2.0e0));
            workcurlbasis[1][0] = (-(3.0e0 * f64::sqrt(2.0e0)) - f64::sqrt(2.0e0) * u + f64::sqrt(6.0e0) * v + f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[1][1] = (f64::sqrt(6.0e0) + 3.0e0 * f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v - 3.0e0 * w) / 16.0e0;
            workcurlbasis[1][2] = w / (4.0e0 * f64::sqrt(2.0e0));
            workbasis[2][0] = w * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0));
            workbasis[2][1] = -(w * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (48.0e0 * f64::sqrt(2.0e0)));
            workbasis[2][2] = -((f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v) * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / 96.0e0);
            workcurlbasis[2][0] = (-(3.0e0 * f64::sqrt(2.0e0)) + f64::sqrt(2.0e0) * u + f64::sqrt(6.0e0) * v + f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[2][1] = (-f64::sqrt(6.0e0) + 3.0e0 * f64::sqrt(6.0e0) * u + f64::sqrt(2.0e0) * v + 3.0e0 * w) / 16.0e0;
            workcurlbasis[2][2] = -(w / (4.0e0 * f64::sqrt(2.0e0)));
        }
        3 => {
            workbasis[0][0] = w * (-(2.0e0 * f64::sqrt(2.0e0) * v) + w) / 16.0e0;
            workbasis[0][1] = w * (4.0e0 * f64::sqrt(3.0e0) + 4.0e0 * f64::sqrt(3.0e0) * u - 3.0e0 * f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0));
            workbasis[0][2] = -((1.0e0 + u - f64::sqrt(3.0e0) * v) * w / 16.0e0);
            workcurlbasis[0][0] = (-(2.0e0 * f64::sqrt(2.0e0)) - 2.0e0 * f64::sqrt(2.0e0) * u + 3.0e0 * f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[0][1] = (-(2.0e0 * f64::sqrt(2.0e0) * v) + 3.0e0 * w) / 16.0e0;
            workcurlbasis[0][2] = w / (2.0e0 * f64::sqrt(2.0e0));
            workbasis[1][0] = w * (-(2.0e0 * f64::sqrt(2.0e0) * v) + w) / 16.0e0;
            workbasis[1][1] = -(w * (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0)));
            workbasis[1][2] = -((f64::sqrt(6.0e0) + f64::sqrt(6.0e0) * u - f64::sqrt(2.0e0) * v) * (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (32.0e0 * f64::sqrt(3.0e0)));
            workcurlbasis[1][0] = (2.0e0 * f64::sqrt(2.0e0) + 2.0e0 * f64::sqrt(2.0e0) * u - 2.0e0 * f64::sqrt(6.0e0) * v + f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[1][1] = (-(4.0e0 * f64::sqrt(2.0e0) * v) + 3.0e0 * w) / 16.0e0;
            workcurlbasis[1][2] = w / (4.0e0 * f64::sqrt(2.0e0));
            workbasis[2][0] = 0.0e0;
            workbasis[2][1] = w * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (24.0e0 * f64::sqrt(2.0e0));
            workbasis[2][2] = -(v * (-6.0e0 - 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (24.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[2][0] = (2.0e0 * f64::sqrt(2.0e0) + 2.0e0 * f64::sqrt(2.0e0) * u - f64::sqrt(6.0e0) * v - f64::sqrt(3.0e0) * w) / 8.0e0;
            workcurlbasis[2][1] = -(v / (4.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[2][2] = -(w / (4.0e0 * f64::sqrt(2.0e0)));
        }
        4 => {
            workbasis[0][0] = w * (-(2.0e0 * f64::sqrt(2.0e0) * v) + w) / 16.0e0;
            workbasis[0][1] = w * (-(4.0e0 * f64::sqrt(3.0e0)) + 4.0e0 * f64::sqrt(3.0e0) * u + 3.0e0 * f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0));
            workbasis[0][2] = -((-1.0e0 + u + f64::sqrt(3.0e0) * v) * w / 16.0e0);
            workcurlbasis[0][0] = (2.0e0 * f64::sqrt(2.0e0) - 2.0e0 * f64::sqrt(2.0e0) * u - 3.0e0 * f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[0][1] = (-(2.0e0 * f64::sqrt(2.0e0) * v) + 3.0e0 * w) / 16.0e0;
            workcurlbasis[0][2] = w / (2.0e0 * f64::sqrt(2.0e0));
            workbasis[1][0] = 0.0e0;
            workbasis[1][1] = w * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (24.0e0 * f64::sqrt(2.0e0));
            workbasis[1][2] = -(v * (-6.0e0 + 6.0e0 * u + 2.0e0 * f64::sqrt(3.0e0) * v + f64::sqrt(6.0e0) * w) / (24.0e0 * f64::sqrt(2.0e0)));
            workcurlbasis[1][0] = (2.0e0 * f64::sqrt(2.0e0) - 2.0e0 * f64::sqrt(2.0e0) * u - f64::sqrt(6.0e0) * v - f64::sqrt(3.0e0) * w) / 8.0e0;
            workcurlbasis[1][1] = v / (4.0e0 * f64::sqrt(2.0e0));
            workcurlbasis[1][2] = w / (4.0e0 * f64::sqrt(2.0e0));
            workbasis[2][0] = (2.0e0 * f64::sqrt(2.0e0) * v - w) * w / 16.0e0;
            workbasis[2][1] = -(w * (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0)));
            workbasis[2][2] = (-f64::sqrt(3.0e0) + f64::sqrt(3.0e0) * u + v) * (-(4.0e0 * v) + f64::sqrt(2.0e0) * w) / (16.0e0 * f64::sqrt(6.0e0));
            workcurlbasis[2][0] = (2.0e0 * f64::sqrt(2.0e0) - 2.0e0 * f64::sqrt(2.0e0) * u - 2.0e0 * f64::sqrt(6.0e0) * v + f64::sqrt(3.0e0) * w) / 16.0e0;
            workcurlbasis[2][1] = (4.0e0 * f64::sqrt(2.0e0) * v - 3.0e0 * w) / 16.0e0;
            workcurlbasis[2][2] = -(w / (4.0e0 * f64::sqrt(2.0e0)));
        }
        _ => {
            panic!("weighted_whitney_forms: wrong face index");
        }
    }
}

// ---------------------------------------------------------------------------
// Result type and element data
// ---------------------------------------------------------------------------

/// Evaluated edge basis at one point.
#[derive(Clone, Debug, Default)]
pub struct EdgeEval {
    /// Number of H(curl) basis functions (`DOFs` / `SIZE(EdgeMap,1)`).
    pub dofs: usize,
    /// The integration metric returned by Elmer: `|det F|` (or `sqrt(det G)`
    /// for surface elements) for [`edge_element_info`], `DetJ` for
    /// [`get_edge_basis`].
    pub det: f64,
    /// `EdgeBasis(1:dofs,:)`.
    pub edge: Vec<[f64; 3]>,
    /// `RotBasis(1:dofs,:)` / `RotWBasis`.
    pub curl: Vec<[f64; 3]>,
}

fn supported(code: u16) -> bool {
    matches!(code, 303 | 404 | 504 | 605 | 706 | 808)
}

fn element_dim(code: u16) -> usize {
    match code / 100 {
        3 | 4 => 2,
        _ => 3,
    }
}

/// `Element % TYPE % NodeU/NodeV/NodeW` of the standard Elmer reference
/// elements (elements.def).
fn std_node_uvw(code: u16) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    match code {
        303 => (vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0], vec![0.0; 3]),
        404 => (
            vec![-1.0, 1.0, 1.0, -1.0],
            vec![-1.0, -1.0, 1.0, 1.0],
            vec![0.0; 4],
        ),
        504 => (
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ),
        605 => (
            vec![-1.0, 1.0, 1.0, -1.0, 0.0],
            vec![-1.0, -1.0, 1.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 1.0],
        ),
        706 => (
            vec![0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            vec![-1.0, -1.0, -1.0, 1.0, 1.0, 1.0],
        ),
        808 => (
            vec![-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
            vec![-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0],
            vec![-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0],
        ),
        _ => unreachable!(),
    }
}

/// `GetEdgeMap(family)` (ElementBasis.F90), 1-based local nodes.
fn get_edge_map(family: u16) -> &'static [[usize; 2]] {
    match family {
        3 => &[[1, 2], [2, 3], [3, 1]],
        4 => &[[1, 2], [2, 3], [3, 4], [4, 1]],
        5 => &[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]],
        6 => &[[1, 2], [2, 3], [4, 3], [1, 4], [1, 5], [2, 5], [3, 5], [4, 5]],
        7 => &[[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4], [1, 4], [2, 5], [3, 6]],
        8 => &[
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
        ],
        _ => &[],
    }
}

fn sum_xd(x: &[f64], d: &[[f64; 3]], i: usize, n: usize) -> f64 {
    let mut s = 0.0;
    for q in 0..n {
        s += x[q] * d[q][i];
    }
    s
}

/// `InvertMatrix3x3` (LinearAlgebra.F90); `gi` is the inverse of `g`.
fn invert_matrix_3x3(g: &[[f64; 3]; 3], det_g: f64) -> [[f64; 3]; 3] {
    // Fortran quirk: `s = 1.0 / DetG` with a default-REAL literal numerator.
    let s = 1.0 / det_g;
    let mut gi = [[0.0f64; 3]; 3];
    gi[0][0] = s * (g[1][1] * g[2][2] - g[2][1] * g[1][2]);
    gi[1][0] = -s * (g[1][0] * g[2][2] - g[2][0] * g[1][2]);
    gi[2][0] = s * (g[1][0] * g[2][1] - g[2][0] * g[1][1]);
    gi[0][1] = -s * (g[0][1] * g[2][2] - g[2][1] * g[0][2]);
    gi[1][1] = s * (g[0][0] * g[2][2] - g[2][0] * g[0][2]);
    gi[2][1] = -s * (g[0][0] * g[2][1] - g[2][0] * g[0][1]);
    gi[0][2] = s * (g[0][1] * g[1][2] - g[1][1] * g[0][2]);
    gi[1][2] = -s * (g[0][0] * g[1][2] - g[1][0] * g[0][2]);
    gi[2][2] = s * (g[0][0] * g[1][1] - g[1][0] * g[0][1]);
    gi
}

/// `ElementMetric` (ElemInfo.F90) for a Cartesian mesh: returns
/// `(Metric, DetG, LtoGMap)` or `None` for a degenerate element.
fn element_metric(
    dim: usize,
    cdim: usize,
    nodes: &[[f64; 3]],
    n: usize,
    dl: &[[f64; 3]],
) -> Option<([[f64; 3]; 3], f64, [[f64; 3]; 3])> {
    let eps = f64::EPSILON.powi(dim as i32);
    let mut dx = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for c in 0..3 {
            let mut s = 0.0;
            for q in 0..n {
                s += nodes[q][c] * dl[q][i];
            }
            dx[c][i] = s;
        }
    }
    let mut g = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            let mut s = 0.0;
            for k in 0..cdim {
                s += dx[k][i] * dx[k][j];
            }
            g[i][j] = s;
        }
    }
    let mut metric = [[0.0f64; 3]; 3];
    let det_g;
    match dim {
        1 => {
            let d = g[0][0];
            if d <= eps {
                return None;
            }
            metric[0][0] = 1.0 / d;
            det_g = d.sqrt();
        }
        2 => {
            let d = g[0][0] * g[1][1] - g[0][1] * g[1][0];
            if d <= eps {
                return None;
            }
            metric[0][0] = g[1][1] / d;
            metric[0][1] = -g[0][1] / d;
            metric[1][0] = -g[1][0] / d;
            metric[1][1] = g[0][0] / d;
            det_g = d.sqrt();
        }
        _ => {
            let d = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                + g[0][1] * (g[1][2] * g[2][0] - g[1][0] * g[2][2])
                + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
            if d <= eps {
                return None;
            }
            metric = invert_matrix_3x3(&g, d);
            det_g = d.sqrt();
        }
    }
    let mut ltog = [[0.0f64; 3]; 3];
    for i in 0..cdim {
        for j in 0..dim {
            let mut s = 0.0;
            for k in 0..dim {
                s += dx[i][k] * metric[k][j];
            }
            ltog[i][j] = s;
        }
    }
    Some((metric, det_g, ltog))
}

// ---------------------------------------------------------------------------
// EdgeElementInfo (Piola version)
// ---------------------------------------------------------------------------

/// Number of H(curl) basis functions of `EdgeElementInfo` (`DOFs`).
pub fn edge_dofs(etype: ElementType, degree: u8) -> Result<usize, String> {
    let second = match degree {
        1 => false,
        2 => true,
        _ => return Err(format!("edge basis degree {degree} is not supported")),
    };
    Ok(match (etype.0, second) {
        (303, false) => 3,
        (303, true) => 8,
        (404, false) => 6,
        (404, true) => 12,
        (504, false) => 6,
        (504, true) => 20,
        (605, false) => 10,
        (605, true) => 31,
        (706, false) => 15,
        (706, true) => 36,
        (808, false) => 27,
        (808, true) => 54,
        (c, _) => return Err(format!("edge element type {c} is not supported")),
    })
}

/// Reference-element basis functions `b` and `Curl b` of `EdgeElementInfo`
/// (before the Piola transformation), with edge/face orientations taken from
/// the global node numbers `gindex`.
pub fn edge_basis_reference(
    etype: ElementType,
    degree: u8,
    gindex: &[i64],
    uvw: [f64; 3],
) -> Result<EdgeEval, String> {
    let dofs = edge_dofs(etype, degree)?;
    let n = (etype.0 % 100) as usize;
    if gindex.len() < n {
        return Err("gindex shorter than the node count".into());
    }
    let [u, v, w] = uvw;
    let mut eb: Mat = [[0.0; 3]; 64];
    let mut cb: Mat = [[0.0; 3]; 64];
    let so = degree == 2;
    match etype.0 / 100 {
        3 => tri_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
        4 => quad_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
        5 => tetra_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
        6 => pyramid_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
        7 => prism_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
        _ => brick_family(so, u, v, w, gindex, &mut eb, &mut cb)?,
    }
    Ok(EdgeEval {
        dofs,
        det: 1.0,
        edge: eb[..dofs].to_vec(),
        curl: cb[..dofs].to_vec(),
    })
}

/// p-reference-element nodal basis and local derivatives used by
/// `EdgeElementInfo` for the element map.
fn p_ref_nodal(code: u16, u: f64, v: f64, w: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    let n = (code % 100) as usize;
    let mut b = vec![0.0; n];
    let mut d = vec![[0.0; 3]; n];
    for q in 1..=n {
        match code / 100 {
            3 => {
                b[q - 1] = triangle_nodal_p_basis(q, u, v);
                d[q - 1] = d_triangle_nodal_p_basis(q, u, v);
            }
            4 => {
                b[q - 1] = quad_nodal_p_basis(q, u, v);
                d[q - 1] = d_quad_nodal_p_basis(q, u, v);
            }
            5 => {
                b[q - 1] = tetra_nodal_p_basis(q, u, v, w);
                d[q - 1] = d_tetra_nodal_p_basis(q, u, v, w);
            }
            6 => {
                b[q - 1] = pyramid_nodal_p_basis(q, u, v, w);
                d[q - 1] = d_pyramid_nodal_p_basis(q, u, v, w);
            }
            7 => {
                b[q - 1] = wedge_nodal_p_basis(q, u, v, w);
                d[q - 1] = d_wedge_nodal_p_basis(q, u, v, w);
            }
            _ => {
                b[q - 1] = brick_nodal_p_basis(q, u, v, w);
                d[q - 1] = d_brick_nodal_p_basis(q, u, v, w);
            }
        }
    }
    (b, d)
}

/// `EdgeElementInfo` with `ApplyPiolaTransform = .TRUE.`: H(curl)-conforming
/// Nédélec basis of the first kind, degree `degree` (1 = lowest-order
/// "optimal" family, 2 = second-order hierarchic), on the physical element
/// with nodal coordinates `nodes`. `gindex` are the global node numbers
/// (`Element % NodeIndexes` in a serial run) that fix the orientation, and
/// `cdim` is `CoordinateSystemDimension()`.
///
/// Returns the basis values `B(f(p))`, the spatial curl and the metric
/// `|det F|` (`sqrt(det G)` for a surface element, `cdim > dim`) for
/// integration.
pub fn edge_element_info(
    etype: ElementType,
    nodes: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
    degree: u8,
    cdim: usize,
) -> Result<EdgeEval, String> {
    let code = etype.0;
    if !supported(code) {
        return Err(format!("edge element type {code} is not supported"));
    }
    let n = (code % 100) as usize;
    if nodes.len() < n {
        return Err("fewer nodal coordinates than element nodes".into());
    }
    let dim = element_dim(code);
    if cdim < dim {
        return Err("coordinate dimension smaller than the element dimension".into());
    }
    let mut ev = edge_basis_reference(etype, degree, gindex, uvw)?;
    let dofs = ev.dofs;
    let [u, v, w] = uvw;
    let (_basis, dl) = p_ref_nodal(code, u, v, w);

    // PiolaTransformationData
    let mut lf = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for c in 0..3 {
            lf[c][i] = sum_xd(&nodes.iter().map(|p| p[c]).collect::<Vec<_>>(), &dl, i, n);
        }
    }
    let mut det_f = match dim {
        2 => lf[0][0] * lf[1][1] - lf[0][1] * lf[1][0],
        _ => {
            lf[0][0] * (lf[1][1] * lf[2][2] - lf[1][2] * lf[2][1])
                + lf[0][1] * (lf[1][2] * lf[2][0] - lf[1][0] * lf[2][2])
                + lf[0][2] * (lf[1][0] * lf[2][1] - lf[1][1] * lf[2][0])
        }
    };

    let mut lg = [[0.0f64; 3]; 3];
    if cdim == dim {
        if dim == 2 {
            lg[0][0] = 1.0 / det_f * lf[1][1];
            lg[0][1] = -1.0 / det_f * lf[0][1];
            lg[1][0] = -1.0 / det_f * lf[1][0];
            lg[1][1] = 1.0 / det_f * lf[0][0];
        } else {
            lg = invert_matrix_3x3(&lf, det_f);
        }
        // LG(1:dim,1:dim) = TRANSPOSE( LG(1:dim,1:dim) )
        let t = lg;
        for i in 0..dim {
            for j in 0..dim {
                lg[i][j] = t[j][i];
            }
        }
        for j in 0..dofs {
            let mut b = [0.0f64; 3];
            for k in 0..dim {
                let mut s = 0.0;
                for c in 0..dim {
                    s += lg[k][c] * ev.edge[j][c];
                }
                b[k] = s;
            }
            ev.edge[j][..dim].copy_from_slice(&b[..dim]);
            if dim == 2 {
                ev.curl[j][2] = 1.0 / det_f * ev.curl[j][2];
            } else {
                for k in 0..dim {
                    let mut s = 0.0;
                    for c in 0..dim {
                        s += lf[k][c] * ev.curl[j][c];
                    }
                    b[k] = 1.0 / det_f * s;
                }
                ev.curl[j][..dim].copy_from_slice(&b[..dim]);
            }
        }
        det_f = det_f.abs();
        ev.det = det_f;
    } else {
        // 2-D element embedded in 3-D: LG is the transpose of the
        // pseudoinverse of F (ElementMetric's LtoGMap).
        let (_m, det_j, ltog) = element_metric(dim, cdim, nodes, n, &dl)
            .ok_or_else(|| "degenerate element".to_string())?;
        lg = ltog;
        for j in 0..dofs {
            let mut b = [0.0f64; 3];
            for k in 0..cdim {
                let mut s = 0.0;
                for c in 0..dim {
                    s += lg[k][c] * ev.edge[j][c];
                }
                b[k] = s;
            }
            ev.edge[j][..cdim].copy_from_slice(&b[..cdim]);
            ev.curl[j][2] = 1.0 / det_j * ev.curl[j][2];
        }
        ev.det = det_j;
    }
    Ok(ev)
}

// ---------------------------------------------------------------------------
// ElementInfo (linear nodal part) + GetEdgeBasis
// ---------------------------------------------------------------------------

/// `NodalBasisFunctions` / `NodalFirstDerivatives` of the standard linear
/// Elmer elements (non-p path).
fn std_nodal(code: u16, u: f64, v: f64, w: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    match code {
        303 => (
            vec![1.0 - u - v, u, v],
            vec![[-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ),
        404 => {
            let c = 1.0 / 4.0;
            let us = [-1.0, 1.0, 1.0, -1.0];
            let vs = [-1.0, -1.0, 1.0, 1.0];
            let mut b = vec![0.0; 4];
            let mut d = vec![[0.0; 3]; 4];
            for q in 0..4 {
                b[q] = c * (1.0 + us[q] * u) * (1.0 + vs[q] * v);
                d[q][0] = c * us[q] * (1.0 + vs[q] * v);
                d[q][1] = c * (1.0 + us[q] * u) * vs[q];
            }
            (b, d)
        }
        504 => (
            vec![1.0 - u - v - w, u, v, w],
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ),
        605 => {
            // InterpolateInElement3D / FirstDerivativeIn[UVW]3D, element 605.
            let s = 1.0 / (1.0 - w);
            let b = vec![
                ((1.0 - u) * (1.0 - v) - w + u * v * w * s) / 4.0,
                ((1.0 + u) * (1.0 - v) - w - u * v * w * s) / 4.0,
                ((1.0 + u) * (1.0 + v) - w + u * v * w * s) / 4.0,
                ((1.0 - u) * (1.0 + v) - w - u * v * w * s) / 4.0,
                w,
            ];
            let d = vec![
                [
                    (-(1.0 - v) + v * w * s) / 4.0,
                    (-(1.0 - u) + u * w * s) / 4.0,
                    (-1.0 + u * v * sq(s)) / 4.0,
                ],
                [
                    ((1.0 - v) - v * w * s) / 4.0,
                    (-(1.0 + u) - u * w * s) / 4.0,
                    (-1.0 - u * v * sq(s)) / 4.0,
                ],
                [
                    ((1.0 + v) + v * w * s) / 4.0,
                    ((1.0 + u) + u * w * s) / 4.0,
                    (-1.0 + u * v * sq(s)) / 4.0,
                ],
                [
                    (-(1.0 + v) - v * w * s) / 4.0,
                    ((1.0 - u) - u * w * s) / 4.0,
                    (-1.0 - u * v * sq(s)) / 4.0,
                ],
                [0.0, 0.0, 1.0],
            ];
            (b, d)
        }
        706 => {
            let half = 1.0 / 2.0;
            let tri = [1.0 - u - v, u, v];
            let line = [half * (1.0 - w), half * (1.0 + w)];
            let gradtri = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            let gradline = [-half, half];
            let mut b = vec![0.0; 6];
            let mut d = vec![[0.0; 3]; 6];
            for q in 0..3 {
                b[q] = line[0] * tri[q];
                b[q + 3] = line[1] * tri[q];
                d[q][0] = gradtri[q][0] * line[0];
                d[q][1] = gradtri[q][1] * line[0];
                d[q + 3][0] = gradtri[q][0] * line[1];
                d[q + 3][1] = gradtri[q][1] * line[1];
                d[q][2] = tri[q] * gradline[0];
                d[q + 3][2] = tri[q] * gradline[1];
            }
            (b, d)
        }
        808 => {
            let c = 1.0 / 8.0;
            let us = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
            let vs = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
            let ws = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
            let mut b = vec![0.0; 8];
            let mut d = vec![[0.0; 3]; 8];
            for q in 0..8 {
                b[q] = c * (1.0 + us[q] * u) * (1.0 + vs[q] * v) * (1.0 + ws[q] * w);
                d[q][0] = c * us[q] * (1.0 + vs[q] * v) * (1.0 + ws[q] * w);
                d[q][1] = c * (1.0 + us[q] * u) * vs[q] * (1.0 + ws[q] * w);
                d[q][2] = c * (1.0 + us[q] * u) * (1.0 + vs[q] * v) * ws[q];
            }
            (b, d)
        }
        _ => unreachable!(),
    }
}

/// `ElementInfo(..., EdgeBasis=, RotBasis=)` without the Piola transform:
/// the nodal basis and global derivatives are computed first and
/// `GetEdgeBasis` then builds the lowest-order Whitney-type basis on the
/// physical element. Orientation follows `gindex`. Only the linear element
/// types 303, 404, 504, 605, 706 and 808 are handled.
///
/// `det` of the result is `DetJ` from `ElementMetric`.
pub fn get_edge_basis(
    etype: ElementType,
    nodes: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
    cdim: usize,
) -> Result<EdgeEval, String> {
    let code = etype.0;
    if !supported(code) {
        return Err(format!("edge element type {code} is not supported"));
    }
    let n = (code % 100) as usize;
    if nodes.len() < n || gindex.len() < n {
        return Err("fewer nodal coordinates / indices than element nodes".into());
    }
    let dim = element_dim(code);
    let family = code / 100;
    let [u0, v0, w0] = uvw;
    let (basis_v, dl) = std_nodal(code, u0, v0, w0);
    let (_m, det_j, ltog) =
        element_metric(dim, cdim, nodes, n, &dl).ok_or_else(|| "degenerate element".to_string())?;
    // dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k) * LtoGMap(j,k)
    let mut dbasisdx = [[0.0f64; 3]; 64];
    for k in 0..dim {
        for j in 0..cdim {
            for i in 0..n {
                dbasisdx[i][j] += dl[i][k] * ltog[j][k];
            }
        }
    }
    let mut basis = [0.0f64; 64];
    basis[..n].copy_from_slice(&basis_v);

    // ---- GetEdgeBasis ----
    let (nu, nv, nw) = std_node_uvw(code);
    let (mut u, mut v, mut w) = (0.0f64, 0.0f64, 0.0f64);
    let mut dudx = [[0.0f64; 3]; 3];
    let mut tribase = [0.0f64; 3];
    let mut dtribase = [[0.0f64; 3]; 3];
    let mut g = [[0.0f64; 3]; 3];
    let mut f = [[0.0f64; 3]; 3];
    let mut det_f = 0.0f64;
    let mut edgebasis_l = [[0.0f64; 3]; 8];
    let mut curlbasis_l = [[0.0f64; 3]; 8];
    if matches!(family, 4 | 7 | 8 | 6) {
        let (mut su, mut sv, mut sw) = (0.0, 0.0, 0.0);
        for q in 0..n {
            su += basis[q] * nu[q];
        }
        for q in 0..n {
            sv += basis[q] * nv[q];
        }
        for q in 0..n {
            sw += basis[q] * nw[q];
        }
        u = su;
        v = sv;
        w = sw;
        for c in 0..3 {
            let mut a = 0.0;
            let mut b = 0.0;
            let mut d = 0.0;
            for q in 0..n {
                a += nu[q] * dbasisdx[q][c];
                b += nv[q] * dbasisdx[q][c];
                d += nw[q] * dbasisdx[q][c];
            }
            dudx[0][c] = a;
            dudx[1][c] = b;
            dudx[2][c] = d;
        }
    }
    if family != 6 {
        tribase[0] = 1.0 - u - v;
        tribase[1] = u;
        tribase[2] = v;
        for c in 0..3 {
            dtribase[0][c] = -dudx[0][c] - dudx[1][c];
            dtribase[1][c] = dudx[0][c];
            dtribase[2][c] = dudx[1][c];
        }
    }
    if family == 6 {
        g = dudx;
        let det_g = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
            + g[0][1] * (g[1][2] * g[2][0] - g[1][0] * g[2][2])
            + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
        det_f = 1.0 / det_g;
        f = invert_matrix_3x3(&g, det_g);
        let edgebasis = &mut edgebasis_l;
        let curlbasis = &mut curlbasis_l;
            edgebasis[0][0] = (1.0e0 - v - w) / 4.0e0;
            edgebasis[0][1] = 0.0e0;
            edgebasis[0][2] = u * (-1.0e0 + v + w) / (4.0e0 * (-1.0e0 + w));
            curlbasis[0][0] = u / (4.0e0 * (-1.0e0 + w));
            curlbasis[0][1] = -((-2.0e0 + v + 2.0e0 * w) / (4.0e0 * (-1.0e0 + w)));
            curlbasis[0][2] = 0.25e0;
            edgebasis[1][0] = 0.0e0;
            edgebasis[1][1] = (1.0e0 + u - w) / 4.0e0;
            edgebasis[1][2] = v * (1.0e0 + u - w) / (4.0e0 - 4.0e0 * w);
            curlbasis[1][0] = (2.0e0 + u - 2.0e0 * w) / (4.0e0 - 4.0e0 * w);
            curlbasis[1][1] = v / (4.0e0 * (-1.0e0 + w));
            curlbasis[1][2] = 0.25e0;
            edgebasis[2][0] = (1.0e0 + v - w) / 4.0e0;
            edgebasis[2][1] = 0.0e0;
            edgebasis[2][2] = u * (1.0e0 + v - w) / (4.0e0 - 4.0e0 * w);
            curlbasis[2][0] = u / (4.0e0 - 4.0e0 * w);
            curlbasis[2][1] = (2.0e0 + v - 2.0e0 * w) / (4.0e0 * (-1.0e0 + w));
            curlbasis[2][2] = -0.25e0;
            edgebasis[3][0] = 0.0e0;
            edgebasis[3][1] = (1.0e0 - u - w) / 4.0e0;
            edgebasis[3][2] = v * (-1.0e0 + u + w) / (4.0e0 * (-1.0e0 + w));
            curlbasis[3][0] = (-2.0e0 + u + 2.0e0 * w) / (4.0e0 * (-1.0e0 + w));
            curlbasis[3][1] = v / (4.0e0 - 4.0e0 * w);
            curlbasis[3][2] = -0.25e0;
            edgebasis[4][0] = w * (-1.0e0 + v + w) / (4.0e0 * (-1.0e0 + w));
            edgebasis[4][1] = w * (-1.0e0 + u + w) / (4.0e0 * (-1.0e0 + w));
            edgebasis[4][2] = (-((-1.0e0 + v) * sq(-1.0e0 + w)) + u * (v - sq(-1.0e0 + w) - 2.0e0 * v * w)) / (4.0e0 * sq(-1.0e0 + w));
            curlbasis[4][0] = -((-1.0e0 + u + w) / (2.0e0 * (-1.0e0 + w)));
            curlbasis[4][1] = (-1.0e0 + v + w) / (2.0e0 * (-1.0e0 + w));
            curlbasis[4][2] = 0.0e0;
            edgebasis[5][0] = -(w * (-1.0e0 + v + w) / (4.0e0 * (-1.0e0 + w)));
            edgebasis[5][1] = w * (-1.0e0 - u + w) / (4.0e0 * (-1.0e0 + w));
            edgebasis[5][2] = (-((-1.0e0 + v) * sq(-1.0e0 + w)) + u * (sq(-1.0e0 + w) + v * (-1.0e0 + 2.0e0 * w))) / (4.0e0 * sq(-1.0e0 + w));
            curlbasis[5][0] = (1.0e0 + u - w) / (2.0e0 * (-1.0e0 + w));
            curlbasis[5][1] = -((-1.0e0 + v + w) / (2.0e0 * (-1.0e0 + w)));
            curlbasis[5][2] = 0.0e0;
            edgebasis[6][0] = (1.0e0 + v - w) * w / (4.0e0 * (-1.0e0 + w));
            edgebasis[6][1] = (1.0e0 + u - w) * w / (4.0e0 * (-1.0e0 + w));
            edgebasis[6][2] = ((1.0e0 + v) * sq(-1.0e0 + w) + u * (v + sq(-1.0e0 + w) - 2.0e0 * v * w)) / (4.0e0 * sq(-1.0e0 + w));
            curlbasis[6][0] = (1.0e0 + u - w) / (2.0e0 - 2.0e0 * w);
            curlbasis[6][1] = (1.0e0 + v - w) / (2.0e0 * (-1.0e0 + w));
            curlbasis[6][2] = 0.0e0;
            edgebasis[7][0] = w * (-1.0e0 - v + w) / (4.0e0 * (-1.0e0 + w));
            edgebasis[7][1] = -(w * (-1.0e0 + u + w) / (4.0e0 * (-1.0e0 + w)));
            edgebasis[7][2] = ((1.0e0 + v) * sq(-1.0e0 + w) - u * (v + sq(-1.0e0 + w) - 2.0e0 * v * w)) / (4.0e0 * sq(-1.0e0 + w));
            curlbasis[7][0] = (-1.0e0 + u + w) / (2.0e0 * (-1.0e0 + w));
            curlbasis[7][1] = (1.0e0 + v - w) / (2.0e0 - 2.0e0 * w);
            curlbasis[7][2] = 0.0e0;
    }

    let edgemap = get_edge_map(family);
    let nedges = edgemap.len();
    let mut wbasis = [[0.0f64; 3]; 64];
    let mut rotwbasis = [[0.0f64; 3]; 64];
    let (mut i, mut j, mut k);
    let mut base = 0.0f64;
    let mut dbase = [0.0f64; 3];
    let mut du = [0.0f64; 3];
    let mut tbase = [0.0f64; 3];
    let mut rbase = [0.0f64; 3];
    for ii in 1..=nedges {
        i = ii;
        j = edgemap[i - 1][0];
        k = edgemap[i - 1][1];
        let nj = gindex[j - 1];
        let nk = gindex[k - 1];
        match family {
            3 | 5 => {
                for c in 0..3 {
                    wbasis[i - 1][c] = basis[j - 1] * dbasisdx[k - 1][c] - basis[k - 1] * dbasisdx[j - 1][c];
                }
                rotwbasis[i - 1][0] = 2.0 * (dbasisdx[j - 1][1] * dbasisdx[k - 1][2] - dbasisdx[j - 1][2] * dbasisdx[k - 1][1]);
                rotwbasis[i - 1][1] = 2.0 * (dbasisdx[j - 1][2] * dbasisdx[k - 1][0] - dbasisdx[j - 1][0] * dbasisdx[k - 1][2]);
                rotwbasis[i - 1][2] = 2.0 * (dbasisdx[j - 1][0] * dbasisdx[k - 1][1] - dbasisdx[j - 1][1] * dbasisdx[k - 1][0]);
            }
            6 => {
                // Piola transform of the reference-element functions
                for kk in 0..3 {
                    let mut s = 0.0;
                    for c in 0..3 {
                        s += g[c][kk] * edgebasis_l[i - 1][c];
                    }
                    wbasis[i - 1][kk] = s;
                }
                for kk in 0..3 {
                    let mut s = 0.0;
                    for c in 0..3 {
                        s += f[kk][c] * curlbasis_l[i - 1][c];
                    }
                    rotwbasis[i - 1][kk] = 1.0 / det_f * s;
                }
            }
            7 => {
                match i {
                    1 => {
                        j = 1;
                        k = 2;
                        base = (1.0 - w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[2][s0] / 2.0);
                        }
                    }
                    2 => {
                        j = 2;
                        k = 3;
                        base = (1.0 - w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[2][s0] / 2.0);
                        }
                    }
                    3 => {
                        j = 3;
                        k = 1;
                        base = (1.0 - w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[2][s0] / 2.0);
                        }
                    }
                    4 => {
                        j = 1;
                        k = 2;
                        base = (1.0 + w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    5 => {
                        j = 2;
                        k = 3;
                        base = (1.0 + w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    6 => {
                        j = 3;
                        k = 1;
                        base = (1.0 + w) / 2.0;
                        for s0 in 0..3 {
                            dbase[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    7 => {
                        base = tribase[0];
                        for s0 in 0..3 {
                            dbase[s0] = dtribase[0][s0];
                        }
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    8 => {
                        base = tribase[1];
                        for s0 in 0..3 {
                            dbase[s0] = dtribase[1][s0];
                        }
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    9 => {
                        base = tribase[2];
                        for s0 in 0..3 {
                            dbase[s0] = dtribase[2][s0];
                        }
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0] / 2.0;
                        }
                    }
                    _ => {}
                }
                if i <= 6 {
                    for s0 in 0..3 {
                        tbase[s0] = tribase[j - 1] * dtribase[k - 1][s0] - tribase[k - 1] * dtribase[j - 1][s0];
                    }
                    rbase[0] = 2.0 * base * (dtribase[j - 1][1] * dtribase[k - 1][2] - dtribase[k - 1][1] * dtribase[j - 1][2]) + dbase[1] * tbase[2] - dbase[2] * tbase[1];
                    rbase[1] = 2.0 * base * (dtribase[j - 1][2] * dtribase[k - 1][0] - dtribase[k - 1][2] * dtribase[j - 1][0]) + dbase[2] * tbase[0] - dbase[0] * tbase[2];
                    rbase[2] = 2.0 * base * (dtribase[j - 1][0] * dtribase[k - 1][1] - dtribase[k - 1][0] * dtribase[j - 1][1]) + dbase[0] * tbase[1] - dbase[1] * tbase[0];
                    for s0 in 0..3 {
                        rotwbasis[i - 1][s0] = rbase[s0];
                    }
                    for s0 in 0..3 {
                        wbasis[i - 1][s0] = tbase[s0] * base;
                    }
                } else {
                    for s0 in 0..3 {
                        wbasis[i - 1][s0] = base * du[s0];
                    }
                    rotwbasis[i - 1][0] = dbase[1] * du[2] - dbase[2] * du[1];
                    rotwbasis[i - 1][1] = dbase[2] * du[0] - dbase[0] * du[2];
                    rotwbasis[i - 1][2] = dbase[0] * du[1] - dbase[1] * du[0];
                }
            }
            4 => {
                match i {
                    1 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[0][s0];
                        }
                        base = (1.0 - v) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[1][s0] * (1.0 - w)) - (1.0 - v) * dudx[2][s0];
                        }
                    }
                    2 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[1][s0];
                        }
                        base = (1.0 + u) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[0][s0] * (1.0 - w) - (1.0 + u) * dudx[2][s0];
                        }
                    }
                    3 => {
                        for s0 in 0..3 {
                            du[s0] = -dudx[0][s0];
                        }
                        base = (1.0 + v) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[1][s0] * (1.0 - w) - (1.0 + v) * dudx[2][s0];
                        }
                    }
                    4 => {
                        for s0 in 0..3 {
                            du[s0] = -dudx[1][s0];
                        }
                        base = (1.0 - u) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[0][s0] * (1.0 - w)) - (1.0 - u) * dudx[2][s0];
                        }
                    }
                    _ => {}
                }
                for s0 in 0..3 {
                    wbasis[i - 1][s0] = base * du[s0] / (n as f64);
                }
                rotwbasis[i - 1][0] = (dbase[1] * du[2] - dbase[2] * du[1]) / (n as f64);
                rotwbasis[i - 1][1] = (dbase[2] * du[0] - dbase[0] * du[2]) / (n as f64);
                rotwbasis[i - 1][2] = (dbase[0] * du[1] - dbase[1] * du[0]) / (n as f64);
                for c in 0..3 {
                    wbasis[i - 1][c] = base * du[c] / (n as f64);
                }
                rotwbasis[i - 1][0] = (dbase[1] * du[2] - dbase[2] * du[1]) / (n as f64);
                rotwbasis[i - 1][1] = (dbase[2] * du[0] - dbase[0] * du[2]) / (n as f64);
                rotwbasis[i - 1][2] = (dbase[0] * du[1] - dbase[1] * du[0]) / (n as f64);
            }
            8 => {
                match i {
                    1 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[0][s0];
                        }
                        base = (1.0 - v) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[1][s0] * (1.0 - w)) - (1.0 - v) * dudx[2][s0];
                        }
                    }
                    2 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[1][s0];
                        }
                        base = (1.0 + u) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[0][s0] * (1.0 - w) - (1.0 + u) * dudx[2][s0];
                        }
                    }
                    3 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[0][s0];
                        }
                        base = (1.0 + v) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[1][s0] * (1.0 - w) - (1.0 + v) * dudx[2][s0];
                        }
                    }
                    4 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[1][s0];
                        }
                        base = (1.0 - u) * (1.0 - w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[0][s0] * (1.0 - w)) - (1.0 - u) * dudx[2][s0];
                        }
                    }
                    5 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[0][s0];
                        }
                        base = (1.0 - v) * (1.0 + w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[1][s0] * (1.0 + w)) + (1.0 - v) * dudx[2][s0];
                        }
                    }
                    6 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[1][s0];
                        }
                        base = (1.0 + u) * (1.0 + w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[0][s0] * (1.0 + w) + (1.0 + u) * dudx[2][s0];
                        }
                    }
                    7 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[0][s0];
                        }
                        base = (1.0 + v) * (1.0 + w);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[1][s0] * (1.0 + w) + (1.0 + v) * dudx[2][s0];
                        }
                    }
                    8 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[1][s0];
                        }
                        base = (1.0 - u) * (1.0 + w);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[0][s0] * (1.0 + w)) + (1.0 - u) * dudx[2][s0];
                        }
                    }
                    9 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0];
                        }
                        base = (1.0 - u) * (1.0 - v);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[0][s0] * (1.0 - v)) - (1.0 - u) * dudx[1][s0];
                        }
                    }
                    10 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0];
                        }
                        base = (1.0 + u) * (1.0 - v);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[0][s0] * (1.0 - v) - (1.0 + u) * dudx[1][s0];
                        }
                    }
                    11 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0];
                        }
                        base = (1.0 + u) * (1.0 + v);
                        for s0 in 0..3 {
                            dbase[s0] = dudx[0][s0] * (1.0 + v) + (1.0 + u) * dudx[1][s0];
                        }
                    }
                    12 => {
                        for s0 in 0..3 {
                            du[s0] = dudx[2][s0];
                        }
                        base = (1.0 - u) * (1.0 + v);
                        for s0 in 0..3 {
                            dbase[s0] = -(dudx[0][s0] * (1.0 + v)) + (1.0 - u) * dudx[1][s0];
                        }
                    }
                    _ => {}
                }
                for s0 in 0..3 {
                    wbasis[i - 1][s0] = base * du[s0] / (n as f64);
                }
                rotwbasis[i - 1][0] = (dbase[1] * du[2] - dbase[2] * du[1]) / (n as f64);
                rotwbasis[i - 1][1] = (dbase[2] * du[0] - dbase[0] * du[2]) / (n as f64);
                rotwbasis[i - 1][2] = (dbase[0] * du[1] - dbase[1] * du[0]) / (n as f64);
                for c in 0..3 {
                    wbasis[i - 1][c] = base * du[c] / (n as f64);
                }
                rotwbasis[i - 1][0] = (dbase[1] * du[2] - dbase[2] * du[1]) / (n as f64);
                rotwbasis[i - 1][1] = (dbase[2] * du[0] - dbase[0] * du[2]) / (n as f64);
                rotwbasis[i - 1][2] = (dbase[0] * du[1] - dbase[1] * du[0]) / (n as f64);
            }
            _ => return Err("Edge Basis: not implemented for this element type".into()),
        }
        if nk < nj {
            for c in 0..3 {
                wbasis[i - 1][c] = -wbasis[i - 1][c];
                rotwbasis[i - 1][c] = -rotwbasis[i - 1][c];
            }
        }
    }
    Ok(EdgeEval {
        dofs: nedges,
        det: det_j,
        edge: wbasis[..nedges].to_vec(),
        curl: rotwbasis[..nedges].to_vec(),
    })
}

/// Edge basis on the reference element (`gindex` = local node order, i.e.
/// every edge points from its first to its second node as listed in
/// `GetEdgeMap`, flipped when the second local number is smaller).
///
/// * `piola == true`: the Piola-version reference functions `b` and
///   `Curl b` of `EdgeElementInfo` for `order` 1 or 2 (also the face and
///   interior functions), evaluated on Elmer's p-reference element.
/// * `piola == false`: the older lowest-order basis of `GetEdgeBasis`
///   evaluated on the standard Elmer reference element (identity element
///   map); `order` must be 1.
///
/// Writes up to `out.len()` / `out_curl.len()` functions and returns the
/// number of functions the element has; an error is returned for an
/// unsupported combination or too small output slices.
pub fn edge_basis(
    etype: ElementType,
    order: u8,
    piola: bool,
    uvw: [f64; 3],
    out: &mut [[f64; 3]],
    out_curl: &mut [[f64; 3]],
) -> Result<usize, String> {
    let code = etype.0;
    if !supported(code) {
        return Err(format!("edge element type {code} is not supported"));
    }
    let n = (code % 100) as usize;
    let gindex: Vec<i64> = (1..=n as i64).collect();
    let ev = if piola {
        edge_basis_reference(etype, order, &gindex, uvw)?
    } else {
        if order != 1 {
            return Err("the non-Piola GetEdgeBasis design is lowest-order only".into());
        }
        let (nu, nv, nw) = std_node_uvw(code);
        let nodes: Vec<[f64; 3]> = (0..n).map(|q| [nu[q], nv[q], nw[q]]).collect();
        get_edge_basis(etype, &nodes, &gindex, uvw, element_dim(code))?
    };
    if out.len() < ev.dofs || out_curl.len() < ev.dofs {
        return Err(format!("output slices must hold {} functions", ev.dofs));
    }
    out[..ev.dofs].copy_from_slice(&ev.edge);
    out_curl[..ev.dofs].copy_from_slice(&ev.curl);
    Ok(ev.dofs)
}
