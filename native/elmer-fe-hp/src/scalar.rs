//! The scalar `ElementInfo` path for p elements (ElemInfo.F90:331), as far
//! as it differs from the vectorised `ElementInfoVec` path that the rest of
//! this crate ports.
//!
//! A solver that calls `ElementInfo` gets its nodal functions from
//! `NodalBasisFunctions` / `NodalFirstDerivatives` (the `*NodalPBasisAll`
//! routines, ElementBasis.F90:2148, :2253), its edge functions from
//! PElementBase.F90's scalar `TriangleEdgePBasis` / `dTriangleEdgePBasis`
//! (each edge hand-expanded), and its geometry from `ElementMetric`
//! (ElemInfo.F90:9893), which divides by the metric determinant where
//! `ElementMetricVec` multiplies by its reciprocal. The values agree with
//! H1Basis.F90's to rounding. Matching them bit for bit matters when the
//! assembled matrix feeds an incomplete factorisation: a round-off entry
//! that one path leaves nonzero and the other cancels to zero changes
//! ILU(0)'s sparsity pattern (StokesPFEM).
//!
//! Ported so far, for degrees up to 5 (the `varPhi` / `dVarPhi` / `dPhi`
//! cases whose powers are `x**2` and `x**4`, which every gfortran
//! expansion of `**` evaluates the same way):
//! * triangles (303): nodal and edge functions;
//! * quadrilaterals (404, full scheme): nodal, edge and interior bubble
//!   functions (`QuadEdgePBasis`, `QuadBubblePBasis` without
//!   `localNumbers`, PElementBase.F90:915, :981, :1168, :1219).
//!
//! Triangle bubbles, face-directed quad bubbles (a 2D element on a 3D
//! face), serendipity quads and higher degrees keep the H1Basis values.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

/// `TriangleNodalPBasis(node, u, v)` (PElementBase.F90:1493), the form the
/// scalar edge functions use.
fn l(node: usize, u: f64, v: f64) -> f64 {
    let s3 = 3f64.sqrt();
    match node {
        1 => (1.0 - u - v / s3) / 2.0,
        2 => (1.0 + u - v / s3) / 2.0,
        _ => v / s3,
    }
}

/// `varPhi(i, x)` (PElementBase.F90:8600), cases 2-5.
fn var_phi(i: usize, x: f64) -> Option<f64> {
    Some(match i {
        2 => -(6f64.sqrt()),
        3 => -(x * 10f64.sqrt()),
        4 => -((5.0 * (x * x) - 1.0) * 14f64.sqrt() / 4.0),
        5 => -(0.75 * (7.0 * (x * x) - 3.0) * x * 2f64.sqrt()),
        _ => return None,
    })
}

/// `dVarPhi(i, x)` (PElementBase.F90:8709), cases 2-5.
fn d_var_phi(i: usize, x: f64) -> Option<f64> {
    Some(match i {
        2 => 0.0,
        3 => -(10f64.sqrt()),
        4 => -(2.5 * x * 14f64.sqrt()),
        5 => -(15.75 * (x * x) * 2f64.sqrt()) + 2.25 * 2f64.sqrt(),
        _ => return None,
    })
}

/// `TriangleNodalPBasisAll` and `dTriangleNodalPBasisAll`
/// (PElementBase.F90:1531, :1629) into `val[0..3]`, `grad[0..3]`.
pub(crate) fn triangle_nodal(u: f64, v: f64, val: &mut [f64], grad: &mut [[f64; 3]]) {
    let half = 1.0 / 2.0;
    let c3 = 1.0 / 3f64.sqrt();
    let c6 = 3f64.sqrt() / 6.0;
    val[0] = half * (1.0 - u - c3 * v);
    val[1] = half * (1.0 + u - c3 * v);
    val[2] = c3 * v;
    grad[0] = [-half, -c6, 0.0];
    grad[1] = [half, -c6, 0.0];
    grad[2] = [0.0, 2.0 * c6, 0.0];
}

/// `TriangleEdgePBasis(edge, i, u, v, invert)` and `dTriangleEdgePBasis`
/// (PElementBase.F90:1660, :1740) for `edge` in 1..=3 and `i` = degree.
/// `None` beyond the ported degrees.
pub(crate) fn triangle_edge(
    edge: usize,
    i: usize,
    u: f64,
    v: f64,
    invert: bool,
) -> Option<(f64, [f64; 2])> {
    let s3 = 3f64.sqrt();
    let (a, b, c) = (s3 / 6.0, s3 / 3.0, s3 / 2.0);
    Some(match edge {
        1 => {
            let (l1, l2) = (l(1, u, v), l(2, u, v));
            let value = l1 * l2 * var_phi(i, if invert { l1 - l2 } else { l2 - l1 })?;
            let x = if invert { -u } else { u };
            let vphi = var_phi(i, x)?;
            let dv = d_var_phi(i, x)?;
            let g0 = if invert {
                -(0.5 * l2 * vphi) + 0.5 * l1 * vphi - l1 * l2 * dv
            } else {
                -(0.5 * l2 * vphi) + 0.5 * l1 * vphi + l1 * l2 * dv
            };
            let g1 = -(a * l2 * vphi) - a * l1 * vphi;
            (value, [g0, g1])
        }
        2 => {
            let (l2, l3) = (l(2, u, v), l(3, u, v));
            let x = if invert { l2 - l3 } else { l3 - l2 };
            let vphi = var_phi(i, x)?;
            let dv = d_var_phi(i, x)?;
            let value = l2 * l3 * vphi;
            if invert {
                (
                    value,
                    [
                        0.5 * l3 * vphi + 0.5 * l2 * l3 * dv,
                        -(a * l3 * vphi) + b * l2 * vphi - c * l2 * l3 * dv,
                    ],
                )
            } else {
                (
                    value,
                    [
                        0.5 * l3 * vphi - 0.5 * l2 * l3 * dv,
                        -(a * l3 * vphi) + b * l2 * vphi + c * l2 * l3 * dv,
                    ],
                )
            }
        }
        _ => {
            let (l1, l3) = (l(1, u, v), l(3, u, v));
            let x = if invert { l3 - l1 } else { l1 - l3 };
            let vphi = var_phi(i, x)?;
            let dv = d_var_phi(i, x)?;
            let value = l1 * l3 * vphi;
            if invert {
                (
                    value,
                    [
                        -(0.5 * l3 * vphi) + 0.5 * l1 * l3 * dv,
                        -(a * l3 * vphi) + b * l1 * vphi + c * l1 * l3 * dv,
                    ],
                )
            } else {
                (
                    value,
                    [
                        -(0.5 * l3 * vphi) - 0.5 * l1 * l3 * dv,
                        -(a * l3 * vphi) + b * l1 * vphi - c * l1 * l3 * dv,
                    ],
                )
            }
        }
    })
}

/// `QuadNodalPBasis(node, u, v)` and `dQuadNodalPBasis` (PElementBase.F90:327,
/// :384), the forms the quad edge functions use.
fn quad_n(node: usize, u: f64, v: f64) -> (f64, [f64; 2]) {
    match node {
        1 => (
            (1.0 - u) * (1.0 - v) / 4.0,
            [-((1.0 - v) / 4.0), -((1.0 - u) / 4.0)],
        ),
        2 => (
            (1.0 + u) * (1.0 - v) / 4.0,
            [(1.0 - v) / 4.0, -((1.0 + u) / 4.0)],
        ),
        3 => (
            (1.0 + u) * (1.0 + v) / 4.0,
            [(1.0 + v) / 4.0, (1.0 + u) / 4.0],
        ),
        _ => (
            (1.0 - u) * (1.0 + v) / 4.0,
            [-((1.0 + v) / 4.0), (1.0 - u) / 4.0],
        ),
    }
}

/// `QuadL(which, u, v)` and `dQuadL`.
fn quad_l(node: usize, u: f64, v: f64) -> (f64, [f64; 2]) {
    match node {
        1 => ((2.0 - u - v) / 2.0, [-0.5, -0.5]),
        2 => ((2.0 + u - v) / 2.0, [0.5, -0.5]),
        3 => ((2.0 + u + v) / 2.0, [0.5, 0.5]),
        _ => ((2.0 - u + v) / 2.0, [-0.5, 0.5]),
    }
}

/// `QuadEdgePBasis(edge, i, u, v, invert)` and `dQuadEdgePBasis`, edge in
/// 1..=4 on `getQuadEdgeMap` ((1,2), (2,3), (4,3), (1,4)).
pub(crate) fn quad_edge(
    edge: usize,
    i: usize,
    u: f64,
    v: f64,
    invert: bool,
) -> Option<(f64, [f64; 2])> {
    let [a, b] = [[1, 2], [2, 3], [4, 3], [1, 4]][edge - 1];
    let ((na, dna), (nb, dnb)) = (quad_n(a, u, v), quad_n(b, u, v));
    let ((la, dla), (lb, dlb)) = (quad_l(a, u, v), quad_l(b, u, v));
    let (mut par, mut dpar) = (lb - la, [dlb[0] - dla[0], dlb[1] - dla[1]]);
    if invert {
        par = -par;
        dpar = [-dpar[0], -dpar[1]];
    }
    let vphi = var_phi(i, par)?;
    let dv = d_var_phi(i, par)?;
    let dvphi = [dv * dpar[0], dv * dpar[1]];
    let g = |k: usize| dna[k] * nb * vphi + na * dnb[k] * vphi + na * nb * dvphi[k];
    Some((na * nb * vphi, [g(0), g(1)]))
}

/// `Phi(i, x) = varPhi(i, x) (1 - x**2) / 4` (PElementBase.F90:8329) and
/// `dPhi` (:8380), cases 2-5.
fn phi(i: usize, x: f64) -> Option<f64> {
    Some(var_phi(i, x)? * (1.0 - x * x) / 4.0)
}

fn d_phi(i: usize, x: f64) -> Option<f64> {
    let x2 = x * x;
    Some(match i {
        2 => 6f64.sqrt() * x / 2.0,
        3 => 10f64.sqrt() * (3.0 * x2 - 1.0) / 4.0,
        4 => 14f64.sqrt() * x * (5.0 * x2 - 3.0) / 4.0,
        5 => 3.0 / 16.0 * 2f64.sqrt() * (35.0 * (x2 * x2) - 30.0 * x2 + 3.0),
        _ => return None,
    })
}

/// `QuadBubblePBasis(i, j, u, v)` and `dQuadBubblePBasis` without
/// `localNumbers`: `Phi(i+2, u) Phi(j+2, v)`.
pub(crate) fn quad_bubble(i: usize, j: usize, u: f64, v: f64) -> Option<(f64, [f64; 2])> {
    let (pu, pv) = (phi(i + 2, u)?, phi(j + 2, v)?);
    Some((pu * pv, [d_phi(i + 2, u)? * pv, pu * d_phi(j + 2, v)?]))
}

/// `ElementMetric` (ElemInfo.F90:9893) for a 1D or 2D element in a `cdim`
/// coordinate system: `dx` from the first `n` functions, then `detJ`
/// and `LtoGMap`. `Err` for a degenerate element.
pub(crate) fn element_metric(
    dim: usize,
    cdim: usize,
    xyz: &[[f64; 3]],
    dl: &[[f64; 3]],
) -> Result<(f64, [[f64; 3]; 3]), String> {
    let n = xyz.len().min(dl.len());
    let mut dx = [[0.0f64; 3]; 3];
    for i in 0..dim {
        for k in 0..3 {
            let mut s = 0.0;
            for m in 0..n {
                s += xyz[m][k] * dl[m][i];
            }
            dx[k][i] = s;
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
    let eps = f64::EPSILON.powi(dim as i32);
    let mut metric = [[0.0f64; 3]; 3];
    let detg = match dim {
        1 => {
            let d = g[0][0];
            if d <= eps {
                return Err("degenerate 1D element".into());
            }
            metric[0][0] = 1.0 / d;
            d.sqrt()
        }
        2 => {
            let d = g[0][0] * g[1][1] - g[0][1] * g[1][0];
            if d <= eps {
                return Err("degenerate 2D element".into());
            }
            metric[0][0] = g[1][1] / d;
            metric[0][1] = -g[0][1] / d;
            metric[1][0] = -g[1][0] / d;
            metric[1][1] = g[0][0] / d;
            d.sqrt()
        }
        _ => {
            return Err(format!(
                "scalar ElementMetric of a {dim}D element is not ported"
            ))
        }
    };
    let mut l2g = [[0.0f64; 3]; 3];
    for i in 0..cdim {
        for j in 0..dim {
            let mut s = 0.0;
            for k in 0..dim {
                s += dx[i][k] * metric[k][j];
            }
            l2g[i][j] = s;
        }
    }
    Ok((detg, l2g))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The hand-expanded scalar gradients are the derivatives of the values.
    #[test]
    fn edge_gradients_differentiate_the_values() {
        let (u, v, h) = (0.13, 0.41, 1e-6);
        for edge in 1..=3 {
            for i in 2..=5 {
                for invert in [false, true] {
                    let (_, g) = triangle_edge(edge, i, u, v, invert).unwrap();
                    let f = |u: f64, v: f64| triangle_edge(edge, i, u, v, invert).unwrap().0;
                    let du = (f(u + h, v) - f(u - h, v)) / (2.0 * h);
                    let dv = (f(u, v + h) - f(u, v - h)) / (2.0 * h);
                    assert!(
                        (g[0] - du).abs() < 1e-7 && (g[1] - dv).abs() < 1e-7,
                        "edge {edge} i {i} inv {invert}: {g:?} vs {du} {dv}"
                    );
                }
            }
        }
        assert!(triangle_edge(1, 6, u, v, false).is_none());
    }

    #[test]
    fn quad_gradients_differentiate_the_values() {
        let (u, v, h) = (0.13, -0.41, 1e-6);
        let check = |f: &dyn Fn(f64, f64) -> (f64, [f64; 2]), what: &str| {
            let (_, g) = f(u, v);
            let du = (f(u + h, v).0 - f(u - h, v).0) / (2.0 * h);
            let dv = (f(u, v + h).0 - f(u, v - h).0) / (2.0 * h);
            assert!(
                (g[0] - du).abs() < 1e-7 && (g[1] - dv).abs() < 1e-7,
                "{what}: {g:?} vs {du} {dv}"
            );
        };
        for edge in 1..=4 {
            for i in 2..=5 {
                for invert in [false, true] {
                    check(
                        &|u, v| quad_edge(edge, i, u, v, invert).unwrap(),
                        &format!("edge {edge} i {i} {invert}"),
                    );
                }
            }
        }
        for i in 0..=3 {
            for j in 0..=3 {
                check(
                    &|u, v| quad_bubble(i, j, u, v).unwrap(),
                    &format!("bubble {i} {j}"),
                );
            }
        }
        assert!(quad_bubble(4, 0, u, v).is_none());
    }
}
