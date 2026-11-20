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
//! Ported so far, for degrees up to 9 (`varPhi` / `dVarPhi` / `dPhi`
//! cases 2-9). Their integer powers `x**n` are evaluated the way the
//! release build (gfortran `-O3`) does: GCC's `powi_table` addition chain
//! ([`pow`]). The debug build calls libgcc's square-and-multiply
//! `__powidf2` instead, which rounds differently from n = 5 on. Up to
//! degree 5 the two agree, which is what the debug-captured per-function
//! fixtures check.
//! * triangles (303): nodal and edge functions;
//! * quadrilaterals (404): nodal, edge and interior bubble functions, full
//!   scheme (`QuadEdgePBasis`, `QuadBubblePBasis` without `localNumbers`,
//!   PElementBase.F90:915, :981, :1168, :1219) and serendipity
//!   (`SD_QuadEdgePBasis`, `SD_QuadBubblePBasis`, :495, :657, :738, :786).
//!
//! * triangle bubbles of a 2D element (`TriangleBubblePBasis`,
//!   `dTriangleBubblePBasis`, PElementBase.F90).
//!
//! Face-directed bubbles (a 2D element on a 3D face) and degrees above 9
//! keep the H1Basis values.

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

/// `x**n` as gfortran `-O3` expands it: GCC's `powi_table` chain,
/// `x^n = x^(n - t[n]) * x^t[n]` (tree-ssa-math-opts.cc), for n up to 10.
fn pow(x: f64, n: usize) -> f64 {
    const T: [usize; 11] = [0, 0, 1, 2, 2, 3, 3, 4, 4, 6, 5];
    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        _ => pow(x, n - T[n]) * pow(x, T[n]),
    }
}

/// `varPhi(i, x)` (PElementBase.F90:8600), cases 2-9.
fn var_phi(i: usize, x: f64) -> Option<f64> {
    let (x2, x4, x6) = (pow(x, 2), pow(x, 4), pow(x, 6));
    Some(match i {
        2 => -(6f64.sqrt()),
        3 => -(x * 10f64.sqrt()),
        4 => -((5.0 * x2 - 1.0) * 14f64.sqrt() / 4.0),
        5 => -(0.75 * (7.0 * x2 - 3.0) * x * 2f64.sqrt()),
        6 => -((21.0 * x4 - 14.0 * x2 + 1.0) * 22f64.sqrt() / 8.0),
        7 => -((33.0 * x4 - 30.0 * x2 + 5.0) * x * 26f64.sqrt() / 8.0),
        8 => -((429.0 * x6 - 495.0 * x4 + 135.0 * x2 - 5.0) * 30f64.sqrt() / 64.0),
        9 => -((715.0 * x6 - 1001.0 * x4 + 385.0 * x2 - 35.0) * x * 34f64.sqrt() / 64.0),
        _ => return None,
    })
}

/// `dVarPhi(i, x)` (PElementBase.F90:8709), cases 2-9.
fn d_var_phi(i: usize, x: f64) -> Option<f64> {
    let (x2, x4, x6) = (pow(x, 2), pow(x, 4), pow(x, 6));
    let (s26, s34) = (26f64.sqrt(), 34f64.sqrt());
    Some(match i {
        2 => 0.0,
        3 => -(10f64.sqrt()),
        4 => -(2.5 * x * 14f64.sqrt()),
        5 => -(15.75 * x2 * 2f64.sqrt()) + 2.25 * 2f64.sqrt(),
        6 => -(3.5 * x * (3.0 * x2 - 1.0) * 22f64.sqrt()),
        7 => -(20.625 * x4 * s26) + 11.25 * x2 * s26 - 0.625 * s26,
        8 => -(0.28125 * x * (143.0 * x4 - 110.0 * x2 + 15.0) * 30f64.sqrt()),
        9 => -(78.203125 * x6 * s34) + 78.203125 * x4 * s34 - 18.046875 * x2 * s34 + 0.546875 * s34,
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
/// `dPhi` (:8380), cases 2-9.
fn phi(i: usize, x: f64) -> Option<f64> {
    Some(var_phi(i, x)? * (1.0 - x * x) / 4.0)
}

fn d_phi(i: usize, x: f64) -> Option<f64> {
    let (x2, x4, x6, x8) = (pow(x, 2), pow(x, 4), pow(x, 6), pow(x, 8));
    Some(match i {
        2 => 6f64.sqrt() * x / 2.0,
        3 => 10f64.sqrt() * (3.0 * x2 - 1.0) / 4.0,
        4 => 14f64.sqrt() * x * (5.0 * x2 - 3.0) / 4.0,
        5 => 3.0 / 16.0 * 2f64.sqrt() * (35.0 * x4 - 30.0 * x2 + 3.0),
        6 => 22f64.sqrt() * x * (63.0 * x4 - 70.0 * x2 + 15.0) / 16.0,
        7 => 26f64.sqrt() * (231.0 * x6 - 315.0 * x4 + 105.0 * x2 - 5.0) / 32.0,
        8 => 30f64.sqrt() * x * (429.0 * x6 - 693.0 * x4 + 315.0 * x2 - 35.0) / 32.0,
        9 => 34f64.sqrt() * (6435.0 * x8 - 12012.0 * x6 + 6930.0 * x4 - 1260.0 * x2 + 35.0) / 256.0,
        _ => return None,
    })
}

/// `QuadBubblePBasis(i, j, u, v)` and `dQuadBubblePBasis` without
/// `localNumbers`: `Phi(i+2, u) Phi(j+2, v)`.
pub(crate) fn quad_bubble(i: usize, j: usize, u: f64, v: f64) -> Option<(f64, [f64; 2])> {
    let (pu, pv) = (phi(i + 2, u)?, phi(j + 2, v)?);
    Some((pu * pv, [d_phi(i + 2, u)? * pv, pu * d_phi(j + 2, v)?]))
}

/// `SD_QuadEdgePBasis(edge, i, u, v, invert)` and `SD_dQuadEdgePBasis`.
pub fn sd_quad_edge(
    edge: usize,
    i: usize,
    u: f64,
    v: f64,
    invert: bool,
) -> Option<(f64, [f64; 2])> {
    // The edge's variable runs along u (edges 1, 3) or v (2, 4), reversed
    // when inverted; `c` is the cross factor 1/2 (1 -/+ v) or 1/2 (1 -/+ u).
    let s = if invert { -1.0 } else { 1.0 };
    Some(match edge {
        1 => {
            let x = s * u;
            (
                0.5 * (1.0 - v) * phi(i, x)?,
                [s * (0.5 * (1.0 - v) * d_phi(i, x)?), -(0.5 * phi(i, x)?)],
            )
        }
        2 => {
            let x = s * v;
            (
                0.5 * (1.0 + u) * phi(i, x)?,
                [0.5 * phi(i, x)?, s * (0.5 * (u + 1.0) * d_phi(i, x)?)],
            )
        }
        3 => {
            let x = s * u;
            (
                0.5 * (1.0 + v) * phi(i, x)?,
                [s * (0.5 * (1.0 + v) * d_phi(i, x)?), 0.5 * phi(i, x)?],
            )
        }
        _ => {
            let x = s * v;
            (
                0.5 * (1.0 - u) * phi(i, x)?,
                [-(0.5 * phi(i, x)?), s * (0.5 * (1.0 - u) * d_phi(i, x)?)],
            )
        }
    })
}

/// `SD_QuadBubblePBasis(i, j, u, v)` and `SD_dQuadBubblePBasis` without
/// `localNumbers`: `Phi(i, u) Phi(j, v)`.
pub(crate) fn sd_quad_bubble(i: usize, j: usize, u: f64, v: f64) -> Option<(f64, [f64; 2])> {
    let (pu, pv) = (phi(i, u)?, phi(j, v)?);
    Some((pu * pv, [d_phi(i, u)? * pv, pu * d_phi(j, v)?]))
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

/// `toExp(x, n)` (PElementBase.F90): 1 for `n == 0`, 0 for `x == 0`, else
/// `x**n` (a negative `n` as `1 / x**-n`).
fn to_exp(x: f64, n: i32) -> f64 {
    if n == 0 {
        1.0
    } else if x == 0.0 {
        0.0
    } else if n < 0 {
        1.0 / pow(x, (-n) as usize)
    } else {
        pow(x, n as usize)
    }
}

/// `TriangleBubblePBasis(i, j, u, v)` and `dTriangleBubblePBasis` without
/// `localNumbers` (a 2D element's bubble).
pub(crate) fn triangle_bubble(i: usize, j: usize, u: f64, v: f64) -> (f64, [f64; 2]) {
    let (la, lb, lc) = (l(1, u, v), l(2, u, v), l(3, u, v));
    let s3 = 3f64.sqrt();
    let dla = [-1.0 / 2.0, -s3 / 6.0];
    let dlb = [1.0 / 2.0, -s3 / 6.0];
    let dlc = [0.0, s3 / 3.0];
    let value = la * lb * lc * pow(lb - la, i) * pow(2.0 * lc - 1.0, j);
    let (ji, ni) = (i as i32, j as i32);
    let lb_laj = to_exp(lb - la, ji);
    let lc_1n = to_exp(2.0 * lc - 1.0, ni);
    let mut grad = [0.0; 2];
    for k in 0..2 {
        grad[k] = dla[k] * lb * lc * lb_laj * lc_1n
            + la * dlb[k] * lc * lb_laj * lc_1n
            + la * lb * dlc[k] * lb_laj * lc_1n
            + la * lb * lc * f64::from(ji) * to_exp(lb - la, ji - 1) * (dlb[k] - dla[k]) * lc_1n
            + la * lb * lc * lb_laj * f64::from(ni) * to_exp(2.0 * lc - 1.0, ni - 1) * (2.0 * dlc[k]);
    }
    (value, grad)
}

#[cfg(test)]
mod tests {
    /// The scalar triangle bubble is `H1Basis`'s away from `Lb == La` and
    /// stays finite on it (`toExp`), where `0 * (Lb-La)**(-1)` is NaN.
    #[test]
    fn triangle_bubble_is_finite_on_the_symmetry_line() {
        let (v0, g0) = triangle_bubble(0, 0, 0.0, 0.1);
        assert!(v0.is_finite() && g0.iter().all(|x| x.is_finite()), "{v0} {g0:?}");
        for &(u, v) in &[(0.3, 0.4), (-0.2, 0.9)] {
            let mut val = vec![0.0; 1];
            let mut grad = vec![[0.0; 3]; 1];
            let mut nb = 0;
            crate::triangle::triangle_bubble_p(u, v, 3, &mut val, &mut nb, None);
            nb = 0;
            crate::triangle::d_triangle_bubble_p(u, v, 3, &mut grad, &mut nb, None);
            let (b, g) = triangle_bubble(0, 0, u, v);
            assert!((b - val[0]).abs() < 1e-15, "{b} {}", val[0]);
            assert!((g[0] - grad[0][0]).abs() < 1e-14 && (g[1] - grad[0][1]).abs() < 1e-14, "{g:?} {:?}", grad[0]);
        }
    }

    use super::*;

    /// The hand-expanded scalar gradients are the derivatives of the values.
    #[test]
    fn edge_gradients_differentiate_the_values() {
        let (u, v, h) = (0.13, 0.41, 1e-6);
        for edge in 1..=3 {
            for i in 2..=9 {
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
        assert!(triangle_edge(1, 10, u, v, false).is_none());
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
            for i in 2..=9 {
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
        assert!(quad_bubble(8, 0, u, v).is_none());
        for edge in 1..=4 {
            for i in 2..=9 {
                for invert in [false, true] {
                    check(
                        &|u, v| sd_quad_edge(edge, i, u, v, invert).unwrap(),
                        &format!("sd edge {edge} i {i} {invert}"),
                    );
                }
            }
        }
        for i in 2..=3 {
            for j in 2..=3 {
                check(
                    &|u, v| sd_quad_bubble(i, j, u, v).unwrap(),
                    &format!("sd bubble {i} {j}"),
                );
            }
        }
    }
}
