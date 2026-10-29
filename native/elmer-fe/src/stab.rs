//! Stabilization parameters: `ElementDiameter` (ElementBasis.F90:3051) and
//! `StabParam` (ElemInfo.F90:152), the `hK` and `mK` of the convection and
//! flow stabilizations, as `MeshStabParams` (MeshBasics.F90:2318) sets them
//! per element.

use elmer_core::{BasisEval, ElementType};

use crate::element_type::{element_type, ElementTypeDef};
use crate::error::FeError;
use crate::info::element_info;
use crate::metric::MetricStatus;
use crate::quadrature::gauss_points;

/// Elmer's `AEPS = 10 * EPSILON(1.0_dp)` (Types.F90:71).
const AEPS: f64 = 10.0 * f64::EPSILON;

/// `GetEdgeMap` (ElementBasis.F90:2937), 1-based node pairs.
fn edge_map(family: u16) -> &'static [[usize; 2]] {
    match family {
        2 => &[[1, 2]],
        3 => &[[1, 2], [2, 3], [3, 1]],
        4 => &[[1, 2], [2, 3], [3, 4], [4, 1]],
        5 => &[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]],
        6 => &[
            [1, 2],
            [2, 3],
            [4, 3],
            [1, 4],
            [1, 5],
            [2, 5],
            [3, 5],
            [4, 5],
        ],
        7 => &[
            [1, 2],
            [2, 3],
            [3, 1],
            [4, 5],
            [5, 6],
            [6, 4],
            [1, 4],
            [2, 5],
            [3, 6],
        ],
        _ => &[
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
    }
}

/// `ElementDiameter( elm, nodes, UseLongEdge )`.
///
/// Elmer keeps `LongEdge` in a SAVEd local, so a call without the argument
/// reuses the last value passed; here the flag is always explicit.
pub fn element_diameter(etype: ElementType, xyz: &[[f64; 3]], use_long_edge: bool) -> f64 {
    let (x, y, z) = (
        |i: usize| xyz[i - 1][0],
        |i: usize| xyz[i - 1][1],
        |i: usize| xyz[i - 1][2],
    );
    let hk = match etype.family() {
        1 => 0.0,
        3 => {
            let j11 = x(2) - x(1);
            let j12 = y(2) - y(1);
            let j13 = z(2) - z(1);
            let j21 = x(3) - x(1);
            let j22 = y(3) - y(1);
            let j23 = z(3) - z(1);
            let g11 = j11 * j11 + j12 * j12 + j13 * j13;
            let g12 = j11 * j21 + j12 * j22 + j13 * j23;
            let g22 = j21 * j21 + j22 * j22 + j23 * j23;
            let a = (g11 * g22 - g12 * g12).sqrt() / 2.0;

            let cx = (x(1) + x(2) + x(3)) / 3.0;
            let cy = (y(1) + y(2) + y(3)) / 3.0;
            let cz = (z(1) + z(2) + z(3)) / 3.0;

            let sq = |d: f64| d * d;
            let mut s = sq(x(1) - cx) + sq(y(1) - cy) + sq(z(1) - cz);
            s = s + sq(x(2) - cx) + sq(y(2) - cy) + sq(z(2) - cz);
            s = s + sq(x(3) - cx) + sq(y(3) - cy) + sq(z(3) - cz);

            16.0 * a * a / (3.0 * s)
        }
        4 => {
            let sq = |d: f64| d * d;
            let cx = sq(x(2) - x(1)) + sq(y(2) - y(1)) + sq(z(2) - z(1));
            let cy = sq(x(4) - x(1)) + sq(y(4) - y(1)) + sq(z(4) - z(1));
            2.0 * cx * cy / (cx + cy)
        }
        f => {
            let mut hk = if use_long_edge { -f64::MAX } else { f64::MAX };
            for &[j, k] in edge_map(f) {
                let x0 = x(j) - x(k);
                let y0 = y(j) - y(k);
                let z0 = z(j) - z(k);
                let d = x0 * x0 + y0 * y0 + z0 * z0;
                // MAX/MIN( hK, d ): the second argument wins only when strictly better.
                if use_long_edge {
                    if d > hk {
                        hk = d;
                    }
                } else if d < hk {
                    hk = d;
                }
            }
            hk
        }
    };
    hk.sqrt()
}

/// The `L` and `G` matrices of `StabParam` (order `n-1`, row-major), or
/// `None` when the element is degenerate.
fn stab_matrices(
    def: &ElementTypeDef,
    xyz: &[[f64; 3]],
    cdim: usize,
) -> Result<(Vec<f64>, Vec<f64>), FeError> {
    let n = def.n_nodes;
    let m = n - 1;
    let mut status = MetricStatus::Unchecked;
    let mut ev = BasisEval::default();
    // dnodal[p][k][i]: d(phi_p)/dx_i at node k.
    let mut dnodal = vec![vec![[0.0f64; 3]; n]; n];
    for k in 0..n {
        let uvw = [def.node_u[k], def.node_v[k], def.node_w[k]];
        element_info(def, xyz, uvw, cdim, &mut status, &mut ev)?;
        for (dn, d) in dnodal.iter_mut().zip(&ev.dbasis_dx) {
            dn[k] = *d;
        }
    }
    let dim = cdim;
    let rule = gauss_points(ElementType(def.code), None)?;
    let mut l = vec![0.0; m * m];
    let mut g = vec![0.0; m * m];
    for t in 0..rule.len() {
        let uvw = [rule.u[t], rule.v[t], rule.w[t]];
        element_info(def, xyz, uvw, cdim, &mut status, &mut ev)?;
        let s = ev.det_j * rule.s[t];
        let db = &ev.dbasis_dx;
        for p in 1..n {
            for q in 1..n {
                let mut ddp = [0.0; 3];
                let mut ddq = [0.0; 3];
                for i in 0..dim {
                    g[(p - 1) * m + q - 1] += s * db[p][i] * db[q][i];
                    let mut sp = 0.0;
                    let mut sq = 0.0;
                    for (dn, d) in dnodal[p].iter().zip(db) {
                        sp += dn[i] * d[i];
                    }
                    for (dn, d) in dnodal[q].iter().zip(db) {
                        sq += dn[i] * d[i];
                    }
                    ddp[i] += sp;
                    ddq[i] += sq;
                }
                let sum = |d: [f64; 3]| 0.0 + d[0] + d[1] + d[2];
                l[(p - 1) * m + q - 1] += s * sum(ddp) * sum(ddq);
            }
        }
    }
    Ok((l, g))
}

/// Largest `lambda` of `L x = lambda G x` with `G` symmetric positive
/// definite (the `DSYGV` call of `StabParam`): Cholesky `G = U^T U`, then
/// cyclic Jacobi on `U^-T L U^-1`. Agrees with LAPACK to rounding, not bit
/// for bit.
fn largest_generalized_eigenvalue(l: &[f64], g: &[f64], m: usize) -> f64 {
    // Upper Cholesky factor, row-major u[i*m+j] for j >= i.
    let mut u = vec![0.0; m * m];
    for j in 0..m {
        let mut d = g[j * m + j];
        for k in 0..j {
            d -= u[k * m + j] * u[k * m + j];
        }
        let d = d.sqrt();
        u[j * m + j] = d;
        for i in j + 1..m {
            let mut s = g[j * m + i];
            for k in 0..j {
                s -= u[k * m + j] * u[k * m + i];
            }
            u[j * m + i] = s / d;
        }
    }
    // C = U^-T L U^-1: solve U^T Y = L, then C^T = U^-T Y^T.
    let solve_ut = |b: &mut [f64]| {
        for i in 0..m {
            let mut s = b[i];
            for k in 0..i {
                s -= u[k * m + i] * b[k];
            }
            b[i] = s / u[i * m + i];
        }
    };
    let mut y = vec![0.0; m * m];
    for c in 0..m {
        let mut col: Vec<f64> = (0..m).map(|r| l[r * m + c]).collect();
        solve_ut(&mut col);
        for r in 0..m {
            y[r * m + c] = col[r];
        }
    }
    let mut a = vec![0.0; m * m];
    for r in 0..m {
        let mut row: Vec<f64> = y[r * m..r * m + m].to_vec();
        solve_ut(&mut row);
        for c in 0..m {
            a[c * m + r] = row[c];
        }
    }
    for i in 0..m {
        for j in i + 1..m {
            let s = 0.5 * (a[i * m + j] + a[j * m + i]);
            a[i * m + j] = s;
            a[j * m + i] = s;
        }
    }
    jacobi_eigenvalues(&mut a, m)
        .into_iter()
        .fold(f64::NEG_INFINITY, f64::max)
}

fn jacobi_eigenvalues(a: &mut [f64], m: usize) -> Vec<f64> {
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..m {
            for j in i + 1..m {
                off += a[i * m + j] * a[i * m + j];
            }
        }
        let mut diag = 0.0;
        for i in 0..m {
            diag += a[i * m + i] * a[i * m + i];
        }
        if off <= 1.0e-36 * diag {
            break;
        }
        for p in 0..m {
            for q in p + 1..m {
                let apq = a[p * m + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * m + q] - a[p * m + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..m {
                    let akp = a[k * m + p];
                    let akq = a[k * m + q];
                    a[k * m + p] = c * akp - s * akq;
                    a[k * m + q] = s * akp + c * akq;
                }
                for k in 0..m {
                    let apk = a[p * m + k];
                    let aqk = a[q * m + k];
                    a[p * m + k] = c * apk - s * aqk;
                    a[q * m + k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..m).map(|i| a[i * m + i]).collect()
}

/// The eigenvalue part of `StabParam`: `Some(lambda_max)`, or `None` where
/// Elmer falls back to `mK = 1/3` (all of `L` below `AEPS`, or `lambda_max`
/// below `10*AEPS`).
fn stab_eigenvalue(
    def: &ElementTypeDef,
    xyz: &[[f64; 3]],
    cdim: usize,
) -> Result<Option<f64>, FeError> {
    let (l, g) = stab_matrices(def, xyz, cdim)?;
    if l.iter().all(|x| x.abs() < AEPS) {
        return Ok(None);
    }
    let mk = largest_generalized_eigenvalue(&l, &g, def.n_nodes - 1);
    Ok(if mk < 10.0 * AEPS { None } else { Some(mk) })
}

/// The type-level `StabilizationMK` that `GetElementType` computes on the
/// reference element with `CurrentModel % Dimension` set to the element
/// dimension (the elements.def value is discarded on reading).
pub fn type_stabilization_mk(etype: ElementType) -> Result<f64, FeError> {
    let def = element_type(etype.0)?;
    if let Some(&mk) = def.type_mk.get() {
        return Ok(mk);
    }
    let xyz: Vec<[f64; 3]> = (0..def.n_nodes)
        .map(|k| [def.node_u[k], def.node_v[k], def.node_w[k]])
        .collect();
    let mk = stab_param_mk(def, &xyz, def.dim)?;
    Ok(*def.type_mk.get_or_init(|| mk))
}

/// `StabParam( Element, Nodes, n, mK )` (no `hK`): the `mK` that
/// `GetElementType` stores per type, here for any nodal coordinates.
pub fn stab_param_mk_only(
    etype: ElementType,
    xyz: &[[f64; 3]],
    cdim: usize,
) -> Result<f64, FeError> {
    let def = element_type(etype.0)?;
    stab_param_mk(def, xyz, cdim)
}

fn stab_param_mk(def: &ElementTypeDef, xyz: &[[f64; 3]], cdim: usize) -> Result<f64, FeError> {
    if def.basis_degree <= 1 {
        return Ok(linear_mk(def.code));
    }
    Ok(match stab_eigenvalue(def, xyz, cdim)? {
        None => 1.0 / 3.0,
        Some(mut mk) => {
            if matches!(def.family(), 2 | 4 | 8) {
                mk *= 4.0;
            }
            (1.0f64 / 3.0).min(2.0 / mk)
        }
    })
}

/// The constant `mK` of the linear elements; Elmer leaves `mK` untouched
/// for the other degree-one types, which start at zero.
fn linear_mk(code: u16) -> f64 {
    match code {
        202 | 303 | 404 | 504 | 605 | 706 => 1.0 / 3.0,
        808 => 1.0 / 6.0,
        _ => 0.0,
    }
}

/// `StabParam( Element, Nodes, n, mK, hK, UseLongEdge )` for one element:
/// returns `(mK, hK)` as `MeshStabParams` stores them in
/// `Element % StabilizationMK` and `Element % hK`. `cdim` is
/// `CoordinateSystemDimension()`.
pub fn stab_param(
    etype: ElementType,
    xyz: &[[f64; 3]],
    cdim: usize,
    use_long_edge: bool,
) -> Result<(f64, f64), FeError> {
    let def = element_type(etype.0)?;
    if xyz.len() != def.n_nodes {
        return Err(FeError::NodeCount {
            etype: etype.0,
            expected: def.n_nodes,
            got: xyz.len(),
        });
    }
    if def.basis_degree <= 1 {
        return Ok((
            linear_mk(def.code),
            element_diameter(etype, xyz, use_long_edge),
        ));
    }
    match stab_eigenvalue(def, xyz, cdim)? {
        None => Ok((1.0 / 3.0, element_diameter(etype, xyz, use_long_edge))),
        Some(mk) => {
            let type_mk = type_stabilization_mk(etype)?;
            let hk = (2.0 / (mk * type_mk)).sqrt();
            Ok(((1.0f64 / 3.0).min(type_mk), hk))
        }
    }
}
