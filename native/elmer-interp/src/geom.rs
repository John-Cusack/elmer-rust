//! Point-in-element test and global-to-local mapping.
//!
//! Ports `GlobalToLocal` (`ElementGeometry.F90`) and `PointInElement`
//! (`Interpolation.F90`), including the Newton relaxation schedule and the
//! per-family distance measures.

use crate::shape::Shape;

/// Machine epsilon of `REAL(KIND=dp)`, `EPSILON(1.0_dp)`.
const EPS: f64 = f64::EPSILON;

/// Tolerances of one `PointInElement` call. `None` fields take the Fortran
/// defaults: global `1e-4`, local `1e-10`, numeric `EPSILON(dp)`.
#[derive(Copy, Clone, Debug, Default)]
pub struct PointTol {
    pub global: Option<f64>,
    pub local: Option<f64>,
    pub numeric: Option<f64>,
}

/// Result of a successful or failed element probe.
#[derive(Copy, Clone, Debug)]
pub struct Probe {
    pub inside: bool,
    /// Local coordinates (after the optional Piola remapping).
    pub uvw: [f64; 3],
}

/// Solve a 2x2 system as `SolveLinSys2x2`; a singular matrix leaves `x` zero.
fn solve2(a: [[f64; 2]; 2], b: [f64; 2]) -> [f64; 2] {
    let mut det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det == 0.0 {
        return [0.0; 2];
    }
    det = 1.0 / det;
    [
        det * (a[1][1] * b[0] - a[0][1] * b[1]),
        det * (a[0][0] * b[1] - a[1][0] * b[0]),
    ]
}

/// Solve a 3x3 system exactly as `SolveLinSys3x3` (pivot on the first row).
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let mut x = [0.0; 3];
    if a[0][0].abs() > a[0][1].abs() && a[0][0].abs() > a[0][2].abs() {
        let q = 1.0 / a[0][0];
        let s = q * a[1][0];
        let t = q * a[2][0];
        let c = [
            [a[1][1] - s * a[0][1], a[1][2] - s * a[0][2]],
            [a[2][1] - t * a[0][1], a[2][2] - t * a[0][2]],
        ];
        let y = solve2(c, [b[1] - s * b[0], b[2] - t * b[0]]);
        x[1] = y[0];
        x[2] = y[1];
        x[0] = q * (b[0] - a[0][1] * x[1] - a[0][2] * x[2]);
    } else if a[0][1].abs() > a[0][2].abs() {
        let q = 1.0 / a[0][1];
        let s = q * a[1][1];
        let t = q * a[2][1];
        let c = [
            [a[1][0] - s * a[0][0], a[1][2] - s * a[0][2]],
            [a[2][0] - t * a[0][0], a[2][2] - t * a[0][2]],
        ];
        let y = solve2(c, [b[1] - s * b[0], b[2] - t * b[0]]);
        x[0] = y[0];
        x[2] = y[1];
        x[1] = q * (b[0] - a[0][0] * x[0] - a[0][2] * x[2]);
    } else {
        let q = 1.0 / a[0][2];
        let s = q * a[1][2];
        let t = q * a[2][2];
        let c = [
            [a[1][0] - s * a[0][0], a[1][1] - s * a[0][1]],
            [a[2][0] - t * a[0][0], a[2][1] - t * a[0][1]],
        ];
        let y = solve2(c, [b[1] - s * b[0], b[2] - t * b[0]]);
        x[0] = y[0];
        x[1] = y[1];
        x[2] = q * (b[0] - a[0][0] * x[0] - a[0][1] * x[1]);
    }
    x
}

/// `GlobalToLocal`: Newton iteration for the local coordinates of `p` in an
/// element given by its nodal coordinates. `coord_dim` is
/// `CoordinateSystemDimension()` and only matters for surface elements.
pub fn global_to_local(
    shape: &Shape,
    nodes: &[[f64; 3]],
    p: [f64; 3],
    coord_dim: usize,
) -> [f64; 3] {
    const MAXITER: usize = 50;
    let mut uvw = [0.0f64; 3];
    let edim = shape.dim();
    if edim == 0 {
        return uvw;
    }
    let n = shape.n_nodes();
    let nodes = &nodes[..n];
    let ext = |d: usize| {
        let mx = nodes.iter().map(|x| x[d]).fold(f64::NEG_INFINITY, f64::max);
        let mn = nodes.iter().map(|x| x[d]).fold(f64::INFINITY, f64::min);
        mx - mn
    };
    let scl = ext(0) + ext(1) + ext(2);
    let acc = EPS * scl * scl;

    let mut basis = vec![0.0; n];
    let mut dbasis = vec![[0.0; 3]; n];
    let mut delta = [0.0f64; 3];
    for i in 1..=MAXITER {
        shape.basis(uvw, &mut basis);
        let mut res = [0.0; 3];
        for d in 0..3 {
            res[d] = basis.iter().zip(nodes).map(|(b, x)| b * x[d]).sum::<f64>() - p[d];
        }
        let err = res[0] * res[0] + res[1] * res[1] + res[2] * res[2];
        if err < acc {
            break;
        }
        let prev = delta;
        delta = [0.0; 3];
        shape.dbasis(uvw, &mut dbasis);
        // J[r][c] = sum_k dN_k/d(c) * x_k[r]
        let jac = |r: usize, c: usize| -> f64 {
            dbasis.iter().zip(nodes).map(|(g, x)| g[c] * x[r]).sum()
        };
        match edim {
            1 => {
                let j = [jac(0, 0), jac(1, 0), jac(2, 0)];
                let det = j[0] * j[0] + j[1] * j[1] + j[2] * j[2];
                delta[0] = (res[0] * j[0] + res[1] * j[1] + res[2] * j[2]) / det;
            }
            2 => {
                let (mut r, mut s) = (res[0], res[1]);
                let mut j = [[jac(0, 0), jac(0, 1)], [jac(1, 0), jac(1, 1)]];
                if coord_dim == 3 {
                    let j3 = [
                        [j[0][0], j[0][1]],
                        [j[1][0], j[1][1]],
                        [jac(2, 0), jac(2, 1)],
                    ];
                    r = j3[0][0] * res[0] + j3[1][0] * res[1] + j3[2][0] * res[2];
                    s = j3[0][1] * res[0] + j3[1][1] * res[1] + j3[2][1] * res[2];
                    for a in 0..2 {
                        for b in 0..2 {
                            j[a][b] = (0..3).map(|k| j3[k][a] * j3[k][b]).sum();
                        }
                    }
                }
                let d = solve2(j, [r, s]);
                delta[0] = d[0];
                delta[1] = d[1];
            }
            _ => {
                let j = [
                    [jac(0, 0), jac(0, 1), jac(0, 2)],
                    [jac(1, 0), jac(1, 1), jac(1, 2)],
                    [jac(2, 0), jac(2, 1), jac(2, 2)],
                ];
                delta = solve3(j, res);
            }
        }
        if i > 10 {
            let change: f64 = (0..3).map(|d| (delta[d] - prev[d]).abs()).sum();
            // Same correction proposed repeatedly: the point is off-element.
            if change < EPS {
                break;
            }
            if i > 20 && change < 1.0e-8 {
                break;
            }
            for d in &mut delta {
                *d *= 0.5;
            }
        }
        for d in 0..3 {
            uvw[d] -= delta[d];
        }
    }
    uvw
}

/// Distance of local coordinates from the reference element, summed over
/// the violated bounds of each family (`PointInElement`).
fn local_distance(family: u16, [u, v, w]: [f64; 3]) -> f64 {
    let ext = |x: f64| (x - 1.0).max((-x - 1.0).max(0.0));
    let simplex = |x: f64, y: f64| (-x).max(0.0) + (-y).max(0.0) + (x + y - 1.0).max(0.0);
    match family {
        1 => u,
        2 => ext(u),
        3 => simplex(u, v),
        4 => ext(u) + ext(v),
        5 => (-u).max(0.0) + (-v).max(0.0) + (-w).max(0.0) + (u + v + w - 1.0).max(0.0),
        // The Fortran has no pyramid case (its `sumdist` is left undefined).
        // Use the pyramid's own bounds: 0 <= w <= 1, |u|,|v| <= 1 - w.
        6 => {
            (-w).max(0.0)
                + (w - 1.0).max(0.0)
                + (u.abs() - (1.0 - w)).max(0.0)
                + (v.abs() - (1.0 - w)).max(0.0)
        }
        7 => simplex(u, v) + ext(w),
        8 => ext(u) + ext(v) + ext(w),
        _ => f64::INFINITY,
    }
}

/// `PointInElement`. `nodes` are the element's nodal coordinates (already
/// flattened along a passive coordinate by the caller); `piola` applies the
/// edge-basis remapping of local coordinates (`EdgeBasis=.TRUE.`).
pub fn point_in_element(
    shape: &Shape,
    nodes: &[[f64; 3]],
    p: [f64; 3],
    tol: PointTol,
    piola: bool,
    coord_dim: usize,
) -> Probe {
    let n = shape.n_nodes();
    let eps0 = tol.numeric.unwrap_or(EPS);
    let eps1 = tol.global.unwrap_or(1.0e-4);
    let eps2 = tol.local.unwrap_or(1.0e-10);
    let miss = Probe {
        inside: false,
        uvw: [0.0; 3],
    };
    if eps1 >= 0.0 {
        for d in 0..3 {
            let mn = nodes[..n]
                .iter()
                .map(|x| x[d])
                .fold(f64::INFINITY, f64::min);
            let mx = nodes[..n]
                .iter()
                .map(|x| x[d])
                .fold(f64::NEG_INFINITY, f64::max);
            let dist = (p[d] - mx).max(0.0).max(mn - p[d]);
            if dist > eps0 + eps1 * (mx - mn) {
                return miss;
            }
        }
    }
    let mut uvw = global_to_local(shape, nodes, p, coord_dim);
    let family = shape.code() / 100;
    let inside = local_distance(family, uvw) < eps0 + eps2;
    if piola {
        let [u, v, w] = uvw;
        match family {
            3 | 7 => {
                uvw[0] = 2.0 * u + v - 1.0;
                uvw[1] = 3f64.sqrt() * v;
            }
            5 => {
                uvw[0] = 2.0 * u + v + w - 1.0;
                uvw[1] = 3f64.sqrt() * v + w / 3f64.sqrt();
                uvw[2] = 2.0 * (2.0f64 / 3.0).sqrt() * w;
            }
            6 => uvw[2] = 2f64.sqrt() * w,
            _ => {}
        }
    }
    Probe { inside, uvw }
}
