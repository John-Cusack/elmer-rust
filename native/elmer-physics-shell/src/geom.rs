//! The shell geometry of ShellSolver.F90 for quadrilateral background
//! elements: the curved edges (`HermiteForm`), the blending surface
//! (`BlendingSurfaceInfo`), the lines-of-curvature frame and its Taylor
//! polynomial (`LinesOfCurvatureFrame`), the reparametrised patch
//! (`LinesOfCurvaturePatch`) and the surface basis of the patch
//! (`SurfaceBasisVectors`).
//!
//! Every sum is taken in the Fortran's order (`SUM`, `DOT_PRODUCT` from
//! zero, array expressions left to right).

use elmer_fe::basis::nodal_basis;
use elmer_fe::element_type;

use crate::ShellError;

pub(crate) type V3 = [f64; 3];

/// `GeometryEpsilon`, `GeometryMaxIters` and `ParabolicDelta`.
const GEOMETRY_EPSILON: f64 = 1.0e-6;
const GEOMETRY_MAX_ITERS: usize = 50;
const PARABOLIC_DELTA: f64 = 1.0e-8;
/// `MaxPatchNodes`.
pub(crate) const MAX_PATCH_NODES: usize = 16;
/// `EPSILON(1.0)`: single precision.
const EPS_SP: f64 = f32::EPSILON as f64;

/// `DOT_PRODUCT( a, b )`.
pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `SQRT( SUM( a**2 ) )`.
pub(crate) fn norm(a: &[f64]) -> f64 {
    let mut s = 0.0;
    for x in a {
        s += x * x;
    }
    s.sqrt()
}

/// `CrossProduct( v1, v2 )`.
pub(crate) fn cross(v1: V3, v2: V3) -> V3 {
    [
        v1[1] * v2[2] - v1[2] * v2[1],
        -v1[0] * v2[2] + v1[2] * v2[0],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
}

/// `c * v` elementwise.
pub(crate) fn scale(c: f64, v: V3) -> V3 {
    [c * v[0], c * v[1], c * v[2]]
}

/// `SolveLinSys2x2( A, x, b )` (`A(i,j)` = `a[i][j]`).
pub(crate) fn solve2x2(a: [[f64; 2]; 2], b: [f64; 2]) -> Result<[f64; 2], ShellError> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det == 0.0 {
        return Err(ShellError::Fatal("SolveLinSys2x2: Singular matrix, sorry!".into()));
    }
    let det = 1.0 / det;
    Ok([
        det * (a[1][1] * b[0] - a[0][1] * b[1]),
        det * (a[0][0] * b[1] - a[1][0] * b[0]),
    ])
}

/// `QuadNodalPBasis( q, u, v )` and `dQuadNodalPBasis( q, u, v )`.
pub(crate) fn quad_p_basis(u: f64, v: f64) -> ([f64; 4], [[f64; 2]; 4]) {
    (
        [
            (1.0 - u) * (1.0 - v) / 4.0,
            (1.0 + u) * (1.0 - v) / 4.0,
            (1.0 + u) * (1.0 + v) / 4.0,
            (1.0 - u) * (1.0 + v) / 4.0,
        ],
        [
            [-(1.0 - v) / 4.0, -(1.0 - u) / 4.0],
            [(1.0 - v) / 4.0, -(1.0 + u) / 4.0],
            [(1.0 + v) / 4.0, (1.0 + u) / 4.0],
            [-(1.0 + v) / 4.0, (1.0 - u) / 4.0],
        ],
    )
}

/// `HermiteForm( X1, X2, d1, d2, cpars )`: the end tangents of the cubic
/// Hermite curve of an edge, normal to the end directors.
pub(crate) fn hermite_form(x1: V3, x2: V3, d1: V3, d2: V3) -> [f64; 6] {
    let v21 = [x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]];
    let h = norm(&v21);
    let b = scale(1.0 / h, v21);
    let tangent = |d: V3| {
        let c = dot(&v21, &d);
        let t = [v21[0] - c * d[0], v21[1] - c * d[1], v21[2] - c * d[2]];
        let t = scale(1.0 / norm(&t), t);
        scale(h / dot(&t, &b), t)
    };
    let (t1, t2) = (tangent(d1), tangent(d2));
    [t1[0], t1[1], t1[2], t2[0], t2[1], t2[2]]
}

/// `HermiteBasis( s, h, ..., 2 )`: the two-node cubic Hermite basis and its
/// first and second derivatives.
fn hermite_basis(s: f64, h: f64) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let s2 = s.powi(2);
    let s3 = s.powi(3);
    (
        [
            (2.0 - 3.0 * s + s3) / 4.0,
            (2.0 + 3.0 * s - s3) / 4.0,
            h / 8.0 - (h * s) / 8.0 - (h * s2) / 8.0 + (h * s3) / 8.0,
            -h / 8.0 - (h * s) / 8.0 + (h * s2) / 8.0 + (h * s3) / 8.0,
        ],
        [
            (-3.0 + 3.0 * s2) / 4.0,
            (3.0 - 3.0 * s2) / 4.0,
            -h / 8.0 - (h * s) / 4.0 + (3.0 * h * s2) / 8.0,
            -h / 8.0 + (h * s) / 4.0 + (3.0 * h * s2) / 8.0,
        ],
        [
            (3.0 * s) / 2.0,
            (-3.0 * s) / 2.0,
            (h * (-1.0 + 3.0 * s)) / 4.0,
            (h + 3.0 * h * s) / 4.0,
        ],
    )
}

/// What `BlendingSurfaceInfo` returns: the point, the covariant basis, the
/// unit normal, the metric and the second fundamental form.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SurfacePoint {
    pub x: V3,
    pub a1: V3,
    pub a2: V3,
    #[allow(dead_code)]
    pub a3: V3,
    pub a: [[f64; 2]; 2],
    pub deta: f64,
    pub b: [[f64; 2]; 2],
}

/// `BlendingSurfaceInfo( Element, Nodes, u, v, ... )` of a four-node
/// quadrilateral: the bilinear map corrected by the Hermite curves of its
/// edges (`edge`, the six parameters of each of the four edges).
pub(crate) fn blending_surface_info(nodes: &[V3; 4], edge: &[[f64; 6]; 4], u: f64, v: f64) -> SurfacePoint {
    let (basis, dbasis) = quad_p_basis(u, v);
    let mut x = [0.0f64; 3];
    let mut a1 = [0.0f64; 3];
    let mut a2 = [0.0f64; 3];
    for k in 0..3 {
        let (mut s, mut s1, mut s2) = (0.0, 0.0, 0.0);
        for q in 0..4 {
            s += nodes[q][k] * basis[q];
        }
        for q in 0..4 {
            s1 += dbasis[q][0] * nodes[q][k];
        }
        for q in 0..4 {
            s2 += dbasis[q][1] * nodes[q][k];
        }
        x[k] = s;
        a1[k] = s1;
        a2[k] = s2;
    }
    let mut d1a1 = [0.0f64; 3];
    let mut d2a1 = [0.0f64; 3];
    for k in 0..3 {
        d2a1[k] = 1.0 / 4.0 * (nodes[0][k] - nodes[1][k] + nodes[2][k] - nodes[3][k]);
    }
    let mut d2a2 = [0.0f64; 3];
    let h = 2.0;
    for e in 0..4 {
        let d = edge[e];
        let (r1, r2, s, t) = match e {
            0 => (nodes[0], nodes[1], u, v),
            1 => (nodes[1], nodes[2], v, u),
            2 => (nodes[3], nodes[2], u, v),
            _ => (nodes[0], nodes[3], v, u),
        };
        let (hb, dhb, ddhb) = hermite_basis(s, h);
        let h1 = 0.5 * (1.0 - t);
        let h2 = 0.5 * (1.0 + t);
        let (dh1, dh2, ddh1, ddh2) = (-0.5, 0.5, 0.0, 0.0);
        let l1 = 0.5 * (1.0 - s);
        let l2 = 0.5 * (1.0 + s);
        let (dl1, dl2) = (-0.5, 0.5);
        let mut c = [0.0f64; 3];
        let mut dc = [0.0f64; 3];
        let mut ddc = [0.0f64; 3];
        for k in 0..3 {
            c[k] = r1[k] * hb[0] + r2[k] * hb[1] + d[k] * 0.5 * hb[2] + d[3 + k] * 0.5 * hb[3];
            dc[k] = r1[k] * dhb[0] + r2[k] * dhb[1] + d[k] * 0.5 * dhb[2] + d[3 + k] * 0.5 * dhb[3];
            ddc[k] = r1[k] * ddhb[0] + r2[k] * ddhb[1] + d[k] * 0.5 * ddhb[2] + d[3 + k] * 0.5 * ddhb[3];
        }
        for k in 0..3 {
            let b12 = c[k] - l1 * r1[k] - l2 * r2[k];
            let db12 = dc[k] - dl1 * r1[k] - dl2 * r2[k];
            match e {
                0 => {
                    x[k] += h1 * b12;
                    a1[k] += h1 * db12;
                    a2[k] += dh1 * b12;
                    d1a1[k] += h1 * ddc[k];
                    d2a1[k] += dh1 * db12;
                    d2a2[k] += ddh1 * b12;
                }
                1 => {
                    x[k] += h2 * b12;
                    a1[k] += dh2 * b12;
                    a2[k] += h2 * db12;
                    d1a1[k] += ddh2 * b12;
                    d2a1[k] += dh2 * db12;
                    d2a2[k] += h2 * ddc[k];
                }
                2 => {
                    x[k] += h2 * b12;
                    a1[k] += h2 * db12;
                    a2[k] += dh2 * b12;
                    d1a1[k] += h2 * ddc[k];
                    d2a1[k] += dh2 * db12;
                    d2a2[k] += ddh2 * b12;
                }
                _ => {
                    x[k] += h1 * b12;
                    a1[k] += dh1 * b12;
                    a2[k] += h1 * db12;
                    d1a1[k] += ddh1 * b12;
                    d2a1[k] += dh1 * db12;
                    d2a2[k] += h1 * ddc[k];
                }
            }
        }
    }
    let a11 = dot(&a1, &a1);
    let a12 = dot(&a1, &a2);
    let a22 = dot(&a2, &a2);
    let deta = a11 * a22 - a12 * a12;
    let a3 = cross(a1, a2);
    let a3 = scale(1.0 / norm(&a3), a3);
    let b11 = dot(&a3, &d1a1);
    let b12 = dot(&a3, &d2a1);
    let b22 = dot(&a3, &d2a2);
    SurfacePoint {
        x,
        a1,
        a2,
        a3,
        a: [[a11, a12], [a12, a22]],
        deta,
        b: [[b11, b12], [b12, b22]],
    }
}

/// The elementwise geometry `LinesOfCurvatureFrame( ..., SaveProperties )`
/// and `LinesOfCurvaturePatch` store (`element frame`, `taylor
/// parameters`, `planar point`, `umbilical point`, `patch nodes`).
#[derive(Clone, Debug)]
pub(crate) struct Frame {
    pub e1: V3,
    pub e2: V3,
    pub e3: V3,
    /// The frame origin (`XGlob` of the patch; not needed by the kernel).
    #[allow(dead_code)]
    pub o: V3,
    pub taylor: [f64; 6],
    pub planar: bool,
    pub umbilical: bool,
    /// The 16 patch nodes `(y1, y2)` of the 416 geometry element.
    pub patch: [[f64; 2]; MAX_PATCH_NODES],
}

/// The node coordinates of the 416 geometry element as
/// `LinesOfCurvatureFrame` sets them (`NodeU(5:16)`, `NodeV(5:16)`).
fn geometry_nodes() -> Result<([f64; 16], [f64; 16]), ShellError> {
    let def = element_type(416).map_err(|e| ShellError::Fatal(e.to_string()))?;
    let mut nu = [0.0f64; 16];
    let mut nv = [0.0f64; 16];
    for i in 0..4 {
        nu[i] = def.node_u[i];
        nv[i] = def.node_v[i];
    }
    let t = 1.0 / 3.0;
    let tail_u = [-t, t, 1.0, 1.0, t, -t, -1.0, -1.0, -t, t, t, -t];
    let tail_v = [-1.0, -1.0, -t, t, 1.0, 1.0, t, -t, -t, -t, t, t];
    nu[4..].copy_from_slice(&tail_u);
    nv[4..].copy_from_slice(&tail_v);
    Ok((nu, nv))
}

/// `NodalBasisFunctions2D( GBasis, GElement, u, v )` of the 416 element.
pub(crate) fn geometry_basis(u: f64, v: f64) -> Result<[f64; 16], ShellError> {
    let def = element_type(416).map_err(|e| ShellError::Fatal(e.to_string()))?;
    let mut b = [0.0f64; 16];
    nodal_basis(def, u, v, 0.0, &mut b);
    Ok(b)
}

/// `LinesOfCurvatureFrame( Element, TaylorApproximation=.TRUE.,
/// LagrangeNodes=..., d=d, PlanarSurface=plate, ..., SaveProperties=.TRUE. )`
/// at the element centre, then `LinesOfCurvaturePatch`. `spherical` is the
/// body's `Spherical Body`.
pub(crate) fn element_frame(
    nodes: &[V3; 4],
    edge: &[[f64; 6]; 4],
    d: V3,
    plate: bool,
    spherical: bool,
) -> Result<Frame, ShellError> {
    let (u, v) = (0.0, 0.0);
    let sp = blending_surface_info(nodes, edge, u, v);
    let (a1, a2, x0, a, b, deta) = (sp.a1, sp.a2, sp.x, sp.a, sp.b, sp.deta);
    let ca = [
        [1.0 / deta * a[1][1], -1.0 / deta * a[0][1]],
        [-1.0 / deta * a[1][0], 1.0 / deta * a[0][0]],
    ];
    let mut dual1 = [0.0f64; 3];
    let mut dual2 = [0.0f64; 3];
    for k in 0..3 {
        dual1[k] = ca[0][0] * a1[k] + ca[0][1] * a2[k];
        dual2[k] = ca[1][0] * a1[k] + ca[1][1] * a2[k];
    }
    // c = MATMUL( b, ContravA )
    let mut c = [[0.0f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            let mut s = 0.0;
            for k in 0..2 {
                s += b[i][k] * ca[k][j];
            }
            c[i][j] = s;
        }
    }
    let detc = c[0][0] * c[1][1] - c[0][1] * c[1][0];
    let trc = c[0][0] + c[1][1];
    let mut disc = trc.powi(2) - 4.0 * detc;
    if disc < 0.0 {
        disc = 0.0;
    }
    let mut umbilical = spherical;
    let (lambda1, lambda2) = if trc > 0.0 {
        (0.5 * (trc - disc.sqrt()), 0.5 * (trc + disc.sqrt()))
    } else {
        (0.5 * (trc + disc.sqrt()), 0.5 * (trc - disc.sqrt()))
    };
    let lambda_max = lambda1.abs().max(lambda2.abs());
    let mut planar = lambda1.abs() < EPS_SP && lambda2.abs() < EPS_SP;
    if plate && !planar {
        planar = true;
    }
    if !planar && !umbilical {
        let mut rk = 0.0f64;
        for (i, j) in [(2usize, 0usize), (3, 1)] {
            let s = ((nodes[i][0] - nodes[j][0]).powi(2)
                + (nodes[i][1] - nodes[j][1]).powi(2)
                + (nodes[i][2] - nodes[j][2]).powi(2))
            .sqrt();
            rk = rk.max(s);
        }
        umbilical = (lambda1 - lambda2).abs() < 0.5 * rk * lambda_max.powi(2);
    }
    let (g1, mut g2) = if planar || umbilical {
        let g1 = scale(1.0 / norm(&a1), a1);
        let c = dot(&a2, &g1);
        let g2 = [a2[0] - c * g1[0], a2[1] - c * g1[1], a2[2] - c * g1[2]];
        let g2 = scale(1.0 / norm(&g2), g2);
        (g1, g2)
    } else {
        let pdir = |lambda: f64| {
            let mut m = [
                [c[0][0] - lambda * 1.0, c[0][1] - lambda * 0.0],
                [c[1][0] - lambda * 0.0, c[1][1] - lambda * 1.0],
            ];
            let i = if (m[0][0].abs() + m[0][1].abs()) > (m[1][0].abs() + m[1][1].abs()) {
                0
            } else {
                1
            };
            let sc = m[i][0].abs() + m[i][1].abs();
            m[i] = [1.0 / sc * m[i][0], 1.0 / sc * m[i][1]];
            let p = [m[i][1], -m[i][0]];
            let nr = norm(&p);
            [1.0 / nr * p[0], 1.0 / nr * p[1]]
        };
        let p1 = pdir(lambda1);
        let p2 = pdir(lambda2);
        let mut g1 = [0.0f64; 3];
        let mut g2 = [0.0f64; 3];
        for k in 0..3 {
            g1[k] = p1[0] * dual1[k] + p1[1] * dual2[k];
            g2[k] = p2[0] * dual1[k] + p2[1] * dual2[k];
        }
        (scale(1.0 / norm(&g1), g1), scale(1.0 / norm(&g2), g2))
    };
    if dot(&d, &cross(g1, g2)) < 0.0 {
        g2 = scale(-1.0, g2);
    }
    let g3 = cross(g1, g2);

    // The Lagrange nodes of the geometry element in the frame.
    let (gu, gv) = geometry_nodes()?;
    let mut ln = [[0.0f64; 3]; MAX_PATCH_NODES];
    for j in 0..MAX_PATCH_NODES {
        let mut p = if j >= 4 {
            blending_surface_info(nodes, edge, gu[j], gv[j]).x
        } else {
            nodes[j]
        };
        for k in 0..3 {
            p[k] -= x0[k];
        }
        ln[j] = [dot(&p, &g1), dot(&p, &g2), dot(&p, &g3)];
    }

    // The Taylor polynomial from a 4x4 subgrid.
    let mut rk = f64::MAX;
    for (xi, eta) in [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)] {
        let gb = geometry_basis(xi, eta)?;
        let (mut p0, mut p1) = (0.0, 0.0);
        for j in 0..16 {
            p0 += ln[j][0] * gb[j];
        }
        for j in 0..16 {
            p1 += ln[j][1] * gb[j];
        }
        rk = rk.min(norm(&[p0, p1]));
    }
    let hk = rk / 4.0;
    let mut x1 = [0.0f64; 4];
    x1[0] = -hk / 2.0;
    x1[1] = x1[0] + hk / 3.0;
    x1[2] = x1[1] + hk / 3.0;
    x1[3] = hk / 2.0;
    let p0 = [dot(&x0, &g1), dot(&x0, &g2)];
    let mut z = [0.0f64; 16];
    let mut grid = 0;
    for j in 0..4 {
        let mut vk = -1.0 / 6.0 + j as f64 / 9.0;
        for i in 0..4 {
            let target = [x1[i], x1[j]];
            // `vk` carries over from the previous point of the row.
            let mut uk = -1.0 / 6.0 + i as f64 / 9.0;
            let mut converged = false;
            let mut pk = [0.0f64; 3];
            for _ in 0..GEOMETRY_MAX_ITERS {
                let s = blending_surface_info(nodes, edge, uk, vk);
                pk = s.x;
                let r = [p0[0] + target[0] - dot(&pk, &g1), p0[1] + target[1] - dot(&pk, &g2)];
                if norm(&r) < GEOMETRY_EPSILON * hk {
                    converged = true;
                    break;
                }
                let m = [
                    [dot(&s.a1, &g1), dot(&s.a2, &g1)],
                    [dot(&s.a1, &g2), dot(&s.a2, &g2)],
                ];
                let delta = solve2x2(m, r)?;
                uk += delta[0];
                vk += delta[1];
            }
            if !converged {
                return Err(ShellError::Fatal(
                    "LinesOfCurvatureFrame: Cannot create a subgrid: elements may not be shape-regular".into(),
                ));
            }
            if uk.abs() > 1.0 + 1.0e-3 || vk.abs() > 1.0 + 1.0e-3 {
                return Err(ShellError::Fatal(
                    "LinesOfCurvatureFrame: Cannot create a subgrid: a subgrid point outside element".into(),
                ));
            }
            let q = [pk[0] - x0[0], pk[1] - x0[1], pk[2] - x0[2]];
            z[grid] = dot(&q, &g3);
            grid += 1;
        }
    }
    let taylor = taylor_parameters(&z, hk, hk);
    let mut patch = [[0.0f64; 2]; MAX_PATCH_NODES];
    let umbilical_point = !planar && umbilical;
    if planar {
        for j in 0..MAX_PATCH_NODES {
            patch[j] = [ln[j][0], ln[j][1]];
        }
    } else {
        patch = curvature_patch(&ln, &taylor, umbilical_point)?;
    }
    Ok(Frame {
        e1: g1,
        e2: g2,
        e3: g3,
        o: x0,
        taylor,
        planar,
        umbilical: umbilical_point,
        patch,
    })
}

/// The Taylor coefficients `A`, `B`, `D`, `E`, `F`, `G` of the subgrid
/// heights `z` (LinesOfCurvatureFrame).
fn taylor_parameters(z0: &[f64; 16], h1: f64, h2: f64) -> [f64; 6] {
    // 1-based as in the Fortran.
    let z = |i: usize| z0[i - 1];
    let apar = (-9.0
        * (z(1) + 9.0 * z(10) + 9.0 * z(11) - 9.0 * z(12) + z(13) - z(14) - z(15) + z(16) - z(2) - z(3) + z(4)
            - 9.0 * z(5)
            + 9.0 * z(6)
            + 9.0 * z(7)
            - 9.0 * z(8)
            - 9.0 * z(9)))
        / (32.0 * h1.powi(2));
    let bpar = (-9.0
        * (z(1) + 9.0 * z(10) + 9.0 * z(11) - z(12) + z(13) - 9.0 * z(14) - 9.0 * z(15) + z(16) - 9.0 * z(2)
            - 9.0 * z(3)
            + z(4)
            - z(5)
            + 9.0 * z(6)
            + 9.0 * z(7)
            - z(8)
            - z(9)))
        / (32.0 * h2.powi(2));
    let dpar = (27.0
        * (z(1) + 27.0 * z(10) - 27.0 * z(11) + 9.0 * z(12) + z(13) - 3.0 * z(14) + 3.0 * z(15) - z(16)
            - 3.0 * z(2)
            + 3.0 * z(3)
            - z(4)
            - 9.0 * z(5)
            + 27.0 * z(6)
            - 27.0 * z(7)
            + 9.0 * z(8)
            - 9.0 * z(9)))
        / (16.0 * h1.powi(3));
    let epar = (9.0
        * (z(1) - 27.0 * z(10) - 27.0 * z(11) + 27.0 * z(12) - z(13) + z(14) + z(15) - z(16) - z(2) - z(3)
            + z(4)
            - 27.0 * z(5)
            + 27.0 * z(6)
            + 27.0 * z(7)
            - 27.0 * z(8)
            + 27.0 * z(9)))
        / (16.0 * h1.powi(2) * h2);
    let fpar = (9.0
        * (z(1) + 27.0 * z(10) - 27.0 * z(11) + z(12) + z(13) - 27.0 * z(14) + 27.0 * z(15) - z(16)
            - 27.0 * z(2)
            + 27.0 * z(3)
            - z(4)
            - z(5)
            + 27.0 * z(6)
            - 27.0 * z(7)
            + z(8)
            - z(9)))
        / (16.0 * h1 * h2.powi(2));
    let gpar = (27.0
        * (z(1) - 27.0 * z(10) - 27.0 * z(11) + 3.0 * z(12) - z(13) + 9.0 * z(14) + 9.0 * z(15) - z(16)
            - 9.0 * z(2)
            - 9.0 * z(3)
            + z(4)
            - 3.0 * z(5)
            + 27.0 * z(6)
            + 27.0 * z(7)
            - 3.0 * z(8)
            + 3.0 * z(9)))
        / (16.0 * h2.powi(3));
    [apar, bpar, dpar, epar, fpar, gpar]
}

/// The reparametrisation constants of `LinesOfCurvaturePatch` and
/// `SurfaceBasisVectors`: `c1, c2, c4, c5, c6, c7, b5, b7`.
fn patch_constants(t: &[f64; 6], umbilical: bool) -> [f64; 8] {
    let (a, b, e, f) = (t[0], t[1], t[3], t[4]);
    let parabolic = !umbilical && a.abs() / b.abs() < PARABOLIC_DELTA;
    let (c1, c2) = if umbilical { (0.0, 0.0) } else { (f / (b - a), e / (b - a)) };
    let c4 = if parabolic { 0.0 } else { -0.5 * (a.powi(2) + a * b + b.powi(2)) };
    let c5 = -2.0 * c1 * c2;
    let c6 = -a.powi(2) - a * b - c4;
    let c7 = -c5;
    let b7 = -b.powi(2) + c6;
    let b5 = -a * b - c6;
    [c1, c2, c4, c5, c6, c7, b5, b7]
}

/// `LinesOfCurvaturePatch` of a non-planar element: the patch nodes whose
/// reparametrised images are the Lagrange nodes.
fn curvature_patch(
    ln: &[[f64; 3]; MAX_PATCH_NODES],
    taylor: &[f64; 6],
    umbilical: bool,
) -> Result<[[f64; 2]; MAX_PATCH_NODES], ShellError> {
    let [c1, c2, c4, c5, c6, c7, b5, b7] = patch_constants(taylor, umbilical);
    let hk = ((ln[2][0] - ln[0][0]).powi(2) + (ln[2][1] - ln[0][1]).powi(2))
        .sqrt()
        .max(((ln[1][0] - ln[3][0]).powi(2) + (ln[1][1] - ln[3][1]).powi(2)).sqrt());
    let mut patch = [[0.0f64; 2]; MAX_PATCH_NODES];
    for i in 0..MAX_PATCH_NODES {
        let (mut y1, mut y2) = (ln[i][0], ln[i][1]);
        let mut converged = false;
        for _ in 0..GEOMETRY_MAX_ITERS {
            let x = y1 - c1 * y1.powi(2) / 2.0 + c2 * y1 * y2 + c1 * y2.powi(2) / 2.0 + c4 * y1.powi(3) / 6.0
                + c5 * y1.powi(2) * y2 / 2.0
                + c6 * y1 * y2.powi(2) / 2.0
                + c7 * y2.powi(3) / 6.0;
            let y = y2 - c2 * y1.powi(2) / 2.0 - c1 * y1 * y2 + c2 * y2.powi(2) / 2.0 - c5 * y1.powi(3) / 6.0
                + b5 * y1.powi(2) * y2 / 2.0
                - c7 * y1 * y2.powi(2) / 2.0
                + b7 * y2.powi(3) / 6.0;
            let r = [ln[i][0] - x, ln[i][1] - y];
            if norm(&r) < GEOMETRY_EPSILON * hk {
                converged = true;
                break;
            }
            let m = [
                [
                    1.0 - c1 * y1 + c2 * y2 + c4 * y1.powi(2) / 2.0 + c5 * y1 * y2 + c6 * y2.powi(2) / 2.0,
                    c2 * y1 + c1 * y2 + c5 * y1.powi(2) / 2.0 + c6 * y1 * y2 + c7 * y2.powi(2) / 2.0,
                ],
                [
                    -c2 * y1 - c1 * y2 - c5 * y1.powi(2) / 2.0 + b5 * y1 * y2 - c7 * y2.powi(2) / 2.0,
                    1.0 - c1 * y1 + c2 * y2 + b5 * y1.powi(2) / 2.0 - c7 * y1 * y2 + b7 * y2.powi(2) / 2.0,
                ],
            ];
            let delta = solve2x2(m, r)?;
            y1 += delta[0];
            y2 += delta[1];
        }
        if !converged {
            return Err(ShellError::Fatal("LinesOfCurvaturePatch: Nonlinear iteration fails".into()));
        }
        patch[i] = [y1, y2];
    }
    Ok(patch)
}

/// What `SurfaceBasisVectors` gives at a patch point.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Basis3 {
    pub a1: V3,
    pub a2: V3,
    pub a3: V3,
    pub a11: f64,
    pub a22: f64,
    pub sqrt_det_a: f64,
    pub b11: f64,
    pub b22: f64,
    pub c111: f64,
    pub c112: f64,
    pub c221: f64,
    pub c222: f64,
    pub c211: f64,
    pub c212: f64,
    pub dual1: V3,
    pub dual2: V3,
}

/// `SurfaceBasisVectors( x, y, TaylorParams, e1, e2, e3, o, ... )` with the
/// dual basis.
pub(crate) fn surface_basis(x: f64, y: f64, f: &Frame) -> Basis3 {
    if f.planar {
        return Basis3 {
            a1: f.e1,
            a2: f.e2,
            a3: f.e3,
            a11: 1.0,
            a22: 1.0,
            sqrt_det_a: 1.0,
            b11: 0.0,
            b22: 0.0,
            c111: 0.0,
            c112: 0.0,
            c221: 0.0,
            c222: 0.0,
            c211: 0.0,
            c212: 0.0,
            dual1: f.e1,
            dual2: f.e2,
        };
    }
    let t = &f.taylor;
    let (a, b, d, e, ff, g) = (t[0], t[1], t[2], t[3], t[4], t[5]);
    let [c1, c2, c4, c5, c6, c7, b5, b7] = patch_constants(t, f.umbilical);
    let k11 = 1.0 - c1 * x + c2 * y + c4 * x.powi(2) / 2.0 + c5 * x * y + c6 * y.powi(2) / 2.0;
    let k12 = -c2 * x - c1 * y - c5 * x.powi(2) / 2.0 + b5 * x * y - c7 * y.powi(2) / 2.0;
    let k13 = a * x
        + (-3.0 / 2.0 * a * c1 + d / 2.0) * x.powi(2)
        + (2.0 * a * c2 - b * c2 + e) * y * x
        + (a * c1 / 2.0 - b * c1 + ff / 2.0) * y.powi(2);
    let k21 = c2 * x + c1 * y + c5 * x.powi(2) / 2.0 + c6 * x * y + c7 * y.powi(2) / 2.0;
    let k22 = 1.0 - c1 * x + c2 * y + b5 * x.powi(2) / 2.0 - c7 * x * y + b7 * y.powi(2) / 2.0;
    let k23 = b * y
        + (a * c2 - b * c2 / 2.0 + e / 2.0) * x.powi(2)
        + (a * c1 - 2.0 * b * c1 + ff) * y * x
        + (3.0 / 2.0 * b * c2 + g / 2.0) * y.powi(2);
    let mut a1 = [0.0f64; 3];
    let mut a2 = [0.0f64; 3];
    for k in 0..3 {
        a1[k] = k11 * f.e1[k] + k12 * f.e2[k] + k13 * f.e3[k];
        a2[k] = k21 * f.e1[k] + k22 * f.e2[k] + k23 * f.e3[k];
    }
    let a3 = cross(a1, a2);
    let a3 = scale(1.0 / norm(&a3), a3);
    let a11 = dot(&a1, &a1);
    let a22 = dot(&a2, &a2);
    let sqrt_det_a = (a11 * a22).sqrt();
    let (b11, b22) = if f.umbilical {
        (a, b)
    } else {
        (
            a + (-2.0 * a * c1 + d) * x + (2.0 * a * c2 + e) * y,
            b + (-2.0 * b * c1 + ff) * x + (2.0 * b * c2 + g) * y,
        )
    };
    let d1 = a.powi(2) + c4 - c1.powi(2) + c2.powi(2);
    let d2 = a.powi(2) + a * b + c4 - c1.powi(2) + c2.powi(2);
    let c111 = -c1 + d1 * x;
    let c222 = c2 - d2 * y;
    Basis3 {
        a1,
        a2,
        a3,
        a11,
        a22,
        sqrt_det_a,
        b11,
        b22,
        c111,
        c112: -c222,
        c221: -c111,
        c222,
        c211: c222,
        c212: c111,
        dual1: scale(1.0 / a11, a1),
        dual2: scale(1.0 / a22, a2),
    }
}
