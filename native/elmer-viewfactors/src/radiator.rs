//! Radiator factors of point sources: the share of the power of an isotropic
//! point source that reaches each radiating surface element. A port of
//! `RadiatorFactors3D` of `fem/src/view3d` (`ViewFactors.c`,
//! `BiLinearUtil.c`, `TriangleUtil.c`, `RayTrace.c`) for triangles (303) and
//! bilinear quadrilaterals (404).
//!
//! Every element is integrated by adaptive bisection: an element whose
//! factor is below `factor_eps`, or whose area is below `area_eps`, gets
//! `n_rays` rays towards the source; the elements the rays pass through
//! unblocked count, the others are cut again (odd levels are always cut).
//! Blocking is tested against the radiating elements themselves. Line
//! radiators (`Radiator Coordinates(n, 6)`, a cylinder shading itself) are
//! not ported.

use crate::Error;

/// A radiating surface element: 3 or 4 corner points in Elmer node order and
/// the unit normal that points where the element receives radiation from
/// (`Normals` of the Fortran, after `DirectedNormalVector`).
#[derive(Clone, Debug)]
pub struct RadiatorSurface {
    pub nodes: Vec<[f64; 3]>,
    pub normal: [f64; 3],
}

/// The view factor keywords as `GetCartParameters` reads them, with the
/// names of the C side. Note that the Fortran passes `Viewfactor Area
/// Tolerance` as the C `Feps` and `Viewfactor Factor Tolerance` as `Aeps`:
/// [`RadiatorParams::from_keywords`] keeps that crossing.
#[derive(Copy, Clone, Debug)]
pub struct RadiatorParams {
    /// `FactorEPS` (C): the factor below which an element is not cut again.
    pub factor_eps: f64,
    /// `AreaEPS` (C): the area below which an element is not cut again.
    pub area_eps: f64,
    /// `RayEPS`: the ray parameter margin of the blocking test.
    pub ray_eps: f64,
    pub n_rays: usize,
    /// Gauss points: per direction of the quad rule (`Viewfactor Quad
    /// Integration Points`, a square), and of the triangle rule.
    pub quad_points: usize,
    pub tri_points: usize,
    /// `Viewfactor Line Integration Points`; also the 1D rule that measures
    /// the edges of a quad.
    pub line_points: usize,
}

impl RadiatorParams {
    /// From `Viewfactor Area Tolerance` (default 0.1), `Viewfactor Factor
    /// Tolerance` (0.01), `Viewfactor Number Of Rays` (1) and the point
    /// counts, as `GetCartParameters` and the call of `RadiatorFactors3D`
    /// hand them over. The ray tolerance keyword is misspelt in the Fortran
    /// and never found, so it is always 1e-5.
    pub fn from_keywords(
        area_tol: Option<f64>,
        factor_tol: Option<f64>,
        n_rays: Option<usize>,
        line_points: Option<usize>,
        quad_points: Option<usize>,
        tri_points: Option<usize>,
    ) -> Self {
        let fortran_area_eps = area_tol.unwrap_or(1.0e-1);
        let fortran_fact_eps = factor_tol.unwrap_or(1.0e-2);
        Self {
            factor_eps: fortran_area_eps,
            area_eps: fortran_fact_eps,
            ray_eps: 1.0e-5,
            n_rays: n_rays.unwrap_or(1),
            quad_points: quad_points.unwrap_or(4),
            tri_points: tri_points.unwrap_or(3),
            line_points: line_points.unwrap_or(2),
        }
    }
}

/// `IPoints1d`: the 1D Gauss points on [0, 1] with their 14-digit values.
fn gauss_1d(n: usize) -> Result<(Vec<f64>, Vec<f64>), Error> {
    Ok(match n {
        1 => (vec![0.5], vec![1.0]),
        2 => (vec![0.78867513459481, 0.21132486540519], vec![0.5, 0.5]),
        3 => (
            vec![0.88729833462074, 0.11270166537926, 0.5],
            vec![0.27777777777778, 0.27777777777778, 0.44444444444444],
        ),
        4 => (
            vec![
                0.66999052179243,
                0.33000947820757,
                0.93056815579703,
                0.06943184420297,
            ],
            vec![
                0.32607257743127,
                0.32607257743127,
                0.17392742256873,
                0.17392742256873,
            ],
        ),
        _ => {
            return Err(Error::Invalid(format!(
                "{n} Gauss points per direction are not tabulated"
            )))
        }
    })
}

/// The triangle rule of `FillIPointArrays`.
/// `(u, v, weight)` columns of the triangle rule.
type TriRule = (Vec<f64>, Vec<f64>, Vec<f64>);

fn gauss_tri(n: usize) -> Result<TriRule, Error> {
    Ok(match n {
        1 => (vec![0.33333333333333], vec![0.33333333333333], vec![0.5]),
        3 => (
            vec![0.16666666666667, 0.66666666666667, 0.16666666666667],
            vec![0.16666666666667, 0.16666666666667, 0.66666666666667],
            vec![0.16666666666667; 3],
        ),
        _ => {
            return Err(Error::Invalid(format!(
                "{n} triangle Gauss points are not tabulated"
            )))
        }
    })
}

/// `PolyFactors`: coefficients of the coordinates (rows 0-2) and of the
/// normal (rows 3-5) in `1, u, v, uv` (quads) or `1, u, v` (triangles).
type Poly = [[f64; 4]; 6];

#[derive(Copy, Clone, PartialEq, Eq)]
enum Kind {
    Tri,
    Quad,
}

fn poly_of(s: &RadiatorSurface) -> Result<(Kind, Poly), Error> {
    let mut p: Poly = [[0.0; 4]; 6];
    let kind = match s.nodes.len() {
        3 => {
            // ShapeFunctionMatrix3: node k, coefficient j.
            const M: [[f64; 3]; 3] = [[1., -1., -1.], [0., 1., 0.], [0., 0., 1.]];
            for j in 0..3 {
                for k in 0..3 {
                    for n in 0..3 {
                        p[n][j] += M[k][j] * s.nodes[k][n];
                        p[n + 3][j] += M[k][j] * s.normal[n];
                    }
                }
            }
            Kind::Tri
        }
        4 => {
            const M: [[f64; 4]; 4] = [
                [1., -1., -1., 1.],
                [0., 1., 0., -1.],
                [0., 0., 0., 1.],
                [0., 0., 1., -1.],
            ];
            for j in 0..4 {
                for k in 0..4 {
                    for n in 0..3 {
                        p[n][j] += M[k][j] * s.nodes[k][n];
                        p[n + 3][j] += M[k][j] * s.normal[n];
                    }
                }
            }
            Kind::Quad
        }
        n => return Err(Error::Invalid(format!("a radiating element of {n} nodes"))),
    };
    Ok((kind, p))
}

fn bilin(u: f64, v: f64, n: &[f64; 4]) -> f64 {
    (n[3] * u + n[2]) * v + n[1] * u + n[0]
}

fn tri(u: f64, v: f64, n: &[f64; 4]) -> f64 {
    n[2] * v + n[1] * u + n[0]
}

/// A node of the subdivision tree of one element, rebuilt for every source.
struct Node {
    kind: Kind,
    poly: Poly,
    area: f64,
}

struct Rules {
    u1: Vec<f64>,
    s1: Vec<f64>,
    uq: Vec<f64>,
    vq: Vec<f64>,
    sq: Vec<f64>,
    ut: Vec<f64>,
    vt: Vec<f64>,
    st: Vec<f64>,
}

impl Rules {
    fn new(p: &RadiatorParams) -> Result<Self, Error> {
        // The quad rule is the product of a square-root-of-N point 1D rule.
        let n1d = ((p.quad_points as f64).sqrt() + 0.5) as usize;
        let (uq1, sq1) = gauss_1d(n1d)?;
        let (mut uq, mut vq, mut sq) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n1d {
            for j in 0..n1d {
                uq.push(uq1[i]);
                vq.push(uq1[j]);
                sq.push(sq1[i] * sq1[j]);
            }
        }
        let (u1, s1) = gauss_1d(p.line_points)?;
        let (ut, vt, st) = gauss_tri(p.tri_points)?;
        Ok(Self {
            u1,
            s1,
            uq,
            vq,
            sq,
            ut,
            vt,
            st,
        })
    }
}

fn eof_a(u: f64, v: f64, p: &Poly, kind: Kind) -> f64 {
    let d = |n: &[f64; 4]| match kind {
        Kind::Quad => (n[3] * v + n[1], n[3] * u + n[2]),
        Kind::Tri => (n[1], n[2]),
    };
    let (dxu, dxv) = d(&p[0]);
    let (dyu, dyv) = d(&p[1]);
    let (dzu, dzv) = d(&p[2]);
    let auu = dxu * dxu + dyu * dyu + dzu * dzu;
    let auv = dxu * dxv + dyu * dyv + dzu * dzv;
    let avv = dxv * dxv + dyv * dyv + dzv * dzv;
    (auu * avv - auv * auv).sqrt()
}

fn eof_l(u: f64, v: f64, p: &Poly, u_not_v: bool) -> f64 {
    let (dxu, dxv) = (p[0][3] * v + p[0][1], p[0][3] * u + p[0][2]);
    let (dyu, dyv) = (p[1][3] * v + p[1][1], p[1][3] * u + p[1][2]);
    let (dzu, dzv) = (p[2][3] * v + p[2][1], p[2][3] * u + p[2][2]);
    if u_not_v {
        (dxu * dxu + dyu * dyu + dzu * dzu).sqrt()
    } else {
        (dxv * dxv + dyv * dyv + dzv * dzv).sqrt()
    }
}

fn area_of(node_kind: Kind, p: &Poly, r: &Rules) -> f64 {
    match node_kind {
        Kind::Quad => (0..r.uq.len())
            .map(|i| r.sq[i] * eof_a(r.uq[i], r.vq[i], p, Kind::Quad))
            .sum(),
        Kind::Tri => (0..r.ut.len())
            .map(|i| r.st[i] * eof_a(r.ut[i], r.vt[i], p, Kind::Tri))
            .sum(),
    }
}

fn quad_length(p: &Poly, u_not_v: bool, r: &Rules) -> f64 {
    let mut len = 0.0;
    for j in 0..2 {
        for i in 0..r.u1.len() {
            let e = if u_not_v {
                eof_l(r.u1[i], j as f64, p, true)
            } else {
                eof_l(j as f64, r.u1[i], p, false)
            };
            len += r.s1[i] * e;
        }
    }
    len
}

const SHAPE4: [[f64; 4]; 4] = [
    [1., -1., -1., 1.],
    [0., 1., 0., -1.],
    [0., 0., 0., 1.],
    [0., 0., 1., -1.],
];

/// `BiLinearSubdivide`: cut in the longer of the parameter directions.
fn subdivide_quad(p: &Poly, r: &Rules) -> (Poly, Poly) {
    let half_u = quad_length(p, true, r) > quad_length(p, false, r);
    type Corners = [(f64, f64); 4];
    let (lc, rc): (Corners, Corners) = if half_u {
        (
            [(0.0, 0.0), (0.5, 0.0), (0.5, 1.0), (0.0, 1.0)],
            [(0.5, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.0)],
        )
    } else {
        (
            [(0.0, 0.0), (1.0, 0.0), (1.0, 0.5), (0.0, 0.5)],
            [(0.0, 0.5), (1.0, 0.5), (1.0, 1.0), (0.0, 1.0)],
        )
    };
    let (mut l, mut rt): (Poly, Poly) = ([[0.0; 4]; 6], [[0.0; 4]; 6]);
    for row in 0..6 {
        let xl: Vec<f64> = lc.iter().map(|&(u, v)| bilin(u, v, &p[row])).collect();
        let xr: Vec<f64> = rc.iter().map(|&(u, v)| bilin(u, v, &p[row])).collect();
        for i in 0..4 {
            for j in 0..4 {
                l[row][i] += SHAPE4[j][i] * xl[j];
                rt[row][i] += SHAPE4[j][i] * xr[j];
            }
        }
    }
    (l, rt)
}

/// `TriangleSubdivide`: cut the longest edge.
fn subdivide_tri(p: &Poly) -> (Poly, Poly) {
    let at =
        |u: f64, v: f64| -> [f64; 3] { [tri(u, v, &p[0]), tri(u, v, &p[1]), tri(u, v, &p[2])] };
    let (a0, a1, a2) = (at(0.0, 0.0), at(1.0, 0.0), at(0.0, 1.0));
    let d2 = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| (b[i] - a[i]).powi(2)).sum::<f64>();
    let (l1, l2, l3) = (d2(a0, a1), d2(a0, a2), d2(a1, a2));
    let mut l: Poly = [[0.0; 4]; 6];
    let mut r: Poly = [[0.0; 4]; 6];
    let sub = |a: [f64; 3], b: [f64; 3]| [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let set = |poly: &mut Poly, o: [f64; 3], e1: [f64; 3], e2: [f64; 3]| {
        for n in 0..3 {
            poly[n][0] = o[n];
            poly[n][1] = e1[n];
            poly[n][2] = e2[n];
        }
    };
    if l1 > l2 && l1 > l3 {
        let n = at(0.5, 0.0);
        set(&mut l, a0, sub(a0, n), sub(a0, a2));
        set(&mut r, n, sub(n, a1), sub(n, a2));
    } else if l2 > l3 {
        let n = at(0.0, 0.5);
        set(&mut l, a0, sub(a0, a1), sub(a0, n));
        set(&mut r, n, sub(n, a1), sub(n, a2));
    } else {
        let n = at(0.5, 0.5);
        set(&mut l, a0, sub(a0, a1), sub(a0, n));
        set(&mut r, a0, sub(a0, n), sub(a0, a2));
    }
    l[3..].copy_from_slice(&p[3..]);
    r[3..].copy_from_slice(&p[3..]);
    (l, r)
}

/// The PCG32 stream of `vrand()` (thread 0, rank 0).
struct Pcg {
    state: u64,
    inc: u64,
}

impl Pcg {
    fn new() -> Self {
        Self {
            state: 0x853c49e6748fea9b,
            inc: 0xda3e39cb94b95bdb,
        }
    }

    fn next(&mut self) -> f64 {
        let old = self.state;
        self.state = old
            .wrapping_mul(6364136223846793005)
            .wrapping_add(self.inc | 1);
        let x = (((old >> 18) ^ old) >> 27) as u32;
        let r = (old >> 59) as u32;
        let out = (x >> r) | (x << (r.wrapping_neg() & 31));
        out as f64 / u32::MAX as f64
    }
}

/// `RayHitTriangle`.
fn ray_hit_tri(p: &Poly, f: [f64; 3], d: [f64; 3], eps: f64) -> bool {
    let (x, y, z) = (&p[0], &p[1], &p[2]);
    let (a11, a12, a13) = (x[1], x[2], -d[0]);
    let (a21, a22, a23) = (y[1], y[2], -d[1]);
    let (a31, a32, a33) = (z[1], z[2], -d[2]);
    let mut det = a11 * (a22 * a33 - a32 * a23)
        + a12 * (a31 * a23 - a21 * a33)
        + a13 * (a21 * a32 - a31 * a22);
    if det.abs() < 1.0e-9 {
        return false;
    }
    det = 1.0 / det;
    let (px, py, pz) = (f[0] - x[0], f[1] - y[0], f[2] - z[0]);
    let (b31, b32, b33) = (
        a21 * a32 - a22 * a31,
        a12 * a31 - a11 * a32,
        a11 * a22 - a21 * a12,
    );
    let t = det * (b31 * px + b32 * py + b33 * pz);
    if t < eps || t > 1.0 - eps {
        return false;
    }
    let (b11, b12, b13) = (
        a22 * a33 - a23 * a32,
        a32 * a13 - a12 * a33,
        a12 * a23 - a22 * a13,
    );
    let u = det * (b11 * px + b12 * py + b13 * pz);
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let (b21, b22, b23) = (
        a23 * a31 - a21 * a33,
        a11 * a33 - a31 * a13,
        a21 * a13 - a11 * a23,
    );
    let v = det * (b21 * px + b22 * py + b23 * pz);
    if !(0.0..=1.0).contains(&v) {
        return false;
    }
    u + v <= 1.0
}

/// `RayHitBiLinear`.
fn ray_hit_quad(p: &Poly, f: [f64; 3], d: [f64; 3], eps: f64) -> bool {
    let (x, y, z) = (&p[0], &p[1], &p[2]);
    let (dx, dy, dz) = (d[0], d[1], d[2]);
    let a0 = -dy * x[1] + dx * y[1];
    let b0 = -dy * x[2] + dx * y[2];
    let c0 = -dy * x[3] + dx * y[3];
    let d0 = -dy * (f[0] - x[0]) + dx * (f[1] - y[0]);
    let a1 = -dy * z[1] + dz * y[1];
    let b1 = -dy * z[2] + dz * y[2];
    let c1 = -dy * z[3] + dz * y[3];
    let d1 = -dy * (f[2] - z[0]) + dz * (f[1] - y[0]);

    let (td, tf, tx) = if dx.abs() > dy.abs() && dx.abs() > dz.abs() {
        (dx, f[0], x)
    } else if dy.abs() > dz.abs() {
        (dy, f[1], y)
    } else {
        (dz, f[2], z)
    };
    let hit_at = |u: f64, v: f64| -> bool {
        let t = (bilin(u, v, tx) - tf) / td;
        t > eps && t < 1.0 - eps
    };
    let unit = |s: f64| (0.0..=1.0).contains(&s);

    if c0.abs() < 1.0e-12 && c1.abs() < 1.0e-12 {
        let a = a0 * b1 - a1 * b0;
        if a.abs() < 1.0e-12 {
            return false;
        }
        let a = 1.0 / a;
        let u = a * (b1 * d0 - b0 * d1);
        if unit(u) {
            let v = a * (a0 * d1 - a1 * d0);
            if unit(v) {
                return hit_at(u, v);
            }
        }
        return false;
    }
    if c0.abs() < 1.0e-12 && a0.abs() < 1.0e-12 {
        let v = d0 / b0;
        if unit(v) {
            let u = (d1 - b1 * v) / (a1 + c1 * v);
            if unit(u) {
                return hit_at(u, v);
            }
        }
        return false;
    }
    if c1.abs() < 1.0e-12 && a1.abs() < 1.0e-12 {
        let v = d1 / b1;
        if unit(v) {
            let u = (d0 - b0 * v) / (a0 + c0 * v);
            if unit(u) {
                return hit_at(u, v);
            }
        }
        return false;
    }
    if c0.abs() < 1.0e-12 && b0.abs() < 1.0e-12 {
        let u = d0 / a0;
        if unit(u) {
            let v = (d1 - a1 * u) / (b1 + c1 * u);
            if unit(v) {
                return hit_at(u, v);
            }
        }
        return false;
    }
    if c1.abs() < 1.0e-12 && b1.abs() < 1.0e-12 {
        let u = d1 / a1;
        if unit(u) {
            let v = (d0 - a0 * u) / (b0 + c0 * u);
            if unit(v) {
                return hit_at(u, v);
            }
        }
        return false;
    }

    let a = c0 * b1 - b0 * c1;
    let b = d0 * c1 - b0 * a1 + a0 * b1 - c0 * d1;
    let c = d0 * a1 - a0 * d1;
    let solve_u = |v: f64| -> f64 {
        if (a0 + c0 * v).abs() < 1.0e-12 {
            (d1 - b1 * v) / (a1 + c1 * v)
        } else {
            (d0 - b0 * v) / (a0 + c0 * v)
        }
    };

    if a.abs() < 1.0e-12 {
        if b.abs() < 1.0e-12 {
            return false;
        }
        let v = -c / b;
        if !unit(v) {
            return false;
        }
        let u = solve_u(v);
        if !unit(u) {
            return false;
        }
        return hit_at(u, v);
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return false;
    }
    let disc = disc.sqrt();
    let v1 = if b > 0.0 {
        -2.0 * c / (b + disc)
    } else {
        (-b + disc) / (2.0 * a)
    };
    if unit(v1) {
        let u = solve_u(v1);
        if unit(u) && hit_at(u, v1) {
            return true;
        }
    }
    let v2 = if b > 0.0 {
        -(b + disc) / (2.0 * a)
    } else {
        2.0 * c / (-b + disc)
    };
    if unit(v2) {
        let u = solve_u(v2);
        if unit(u) && hit_at(u, v2) {
            return true;
        }
    }
    false
}

struct Ctx<'a> {
    params: RadiatorParams,
    rules: &'a Rules,
    elements: &'a [(Kind, Poly)],
    rng: &'a mut Pcg,
}

impl Ctx<'_> {
    fn blocked(&self, f: [f64; 3], d: [f64; 3]) -> bool {
        self.elements.iter().any(|(k, p)| match k {
            Kind::Tri => ray_hit_tri(p, f, d, self.params.ray_eps),
            Kind::Quad => ray_hit_quad(p, f, d, self.params.ray_eps),
        })
    }

    /// `Linear/BiLinear/TriangleIntegrateDiffToArea` for a point source:
    /// `integral cosB / R^2 dA` over the node's element.
    fn diff_to_area(&self, n: &Node, src: [f64; 3]) -> f64 {
        let r = self.rules;
        let p = &n.poly;
        match n.kind {
            Kind::Quad => {
                let mut f = 0.0;
                for i in 0..r.uq.len() {
                    let (u, v) = (r.uq[i], r.vq[i]);
                    let mut nt = [bilin(u, v, &p[3]), bilin(u, v, &p[4]), bilin(u, v, &p[5])];
                    let mut m = nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2];
                    if (1.0 - m).abs() > 1.0e-8 {
                        m = 1.0 / m.sqrt();
                        nt.iter_mut().for_each(|c| *c *= m);
                    }
                    let d = [
                        bilin(u, v, &p[0]) - src[0],
                        bilin(u, v, &p[1]) - src[1],
                        bilin(u, v, &p[2]) - src[2],
                    ];
                    let rr = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    let cos_b = (-d[0] * nt[0] - d[1] * nt[1] - d[2] * nt[2]) / rr;
                    if cos_b < 1.0e-8 {
                        continue;
                    }
                    f += eof_a(u, v, p, Kind::Quad) * cos_b * r.sq[i] / (rr * rr);
                }
                f
            }
            Kind::Tri => {
                let (u, v) = (1.0 / 3.0, 1.0 / 3.0);
                let mut nt = [tri(u, v, &p[3]), tri(u, v, &p[4]), tri(u, v, &p[5])];
                let m = 1.0 / (nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2]).sqrt();
                nt.iter_mut().for_each(|c| *c *= m);
                let mut f = 0.0;
                for i in 0..r.ut.len() {
                    let (u, v) = (r.ut[i], r.vt[i]);
                    let d = [
                        tri(u, v, &p[0]) - src[0],
                        tri(u, v, &p[1]) - src[1],
                        tri(u, v, &p[2]) - src[2],
                    ];
                    let rr = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    let cos_b = (-d[0] * nt[0] - d[1] * nt[1] - d[2] * nt[2]) / rr;
                    if cos_b < 1.0e-8 {
                        continue;
                    }
                    f += 2.0 * n.area * cos_b * r.st[i] / (rr * rr);
                }
                f
            }
        }
    }

    /// `*ComputeRadiatorFactors` and `ComputeViewFactorValue` in one: the
    /// sum over the leaves of `area * view factor` below `node`.
    fn factor(&mut self, node: &Node, level: usize, src: [f64; 3]) -> f64 {
        let pr = self.params;
        let mut subdivide = level & 1 == 1;
        if !subdivide {
            let fa = self.diff_to_area(node, src);
            if fa < 1.0e-10 {
                return 0.0;
            }
            if fa < pr.factor_eps || node.area < pr.area_eps {
                let n_rays = pr.n_rays as f64;
                let mut hit = n_rays;
                for _ in 0..pr.n_rays {
                    let (fx, d) = self.sample_ray(node, src);
                    if self.blocked(fx, d) {
                        hit -= 1.0;
                    }
                }
                if hit == 0.0 {
                    return 0.0;
                }
                if hit == n_rays || fa < pr.factor_eps / 2.0 || node.area < pr.area_eps / 2.0 {
                    let f = hit * fa / n_rays;
                    let vf = f / node.area / (4.0 * std::f64::consts::PI);
                    return node.area * vf;
                }
            }
            subdivide = true;
        }
        debug_assert!(subdivide);
        let (lp, rp) = match node.kind {
            Kind::Quad => subdivide_quad(&node.poly, self.rules),
            Kind::Tri => subdivide_tri(&node.poly),
        };
        let left = Node {
            kind: node.kind,
            area: area_of(node.kind, &lp, self.rules),
            poly: lp,
        };
        let right = Node {
            kind: node.kind,
            area: area_of(node.kind, &rp, self.rules),
            poly: rp,
        };
        self.factor(&left, level + 1, src) + self.factor(&right, level + 1, src)
    }

    /// A random point of the node's element and the ray from it to the
    /// source.
    fn sample_ray(&mut self, node: &Node, src: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let p = &node.poly;
        let (mut u, mut v) = (self.rng.next(), self.rng.next());
        let f = match node.kind {
            Kind::Tri => {
                while u + v > 1.0 {
                    u = 1.0 - u;
                    v = 1.0 - v;
                }
                [tri(u, v, &p[0]), tri(u, v, &p[1]), tri(u, v, &p[2])]
            }
            Kind::Quad => [bilin(u, v, &p[0]), bilin(u, v, &p[1]), bilin(u, v, &p[2])],
        };
        (f, [src[0] - f[0], src[1] - f[1], src[2] - f[2]])
    }
}

/// The raw radiator factors, row-major `radiators.len() * surfaces.len()`:
/// for each source the fraction of its power the element receives, before
/// the normalisation to one that `ViewFactors` applies to a closed set.
pub fn radiator_factors(
    surfaces: &[RadiatorSurface],
    radiators: &[[f64; 3]],
    params: &RadiatorParams,
) -> Result<Vec<f64>, Error> {
    let rules = Rules::new(params)?;
    let elements: Vec<(Kind, Poly)> = surfaces.iter().map(poly_of).collect::<Result<_, _>>()?;
    let roots: Vec<Node> = elements
        .iter()
        .map(|(k, p)| Node {
            kind: *k,
            poly: *p,
            area: area_of(*k, p, &rules),
        })
        .collect();
    let mut rng = Pcg::new();
    let mut ctx = Ctx {
        params: *params,
        rules: &rules,
        elements: &elements,
        rng: &mut rng,
    };
    let n = surfaces.len();
    let mut out = vec![0.0; radiators.len() * n];
    for (i, src) in radiators.iter().enumerate() {
        for (j, root) in roots.iter().enumerate() {
            if root.area < 1.0e-10 {
                continue;
            }
            out[i * n + j] = ctx.factor(root, 0, *src);
        }
    }
    Ok(out)
}

/// The areas `ElementArea` gives the elements (the rule of the Fortran
/// would integrate the same planar patch exactly).
pub fn surface_areas(
    surfaces: &[RadiatorSurface],
    params: &RadiatorParams,
) -> Result<Vec<f64>, Error> {
    let rules = Rules::new(params)?;
    surfaces
        .iter()
        .map(|s| poly_of(s).map(|(k, p)| area_of(k, &p, &rules)))
        .collect()
}
