//! The integration rules `GaussPoints` (Integration.F90:2882) picks for
//! active p elements that `elmer-fe` does not provide:
//!
//! * `GaussPointsQuad(np, PMethod=.TRUE.)` (Integration.F90:1850): the
//!   economic quadrilateral rules for serendipity p elements, 8 to 60
//!   points (tables `UPQuad*`, `VPQuad*`, `SPQuad*`, Integration.F90:425);
//! * triangles on the p reference element (vertices (-1, 0), (1, 0) and
//!   (0, sqrt 3)): `TriangleSimplexRulePoints` (Integration.F90:1653) picks
//!   a tabulated rule, mapped by `GaussPointsTriangle(n,
//!   PReferenceElement=.TRUE.)` (Integration.F90:1742), and otherwise
//!   `GaussPointsPTriangle` (Integration.F90:1689) collapses a quad rule;
//! * wedges on the p reference element: `WedgeRulePoints`
//!   (Integration.F90:2343), `GaussPointsWedge2` (triangle x line,
//!   Integration.F90:2406), `GaussPointsWedgeEconomic`
//!   (Integration.F90:2610) and `GaussPointsPWedge` (collapsed brick,
//!   Integration.F90:2242), dispatched as `GaussPoints` does.

use elmer_fe::quadrature::{
    gauss_points_1d, gauss_points_brick, gauss_points_quad, gauss_points_tetra,
    gauss_points_triangle, gauss_points_wedge_economic, GaussRule,
};

use crate::Error;

fn fe(r: Result<GaussRule, elmer_fe::FeError>) -> Result<GaussRule, Error> {
    r.map_err(|err| Error::Unsupported(err.to_string()))
}

fn table(u: &[f64], v: &[f64], s: &[f64]) -> GaussRule {
    GaussRule {
        u: u.to_vec(),
        v: v.to_vec(),
        w: vec![0.0; u.len()],
        s: s.to_vec(),
    }
}

/// `GaussPointsQuad(np, PMethod)`: with `economic` and `4 < np <= 60` one
/// of the tabulated rules (other counts in that range are Elmer's
/// `Fatal`), else the tensor rule.
pub fn gauss_points_quad_p(np: usize, economic: bool) -> Result<GaussRule, Error> {
    if !(economic && np > 4 && np <= 60) {
        return fe(gauss_points_quad(np));
    }
    Ok(match np {
        8 => table(&U_PQUAD_8, &V_PQUAD_8, &S_PQUAD_8),
        12 => table(&U_PQUAD_12, &V_PQUAD_12, &S_PQUAD_12),
        20 => table(&U_PQUAD_20, &V_PQUAD_20, &S_PQUAD_20),
        25 => table(&U_PQUAD_25, &V_PQUAD_25, &S_PQUAD_25),
        36 => table(&U_PQUAD_36, &V_PQUAD_36, &S_PQUAD_36),
        45 => table(&U_PQUAD_45, &V_PQUAD_45, &S_PQUAD_45),
        60 => table(&U_PQUAD_60, &V_PQUAD_60, &S_PQUAD_60),
        _ => {
            return Err(Error::BadDefinition(format!(
                "GaussPointsQuad: Invalid number of points for p-quadrature: {np}"
            )))
        }
    })
}

/// `TriangleSimplexRulePoints(np)`: the tabulated simplex rule exact for
/// the degree a `np`-point tensor rule integrates, 0 when none is.
pub fn triangle_simplex_rule_points(np: usize) -> usize {
    let maxp = (np as f64).sqrt().round() as i64;
    let deg = 2 * (maxp - 1).max(0) + 1;
    match deg {
        i64::MIN..=1 => 1,
        2 => 3,
        3 => 4,
        4 => 6,
        5 => 7,
        6 => 11,
        7 => 12,
        8 => 17,
        9 => 20,
        _ => 0,
    }
}

/// `WedgeRulePoints(np)`: the tabulated wedge rule for the degree an
/// `np`-point tensor rule integrates (a triangle rule times a line rule, or
/// a smaller economic rule), 0 when none is.
pub fn wedge_rule_points(np: usize) -> usize {
    let maxp = (np as f64).powf(1.0 / 3.0).round() as i64;
    let deg = 2 * (maxp - 1).max(0) + 1;
    let mtri = match deg {
        i64::MIN..=1 => 1,
        2 => 3,
        3 => 4,
        4 => 6,
        5 => 7,
        6 => 11,
        7 => 12,
        8 => 17,
        9 => 20,
        _ => 0,
    };
    if mtri == 0 {
        return 0;
    }
    let nseg = ((deg + 1) / 2) as usize;
    let m = mtri * nseg;
    let mecon = match deg {
        i64::MIN..=2 => 4,
        3 => 7,
        4 => 11,
        5 => 15,
        _ => 0,
    };
    if mecon > 0 {
        m.min(mecon)
    } else {
        m
    }
}

/// `GaussPointsPWedge(n)`: the `n`-point brick rule collapsed onto the p
/// reference wedge.
pub fn gauss_points_p_wedge(n: usize) -> Result<GaussRule, Error> {
    let mut p = fe(gauss_points_brick(n))?;
    let sqrt3 = 3.0f64.sqrt();
    for i in 0..p.len() {
        let (uh, vh, sh) = (p.u[i], p.v[i], p.s[i]);
        p.u[i] = 1.0 / 2.0 * (uh - uh * vh);
        p.v[i] = sqrt3 / 2.0 * (1.0 + vh);
        p.s[i] = sh * sqrt3 * (1.0 - vh) / 4.0;
    }
    Ok(p)
}

/// `GaussPointsWedge2(m, n, PReferenceElement=.TRUE.)`: an `m`-point
/// triangle rule times an `n`-point line rule on the p reference wedge, or
/// the collapsed `n^3` brick rule when no `m`-point triangle rule exists.
pub fn gauss_points_wedge2_p(m: usize, n: usize) -> Result<GaussRule, Error> {
    if !triangle_rule_tabulated(m) {
        return gauss_points_p_wedge(n * n * n);
    }
    // The raw triangle tables: elmer-fe halves the weights of the rules up
    // to 7 points (exactly, a power of two).
    let t = fe(gauss_points_triangle(m))?;
    let line = fe(gauss_points_1d(n))?;
    let sqrt3 = 3.0f64.sqrt();
    let mut p = GaussRule::default();
    for i in 0..m {
        let (u, v) = (t.u[i], t.v[i]);
        let uq = -1.0 + 2.0 * u + v;
        let vq = sqrt3 * v;
        let sq = if m <= 7 {
            sqrt3 * (t.s[i] * 2.0)
        } else {
            sqrt3 * 2.0 * t.s[i]
        };
        for j in 0..n {
            p.u.push(uq);
            p.v.push(vq);
            p.w.push(line.u[j]);
            p.s.push(sq * line.s[j]);
        }
    }
    Ok(p)
}

/// `GaussPointsWedgeEconomic(n, PReferenceElement=.TRUE.)`.
pub fn gauss_points_wedge_economic_p(n: usize) -> Result<GaussRule, Error> {
    // elmer-fe maps the raw tables with (x + 1) / 2 and s / 4 (exact).
    let mut p = fe(gauss_points_wedge_economic(n))?;
    let sqrt3 = 3.0f64.sqrt();
    for i in 0..p.len() {
        let (uq, vq, sq) = (p.u[i], p.v[i], p.s[i] * 4.0);
        p.u[i] = -1.0 + 2.0 * uq + vq;
        p.v[i] = sqrt3 * vq;
        p.s[i] = sqrt3 * sq / 2.0;
    }
    Ok(p)
}

/// The rule `GaussPoints` gives a p wedge with `n` points: `np_given` is
/// `GaussPoints(Element, np)` (the tabulated rules by `n` itself),
/// otherwise `n = PDefs % GaussPoints` goes through `WedgeRulePoints`.
pub fn p_wedge_rule(n: usize, np_given: bool) -> Result<GaussRule, Error> {
    let mut n = n;
    let mut tabulated = np_given;
    if !tabulated {
        let m = wedge_rule_points(n);
        if m > 0 {
            n = m;
            tabulated = true;
        }
    }
    if tabulated {
        let nseg = match n {
            1..=3 => n,
            6 | 8 => 2,
            12 | 18 | 21 => 3,
            28 | 44 | 48 => 4,
            85 | 100 => 5,
            _ => 0,
        };
        if let Some(ntri) = n.checked_div(nseg) {
            return gauss_points_wedge2_p(ntri, nseg);
        }
        if [4, 5, 7, 10, 11, 14, 15, 16, 24].contains(&n) {
            return gauss_points_wedge_economic_p(n);
        }
    }
    gauss_points_p_wedge(n)
}

/// `GaussPointsPBrick(nx, ny, nz)`: the tensor rule with `u` outermost.
fn gauss_points_p_brick(nx: usize, ny: usize, nz: usize) -> Result<GaussRule, Error> {
    let (x, y, z) = (
        fe(gauss_points_1d(nx))?,
        fe(gauss_points_1d(ny))?,
        fe(gauss_points_1d(nz))?,
    );
    let mut p = GaussRule::default();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                p.u.push(x.u[i]);
                p.v.push(y.u[j]);
                p.w.push(z.u[k]);
                p.s.push(x.s[i] * y.s[j] * z.s[k]);
            }
        }
    }
    Ok(p)
}

/// `n = DBLE(np)**(1.0D0/3.0D0) + 0.5D0`, truncated.
fn cube_root_points(np: usize) -> usize {
    ((np as f64).powf(1.0 / 3.0) + 0.5) as usize
}

/// `TetraSimplexRulePoints(np)`: the tabulated tetra rule for the degree an
/// `np`-point tensor rule integrates, 0 when none is.
pub fn tetra_simplex_rule_points(np: usize) -> usize {
    let maxp = (np as f64).powf(1.0 / 3.0).round() as i64;
    let deg = 2 * (maxp - 1).max(0) + 1;
    match deg {
        i64::MIN..=1 => 1,
        2 => 4,
        3 => 5,
        4 => 11,
        5 | 6 => 24,
        _ => 0,
    }
}

/// `TetraRuleTabulated(n)`.
pub fn tetra_rule_tabulated(n: usize) -> bool {
    [1, 4, 5, 11, 24].contains(&n)
}

/// `GaussPointsTetra(n, PReferenceElement=.TRUE.)`.
pub fn gauss_points_tetra_p(n: usize) -> Result<GaussRule, Error> {
    let mut p = fe(gauss_points_tetra(n))?;
    for i in 0..p.len() {
        let (uq, vq, wq, sq) = (p.u[i], p.v[i], p.w[i], p.s[i]);
        p.u[i] = -1.0 + 2.0 * uq + vq + wq;
        p.v[i] = 3.0f64.sqrt() * vq + 1.0 / 3.0f64.sqrt() * wq;
        p.w[i] = 8.0f64.sqrt() / 3.0f64.sqrt() * wq;
        p.s[i] = 8.0f64.sqrt() * 2.0 * sq;
    }
    Ok(p)
}

/// `GaussPointsPTetra(np)`: a collapsed `n x n x (n+1)` brick rule on the
/// p reference tetrahedron.
pub fn gauss_points_p_tetra(np: usize) -> Result<GaussRule, Error> {
    let n = cube_root_points(np);
    let mut p = gauss_points_p_brick(n, n, n + 1)?;
    for i in 0..p.len() {
        let (uh, vh, wh, sh) = (p.u[i], p.v[i], p.w[i], p.s[i]);
        p.u[i] = 1.0 / 4.0 * (uh - uh * vh - uh * wh + uh * vh * wh);
        p.v[i] = 3.0f64.sqrt() / 4.0 * (5.0 / 3.0 + vh - wh / 3.0 - vh * wh);
        p.w[i] = 6.0f64.sqrt() / 3.0 * (1.0 + wh);
        p.s[i] = -sh * 2.0f64.sqrt() / 16.0 * (1.0 - vh - wh + vh * wh) * (-1.0 + wh);
    }
    Ok(p)
}

/// `GaussPointsPPyramid(np)`: a collapsed `n x n x (n+1)` brick rule on the
/// p reference pyramid.
pub fn gauss_points_p_pyramid(np: usize) -> Result<GaussRule, Error> {
    let n = cube_root_points(np);
    let mut p = gauss_points_p_brick(n, n, n + 1)?;
    for i in 0..p.len() {
        let (uh, vh, wh, sh) = (p.u[i], p.v[i], p.w[i], p.s[i]);
        p.u[i] = 1.0 / 2.0 * uh * (1.0 - wh);
        p.v[i] = 1.0 / 2.0 * vh * (1.0 - wh);
        p.w[i] = 2.0f64.sqrt() / 2.0 * (1.0 + wh);
        p.s[i] = sh * 2.0f64.sqrt() / 8.0 * (-1.0 + wh).powi(2);
    }
    Ok(p)
}

/// The rule `GaussPoints` gives a p tetrahedron with `n` points
/// (`np_given`: `GaussPoints(Element, np)`).
pub fn p_tetra_rule(n: usize, np_given: bool) -> Result<GaussRule, Error> {
    let m = if np_given {
        if tetra_rule_tabulated(n) {
            n
        } else {
            0
        }
    } else {
        tetra_simplex_rule_points(n)
    };
    if m > 0 {
        gauss_points_tetra_p(m)
    } else {
        gauss_points_p_tetra(n)
    }
}

/// `TriangleRuleTabulated(n)`: whether `n` is one of the tabulated
/// triangle rules.
pub fn triangle_rule_tabulated(n: usize) -> bool {
    [1, 3, 4, 6, 7, 11, 12, 17, 20].contains(&n)
}

/// `GaussPointsTriangle(n, PReferenceElement=.TRUE.)`.
pub fn gauss_points_triangle_p(n: usize) -> Result<GaussRule, Error> {
    let mut p = fe(gauss_points_triangle(n))?;
    let sqrt3 = 3.0f64.sqrt();
    for i in 0..p.len() {
        let (uq, vq, sq) = (p.u[i], p.v[i], p.s[i]);
        p.u[i] = -1.0 + 2.0 * uq + vq;
        p.v[i] = sqrt3 * vq;
        p.s[i] = sqrt3 * 2.0 * sq;
    }
    Ok(p)
}

/// `GaussPointsPTriangle(n)`: the `n`-point quad rule collapsed onto the p
/// reference triangle.
pub fn gauss_points_p_triangle(n: usize) -> Result<GaussRule, Error> {
    let mut p = fe(gauss_points_quad(n))?;
    let sqrt3 = 3.0f64.sqrt();
    for i in 0..p.len() {
        let (uq, vq, sq) = (p.u[i], p.v[i], p.s[i]);
        p.u[i] = 1.0 / 2.0 * (uq - uq * vq);
        p.v[i] = sqrt3 / 2.0 * (1.0 + vq);
        p.s[i] = -sqrt3 / 4.0 * (-1.0 + vq) * sq;
    }
    Ok(p)
}

/// The rule `GaussPoints(Element)` gives a p triangle with `PDefs %
/// GaussPoints = n`.
pub fn p_triangle_rule(n: usize) -> Result<GaussRule, Error> {
    match triangle_simplex_rule_points(n) {
        0 => gauss_points_p_triangle(n),
        m => gauss_points_triangle_p(m),
    }
}

// Generated from Integration.F90 (UPQuad*/VPQuad*/SPQuad*) by gen_pquad.py; do not edit.
const U_PQUAD_8: [f64; 8] = [
    0.683130051063973e0,
    -0.683130051063973e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.881917103688197e0,
    0.881917103688197e0,
    -0.881917103688197e0,
    -0.881917103688197e0,
];
const V_PQUAD_8: [f64; 8] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.683130051063973e0,
    -0.683130051063973e0,
    0.881917103688197e0,
    -0.881917103688197e0,
    0.881917103688197e0,
    -0.881917103688197e0,
];
const S_PQUAD_8: [f64; 8] = [
    0.816326530612245e0,
    0.816326530612245e0,
    0.816326530612245e0,
    0.816326530612245e0,
    0.183673469387755e0,
    0.183673469387755e0,
    0.183673469387755e0,
    0.183673469387755e0,
];
const U_PQUAD_12: [f64; 12] = [
    0.925820099772551e0,
    -0.925820099772551e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.805979782918599e0,
    0.805979782918599e0,
    -0.805979782918599e0,
    -0.805979782918599e0,
    0.380554433208316e0,
    0.380554433208316e0,
    -0.380554433208316e0,
    -0.380554433208316e0,
];
const V_PQUAD_12: [f64; 12] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.925820099772551e0,
    -0.925820099772551e0,
    0.805979782918599e0,
    -0.805979782918599e0,
    0.805979782918599e0,
    -0.805979782918599e0,
    0.380554433208316e0,
    -0.380554433208316e0,
    0.380554433208316e0,
    -0.380554433208316e0,
];
const S_PQUAD_12: [f64; 12] = [
    0.241975308641975e0,
    0.241975308641975e0,
    0.241975308641975e0,
    0.241975308641975e0,
    0.237431774690630e0,
    0.237431774690630e0,
    0.237431774690630e0,
    0.237431774690630e0,
    0.520592916667394e0,
    0.520592916667394e0,
    0.520592916667394e0,
    0.520592916667394e0,
];
const U_PQUAD_20: [f64; 20] = [
    0.112122576386656e1,
    -0.112122576386656e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.451773049920657e0,
    -0.451773049920657e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.891849420851512e0,
    0.891849420851512e0,
    -0.891849420851512e0,
    -0.891849420851512e0,
    0.824396370749276e0,
    0.824396370749276e0,
    -0.824396370749276e0,
    -0.824396370749276e0,
    0.411623426336542e0,
    0.411623426336542e0,
    -0.411623426336542e0,
    -0.411623426336542e0,
];
const V_PQUAD_20: [f64; 20] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.112122576386656e1,
    -0.112122576386656e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.451773049920657e0,
    -0.451773049920657e0,
    0.891849420851512e0,
    -0.891849420851512e0,
    0.891849420851512e0,
    -0.891849420851512e0,
    0.411623426336542e0,
    -0.411623426336542e0,
    0.411623426336542e0,
    -0.411623426336542e0,
    0.824396370749276e0,
    -0.824396370749276e0,
    0.824396370749276e0,
    -0.824396370749276e0,
];
const S_PQUAD_20: [f64; 20] = [
    0.184758425074910e-1,
    0.184758425074910e-1,
    0.184758425074910e-1,
    0.184758425074910e-1,
    0.390052939160735e0,
    0.390052939160735e0,
    0.390052939160735e0,
    0.390052939160735e0,
    0.830951780264820e-1,
    0.830951780264820e-1,
    0.830951780264820e-1,
    0.830951780264820e-1,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
    0.254188020152646e0,
];
const U_PQUAD_25: [f64; 25] = [
    0.000000000000000e0,
    0.104440291540981e1,
    -0.104440291540981e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.769799068396649e0,
    -0.769799068396649e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.935787012440540e0,
    0.935787012440540e0,
    -0.935787012440540e0,
    -0.935787012440540e0,
    0.413491953449114e0,
    0.413491953449114e0,
    -0.413491953449114e0,
    -0.413491953449114e0,
    0.883025508525690e0,
    0.883025508525690e0,
    -0.883025508525690e0,
    -0.883025508525690e0,
    0.575653595840465e0,
    0.575653595840465e0,
    -0.575653595840465e0,
    -0.575653595840465e0,
];
const V_PQUAD_25: [f64; 25] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.104440291540981e1,
    -0.104440291540981e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.769799068396649e0,
    -0.769799068396649e0,
    0.935787012440540e0,
    -0.935787012440540e0,
    0.935787012440540e0,
    -0.935787012440540e0,
    0.413491953449114e0,
    -0.413491953449114e0,
    0.413491953449114e0,
    -0.413491953449114e0,
    0.575653595840465e0,
    -0.575653595840465e0,
    0.575653595840465e0,
    -0.575653595840465e0,
    0.883025508525690e0,
    -0.883025508525690e0,
    0.883025508525690e0,
    -0.883025508525690e0,
];
const S_PQUAD_25: [f64; 25] = [
    0.365379525585903e0,
    0.277561655642040e-1,
    0.277561655642040e-1,
    0.277561655642040e-1,
    0.277561655642040e-1,
    0.244272057751754e0,
    0.244272057751754e0,
    0.244272057751754e0,
    0.244272057751754e0,
    0.342651038512290e-1,
    0.342651038512290e-1,
    0.342651038512290e-1,
    0.342651038512290e-1,
    0.308993036133713e0,
    0.308993036133713e0,
    0.308993036133713e0,
    0.308993036133713e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
    0.146684377651312e0,
];
const U_PQUAD_36: [f64; 36] = [
    0.108605615857397e1,
    -0.108605615857397e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.658208197042585e0,
    -0.658208197042585e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.100130060299173e1,
    0.100130060299173e1,
    -0.100130060299173e1,
    -0.100130060299173e1,
    0.584636168775946e0,
    0.584636168775946e0,
    -0.584636168775946e0,
    -0.584636168775946e0,
    0.246795612720261e0,
    0.246795612720261e0,
    -0.246795612720261e0,
    -0.246795612720261e0,
    0.900258815287201e0,
    0.900258815287201e0,
    -0.900258815287201e0,
    -0.900258815287201e0,
    0.304720678579870e0,
    0.304720678579870e0,
    -0.304720678579870e0,
    -0.304720678579870e0,
    0.929866705560780e0,
    0.929866705560780e0,
    -0.929866705560780e0,
    -0.929866705560780e0,
    0.745052720131169e0,
    0.745052720131169e0,
    -0.745052720131169e0,
    -0.745052720131169e0,
];
const V_PQUAD_36: [f64; 36] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.108605615857397e1,
    -0.108605615857397e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.658208197042585e0,
    -0.658208197042585e0,
    0.100130060299173e1,
    -0.100130060299173e1,
    0.100130060299173e1,
    -0.100130060299173e1,
    0.584636168775946e0,
    -0.584636168775946e0,
    0.584636168775946e0,
    -0.584636168775946e0,
    0.246795612720261e0,
    -0.246795612720261e0,
    0.246795612720261e0,
    -0.246795612720261e0,
    0.304720678579870e0,
    -0.304720678579870e0,
    0.304720678579870e0,
    -0.304720678579870e0,
    0.900258815287201e0,
    -0.900258815287201e0,
    0.900258815287201e0,
    -0.900258815287201e0,
    0.745052720131169e0,
    -0.745052720131169e0,
    0.745052720131169e0,
    -0.745052720131169e0,
    0.929866705560780e0,
    -0.929866705560780e0,
    0.929866705560780e0,
    -0.929866705560780e0,
];
const S_PQUAD_36: [f64; 36] = [
    0.565616969376400e-2,
    0.565616969376400e-2,
    0.565616969376400e-2,
    0.565616969376400e-2,
    0.192443967470396e0,
    0.192443967470396e0,
    0.192443967470396e0,
    0.192443967470396e0,
    0.516683297977300e-2,
    0.516683297977300e-2,
    0.516683297977300e-2,
    0.516683297977300e-2,
    0.200302559622138e0,
    0.200302559622138e0,
    0.200302559622138e0,
    0.200302559622138e0,
    0.228125175912536e0,
    0.228125175912536e0,
    0.228125175912536e0,
    0.228125175912536e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.117496926974491e0,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
    0.666557701862050e-1,
];
const U_PQUAD_45: [f64; 45] = [
    0.000000000000000e0,
    0.102731435771937e1,
    -0.102731435771937e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.856766776147643e0,
    -0.856766776147643e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.327332998189723e0,
    -0.327332998189723e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.967223740028505e0,
    0.967223740028505e0,
    -0.967223740028505e0,
    -0.967223740028505e0,
    0.732168901749711e0,
    0.732168901749711e0,
    -0.732168901749711e0,
    -0.732168901749711e0,
    0.621974427996805e0,
    0.621974427996805e0,
    -0.621974427996805e0,
    -0.621974427996805e0,
    0.321696694921009e0,
    0.321696694921009e0,
    -0.321696694921009e0,
    -0.321696694921009e0,
    0.928618480068352e0,
    0.928618480068352e0,
    -0.928618480068352e0,
    -0.928618480068352e0,
    0.455124178121179e0,
    0.455124178121179e0,
    -0.455124178121179e0,
    -0.455124178121179e0,
    0.960457474887516e0,
    0.960457474887516e0,
    -0.960457474887516e0,
    -0.960457474887516e0,
    0.809863684081217e0,
    0.809863684081217e0,
    -0.809863684081217e0,
    -0.809863684081217e0,
];
const V_PQUAD_45: [f64; 45] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.102731435771937e1,
    -0.102731435771937e1,
    0.000000000000000e0,
    0.000000000000000e0,
    0.856766776147643e0,
    -0.856766776147643e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.327332998189723e0,
    -0.327332998189723e0,
    0.967223740028505e0,
    -0.967223740028505e0,
    0.967223740028505e0,
    -0.967223740028505e0,
    0.732168901749711e0,
    -0.732168901749711e0,
    0.732168901749711e0,
    -0.732168901749711e0,
    0.321696694921009e0,
    -0.321696694921009e0,
    0.321696694921009e0,
    -0.321696694921009e0,
    0.621974427996805e0,
    -0.621974427996805e0,
    0.621974427996805e0,
    -0.621974427996805e0,
    0.455124178121179e0,
    -0.455124178121179e0,
    0.455124178121179e0,
    -0.455124178121179e0,
    0.928618480068352e0,
    -0.928618480068352e0,
    0.928618480068352e0,
    -0.928618480068352e0,
    0.809863684081217e0,
    -0.809863684081217e0,
    0.809863684081217e0,
    -0.809863684081217e0,
    0.960457474887516e0,
    -0.960457474887516e0,
    0.960457474887516e0,
    -0.960457474887516e0,
];
const S_PQUAD_45: [f64; 45] = [
    -0.176897982720700e-2,
    0.128167266175120e-1,
    0.128167266175120e-1,
    0.128167266175120e-1,
    0.128167266175120e-1,
    0.119897873101347e0,
    0.119897873101347e0,
    0.119897873101347e0,
    0.119897873101347e0,
    0.210885452208801e0,
    0.210885452208801e0,
    0.210885452208801e0,
    0.210885452208801e0,
    0.639272012821500e-2,
    0.639272012821500e-2,
    0.639272012821500e-2,
    0.639272012821500e-2,
    0.104415680788580e0,
    0.104415680788580e0,
    0.104415680788580e0,
    0.104415680788580e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.168053047203816e0,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.761696944522940e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
    0.287941544000640e-1,
];
const U_PQUAD_60: [f64; 60] = [
    0.989353074512600e0,
    -0.989353074512600e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.376285207157973e0,
    -0.376285207157973e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.978848279262233e0,
    0.978848279262233e0,
    -0.978848279262233e0,
    -0.978848279262233e0,
    0.885794729164116e0,
    0.885794729164116e0,
    -0.885794729164116e0,
    -0.885794729164116e0,
    0.171756123838348e0,
    0.171756123838348e0,
    -0.171756123838348e0,
    -0.171756123838348e0,
    0.590499273806002e0,
    0.590499273806002e0,
    -0.590499273806002e0,
    -0.590499273806002e0,
    0.319505036634574e0,
    0.319505036634574e0,
    -0.319505036634574e0,
    -0.319505036634574e0,
    0.799079131916863e0,
    0.799079131916863e0,
    -0.799079131916863e0,
    -0.799079131916863e0,
    0.597972451929457e0,
    0.597972451929457e0,
    -0.597972451929457e0,
    -0.597972451929457e0,
    0.803743962958745e0,
    0.803743962958745e0,
    -0.803743962958745e0,
    -0.803743962958745e0,
    0.583444817765510e-1,
    0.583444817765510e-1,
    -0.583444817765510e-1,
    -0.583444817765510e-1,
    0.936506276127495e0,
    0.936506276127495e0,
    -0.936506276127495e0,
    -0.936506276127495e0,
    0.347386316166203e0,
    0.347386316166203e0,
    -0.347386316166203e0,
    -0.347386316166203e0,
    0.981321179805452e0,
    0.981321179805452e0,
    -0.981321179805452e0,
    -0.981321179805452e0,
    0.706000287798646e0,
    0.706000287798646e0,
    -0.706000287798646e0,
    -0.706000287798646e0,
];
const V_PQUAD_60: [f64; 60] = [
    0.000000000000000e0,
    0.000000000000000e0,
    0.989353074512600e0,
    -0.989353074512600e0,
    0.000000000000000e0,
    0.000000000000000e0,
    0.376285207157973e0,
    -0.376285207157973e0,
    0.978848279262233e0,
    -0.978848279262233e0,
    0.978848279262233e0,
    -0.978848279262233e0,
    0.885794729164116e0,
    -0.885794729164116e0,
    0.885794729164116e0,
    -0.885794729164116e0,
    0.171756123838348e0,
    -0.171756123838348e0,
    0.171756123838348e0,
    -0.171756123838348e0,
    0.319505036634574e0,
    -0.319505036634574e0,
    0.319505036634574e0,
    -0.319505036634574e0,
    0.590499273806002e0,
    -0.590499273806002e0,
    0.590499273806002e0,
    -0.590499273806002e0,
    0.597972451929457e0,
    -0.597972451929457e0,
    0.597972451929457e0,
    -0.597972451929457e0,
    0.799079131916863e0,
    -0.799079131916863e0,
    0.799079131916863e0,
    -0.799079131916863e0,
    0.583444817765510e-1,
    -0.583444817765510e-1,
    0.583444817765510e-1,
    -0.583444817765510e-1,
    0.803743962958745e0,
    -0.803743962958745e0,
    0.803743962958745e0,
    -0.803743962958745e0,
    0.347386316166203e0,
    -0.347386316166203e0,
    0.347386316166203e0,
    -0.347386316166203e0,
    0.936506276127495e0,
    -0.936506276127495e0,
    0.936506276127495e0,
    -0.936506276127495e0,
    0.706000287798646e0,
    -0.706000287798646e0,
    0.706000287798646e0,
    -0.706000287798646e0,
    0.981321179805452e0,
    -0.981321179805452e0,
    0.981321179805452e0,
    -0.981321179805452e0,
];
const S_PQUAD_60: [f64; 60] = [
    0.206149159199910e-1,
    0.206149159199910e-1,
    0.206149159199910e-1,
    0.206149159199910e-1,
    0.128025716179910e0,
    0.128025716179910e0,
    0.128025716179910e0,
    0.128025716179910e0,
    0.551173953403200e-2,
    0.551173953403200e-2,
    0.551173953403200e-2,
    0.551173953403200e-2,
    0.392077124571420e-1,
    0.392077124571420e-1,
    0.392077124571420e-1,
    0.392077124571420e-1,
    0.763969450798630e-1,
    0.763969450798630e-1,
    0.763969450798630e-1,
    0.763969450798630e-1,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.141513729949972e0,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.839032793637980e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.603941636496850e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.573877529692130e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
    0.219225594818640e-1,
];

#[cfg(test)]
mod tests {
    use super::*;

    /// Integral of `x^a y^b` over a rule.
    fn moment(r: &GaussRule, a: i32, b: i32) -> f64 {
        (0..r.len())
            .map(|t| r.s[t] * r.u[t].powi(a) * r.v[t].powi(b))
            .sum()
    }

    #[test]
    fn economic_quad_rules_integrate_their_degree() {
        // Exact integrals over [-1, 1]^2 of x^a y^b with a + b <= degree.
        for (np, deg) in [
            (8, 5),
            (12, 7),
            (20, 9),
            (25, 11),
            (36, 13),
            (45, 15),
            (60, 17),
        ] {
            let r = gauss_points_quad_p(np, true).unwrap();
            assert_eq!(r.len(), np);
            for a in 0..=deg {
                for b in 0..=deg - a {
                    let ex = |k: i32| {
                        if k % 2 == 1 {
                            0.0
                        } else {
                            2.0 / f64::from(k + 1)
                        }
                    };
                    let want = ex(a) * ex(b);
                    let got = moment(&r, a, b);
                    // Elmer's 36-point table carries weights summing to
                    // 4.0000004; the others are exact to the printed digits.
                    let tol = if np == 36 { 1e-6 } else { 1e-10 };
                    assert!(
                        (got - want).abs() < tol,
                        "np {np}: x^{a} y^{b}: {got} vs {want}"
                    );
                }
            }
        }
        assert_eq!(gauss_points_quad_p(16, false).unwrap().len(), 16);
        assert_eq!(gauss_points_quad_p(4, true).unwrap().len(), 4);
        assert!(gauss_points_quad_p(16, true).is_err());
    }

    /// Integral of `x^a y^b z^c` over a 3D rule.
    fn moment3(r: &GaussRule, a: i32, b: i32, c: i32) -> f64 {
        (0..r.len())
            .map(|t| r.s[t] * r.u[t].powi(a) * r.v[t].powi(b) * r.w[t].powi(c))
            .sum()
    }

    #[test]
    fn p_wedge_rules_integrate_on_the_p_reference_wedge() {
        // p reference wedge: the p triangle (area sqrt 3, centroid
        // (0, 1/sqrt 3)) times [-1, 1].
        let vol = 2.0 * 3.0f64.sqrt();
        let counts = [
            64, 44, 28, 27, 24, 21, 18, 16, 15, 14, 12, 11, 10, 8, 7, 6, 5,
        ];
        for np in counts {
            let r = p_wedge_rule(np, true).unwrap();
            assert!(
                (moment3(&r, 0, 0, 0) - vol).abs() < 1e-12,
                "np {np}: volume"
            );
            assert!(moment3(&r, 1, 0, 0).abs() < 1e-12, "np {np}: x");
            assert!((moment3(&r, 0, 1, 0) - 2.0).abs() < 1e-12, "np {np}: y");
            assert!(moment3(&r, 0, 0, 1).abs() < 1e-12, "np {np}: z");
            // Every rule of 6 points or more integrates z^2 exactly.
            if np >= 6 {
                assert!(
                    (moment3(&r, 0, 0, 2) - vol / 3.0).abs() < 1e-12,
                    "np {np}: z^2"
                );
            }
        }
        // Without np, PDefs % GaussPoints goes through WedgeRulePoints.
        // 27 points: degree 5, min(7 x 3, economic 15); 64: degree 7,
        // 12 x 4; 125: degree 9, 20 x 5.
        assert_eq!(wedge_rule_points(27), 15);
        assert_eq!(wedge_rule_points(64), 48);
        assert_eq!(wedge_rule_points(125), 100);
        assert_eq!(p_wedge_rule(27, false).unwrap().len(), 15);
        assert_eq!(p_wedge_rule(125, false).unwrap().len(), 100);
        assert_eq!(
            p_wedge_rule(13 * 13 * 13, false).unwrap().len(),
            13 * 13 * 13
        );
    }

    #[test]
    fn p_tetra_and_pyramid_rules_cover_their_reference_elements() {
        // p tetrahedron: the p triangle base (area sqrt 3) under the apex
        // (0, 1/sqrt 3, sqrt(8/3)); volume 2 sqrt 2 / 3, centroid height
        // sqrt(8/3) / 4.
        let vol = 2.0 * 2.0f64.sqrt() / 3.0;
        for (n, np_given) in [
            (1, true),
            (4, true),
            (5, true),
            (11, true),
            (24, true),
            (27, true),
            (8, false),
            (27, false),
            (64, false),
            (343, false),
        ] {
            let r = p_tetra_rule(n, np_given).unwrap();
            assert!(
                (moment3(&r, 0, 0, 0) - vol).abs() < 1e-13,
                "tetra {n} {np_given}: volume"
            );
            assert!(moment3(&r, 1, 0, 0).abs() < 1e-13, "tetra {n}: x");
            let zc = (8.0f64 / 3.0).sqrt() / 4.0;
            assert!(
                (moment3(&r, 0, 0, 1) - vol * zc).abs() < 1e-13,
                "tetra {n}: z"
            );
        }
        assert_eq!(tetra_simplex_rule_points(27), 24);
        assert_eq!(p_tetra_rule(343, false).unwrap().len(), 7 * 7 * 8);
        // p pyramid: base [-1, 1]^2, apex at height sqrt 2; volume
        // 4 sqrt 2 / 3, centroid height sqrt 2 / 4.
        let vol = 4.0 * 2.0f64.sqrt() / 3.0;
        for n in [1, 8, 27, 64] {
            let r = gauss_points_p_pyramid(n).unwrap();
            assert!(
                (moment3(&r, 0, 0, 0) - vol).abs() < 1e-13,
                "pyramid {n}: volume"
            );
            assert!(
                (moment3(&r, 0, 0, 1) - vol * 2.0f64.sqrt() / 4.0).abs() < 1e-13,
                "pyramid {n}: z"
            );
        }
    }

    #[test]
    fn p_triangle_rules_cover_the_p_reference_triangle() {
        // Area sqrt 3; centroid (0, 1 / sqrt 3).
        let area = 3.0f64.sqrt();
        for n in [1, 4, 9, 16, 25, 36, 49, 64, 81, 100] {
            let r = p_triangle_rule(n).unwrap();
            assert!((moment(&r, 0, 0) - area).abs() < 1e-13, "n {n}");
            assert!(moment(&r, 1, 0).abs() < 1e-13, "n {n}");
            assert!((moment(&r, 0, 1) - 1.0).abs() < 1e-13, "n {n}");
        }
        // NINT(SQRT(np)) points per direction: degree 2 (maxp - 1) + 1.
        assert_eq!(triangle_simplex_rule_points(4), 4);
        assert_eq!(triangle_simplex_rule_points(9), 7);
        assert_eq!(triangle_simplex_rule_points(25), 20);
        assert_eq!(triangle_simplex_rule_points(36), 0);
        assert_eq!(p_triangle_rule(36).unwrap().len(), 36);
    }
}
