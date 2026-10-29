//! Gauss integration rules of Integration.F90 for non-p elements: the
//! per-family rule builders and the `GaussPoints` dispatch (Integration.F90:2882).

use elmer_core::{ElementType, QuadPoint, Quadrature};

use crate::element_type::{element_type, powi};
use crate::error::FeError;
use crate::gauss1d::{points_1d, weights_1d, MAXN};
use crate::tables::*;

/// A rule as `GaussIntegrationPoints_t` holds it: point `t` is
/// `(u[t], v[t], w[t])` with weight `s[t]`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GaussRule {
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub w: Vec<f64>,
    pub s: Vec<f64>,
}

impl GaussRule {
    fn with_len(n: usize) -> Self {
        GaussRule {
            u: vec![0.0; n],
            v: vec![0.0; n],
            w: vec![0.0; n],
            s: vec![0.0; n],
        }
    }

    fn from_tables(u: &[f64], v: &[f64], w: Option<&[f64]>, s: &[f64], sdiv: f64) -> Self {
        GaussRule {
            u: u.to_vec(),
            v: v.to_vec(),
            w: w.map_or_else(|| vec![0.0; u.len()], <[f64]>::to_vec),
            s: s.iter().map(|x| x / sdiv).collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.s.len()
    }

    pub fn is_empty(&self) -> bool {
        self.s.is_empty()
    }

    pub fn points(&self) -> Vec<QuadPoint> {
        (0..self.len())
            .map(|t| QuadPoint {
                uvw: [self.u[t], self.v[t], self.w[t]],
                weight: self.s[t],
            })
            .collect()
    }
}

fn check_n(family: u16, n: usize, np: usize) -> Result<(), FeError> {
    if !(1..=MAXN).contains(&n) {
        Err(FeError::InvalidPoints { family, np })
    } else {
        Ok(())
    }
}

/// `n = REAL(np)**(1.0D0/3.0D0) + 0.5D0`: single-precision `REAL(np)`
/// widened, libm `pow`, truncation toward zero.
fn cube_root_count(np: usize) -> usize {
    let x = (np as f32) as f64;
    (x.powf(1.0 / 3.0) + 0.5) as usize
}

/// GaussPoints0D.
pub fn gauss_points_0d() -> GaussRule {
    GaussRule {
        u: vec![0.0],
        v: vec![0.0],
        w: vec![0.0],
        s: vec![1.0],
    }
}

/// GaussPoints1D.
pub fn gauss_points_1d(n: usize) -> Result<GaussRule, FeError> {
    check_n(2, n, n)?;
    Ok(GaussRule {
        u: points_1d(n).to_vec(),
        v: vec![0.0; n],
        w: vec![0.0; n],
        s: weights_1d(n).to_vec(),
    })
}

/// GaussPointsTriangle with `PReferenceElement` absent.
pub fn gauss_points_triangle(n: usize) -> Result<GaussRule, FeError> {
    Ok(match n {
        1 => GaussRule::from_tables(&UTRIANGLE_1P, &VTRIANGLE_1P, None, &STRIANGLE_1P, 2.0),
        3 => GaussRule::from_tables(&UTRIANGLE_3P, &VTRIANGLE_3P, None, &STRIANGLE_3P, 2.0),
        4 => GaussRule::from_tables(&UTRIANGLE_4P, &VTRIANGLE_4P, None, &STRIANGLE_4P, 2.0),
        6 => GaussRule::from_tables(&UTRIANGLE_6P, &VTRIANGLE_6P, None, &STRIANGLE_6P, 2.0),
        7 => GaussRule::from_tables(&UTRIANGLE_7P, &VTRIANGLE_7P, None, &STRIANGLE_7P, 2.0),
        11 => GaussRule::from_tables(&UTRIANGLE_11P, &VTRIANGLE_11P, None, &STRIANGLE_11P, 1.0),
        12 => GaussRule::from_tables(&UTRIANGLE_12P, &VTRIANGLE_12P, None, &STRIANGLE_12P, 1.0),
        17 => GaussRule::from_tables(&UTRIANGLE_17P, &VTRIANGLE_17P, None, &STRIANGLE_17P, 1.0),
        20 => GaussRule::from_tables(&UTRIANGLE_20P, &VTRIANGLE_20P, None, &STRIANGLE_20P, 1.0),
        _ => {
            // Collapsed quadrilateral rule, renormalised to the area 1/2.
            let mut p =
                gauss_points_quad(n).map_err(|_| FeError::InvalidPoints { family: 3, np: n })?;
            for i in 0..p.len() {
                p.v[i] = (p.v[i] + 1.0) / 2.0;
                p.u[i] = (p.u[i] + 1.0) / 2.0 * (1.0 - p.v[i]);
                p.s[i] *= 1.0 - p.v[i];
            }
            let sum = fortran_sum(&p.s);
            for s in &mut p.s {
                *s = 0.5 * *s / sum;
            }
            p
        }
    })
}

/// GaussPointsQuad with `PMethod` absent (the tensor-product rule).
pub fn gauss_points_quad(np: usize) -> Result<GaussRule, FeError> {
    // n = SQRT( REAL(np) ) + 0.5, all in default (single) precision.
    let n = ((np as f32).sqrt() + 0.5f32) as usize;
    check_n(4, n, np)?;
    let (pts, wts) = (points_1d(n), weights_1d(n));
    let mut p = GaussRule::with_len(n * n);
    let mut t = 0;
    for i in 0..n {
        for j in 0..n {
            p.u[t] = pts[j];
            p.v[t] = pts[i];
            p.s[t] = wts[i] * wts[j];
            t += 1;
        }
    }
    Ok(p)
}

/// GaussPointsTetra with `PReferenceElement` absent.
pub fn gauss_points_tetra(n: usize) -> Result<GaussRule, FeError> {
    let tab =
        |u: &[f64], v: &[f64], w: &[f64], s: &[f64]| GaussRule::from_tables(u, v, Some(w), s, 6.0);
    Ok(match n {
        1 => tab(&UTETRA_1P, &VTETRA_1P, &WTETRA_1P, &STETRA_1P),
        4 => tab(&UTETRA_4P, &VTETRA_4P, &WTETRA_4P, &STETRA_4P),
        5 => tab(&UTETRA_5P, &VTETRA_5P, &WTETRA_5P, &STETRA_5P),
        11 => tab(&UTETRA_11P, &VTETRA_11P, &WTETRA_11P, &STETRA_11P),
        24 => tab(&UTETRA_24P, &VTETRA_24P, &WTETRA_24P, &STETRA_24P),
        _ => {
            // Collapsed brick rule (Duffy-type map).
            let mut p =
                gauss_points_brick(n).map_err(|_| FeError::InvalidPoints { family: 5, np: n })?;
            for i in 0..p.len() {
                let mut scale = 0.5;
                p.u[i] = (p.u[i] + 1.0) * scale;
                p.v[i] = (p.v[i] + 1.0) * scale;
                p.w[i] = (p.w[i] + 1.0) * scale;
                p.s[i] *= powi(scale, 3);

                scale = 1.0 - p.w[i];
                p.u[i] *= scale;
                p.v[i] *= scale;
                p.s[i] *= powi(scale, 2);

                scale = 1.0 - p.v[i] / scale;
                p.u[i] *= scale;
                p.s[i] *= scale;
            }
            p
        }
    })
}

/// Shared tensor loop of the brick, wedge and pyramid rules: `u` runs
/// fastest, then `v`, then `w`.
fn tensor3(n: usize) -> GaussRule {
    let (pts, wts) = (points_1d(n), weights_1d(n));
    let mut p = GaussRule::with_len(n * n * n);
    let mut t = 0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                p.u[t] = pts[k];
                p.v[t] = pts[j];
                p.w[t] = pts[i];
                p.s[t] = wts[i] * wts[j] * wts[k];
                t += 1;
            }
        }
    }
    p
}

/// GaussPointsPyramid.
pub fn gauss_points_pyramid(np: usize) -> Result<GaussRule, FeError> {
    let n = cube_root_count(np);
    check_n(6, n, np)?;
    let mut p = tensor3(n);
    for t in 0..p.len() {
        p.w[t] = (p.w[t] + 1.0) / 2.0;
        p.u[t] *= 1.0 - p.w[t];
        p.v[t] *= 1.0 - p.w[t];
        p.s[t] = p.s[t] * powi(1.0 - p.w[t], 2) / 2.0;
    }
    Ok(p)
}

/// GaussPointsWedge (the collapsed tensor rule).
pub fn gauss_points_wedge(np: usize) -> Result<GaussRule, FeError> {
    let n = cube_root_count(np);
    check_n(7, n, np)?;
    let mut p = tensor3(n);
    for i in 0..p.len() {
        p.v[i] = (p.v[i] + 1.0) / 2.0;
        p.u[i] = (p.u[i] + 1.0) / 2.0 * (1.0 - p.v[i]);
        p.s[i] = p.s[i] * (1.0 - p.v[i]) / 4.0;
    }
    Ok(p)
}

/// GaussPointsWedge2 with `PReferenceElement` false: a tabulated `m`-point
/// triangle rule times the `n`-point segment rule.
pub fn gauss_points_wedge2(m: usize, n: usize) -> Result<GaussRule, FeError> {
    let (uq, vq, sq): (&[f64], &[f64], Vec<f64>) = match m {
        1 => (&UTRIANGLE_1P, &VTRIANGLE_1P, halve(&STRIANGLE_1P)),
        3 => (&UTRIANGLE_3P, &VTRIANGLE_3P, halve(&STRIANGLE_3P)),
        4 => (&UTRIANGLE_4P, &VTRIANGLE_4P, halve(&STRIANGLE_4P)),
        6 => (&UTRIANGLE_6P, &VTRIANGLE_6P, halve(&STRIANGLE_6P)),
        7 => (&UTRIANGLE_7P, &VTRIANGLE_7P, halve(&STRIANGLE_7P)),
        11 => (&UTRIANGLE_11P, &VTRIANGLE_11P, STRIANGLE_11P.to_vec()),
        12 => (&UTRIANGLE_12P, &VTRIANGLE_12P, STRIANGLE_12P.to_vec()),
        17 => (&UTRIANGLE_17P, &VTRIANGLE_17P, STRIANGLE_17P.to_vec()),
        20 => (&UTRIANGLE_20P, &VTRIANGLE_20P, STRIANGLE_20P.to_vec()),
        _ => {
            // Any other triangle count: the collapsed n^3 rule, as
            // GaussPointsWedge but keyed on the segment count.
            check_n(7, n, m * n)?;
            let mut p = tensor3(n);
            for i in 0..p.len() {
                p.v[i] = (p.v[i] + 1.0) / 2.0;
                p.u[i] = (p.u[i] + 1.0) / 2.0 * (1.0 - p.v[i]);
                p.s[i] = p.s[i] * (1.0 - p.v[i]) / 4.0;
            }
            return Ok(p);
        }
    };
    check_n(7, n, m * n)?;
    let (pts, wts) = (points_1d(n), weights_1d(n));
    let mut p = GaussRule::with_len(m * n);
    let mut t = 0;
    for i in 0..m {
        for j in 0..n {
            p.u[t] = uq[i];
            p.v[t] = vq[i];
            p.w[t] = pts[j];
            p.s[t] = sq[i] * wts[j];
            t += 1;
        }
    }
    Ok(p)
}

fn halve(s: &[f64]) -> Vec<f64> {
    s.iter().map(|x| x / 2.0).collect()
}

/// GaussPointsWedgeEconomic with `PReferenceElement` false.
pub fn gauss_points_wedge_economic(n: usize) -> Result<GaussRule, FeError> {
    let (u, v, w, s): (&[f64], &[f64], &[f64], &[f64]) = match n {
        4 => (&UWEDGE_4P, &VWEDGE_4P, &WWEDGE_4P, &SWEDGE_4P),
        5 => (&UWEDGE_5P, &VWEDGE_5P, &WWEDGE_5P, &SWEDGE_5P),
        6 => (&UWEDGE_6P, &VWEDGE_6P, &WWEDGE_6P, &SWEDGE_6P),
        7 => (&UWEDGE_7P, &VWEDGE_7P, &WWEDGE_7P, &SWEDGE_7P),
        10 => (&UWEDGE_10P, &VWEDGE_10P, &WWEDGE_10P, &SWEDGE_10P),
        11 => (&UWEDGE_11P, &VWEDGE_11P, &WWEDGE_11P, &SWEDGE_11P),
        14 => (&UWEDGE_14P, &VWEDGE_14P, &WWEDGE_14P, &SWEDGE_14P),
        15 => (&UWEDGE_15P, &VWEDGE_15P, &WWEDGE_15P, &SWEDGE_15P),
        16 => (&UWEDGE_16P, &VWEDGE_16P, &WWEDGE_16P, &SWEDGE_16P),
        24 => (&UWEDGE_24P, &VWEDGE_24P, &WWEDGE_24P, &SWEDGE_24P),
        _ => return Err(FeError::InvalidPoints { family: 7, np: n }),
    };
    Ok(GaussRule {
        u: u.iter().map(|x| (x + 1.0) / 2.0).collect(),
        v: v.iter().map(|x| (x + 1.0) / 2.0).collect(),
        w: w.to_vec(),
        s: s.iter().map(|x| x / 4.0).collect(),
    })
}

/// GaussPointsBrick.
pub fn gauss_points_brick(np: usize) -> Result<GaussRule, FeError> {
    let n = match np {
        8 => 2,
        27 => 3,
        64 => 4,
        _ => {
            let n = cube_root_count(np);
            check_n(8, n, np)?;
            n
        }
    };
    Ok(tensor3(n))
}

/// gfortran's inline `SUM`: a left-to-right accumulation from zero.
pub(crate) fn fortran_sum(x: &[f64]) -> f64 {
    x.iter().fold(0.0, |s, v| s + v)
}

/// The wedge `np` dispatch of `GaussPoints` when a count is named: the
/// segment x triangle rules, the economical rules, else the collapsed rule.
fn wedge_named(n: usize) -> Result<GaussRule, FeError> {
    let nseg = match n {
        1..=3 => n,
        6 | 8 => 2,
        12 | 18 | 21 => 3,
        28 | 44 | 48 => 4,
        85 | 100 => 5,
        _ => 0,
    };
    if let Some(ntri) = n.checked_div(nseg) {
        return gauss_points_wedge2(ntri, nseg);
    }
    if matches!(n, 4 | 5 | 7 | 10 | 11 | 14 | 15 | 16 | 24) {
        return gauss_points_wedge_economic(n);
    }
    gauss_points_wedge(n)
}

/// `GaussPoints(elm, np)` for a non-p element: `np = None` is the element
/// type's default rule.
pub fn gauss_points(etype: ElementType, np: Option<usize>) -> Result<GaussRule, FeError> {
    let def = element_type(etype.0)?;
    let n = np.unwrap_or(def.gauss_points);
    family_rule(def.family(), n, np.is_some())
}

/// `GaussPoints(elm, RelOrder=rel)` for a non-p element: -1, 0 and +1 pick
/// `GaussPoints0`, `GaussPoints` and `GaussPoints2` of elements.def.
pub fn gauss_points_rel(etype: ElementType, rel: i32) -> Result<GaussRule, FeError> {
    let def = element_type(etype.0)?;
    let n = match rel {
        0 => def.gauss_points,
        1 => def.gauss_points2,
        -1 => def.gauss_points0,
        r => return Err(FeError::InvalidRelOrder(r)),
    };
    family_rule(def.family(), n, false)
}

fn family_rule(family: u16, n: usize, named: bool) -> Result<GaussRule, FeError> {
    match family {
        1 => Ok(gauss_points_0d()),
        2 => gauss_points_1d(n),
        3 => gauss_points_triangle(n),
        4 => gauss_points_quad(n),
        5 => gauss_points_tetra(n),
        6 => gauss_points_pyramid(n),
        7 if named => wedge_named(n),
        7 => gauss_points_wedge(n),
        8 => gauss_points_brick(n),
        f => Err(FeError::InvalidPoints { family: f, np: n }),
    }
}

/// [`Quadrature`] over Elmer's non-p Gauss rules.
#[derive(Copy, Clone, Debug, Default)]
pub struct ElmerQuadrature;

impl Quadrature for ElmerQuadrature {
    /// Panics where Elmer calls `Fatal`: an unknown element type or a point
    /// count with no rule. Use [`gauss_points`] to handle those as errors.
    fn gauss_points(&self, etype: ElementType, np: Option<usize>) -> Vec<QuadPoint> {
        match gauss_points(etype, np) {
            Ok(r) => r.points(),
            Err(e) => panic!("GaussPoints: {e}"),
        }
    }
}
