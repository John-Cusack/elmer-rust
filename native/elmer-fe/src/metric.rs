//! Element metric: `ElementMetric`, `CheckMetric` and `ElementMetricQP`
//! (ElemInfo.F90:9793-10186) and `InvertMatrix3x3` (LinearAlgebra.F90:830).
//!
//! The pristine build defines `HAVE_QP`: the first `ElementInfo` call on an
//! element compares the metric determinant in double and quad precision and
//! fixes the element's `Status`; status 2 elements take the quad-precision
//! metric from then on, and a double-precision metric that fails the
//! degeneracy test retries in quad precision. Quad precision is emulated with
//! double-double ([`crate::dd`]).

use crate::dd::Dd;
use crate::error::FeError;

/// `Element % Status`: 0 = not yet checked, 1 = double precision is
/// accurate, 2 = use the quad-precision metric.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum MetricStatus {
    #[default]
    Unchecked,
    Double,
    Quad,
}

/// Result of `ElementMetric`: the contravariant metric, `sqrt(det G)` and the
/// local-to-global gradient map `LtoGMap(i,j)` (`ltog[i][j]`).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Metric {
    pub metric: [[f64; 3]; 3],
    pub det_g: f64,
    pub ltog: [[f64; 3]; 3],
}

/// `dx(c,i) = SUM( x_c(1:n) * dLBasisdx(1:n,i) )` for `i = 1..dim`.
fn covariant(xyz: &[[f64; 3]], dl: &[[f64; 3]], dim: usize) -> [[f64; 3]; 3] {
    let mut dx = [[0.0; 3]; 3];
    for i in 0..dim {
        for c in 0..3 {
            let mut s = 0.0;
            for (x, d) in xyz.iter().zip(dl) {
                s += x[c] * d[i];
            }
            dx[c][i] = s;
        }
    }
    dx
}

/// `CheckMetric`: does the double-precision determinant agree with the
/// quad-precision one to a relative 1e-6?
pub fn check_metric(xyz: &[[f64; 3]], dl: &[[f64; 3]], dim: usize, cdim: usize) -> bool {
    let dx = covariant(xyz, dl, dim);
    let mut g = [[0.0; 3]; 3];
    let mut gq = [[Dd::ZERO; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            let (mut s, mut sq) = (0.0, Dd::ZERO);
            for row in dx.iter().take(cdim) {
                s += row[i] * row[j];
                sq = sq.add(Dd::from(row[i]).mul(Dd::from(row[j])));
            }
            g[i][j] = s;
            gq[i][j] = sq;
        }
    }
    let (dp, qp) = match dim {
        1 => (g[0][0], gq[0][0]),
        2 => (
            g[0][0] * g[1][1] - g[0][1] * g[1][0],
            gq[0][0].mul(gq[1][1]).sub(gq[0][1].mul(gq[1][0])),
        ),
        _ => (det3(&g), det3_dd(&gq)),
    };
    Dd::from(dp)
        .sub(qp)
        .abs()
        .le(Dd::from(1.0e-6).mul(qp.abs()))
}

fn det3(g: &[[f64; 3]; 3]) -> f64 {
    g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
        + g[0][1] * (g[1][2] * g[2][0] - g[1][0] * g[2][2])
        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0])
}

fn det3_dd(g: &[[Dd; 3]; 3]) -> Dd {
    let a = g[0][0].mul(g[1][1].mul(g[2][2]).sub(g[1][2].mul(g[2][1])));
    let b = g[0][1].mul(g[1][2].mul(g[2][0]).sub(g[1][0].mul(g[2][2])));
    let c = g[0][2].mul(g[1][0].mul(g[2][1]).sub(g[1][1].mul(g[2][0])));
    a.add(b).add(c)
}

/// `InvertMatrix3x3( G, GI, detG )`.
fn invert3(g: &[[f64; 3]; 3], det_g: f64) -> [[f64; 3]; 3] {
    let s = 1.0 / det_g;
    let mut gi = [[0.0; 3]; 3];
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

/// `InvertMatrix3x3QP`, each entry rounded to f64 as `Metric = GI` does.
fn invert3_dd(g: &[[Dd; 3]; 3], det_g: Dd) -> [[f64; 3]; 3] {
    let s = Dd::from(1.0).div(det_g);
    let m = |a: Dd, b: Dd, c: Dd, d: Dd| a.mul(b).sub(c.mul(d));
    let mut gi = [[0.0; 3]; 3];
    gi[0][0] = s.mul(m(g[1][1], g[2][2], g[2][1], g[1][2])).to_f64();
    gi[1][0] = s.neg().mul(m(g[1][0], g[2][2], g[2][0], g[1][2])).to_f64();
    gi[2][0] = s.mul(m(g[1][0], g[2][1], g[2][0], g[1][1])).to_f64();
    gi[0][1] = s.neg().mul(m(g[0][1], g[2][2], g[2][1], g[0][2])).to_f64();
    gi[1][1] = s.mul(m(g[0][0], g[2][2], g[2][0], g[0][2])).to_f64();
    gi[2][1] = s.neg().mul(m(g[0][0], g[2][1], g[2][0], g[0][1])).to_f64();
    gi[0][2] = s.mul(m(g[0][1], g[1][2], g[1][1], g[0][2])).to_f64();
    gi[1][2] = s.neg().mul(m(g[0][0], g[1][2], g[1][0], g[0][2])).to_f64();
    gi[2][2] = s.mul(m(g[0][0], g[1][1], g[1][0], g[0][1])).to_f64();
    gi
}

/// The double-precision body of `ElementMetric`; `Err(detg)` on the
/// degeneracy test.
fn metric_double(dx: &[[f64; 3]; 3], dim: usize, cdim: usize) -> Result<Metric, f64> {
    let eps = crate::element_type::powi(f64::EPSILON, dim as u32);
    let mut g = [[0.0; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            let mut s = 0.0;
            for row in dx.iter().take(cdim) {
                s += row[i] * row[j];
            }
            g[i][j] = s;
        }
    }
    let mut m = Metric::default();
    match dim {
        1 => {
            let det = g[0][0];
            if det <= eps {
                return Err(det);
            }
            m.metric[0][0] = 1.0 / det;
            m.det_g = det.sqrt();
        }
        2 => {
            let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
            if det <= eps {
                return Err(det);
            }
            m.metric[0][0] = g[1][1] / det;
            m.metric[0][1] = -g[0][1] / det;
            m.metric[1][0] = -g[1][0] / det;
            m.metric[1][1] = g[0][0] / det;
            m.det_g = det.sqrt();
        }
        _ => {
            let det = det3(&g);
            if det <= eps {
                return Err(det);
            }
            m.metric = invert3(&g, det);
            m.det_g = det.sqrt();
        }
    }
    ltog(dx, &mut m, dim, cdim);
    Ok(m)
}

// Index loops keep the Fortran summation order explicit.
#[allow(clippy::needless_range_loop)]
fn ltog(dx: &[[f64; 3]; 3], m: &mut Metric, dim: usize, cdim: usize) {
    for (row, out) in dx.iter().zip(m.ltog.iter_mut()).take(cdim) {
        for j in 0..dim {
            let mut s = 0.0;
            for k in 0..dim {
                s += row[k] * m.metric[k][j];
            }
            out[j] = s;
        }
    }
}

/// `ElementMetricQP`; `None` when `DetG <= TINY`.
#[allow(clippy::needless_range_loop)]
fn metric_quad(dx: &[[f64; 3]; 3], dim: usize, cdim: usize) -> Option<Metric> {
    let dq = dx.map(|r| r.map(Dd::from));
    let mut g = [[Dd::ZERO; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            let mut s = Dd::ZERO;
            for row in dq.iter().take(cdim) {
                s = s.add(row[i].mul(row[j]));
            }
            g[i][j] = s;
        }
    }
    let tiny = Dd::from(f64::MIN_POSITIVE);
    let mut m = Metric::default();
    let det = match dim {
        1 => {
            let det = g[0][0];
            if det.le(tiny) {
                return None;
            }
            m.metric[0][0] = Dd::from(1.0).div(det).to_f64();
            det
        }
        2 => {
            let det = g[0][0].mul(g[1][1]).sub(g[0][1].mul(g[1][0]));
            if det.le(tiny) {
                return None;
            }
            m.metric[0][0] = g[1][1].div(det).to_f64();
            m.metric[0][1] = g[0][1].neg().div(det).to_f64();
            m.metric[1][0] = g[1][0].neg().div(det).to_f64();
            m.metric[1][1] = g[0][0].div(det).to_f64();
            det
        }
        _ => {
            let det = det3_dd(&g);
            if det.le(tiny) {
                return None;
            }
            m.metric = invert3_dd(&g, det);
            det
        }
    };
    m.det_g = det.sqrt().to_f64();
    for (row, out) in dq.iter().zip(m.ltog.iter_mut()).take(cdim) {
        for j in 0..dim {
            let mut s = Dd::ZERO;
            for k in 0..dim {
                s = s.add(row[k].mul(Dd::from(m.metric[k][j])));
            }
            out[j] = s.to_f64();
        }
    }
    Some(m)
}

/// `ElementMetric` with the `HAVE_QP` fallbacks; `status` must already be
/// checked (see [`check_metric`]).
pub fn element_metric(
    xyz: &[[f64; 3]],
    dl: &[[f64; 3]],
    dim: usize,
    cdim: usize,
    status: MetricStatus,
) -> Result<Metric, FeError> {
    let dx = covariant(xyz, dl, dim);
    if status == MetricStatus::Quad {
        return metric_quad(&dx, dim, cdim).ok_or(FeError::Degenerate { dim, detg: 0.0 });
    }
    match metric_double(&dx, dim, cdim) {
        Ok(m) => Ok(m),
        Err(detg) => metric_quad(&dx, dim, cdim).ok_or(FeError::Degenerate { dim, detg }),
    }
}
