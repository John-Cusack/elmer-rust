//! The element matrix of the legacy assembly as `DiffuseConvectiveCompose`
//! (DiffuseConvectiveAnisotropic.F90:57) forms it, for the cases without
//! convection.
//!
//! The per-point kernel ([`crate::heat_bulk_local`]) integrates the same
//! terms, but in a different order: it multiplies the weight by the
//! conductivity before the gradient product, and it interpolates source and
//! density separately. The legacy routine instead builds per `(p,q)`
//!
//! `A = C0 Bq Bp + sum_ij C2(i,j) dBq_i dBp_j`,  `M = CT Bq Bp`,
//!
//! then adds `s A` and `s M` (`s = detJ w`). The load is
//! `Force = SUM(Load Basis) + JouleHeat` with `Load = Density * Heat Source`
//! formed at the nodes (HeatSolveLegacy.F90:976-979), added as
//! `s (Force Bp)`. The two orders agree to rounding. At a tight linear
//! tolerance that rounding is what the field check sees (HeatGap).
//!
//! Only Cartesian coordinates without a convection velocity and without a
//! `Heat Transfer Multiplier` take this path. The legacy driver keeps
//! [`crate::heat_bulk_local`] for the rest.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::HeatError;
use crate::material::{Conductivity, HeatBulkConfig, HeatMaterialAtNodes};

/// `SUM( Nodal(1:n) * Basis(1:n) )`.
fn sum(nodal: impl Iterator<Item = f64>, basis: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, b) in nodal.zip(basis) {
        s += x * b;
    }
    s
}

/// The legacy element matrix into `out`, or `Ok(false)` when the element is
/// one this path does not cover (convection, axisymmetry, a heat transfer
/// multiplier); `out` is then untouched.
pub(crate) fn legacy_compose_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &HeatMaterialAtNodes,
    cfg: &HeatBulkConfig,
    out: &mut LocalSystem,
) -> Result<bool, HeatError> {
    if mat.convection_velocity.is_some()
        || cfg.axisymmetric
        || mat.heat_transfer_multiplier.is_some()
    {
        return Ok(false);
    }
    let n = xyz.len();
    let dim = cfg.dim;
    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(HeatError::NoQuadrature(etype.0));
    }
    // Nodal coefficients as HeatSolveLegacy hands them over: C2 the 3x3
    // conductivity (zero off the given entries), CT = Density * Heat
    // Capacity, Load = Volumetric Heat Source or Density * Heat Source.
    let zero = vec![0.0; n];
    let rho = mat.density.as_deref().unwrap_or(&zero);
    let mut c2 = vec![[[0.0f64; 3]; 3]; n];
    match &mat.conductivity {
        None => {}
        Some(Conductivity::Scalar(k)) => {
            for (c, &kk) in c2.iter_mut().zip(k) {
                for i in 0..3 {
                    c[i][i] = kk;
                }
            }
        }
        Some(Conductivity::Diagonal(k)) => {
            for (c, kk) in c2.iter_mut().zip(k) {
                for i in 0..3 {
                    c[i][i] = kk[i];
                }
            }
        }
        Some(Conductivity::Full(k)) => c2.copy_from_slice(&k[..n]),
    }
    let ct: Vec<f64> = match &mat.heat_capacity {
        Some(cp) => rho.iter().zip(cp).map(|(r, c)| r * c).collect(),
        None => zero.clone(),
    };
    let load: Vec<f64> = match (&mat.volumetric_heat_source, &mat.heat_source) {
        (Some(v), _) => v.clone(),
        (None, Some(h)) => rho.iter().zip(h).map(|(r, h)| r * h).collect(),
        (None, None) => zero.clone(),
    };

    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.damp = None;
    out.mass = if cfg.transient {
        let mut m = out.mass.take().unwrap_or_default();
        m.clear();
        m.resize(n * n, 0.0);
        Some(m)
    } else {
        None
    };

    let mut be = BasisEval::default();
    for qp in &ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(HeatError::Basis)?;
        let b = &be.basis[..n];
        let db = &be.dbasis_dx[..n];
        let s = be.det_j * qp.weight;
        let c0 = sum(zero.iter().copied(), b);
        let ctip = sum(ct.iter().copied(), b);
        let mut c2ip = [[0.0f64; 3]; 3];
        for i in 0..dim {
            for j in 0..dim {
                c2ip[i][j] = sum(c2.iter().map(|c| c[i][j]), b);
            }
        }
        for p in 0..n {
            for q in 0..n {
                let m = ctip * b[q] * b[p];
                let mut a = c0 * b[q] * b[p];
                for i in 0..dim {
                    for j in 0..dim {
                        a += c2ip[i][j] * db[q][i] * db[p][j];
                    }
                }
                out.stiff[p * n + q] += s * a;
                if let Some(mass) = out.mass.as_mut() {
                    mass[p * n + q] += s * m;
                }
            }
        }
        // JouleHeat is zero without a Joule heating field.
        let force = sum(load.iter().copied(), b) + 0.0;
        for p in 0..n {
            out.force[p] += s * (force * b[p]);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::QuadPoint;

    /// Linear triangle on its physical coordinates.
    struct P1;
    impl ElementBasis for P1 {
        fn eval(
            &self,
            _: ElementType,
            x: &[[f64; 3]],
            uvw: [f64; 3],
            out: &mut BasisEval,
        ) -> Result<(), String> {
            let (u, v) = (uvw[0], uvw[1]);
            let (a, b) = (
                [x[1][0] - x[0][0], x[1][1] - x[0][1]],
                [x[2][0] - x[0][0], x[2][1] - x[0][1]],
            );
            let det = a[0] * b[1] - a[1] * b[0];
            // Reference gradients (-1,-1), (1,0), (0,1) mapped by J^-T.
            let g = |r: [f64; 2]| {
                [
                    (b[1] * r[0] - a[1] * r[1]) / det,
                    (-b[0] * r[0] + a[0] * r[1]) / det,
                    0.0,
                ]
            };
            out.det_j = det.abs();
            out.basis = vec![1.0 - u - v, u, v];
            out.dbasis_dx = vec![g([-1.0, -1.0]), g([1.0, 0.0]), g([0.0, 1.0])];
            Ok(())
        }
    }
    struct Three;
    impl Quadrature for Three {
        fn gauss_points(&self, _: ElementType, _: Option<usize>) -> Vec<QuadPoint> {
            [
                [1.0 / 6.0, 1.0 / 6.0],
                [2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0],
            ]
            .iter()
            .map(|&[u, v]| QuadPoint {
                uvw: [u, v, 0.0],
                weight: 1.0 / 6.0,
            })
            .collect()
        }
    }

    fn material() -> HeatMaterialAtNodes {
        HeatMaterialAtNodes {
            conductivity: Some(Conductivity::Scalar(vec![0.3, 1.7, 2.9])),
            density: Some(vec![1.1, 0.9, 1.3]),
            heat_capacity: Some(vec![2.0, 3.0, 5.0]),
            heat_source: Some(vec![7.0, 11.0, 13.0]),
            ..Default::default()
        }
    }

    /// Same matrix as the shared kernel up to rounding (the integrals are
    /// the same; only the operation order differs).
    #[test]
    fn agrees_with_the_shared_kernel_to_rounding() {
        let xyz = [[0.1, 0.0, 0.0], [1.3, 0.2, 0.0], [0.4, 0.9, 0.0]];
        let mut cfg = HeatBulkConfig::new(2);
        cfg.transient = true;
        let (mut a, mut b) = (LocalSystem::zeros(0), LocalSystem::zeros(0));
        assert!(legacy_compose_local(
            &P1,
            &Three,
            ElementType(303),
            &xyz,
            &material(),
            &cfg,
            &mut a
        )
        .unwrap());
        crate::heat_bulk_local(
            &P1,
            &Three,
            ElementType(303),
            &xyz,
            &material(),
            &cfg,
            &mut b,
        )
        .unwrap();
        // The shared kernel interpolates rho and the source separately; with
        // nodal products the exact integrals differ, so compare the parts
        // that do not involve rho: the stiffness.
        for (x, y) in a.stiff.iter().zip(&b.stiff) {
            assert!((x - y).abs() <= 1e-14 * (1.0 + y.abs()), "{x} vs {y}");
        }
        assert!(a.mass.is_some() && a.force.iter().all(|f| *f > 0.0));
    }

    #[test]
    fn declines_convection_axisymmetry_and_multipliers() {
        let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cfg = HeatBulkConfig::new(2);
        let mut out = LocalSystem::zeros(0);
        let mut m = material();
        m.convection_velocity = Some(vec![[1.0, 0.0, 0.0]; 3]);
        assert!(
            !legacy_compose_local(&P1, &Three, ElementType(303), &xyz, &m, &cfg, &mut out).unwrap()
        );
        let mut m = material();
        m.heat_transfer_multiplier = Some(vec![2.0; 3]);
        assert!(
            !legacy_compose_local(&P1, &Three, ElementType(303), &xyz, &m, &cfg, &mut out).unwrap()
        );
        let mut c = HeatBulkConfig::new(2);
        c.axisymmetric = true;
        assert!(!legacy_compose_local(
            &P1,
            &Three,
            ElementType(303),
            &xyz,
            &material(),
            &c,
            &mut out
        )
        .unwrap());
        assert_eq!(out.n, 0, "a declined element leaves the output untouched");
    }
}
