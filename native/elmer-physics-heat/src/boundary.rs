use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::HeatError;
use crate::material::{check_opt, interp, HeatBcAtNodes, HeatBoundaryConfig};

/// Boundary element matrices of `LocalMatrixBC` (HeatSolve.F90:2399) for the
/// `Heat Flux` and Robin (`Heat Transfer Coefficient`, `External
/// Temperature`) conditions.
///
/// `out` is reset to `n x n` with no mass: `force` gets `w * q * phi` for the
/// flux and `w * c * Text * phi` for the Robin source, `stiff` gets
/// `w * c * phi_p * phi_q`. The weight is the quadrature weight times the
/// boundary Jacobian, times `r` when axisymmetric.
///
/// Radiation, `Farfield Temperature`, `Infinity BC`, `Phase Change` and
/// radiator conditions are not implemented (see [`crate::check_supported`]).
pub fn heat_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &HeatBcAtNodes,
    cfg: &HeatBoundaryConfig,
    out: &mut LocalSystem,
) -> Result<(), HeatError> {
    let n = xyz.len();
    check_opt("heat flux", n, &bc.heat_flux)?;
    check_opt(
        "heat transfer coefficient",
        n,
        &bc.heat_transfer_coefficient,
    )?;
    check_opt("external temperature", n, &bc.external_temperature)?;

    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.mass = None;
    out.damp = None;

    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(HeatError::NoQuadrature(etype.0));
    }

    let mut be = BasisEval::default();
    for qp in &ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(HeatError::Basis)?;
        let phi = &be.basis[..n];

        let mut w = qp.weight * be.det_j;
        if cfg.axisymmetric {
            w *= phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>();
        }

        if let Some(q) = &bc.heat_flux {
            let f = interp(phi, q);
            for p in 0..n {
                out.force[p] += w * f * phi[p];
            }
        }

        if let Some(h) = &bc.heat_transfer_coefficient {
            let c = interp(phi, h);
            let text = bc
                .external_temperature
                .as_deref()
                .map_or(0.0, |t| interp(phi, t));
            for p in 0..n {
                for q in 0..n {
                    out.stiff[p * n + q] += w * c * phi[q] * phi[p];
                }
                out.force[p] += w * c * text * phi[p];
            }
        }
    }
    Ok(())
}
