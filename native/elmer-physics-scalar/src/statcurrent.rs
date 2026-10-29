//! StatCurrentSolve element kernels: `LocalMatrix` and `LocalMatrixBC` of
//! `fem/src/modules/StatCurrentSolve.F90` for the real potential.

use std::f64::consts::PI;

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::common::{check_len, check_opt, interp};
use crate::error::ScalarError;
use crate::potential::Coefficient;
use crate::statelec::PHYSICAL_EPS0;

/// Element-level switches for the StatCurrentSolve kernels.
#[derive(Clone, Debug, PartialEq)]
pub struct StatCurrentConfig {
    /// Coordinate-system dimension (`CoordinateSystemDimension()`), 1 to 3.
    pub dim: usize,
    /// Axisymmetric or cylindric coordinates: the measure becomes
    /// `2 pi r dr dz` with `r` the first coordinate.
    pub axisymmetric: bool,
    /// Build the displacement-current mass matrix (`Transient`).
    pub transient: bool,
    /// `Permittivity Of Vacuum`, multiplying `Relative Permittivity` in the
    /// mass matrix.
    pub eps0: f64,
    /// Requested Gauss points, or `None` for the element default.
    pub quad_points: Option<usize>,
}

impl StatCurrentConfig {
    /// Steady, Cartesian, default quadrature.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            axisymmetric: false,
            transient: false,
            eps0: PHYSICAL_EPS0,
            quad_points: None,
        }
    }
}

/// Bulk properties at the element nodes.
#[derive(Clone, Debug, PartialEq)]
pub struct StatCurrentAtNodes {
    /// `Electric Conductivity` (material).
    pub conductivity: Coefficient,
    /// `Relative Permittivity` (material), used only for the mass matrix.
    pub relative_permittivity: Option<Vec<f64>>,
    /// `Current Source` (body force).
    pub current_source: Option<Vec<f64>>,
}

/// Boundary conditions at the nodes of a boundary element.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StatCurrentBcAtNodes {
    /// `Current Density`: the given flux.
    pub current_density: Option<Vec<f64>>,
    /// `External Conductivity`: the Robin coefficient; switches Robin on.
    pub external_conductivity: Option<Vec<f64>>,
    /// `External Potential`.
    pub external_potential: Option<Vec<f64>>,
}

fn weight(axisymmetric: bool, w: f64, phi: &[f64], xyz: &[[f64; 3]]) -> f64 {
    if axisymmetric {
        w * 2.0 * PI * phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>()
    } else {
        w
    }
}

/// Bulk element matrices of `LocalMatrix` (StatCurrentSolve.F90:596) with
/// `nd = n` nodes: `stiff = int grad phi_p . sigma grad phi_q`, `force = int
/// source phi_p`, and when transient `mass = int eps0 eps_r grad phi_p . grad
/// phi_q` (the Fortran's displacement-current term).
pub fn stat_current_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &StatCurrentAtNodes,
    cfg: &StatCurrentConfig,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("electric conductivity", n, mat.conductivity.len())?;
    check_opt("relative permittivity", n, &mat.relative_permittivity)?;
    check_opt("current source", n, &mat.current_source)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!(
            "dim {} not in 1..=3",
            cfg.dim
        )));
    }
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.damp = None;
    out.mass = cfg.transient.then(|| vec![0.0; n * n]);

    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(ScalarError::Basis)?;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];
        let w = weight(cfg.axisymmetric, qp.weight * be.det_j, phi, xyz);

        for p in 0..n {
            for q in 0..n {
                out.stiff[p * n + q] += w * mat.conductivity.stiffness(phi, dphi, p, q, cfg.dim);
            }
        }
        if let (Some(mass), Some(eps)) = (out.mass.as_mut(), &mat.relative_permittivity) {
            let e = cfg.eps0 * interp(phi, eps);
            for p in 0..n {
                for q in 0..n {
                    let s: f64 = (0..3).map(|i| dphi[p][i] * dphi[q][i]).sum();
                    mass[p * n + q] += w * e * s;
                }
            }
        }
        if let Some(src) = &mat.current_source {
            let f = interp(phi, src);
            for p in 0..n {
                out.force[p] += w * f * phi[p];
            }
        }
    }
    Ok(())
}

/// Boundary element matrices of `LocalMatrixBC` (StatCurrentSolve.F90:721)
/// for the given current density and the Robin condition `C (u - u_ext)`:
/// `stiff += w C phi_q phi_p`, `force += w (J + C u_ext) phi_p`. `Farfield
/// Potential` is not here: it needs the outward normal, which the element
/// basis does not provide.
pub fn stat_current_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &StatCurrentBcAtNodes,
    cfg: &StatCurrentConfig,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_opt("current density", n, &bc.current_density)?;
    check_opt("external conductivity", n, &bc.external_conductivity)?;
    check_opt("external potential", n, &bc.external_potential)?;
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.mass = None;
    out.damp = None;

    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(ScalarError::Basis)?;
        let phi = &be.basis[..n];
        let w = weight(cfg.axisymmetric, qp.weight * be.det_j, phi, xyz);

        if let Some(j) = &bc.current_density {
            let f = interp(phi, j);
            for p in 0..n {
                out.force[p] += w * f * phi[p];
            }
        }
        if let Some(c) = &bc.external_conductivity {
            let c = interp(phi, c);
            let ext = bc
                .external_potential
                .as_deref()
                .map_or(0.0, |e| interp(phi, e));
            for p in 0..n {
                for q in 0..n {
                    out.stiff[p * n + q] += w * c * phi[q] * phi[p];
                }
                out.force[p] += w * c * ext * phi[p];
            }
        }
    }
    Ok(())
}
