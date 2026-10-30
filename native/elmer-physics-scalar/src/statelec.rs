//! StatElecSolve element kernels: `LocalMatrix`, `LocalMatrixBC` and
//! `LocalPostAssembly` of `fem/src/modules/StatElecSolve.F90` for the real
//! (single potential) path.

use std::f64::consts::PI;

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint, Quadrature};

use crate::common::{check_len, check_opt, interp};
use crate::error::ScalarError;

/// `PhysicalEps0`: the permittivity of vacuum StatElecSolve uses when
/// `Constants` has no `Permittivity Of Vacuum` (StatElecSolve.F90:49).
pub const PHYSICAL_EPS0: f64 = 8.854187817e-12;

/// Element-level switches for the StatElecSolve kernels.
#[derive(Clone, Debug, PartialEq)]
pub struct StatElecConfig {
    /// Coordinate-system dimension (`CoordinateSystemDimension()`), 1 to 3.
    pub dim: usize,
    /// Axisymmetric or cylindric coordinates: the measure becomes
    /// `2 pi r dr dz` with `r` the first coordinate.
    pub axisymmetric: bool,
    /// `Permittivity Of Vacuum`, multiplying the material's relative value.
    pub eps0: f64,
    /// Requested Gauss points, or `None` for the element default.
    pub quad_points: Option<usize>,
}

impl StatElecConfig {
    /// Cartesian, default quadrature, the physical vacuum permittivity.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            axisymmetric: false,
            eps0: PHYSICAL_EPS0,
            quad_points: None,
        }
    }
}

/// Bulk properties at the element nodes. A missing keyword is a zero array,
/// as `ListGetElementReal` returns it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StatElecAtNodes {
    /// `Relative Permittivity` (material).
    pub relative_permittivity: Vec<f64>,
    /// `Charge Density` (body force).
    pub charge_density: Option<Vec<f64>>,
}

impl StatElecAtNodes {
    /// The same values at `nodes`, given at `mesh_nodes` (the same set in
    /// another order).
    pub fn reordered(&self, mesh_nodes: &[u32], nodes: &[u32]) -> Self {
        let map: Vec<usize> = nodes
            .iter()
            .map(|n| {
                mesh_nodes
                    .iter()
                    .position(|m| m == n)
                    .expect("same node set")
            })
            .collect();
        let f = |v: &[f64]| map.iter().map(|&i| v[i]).collect();
        Self {
            relative_permittivity: f(&self.relative_permittivity),
            charge_density: self.charge_density.as_deref().map(f),
        }
    }
}

/// Boundary conditions at the nodes of a boundary element.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StatElecBcAtNodes {
    /// `Electric Flux`, else `Surface Charge Density`: the given `D.n`.
    pub electric_flux: Option<Vec<f64>>,
    /// `Layer Relative Permittivity`; switches the layer condition on.
    pub layer_relative_permittivity: Option<Vec<f64>>,
    /// `Layer Thickness`, required with the layer permittivity.
    pub layer_thickness: Option<Vec<f64>>,
    /// `Electrode Potential` of the layer.
    pub electrode_potential: Option<Vec<f64>>,
    /// `Layer Charge Density`.
    pub layer_charge_density: Option<Vec<f64>>,
}

fn reset(out: &mut LocalSystem, n: usize) {
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.mass = None;
    out.damp = None;
}

fn weight(cfg: &StatElecConfig, w: f64, phi: &[f64], xyz: &[[f64; 3]]) -> f64 {
    if cfg.axisymmetric {
        w * 2.0 * PI * phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>()
    } else {
        w
    }
}

impl StatElecBcAtNodes {
    /// The same values at `nodes`, given at `mesh_nodes` (the same set in
    /// another order).
    pub fn reordered(&self, mesh_nodes: &[u32], nodes: &[u32]) -> Self {
        let map: Vec<usize> = nodes
            .iter()
            .map(|n| {
                mesh_nodes
                    .iter()
                    .position(|m| m == n)
                    .expect("same node set")
            })
            .collect();
        let f = |v: &Option<Vec<f64>>| v.as_ref().map(|v| map.iter().map(|&i| v[i]).collect());
        Self {
            electric_flux: f(&self.electric_flux),
            layer_relative_permittivity: f(&self.layer_relative_permittivity),
            layer_thickness: f(&self.layer_thickness),
            electrode_potential: f(&self.electrode_potential),
            layer_charge_density: f(&self.layer_charge_density),
        }
    }
}

/// Bulk element matrices of `LocalMatrix` (StatElecSolve.F90:748) for an
/// element with `nd = n` nodes: `stiff[p*n + q] = int eps0 eps_r grad(phi_q) .
/// grad(phi_p)` and `force[p] = int rho phi_p`. Gradients use all three
/// columns of `dbasis_dx`, as the Fortran's `MATMUL`.
pub fn stat_elec_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &StatElecAtNodes,
    cfg: &StatElecConfig,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("relative permittivity", n, mat.relative_permittivity.len())?;
    check_opt("charge density", n, &mat.charge_density)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!(
            "dim {} not in 1..=3",
            cfg.dim
        )));
    }
    let ips = quad.gauss_points(etype, cfg.quad_points);
    bulk_with(basis, &ips, etype, xyz, mat, cfg, n, out)
}

/// [`stat_elec_bulk_local`] for `nd >= n` basis functions (the nodal ones
/// first, then a p element's edge, face and bubble functions) on the given
/// rule: `LocalMatrix(Element, n, nd, ...)` before `CondensateP`, material
/// values interpolated with the nodal functions only. `out` is `nd x nd`.
#[allow(clippy::too_many_arguments)]
pub fn stat_elec_bulk_local_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &StatElecAtNodes,
    cfg: &StatElecConfig,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("relative permittivity", n, mat.relative_permittivity.len())?;
    check_opt("charge density", n, &mat.charge_density)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!(
            "dim {} not in 1..=3",
            cfg.dim
        )));
    }
    if nd < n {
        return Err(ScalarError::Invalid(format!("{nd} dofs for {n} nodes")));
    }
    bulk_with(basis, ips, etype, xyz, mat, cfg, nd, out)
}

#[allow(clippy::too_many_arguments)]
fn bulk_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &StatElecAtNodes,
    cfg: &StatElecConfig,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    reset(out, nd);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(ScalarError::Basis)?;
        if be.basis.len() < nd {
            return Err(ScalarError::Invalid(format!(
                "{} basis functions for {nd} dofs",
                be.basis.len()
            )));
        }
        let phi = &be.basis[..nd];
        let dphi = &be.dbasis_dx[..nd];
        let nphi = &phi[..n];
        let w = weight(cfg, qp.weight * be.det_j, nphi, xyz);

        let eps = cfg.eps0 * interp(nphi, &mat.relative_permittivity);
        for p in 0..nd {
            for q in 0..nd {
                let s: f64 = (0..3).map(|i| dphi[p][i] * dphi[q][i]).sum();
                out.stiff[p * nd + q] += w * eps * s;
            }
        }
        if let Some(rho) = &mat.charge_density {
            let f = interp(nphi, rho);
            for p in 0..nd {
                out.force[p] += w * f * phi[p];
            }
        }
    }
    Ok(())
}

/// Boundary element matrices of `LocalMatrixBC` (StatElecSolve.F90:867) for
/// the given flux and the electrostatic layer: with `alpha = eps_l / h` and
/// `beta = flux + alpha V + rho h / (2 eps0)`, `stiff += alpha phi_q phi_p`
/// and `force += beta phi_p`. `Electric Infinity BC` is not here: it needs
/// the outward normal, which the element basis does not provide.
pub fn stat_elec_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &StatElecBcAtNodes,
    cfg: &StatElecConfig,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_opt("electric flux", n, &bc.electric_flux)?;
    check_opt(
        "layer relative permittivity",
        n,
        &bc.layer_relative_permittivity,
    )?;
    check_opt("layer thickness", n, &bc.layer_thickness)?;
    check_opt("electrode potential", n, &bc.electrode_potential)?;
    check_opt("layer charge density", n, &bc.layer_charge_density)?;
    let ips = quad.gauss_points(etype, cfg.quad_points);
    boundary_with(basis, &ips, etype, xyz, bc, cfg, n, out)
}

/// [`stat_elec_boundary_local`] for `nd >= n` basis functions on the given
/// rule (`LocalMatrixBC(Element, n, nd, ...)`); `out` is `nd x nd`.
#[allow(clippy::too_many_arguments)]
pub fn stat_elec_boundary_local_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &StatElecBcAtNodes,
    cfg: &StatElecConfig,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_opt("electric flux", n, &bc.electric_flux)?;
    check_opt(
        "layer relative permittivity",
        n,
        &bc.layer_relative_permittivity,
    )?;
    check_opt("layer thickness", n, &bc.layer_thickness)?;
    check_opt("electrode potential", n, &bc.electrode_potential)?;
    check_opt("layer charge density", n, &bc.layer_charge_density)?;
    if nd < n {
        return Err(ScalarError::Invalid(format!("{nd} dofs for {n} nodes")));
    }
    boundary_with(basis, ips, etype, xyz, bc, cfg, nd, out)
}

#[allow(clippy::too_many_arguments)]
fn boundary_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &StatElecBcAtNodes,
    cfg: &StatElecConfig,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    reset(out, nd);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(ScalarError::Basis)?;
        if be.basis.len() < nd {
            return Err(ScalarError::Invalid(format!(
                "{} basis functions for {nd} dofs",
                be.basis.len()
            )));
        }
        let all = &be.basis[..nd];
        let phi = &all[..n];
        let w = weight(cfg, qp.weight * be.det_j, phi, xyz);

        let mut alpha = 0.0;
        let mut beta = bc.electric_flux.as_deref().map_or(0.0, |f| interp(phi, f));
        let mut got_some = bc.electric_flux.is_some();
        if let Some(eps_l) = &bc.layer_relative_permittivity {
            got_some = true;
            let h = bc.layer_thickness.as_deref().ok_or_else(|| {
                ScalarError::Invalid("Charge > Layer thickness < not given!".into())
            })?;
            let h = interp(phi, h);
            let v = bc
                .electrode_potential
                .as_deref()
                .map_or(0.0, |v| interp(phi, v));
            let rho = bc
                .layer_charge_density
                .as_deref()
                .map_or(0.0, |r| interp(phi, r));
            alpha = interp(phi, eps_l) / h;
            beta += alpha * v + 0.5 * rho * h / cfg.eps0;
        }
        if !got_some {
            continue;
        }
        for p in 0..nd {
            for q in 0..nd {
                out.stiff[p * nd + q] += w * alpha * all[q] * all[p];
            }
            out.force[p] += w * beta * all[p];
        }
    }
    Ok(())
}
