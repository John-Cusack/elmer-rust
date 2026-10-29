//! StatElecSolve element kernels: `LocalMatrix`, `LocalMatrixBC` and
//! `LocalPostAssembly` of `fem/src/modules/StatElecSolve.F90` for the real
//! (single potential) path.

use std::f64::consts::PI;

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

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
        return Err(ScalarError::Invalid(format!("dim {} not in 1..=3", cfg.dim)));
    }
    reset(out, n);
    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        basis.eval(etype, xyz, qp.uvw, &mut be).map_err(ScalarError::Basis)?;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];
        let w = weight(cfg, qp.weight * be.det_j, phi, xyz);

        let eps = cfg.eps0 * interp(phi, &mat.relative_permittivity);
        for p in 0..n {
            for q in 0..n {
                let s: f64 = (0..3).map(|i| dphi[p][i] * dphi[q][i]).sum();
                out.stiff[p * n + q] += w * eps * s;
            }
        }
        if let Some(rho) = &mat.charge_density {
            let f = interp(phi, rho);
            for p in 0..n {
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
    check_opt("layer relative permittivity", n, &bc.layer_relative_permittivity)?;
    check_opt("layer thickness", n, &bc.layer_thickness)?;
    check_opt("electrode potential", n, &bc.electrode_potential)?;
    check_opt("layer charge density", n, &bc.layer_charge_density)?;
    reset(out, n);
    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        basis.eval(etype, xyz, qp.uvw, &mut be).map_err(ScalarError::Basis)?;
        let phi = &be.basis[..n];
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
            let v = bc.electrode_potential.as_deref().map_or(0.0, |v| interp(phi, v));
            let rho = bc.layer_charge_density.as_deref().map_or(0.0, |r| interp(phi, r));
            alpha = interp(phi, eps_l) / h;
            beta += alpha * v + 0.5 * rho * h / cfg.eps0;
        }
        if !got_some {
            continue;
        }
        for p in 0..n {
            for q in 0..n {
                out.stiff[p * n + q] += w * alpha * phi[q] * phi[p];
            }
            out.force[p] += w * beta * phi[p];
        }
    }
    Ok(())
}

/// One element's contribution to the postprocessing fields
/// (`LocalPostAssembly`, StatElecSolve.F90:1206). The `*_rows` are the
/// right-hand sides `int f phi_i` of the nodal projection of each field
/// against `mass`.
#[derive(Clone, Debug, Default)]
pub struct StatElecPostLocal {
    pub n: usize,
    /// Consistent mass matrix, row-major `n x n`.
    pub mass: Vec<f64>,
    /// `int phi_i`: the nodal weights.
    pub weights: Vec<f64>,
    /// `int E.D phi_i`, the energy density times two.
    pub energy: Vec<f64>,
    /// `int -eps0 eps_r grad(V)_j phi_i` for j in 0..dim (rest zero).
    pub flux: [Vec<f64>; 3],
    /// `int -grad(V)_j phi_i` for j in 0..dim (rest zero).
    pub field: [Vec<f64>; 3],
    /// `int dV`.
    pub volume: f64,
    /// `int E.D dV`.
    pub energy_integral: f64,
    /// `int V dV`.
    pub potential_integral: f64,
}

/// The postprocessing integrals of one element for the potential `pot` at its
/// nodes and the relative permittivity `eps_r` at its nodes.
pub fn stat_elec_post_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    pot: &[f64],
    eps_r: &[f64],
    cfg: &StatElecConfig,
    out: &mut StatElecPostLocal,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("potential", n, pot.len())?;
    check_len("relative permittivity", n, eps_r.len())?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!("dim {} not in 1..=3", cfg.dim)));
    }
    let dim = cfg.dim;
    out.n = n;
    out.mass.clear();
    out.mass.resize(n * n, 0.0);
    for v in [&mut out.weights, &mut out.energy] {
        v.clear();
        v.resize(n, 0.0);
    }
    for v in out.flux.iter_mut().chain(out.field.iter_mut()) {
        v.clear();
        v.resize(n, 0.0);
    }
    out.volume = 0.0;
    out.energy_integral = 0.0;
    out.potential_integral = 0.0;

    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(ScalarError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        basis.eval(etype, xyz, qp.uvw, &mut be).map_err(ScalarError::Basis)?;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];
        let w = weight(cfg, qp.weight * be.det_j, phi, xyz);

        for i in 0..n {
            for j in 0..n {
                out.mass[i * n + j] += w * phi[i] * phi[j];
            }
            out.weights[i] += w * phi[i];
        }
        let eps = cfg.eps0 * interp(phi, eps_r);
        let mut grad = [0.0; 3];
        for (j, g) in grad.iter_mut().enumerate().take(dim) {
            *g = (0..n).map(|k| dphi[k][j] * pot[k]).sum();
        }
        let heat: f64 = (0..dim).map(|j| grad[j] * eps * grad[j]).sum();
        for j in 0..dim {
            for i in 0..n {
                out.field[j][i] -= grad[j] * w * phi[i];
                out.flux[j][i] -= eps * grad[j] * w * phi[i];
            }
        }
        for i in 0..n {
            out.energy[i] += heat * w * phi[i];
        }
        out.volume += w;
        out.energy_integral += w * heat;
        out.potential_integral += w * interp(phi, pot);
    }
    Ok(())
}

/// `LUdecomp` of the small element mass matrix (partial pivoting), for
/// repeated `LUSolve`s against different right-hand sides.
pub struct DenseLu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
}

impl DenseLu {
    /// Factors the row-major `n x n` matrix `a`.
    pub fn new(mut a: Vec<f64>, n: usize) -> Result<Self, ScalarError> {
        let mut piv = vec![0usize; n];
        for c in 0..n {
            let p = (c..n)
                .max_by(|&x, &y| a[x * n + c].abs().total_cmp(&a[y * n + c].abs()))
                .expect("non-empty range");
            if a[p * n + c] == 0.0 {
                return Err(ScalarError::Invalid("LU-decomposition fails".into()));
            }
            piv[c] = p;
            if p != c {
                for j in 0..n {
                    a.swap(c * n + j, p * n + j);
                }
            }
            for i in c + 1..n {
                let f = a[i * n + c] / a[c * n + c];
                a[i * n + c] = f;
                for j in c + 1..n {
                    a[i * n + j] -= f * a[c * n + j];
                }
            }
        }
        Ok(Self { n, lu: a, piv })
    }

    /// Solves `A x = b` in place.
    pub fn solve(&self, x: &mut [f64]) {
        let n = self.n;
        for c in 0..n {
            x.swap(c, self.piv[c]);
            for i in c + 1..n {
                x[i] -= self.lu[i * n + c] * x[c];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= self.lu[i * n + j] * x[j];
            }
            x[i] /= self.lu[i * n + i];
        }
    }
}
