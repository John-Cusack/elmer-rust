use elmer_core::{BasisEval, ElementBasis, ElementType, Quadrature};

use crate::error::FluxError;

/// Which fields `FluxSolver` projects and how the element is integrated.
///
/// The right-hand sides come out in Elmer's order: `Flux` (`dim` rows), `Grad`
/// (`dim` rows), `Flux_mag`, `Grad_mag`; only the requested ones are present.
#[derive(Clone, Debug, PartialEq)]
pub struct FluxConfig {
    /// Coordinate system dimension, 2 or 3.
    pub dim: usize,
    /// `Coordinate System = Axi Symmetric` or `Cylindric Symmetric`: the
    /// integration weight is multiplied by the radius `x`.
    pub axisymmetric: bool,
    pub calculate_flux: bool,
    pub calculate_flux_magnitude: bool,
    pub calculate_grad: bool,
    pub calculate_grad_magnitude: bool,
    /// Gauss points per element; `None` is the element type's default
    /// (`GaussPoints(Element)`).
    pub quad_points: Option<usize>,
}

impl FluxConfig {
    /// Nothing requested yet.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            axisymmetric: false,
            calculate_flux: false,
            calculate_flux_magnitude: false,
            calculate_grad: false,
            calculate_grad_magnitude: false,
            quad_points: None,
        }
    }

    /// Number of projected scalar fields, `Dofs` in FluxSolver.F90:96.
    pub fn dofs(&self) -> usize {
        let mut d = 0;
        if self.calculate_flux {
            d += self.dim;
        }
        if self.calculate_flux_magnitude {
            d += 1;
        }
        if self.calculate_grad {
            d += self.dim;
        }
        if self.calculate_grad_magnitude {
            d += 1;
        }
        d
    }

    /// Index of the first magnitude right-hand side (`firstmag - 1`).
    pub fn first_magnitude(&self) -> usize {
        let mut i = 0;
        if self.calculate_flux {
            i += self.dim;
        }
        if self.calculate_grad {
            i += self.dim;
        }
        i
    }

    fn needs_flux(&self) -> bool {
        self.calculate_flux || self.calculate_flux_magnitude
    }
}

/// The `Flux Coefficient` at the element nodes, by the tensor rank
/// `GetTensorRank` finds in the material array (FluxSolver.F90:553).
#[derive(Clone, Debug, PartialEq)]
pub enum FluxCoefficient {
    /// The keyword is absent: the flux is the gradient (rank 0, `C = I`).
    Identity,
    /// A scalar (rank 1) times the identity.
    Scalar(Vec<f64>),
    /// The diagonal of `C` (rank 2); entries past `dim` are ignored.
    Diagonal(Vec<[f64; 3]>),
    /// A full tensor (rank 3), `[row][col]`.
    Full(Vec<[[f64; 3]; 3]>),
}

/// One element's projection system: the mass matrix `STIFF` and one
/// right-hand side per projected field.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FluxLocal {
    /// Number of element nodes.
    pub n: usize,
    /// Row-major `n x n`, `mass[p*n + q] = sum w phi_q phi_p`.
    pub mass: Vec<f64>,
    /// `force[field][p]`, `cfg.dofs()` rows of `n`.
    pub force: Vec<Vec<f64>>,
}

fn interp(phi: &[f64], v: &[f64]) -> f64 {
    phi.iter().zip(v).map(|(a, b)| a * b).sum()
}

fn check_len(what: &'static str, n: usize, got: usize) -> Result<(), FluxError> {
    if got == n {
        Ok(())
    } else {
        Err(FluxError::NodalLength {
            what,
            expected: n,
            got,
        })
    }
}

/// Bulk element matrices of `FluxSolver`'s `BulkAssembly` (FluxSolver.F90:289)
/// for an element with `nd = n` nodes (no bubble or p dofs).
///
/// `potential` holds the target variable at the element nodes
/// (`LocalPotential`); the gradient at an integration point is
/// `sum_p potential[p] * dphi[p]`. The flux is `-C grad` and the stored right-hand
/// side is `-w (C grad)_i phi`, as in the Fortran, so the projected `Flux`
/// field is the physical flux. Gradient and magnitude right-hand sides carry
/// `+w`.
///
/// `out` is reset to `n x n` and `cfg.dofs()` rows.
#[allow(clippy::too_many_arguments)] // the lane's kernel signature: basis, quadrature, geometry, inputs, config, output
pub fn flux_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    potential: &[f64],
    coeff: &FluxCoefficient,
    cfg: &FluxConfig,
    out: &mut FluxLocal,
) -> Result<(), FluxError> {
    let n = xyz.len();
    check_len("potential", n, potential.len())?;
    match coeff {
        FluxCoefficient::Identity => {}
        FluxCoefficient::Scalar(v) => check_len("flux coefficient", n, v.len())?,
        FluxCoefficient::Diagonal(v) => check_len("flux coefficient", n, v.len())?,
        FluxCoefficient::Full(v) => check_len("flux coefficient", n, v.len())?,
    }
    if !(2..=3).contains(&cfg.dim) {
        return Err(FluxError::Invalid(format!(
            "dim {} not in 2..=3 (flux computation makes sense only in 2D and 3D)",
            cfg.dim
        )));
    }
    let dim = cfg.dim;
    let dofs = cfg.dofs();

    out.n = n;
    out.mass.clear();
    out.mass.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(dofs, vec![0.0; n]);

    let ips = quad.gauss_points(etype, cfg.quad_points);
    if ips.is_empty() {
        return Err(FluxError::NoQuadrature(etype.0));
    }

    let mut be = BasisEval::default();
    for qp in &ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(FluxError::Basis)?;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];

        let mut w = qp.weight * be.det_j;
        if cfg.axisymmetric {
            w *= phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>();
        }

        for p in 0..n {
            for q in 0..n {
                out.mass[p * n + q] += w * phi[q] * phi[p];
            }
        }

        let mut grad = [0.0; 3];
        for (i, g) in grad.iter_mut().enumerate().take(dim) {
            *g = potential.iter().zip(dphi).map(|(u, d)| u * d[i]).sum();
        }

        let mut flux = [0.0; 3];
        if cfg.needs_flux() {
            let mut c = [[0.0; 3]; 3];
            for (i, row) in c.iter_mut().enumerate().take(dim) {
                row[i] = 1.0;
            }
            match coeff {
                FluxCoefficient::Identity => {}
                FluxCoefficient::Scalar(k) => {
                    let kt = interp(phi, k);
                    for (i, row) in c.iter_mut().enumerate().take(dim) {
                        row[i] = kt;
                    }
                }
                FluxCoefficient::Diagonal(k) => {
                    for (i, row) in c.iter_mut().enumerate().take(dim) {
                        row[i] = phi.iter().zip(k).map(|(b, kn)| b * kn[i]).sum();
                    }
                }
                FluxCoefficient::Full(k) => {
                    for (i, row) in c.iter_mut().enumerate().take(dim) {
                        for (j, cij) in row.iter_mut().enumerate().take(dim) {
                            *cij = phi.iter().zip(k).map(|(b, kn)| b * kn[i][j]).sum();
                        }
                    }
                }
            }
            for i in 0..dim {
                flux[i] = (0..dim).map(|j| c[i][j] * grad[j]).sum();
            }
        }

        let mut k = 0;
        if cfg.calculate_flux {
            for i in 0..dim {
                let c = w * flux[i];
                for p in 0..n {
                    out.force[i][p] -= c * phi[p];
                }
            }
            k += dim;
        }
        if cfg.calculate_grad {
            for i in 0..dim {
                let c = w * grad[i];
                for p in 0..n {
                    out.force[k + i][p] += c * phi[p];
                }
            }
            k += dim;
        }
        if cfg.calculate_flux_magnitude {
            let c = w * flux.iter().map(|f| f * f).sum::<f64>().sqrt();
            for p in 0..n {
                out.force[k][p] += c * phi[p];
            }
            k += 1;
        }
        if cfg.calculate_grad_magnitude {
            let c = w * grad.iter().map(|g| g * g).sum::<f64>().sqrt();
            for p in 0..n {
                out.force[k][p] += c * phi[p];
            }
        }
    }
    Ok(())
}
