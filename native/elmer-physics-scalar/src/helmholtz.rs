//! HelmholtzSolve element kernels: `LocalMatrix` and `LocalMatrixBoundary` of
//! `fem/src/modules/HelmholtzSolve.F90` without bubbles, convection and the
//! preconditioning matrix.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};
use num_complex::Complex64;

use crate::common::{check_len, interp};
use crate::error::ScalarError;

/// `AEPS` (Types.F90): the impedance below which the admittance is zero.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// Element-level switches for the Helmholtz kernels.
#[derive(Clone, Debug, PartialEq)]
pub struct HelmholtzConfig {
    /// Coordinate-system dimension (`CoordinateSystemDimension()`), 1 to 3.
    pub dim: usize,
    /// Axisymmetric or cylindric coordinates: the measure becomes `r dr dz`
    /// (`SqrtMetric`) with `r` the first coordinate.
    pub axisymmetric: bool,
    /// Angular frequency in rad/s.
    pub angular_frequency: f64,
    /// Requested Gauss points, or `None` for the element default.
    pub quad_points: Option<usize>,
}

impl HelmholtzConfig {
    pub fn new(dim: usize, angular_frequency: f64) -> Self {
        Self {
            dim,
            axisymmetric: false,
            angular_frequency,
            quad_points: None,
        }
    }
}

/// Bulk properties at the element nodes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HelmholtzBulkAtNodes {
    pub sound_speed: Vec<f64>,
    pub density: Vec<f64>,
    /// `Sound Damping`; absent is zero.
    pub damping: Option<Vec<f64>>,
    /// `Pressure Source 1` + i `Pressure Source 2`; absent is zero.
    pub load: Option<Vec<Complex64>>,
}

/// Boundary properties at the nodes of a boundary element.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HelmholtzBcAtNodes {
    /// The density of the parent element's material.
    pub density: Vec<f64>,
    /// `Wave Impedance 1` + i `Wave Impedance 2`; absent or zero is a
    /// vanishing admittance.
    pub impedance: Option<Vec<Complex64>>,
    /// `Wave Flux 1` + i `Wave Flux 2`.
    pub wave_flux: Option<Vec<Complex64>>,
}

/// A complex element system, row-major `n x n` with the test function on the
/// row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HelmholtzLocal {
    pub n: usize,
    pub stiff: Vec<Complex64>,
    pub force: Vec<Complex64>,
}

impl HelmholtzLocal {
    fn reset(&mut self, n: usize) {
        self.n = n;
        self.stiff.clear();
        self.stiff.resize(n * n, Complex64::new(0.0, 0.0));
        self.force.clear();
        self.force.resize(n, Complex64::new(0.0, 0.0));
    }

    /// The equivalent real system with each node's real and imaginary parts
    /// interleaved, as `DefaultUpdateEquations` stores a complex local system
    /// in a matrix with two dofs per node: `[[Re, -Im], [Im, Re]]` per entry.
    pub fn to_real(&self, out: &mut LocalSystem) {
        let n = self.n;
        let m = 2 * n;
        out.n = m;
        out.stiff.clear();
        out.stiff.resize(m * m, 0.0);
        out.force.clear();
        out.force.resize(m, 0.0);
        out.mass = None;
        out.damp = None;
        for p in 0..n {
            for q in 0..n {
                let a = self.stiff[p * n + q];
                out.stiff[2 * p * m + 2 * q] = a.re;
                out.stiff[2 * p * m + 2 * q + 1] = -a.im;
                out.stiff[(2 * p + 1) * m + 2 * q] = a.im;
                out.stiff[(2 * p + 1) * m + 2 * q + 1] = a.re;
            }
            out.force[2 * p] = self.force[p].re;
            out.force[2 * p + 1] = self.force[p].im;
        }
    }
}

fn check_complex(
    what: &'static str,
    n: usize,
    v: &Option<Vec<Complex64>>,
) -> Result<(), ScalarError> {
    v.as_ref().map_or(Ok(()), |v| check_len(what, n, v.len()))
}

fn interp_c(phi: &[f64], nodal: &[Complex64]) -> Complex64 {
    phi.iter().zip(nodal).map(|(b, v)| *v * *b).sum()
}

fn measure(cfg: &HelmholtzConfig, w: f64, phi: &[f64], xyz: &[[f64; 3]]) -> f64 {
    if cfg.axisymmetric {
        w * phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>()
    } else {
        w
    }
}

/// Bulk element system of `LocalMatrix` (HelmholtzSolve.F90:497) for an
/// element with `nd = n` nodes and no bubbles: with `k = omega / c`,
/// `stiff[p*n + q] = int ((-k^2 + i k d) phi_q phi_p + grad phi_q . grad phi_p) / rho`
/// and `force[p] = int phi_p (f1 + i f2) / rho`.
pub fn helmholtz_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &HelmholtzBulkAtNodes,
    cfg: &HelmholtzConfig,
    out: &mut HelmholtzLocal,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("sound speed", n, mat.sound_speed.len())?;
    check_len("density", n, mat.density.len())?;
    if let Some(d) = &mat.damping {
        check_len("sound damping", n, d.len())?;
    }
    check_complex("pressure source", n, &mat.load)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!(
            "dim {} not in 1..=3",
            cfg.dim
        )));
    }
    out.reset(n);
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
        let s = measure(cfg, qp.weight * be.det_j, phi, xyz);

        let wave_number = cfg.angular_frequency / interp(phi, &mat.sound_speed);
        let d = wave_number * mat.damping.as_deref().map_or(0.0, |d| interp(phi, d));
        let m = -wave_number * wave_number;
        let load = mat
            .load
            .as_deref()
            .map_or(Complex64::new(0.0, 0.0), |l| interp_c(phi, l));
        let rho = interp(phi, &mat.density);

        for p in 0..n {
            for q in 0..n {
                let grad: f64 = (0..cfg.dim).map(|i| dphi[q][i] * dphi[p][i]).sum();
                let a = Complex64::new(m, d) * (phi[q] * phi[p]) + grad;
                out.stiff[p * n + q] += a * (s / rho);
            }
            out.force[p] += load * (s * phi[p] / rho);
        }
    }
    Ok(())
}

/// Boundary element system of `LocalMatrixBoundary` (HelmholtzSolve.F90:640)
/// without convection: `stiff[p*n + q] = int Y phi_q phi_p / rho` with the
/// admittance `Y = i omega / (Z1 + i Z2)` and `force[p] = int phi_p (q1 + i q2) / rho`.
pub fn helmholtz_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &HelmholtzBcAtNodes,
    cfg: &HelmholtzConfig,
    out: &mut HelmholtzLocal,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("density", n, bc.density.len())?;
    check_complex("wave impedance", n, &bc.impedance)?;
    check_complex("wave flux", n, &bc.wave_flux)?;
    out.reset(n);
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
        let s = measure(cfg, qp.weight * be.det_j, phi, xyz);

        let z = bc
            .impedance
            .as_deref()
            .map_or(Complex64::new(0.0, 0.0), |z| interp_c(phi, z));
        let admittance = if z.re.abs() < AEPS && z.im.abs() < AEPS {
            Complex64::new(0.0, 0.0)
        } else {
            Complex64::new(0.0, cfg.angular_frequency) / z
        };
        let rho = interp(phi, &bc.density);
        let load = bc
            .wave_flux
            .as_deref()
            .map_or(Complex64::new(0.0, 0.0), |l| interp_c(phi, l));
        for p in 0..n {
            for q in 0..n {
                out.stiff[p * n + q] += admittance * (s * phi[q] * phi[p] / rho);
            }
            out.force[p] += load * (s * phi[p] / rho);
        }
    }
    Ok(())
}
