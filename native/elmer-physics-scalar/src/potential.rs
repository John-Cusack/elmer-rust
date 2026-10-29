//! The potential-field recovery shared by StatElecSolve and StatCurrentSolve:
//! `LocalPostAssembly` of both modules (StatElecSolve.F90:1206,
//! StatCurrentSolve.F90:1028) and their dense element solve.

use std::f64::consts::PI;

use elmer_core::{BasisEval, ElementBasis, ElementType, Quadrature};

use crate::common::{check_len, interp};
use crate::error::ScalarError;

/// A material coefficient at the element nodes, by the shape of its keyword
/// (`ListGetElementReal`'s `Rdim`): one value, a diagonal or a full tensor.
#[derive(Clone, Debug, PartialEq)]
pub enum Coefficient {
    Scalar(Vec<f64>),
    Diagonal(Vec<[f64; 3]>),
    Full(Vec<[[f64; 3]; 3]>),
}

impl Coefficient {
    pub fn len(&self) -> usize {
        match self {
            Coefficient::Scalar(v) => v.len(),
            Coefficient::Diagonal(v) => v.len(),
            Coefficient::Full(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `C grad`, restricted to the first `dim` components, at the point with
    /// basis values `phi`.
    pub(crate) fn times(&self, phi: &[f64], grad: &[f64; 3], dim: usize) -> [f64; 3] {
        let mut out = [0.0; 3];
        match self {
            Coefficient::Scalar(c) => {
                let c = interp(phi, c);
                for i in 0..dim {
                    out[i] = c * grad[i];
                }
            }
            Coefficient::Diagonal(c) => {
                for i in 0..dim {
                    out[i] = phi.iter().zip(c).map(|(b, c)| b * c[i]).sum::<f64>() * grad[i];
                }
            }
            Coefficient::Full(c) => {
                for i in 0..dim {
                    for j in 0..dim {
                        out[i] +=
                            phi.iter().zip(c).map(|(b, c)| b * c[i][j]).sum::<f64>() * grad[j];
                    }
                }
            }
        }
        out
    }

    /// `int`-ready stiffness integrand `grad phi_p . C grad phi_q`: the first
    /// `dim` components for a tensor, all three for a scalar (`MATMUL`).
    pub(crate) fn stiffness(
        &self,
        phi: &[f64],
        dphi: &[[f64; 3]],
        p: usize,
        q: usize,
        dim: usize,
    ) -> f64 {
        match self {
            Coefficient::Scalar(c) => {
                let s: f64 = (0..3).map(|i| dphi[p][i] * dphi[q][i]).sum();
                interp(phi, c) * s
            }
            Coefficient::Diagonal(c) => (0..dim)
                .map(|i| {
                    phi.iter().zip(c).map(|(b, c)| b * c[i]).sum::<f64>() * dphi[p][i] * dphi[q][i]
                })
                .sum(),
            Coefficient::Full(c) => {
                let mut a = 0.0;
                for i in 0..dim {
                    for j in 0..dim {
                        a += phi.iter().zip(c).map(|(b, c)| b * c[i][j]).sum::<f64>()
                            * dphi[p][i]
                            * dphi[q][j];
                    }
                }
                a
            }
        }
    }
}

/// Element-level switches for [`potential_post_local`].
#[derive(Clone, Debug, PartialEq)]
pub struct PostConfig {
    pub dim: usize,
    /// The measure becomes `2 pi r dr dz`.
    pub axisymmetric: bool,
    pub quad_points: Option<usize>,
}

/// One element's contribution to the postprocessing fields. The rows are the
/// right-hand sides `int f phi_i` of the nodal projection of each field
/// against `mass`.
#[derive(Clone, Debug, Default)]
pub struct PotentialPostLocal {
    pub n: usize,
    /// Consistent mass matrix, row-major `n x n`.
    pub mass: Vec<f64>,
    /// `int phi_i`: the nodal weights.
    pub weights: Vec<f64>,
    /// `int E.J phi_i`: twice the energy density, or the Joule heating.
    pub heat: Vec<f64>,
    /// `int -C grad(V)_j phi_i` for j in 0..dim (rest zero).
    pub flux: [Vec<f64>; 3],
    /// `int -grad(V)_j phi_i` for j in 0..dim (rest zero).
    pub field: [Vec<f64>; 3],
    /// `int dV`.
    pub volume: f64,
    /// `int E.J dV`.
    pub heat_integral: f64,
    /// `int V dV`.
    pub potential_integral: f64,
}

/// The postprocessing integrals of one element for the potential `pot` at its
/// nodes and the coefficient `coef` (relative permittivity or electric
/// conductivity) times `scale` (the permittivity of vacuum, or 1).
#[allow(clippy::too_many_arguments)]
pub fn potential_post_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    pot: &[f64],
    coef: &Coefficient,
    scale: f64,
    cfg: &PostConfig,
    out: &mut PotentialPostLocal,
) -> Result<(), ScalarError> {
    let n = xyz.len();
    check_len("potential", n, pot.len())?;
    check_len("coefficient", n, coef.len())?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(ScalarError::Invalid(format!(
            "dim {} not in 1..=3",
            cfg.dim
        )));
    }
    let dim = cfg.dim;
    out.n = n;
    out.mass.clear();
    out.mass.resize(n * n, 0.0);
    for v in [&mut out.weights, &mut out.heat] {
        v.clear();
        v.resize(n, 0.0);
    }
    for v in out.flux.iter_mut().chain(out.field.iter_mut()) {
        v.clear();
        v.resize(n, 0.0);
    }
    out.volume = 0.0;
    out.heat_integral = 0.0;
    out.potential_integral = 0.0;

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
        let mut w = qp.weight * be.det_j;
        if cfg.axisymmetric {
            w *= 2.0 * PI * phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>();
        }

        for i in 0..n {
            for j in 0..n {
                out.mass[i * n + j] += w * phi[i] * phi[j];
            }
            out.weights[i] += w * phi[i];
        }
        let mut grad = [0.0; 3];
        for (j, g) in grad.iter_mut().enumerate().take(dim) {
            *g = (0..n).map(|k| dphi[k][j] * pot[k]).sum();
        }
        let mut cgrad = coef.times(phi, &grad, dim);
        cgrad.iter_mut().for_each(|c| *c *= scale);
        let heat: f64 = (0..dim).map(|j| grad[j] * cgrad[j]).sum();
        for j in 0..dim {
            for i in 0..n {
                out.field[j][i] -= grad[j] * w * phi[i];
                out.flux[j][i] -= cgrad[j] * w * phi[i];
            }
        }
        for i in 0..n {
            out.heat[i] += heat * w * phi[i];
        }
        out.volume += w;
        out.heat_integral += w * heat;
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
        // LAPACK order: all row interchanges first, then the unit lower factor
        // (whose rows carry the later interchanges).
        for (c, &p) in self.piv.iter().enumerate() {
            x.swap(c, p);
        }
        for c in 0..n {
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
