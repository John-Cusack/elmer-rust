//! `DiffuseConvectiveCompose` and `DiffuseConvectiveBoundary`
//! (AdvectionDiffusion.F90:1135, :1668): the Cartesian element matrices of
//! `AdvectionDiffusionSolver`.

use elmer_core::{BasisEval, ElementType, LocalSystem};

use crate::error::{check_len, interp, AdvDiffError};
use crate::fe::{bubble_gauss_points, eval_with_bubbles, Fe};

/// Elmer's `AEPS`, the threshold below which a Soret or Nernst-Planck
/// coefficient counts as absent.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// How the convection term is stabilised.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Stabilization {
    /// Plain Galerkin.
    None,
    /// Residual-free bubbles (`Bubbles = True`, the Elmer default): the
    /// element gains one bubble function per node, condensed out by
    /// [`crate::condensate`] after the time terms. Active only where the
    /// element has convection.
    Bubbles,
    /// SUPG (`Stabilize = True`) with the element size `Element % hK` and
    /// `Element % StabilizationMK` (`StabParam`, Stabilize.F90), which the
    /// element layer provides. Needs [`Fe::hessian`].
    Supg { h_k: f64, m_k: f64 },
}

/// Element-level switches of [`advdiff_bulk_local`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AdvDiffConfig {
    /// `CoordinateSystemDimension()`, 1 to 3.
    pub dim: usize,
    pub stabilization: Stabilization,
}

impl AdvDiffConfig {
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            stabilization: Stabilization::None,
        }
    }
}

/// Nodal coefficient values of one bulk element, in the notation of the
/// Fortran (the users' guide's `C_T`, `C_1`, `C_2`):
/// relative mass units give `ct = c1 = density` and `c2 = density * D`;
/// absolute mass units give `ct = 1`, `c1 = 1` (when there is convection or
/// migration) and `c2 = D`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdvDiffNodal {
    /// Time derivative coefficient `C_T`.
    pub ct: Vec<f64>,
    /// Convection coefficient `C_1`; all zero turns convection off.
    pub c1: Vec<f64>,
    /// Diffusion tensor `C_2`.
    pub c2: Vec<[[f64; 3]; 3]>,
    /// Convection velocity `U` (`Ux, Uy, Uz` of the Fortran).
    pub velocity: Vec<[f64; 3]>,
    /// Mesh velocity, subtracted from the convection velocity; empty for none.
    pub mesh_velocity: Vec<[f64; 3]>,
    /// Source `Load`, the `<name> Diffusion Source`.
    pub load: Vec<f64>,
    /// `NodalC0` for the Soret term: `density` in relative mass units, `1` in
    /// absolute. The zero-order term of the Fortran is always switched off
    /// (`NodalC0 = 0` after this is copied), so this is its only use.
    pub c0: Vec<f64>,
    /// `<name> Soret Diffusivity`; empty for none.
    pub soret_diffusivity: Vec<f64>,
    /// Temperature, for the Soret flux `grad T`.
    pub temperature: Vec<f64>,
    /// Nernst-Einstein constant `NEConst`; empty for no migration.
    pub ne_const: Vec<f64>,
    /// Electric potential, for the migration velocity `grad phi`.
    pub potential: Vec<f64>,
}

/// Bulk element matrices of `DiffuseConvectiveCompose` for a scalar species
/// concentration in Cartesian coordinates.
///
/// `out` is reset to `nb x nb` (`nb = n`, or `2n` when [`Stabilization::Bubbles`]
/// is in effect and the element has convection): `stiff` holds diffusion,
/// convection, migration and stabilisation, `mass` the `C_T` mass matrix
/// (stabilised for SUPG), `force` the source and Soret terms. Matrices are
/// row-major with the test function on the row.
///
/// Not implemented, as in [`crate::check_supported`]: compressible flow (the
/// `div u` zero-order term), which the preflight refuses.
///
/// If a point has no convection (`C1 = 0` there) under SUPG, the stabilisation
/// parameter and residual of the previous point stay in force for that point's
/// right-hand side, as the Fortran's variables do.
pub fn advdiff_bulk_local(
    fe: &Fe<'_>,
    etype: ElementType,
    xyz: &[[f64; 3]],
    nodal: &AdvDiffNodal,
    cfg: &AdvDiffConfig,
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    check_nodal(n, nodal)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(AdvDiffError::Invalid(format!("dim {} not in 1..=3", cfg.dim)));
    }
    let dim = cfg.dim;

    let any_convection = nodal.c1.iter().any(|&c| c != 0.0);
    let supg = match cfg.stabilization {
        Stabilization::Supg { h_k, m_k } => Some((h_k, m_k)),
        _ => None,
    };
    let bubbles = any_convection && cfg.stabilization == Stabilization::Bubbles;
    let nb = if bubbles { 2 * n } else { n };

    let thermal = nodal.soret_diffusivity.iter().any(|d| d.abs() > AEPS);
    let migration = nodal.ne_const.iter().any(|d| d.abs() > AEPS);
    let stabilized = supg.is_some() && (any_convection || migration);

    out.n = nb;
    out.stiff.clear();
    out.stiff.resize(nb * nb, 0.0);
    out.force.clear();
    out.force.resize(nb, 0.0);
    out.damp = None;
    let mut mass = out.mass.take().unwrap_or_default();
    mass.clear();
    mass.resize(nb * nb, 0.0);

    let np = if bubbles {
        Some(bubble_gauss_points(etype).ok_or(AdvDiffError::NoQuadrature(etype.0))?)
    } else {
        None
    };
    let ips = fe.quad.gauss_points(etype, np);
    if ips.is_empty() {
        return Err(AdvDiffError::NoQuadrature(etype.0));
    }

    let mut be = BasisEval::default();
    let mut hess = Vec::new();
    let (mut tau, mut su) = (0.0, vec![0.0; n]);
    let mut velo = [0.0; 3];

    for qp in &ips {
        if bubbles {
            eval_with_bubbles(fe.basis, etype, xyz, qp.uvw, dim, &mut be)?;
        } else {
            fe.basis
                .eval(etype, xyz, qp.uvw, &mut be)
                .map_err(AdvDiffError::Basis)?;
        }
        let phi = &be.basis[..];
        let dphi = &be.dbasis_dx[..];
        let s = be.det_j * qp.weight;

        let c1 = interp(&phi[..n], &nodal.c1);
        let ct = interp(&phi[..n], &nodal.ct);
        let mut c2 = [[0.0; 3]; 3];
        for (i, row) in c2.iter_mut().enumerate().take(dim) {
            for (j, c) in row.iter_mut().enumerate().take(dim) {
                *c = (0..n).map(|m| nodal.c2[m][i][j] * phi[m]).sum();
            }
        }

        let convection = c1 != 0.0;
        if convection {
            for (i, v) in velo.iter_mut().enumerate().take(dim.max(2)) {
                *v = (0..n)
                    .map(|m| {
                        let mv = nodal.mesh_velocity.get(m).map_or(0.0, |x| x[i]);
                        (nodal.velocity[m][i] - mv) * phi[m]
                    })
                    .sum();
            }
            if migration {
                let ne = interp(&phi[..n], &nodal.ne_const);
                let mut grad_phi = [0.0; 3];
                for (i, g) in grad_phi.iter_mut().enumerate().take(dim) {
                    *g = (0..n).map(|m| dphi[m][i] * nodal.potential[m]).sum();
                }
                for i in 0..dim {
                    for j in 0..dim {
                        velo[i] -= c2[i][j] * ne * grad_phi[j];
                    }
                }
            }

            if let (true, Some((h_k, m_k))) = (stabilized, supg) {
                let vnorm = velo[..dim].iter().map(|v| v * v).sum::<f64>().sqrt();
                let pe = (m_k * h_k * c1 * vnorm / (2.0 * c2[0][0].abs())).min(1.0);
                tau = if vnorm != 0.0 {
                    h_k * pe / (2.0 * c1 * vnorm)
                } else {
                    0.0
                };
                let hessian = fe.hessian.ok_or_else(|| AdvDiffError::Unsupported {
                    section: "solver",
                    keyword: "Stabilize".into(),
                    reason: "SUPG needs second derivatives of the basis (ddBasisddx)",
                })?;
                hessian
                    .eval_hessian(etype, xyz, qp.uvw, &mut hess)
                    .map_err(AdvDiffError::Basis)?;
                let mut dc2dx = [[[0.0; 3]; 3]; 3];
                for (i, a) in dc2dx.iter_mut().enumerate().take(dim) {
                    for (j, b) in a.iter_mut().enumerate().take(dim) {
                        for (k, c) in b.iter_mut().enumerate().take(dim) {
                            *c = (0..n).map(|m| nodal.c2[m][i][j] * dphi[m][k]).sum();
                        }
                    }
                }
                for p in 0..n {
                    let mut r = 0.0;
                    for i in 0..dim {
                        r += c1 * dphi[p][i] * velo[i];
                        for j in 0..dim {
                            r -= c2[i][j] * hess[p][i][j];
                            r -= dc2dx[i][j][j] * dphi[p][i];
                        }
                    }
                    su[p] = r;
                }
            }
        }

        for p in 0..nb {
            for q in 0..nb {
                let mut m = ct * phi[q] * phi[p];
                let mut a = 0.0;
                for i in 0..dim {
                    for j in 0..dim {
                        a += c2[i][j] * dphi[q][i] * dphi[p][j];
                    }
                }
                if convection {
                    for i in 0..dim {
                        a += c1 * velo[i] * dphi[q][i] * phi[p];
                    }
                    if stabilized {
                        a += tau * su[q] * su[p];
                        m += tau * ct * phi[q] * su[p];
                    }
                }
                out.stiff[p * nb + q] += s * a;
                mass[p * nb + q] += s * m;
            }
        }

        let force = interp(&phi[..n], &nodal.load);
        for p in 0..nb {
            let mut load = phi[p];
            if stabilized {
                load += tau * su[p];
            }
            out.force[p] += s * force * load;
        }

        if thermal {
            let c_thermal = interp(&phi[..n], &nodal.c0);
            let sor_d = interp(&phi[..n], &nodal.soret_diffusivity);
            let mut grad_t = [0.0; 3];
            for (i, g) in grad_t.iter_mut().enumerate().take(dim) {
                *g = (0..n).map(|m| dphi[m][i] * nodal.temperature[m]).sum();
            }
            for p in 0..nb {
                let load = if stabilized { tau * su[p] } else { 1.0 };
                let soret: f64 = c_thermal
                    * sor_d
                    * (0..dim).map(|i| grad_t[i] * dphi[p][i]).sum::<f64>();
                out.force[p] -= s * soret * load;
            }
        }
    }
    out.mass = Some(mass);
    Ok(())
}

/// Boundary element matrices of `DiffuseConvectiveBoundary`
/// (AdvectionDiffusion.F90:1668): `stiff = int alpha N_p N_q`, `force =
/// int load N_q`, where `alpha` is `Mass Transfer Coefficient` and `load` is
/// `alpha * External Concentration + <name> Flux`.
///
/// `out` is reset to `n x n` with no mass.
pub fn advdiff_boundary_local(
    fe: &Fe<'_>,
    etype: ElementType,
    xyz: &[[f64; 3]],
    alpha: &[f64],
    load: &[f64],
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    check_len("mass transfer coefficient", n, alpha.len())?;
    check_len("boundary load", n, load.len())?;
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.mass = None;
    out.damp = None;

    let ips = fe.quad.gauss_points(etype, None);
    if ips.is_empty() {
        return Err(AdvDiffError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in &ips {
        fe.basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(AdvDiffError::Basis)?;
        let phi = &be.basis[..n];
        let s = be.det_j * qp.weight;
        let force = interp(phi, load);
        let a = interp(phi, alpha);
        for p in 0..n {
            for q in 0..n {
                out.stiff[p * n + q] += s * a * phi[q] * phi[p];
            }
        }
        for q in 0..n {
            out.force[q] += s * phi[q] * force;
        }
    }
    Ok(())
}

fn check_nodal(n: usize, a: &AdvDiffNodal) -> Result<(), AdvDiffError> {
    check_len("time coefficient", n, a.ct.len())?;
    check_len("convection coefficient", n, a.c1.len())?;
    check_len("diffusion tensor", n, a.c2.len())?;
    check_len("load", n, a.load.len())?;
    check_len("convection velocity", n, a.velocity.len())?;
    if !a.mesh_velocity.is_empty() {
        check_len("mesh velocity", n, a.mesh_velocity.len())?;
    }
    if !a.soret_diffusivity.is_empty() {
        check_len("soret diffusivity", n, a.soret_diffusivity.len())?;
        check_len("soret c0", n, a.c0.len())?;
        check_len("temperature", n, a.temperature.len())?;
    }
    if !a.ne_const.is_empty() {
        check_len("nernst-einstein constant", n, a.ne_const.len())?;
        check_len("electric potential", n, a.potential.len())?;
    }
    Ok(())
}
