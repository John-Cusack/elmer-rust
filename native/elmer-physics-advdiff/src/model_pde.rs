//! `LocalMatrix` and `LocalMatrixBC` of `AdvDiffSolver` in ModelPDE.F90: the
//! prototype advection-diffusion-reaction equation
//! `rho du/dt + div(C a u) - div(D grad u) + R u = f` of the Elmer tutorials,
//! for one scalar degree of freedom per node.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint};

use crate::error::{check_len, interp, AdvDiffError};
use crate::fe::Fe;

/// Nodal material and body force values of one bulk element, from the
/// Material keywords `Diffusion Coefficient`, `Convection Coefficient`,
/// `Reaction Coefficient`, `Time Derivative Coefficient` and
/// `Convection Velocity <i>`, and the body force `Field Source`. A keyword the
/// section lacks reads as zero, as `GetReal` does.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ModelPdeNodal {
    pub diffusion: Vec<f64>,
    pub convection: Vec<f64>,
    pub reaction: Vec<f64>,
    pub time: Vec<f64>,
    pub velocity: Vec<[f64; 3]>,
    pub source: Vec<f64>,
}

/// Boundary condition values at the nodes of a boundary element: `Field Flux`,
/// `Robin Coefficient` and `External Field`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ModelPdeBcNodal {
    pub flux: Vec<f64>,
    pub robin: Vec<f64>,
    pub external: Vec<f64>,
}

/// Bulk element matrices of `LocalMatrix` for an element with `nd = n` nodes
/// (no p-bubble degrees of freedom).
///
/// `out` is reset to `n x n`: `stiff` holds diffusion, convection and
/// reaction, `mass` the time derivative term, `force` the source. Matrices are
/// row-major with the test function on the row. `dim` is
/// `CoordinateSystemDimension()`; `quad_points` the requested Gauss point
/// count (`None` for the element default).
pub fn model_pde_bulk_local(
    fe: &Fe<'_>,
    etype: ElementType,
    xyz: &[[f64; 3]],
    nodal: &ModelPdeNodal,
    dim: usize,
    quad_points: Option<usize>,
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    let ips = fe.quad.gauss_points(etype, quad_points);
    model_pde_bulk_local_with(fe.basis, &ips, etype, xyz, nodal, dim, n, out)
}

/// [`model_pde_bulk_local`] for `nd >= n` basis functions (the nodal ones
/// first, then edge, face and bubble functions of a p element) on the
/// given rule: `LocalMatrix(Element, n, nd)`, coefficients interpolated
/// with the nodal functions only. `out` is reset to `nd x nd`; bubbles are
/// not condensed here.
#[allow(clippy::too_many_arguments)]
pub fn model_pde_bulk_local_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    nodal: &ModelPdeNodal,
    dim: usize,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    check_len("diffusion coefficient", n, nodal.diffusion.len())?;
    check_len("convection coefficient", n, nodal.convection.len())?;
    check_len("reaction coefficient", n, nodal.reaction.len())?;
    check_len("time derivative coefficient", n, nodal.time.len())?;
    check_len("convection velocity", n, nodal.velocity.len())?;
    check_len("field source", n, nodal.source.len())?;
    if !(1..=3).contains(&dim) {
        return Err(AdvDiffError::Invalid(format!("dim {dim} not in 1..=3")));
    }
    if nd < n {
        return Err(AdvDiffError::Invalid(format!("{nd} dofs for {n} nodes")));
    }

    out.n = nd;
    out.stiff.clear();
    out.stiff.resize(nd * nd, 0.0);
    out.force.clear();
    out.force.resize(nd, 0.0);
    out.damp = None;
    let mut mass = out.mass.take().unwrap_or_default();
    mass.clear();
    mass.resize(nd * nd, 0.0);

    if ips.is_empty() {
        return Err(AdvDiffError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(AdvDiffError::Basis)?;
        if be.basis.len() < nd {
            return Err(AdvDiffError::Invalid(format!(
                "{} basis functions for {nd} dofs",
                be.basis.len()
            )));
        }
        let phi = &be.basis[..nd];
        let dphi = &be.dbasis_dx[..nd];

        let nphi = &phi[..n];
        let load = interp(nphi, &nodal.source);
        let rho = interp(nphi, &nodal.time);
        let mut a = [0.0; 3];
        for (i, ai) in a.iter_mut().enumerate().take(dim) {
            *ai = nphi
                .iter()
                .zip(&nodal.velocity)
                .map(|(b, v)| b * v[i])
                .sum();
        }
        let d = interp(nphi, &nodal.diffusion);
        let c = interp(nphi, &nodal.convection);
        let r = interp(nphi, &nodal.reaction);
        let w = qp.weight * be.det_j;

        for p in 0..nd {
            for q in 0..nd {
                let grad: f64 = (0..3).map(|i| dphi[p][i] * dphi[q][i]).sum();
                let adv: f64 = (0..dim).map(|i| a[i] * dphi[q][i]).sum();
                let st = &mut out.stiff[p * nd + q];
                *st += w * d * grad;
                *st += w * c * adv * phi[p];
                *st += w * r * phi[q] * phi[p];
                mass[p * nd + q] += w * rho * phi[q] * phi[p];
            }
            out.force[p] += w * load * phi[p];
        }
    }
    out.mass = Some(mass);
    Ok(())
}

/// Boundary element matrices of `LocalMatrixBC`: `stiff = int c N_p N_q` and
/// `force = int (flux + c u_ext) N_p`, with `c` the `Robin Coefficient`.
///
/// `out` is reset to `n x n` with no mass.
pub fn model_pde_boundary_local(
    fe: &Fe<'_>,
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &ModelPdeBcNodal,
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let ips = fe.quad.gauss_points(etype, None);
    model_pde_boundary_local_with(fe.basis, &ips, etype, xyz, bc, xyz.len(), out)
}

/// [`model_pde_boundary_local`] for `nd >= n` basis functions on the given
/// rule (`LocalMatrixBC(Element, n, nd)`); `out` is reset to `nd x nd`.
pub fn model_pde_boundary_local_with(
    basis: &dyn ElementBasis,
    ips: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    bc: &ModelPdeBcNodal,
    nd: usize,
    out: &mut LocalSystem,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    check_len("field flux", n, bc.flux.len())?;
    check_len("robin coefficient", n, bc.robin.len())?;
    check_len("external field", n, bc.external.len())?;
    if nd < n {
        return Err(AdvDiffError::Invalid(format!("{nd} dofs for {n} nodes")));
    }
    out.n = nd;
    out.stiff.clear();
    out.stiff.resize(nd * nd, 0.0);
    out.force.clear();
    out.force.resize(nd, 0.0);
    out.mass = None;
    out.damp = None;

    if ips.is_empty() {
        return Err(AdvDiffError::NoQuadrature(etype.0));
    }
    let mut be = BasisEval::default();
    for qp in ips {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(AdvDiffError::Basis)?;
        if be.basis.len() < nd {
            return Err(AdvDiffError::Invalid(format!(
                "{} basis functions for {nd} dofs",
                be.basis.len()
            )));
        }
        let phi = &be.basis[..nd];
        let nphi = &phi[..n];
        let w = qp.weight * be.det_j;
        let f = interp(nphi, &bc.flux);
        let c = interp(nphi, &bc.robin);
        let ext = interp(nphi, &bc.external);
        for p in 0..nd {
            for q in 0..nd {
                out.stiff[p * nd + q] += w * c * phi[q] * phi[p];
            }
            out.force[p] += w * (f + c * ext) * phi[p];
        }
    }
    Ok(())
}
