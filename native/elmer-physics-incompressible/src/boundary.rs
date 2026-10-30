//! `LocalBoundaryMatrix` (IncompressibleNS.F90:1749) for the boundary terms
//! this crate ports: `Surface Traction i`, `Normal Surface Traction` /
//! `External Pressure` along the outward normal, and `Slip Coefficient i`.
//!
//! The element integrals are those of `NavierStokesBoundary`, which
//! `elmer-physics-flow` already carries bit for bit, so this delegates to
//! `flow_boundary_local` with no mass flux; on p elements
//! [`incompressible_boundary_local_with`] integrates them over all the
//! boundary element's functions.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint, Quadrature};
use elmer_physics_flow::{flow_boundary_local, outward_normal, FlowBoundaryInput, FlowError};

use crate::error::IncompressibleError;

/// What `LocalBoundaryMatrix` reads for one boundary element, at its nodes.
pub struct IncompressibleBoundaryInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub dim: usize,
    /// `Surface Traction 1..3`.
    pub traction: &'a [[f64; 3]],
    /// `Normal Surface Traction`, else `External Pressure`.
    pub ext_pressure: &'a [f64],
    /// `Slip Coefficient 1..3`.
    pub slip: &'a [[f64; 3]],
    /// A point inside the parent element, orienting the normal outward.
    pub parent_center: Option<[f64; 3]>,
    pub np: Option<usize>,
    /// `Coordinate System = Axi Symmetric`: the measure is r ds
    /// ([`incompressible_boundary_local_with`] only).
    pub axisymmetric: bool,
}

fn convert(e: FlowError) -> IncompressibleError {
    match e {
        FlowError::NodalLength { what, expected, got } => IncompressibleError::NodalLength { what, expected, got },
        FlowError::Basis(s) => IncompressibleError::Basis(s),
        FlowError::NoQuadrature(t) => IncompressibleError::NoQuadrature(t),
        FlowError::Unsupported(s) => IncompressibleError::Unsupported(s),
    }
}

/// The local boundary matrix and force vector of one boundary element.
pub fn incompressible_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &IncompressibleBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), IncompressibleError> {
    if inp.axisymmetric {
        let points = quad.gauss_points(inp.etype, inp.np);
        return incompressible_boundary_local_with(basis, &points, inp.xyz.len(), inp, out);
    }
    let zero = vec![0.0; inp.xyz.len()];
    let f = FlowBoundaryInput {
        etype: inp.etype,
        xyz: inp.xyz,
        dim: inp.dim,
        traction: inp.traction,
        ext_pressure: inp.ext_pressure,
        mass_flux: &zero,
        slip: inp.slip,
        parent_center: inp.parent_center,
        np: inp.np,
    };
    flow_boundary_local(basis, quad, &f, out).map_err(convert)
}

/// [`incompressible_boundary_local`] for `nd >= n` basis functions (a p
/// boundary element: its nodes, then its edge and face functions) on the
/// given rule: `LocalBoundaryMatrix(Element, n, nd)` without normal-
/// tangential conditions. The loads and slip are interpolated with the
/// nodal functions. `out` is resized to `nd * (dim + 1)`.
pub fn incompressible_boundary_local_with(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nd: usize,
    inp: &IncompressibleBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), IncompressibleError> {
    let n = inp.xyz.len();
    let dim = inp.dim;
    let c = dim + 1;
    for (what, got) in [
        ("element nodes", inp.etype.n_nodes()),
        ("traction", inp.traction.len()),
        ("external pressure", inp.ext_pressure.len()),
        ("slip coefficient", inp.slip.len()),
    ] {
        if got != n {
            return Err(IncompressibleError::NodalLength { what, expected: n, got });
        }
    }
    if nd < n {
        return Err(IncompressibleError::Unsupported(format!("{nd} basis functions for {n} nodes")));
    }
    let m = nd * c;
    *out = LocalSystem::zeros(m);
    let has_normal_force = inp.ext_pressure.iter().any(|&p| p != 0.0);
    let has_slip = inp.slip.iter().any(|s| s.iter().any(|&x| x != 0.0));
    let parent_center = match (has_normal_force, inp.parent_center) {
        (true, None) => {
            return Err(IncompressibleError::Unsupported(
                "a normal force on a boundary element with no parent (normal orientation)".into(),
            ))
        }
        (_, pc) => pc,
    };
    if points.is_empty() {
        return Err(IncompressibleError::NoQuadrature(inp.etype.0));
    }
    let mut be = BasisEval::default();
    for qp in points {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(IncompressibleError::Basis)?;
        if be.basis.len() < nd {
            return Err(IncompressibleError::Unsupported(format!("{} basis functions for {nd} dofs", be.basis.len())));
        }
        let phi = &be.basis[..nd];
        let mut s = be.det_j * qp.weight;
        if inp.axisymmetric {
            s *= (0..n).map(|k| phi[k] * inp.xyz[k][0]).sum::<f64>();
        }
        let mut force = [0.0; 3];
        for (i, f) in force.iter_mut().enumerate().take(dim) {
            *f = (0..n).map(|k| inp.traction[k][i] * phi[k]).sum();
        }
        if let Some(pc) = parent_center.filter(|_| has_normal_force) {
            let normal = outward_normal(inp.etype, inp.xyz, qp.uvw, pc).map_err(convert)?;
            let alpha: f64 = (0..n).map(|k| inp.ext_pressure[k] * phi[k]).sum();
            for i in 0..dim {
                force[i] += alpha * normal[i];
            }
        }
        if has_slip {
            for i in 0..dim {
                let slip: f64 = (0..n).map(|k| inp.slip[k][i] * phi[k]).sum();
                for p in 0..nd {
                    for q in 0..nd {
                        out.stiff[(p * c + i) * m + q * c + i] += s * slip * phi[q] * phi[p];
                    }
                }
            }
        }
        for q in 0..nd {
            for (i, f) in force.iter().enumerate().take(dim) {
                out.force[q * c + i] += s * phi[q] * f;
            }
        }
    }
    Ok(())
}
