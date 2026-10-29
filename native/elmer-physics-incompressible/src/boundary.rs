//! `LocalBoundaryMatrix` (IncompressibleNS.F90:1749) for the boundary terms
//! this crate ports: `Surface Traction i`, `Normal Surface Traction` /
//! `External Pressure` along the outward normal, and `Slip Coefficient i`.
//!
//! The element integrals are those of `NavierStokesBoundary`, which
//! `elmer-physics-flow` already carries bit for bit, so this delegates to
//! `flow_boundary_local` with no mass flux.

use elmer_core::{ElementBasis, ElementType, LocalSystem, Quadrature};
use elmer_physics_flow::{flow_boundary_local, FlowBoundaryInput, FlowError};

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
