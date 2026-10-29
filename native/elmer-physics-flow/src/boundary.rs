//! `NavierStokesBoundary` (NavierStokes.F90:1320) in cartesian coordinates:
//! the traction boundary terms FlowSolve assembles on every boundary element
//! with `Flow Force BC` (default on).
//!
//! Ported: a force in the coordinate directions (`Pressure i`), a normal
//! force (`External Pressure` / `Normal Pressure`, along the outward normal),
//! a `Mass Flux` into the continuity equation and the slip resistance
//! (`Slip Coefficient i`) in coordinate directions. Normal-tangential
//! coordinates, surface tension and its tangential derivative are not
//! ported; the solver refuses them.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::FlowError;

/// What `NavierStokesBoundary` reads for one boundary element, at its nodes.
pub struct FlowBoundaryInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// Coordinate system dimension of the bulk.
    pub dim: usize,
    /// `Pressure 1..3`: traction in the coordinate directions.
    pub traction: &'a [[f64; 3]],
    /// `External Pressure`: traction along the outward normal.
    pub ext_pressure: &'a [f64],
    /// `Mass Flux`.
    pub mass_flux: &'a [f64],
    /// `Slip Coefficient 1..3`.
    pub slip: &'a [[f64; 3]],
    /// A point inside the parent element (its centre), used to orient the
    /// normal outward (`CheckNormalDirection`).
    pub parent_center: Option<[f64; 3]>,
    pub np: Option<usize>,
}

/// The mean of the element nodes: the reference-element centre of a linear
/// element (`NodalBasisFunctions` at the centre in `CheckNormalDirection`).
pub fn element_center(xyz: &[[f64; 3]]) -> [f64; 3] {
    let n = xyz.len() as f64;
    std::array::from_fn(|a| xyz.iter().map(|x| x[a]).sum::<f64>() / n)
}

/// Unit normal of a line (in 2D) or face (in 3D) at reference point `uvw`
/// (`NormalVector` without the orientation check), from the parametric
/// derivatives of the geometry.
fn raw_normal(etype: ElementType, xyz: &[[f64; 3]], uvw: [f64; 3]) -> Result<[f64; 3], FlowError> {
    let (u, v) = (uvw[0], uvw[1]);
    // dx/du, dx/dv of the linear boundary element maps.
    let (du, dv): ([f64; 3], [f64; 3]) = match etype.0 {
        202 => (std::array::from_fn(|a| (xyz[1][a] - xyz[0][a]) / 2.0), [0.0; 3]),
        303 => (
            std::array::from_fn(|a| xyz[1][a] - xyz[0][a]),
            std::array::from_fn(|a| xyz[2][a] - xyz[0][a]),
        ),
        404 => (
            std::array::from_fn(|a| {
                0.25 * ((1.0 - v) * (xyz[1][a] - xyz[0][a]) + (1.0 + v) * (xyz[2][a] - xyz[3][a]))
            }),
            std::array::from_fn(|a| {
                0.25 * ((1.0 - u) * (xyz[3][a] - xyz[0][a]) + (1.0 + u) * (xyz[2][a] - xyz[1][a]))
            }),
        ),
        t => return Err(FlowError::Unsupported(format!("boundary normal of element type {t}"))),
    };
    let nrm = if etype.family() == 2 {
        [-du[1], du[0], 0.0]
    } else {
        [
            du[1] * dv[2] - dv[1] * du[2],
            dv[0] * du[2] - du[0] * dv[2],
            du[0] * dv[1] - dv[0] * du[1],
        ]
    };
    let len = nrm.iter().map(|x| x * x).sum::<f64>().sqrt();
    if len <= 0.0 {
        return Ok([0.0; 3]);
    }
    Ok(nrm.map(|x| x / len))
}

/// `NormalVector(Element, Nodes, u, v, .TRUE.)`: the unit normal pointing
/// away from `parent_center`.
pub fn outward_normal(
    etype: ElementType,
    xyz: &[[f64; 3]],
    uvw: [f64; 3],
    parent_center: [f64; 3],
) -> Result<[f64; 3], FlowError> {
    let nrm = raw_normal(etype, xyz, uvw)?;
    let c = element_center(xyz);
    let d: f64 = (0..3).map(|a| (parent_center[a] - c[a]) * nrm[a]).sum();
    Ok(if d > 0.0 { nrm.map(|x| -x) } else { nrm })
}

/// The local boundary matrix and force vector of one boundary element.
/// `out` is resized to `n * (dim + 1)`; there is no mass matrix.
pub fn flow_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &FlowBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), FlowError> {
    let n = inp.xyz.len();
    let dim = inp.dim;
    let c = dim + 1;
    for (what, got) in [
        ("element nodes", inp.etype.n_nodes()),
        ("traction", inp.traction.len()),
        ("external pressure", inp.ext_pressure.len()),
        ("mass flux", inp.mass_flux.len()),
        ("slip coefficient", inp.slip.len()),
    ] {
        if got != n {
            return Err(FlowError::NodalLength { what, expected: n, got });
        }
    }
    let nd = n * c;
    *out = LocalSystem::zeros(nd);

    let has_normal_force = inp.ext_pressure.iter().any(|&p| p != 0.0);
    let has_slip = inp.slip.iter().any(|s| s.iter().any(|&x| x != 0.0));
    let parent_center = match (has_normal_force, inp.parent_center) {
        (true, None) => {
            return Err(FlowError::Unsupported(
                "a normal force on a boundary element with no parent (normal orientation)".into(),
            ))
        }
        (_, pc) => pc,
    };

    let mut be = BasisEval::default();
    let points = quad.gauss_points(inp.etype, inp.np);
    if points.is_empty() {
        return Err(FlowError::NoQuadrature(inp.etype.0));
    }
    for qp in &points {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(FlowError::Basis)?;
        let s = be.det_j * qp.weight;
        let phi = &be.basis[..n];

        let mut force = [0.0; 3];
        for (i, f) in force.iter_mut().enumerate().take(dim) {
            *f = (0..n).map(|k| inp.traction[k][i] * phi[k]).sum();
        }
        if let Some(pc) = parent_center.filter(|_| has_normal_force) {
            let normal = outward_normal(inp.etype, inp.xyz, qp.uvw, pc)?;
            let alpha: f64 = (0..n).map(|k| inp.ext_pressure[k] * phi[k]).sum();
            for i in 0..dim {
                force[i] += alpha * normal[i];
            }
        }
        let mass_flux: f64 = (0..n).map(|k| inp.mass_flux[k] * phi[k]).sum();

        if has_slip {
            for i in 0..dim {
                let slip: f64 = (0..n).map(|k| inp.slip[k][i] * phi[k]).sum();
                for p in 0..n {
                    for q in 0..n {
                        out.stiff[(p * c + i) * nd + q * c + i] += s * slip * phi[q] * phi[p];
                    }
                }
            }
        }
        for q in 0..n {
            for (i, f) in force.iter().enumerate().take(dim) {
                out.force[q * c + i] += s * phi[q] * f;
            }
            out.force[q * c + dim] += s * mass_flux * phi[q];
        }
    }
    Ok(())
}
