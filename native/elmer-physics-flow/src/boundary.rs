//! `NavierStokesBoundary` (NavierStokes.F90:1320) in cartesian coordinates:
//! the traction boundary terms FlowSolve assembles on every boundary element
//! with `Flow Force BC` (default on).
//!
//! Ported: a force in the coordinate directions (`Pressure i`), a normal
//! force (`External Pressure` / `Normal Pressure`, along the outward normal),
//! a `Mass Flux` into the continuity equation and the slip resistance
//! (`Slip Coefficient i`), in coordinate or normal-tangential directions,
//! and surface tension: `Surface Tension Coefficient` (`Alpha`, the term
//! `-Alpha dBasisdx`) and the tangential force of its variation along the
//! boundary (`Beta`, the Marangoni force `grad(Beta) Basis`).

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
    /// `Normal-Tangential Velocity`: the forces and slip coefficients are
    /// given in normal-tangential components (the normal force adds to the
    /// first), assembled through the element normal and tangents.
    pub normal_tangential: bool,
    /// Surface tension, `(Alpha, Beta)` at the nodes: `Surface Tension
    /// Coefficient`, and the same scaled by `1 - Surface Tension Expansion
    /// Coefficient * Temperature` (FlowSolve.F90:1279-1302).
    pub surface_tension: Option<(&'a [f64], &'a [f64])>,
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

/// `TangentDirections(Normal, Tangent1, Tangent2)`.
fn tangent_directions(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (n1, n2, n3) = (n[0].abs(), n[1].abs(), n[2].abs());
    let mut t1 = if n1 <= n3 && n2 <= n3 { [0.0, -n[2], n[1]] } else { [-n[1], n[0], 0.0] };
    let s = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    t1 = t1.map(|x| x / s);
    let mut t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];
    let s = (t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]).sqrt();
    t2 = t2.map(|x| x / s);
    (t1, t2)
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
    let nt = inp.normal_tangential;
    let parent_center = match (has_normal_force || nt, inp.parent_center) {
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
        let dphi = &be.dbasis_dx[..n];
        // TangentForce(i) = SUM(NodalBeta * dBasisdx(:,i)), Alpha at the point.
        let (alpha, tangent_force) = match inp.surface_tension {
            Some((a, b)) => (
                (0..n).fold(0.0, |acc, k| acc + a[k] * phi[k]),
                std::array::from_fn::<f64, 3, _>(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * dphi[k][i])),
            ),
            None => (0.0, [0.0; 3]),
        };

        let mut force = [0.0; 3];
        for (i, f) in force.iter_mut().enumerate().take(dim) {
            *f = (0..n).map(|k| inp.traction[k][i] * phi[k]).sum();
        }
        // NavierStokesBoundary: the normal, and in n-t form the tangents.
        let mut vects = [[0.0; 3]; 3];
        if let Some(pc) = parent_center.filter(|_| has_normal_force || nt) {
            let normal = outward_normal(inp.etype, inp.xyz, qp.uvw, pc)?;
            let alpha: f64 = (0..n).map(|k| inp.ext_pressure[k] * phi[k]).sum();
            if nt {
                force[0] += alpha;
            } else {
                for i in 0..dim {
                    force[i] += alpha * normal[i];
                }
            }
            vects[0] = normal;
            if inp.etype.family() == 2 {
                vects[1] = [normal[1], -normal[0], 0.0];
            } else {
                let (t1, t2) = tangent_directions(normal);
                vects[1] = t1;
                vects[2] = t2;
            }
        }
        let mass_flux: f64 = (0..n).map(|k| inp.mass_flux[k] * phi[k]).sum();

        if has_slip && nt {
            for p in 0..n {
                for q in 0..n {
                    for (i, vect) in vects.iter().enumerate().take(dim) {
                        let slip: f64 = (0..n).map(|k| inp.slip[k][i] * phi[k]).sum();
                        for j in 0..dim {
                            for k in 0..dim {
                                out.stiff[(p * c + j) * nd + q * c + k] += s * slip * phi[q] * phi[p] * vect[j] * vect[k];
                            }
                        }
                    }
                }
            }
        } else if has_slip {
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
                if nt {
                    for j in 0..dim {
                        out.force[q * c + j] += s * phi[q] * f * vects[i][j];
                    }
                } else {
                    out.force[q * c + i] += s * phi[q] * f;
                }
                if inp.surface_tension.is_some() {
                    out.force[q * c + i] -= s * alpha * dphi[q][i];
                    out.force[q * c + i] += s * tangent_force[i] * phi[q];
                }
            }
            out.force[q * c + dim] += s * mass_flux * phi[q];
        }
    }
    Ok(())
}
