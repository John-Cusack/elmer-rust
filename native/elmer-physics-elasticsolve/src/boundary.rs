//! Boundary element system of `ElasticSolve.F90` `LocalBoundaryMatrix` (`:5318`):
//! surface traction, normal traction (pressure) and spring conditions, with the
//! cofactor terms of `Large Deflection` and `Pseudo-Traction`.
//!
//! Not covered: FSI traction, lumping loads, the mixed formulation and the
//! `Normal-Tangential` rotation (the assembler's job).

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};
use elmer_physics_elastic::{CoordSystem, Spring};

use crate::error::ElasticError;
use crate::tensor::{self as t3, identity, T3, ZERO};

/// The parent element the boundary element's displacement gradient is read from
/// (`GetElementNodes( ParentNodes, ParentElement )`).
#[derive(Clone, Debug)]
pub struct ParentInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// Parent displacement at its nodes, `cdim` components used.
    pub disp: &'a [[f64; 3]],
    /// Local coordinates of each boundary node in the parent, in the boundary
    /// element's node order (`Parent % TYPE % NodeU/V/W` of the matching node).
    pub node_uvw: &'a [[f64; 3]],
    /// Centre of the parent: the normal points away from it
    /// (`CheckNormalDirectionParent`).
    pub center: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct BoundaryInput<'a> {
    /// Boundary element type (202, 303, 404, ...).
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub cdim: usize,
    pub coord: CoordSystem,
    pub large_deflection: bool,
    /// `Pseudo-Traction`: the traction is a force per reference area, so it is
    /// not scaled by the deformed-to-reference area ratio.
    pub pseudo_traction: bool,
    /// `Surface Traction i` / `Force i`, per node.
    pub traction: &'a [[f64; 3]],
    /// `Normal Surface Traction` / `Normal Force`, per node.
    pub normal_traction: &'a [f64],
    pub spring: Spring<'a>,
    pub parent: ParentInput<'a>,
    /// Gauss points per element; `None` is the Elmer default rule.
    pub np: Option<usize>,
    /// `DOFs` per basis function in the local system: `cdim`, or `cdim + 1`
    /// with the pressure of `Incompressible` (its rows get no load).
    pub stride: usize,
}

/// Unit normal at a point from the tangential gradients the basis returns:
/// `sum_k x_k (x) grad N_k` is the projector onto the tangent space.
fn normal_from_projector(xyz: &[[f64; 3]], dbasis_dx: &[[f64; 3]], cdim: usize) -> Result<[f64; 3], ElasticError> {
    let mut proj = [[0.0; 3]; 3];
    for (x, d) in xyz.iter().zip(dbasis_dx) {
        for a in 0..3 {
            for b in 0..3 {
                proj[a][b] += x[a] * d[b];
            }
        }
    }
    let mut nn = [[0.0; 3]; 3];
    for a in 0..cdim {
        for b in 0..cdim {
            nn[a][b] = if a == b { 1.0 } else { 0.0 } - proj[a][b];
        }
    }
    // For an element of dimension cdim - 1, nn is the rank-one n n^T.
    let c = (0..cdim).max_by(|&a, &b| nn[a][a].partial_cmp(&nn[b][b]).unwrap()).unwrap();
    let d = nn[c][c];
    if d <= 0.0 {
        return Err(ElasticError::Basis("degenerate boundary element".into()));
    }
    let inv = 1.0 / d.sqrt();
    let mut n = [0.0; 3];
    for a in 0..cdim {
        n[a] = nn[a][c] * inv;
    }
    Ok(n)
}

/// The unit normal of a boundary element at a point, unoriented: on a
/// linear line in 2-D `NormalVector`'s own arithmetic (ElementGeometry.F90:
/// `(-dydu, dxdu) / SQRT(dxdu**2 + dydu**2)` with `dN/du = (-1/2, 1/2)`),
/// otherwise [`normal_from_projector`].
fn boundary_normal(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    dbasis_dx: &[[f64; 3]],
    cdim: usize,
) -> Result<[f64; 3], ElasticError> {
    if cdim == 2 && etype.0 == 202 {
        let du = [-0.5, 0.5];
        let (mut dxdu, mut dydu) = (0.0, 0.0);
        for k in 0..2 {
            dxdu += du[k] * xyz[k][0];
        }
        for k in 0..2 {
            dydu += du[k] * xyz[k][1];
        }
        let det = dxdu * dxdu + dydu * dydu;
        if det <= 0.0 {
            return Err(ElasticError::Basis("degenerate boundary element".into()));
        }
        let det = 1.0 / det.sqrt();
        return Ok([-dydu * det, dxdu * det, 0.0]);
    }
    normal_from_projector(xyz, dbasis_dx, cdim)
}

/// Boundary element force and spring matrix. Overwrites `out`.
pub fn elastic_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &BoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ElasticError> {
    let points = quad.gauss_points(inp.etype, inp.np);
    elastic_boundary_local_p(basis, &points, inp.etype.n_nodes(), inp, out)
}

/// [`elastic_boundary_local`] over the `ntot >= n` basis functions of a p
/// boundary element on the given rule (`LocalBoundaryMatrix` with `ntot`):
/// the loads are interpolated with the nodal functions and act on every
/// basis function. The parent's deformation (`Large Deflection`) is not
/// available on p elements.
pub fn elastic_boundary_local_p(
    basis: &dyn ElementBasis,
    points: &[elmer_core::QuadPoint],
    ntot: usize,
    inp: &BoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ElasticError> {
    let n = inp.etype.n_nodes();
    if ntot > n && inp.large_deflection {
        return Err(ElasticError::Unsupported("Large Deflection boundary loads on p elements".into()));
    }
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    let edim = match inp.etype.family() {
        2 => 1,
        3 | 4 => 2,
        f => return Err(ElasticError::Unsupported(format!("boundary element family {f}"))),
    };
    if edim + 1 != cdim {
        return Err(ElasticError::Unsupported(format!("{edim}D boundary element in a {cdim}D mesh")));
    }
    if inp.xyz.len() != n || inp.traction.len() != n || inp.normal_traction.len() != n || inp.parent.node_uvw.len() != n {
        return Err(ElasticError::Input(format!("element {} needs {n} nodal values", inp.etype.0)));
    }
    let dofs = inp.stride.max(cdim);
    let ndof = ntot * dofs;
    *out = LocalSystem::zeros(ndof);
    let eye = identity(cdim);

    // Orientation is a property of the element: decided at its centre.
    let centre_uvw = if inp.etype.family() == 3 { [1.0 / 3.0, 1.0 / 3.0, 0.0] } else { [0.0; 3] };
    let mut be = BasisEval::default();
    basis.eval(inp.etype, inp.xyz, centre_uvw, &mut be).map_err(ElasticError::Basis)?;
    let centre: [f64; 3] = std::array::from_fn(|a| (0..n).map(|k| be.basis[k] * inp.xyz[k][a]).sum());
    let n_centre = boundary_normal(inp.etype, inp.xyz, &be.dbasis_dx, cdim)?;
    let flip = (0..3).map(|a| (inp.parent.center[a] - centre[a]) * n_centre[a]).sum::<f64>() > 0.0;
    let sign = if flip { -1.0 } else { 1.0 };

    let mut pbe = BasisEval::default();
    for qp in points {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(ElasticError::Basis)?;
        let bas = &be.basis;
        let mut s = be.det_j * qp.weight;
        let mut r = 0.0;
        if axi {
            r = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            s *= r;
        }

        let mut ref_normal = boundary_normal(inp.etype, inp.xyz, &be.dbasis_dx, cdim)?;
        for v in &mut ref_normal {
            *v *= sign;
        }

        // Cofactor of the deformation gradient at the point, from the parent.
        let cof: T3 = if inp.large_deflection {
            let puvw: [f64; 3] =
                std::array::from_fn(|a| (0..n).map(|k| bas[k] * inp.parent.node_uvw[k][a]).sum());
            basis.eval(inp.parent.etype, inp.parent.xyz, puvw, &mut pbe).map_err(ElasticError::Basis)?;
            let np = inp.parent.xyz.len();
            let mut grad = ZERO;
            for a in 0..cdim {
                for b in 0..cdim {
                    grad[a][b] = (0..np).map(|k| inp.parent.disp[k][a] * pbe.dbasis_dx[k][b]).sum();
                }
            }
            let defg = t3::add(&eye, &grad);
            let mut det = t3::det(&defg, cdim);
            if axi {
                det *= 1.0 + (0..np).map(|k| pbe.basis[k] * inp.parent.disp[k][0]).sum::<f64>() / r;
            }
            let inv = t3::invert(&defg, cdim)
                .ok_or_else(|| ElasticError::Numerical("singular deformation gradient at the boundary".into()))?;
            t3::scale(&t3::transpose(&inv), det)
        } else {
            eye
        };
        let normal: [f64; 3] =
            std::array::from_fn(|a| if a < cdim { (0..cdim).map(|b| cof[a][b] * ref_normal[b]).sum() } else { 0.0 });
        let metric = (0..cdim).map(|a| normal[a] * normal[a]).sum::<f64>().sqrt();

        let beta: f64 = (0..n).map(|k| inp.normal_traction[k] * bas[k]).sum();
        let pressure_force: [f64; 3] = std::array::from_fn(|a| beta * normal[a]);
        let traction: [f64; 3] =
            std::array::from_fn(|a| if a < cdim { (0..n).map(|k| inp.traction[k][a] * bas[k]).sum() } else { 0.0 });
        let scale = if inp.pseudo_traction { 1.0 } else { metric };
        // LocalBoundaryMatrix adds the normal traction, then the traction:
        // `s * Basis(q) * Force(i)` and `s * Basis(q) * Force(i) * MetricTerm`.
        for q in 0..ntot {
            for i in 0..cdim {
                out.force[q * dofs + i] += s * bas[q] * pressure_force[i];
            }
        }
        for q in 0..ntot {
            for i in 0..cdim {
                out.force[q * dofs + i] += if inp.pseudo_traction {
                    s * bas[q] * traction[i]
                } else {
                    s * bas[q] * traction[i] * scale
                };
            }
        }

        match &inp.spring {
            Spring::None => {}
            Spring::Normal(k) => {
                let kc: f64 = (0..n).map(|m| bas[m] * k[m]).sum();
                for p in 0..ntot {
                    for i in 0..cdim {
                        for q in 0..ntot {
                            for j in 0..cdim {
                                out.stiff[(p * dofs + i) * ndof + q * dofs + j] +=
                                    kc * bas[q] * ref_normal[j] * bas[p] * ref_normal[i] * s;
                            }
                        }
                    }
                }
            }
            Spring::Tensor(k) => {
                let mut kc = [[0.0; 3]; 3];
                for (i, row) in kc.iter_mut().enumerate().take(cdim) {
                    for (j, v) in row.iter_mut().enumerate().take(cdim) {
                        *v = (0..n).map(|m| bas[m] * k[m][i][j]).sum();
                    }
                }
                for p in 0..ntot {
                    for i in 0..cdim {
                        for q in 0..ntot {
                            for j in 0..cdim {
                                out.stiff[(p * dofs + i) * ndof + q * dofs + j] += kc[i][j] * bas[q] * bas[p] * s;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

pub use elmer_physics_elastic::nodes::node_uvw;
