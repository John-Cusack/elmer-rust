//! Boundary element matrix and force of the linear small-strain path of
//! `ElasticSolve.F90` `LocalBoundaryMatrix` (`:5318`): surface traction,
//! normal traction (pressure) and spring conditions.
//!
//! Not covered: FSI traction, large-deflection cofactor terms, lumping loads
//! and `Normal-Tangential` rotation (the assembler's job).

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint, Quadrature};

use crate::bulk::CoordSystem;
use crate::error::StressError;

/// Spring condition of a boundary, per element node.
#[derive(Clone, Debug)]
pub enum Spring<'a> {
    None,
    /// `Spring`: `k (u.n) n` against the reference normal.
    Normal(&'a [f64]),
    /// `Spring i` / `Spring ij`: `K(i,j)` per node.
    Tensor(&'a [[[f64; 3]; 3]]),
}

#[derive(Clone, Debug)]
pub struct StressBoundaryInput<'a> {
    /// Boundary element type (202, 303, 404, ...).
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub cdim: usize,
    pub coord: CoordSystem,
    /// `Surface Traction i` / `Force i`, per node.
    pub traction: &'a [[f64; 3]],
    /// `Normal Surface Traction` / `Normal Force`, per node (positive along the
    /// outward normal).
    pub normal_traction: &'a [f64],
    pub spring: Spring<'a>,
    /// Centre of the parent element Elmer chose for the boundary (left parent,
    /// else right), as `NodalBasisFunctions` at `u = 1/3 | 0` gives it. The
    /// normal is turned to point away from it (`CheckNormalDirectionParent`).
    pub parent_center: [f64; 3],
    /// Gauss points per element; `None` is the Elmer default rule.
    pub np: Option<usize>,
}

/// Unit normal at a point from the tangential gradients the basis returns:
/// `sum_k x_k (x) grad N_k` is the projector onto the tangent space.
fn normal_from_projector(
    xyz: &[[f64; 3]],
    dbasis_dx: &[[f64; 3]],
    cdim: usize,
    edim: usize,
) -> Result<[f64; 3], StressError> {
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
    // edim == cdim - 1: nn is the rank-one n n^T.
    debug_assert_eq!(edim + 1, cdim);
    let c = (0..cdim)
        .max_by(|&a, &b| nn[a][a].partial_cmp(&nn[b][b]).unwrap())
        .unwrap();
    let d = nn[c][c];
    if d <= 0.0 {
        return Err(StressError::Basis("degenerate boundary element".into()));
    }
    let inv = 1.0 / d.sqrt();
    let mut n = [0.0; 3];
    for a in 0..cdim {
        n[a] = nn[a][c] * inv;
    }
    Ok(n)
}

/// Boundary element force and spring matrix. Overwrites `out`.
pub fn stress_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &StressBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    let points = quad.gauss_points(inp.etype, inp.np);
    boundary_impl(basis, &points, inp.etype.n_nodes(), inp, out)
}

/// [`stress_boundary_local`] for `nd >= n` basis functions per component on
/// the given rule (`StressBoundary(..., n, ntot, ...)`): loads, springs and
/// the normal come from the nodal functions; `out` is `nd * cdim` square.
pub fn stress_boundary_local_p(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nd: usize,
    inp: &StressBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    if nd < inp.etype.n_nodes() {
        return Err(StressError::Input(format!("{nd} basis functions for {} nodes", inp.etype.n_nodes())));
    }
    boundary_impl(basis, points, nd, inp, out)
}

fn boundary_impl(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nd: usize,
    inp: &StressBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    let n = inp.etype.n_nodes();
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    let edim = match inp.etype.family() {
        2 => 1,
        3 | 4 => 2,
        f => {
            return Err(StressError::Unsupported(format!(
                "boundary element family {f}"
            )))
        }
    };
    if edim + 1 != cdim {
        return Err(StressError::Unsupported(format!(
            "{}D boundary element in a {cdim}D mesh",
            edim
        )));
    }
    if inp.xyz.len() != n || inp.traction.len() != n || inp.normal_traction.len() != n {
        return Err(StressError::Input(format!(
            "element {} needs {n} nodal values",
            inp.etype.0
        )));
    }
    let ndof = nd * cdim;
    *out = LocalSystem::zeros(ndof);

    // Orientation is a property of the element: decided at its centre.
    let centre_uvw = if inp.etype.family() == 3 {
        [1.0 / 3.0, 1.0 / 3.0, 0.0]
    } else {
        [0.0; 3]
    };
    let mut be = BasisEval::default();
    basis
        .eval(inp.etype, inp.xyz, centre_uvw, &mut be)
        .map_err(StressError::Basis)?;
    let centre: [f64; 3] =
        std::array::from_fn(|a| (0..n).map(|k| be.basis[k] * inp.xyz[k][a]).sum());
    let n_centre = normal_from_projector(inp.xyz, &be.dbasis_dx, cdim, edim)?;
    let flip = (0..3)
        .map(|a| (inp.parent_center[a] - centre[a]) * n_centre[a])
        .sum::<f64>()
        > 0.0;
    let sign = if flip { -1.0 } else { 1.0 };

    for qp in points {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(StressError::Basis)?;
        if be.basis.len() < nd {
            return Err(StressError::Input(format!("{} basis functions for {nd} dofs", be.basis.len())));
        }
        let bas = &be.basis;
        let mut s = be.det_j * qp.weight;
        if axi {
            let r: f64 = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            s *= r;
        }
        let mut normal = normal_from_projector(inp.xyz, &be.dbasis_dx, cdim, edim)?;
        for v in &mut normal {
            *v *= sign;
        }

        let beta: f64 = (0..n).map(|k| inp.normal_traction[k] * bas[k]).sum();
        let mut force = [0.0; 3];
        for i in 0..cdim {
            force[i] = beta * normal[i] + (0..n).map(|k| inp.traction[k][i] * bas[k]).sum::<f64>();
        }
        for q in 0..nd {
            for i in 0..cdim {
                out.force[q * cdim + i] += s * bas[q] * force[i];
            }
        }

        match &inp.spring {
            Spring::None => {}
            Spring::Normal(k) => {
                let kc: f64 = (0..n).map(|m| bas[m] * k[m]).sum();
                for p in 0..nd {
                    for i in 0..cdim {
                        for q in 0..nd {
                            for j in 0..cdim {
                                out.stiff[(p * cdim + i) * ndof + q * cdim + j] +=
                                    kc * bas[q] * normal[j] * bas[p] * normal[i] * s;
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
                for p in 0..nd {
                    for i in 0..cdim {
                        for q in 0..nd {
                            for j in 0..cdim {
                                out.stiff[(p * cdim + i) * ndof + q * cdim + j] +=
                                    kc[i][j] * bas[q] * bas[p] * s;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
