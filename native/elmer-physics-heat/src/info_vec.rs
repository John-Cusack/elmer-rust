//! `ElementInfoVec` (ElemInfo.F90:1686) for a nodal element, one point at a
//! time: the nodal basis and its local derivatives (`NodalBasisFunctions`,
//! `NodalFirstDerivatives`), then `ElementMetricVec` and
//! `ElementInfoVec_ElementBasisToGlobal` in their own arithmetic.
//!
//! The scalar `ElementInfo` (elmer-fe) divides the metric by its
//! determinant and can retry in quad precision. `ElementMetricVec`
//! multiplies by the reciprocal and stays in double precision, so the two
//! round differently. `LocalMatrixVec` integrates with the latter.
//!
//! This is the path `ElementInfoVec` takes for at most
//! `VECTOR_SMALL_THRESH` (9) points. With more points on an element whose
//! code is at least its node count (303, 404, 504, 706, 808, ...), the
//! Fortran evaluates the basis with the `H1Basis_*Nodal` routines instead,
//! which this module does not reproduce.

use elmer_core::{BasisEval, ElementBasis, ElementType};
use elmer_fe::basis::{nodal_basis, nodal_first_derivatives};
use elmer_fe::element_type;

/// [`ElementBasis`] with `ElementInfoVec`'s metric. `cdim` plays
/// `CoordinateSystemDimension()`.
#[derive(Copy, Clone, Debug)]
pub struct ElementInfoVec {
    pub cdim: usize,
}

impl ElementBasis for ElementInfoVec {
    fn eval(
        &self,
        etype: ElementType,
        nodal_xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let def = element_type(etype.0).map_err(|e| e.to_string())?;
        let n = def.n_nodes;
        if nodal_xyz.len() != n {
            return Err(format!(
                "element {}: {} nodal coordinates for {n} nodes",
                etype.0,
                nodal_xyz.len()
            ));
        }
        out.basis.clear();
        out.basis.resize(n, 0.0);
        out.dbasis_dx.clear();
        out.dbasis_dx.resize(n, [0.0; 3]);
        if def.code == 101 {
            out.det_j = 1.0;
            out.basis[0] = 1.0;
            return Ok(());
        }
        let [u, v, w] = uvw;
        let mut dl = vec![[0.0; 3]; n];
        nodal_basis(def, u, v, w, &mut out.basis);
        nodal_first_derivatives(def, u, v, w, &mut dl);
        let (det_j, ltog) = element_metric_vec(nodal_xyz, &dl, def.dim, self.cdim)
            .ok_or_else(|| format!("ElementMetricVec: degenerate {}D element", def.dim))?;
        out.det_j = det_j;
        // ElementInfoVec_ElementBasisToGlobal.
        for (d, l) in out.dbasis_dx.iter_mut().zip(&dl) {
            for j in 0..self.cdim {
                d[j] = match def.dim {
                    1 => l[0] * ltog[j][0],
                    2 => l[0] * ltog[j][0] + l[1] * ltog[j][1],
                    _ => l[0] * ltog[j][0] + l[1] * ltog[j][1] + l[2] * ltog[j][2],
                };
            }
        }
        Ok(())
    }
}

/// `ElementMetricVec` at one point: `sqrt(det G)` and `LtoGMap(i, 1:dim)`
/// for `i < cdim`, or `None` for a degenerate point (`DetJ <= TINY`).
fn element_metric_vec(
    xyz: &[[f64; 3]],
    dl: &[[f64; 3]],
    dim: usize,
    cdim: usize,
) -> Option<(f64, [[f64; 3]; 3])> {
    // dx(j, l) = sum_k dLBasisdx(k, l) * xyz(k, j)
    let mut dx = [[0.0f64; 3]; 3];
    for l in 0..dim {
        for j in 0..3 {
            let mut s = 0.0;
            for (d, x) in dl.iter().zip(xyz) {
                s += d[l] * x[j];
            }
            dx[j][l] = s;
        }
    }
    // The upper triangle of G = dx^T dx, packed as Elmer's utind:
    // (1,1) (1,2) (2,2) (1,3) (2,3) (3,3).
    let g = |i: usize, j: usize| -> f64 {
        match cdim {
            1 => dx[0][i] * dx[0][j],
            2 => dx[0][i] * dx[0][j] + dx[1][i] * dx[1][j],
            _ => dx[0][i] * dx[0][j] + dx[1][i] * dx[1][j] + dx[2][i] * dx[2][j],
        }
    };
    let mut m = [0.0f64; 6];
    let det = match dim {
        1 => {
            let d = g(0, 0);
            if d <= f64::MIN_POSITIVE {
                return None;
            }
            m[0] = 1.0 / d;
            d
        }
        2 => {
            let (g1, g2, g3) = (g(0, 0), g(0, 1), g(1, 1));
            let d = g1 * g3 - g2 * g2;
            if d <= f64::MIN_POSITIVE {
                return None;
            }
            let s = 1.0 / d;
            m[0] = s * g3;
            m[1] = -s * g2;
            m[2] = s * g1;
            d
        }
        _ => {
            let (g1, g2, g3, g4, g5, g6) = (g(0, 0), g(0, 1), g(1, 1), g(0, 2), g(1, 2), g(2, 2));
            let d = g1 * (g3 * g6 - g5 * g5) + g2 * (g5 * g4 - g2 * g6) + g4 * (g2 * g5 - g3 * g4);
            if d <= f64::MIN_POSITIVE {
                return None;
            }
            let s = 1.0 / d;
            m[0] = s * (g3 * g6 - g5 * g5);
            m[1] = -s * (g2 * g6 - g4 * g5);
            m[2] = s * (g1 * g6 - g4 * g4);
            m[3] = s * (g2 * g5 - g3 * g4);
            m[4] = -s * (g1 * g5 - g2 * g4);
            m[5] = s * (g1 * g3 - g2 * g2);
            d
        }
    };
    let mut ltog = [[0.0f64; 3]; 3];
    for (i, row) in ltog.iter_mut().enumerate().take(cdim) {
        let x = dx[i];
        match dim {
            1 => row[0] = x[0] * m[0],
            2 => {
                row[0] = x[0] * m[0] + x[1] * m[1];
                row[1] = x[0] * m[1] + x[1] * m[2];
            }
            _ => {
                row[0] = x[0] * m[0] + x[1] * m[1] + x[2] * m[3];
                row[1] = x[0] * m[1] + x[1] * m[2] + x[2] * m[4];
                row[2] = x[0] * m[3] + x[1] * m[4] + x[2] * m[5];
            }
        }
    }
    Some((det.sqrt(), ltog))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same basis, gradients and Jacobian as the scalar `ElementInfo`,
    /// up to rounding, on a distorted quadrilateral and a tetrahedron.
    #[test]
    fn agrees_with_element_info() {
        let cases: [(u16, Vec<[f64; 3]>, [f64; 3]); 3] = [
            (
                404,
                vec![
                    [0.1, 0.0, 0.0],
                    [1.2, 0.1, 0.0],
                    [1.0, 0.9, 0.0],
                    [0.0, 1.1, 0.0],
                ],
                [0.3, -0.2, 0.0],
            ),
            (
                504,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.1, 0.0],
                    [0.2, 1.0, 0.1],
                    [0.1, 0.2, 1.3],
                ],
                [0.2, 0.3, 0.1],
            ),
            (
                202,
                vec![[0.0, 0.0, 0.0], [0.7, 0.4, 0.0]],
                [0.25, 0.0, 0.0],
            ),
        ];
        for (code, xyz, uvw) in cases {
            let cdim = if code == 504 { 3 } else { 2 };
            let et = ElementType(code);
            let (mut a, mut b) = (BasisEval::default(), BasisEval::default());
            ElementInfoVec { cdim }.eval(et, &xyz, uvw, &mut a).unwrap();
            elmer_fe::LagrangeBasis::new(cdim)
                .eval(et, &xyz, uvw, &mut b)
                .unwrap();
            assert!((a.det_j - b.det_j).abs() < 1e-14 * b.det_j.abs(), "{code}");
            for (x, y) in a.basis.iter().zip(&b.basis) {
                assert_eq!(x, y);
            }
            for (x, y) in a.dbasis_dx.iter().zip(&b.dbasis_dx) {
                for k in 0..3 {
                    assert!(
                        (x[k] - y[k]).abs() < 1e-13 * (1.0 + y[k].abs()),
                        "{code}: {x:?} {y:?}"
                    );
                }
            }
        }
    }
}
