//! `ElementInfo` (ElemInfo.F90:331) for Lagrange elements: basis values,
//! global gradients and `detJ = sqrt(det(J^T J))` at one reference point.

use elmer_core::{BasisEval, ElementBasis, ElementType};

use crate::basis::{nodal_basis, nodal_first_derivatives};
use crate::element_type::{element_type, ElementTypeDef};
use crate::error::FeError;
use crate::metric::{check_metric, element_metric, MetricStatus};

/// `ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )`.
///
/// `cdim` is `CoordinateSystemDimension()`, `status` the element's
/// `Element % Status` (fixed by the first call on an element). On success
/// `out.basis` and `out.dbasis_dx` hold `n` entries; gradient components past
/// `cdim` are zero, as Elmer zeroes `dBasisdx` before filling it.
pub fn element_info(
    def: &ElementTypeDef,
    xyz: &[[f64; 3]],
    uvw: [f64; 3],
    cdim: usize,
    status: &mut MetricStatus,
    out: &mut BasisEval,
) -> Result<(), FeError> {
    let n = def.n_nodes;
    if xyz.len() != n {
        return Err(FeError::NodeCount {
            etype: def.code,
            expected: n,
            got: xyz.len(),
        });
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
    if def.code < 200 {
        return Err(FeError::Unsupported(def.code));
    }
    let [u, v, w] = uvw;
    let mut dl = vec![[0.0; 3]; n];
    nodal_basis(def, u, v, w, &mut out.basis);
    nodal_first_derivatives(def, u, v, w, &mut dl);

    let dim = def.dim;
    if *status == MetricStatus::Unchecked {
        *status = if check_metric(xyz, &dl, dim, cdim) {
            MetricStatus::Double
        } else {
            MetricStatus::Quad
        };
    }
    let m = element_metric(xyz, &dl, dim, cdim, *status)?;
    out.det_j = m.det_g;
    for k in 0..dim {
        for j in 0..cdim {
            for (d, l) in out.dbasis_dx.iter_mut().zip(&dl) {
                d[j] += l[k] * m.ltog[j][k];
            }
        }
    }
    Ok(())
}

/// [`ElementBasis`] for the Lagrange element types of elements.def.
///
/// `cdim` plays `CoordinateSystemDimension()`. Each [`ElementBasis::eval`]
/// call is treated as the first `ElementInfo` call on its element, so the
/// double/quad metric choice is made from that point; callers that loop over
/// points of one element and need Elmer's once-per-element choice use
/// [`element_info`] with a [`MetricStatus`] kept per element.
#[derive(Copy, Clone, Debug)]
pub struct LagrangeBasis {
    pub cdim: usize,
}

impl LagrangeBasis {
    pub fn new(cdim: usize) -> Self {
        LagrangeBasis { cdim }
    }
}

impl Default for LagrangeBasis {
    fn default() -> Self {
        LagrangeBasis { cdim: 3 }
    }
}

impl ElementBasis for LagrangeBasis {
    fn eval(
        &self,
        etype: ElementType,
        nodal_xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let def = element_type(etype.0).map_err(|e| e.to_string())?;
        let mut status = MetricStatus::Unchecked;
        element_info(def, nodal_xyz, uvw, self.cdim, &mut status, out).map_err(|e| e.to_string())
    }
}
