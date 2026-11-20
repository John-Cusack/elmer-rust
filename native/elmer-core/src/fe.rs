//! Consumer view of quadrature and basis evaluation; `elmer-fe` and
//! `elmer-fe-hp` implement these.

use crate::mesh::ElementType;

#[derive(Copy, Clone, Debug)]
pub struct QuadPoint {
    pub uvw: [f64; 3],
    pub weight: f64,
}

pub trait Quadrature {
    /// `np`: requested point count, or `None` for the element default.
    fn gauss_points(&self, etype: ElementType, np: Option<usize>) -> Vec<QuadPoint>;

    /// `GaussPoints(Element, RelOrder = rel)` of a nodal element (`rel` of
    /// -1, 0 or 1: `GaussPoints0`, `GaussPoints` or `GaussPoints2` of
    /// elements.def); `None` where the implementation has no such rules.
    fn gauss_points_rel(&self, etype: ElementType, rel: i32) -> Option<Vec<QuadPoint>> {
        (rel == 0).then(|| self.gauss_points(etype, None))
    }
}

#[derive(Clone, Debug, Default)]
pub struct BasisEval {
    pub det_j: f64,
    pub basis: Vec<f64>,
    pub dbasis_dx: Vec<[f64; 3]>,
}

pub trait ElementBasis {
    /// ElementInfo (ElemInfo.F90:331) equivalent for one point.
    fn eval(
        &self,
        etype: ElementType,
        nodal_xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String>;

    /// `ElementInfo(.., ddBasisddx, SecondDerivatives = .TRUE.)`: the
    /// global second derivatives of each basis function at one point.
    fn second_derivatives(
        &self,
        _etype: ElementType,
        _nodal_xyz: &[[f64; 3]],
        _uvw: [f64; 3],
    ) -> Result<Vec<[[f64; 3]; 3]>, String> {
        Err("second derivatives of this basis are not implemented".into())
    }
}
