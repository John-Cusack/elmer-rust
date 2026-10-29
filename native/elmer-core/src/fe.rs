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
}
