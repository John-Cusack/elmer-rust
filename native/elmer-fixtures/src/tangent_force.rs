//! `NormalTangentialDisplacement/TangentForce.f90`: the tangential load of
//! a friction boundary, `Mu = 0.5` times the normal force the previous
//! solve left in `Displacement Loads`, along the outward normal of the
//! current boundary element at its centre (`NormalVector(Element, Nodes,
//! 0, 0)` without the direction check).

use std::sync::Arc;

use elmer_expr::NodeFn;

use crate::ModelView;

/// `NormalVector` of a 2D line element at `u = 0`: `(-dy/du, dx/du)`
/// normalised, `NodalFirstDerivatives` being `-1/2, 1/2 (, 0)` there.
fn line_normal(xyz: &[[f64; 3]]) -> [f64; 3] {
    let d = [-0.5, 0.5, 0.0];
    let mut dxdu = 0.0;
    let mut dydu = 0.0;
    for (k, x) in xyz.iter().enumerate().take(3) {
        dxdu += d[k] * x[0];
        dydu += d[k] * x[1];
    }
    let det = dxdu * dxdu + dydu * dydu;
    if det <= 0.0 {
        return [0.0; 3];
    }
    let det = 1.0 / det.sqrt();
    [-dydu * det, dxdu * det, 0.0]
}

/// `TangentForce(Model, n, t)` bound to `model`.
pub fn tangent_force_fn(model: Arc<dyn ModelView>) -> NodeFn {
    Arc::new(move |node, _t| {
        let xyz = elmer_expr::current_coords();
        if xyz.len() < 2 {
            panic!("TangentForce: called without a current boundary element");
        }
        let normal = line_normal(&xyz);
        let mut force = [0.0; 3];
        for (i, f) in force.iter_mut().enumerate().take(2) {
            *f = model.value(&format!("displacement loads {}", i + 1), node).unwrap_or(0.0);
        }
        let normal_force = -(force[0] * normal[0] + force[1] * normal[1] + force[2] * normal[2]);
        0.5 * normal_force
    })
}
