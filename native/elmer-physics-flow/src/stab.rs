//! Element size `hK` and the stabilization constant `mK` of Franca et al.,
//! as Elmer stores them in `Element % hK` and `Element % StabilizationMK`
//! (`MeshStabParams`, MeshBasics.F90:2370; `StabParam`, ElemInfo.F90:152;
//! `ElementDiameter`, ElementBasis.F90:3051).
//!
//! Only the linear-element branch of `StabParam` is ported: there `mK` is a
//! per-type constant and `hK` the geometric diameter. Higher-order elements
//! get `mK` from a generalized eigenproblem (`DSYGV`) and are refused. These
//! belong with the element geometry (`elmer-fe` / `MeshView`); this crate
//! computes them locally until that lane provides them.

use elmer_core::ElementType;

use crate::error::FlowError;

/// `Element % hK` and `Element % StabilizationMK` of one element.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StabParams {
    pub hk: f64,
    pub mk: f64,
}

/// Node pairs of each edge (Elmer `GetEdgeMap`, 1-based there).
fn edge_map(family: u16) -> &'static [(usize, usize)] {
    match family {
        2 => &[(0, 1)],
        5 => &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
        8 => &[
            (0, 1),
            (1, 2),
            (3, 2),
            (0, 3),
            (4, 5),
            (5, 6),
            (7, 6),
            (4, 7),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ],
        _ => &[],
    }
}

fn sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]) * (a[i] - b[i])).sum()
}

/// `ElementDiameter` (ElementBasis.F90:3051) for linear elements:
/// triangles and quadrilaterals use their closed forms, other families the
/// shortest (or, with `long_edge`, the longest) edge.
pub fn element_diameter(etype: ElementType, xyz: &[[f64; 3]], long_edge: bool) -> Result<f64, FlowError> {
    let family = etype.family();
    if xyz.len() < etype.n_nodes() {
        return Err(FlowError::NodalLength {
            what: "element nodes",
            expected: etype.n_nodes(),
            got: xyz.len(),
        });
    }
    let h2 = match family {
        1 => 0.0,
        3 => {
            let (x, y, z) = (xyz[0], xyz[1], xyz[2]);
            let j1 = [y[0] - x[0], y[1] - x[1], y[2] - x[2]];
            let j2 = [z[0] - x[0], z[1] - x[1], z[2] - x[2]];
            let g11: f64 = j1.iter().map(|v| v * v).sum();
            let g12: f64 = j1.iter().zip(&j2).map(|(a, b)| a * b).sum();
            let g22: f64 = j2.iter().map(|v| v * v).sum();
            let a = (g11 * g22 - g12 * g12).sqrt() / 2.0;
            let c: [f64; 3] = std::array::from_fn(|i| (x[i] + y[i] + z[i]) / 3.0);
            let s = sq(x, c) + sq(y, c) + sq(z, c);
            16.0 * a * a / (3.0 * s)
        }
        4 => {
            let cx = sq(xyz[1], xyz[0]);
            let cy = sq(xyz[3], xyz[0]);
            2.0 * cx * cy / (cx + cy)
        }
        f => {
            let edges = edge_map(f);
            if edges.is_empty() {
                return Err(FlowError::Unsupported(format!("element diameter of element type {}", etype.0)));
            }
            let lens = edges.iter().map(|&(a, b)| sq(xyz[a], xyz[b]));
            if long_edge {
                lens.fold(f64::MIN, f64::max)
            } else {
                lens.fold(f64::MAX, f64::min)
            }
        }
    };
    Ok(h2.sqrt())
}

/// `StabParam` with `hK` present, for elements of basis degree one.
pub fn stab_params(etype: ElementType, xyz: &[[f64; 3]], long_edge: bool) -> Result<StabParams, FlowError> {
    let mk = match etype.0 {
        202 | 303 | 404 | 504 | 605 | 706 => 1.0 / 3.0,
        808 => 1.0 / 6.0,
        t => {
            return Err(FlowError::Unsupported(format!(
                "stabilization parameter mK of element type {t} (StabParam's DSYGV branch for higher-order elements)"
            )))
        }
    };
    Ok(StabParams {
        hk: element_diameter(etype, xyz, long_edge)?,
        mk,
    })
}

/// Reference coordinates of the element nodes (`Element % TYPE % NodeU/V/W`).
pub(crate) fn reference_nodes(etype: ElementType) -> Result<&'static [[f64; 3]], FlowError> {
    Ok(match etype.0 {
        202 => &[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        303 => &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        404 => &[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]],
        504 => &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        808 => &[
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        t => return Err(FlowError::Unsupported(format!("element type {t}"))),
    })
}
