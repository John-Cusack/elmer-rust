//! `ElementDiameter(Element, Nodes)` (ElementBasis.F90:3051), the
//! `Element % hK` of the pressure stabilization. The same function as
//! `elmer-fe`'s `stab::element_diameter` with `UseLongEdge = .FALSE.`,
//! which this crate does not depend on.

use elmer_core::ElementType;

/// `GetEdgeMap` (ElementBasis.F90:2937), 1-based node pairs.
fn edge_map(family: u16) -> &'static [[usize; 2]] {
    match family {
        5 => &[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]],
        6 => &[[1, 2], [2, 3], [4, 3], [1, 4], [1, 5], [2, 5], [3, 5], [4, 5]],
        7 => &[[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4], [1, 4], [2, 5], [3, 6]],
        _ => &[[1, 2], [2, 3], [4, 3], [1, 4], [5, 6], [6, 7], [8, 7], [5, 8], [1, 5], [2, 6], [3, 7], [4, 8]],
    }
}

/// `hK` of an element from its node coordinates.
pub(crate) fn element_diameter(etype: ElementType, xyz: &[[f64; 3]]) -> f64 {
    let (x, y, z) = (|i: usize| xyz[i - 1][0], |i: usize| xyz[i - 1][1], |i: usize| xyz[i - 1][2]);
    let sq = |d: f64| d * d;
    let hk = match etype.family() {
        1 => 0.0,
        2 => sq(x(2) - x(1)) + sq(y(2) - y(1)) + sq(z(2) - z(1)),
        3 => {
            let j11 = x(2) - x(1);
            let j12 = y(2) - y(1);
            let j13 = z(2) - z(1);
            let j21 = x(3) - x(1);
            let j22 = y(3) - y(1);
            let j23 = z(3) - z(1);
            let g11 = j11 * j11 + j12 * j12 + j13 * j13;
            let g12 = j11 * j21 + j12 * j22 + j13 * j23;
            let g22 = j21 * j21 + j22 * j22 + j23 * j23;
            let a = (g11 * g22 - g12 * g12).sqrt() / 2.0;
            let cx = (x(1) + x(2) + x(3)) / 3.0;
            let cy = (y(1) + y(2) + y(3)) / 3.0;
            let cz = (z(1) + z(2) + z(3)) / 3.0;
            let mut s = sq(x(1) - cx) + sq(y(1) - cy) + sq(z(1) - cz);
            s = s + sq(x(2) - cx) + sq(y(2) - cy) + sq(z(2) - cz);
            s = s + sq(x(3) - cx) + sq(y(3) - cy) + sq(z(3) - cz);
            16.0 * a * a / (3.0 * s)
        }
        4 => {
            let cx = sq(x(2) - x(1)) + sq(y(2) - y(1)) + sq(z(2) - z(1));
            let cy = sq(x(4) - x(1)) + sq(y(4) - y(1)) + sq(z(4) - z(1));
            2.0 * cx * cy / (cx + cy)
        }
        f => {
            let mut hk = f64::MAX;
            for &[j, k] in edge_map(f) {
                let d = sq(x(j) - x(k)) + sq(y(j) - y(k)) + sq(z(j) - z(k));
                if d < hk {
                    hk = d;
                }
            }
            hk
        }
    };
    hk.sqrt()
}
