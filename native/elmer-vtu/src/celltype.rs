//! Elmer element code to VTK cell type and node order
//! (`Elmer2VtkElement` and `Elmer2VtkIndexes`, SaveUtils.F90:60-263).

use elmer_core::ElementType;

/// 0-based form of `Order416` (SaveUtils.F90:157).
const ORDER_416: [usize; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 12, 13, 15, 14];
/// 0-based form of `Order820` (SaveUtils.F90:160).
const ORDER_820: [usize; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15,
];
/// 0-based form of `Order827` (SaveUtils.F90:163).
const ORDER_827: [usize; 27] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15, 23, 21, 20, 22, 24, 25,
    26,
];

/// VTK cell type for an Elmer element code, `SaveLinear = .FALSE.`.
///
/// `None` where `Elmer2VtkElement` calls `Fatal` ("Not implemented for
/// elementtype").
fn vtk_type_full(code: u16) -> Option<u8> {
    Some(match code {
        101 => 1,
        202 => 3,
        203 => 21,
        204 => 35,
        205..=209 => 68,
        303 => 5,
        306 => 22,
        310 | 315 | 321 | 328 | 336 | 345 => 69,
        404 => 9,
        408 => 23,
        409 => 28,
        416 | 425 | 436 | 449 | 464 | 481 => 70,
        504 => 10,
        510 => 24,
        520 | 535 => 71,
        605 => 14,
        613 => 27,
        706 => 13,
        715 => 26,
        718 | 740 | 775 => 73,
        808 => 12,
        820 => 25,
        827 => 29,
        864 | 8125 | 8216 | 8343 | 8512 | 8729 => 72,
        _ => return None,
    })
}

/// `Elmer2VtkElement(code, SaveLinear)`: the first-order VTK cell when
/// `save_linear` is set, otherwise the cell of the element's own order.
pub fn vtk_cell_type(etype: ElementType, save_linear: bool) -> Option<u8> {
    let mut t = vtk_type_full(etype.0)?;
    if save_linear {
        t = match t {
            21 | 35 => 3,
            22 | 69 => 5,
            23 | 28 | 70 => 9,
            24 => 10,
            27 => 14,
            26 => 13,
            25 | 29 => 12,
            other => other,
        };
    }
    Some(t)
}

/// Node order of `etype` in VTK, as 0-based Elmer local indices; an empty
/// slice means the order is unchanged. Only 416, 820 and 827 are reordered
/// (and never when saving linear elements).
pub fn vtk_node_order(etype: ElementType) -> &'static [usize] {
    match etype.0 {
        416 => &ORDER_416,
        820 => &ORDER_820,
        827 => &ORDER_827,
        _ => &[],
    }
}

/// Number of nodes written per element: the corner count when saving linear
/// elements (`GetElementCorners`, DefUtils.F90:2016), else all of them.
pub fn written_node_count(etype: ElementType, save_linear: bool) -> usize {
    if save_linear {
        let mut c = (etype.0 / 100) as usize;
        if (5..=7).contains(&c) {
            c -= 1;
        }
        c
    } else {
        etype.n_nodes()
    }
}

/// VTK cell type and node reorder for an Elmer element (the full-order form
/// of `Elmer2VtkElement` + `Elmer2VtkIndexes`). The slice is empty when the
/// node order is unchanged. `None` for element codes Elmer cannot write.
pub fn elmer_to_vtk_type(etype: ElementType) -> Option<(u8, &'static [usize])> {
    Some((vtk_cell_type(etype, false)?, vtk_node_order(etype)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_cells() {
        let et = |c| ElementType(c);
        assert_eq!(elmer_to_vtk_type(et(303)), Some((5, &[][..])));
        assert_eq!(elmer_to_vtk_type(et(808)).unwrap().0, 12);
        assert_eq!(elmer_to_vtk_type(et(605)).unwrap().0, 14);
        assert_eq!(elmer_to_vtk_type(et(706)).unwrap().0, 13);
        assert_eq!(elmer_to_vtk_type(et(101)).unwrap().0, 1);
        assert_eq!(elmer_to_vtk_type(et(999)), None);
    }

    #[test]
    fn reorders_are_permutations() {
        for (code, n) in [(416u16, 16usize), (820, 20), (827, 27)] {
            let mut o = vtk_node_order(ElementType(code)).to_vec();
            assert_eq!(o.len(), n);
            o.sort_unstable();
            assert_eq!(o, (0..n).collect::<Vec<_>>());
        }
        // Spot checks against the Fortran 1-based lists.
        assert_eq!(vtk_node_order(ElementType(416))[8], 9);
        assert_eq!(vtk_node_order(ElementType(820))[12], 16);
        assert_eq!(vtk_node_order(ElementType(827))[20], 23);
    }

    #[test]
    fn linear_collapse_and_corners() {
        let et = |c| ElementType(c);
        assert_eq!(vtk_cell_type(et(510), true), Some(10));
        assert_eq!(vtk_cell_type(et(820), true), Some(12));
        assert_eq!(vtk_cell_type(et(613), true), Some(14));
        assert_eq!(vtk_cell_type(et(715), true), Some(13));
        assert_eq!(vtk_cell_type(et(203), true), Some(3));
        assert_eq!(vtk_cell_type(et(310), true), Some(5));
        assert_eq!(vtk_cell_type(et(303), true), Some(5));
        assert_eq!(written_node_count(et(510), true), 4);
        assert_eq!(written_node_count(et(613), true), 5);
        assert_eq!(written_node_count(et(715), true), 6);
        assert_eq!(written_node_count(et(827), true), 8);
        assert_eq!(written_node_count(et(827), false), 27);
    }
}
