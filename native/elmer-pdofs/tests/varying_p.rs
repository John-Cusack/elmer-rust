//! `fem/tests/p-FEM_with_varying_p` against the pristine captures in
//! `tests/data/varying_p`: the MATC degree `p:%bodywisep` (4 in body 1, 8
//! in body 2), the edge and bubble slots it implies, Perm and the matrix
//! pattern of the 5-component shell system.

mod support;

use elmer_assembly::PermOptions;
use elmer_pdofs::{
    create_p_system, initial_max_defs, DofLayout, ElementDefs, Kind, PMesh, PMeshOptions, PSystem,
    SolverElements,
};
use support::{data_dir, lines, read_perm, TestMesh};

const DEF: &str = "p:%bodywisep";

fn setup() -> (PMesh, ElementDefs, PSystem) {
    let mesh = TestMesh::read(&data_dir("varying_p"), 2);
    let mut max = initial_max_defs();
    let mut defs = ElementDefs::from_solver_element(Some(DEF), 2, false, &mut max).unwrap();
    let none = |_: usize| None;
    let mut pm = {
        let mut solvers = [SolverElements {
            defs: &mut defs,
            element: Some(DEF),
            dg: false,
            body_element: &none,
        }];
        PMesh::new(
            &mesh,
            &mut solvers,
            &mut max,
            PMeshOptions {
                serendipity: false,
                need_edges: None,
                n_bodies: 2,
                body_has_equation: &|_| true,
                // $ function bodywisep(cx) { if (cx(0)>1) {_bodywisep=8;} else {_bodywisep=4;} }
                matc_degree: &mut |f, cx| {
                    assert_eq!(f, "bodywisep(cx)");
                    Ok(if cx[0] > 1.0 { 8.0 } else { 4.0 })
                },
                boundary_body: &|_| 0,
            },
        )
        .unwrap()
    };
    let n_bulk = pm.n_bulk;
    let active = move |e: usize| e < n_bulk;
    let sys = create_p_system(
        &mut pm,
        &defs,
        true,
        false,
        &active,
        PermOptions::default(),
        5,
    );
    (pm, defs, sys)
}

#[test]
fn degrees_and_dof_counts() {
    let (pm, defs, sys) = setup();
    assert_eq!(defs.get(4, 1, Kind::P), 4);
    assert_eq!(defs.get(4, 2, Kind::P), 8);
    assert_eq!(pm.edges.len(), 12);
    assert_eq!((pm.max_edge_dofs, pm.max_bdofs), (7, 49));
    // Edges next to body 2 carry 7 DOFs, the others 3.
    let mut counts: Vec<usize> = pm.edges.iter().map(|e| e.bdofs).collect();
    counts.sort_unstable();
    assert_eq!(counts, [3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 7, 7]);
    // Ndeg = 9 nodes + 12 x 7 edge slots + 4 x 49 bubble slots.
    assert_eq!(sys.perm.perm.len(), 289);
    // 9 + 5 x 3 + 7 x 7 + 2 x 9 + 2 x 49 used.
    assert_eq!(sys.perm.n_dofs, 189);
    assert_eq!(sys.matrix.nrows, 945);
    assert_eq!(sys.matrix.nnz(), 363_425);
    assert_eq!(sys.bandwidth.initial, 187);
    assert_eq!(sys.bandwidth.optimized, Some(187));

    let lay = DofLayout::new(&pm, &defs, true);
    let mut dofs = Vec::new();
    for (e, nd) in [(0, 4 + 4 * 3 + 9), (3, 4 + 4 * 7 + 49)] {
        lay.element_dofs(e, &mut dofs);
        assert_eq!(dofs.len(), nd, "element {e}");
        let spec = lay.element_spec(e).unwrap();
        assert_eq!(spec.n_basis().unwrap(), nd, "element {e}");
    }
    // The body 1 element shared edge uses the first 3 of its 7 slots.
    let shared = pm.elements[0].edges[2];
    assert_eq!(pm.edges[shared as usize].bdofs, 7);
    lay.element_dofs(0, &mut dofs);
    let first = lay.edge_slot(shared, 1);
    assert_eq!(&dofs[4 + 2 * 3..4 + 3 * 3], &[first, first + 1, first + 2]);
}

#[test]
fn perm_and_pattern_match() {
    let (_, _, sys) = setup();
    let elmer = read_perm(&data_dir("varying_p").join("perm.dat"));
    let ours: Vec<u32> = sys
        .perm
        .perm
        .iter()
        .map(|&p| if p == elmer_core::NO_DOF { 0 } else { p + 1 })
        .collect();
    assert_eq!(ours, elmer);

    let a = &sys.matrix;
    let rows = lines(&data_dir("varying_p").join("pattern.dat"));
    assert_eq!(rows.len(), sys.perm.n_dofs);
    for (i, l) in rows.iter().enumerate() {
        let (head, cols) = l.split_once(':').unwrap();
        assert_eq!(head.parse::<usize>().unwrap(), i + 1);
        let mut want = Vec::new();
        for r in cols.split_whitespace() {
            let (lo, hi) = r.split_once('-').unwrap_or((r, r));
            want.extend(lo.parse::<u32>().unwrap() - 1..hi.parse::<u32>().unwrap());
        }
        // Component 0 of slot i, block columns.
        let row = 5 * i;
        let got: Vec<u32> = a.col_idx[a.row_ptr[row]..a.row_ptr[row + 1]]
            .iter()
            .filter(|&&c| c % 5 == 0)
            .map(|&c| c / 5)
            .collect();
        assert_eq!(got, want, "row {}", i + 1);
        assert_eq!(a.row_ptr[row + 1] - a.row_ptr[row], 5 * want.len());
    }
}
