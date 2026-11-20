//! `Discontinuous Galerkin` with `DG Reduced Basis`: two quadrilaterals in
//! two bodies sharing the edge of nodes 2 and 5. Worked out by hand from
//! InitialPermutation (Lists.F90:271) and MakeListMatrix (ElementUtils.F90).

mod support;

use elmer_core::{ElemId, MeshView, ParamMap, Value, NO_DOF};
use elmer_pdofs::PDofsBackend;
use elmer_solver::{DofsBackend, Model};
use support::TestMesh;

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn mesh() -> TestMesh {
    let dir = std::env::temp_dir().join(format!("elmer-pdofs-dg-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("mesh.nodes"),
        "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 2 0 0\n4 -1 0 1 0\n5 -1 1 1 0\n6 -1 2 1 0\n",
    )
    .unwrap();
    std::fs::write(dir.join("mesh.elements"), "1 1 404 1 2 5 4\n2 2 404 2 3 6 5\n").unwrap();
    std::fs::write(dir.join("mesh.boundary"), "1 1 1 2 202 2 5\n2 2 1 0 202 1 4\n").unwrap();
    let m = TestMesh::read(&dir, 2);
    std::fs::remove_dir_all(&dir).ok();
    m
}

#[test]
fn reduced_basis_numbers_each_group_by_node_and_couples_the_interface() {
    let mesh = mesh();
    let material = sec(&[("DG Parent Material", Value::Logical(true))]);
    let body = |b: i32| sec(&[("Equation", Value::Integer(vec![1])), ("Material", Value::Integer(vec![b]))]);
    let model = Model {
        solvers: vec![sec(&[("Equation", Value::Str("heat".into()))])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![body(1), body(1)],
        materials: vec![material],
        boundary_conditions: vec![sec(&[]), sec(&[])],
        ..Model::default()
    };
    let solver = sec(&[
        ("Discontinuous Galerkin", Value::Logical(true)),
        ("DG Reduced Basis", Value::Logical(true)),
        ("DG Reduced Basis Master Bodies", Value::Integer(vec![1])),
        ("DG Indirect Connections", Value::Logical(true)),
    ]);
    let active = [ElemId(0), ElemId(1)];
    let d = PDofsBackend::new().for_solver_id(&mesh, &model, 1, &solver, &active, 1).unwrap().unwrap();
    // Group 1 (body 1) numbers nodes 1, 2, 5, 4; group 2 (body 2) starts
    // again: nodes 2, 3, 6, 5. The Perm has MaxDGDOFs * 4 elements entries.
    let p = d.perm();
    assert_eq!(p.n_slots, 8);
    assert_eq!(p.perm.len(), 16);
    assert_eq!(&p.perm[..8], &[0, 1, 2, 3, 4, 5, 6, 7]);
    assert!(p.perm[8..].iter().all(|&s| s == NO_DOF));
    // The interface joins both elements' dofs.
    let a = d.create_matrix(1);
    assert_eq!(a.nrows, 8);
    assert_eq!(a.values.len(), 64);
    // A boundary element's slots are its parents' (left, then right); its DG
    // indexes are inherited from the left parent (DG Parent Material).
    let mut s = Vec::new();
    d.element_slots(ElemId(2), &mut s);
    assert_eq!(s, [0, 1, 2, 3, 4, 5, 6, 7]);
    let ix = d.dg_layout().unwrap();
    assert_eq!(ix.index[2], [1, 2]);
    assert_eq!(ix.index[3], [0, 3]);
    // The edge between the bodies, its left parent the first to visit it.
    let faces = elmer_solver::dg_interfaces(&mesh).unwrap();
    assert_eq!(faces, [(ElemId(0), ElemId(1), vec![1, 4])]);
    assert_eq!(mesh.n_bulk(), 2);
}
