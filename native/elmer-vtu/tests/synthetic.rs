//! Behaviour the Fortran fixtures do not exercise. Expected values are worked
//! out by hand from the Fortran source, not captured from a run: these check
//! local correctness only, not parity.

mod support;

use std::io::ErrorKind;

use elmer_core::{ElemId, ElementRef, ElementType, FieldKind, MeshView, ParamMap, Value};
use elmer_vtu::{render_pvtu, render_vtu, GeometryIds, VtuOptions};
use support::*;

/// Two triangles (bodies 1 and 2) and one boundary line on node 2-3 (tag 7).
fn mesh() -> TestMesh {
    TestMesh {
        coords: vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [9., 9., 9.]],
        elems: vec![
            (ElementType(303), vec![0, 1, 2], 1),
            (ElementType(303), vec![0, 2, 3], 2),
            (ElementType(202), vec![1, 2], 7),
        ],
        n_bulk: 2,
    }
}

fn ascii_opts(scalars: &[&str]) -> VtuOptions {
    let mut p = ParamMap::new();
    p.insert("Ascii Output", Value::Logical(true));
    for (i, s) in scalars.iter().enumerate() {
        p.insert(&format!("Scalar Field {}", i + 1), Value::Str((*s).into()));
    }
    VtuOptions::from_params(&p)
}

fn text(o: &VtuOptions, m: &TestMesh, f: &TestFields, group: Option<usize>) -> String {
    String::from_utf8(render_vtu(m, f, o, group).unwrap().unwrap()).unwrap()
}

#[test]
fn unreferenced_node_is_dropped_and_connectivity_renumbered() {
    // Node 5 (index 4) belongs to no element, so it is not saved; every
    // other node keeps its position.
    let f = TestFields::new(vec![OwnedField::nodal("t", 1, vec![1., 2., 3., 4., 5.])]);
    let out = text(&ascii_opts(&["t"]), &mesh(), &f, None);
    assert!(out.contains("NumberOfPoints=\"4\" NumberOfCells=\"3\""));
    assert!(!out.contains("9.0000000E+000"));
}

#[test]
fn perm_hides_absent_nodes_as_zero() {
    let mut fld = OwnedField::nodal("t", 1, vec![10., 20.]);
    fld.perm = Some(vec![0, elmer_core::NO_DOF, 1, elmer_core::NO_DOF, elmer_core::NO_DOF]);
    let out = text(&ascii_opts(&["t"]), &mesh(), &TestFields::new(vec![fld]), None);
    let a = parse_vtu(out.as_bytes());
    assert_eq!(a.array("PointData", Some("t")).reals, vec![10., 0., 20., 0.]);
}

#[test]
fn two_component_vector_pads_to_three_and_component_names_resolve() {
    let f = TestFields::new(vec![
        OwnedField::nodal("v 1", 1, vec![1., 2., 3., 4., 5.]),
        OwnedField::nodal("v 2", 1, vec![6., 7., 8., 9., 0.]),
    ]);
    let mut o = ascii_opts(&[]);
    o.vector_fields = vec!["v".into()];
    let a = parse_vtu(text(&o, &mesh(), &f, None).as_bytes());
    let v = a.array("PointData", Some("v"));
    assert_eq!(v.ncomp, 3);
    assert_eq!(v.reals, vec![1., 6., 0., 2., 7., 0., 3., 8., 0., 4., 9., 0.]);
}

#[test]
fn elemental_integration_point_and_boundary_parent_values() {
    // Elemental `e` (slot per element, dofs 1) and integration-point `g`
    // (Perm = offsets 0,2,5 -> 2 and 3 points). The boundary line is element 3
    // and takes its parent's (element 1) elemental value.
    let mut e = OwnedField::nodal("e", 1, vec![5.0, 6.0]);
    e.kind = FieldKind::Elemental;
    e.perm = Some(vec![0, 1]);
    let mut g = OwnedField::nodal("g", 1, vec![1.0, 3.0, 10.0, 20.0, 30.0]);
    g.kind = FieldKind::Ip;
    g.perm = Some(vec![0, 2, 5]);
    let f = TestFields::new(vec![e, g]);

    struct WithParent(TestMesh);
    impl MeshView for WithParent {
        fn dim(&self) -> usize { 3 }
        fn n_nodes(&self) -> usize { self.0.n_nodes() }
        fn coord(&self, n: u32) -> [f64; 3] { self.0.coord(n) }
        fn n_bulk(&self) -> usize { self.0.n_bulk() }
        fn n_boundary(&self) -> usize { self.0.n_boundary() }
        fn element(&self, id: ElemId) -> ElementRef<'_> {
            let mut r = self.0.element(id);
            if id.0 == 2 {
                r.parents = [Some(ElemId(0)), None];
            }
            r
        }
    }

    let mut o = ascii_opts(&["e", "g"]);
    o.scalar_fields = vec!["e".into(), "g".into()];
    let out = render_vtu(&WithParent(mesh()), &f, &o, None).unwrap().unwrap();
    let a = parse_vtu(&out);
    assert_eq!(a.array("CellData", Some("e")).reals, vec![5.0, 6.0, 5.0]);
    // Ip: means 2.0 and 20.0; the boundary element has no points (offset
    // array ends) and gets 0, where Fortran would divide 0 by 0.
    assert_eq!(a.array("CellData", Some("g")).reals, vec![2.0, 20.0, 0.0]);
    // Neither is point data.
    assert!(a.arrays.iter().all(|x| x.section != "PointData"));
}

#[test]
fn geometry_ids_and_part_groups() {
    let mut o = ascii_opts(&[]);
    o.save_geometry_ids = true;
    o.geometry = Some(GeometryIds { n_bodies: 2, bc_tags: vec![7], ..Default::default() });
    let f = TestFields::new(vec![]);
    let a = parse_vtu(text(&o, &mesh(), &f, None).as_bytes());
    // Offset 100 > 2 bodies; BC tag 7 is BC 1 -> 101.
    assert_eq!(a.array("CellData", Some("GeometryIds")).ints, vec![1, 2, 101]);

    // Group 2 = body 2 only: one cell, three nodes (0, 2, 3), renumbered.
    let a = parse_vtu(text(&o, &mesh(), &f, Some(2)).as_bytes());
    assert_eq!((a.npoints, a.ncells), (3, 1));
    assert_eq!(a.array("Cells", Some("connectivity")).ints, vec![0, 1, 2]);
    // Group 3 = the first BC: the line (constraint 7 != 1) is not selected, but
    // Fortran compares the tag to GroupId - NumberOfBodies = 1, so nothing is.
    assert!(render_vtu(&mesh(), &f, &o, Some(3)).unwrap().is_none());

    // Overrides and defaults.
    o.default_bc_id = Some(55);
    o.default_body_id = Some(9);
    let a = parse_vtu(text(&o, &mesh(), &f, None).as_bytes());
    assert_eq!(a.array("CellData", Some("GeometryIds")).ints, vec![9, 9, 55]);

    // Missing model facts are an error, not a guess.
    o.geometry = None;
    assert!(render_vtu(&mesh(), &f, &o, None).is_err());
}

#[test]
fn fixed_mesh_removes_displacement_then_scales() {
    let f = TestFields::new(vec![OwnedField::nodal(
        "Displacement",
        2,
        vec![0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 9.0],
    )]);
    let mut o = ascii_opts(&[]);
    o.fixed_mesh = true;
    o.coordinate_scale = [2.0, 1.0, 1.0];
    o.coordinate_offset = [0.0, 0.0, 1.0];
    let a = parse_vtu(text(&o, &mesh(), &f, None).as_bytes());
    let p = &a.array("Points", None).reals;
    // Node 0: x = 2*(0 - 0.5) = -1; z gets only the offset (2-dof field).
    assert_eq!(&p[..3], &[-1.0, 0.0, 1.0]);
    assert_eq!(&p[3..6], &[2.0, 0.0, 1.0]);
}

#[test]
fn unsupported_keywords_and_dg_fields_are_refused_not_approximated() {
    let f = TestFields::new(vec![]);
    let mut p = ParamMap::new();
    p.insert("Discontinuous Bodies", Value::Logical(true));
    let o = VtuOptions::from_params(&p);
    let e = render_vtu(&mesh(), &f, &o, None).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Unsupported);
    assert!(render_pvtu(&f, &o, &[]).is_err());

    let mut dg = OwnedField::nodal("d", 1, vec![0.0; 6]);
    dg.kind = FieldKind::Dg;
    let o = ascii_opts(&["d"]);
    let e = render_vtu(&mesh(), &TestFields::new(vec![dg]), &o, None).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Unsupported);
}

#[test]
fn nothing_selected_writes_nothing() {
    let mut o = ascii_opts(&[]);
    o.save_bulk_only = true;
    o.save_boundaries_only = true;
    assert!(render_vtu(&mesh(), &TestFields::new(vec![]), &o, None).unwrap().is_none());
}

#[test]
fn missing_variable_is_skipped() {
    let f = TestFields::new(vec![OwnedField::nodal("t", 1, vec![0.0; 5])]);
    let a = parse_vtu(text(&ascii_opts(&["nope", "t"]), &mesh(), &f, None).as_bytes());
    let names: Vec<_> = a.arrays.iter().filter(|x| x.section == "PointData").collect();
    assert_eq!(names.len(), 1);
}
