//! Inputs the operations leave alone or refuse.

use elmer_core::{ElemId, ElementType, ParamMap, Value};
use elmer_mesh_levels::Element;
use elmer_mesh_ops::{create_discont_mesh, increase_element_order, Mesh, MeshOpsError};

/// Two triangles sharing the edge 1-2 and the boundary element on it.
fn two_triangles() -> Mesh {
    let tri = |nodes: [u32; 3], tag| Element {
        etype: ElementType(303),
        nodes: nodes.to_vec(),
        tag,
        parents: [None, None],
    };
    Mesh {
        dim: 2,
        coords: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        elements: vec![
            tri([0, 1, 2], 1),
            tri([1, 3, 2], 2),
            Element {
                etype: ElementType(202),
                nodes: vec![1, 2],
                tag: 1,
                parents: [Some(ElemId(0)), Some(ElemId(1))],
            },
        ],
        n_bulk: 2,
    }
}

fn bc(discontinuous: bool) -> ParamMap {
    let mut bc = ParamMap::new();
    bc.insert("Discontinuous Boundary", Value::Logical(discontinuous));
    bc
}

#[test]
fn no_discontinuous_bc_leaves_the_mesh() {
    let mesh = two_triangles();
    let sim = ParamMap::new();
    assert!(create_discont_mesh(&mesh, &sim, &[bc(false)])
        .unwrap()
        .is_none());
    assert!(create_discont_mesh(&mesh, &sim, &[]).unwrap().is_none());
}

#[test]
fn interface_nodes_of_the_moving_body_are_copied() {
    let mesh = two_triangles();
    let mut b = bc(true);
    b.insert("Discontinuous Target Bodies", Value::Integer(vec![2]));
    let out = create_discont_mesh(&mesh, &ParamMap::new(), &[b])
        .unwrap()
        .unwrap();
    assert!(out.discont_mesh);
    assert_eq!(out.n_discont_nodes, 2);
    assert_eq!(out.perm, vec![0, 1, 2, 0]);
    assert_eq!(out.mesh.coords.len(), 6);
    assert_eq!(out.mesh.coords[4], out.mesh.coords[1]);
    assert_eq!(out.mesh.coords[5], out.mesh.coords[2]);
    assert_eq!(out.mesh.elements[0].nodes, [0, 1, 2]);
    assert_eq!(out.mesh.elements[1].nodes, [4, 3, 5]);
    // The moving parent ends up on the right.
    assert_eq!(
        out.mesh.elements[2].parents,
        [Some(ElemId(0)), Some(ElemId(1))]
    );
}

#[test]
fn the_parent_body_common_to_all_interface_elements_moves() {
    // Left is body 1 and right body 2 on the only interface element.
    let mesh = two_triangles();
    let out = create_discont_mesh(&mesh, &ParamMap::new(), &[bc(true)])
        .unwrap()
        .unwrap();
    assert_eq!(out.mesh.elements[1].nodes, [4, 3, 5]);
}

#[test]
fn doubling_needs_a_target_boundary() {
    let mesh = two_triangles();
    let mut b = bc(true);
    b.insert("Discontinuous BC", Value::Integer(vec![0]));
    let mut b2 = bc(true);
    b2.insert("Mortar BC", Value::Integer(vec![1]));
    // The doubled BC names the target 0, and the other one is never asked.
    let err = create_discont_mesh(&mesh, &ParamMap::new(), &[b, b2]).unwrap_err();
    assert!(matches!(err, MeshOpsError::Invalid(_)), "{err}");
}

#[test]
fn higher_order_input_is_refused() {
    let mut mesh = two_triangles();
    mesh.elements[0].etype = ElementType(306);
    mesh.elements[0].nodes.extend([3, 3, 3]);
    assert!(matches!(
        increase_element_order(&mesh),
        Err(MeshOpsError::Invalid(_))
    ));
}

#[test]
fn pyramids_are_refused() {
    let mut mesh = two_triangles();
    mesh.dim = 3;
    mesh.elements[0].etype = ElementType(605);
    mesh.elements[0].nodes = vec![0, 1, 2, 3, 0];
    assert!(matches!(
        increase_element_order(&mesh),
        Err(MeshOpsError::Unsupported(_))
    ));
}

#[test]
fn edges_are_numbered_once() {
    let out = increase_element_order(&two_triangles()).unwrap();
    // 4 corners and 5 edges; the shared edge and its boundary element have one node.
    assert_eq!(out.coords.len(), 9);
    assert!(out.bulk().iter().all(|e| e.etype.0 == 306));
    assert_eq!(out.boundary()[0].etype.0, 203);
    assert_eq!(
        out.coords[out.boundary()[0].nodes[2] as usize],
        [0.5, 0.5, 0.0]
    );
}
