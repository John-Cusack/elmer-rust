//! Keywords and meshes the port refuses instead of approximating, and the
//! unit-segment divisions on their own.

use elmer_core::{ElemId, ElementType, ParamMap, Value};
use elmer_mesh_levels::{
    division, extrude, mesh_levels, Element, ExtrudeOptions, Mesh, MeshLevelsError,
};

/// One unit square (404) with its four edges as tagged boundaries.
fn square() -> Mesh {
    let el = |etype, nodes: &[u32], tag, parent| Element {
        etype: ElementType(etype),
        nodes: nodes.to_vec(),
        tag,
        parents: [parent, None],
    };
    let p = Some(ElemId(0));
    Mesh {
        dim: 2,
        coords: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        elements: vec![
            el(404, &[0, 1, 2, 3], 1, None),
            el(202, &[0, 1], 1, p),
            el(202, &[1, 2], 2, p),
            el(202, &[2, 3], 1, p),
            el(202, &[3, 0], 2, p),
        ],
        n_bulk: 1,
    }
}

fn params(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn unsupported<T: std::fmt::Debug>(r: Result<T, MeshLevelsError>) {
    assert!(matches!(r, Err(MeshLevelsError::Unsupported(_))), "{r:?}");
}

#[test]
fn split_counts_and_parents() {
    let sim = params(&[("Mesh Levels", Value::Integer(vec![3]))]);
    let m = mesh_levels(&square(), &sim, &[]).unwrap().unwrap();
    assert_eq!((m.coords.len(), m.n_bulk, m.boundary().len()), (25, 16, 16));
    for b in m.boundary() {
        let parent = &m.elements[b.parents[0].unwrap().0 as usize];
        assert!(b.nodes.iter().all(|n| parent.nodes.contains(n)));
    }
    assert!(mesh_levels(&square(), &ParamMap::new(), &[])
        .unwrap()
        .is_none());
}

#[test]
fn mesh_levels_refusals() {
    let levels = ("Mesh Levels", Value::Integer(vec![2]));
    let follow = params(&[("Follow Circle Boundary", Value::Logical(true))]);
    unsupported(mesh_levels(
        &square(),
        &params(std::slice::from_ref(&levels)),
        &[follow],
    ));
    let grading = params(&[levels.clone(), ("Mesh Keep Grading", Value::Logical(true))]);
    unsupported(mesh_levels(&square(), &grading, &[]));
    let quads = params(&[levels.clone(), ("Split Mesh Quads", Value::Logical(true))]);
    unsupported(mesh_levels(&square(), &quads, &[]));
    let mut pyramid = square();
    pyramid.elements[0].etype = ElementType(605);
    unsupported(mesh_levels(&pyramid, &params(&[levels]), &[]));
}

#[test]
fn extrusion_refusals() {
    let levels = ("Extruded Mesh Levels", Value::Integer(vec![3]));
    for (k, v) in [
        ("Extruded Mesh Slices", Value::Logical(true)),
        ("Extruded Mesh Rotational", Value::Logical(true)),
        ("Extruded Mesh End Map", Value::real(1.0)),
    ] {
        unsupported(ExtrudeOptions::from_params(
            &params(&[levels.clone(), (k, v)]),
            2,
            None,
        ));
    }
    let elements = params(&[
        ("Extruded Mesh Ratio", Value::real(2.0)),
        ("Extruded Elements", Value::Integer(vec![2, 3])),
    ]);
    unsupported(ExtrudeOptions::from_params(&elements, 2, None));
    let density = Value::Dependent(elmer_core::Dependent {
        variables: vec!["coordinate 1".into()],
        kind: elmer_core::DependentKind::Matc("1+x".into()),
        rows: 1,
        cols: 1,
    });
    unsupported(ExtrudeOptions::from_params(
        &params(&[levels, ("Extruded Mesh Density", density)]),
        2,
        None,
    ));
    assert!(ExtrudeOptions::from_params(&ParamMap::new(), 2, None)
        .unwrap()
        .is_none());
}

#[test]
fn extruded_square_is_a_brick_column() {
    let sim = params(&[("Extruded Mesh Levels", Value::Integer(vec![3]))]);
    let opts = ExtrudeOptions::from_params(&sim, 2, None).unwrap().unwrap();
    assert_eq!(opts.w, [0.0, 0.5, 1.0]);
    let m = extrude(&square(), &opts).unwrap();
    // 2 bricks, 4 sides x 2 layers, bottom and top.
    assert_eq!((m.coords.len(), m.n_bulk, m.boundary().len()), (12, 2, 10));
    assert_eq!(m.elements[0].nodes, [0, 1, 2, 3, 4, 5, 6, 7]);
    // Bottom gets BC 3 (after the largest side BC 2), top BC 4.
    let planes: Vec<i32> = m.boundary()[8..].iter().map(|e| e.tag).collect();
    assert_eq!(planes, [3, 4]);
}

#[test]
fn geometric_divisions() {
    let w = division::geometric(4, 1.0).unwrap();
    assert_eq!(w, [0.0, 0.25, 0.5, 0.75, 1.0]);
    // Last element twice the first.
    let w = division::geometric(3, 2.0).unwrap();
    let h: Vec<f64> = w.windows(2).map(|p| p[1] - p[0]).collect();
    assert!((h[2] / h[0] - 2.0).abs() < 1e-12);
    // Symmetric about the middle for a negative ratio.
    let w = division::geometric(4, -0.5).unwrap();
    assert!((w[1] - (1.0 - w[3])).abs() < 1e-15 && w[2] == 0.5);
    assert_eq!(division::powi(3.0, 5), 243.0);
    assert_eq!(division::powi(2.0, -2), 0.25);
}
