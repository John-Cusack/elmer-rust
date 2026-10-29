//! Mesh2D against the pristine `Mesh2D` binary: the four files it writes,
//! line by line. `tests/data/mesh2d/<case>` holds the `.mif` input and the
//! `expected/` files of the binary run on it.

use std::path::{Path, PathBuf};

use elmer_mesh_ops::mesh2d;

fn data(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mesh2d")
        .join(case)
}

fn fixture(case: &str) {
    let dir = data(case);
    let mif = std::fs::read_to_string(dir.join("in.mif")).unwrap();
    let got = mesh2d(&mif).unwrap();
    for (file, got) in [
        ("mesh.header", &got.header),
        ("mesh.nodes", &got.nodes),
        ("mesh.elements", &got.elements),
        ("mesh.boundary", &got.boundary),
    ] {
        let want = std::fs::read_to_string(dir.join("expected").join(file)).unwrap();
        let (got, want): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        assert_eq!(got.len(), want.len(), "{case}: {file} line count");
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g, w, "{case}: {file} line {}", i + 1);
        }
    }
}

#[test]
fn l_domain_voronoi_vertex_parabolic() {
    fixture("Ldomain");
}

#[test]
fn corner_moving_front() {
    fixture("kulma");
}

#[test]
fn gmultg_mesh() {
    fixture("GmultgMesh1");
}

#[test]
fn heateq_mesh() {
    fixture("HeateqMesh1");
}

#[test]
fn circle_in_box() {
    fixture("CircleInBox");
}

#[test]
fn unsupported_layers_are_refused() {
    let mif = std::fs::read_to_string(data("kulma").join("in.mif"))
        .unwrap()
        .replace("MovingFront", "QuadGrid");
    assert!(matches!(
        mesh2d(&mif),
        Err(elmer_mesh_ops::MeshOpsError::Unsupported(_))
    ));
}

#[test]
fn malformed_input_is_invalid() {
    assert!(matches!(
        mesh2d("Geometry2D:\n H: x"),
        Err(elmer_mesh_ops::MeshOpsError::Invalid(_))
    ));
}
