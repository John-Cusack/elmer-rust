//! Refined and extruded meshes against the pristine ElmerSolver's own
//! `SaveMesh` output (`WriteMeshToDisk2`): node coordinates as `ES17.10`,
//! element and boundary connectivity, tags and parents, line by line.
//!
//! `tests/data/<case>` holds the input mesh, `expected/` and a `case.sif`
//! reduced to the keywords the mesh step reads. The quick-test cases are
//! named after the test; `TriExtrude` is ExtrudedMesh6 on a triangulated
//! square, and `*Levels` are `Mesh Levels` runs on a hexahedral, wedge
//! (TriExtrude's output) and tetrahedral (IncreaseOrderTet's angle.msh)
//! mesh, since no serial quick test refines a 3D mesh.

mod support;

use std::path::{Path, PathBuf};

fn data(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
}

fn fixture(case: &str) {
    let dir = data(case);
    support::check(&dir, "case.sif", &dir.join("expected"));
}

#[test]
fn extruded_mesh_layers() {
    fixture("ExtrudedMesh");
}

#[test]
fn extruded_mesh_ratio() {
    fixture("ExtrudedMesh6");
}

#[test]
fn extruded_mesh_symmetric_ratio() {
    fixture("ExtrudedMesh7");
}

#[test]
fn extruded_bc_layers() {
    fixture("ExtrudedMesh8");
}

#[test]
fn extruded_preserve_baseline() {
    fixture("ExtrudedMesh9");
}

#[test]
fn extruded_triangles_to_wedges() {
    fixture("TriExtrude");
}

#[test]
fn mesh_levels_quads() {
    fixture("MeshLevelsRestart");
}

#[test]
fn mesh_levels_hexahedra() {
    fixture("Hex808Levels");
}

#[test]
fn mesh_levels_wedges() {
    fixture("Wedge706Levels");
}

#[test]
fn mesh_levels_tetrahedra() {
    fixture("Tet504AngleLevels");
}

/// The larger references, made by the lane's capture scripts outside the
/// repository: `ELMER_MESH_LEVELS_REF=<dir>` with `<dir>/<test>/ref.sif`,
/// its mesh DB and `saved/<test>/`. Skipped without the variable.
#[test]
fn quick_test_references() {
    let Some(root) = std::env::var_os("ELMER_MESH_LEVELS_REF") else {
        eprintln!("ELMER_MESH_LEVELS_REF not set; skipping");
        return;
    };
    let root = PathBuf::from(root);
    for case in [
        "1dtests_MeshLevels",
        "AngleMetisLayer",
        "CoilSolverLoop",
        "SunAngle",
        "TopoOptHeat2DCompMin2",
        "fluxsolver-multimesh",
        "ExtrudedMesh3",
        "ExtrudedMesh4",
        "FixTangentVelo",
        "Tet504Levels",
    ] {
        let dir = root.join(case);
        let saved = dir.join("saved");
        let name = std::fs::read_dir(&saved)
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .file_name();
        support::check(&dir, "ref.sif", &saved.join(name));
        eprintln!("{case}: identical");
    }
}
