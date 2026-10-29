//! `Discontinuous Boundary` against the pristine ElmerSolver's own
//! `SaveMesh` output (`WriteMeshToDisk2`): node coordinates as `ES17.10`,
//! element and boundary connectivity, tags and parents, line by line.
//!
//! `tests/data/<case>` holds the input mesh, `expected/` and a `case.sif`
//! reduced to the keywords the mesh step reads. The cases are the quick
//! tests `DisContBoundary{,3D,Double,DoubleNames,FullAngle}`.

#![allow(dead_code)]

mod support;

use std::path::{Path, PathBuf};

use elmer_mesh_ops::create_discont_mesh;

fn data(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
}

fn fixture(case: &str) {
    let dir = data(case);
    let support::Case {
        simulation,
        bodies,
        bcs,
        mut mesh,
    } = support::load(&dir, "case.sif");
    support::map_body_tags(&mut mesh, &bodies);
    support::map_bc_tags(&mut mesh, &bcs);
    let out = create_discont_mesh(&mesh, &simulation, &bcs)
        .unwrap()
        .expect("a discontinuous mesh");
    support::compare(&out.mesh, &bcs, &dir, &dir.join("expected"));
}

#[test]
fn interface_in_a_square() {
    fixture("DisContBoundary");
}

#[test]
fn interface_in_a_box() {
    fixture("DisContBoundary3D");
}

#[test]
fn doubled_boundary() {
    fixture("DisContBoundaryDouble");
}

#[test]
fn doubled_boundary_by_mesh_names() {
    fixture("DisContBoundaryDoubleNames");
}

#[test]
fn full_angle_tetrahedra() {
    fixture("DisContBoundaryFullAngle");
}
