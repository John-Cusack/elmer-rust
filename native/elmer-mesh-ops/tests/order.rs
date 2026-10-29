//! `Increase Element Order` against the pristine ElmerSolver's own
//! `SaveMesh` output (`WriteMeshToDisk2`): node coordinates as `ES17.10`,
//! element and boundary connectivity, tags and parents, line by line.
//!
//! `tests/data/<case>` holds the input mesh, `expected/` and a `case.sif`
//! reduced to the keywords the mesh step reads. The cases are the quick
//! tests `IncreaseOrder{Tri,Quad,Tet,Prism,Brick}`.

#![allow(dead_code)]

mod support;

use std::path::{Path, PathBuf};

use elmer_core::{Params, Value};
use elmer_mesh_ops::increase_element_order;

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
    assert_eq!(
        simulation.get("increase element order"),
        Some(&Value::Logical(true))
    );
    support::map_body_tags(&mut mesh, &bodies);
    support::map_bc_tags(&mut mesh, &bcs);
    let out = increase_element_order(&mesh).unwrap();
    support::compare(&out, &bcs, &dir, &dir.join("expected"));
}

#[test]
fn triangles_to_306() {
    fixture("IncreaseOrderTri");
}

#[test]
fn quadrilaterals_to_408() {
    fixture("IncreaseOrderQuad");
}

#[test]
fn tetrahedra_to_510() {
    fixture("IncreaseOrderTet");
}

#[test]
fn wedges_to_715() {
    fixture("IncreaseOrderPrism");
}

#[test]
fn hexahedra_to_820() {
    fixture("IncreaseOrderBrick");
}
