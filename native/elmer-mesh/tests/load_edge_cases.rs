//! ASCII mesh DB corners Elmer documents but the 482 quick tests never
//! trigger: halo (`/ partn`) prefixes, permuted node tags, permuted bulk
//! tags with parent remap, and right-only boundary parents.
//!
//! All behaviors mirror `MeshIO.F90` (`ReadNodesFile`, `ReadElementsFile`,
//! `ReadBoundaryFile`, including the right-to-left swap).

use std::fs;
use std::path::PathBuf;

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "elmer_mesh_edge_{}_{}_{}",
        std::process::id(),
        name,
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&d).unwrap();
    d
}

fn write(d: &std::path::Path, header: &str, nodes: &str, elements: &str, boundary: &str) {
    fs::write(d.join("mesh.header"), header).unwrap();
    fs::write(d.join("mesh.nodes"), nodes).unwrap();
    fs::write(d.join("mesh.elements"), elements).unwrap();
    fs::write(d.join("mesh.boundary"), boundary).unwrap();
}

const HEADER: &str = "4 2 1\n2\n303 2\n202 1\n";

#[test]
fn halo_prefix_elements_and_boundary() {
    let d = dir("halo");
    write(
        &d,
        HEADER,
        "1 - 0 0 0\n2 - 1 0 0\n3 - 1 1 0\n4 - 0 1 0\n",
        // `/ partn` prefix: Elmer skips it via the halo offset.
        "1 /0 1 303 1 2 3\n2 /1 1 303 1 3 4\n",
        "1 /0 7 1 0 202 1 2\n",
    );
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    assert_eq!(m.bulk.len(), 2);
    assert_eq!(m.bulk[0].body, 1);
    assert_eq!(m.bulk[0].code, 303);
    assert_eq!(m.bulk[0].nodes, vec![0, 1, 2]);
    assert_eq!(m.boundary[0].constraint, 7);
    assert_eq!(m.boundary[0].parents, [Some(0), None]);
    fs::remove_dir_all(&d).ok();
}

#[test]
fn permuted_node_tags_reorder_by_tag() {
    let d = dir("nodetag");
    write(
        &d,
        HEADER,
        "3 - 1 1 0\n1 - 0 0 0\n4 - 0 1 0\n2 - 1 0 0\n",
        "1 1 303 1 2 3\n2 1 303 1 3 4\n",
        "1 7 1 0 202 1 2\n",
    );
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    // Node with tag t sits at position t-1 (Elmer's PermuteNodeNumbering).
    assert_eq!(m.nodes[0], [0.0, 0.0, 0.0]);
    assert_eq!(m.nodes[1], [1.0, 0.0, 0.0]);
    assert_eq!(m.nodes[2], [1.0, 1.0, 0.0]);
    assert_eq!(m.nodes[3], [0.0, 1.0, 0.0]);
    fs::remove_dir_all(&d).ok();
}

#[test]
fn permuted_bulk_tags_remap_parents() {
    let d = dir("bulktag");
    write(
        &d,
        HEADER,
        "1 - 0 0 0\n2 - 1 0 0\n3 - 1 1 0\n4 - 0 1 0\n",
        // File order tag 2 then tag 1; bulk elements stay in file order.
        "2 1 303 1 2 3\n1 1 303 1 3 4\n",
        // Parent tag 1 lives at file position 1 (Elmer's LocalEPerm).
        "1 7 1 0 202 1 2\n",
    );
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    assert_eq!(m.bulk[0].tag, 2);
    assert_eq!(m.bulk[1].tag, 1);
    assert_eq!(m.boundary[0].parents, [Some(1), None]);
    fs::remove_dir_all(&d).ok();
}

#[test]
fn right_only_parent_becomes_left() {
    let d = dir("swap");
    write(
        &d,
        HEADER,
        "1 - 0 0 0\n2 - 1 0 0\n3 - 1 1 0\n4 - 0 1 0\n",
        "1 1 303 1 2 3\n2 1 303 1 3 4\n",
        // Left 0, right 1: Elmer always wants Left present.
        "1 7 0 1 202 1 2\n",
    );
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    assert_eq!(m.boundary[0].parents, [Some(0), None]);
    fs::remove_dir_all(&d).ok();
}
