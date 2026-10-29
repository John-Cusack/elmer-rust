//! Binary mesh inputs (Elmer's ASCII-then-`.bin`-then-`.sbin` fallback).
//!
//! Layouts straight from `MeshIO.F90`: stream access, so records are packed
//! with no markers - nodes are `(tag:i32, xyz:3xf64|3xf32)`, bulk elements
//! `(tag, partn, body, code):i32` plus 1-based node tags, boundary elements
//! `(tag, partn, bndry, left, right, code):i32` plus 1-based node tags.

use std::fs;
use std::path::PathBuf;

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "elmer_mesh_bin_{}_{}_{}",
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

fn put_i32(v: &mut Vec<u8>, x: i32) {
    v.extend_from_slice(&x.to_ne_bytes());
}

fn put_f64(v: &mut Vec<u8>, x: f64) {
    v.extend_from_slice(&x.to_ne_bytes());
}

fn put_f32(v: &mut Vec<u8>, x: f32) {
    v.extend_from_slice(&x.to_ne_bytes());
}

const HEADER: &str = "4 2 1\n2\n303 2\n202 1\n";
// ASCII content: same topology, all x coords shifted by +10 so the source
// file is distinguishable from the binary one.
const NODES_A: &str = "1 - 10 0 0\n2 - 11 0 0\n3 - 11 1 0\n4 - 10 1 0\n";
const ELEMENTS_A: &str = "1 1 303 1 2 3\n2 1 303 1 3 4\n";
const BOUNDARY_A: &str = "1 7 1 0 202 1 2\n";

fn nodes_bin(single: bool) -> Vec<u8> {
    let mut v = Vec::new();
    for (t, c) in [
        (1, [0.0, 0.0, 0.0]),
        (2, [1.0, 0.0, 0.0]),
        (3, [1.0, 1.0, 0.0]),
        (4, [0.0, 1.0, 0.0]),
    ] {
        put_i32(&mut v, t);
        for x in c {
            if single {
                put_f32(&mut v, x as f32);
            } else {
                put_f64(&mut v, x);
            }
        }
    }
    v
}

fn elements_bin() -> Vec<u8> {
    let mut v = Vec::new();
    for (tag, nodes) in [(1, [1, 2, 3]), (2, [1, 3, 4])] {
        put_i32(&mut v, tag);
        put_i32(&mut v, 0);
        put_i32(&mut v, 1);
        put_i32(&mut v, 303);
        for n in nodes {
            put_i32(&mut v, n);
        }
    }
    v
}

fn boundary_bin() -> Vec<u8> {
    let mut v = Vec::new();
    put_i32(&mut v, 1);
    put_i32(&mut v, 0);
    put_i32(&mut v, 7);
    put_i32(&mut v, 1);
    put_i32(&mut v, 0);
    put_i32(&mut v, 202);
    put_i32(&mut v, 1);
    put_i32(&mut v, 2);
    v
}

fn write_ascii(d: &std::path::Path) {
    fs::write(d.join("mesh.header"), HEADER).unwrap();
    fs::write(d.join("mesh.nodes"), NODES_A).unwrap();
    fs::write(d.join("mesh.elements"), ELEMENTS_A).unwrap();
    fs::write(d.join("mesh.boundary"), BOUNDARY_A).unwrap();
}

fn write_bin(d: &std::path::Path) {
    fs::write(d.join("mesh.header"), HEADER).unwrap();
    fs::write(d.join("mesh.nodes.bin"), nodes_bin(false)).unwrap();
    fs::write(d.join("mesh.elements.bin"), elements_bin()).unwrap();
    fs::write(d.join("mesh.boundary.bin"), boundary_bin()).unwrap();
}

#[test]
fn ascii_wins_over_bin() {
    let d = dir("prefer");
    write_ascii(&d);
    fs::write(d.join("mesh.nodes.bin"), nodes_bin(false)).unwrap();
    fs::write(d.join("mesh.elements.bin"), elements_bin()).unwrap();
    fs::write(d.join("mesh.boundary.bin"), boundary_bin()).unwrap();
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    // ASCII content (x shifted by +10), not the binary content.
    assert_eq!(m.nodes[0], [10.0, 0.0, 0.0]);
    assert_eq!(m.nodes[2], [11.0, 1.0, 0.0]);
    assert_eq!(m.bulk.len(), 2);
    assert_eq!(m.boundary.len(), 1);
    fs::remove_dir_all(&d).ok();
}

#[test]
fn bin_fallback_loads_elmer_layout() {
    let d = dir("fallback");
    write_bin(&d);
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    assert_eq!(m.nodes.len(), 4);
    assert_eq!(m.nodes[1], [1.0, 0.0, 0.0]);
    assert_eq!(m.bulk.len(), 2);
    assert_eq!(m.bulk[0].code, 303);
    assert_eq!(m.bulk[0].body, 1);
    assert_eq!(m.bulk[0].nodes, vec![0, 1, 2]);
    assert_eq!(m.bulk[1].nodes, vec![0, 2, 3]);
    assert_eq!(m.boundary.len(), 1);
    let b = &m.boundary[0];
    assert_eq!(b.code, 202);
    assert_eq!(b.constraint, 7);
    assert_eq!(b.parents, [Some(0), None]);
    assert_eq!(b.nodes, vec![0, 1]);
    fs::remove_dir_all(&d).ok();
}

#[test]
fn sbin_fallback_loads_single_precision() {
    let d = dir("single");
    fs::write(d.join("mesh.header"), HEADER).unwrap();
    fs::write(d.join("mesh.nodes.sbin"), nodes_bin(true)).unwrap();
    fs::write(d.join("mesh.elements.bin"), elements_bin()).unwrap();
    fs::write(d.join("mesh.boundary.bin"), boundary_bin()).unwrap();
    let m = elmer_mesh::load(d.to_str().unwrap()).unwrap();
    assert_eq!(m.nodes.len(), 4);
    for (got, want) in m.nodes.iter().zip([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]) {
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() < 1e-6);
        }
    }
    fs::remove_dir_all(&d).ok();
}
