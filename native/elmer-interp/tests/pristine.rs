//! Mesh-to-mesh interpolation checked against the pristine `ElmerSolver`.
//!
//! Each `fixtures/mesh{2,3}d/case.sif` maps `TempA` from the mesh `src` onto
//! the mesh `dst` with `Mesh2MeshSolver` (`InterpolateMeshToMesh`).
//! `pristine_run.log` is the run's log; `expected_tempa.txt` is the
//! full-precision `TempA` it wrote on `dst` (binary VTU parsed to text; target
//! nodes the Fortran does not find keep their initial value `-1000`) and
//! `src_tempa.txt` is `TempA` on `src` dumped by a second `ResultOutputSolver`
//! in the same run. The sif's MATC polynomial is only evaluated to ~5e-12, so
//! the dumped values, not the polynomial, are the input here.
//!
//! The 2-D pair mixes quads and triangles; the 3-D pair is hexahedra and
//! wedges onto tetrahedra. In both, many `dst` nodes lie just outside `src` and
//! are found only through the Fortran's widening tolerance ladder, and the
//! hexahedron/quad Newton iteration stops on a residual threshold, so agreement
//! to ~1e-14 on which nodes are found and on their values exercises the tree
//! search, the tolerance loop and `GlobalToLocal` together.

mod support;

use std::path::PathBuf;

use elmer_core::NullLogger;
use elmer_interp::{build_projector, interpolate, Extras, FieldData, InterpOptions, Locator};
use support::*;

struct Case {
    src: TestMesh,
    dst: TestMesh,
    src_values: Vec<f64>,
    want: Vec<f64>,
    dim: usize,
}

fn load(name: &str, dim: usize) -> Case {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    let column = |f: &str| -> Vec<f64> {
        std::fs::read_to_string(dir.join(f))
            .unwrap()
            .lines()
            .map(|l| l.parse().unwrap())
            .collect()
    };
    let case = Case {
        src: read_elmer_mesh(&dir.join("src"), dim),
        dst: read_elmer_mesh(&dir.join("dst"), dim),
        src_values: column("src_tempa.txt"),
        want: column("expected_tempa.txt"),
        dim,
    };
    assert_eq!(case.src_values.len(), case.src.coords.len());
    assert_eq!(case.want.len(), case.dst.coords.len());
    case
}

fn options(c: &Case) -> InterpOptions {
    InterpOptions {
        coordinate_system_dim: Some(c.dim),
        ..Default::default()
    }
}

fn check_direct(c: &Case, n_unfound: usize) {
    let old = vec![FieldData::nodal("tempa", c.src_values.clone())];
    let mut new = vec![FieldData::nodal("tempa", vec![-1000.0; c.dst.coords.len()])];
    let loc = Locator::new(&c.src, options(c)).unwrap();
    let rep = interpolate(
        &loc,
        &c.dst,
        &old,
        &mut new,
        &Extras::default(),
        &NullLogger,
    )
    .unwrap();
    // "Points not found: N (found M)" in pristine_run.log.
    assert_eq!(rep.n_not_found, n_unfound);
    for (i, &w) in c.want.iter().enumerate() {
        assert_eq!(rep.found[i], w != -1000.0, "node {i} found-state");
        let got = new[0].values[i];
        if w == -1000.0 {
            assert_eq!(got, -1000.0, "unfound node {i} must keep its value");
        } else {
            assert_close(got, w, 1e-14, &format!("node {i}"));
        }
    }
}

fn check_projector(c: &Case) {
    let loc = Locator::new(&c.src, options(c)).unwrap();
    let p = build_projector(&loc, &c.dst, &NodalBasis, &Extras::default(), &NullLogger).unwrap();
    let old = vec![FieldData::nodal("tempa", c.src_values.clone())];
    let mut new = vec![FieldData::nodal("tempa", vec![0.0; c.dst.coords.len()])];
    p.apply(&old, &mut new);
    for (i, &w) in c.want.iter().enumerate() {
        if w == -1000.0 {
            assert!(!p.found[i]);
            assert_eq!(new[0].values[i], 0.0);
        } else {
            assert_close(new[0].values[i], w, 1e-13, &format!("projector node {i}"));
        }
    }
}

#[test]
fn quads_and_triangles_match_pristine_2d() {
    let c = load("mesh2d", 2);
    check_direct(&c, 27);
    check_projector(&c);
}

#[test]
fn hexes_and_wedges_onto_tets_match_pristine_3d() {
    let c = load("mesh3d", 3);
    check_direct(&c, 25);
    check_projector(&c);
}
