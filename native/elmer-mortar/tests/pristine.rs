//! Periodic projectors checked byte for byte against the pristine
//! `ElmerSolver`'s `Save Projector` files.
//!
//! Each `fixtures/<case>` is an Elmer quick test (`fem/tests/<case>`) whose
//! `test.sif` got `Save Projector = Logical True` next to every `Periodic BC`
//! (and an `Output File`); `mesh/` is its ElmerGrid mesh and `expected/` the
//! `p<bc>.dat` / `p<bc>_rsum.dat` files the serial pristine run wrote. The
//! cases cover the automatic overlay (bounding-box scale and translation),
//! explicit `Periodic BC Translate` / `Scale`, two periodic pairs sharing
//! corner nodes, and a nonconforming interface where target nodes fall inside
//! elements (non-trivial weights from `GlobalToLocal`).

mod support;

use std::path::{Path, PathBuf};

use elmer_core::{MeshView, NullLogger, Params};
use elmer_mortar::{generate_periodic_projectors, ProjectorEnv};
use support::{Case, FileMesh};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// Builds every periodic projector of the case, saves them into a scratch
/// directory and compares each expected file byte for byte.
fn check(name: &str) {
    let dir = fixture(name);
    let case = Case::read(&dir.join("case.sif"));
    let mesh = FileMesh::read(&dir.join("mesh"), 3);
    let mesh = FileMesh::read(&dir.join("mesh"), case.coord_dim(mesh.dim()));
    let bcs: Vec<&dyn Params> = case.bcs.iter().map(|b| b as &dyn Params).collect();
    let bc_of = |t: i32| case.bc_of(t);
    let env = ProjectorEnv {
        mesh: &mesh,
        bcs: &bcs,
        bc_of: &bc_of,
        simulation: &case.simulation,
        solver: case.solvers.first().map(|s| s as &dyn Params),
        coord_dim: mesh.dim(),
        log: &NullLogger,
    };
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("elmer-mortar")
        .join(name);
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let (projectors, stop) = generate_periodic_projectors(&env, &out, 0).unwrap();
    assert!(!stop);
    assert!(
        projectors.iter().any(Option::is_some),
        "{name}: no projector built"
    );

    let mut expected: Vec<PathBuf> = std::fs::read_dir(dir.join("expected"))
        .unwrap()
        .map(|e| e.unwrap().path())
        // The Save Projector files; `case.result` is for `tests/e2e.rs`.
        .filter(|p| p.extension().is_some_and(|x| x == "dat"))
        .collect();
    expected.sort();
    for want in &expected {
        let file = want.file_name().unwrap();
        let got = out.join(file);
        compare(name, want, &got);
    }
    let written = std::fs::read_dir(&out).unwrap().count();
    assert_eq!(written, expected.len(), "{name}: file count");
}

fn compare(name: &str, want: &Path, got: &Path) {
    let w = std::fs::read_to_string(want).unwrap();
    let g = std::fs::read_to_string(got).unwrap_or_else(|e| panic!("{}: {e}", got.display()));
    for (i, (a, b)) in w.lines().zip(g.lines()).enumerate() {
        assert_eq!(
            a,
            b,
            "{name} {}: line {}",
            want.file_name().unwrap().to_string_lossy(),
            i + 1
        );
    }
    assert_eq!(
        w.lines().count(),
        g.lines().count(),
        "{name}: line count of {}",
        want.display()
    );
    assert_eq!(w, g);
}

#[test]
fn periodic1_two_pairs() {
    check("periodic1");
}

#[test]
fn periodic2_default_overlay() {
    check("periodic2");
}

#[test]
fn periodic_offset() {
    check("periodic_offset");
}

#[test]
fn periodic_nonconforming_weights() {
    check("periodic_nonconforming");
}
