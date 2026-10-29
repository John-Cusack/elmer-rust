//! Periodic and mortar projectors checked byte for byte against the
//! pristine `ElmerSolver`'s `Save Projector` files.
//!
//! Each `fixtures/<case>` is an Elmer quick test (`fem/tests/<case>`) whose
//! `.sif` got `Save Projector = Logical True` next to every `Periodic BC` or
//! `Mortar BC` (and an `Output File`); `mesh/` is its ElmerGrid mesh and `expected/` the
//! `p<bc>.dat` / `p<bc>_rsum.dat` files the serial pristine run wrote. The
//! cases cover the automatic overlay (bounding-box scale and translation),
//! explicit `Periodic BC Translate` / `Scale`, two periodic pairs sharing
//! corner nodes, and a nonconforming interface where target nodes fall inside
//! elements (non-trivial weights from `GlobalToLocal`). `mortar_poisson2d`
//! (MortarPoisson2D) has Galerkin projectors between nonmatching meshes,
//! and its `Save Constraint Matrix` output (`cm*.dat`) checks the
//! constraint matrix built from them.

mod support;

use std::path::{Path, PathBuf};

use elmer_core::{ElemId, MeshView, NullLogger, Params};
use elmer_mortar::{
    generate_constraint_matrix, generate_mortar_projectors, generate_periodic_projectors,
    save_constraint_matrix, ConstraintTarget, ProjectorEnv,
};
use support::{Case, FileMesh};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// Which projectors a case builds.
#[derive(Clone, Copy)]
enum Kind {
    /// `GeneratePeriodicProjectors`.
    Periodic,
    /// `GenerateProjectors` for solver 1.
    Mortar,
}

/// Builds every periodic projector of the case, saves them into a scratch
/// directory and compares each expected file byte for byte.
fn check(name: &str) {
    check_kind(name, Kind::Periodic)
}

fn check_kind(name: &str, kind: Kind) {
    let dir = fixture(name);
    let case = Case::read(&dir.join("case.sif"));
    let mesh = FileMesh::read(&dir.join("mesh"), 3);
    let mesh = FileMesh::read(&dir.join("mesh"), case.coord_dim(mesh.bulk_dim()));
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
    let (projectors, stop) = match kind {
        Kind::Periodic => generate_periodic_projectors(&env, &out, 0).unwrap(),
        Kind::Mortar => generate_mortar_projectors(&env, &[], false, &out, 0)
            .unwrap()
            .expect("Apply Mortar BCs"),
    };
    if let Kind::Mortar = kind {
        // Every bulk element is active and `Optimize Bandwidth = False`.
        let perm =
            elmer_assembly::initial_permutation(&mesh, &|e: ElemId| (e.0 as usize) < mesh.n_bulk());
        let solver = &case.solvers[0];
        let t = ConstraintTarget {
            solver,
            name: "temperature",
            dofs: 1,
            perm: &perm.perm,
            n_rows: perm.n_dofs,
        };
        let cm = generate_constraint_matrix(&bcs, &projectors, &t, &NullLogger)
            .unwrap()
            .expect("constraint matrix");
        save_constraint_matrix(&cm, &out.join("cm")).unwrap();
    }
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

#[test]
fn mortar_poisson2d_galerkin() {
    check_kind("mortar_poisson2d", Kind::Mortar);
}
