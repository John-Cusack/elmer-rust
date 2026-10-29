//! `mpirun -np N elmer-mpi-run` on the AngleMetisLayer fixture
//! (`tests/data/angle`) at np = 1..=4: every run must reproduce the pristine
//! serial ElmerSolver norm to 1e-12 relative (`TEST.PASSED[_N]` = 1), also
//! with `Mesh Levels` 2 and 3.
//! Skipped with a message when `mpirun` is not installed. Always
//! `HWLOC_COMPONENTS=-gl`, never more than 4 ranks.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const RUN: &str = env!("CARGO_BIN_EXE_elmer-mpi-run");

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dst = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &dst);
        } else {
            std::fs::copy(e.path(), dst).unwrap();
        }
    }
}

/// A scratch copy of the fixture, so `TEST.PASSED` lands outside the tree.
fn stage(name: &str) -> PathBuf {
    let src = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/angle");
    let dir = std::env::temp_dir().join(format!("elmer-mpi-run-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    copy_dir(&src, &dir);
    dir
}

fn mpirun(np: usize, dir: &Path, sif: &str) -> Option<Output> {
    let mut cmd = Command::new("mpirun");
    cmd.args(["--oversubscribe", "-np", &np.to_string(), RUN, sif])
        .current_dir(dir)
        .env("HWLOC_COMPONENTS", "-gl")
        .env("OMP_NUM_THREADS", "1");
    match cmd.output() {
        Ok(o) => Some(o),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("mpirun not found; skipping np={np}");
            None
        }
        Err(e) => panic!("mpirun: {e}"),
    }
}

fn passes_at(np: usize) {
    let dir = stage(&format!("np{np}"));
    let Some(out) = mpirun(np, &dir, "case.sif") else { return };
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "np={np}:\n{stdout}\n{}", String::from_utf8_lossy(&out.stderr));
    let name = if np > 1 { format!("TEST.PASSED_{np}") } else { "TEST.PASSED".into() };
    let verdict = std::fs::read_to_string(dir.join(&name)).unwrap_or_else(|e| panic!("{name}: {e}\n{stdout}"));
    assert_eq!(verdict.trim(), "1", "np={np}:\n{stdout}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn serial_mesh_at_np1() {
    passes_at(1);
}

#[test]
fn two_partitions() {
    passes_at(2);
}

#[test]
fn three_partitions() {
    passes_at(3);
}

#[test]
fn four_partitions() {
    passes_at(4);
}

/// `case.sif` with `Mesh Levels = levels` and the pristine serial norm of
/// that refined mesh (from its full-precision `Output File`).
fn with_levels(dir: &Path, levels: u32, norm: &str) -> String {
    let sif = std::fs::read_to_string(dir.join("case.sif")).unwrap();
    let sif = sif.replace("Simulation\n", &format!("Simulation\n  Mesh Levels = {levels}\n"));
    let sif = sif.replace("Real 0.16133448050272817", &format!("Real {norm}"));
    let name = format!("levels{levels}.sif");
    std::fs::write(dir.join(&name), sif).unwrap();
    name
}

/// Every rank splits its piece; the new edge and centre nodes get global
/// ids and holders across ranks. Two levels split the first level's new
/// nodes again. Each np must match the serial pristine norm, which the
/// consistent norm only does when every shared new node is counted once.
fn levels_pass_at(levels: u32, norm: &str) {
    for np in 1..=4 {
        let dir = stage(&format!("levels{levels}-np{np}"));
        let sif = with_levels(&dir, levels, norm);
        let Some(out) = mpirun(np, &dir, &sif) else { return };
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(out.status.success(), "np={np}:\n{stdout}\n{}", String::from_utf8_lossy(&out.stderr));
        let name = if np > 1 { format!("TEST.PASSED_{np}") } else { "TEST.PASSED".into() };
        let verdict = std::fs::read_to_string(dir.join(&name)).unwrap_or_else(|e| panic!("{name}: {e}\n{stdout}"));
        assert_eq!(verdict.trim(), "1", "Mesh Levels = {levels}, np={np}:\n{stdout}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

#[test]
fn mesh_levels_2() {
    levels_pass_at(2, "0.1639056784545602");
}

#[test]
fn mesh_levels_3() {
    levels_pass_at(3, "0.16525319400846628");
}

/// What the parallel path cannot run is refused (exit 78) on every rank
/// before anything is written.
#[test]
fn direct_solver_is_refused() {
    let dir = stage("direct");
    let sif = std::fs::read_to_string(dir.join("case.sif")).unwrap();
    std::fs::write(dir.join("direct.sif"), sif.replace("Linear System Solver = Iterative", "Linear System Solver = Direct")).unwrap();
    let Some(out) = mpirun(2, &dir, "direct.sif") else { return };
    assert_eq!(out.status.code(), Some(78), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stderr).contains("Linear System Solver = direct in parallel"));
    assert!(!dir.join("TEST.PASSED_2").exists());
    std::fs::remove_dir_all(&dir).unwrap();
}
