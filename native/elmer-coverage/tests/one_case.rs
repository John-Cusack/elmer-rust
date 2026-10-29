//! `run_one` on two quick-test cases that `elmer-run` already matches, and on
//! one it must refuse. Its own test binary: `run_one` changes the working
//! directory and parses with the process-wide MATC library.

use std::path::{Path, PathBuf};

use elmer_coverage::child::{run_one, Verdict};

fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for e in std::fs::read_dir(src).unwrap() {
        let e = e.unwrap();
        let to = dst.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &to);
        } else {
            std::fs::copy(e.path(), to).unwrap();
        }
    }
}

fn scratch(tag: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("coverage-{tag}"));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

#[test]
fn matches_minimal_example_and_refuses_unknown_solver() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../elmer-run/tests/data/minimal_example");
    let dir = scratch("minimal");
    copy_dir(&data, &dir);

    let report = run_one(&dir, "case.sif", &[], false);
    assert_eq!(report.verdict, Verdict::Match, "{report:?}");
    assert_eq!(report.libs, ["heatsolve"]);
    assert!(report.blockers.is_empty(), "{:?}", report.blockers);
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");

    // The same case with its solver replaced is refused before any output.
    let sif = std::fs::read_to_string(dir.join("case.sif")).unwrap();
    std::fs::write(dir.join("other.sif"), sif.replace("HeatSolve", "NoSuchSolve").replace("HeatSolver", "NoSuchSolver")).unwrap();
    let report = run_one(&dir, "other.sif", &[], false);
    assert_eq!(report.verdict, Verdict::Refused, "{report:?}");
    assert_eq!(report.blockers, ["solver: nosuchsolve nosuchsolver"]);
}

#[test]
fn child_report_round_trips() {
    use elmer_coverage::child::ChildReport;
    let r = ChildReport {
        verdict: Verdict::Mismatch,
        detail: "Solver 1 FAILED".into(),
        blockers: vec!["feature: eigenvalue analysis".into()],
        libs: vec!["heatsolve".into()],
    };
    assert_eq!(ChildReport::parse(&r.render()), Some(r));
}
