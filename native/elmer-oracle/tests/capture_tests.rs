//! End-to-end `capture` through stub binaries: a fake ElmerGrid and a
//! fake ElmerSolver honouring the runner contract (exit 0, write
//! `TEST.PASSED`, print the verdict lines). This exercises plan, stage,
//! execute, verdict-read, and manifest-freeze without a Fortran build.

use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

fn unique_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "elmer-oracle-{}-{}-{}",
        tag,
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_exe(path: &Path, body: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(body.as_bytes()).unwrap();
    std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).unwrap();
}

fn fixture_tree(root: &Path) {
    let tests = root.join("tests");
    let build = root.join("build");
    let t = tests.join("faketest");
    std::fs::create_dir_all(&t).unwrap();
    std::fs::create_dir_all(&build).unwrap();
    std::fs::write(
        t.join("CMakeLists.txt"),
        "INCLUDE(test_macros)\nADD_ELMER_TEST(faketest LABELS quick)\n",
    )
    .unwrap();
    std::fs::write(
        t.join("runtest.cmake"),
        "include(test_macros)\nexecute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh.grd)\nRUN_ELMER_TEST()\n",
    )
    .unwrap();
    std::fs::write(t.join("ELMERSOLVER_STARTINFO"), "case.sif\n").unwrap();
    std::fs::write(
        t.join("case.sif"),
        "Header\n  Mesh DB = \".\" \"m\"\nEnd\n\nSimulation\n  Simulation Type = \"Steady State\"\nEnd\n\nSolver 1\n  Equation = \"Heat Equation\"\nEnd\n",
    )
    .unwrap();
    write_exe(&build.join("bin-grid"), "#!/bin/sh\nexit 0\n");
    let femsrc = build.join("fem").join("src");
    std::fs::create_dir_all(&femsrc).unwrap();
    // The CLI derives tool paths as <build>/fem/src/<Tool>.
    std::fs::copy(build.join("bin-grid"), femsrc.join("ElmerGrid")).unwrap();
    write_exe(
        &femsrc.join("ElmerSolver"),
        "#!/bin/sh\necho 'CompareToReferenceSolution: FakeSolver Norm 1.0 vs 1.0'\necho 'SOLVER TOTAL TIME(CPU,REAL):    1.000    2.000'\necho 1 > TEST.PASSED\n",
    );
    for name in ["ViewFactors", "Radiators"] {
        std::fs::copy(build.join("bin-grid"), femsrc.join(name)).unwrap();
    }
    std::fs::create_dir_all(build.join("meshgen2d").join("src")).unwrap();
}

#[test]
fn capture_freezes_stub_run_as_golden() {
    let root = unique_dir("capture");
    fixture_tree(&root);
    let bin = PathBuf::from(env!("CARGO_BIN_EXE_elmer-oracle"));
    let golden = root.join("golden");
    let manifest = root.join("golden.manifest");
    let out = Command::new(&bin)
        .args([
            "capture",
            &root.join("build").display().to_string(),
            &root.join("tests").display().to_string(),
            "faketest",
            &golden.display().to_string(),
            &manifest.display().to_string(),
            "--o2",
            "--o2m",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "capture failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("passed: true"), "stdout: {stdout}");
    assert!(golden.join("run").join("TEST.PASSED").exists());
    // Injection lands in the staged sif the solver actually read.
    let staged = std::fs::read_to_string(golden.join("run").join("case.sif")).unwrap();
    assert!(staged.contains("Output File = \"case.result\""));
    assert!(staged.contains("Linear System Save = Logical True"));
    let text = std::fs::read_to_string(&manifest).unwrap();
    assert!(text.contains("faketest/TEST.PASSED"));
    assert!(text.contains("faketest/test-stdout.log"));
}
