//! The strangler `exec`s the Fortran solver named by `ELMER_FORTRAN_SOLVER`.
//! A fake solver script stands in for it.

use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("elmer-cli-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Fake solver: prints argv, cwd, a marker env var and the guard, exits 7.
fn fake_solver(dir: &Path) -> PathBuf {
    let p = dir.join("fake_elmersolver");
    std::fs::write(
        &p,
        "#!/bin/sh\necho \"fake args: $*\"\necho \"cwd: $(pwd)\"\n\
         echo \"marker: $MARKER\"\necho \"guard: $ELMER_CLI_DELEGATED\"\nexit 7\n",
    )
    .unwrap();
    std::fs::set_permissions(&p, std::fs::Permissions::from_mode(0o755)).unwrap();
    p
}

fn run(exe: &str, dir: &Path, args: &[&str], envs: &[(&str, &Path)], plain: &[(&str, &str)]) -> Output {
    let mut c = Command::new(exe);
    c.args(args)
        .current_dir(dir)
        .env_remove("ELMER_FORTRAN_SOLVER")
        .env_remove("ELMER_FORTRAN_SOLVER_MPI")
        .env_remove("ELMER_NATIVE")
        .env_remove("ELMER_SHADOW")
        .env_remove("ELMER_CLI_DELEGATED")
        .env_remove("OMPI_COMM_WORLD_SIZE");
    for (k, v) in envs {
        c.env(k, v);
    }
    for (k, v) in plain {
        c.env(k, v);
    }
    c.output().unwrap()
}

fn text(o: &[u8]) -> String {
    String::from_utf8_lossy(o).into_owned()
}

#[test]
fn fortran_first_execs_with_args_cwd_env_and_status() {
    let d = scratch("first");
    let fake = fake_solver(&d);
    let exe = env!("CARGO_BIN_EXE_ElmerSolver");
    let o = run(exe, &d, &["case.sif", "-s"], &[("ELMER_FORTRAN_SOLVER", &fake)], &[("MARKER", "kept")]);
    assert_eq!(o.status.code(), Some(7));
    let out = text(&o.stdout);
    assert_eq!(
        out,
        format!(
            "fake args: case.sif -s\ncwd: {}\nmarker: kept\nguard: 1\n",
            d.canonicalize().unwrap().display()
        )
    );
    assert_eq!(text(&o.stderr), "");
}

#[test]
fn native_refusal_falls_back_without_printing_first() {
    let d = scratch("refuse");
    let fake = fake_solver(&d);
    std::fs::write(d.join("case.sif"), "Header\nEnd\n").unwrap();
    let exe = env!("CARGO_BIN_EXE_ElmerSolver");
    let o = run(exe, &d, &["case.sif"], &[("ELMER_FORTRAN_SOLVER", &fake)], &[("ELMER_NATIVE", "1")]);
    assert_eq!(o.status.code(), Some(7));
    // Nothing from the native side precedes the Fortran output.
    assert!(text(&o.stdout).starts_with("fake args: case.sif\n"));
}

#[test]
fn startinfo_run_with_no_arguments_delegates_unchanged() {
    let d = scratch("startinfo");
    let fake = fake_solver(&d);
    std::fs::write(d.join("ELMERSOLVER_STARTINFO"), "case.sif\n").unwrap();
    std::fs::write(d.join("case.sif"), "Header\nEnd\n").unwrap();
    let o = run(env!("CARGO_BIN_EXE_ElmerSolver"), &d, &[], &[("ELMER_FORTRAN_SOLVER", &fake)], &[("ELMER_NATIVE", "1")]);
    assert_eq!(o.status.code(), Some(7));
    assert!(text(&o.stdout).starts_with("fake args: \n"));
}

#[test]
fn refusal_without_fortran_exits_78_with_reasons_on_stderr() {
    let d = scratch("nofortran");
    std::fs::write(d.join("case.sif"), "Header\nEnd\n").unwrap();
    let o = run(env!("CARGO_BIN_EXE_ElmerSolver"), &d, &["case.sif"], &[], &[]);
    assert_eq!(o.status.code(), Some(78));
    assert_eq!(text(&o.stdout), "");
    let err = text(&o.stderr);
    assert!(err.contains("refused: native engine not available"), "{err}");
    assert!(err.contains("ELMER_FORTRAN_SOLVER"));
}

#[test]
fn missing_input_file_and_startinfo_are_input_errors() {
    let d = scratch("noinput");
    let exe = env!("CARGO_BIN_EXE_ElmerSolver");
    let o = run(exe, &d, &["nope.sif"], &[], &[]);
    assert_eq!(o.status.code(), Some(2));
    assert_eq!(text(&o.stdout), "ERROR:: MAIN: Unable to find input file [nope.sif], can not execute.\n");
    let o = run(exe, &d, &[], &[], &[]);
    assert_eq!(o.status.code(), Some(2));
    assert_eq!(text(&o.stdout), "ERROR:: MAIN: Unable to find ELMERSOLVER_STARTINFO, can not execute.\n");
}

#[test]
fn recursion_into_self_is_refused() {
    let d = scratch("recursion");
    let exe = env!("CARGO_BIN_EXE_ElmerSolver");
    let o = run(exe, &d, &["case.sif"], &[("ELMER_FORTRAN_SOLVER", Path::new(exe))], &[]);
    assert_eq!(o.status.code(), Some(2));
    assert!(text(&o.stderr).contains("would re-enter elmer-cli"));
}

#[test]
fn missing_fortran_path_is_an_input_error() {
    let d = scratch("missing");
    let o = run(env!("CARGO_BIN_EXE_ElmerSolver"), &d, &["case.sif"], &[("ELMER_FORTRAN_SOLVER", Path::new("/nonexistent/x"))], &[]);
    assert_eq!(o.status.code(), Some(2));
    assert!(text(&o.stderr).contains("not an executable file"));
}

#[test]
fn mpi_binary_prefers_its_own_path_and_multi_rank_skips_native() {
    let d = scratch("mpi");
    let fake = fake_solver(&d);
    std::fs::write(d.join("case.sif"), "Header\nEnd\n").unwrap();
    let exe = env!("CARGO_BIN_EXE_ElmerSolver_mpi");
    let o = run(exe, &d, &["case.sif"], &[("ELMER_FORTRAN_SOLVER_MPI", &fake)], &[("OMPI_COMM_WORLD_SIZE", "4"), ("ELMER_NATIVE", "1")]);
    assert_eq!(o.status.code(), Some(7));
    // Multi-rank without a Fortran solver: collective refusal, no native run.
    let o = run(exe, &d, &["case.sif"], &[], &[("OMPI_COMM_WORLD_SIZE", "4"), ("ELMER_NATIVE", "1")]);
    assert_eq!(o.status.code(), Some(78));
    assert!(text(&o.stderr).contains("multi-rank"));
}

#[test]
fn shadow_mode_runs_fortran_as_child_and_keeps_its_status() {
    let d = scratch("shadow");
    let fake = fake_solver(&d);
    let o = run(env!("CARGO_BIN_EXE_ElmerSolver"), &d, &["case.sif"], &[("ELMER_FORTRAN_SOLVER", &fake)], &[("ELMER_SHADOW", "1")]);
    assert_eq!(o.status.code(), Some(7));
    assert!(text(&o.stdout).starts_with("fake args: case.sif\n"));
    assert!(text(&o.stderr).contains("shadow mode"));
}
