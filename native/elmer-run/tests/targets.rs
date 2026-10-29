//! S1 and heat-breadth quick tests end to end from their pristine `.sif`:
//! `elmer-sif` (parse, completion, `Post File` auto-solver) ->
//! `elmer-solver`'s default backends -> `ComputeNorm` against the case's own
//! `Reference Norm` with Elmer's default 1e-5 relative tolerance. Meshes are
//! the ElmerGrid mesh DBs the tests generate, read by the test-local reader.

mod support;

use std::path::{Path, PathBuf};


use elmer_core::SerialComm;
use elmer_run::{run_case, Case, Frontend, Outcome};
use support::case::{scratch_case, sif_frontend};
use support::logger;

/// `Reference Norm Tolerance` default.
const TOLERANCE: f64 = 1.0e-5;

/// Runs `sif` from a scratch copy of `data/<case>`.
fn run(case: &str, sif: &str, tag: &str) -> (Outcome, PathBuf) {
    let dir = scratch_case(case, tag);
    let frontend = sif_frontend();
    let loaded = frontend.load(&dir.join(sif)).unwrap_or_else(|e| panic!("{tag}: {e}"));
    let (log, buf) = logger();
    let case = Case {
        model: loaded.model,
        reloads: loaded.reloads,
        mesh: loaded.mesh.as_ref(),
        case_dir: loaded.case_dir,
        mesh_dir: loaded.mesh_dir,
    };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{tag}: {e}\n{}", buf.text()));
    (out, dir)
}

fn check_norm(tag: &str, out: &Outcome, reference: f64, dir: &Path) {
    let Some(&(_, norm)) = out.norms.iter().find(|(id, _)| *id == 1) else {
        panic!("{tag}: {:?}", out.norms)
    };
    let rel = (norm - reference).abs() / reference;
    println!("{tag}: native norm {norm:.17e}, reference {reference:e}, rel {rel:e}");
    assert!(rel <= TOLERANCE, "{tag}: norm {norm}, rel {rel:e}");
    assert_eq!(out.passed, Some(true));
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
}

#[test]
fn minimal_example() {
    let (out, dir) = run("minimal_example", "case.sif", "sif-minimal");
    check_norm("MinimalExample", &out, 2.88675135e-1, &dir);
    assert_eq!(out.written, vec![dir.join("triangle/case_t0001.vtu")]);
}

#[test]
fn thermal_pipe_2d() {
    let (out, dir) = run("thermal_pipe_2d", "case.sif", "sif-pipe");
    check_norm("ThermalPipe2D", &out, 5.44659991, &dir);
    // Post File with BCs on several tags: GeometryIds from the mapped tags.
    assert_eq!(out.written, vec![dir.join("OnePipe/case_t0001.vtu")]);
}

#[test]
fn coordinate_scaling() {
    let (out, dir) = run("coordinate_scaling", "case.sif", "sif-scaling");
    check_norm("CoordinateScaling", &out, 3.937790364340947e-2, &dir);
}

#[test]
fn wildcardsize() {
    let (out, dir) = run("wildcardsize", "case.sif", "sif-wildcard");
    check_norm("wildcardsize", &out, 3.50844423e-2, &dir);
}

#[test]
fn heateq_bdf2() {
    let (out, dir) = run("heateq_bdf", "TempDist_bdf2.sif", "sif-bdf2");
    check_norm("heateq_bdf2", &out, 5.214225320136149e-2, &dir);
    assert_eq!(out.report.steps, 8);
}

#[test]
fn heateq_bdf3() {
    let (out, dir) = run("heateq_bdf", "TempDist_bdf3.sif", "sif-bdf3");
    check_norm("heateq_bdf3", &out, 5.193298434405153e-2, &dir);
    assert_eq!(out.report.steps, 8);
}
