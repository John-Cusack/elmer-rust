//! S1 quick-test targets (strategy report 5.3) through `run_case`: hand-built
//! case, ElmerGrid mesh DB, `elmer-solver`'s default backends (BiCGStab with
//! ILU, `Linear System Scaling`), `ComputeNorm` against the pristine
//! `Reference Norm` with Elmer's default 1e-5 relative tolerance.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_run::{run_case, Case, Outcome};
use elmer_solver::Model;
use support::case::{self as cases, fe, scratch_case};
use support::logger;
use support::mesh::FileMesh;

/// `Reference Norm Tolerance` default.
const TOLERANCE: f64 = 1.0e-5;

fn run(case: &str, mesh: &str, model: Model, scale: f64) -> (Outcome, std::path::PathBuf) {
    let _ = elmer_physics_heat::install_fe(fe());
    let dir = scratch_case(case, case);
    let mesh_dir = dir.join(mesh);
    let m = FileMesh::read(&mesh_dir, 2).scaled(scale);
    let (log, buf) = logger();
    let out = run_case(Case { model, mesh: &m, case_dir: dir.clone(), mesh_dir }, &log, &SerialComm)
        .unwrap_or_else(|e| panic!("{case}: {e}\n{}", buf.text()));
    (out, dir)
}

fn check_norm(case: &str, out: &Outcome, reference: f64, dir: &Path) {
    let [(1, norm)] = out.norms[..] else {
        panic!("{case}: {:?}", out.norms)
    };
    let rel = (norm - reference).abs() / reference;
    println!("{case}: native norm {norm:.17e}, reference {reference:e}, rel {rel:e}");
    assert!(rel <= TOLERANCE, "{case}: norm {norm}, rel {rel:e}");
    assert_eq!(out.passed, Some(true));
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
}

#[test]
fn thermal_pipe_2d() {
    let (out, dir) = run("thermal_pipe_2d", "OnePipe", cases::thermal_pipe_2d(), 1.0);
    check_norm("ThermalPipe2D", &out, 5.44659991, &dir);
    // Post File with BCs on several tags: GeometryIds from the mapped tags.
    assert_eq!(out.written, vec![dir.join("OnePipe/case_t0001.vtu")]);
}

#[test]
fn coordinate_scaling() {
    let (out, dir) = run("coordinate_scaling", "square", cases::coordinate_scaling(), 0.001);
    check_norm("CoordinateScaling", &out, 3.937790364340947e-2, &dir);
}

#[test]
fn wildcardsize() {
    let (out, dir) = run("wildcardsize", "square", cases::wildcardsize(), 1.0);
    check_norm("wildcardsize", &out, 3.50844423e-2, &dir);
}
