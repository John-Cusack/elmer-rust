//! MinimalExample end to end: hand-built case and mesh -> `elmer-assembly`
//! -> `elmer-direct` band LU -> `ComputeNorm` -> `TEST.PASSED` and the
//! `Post File` VTU, directly and through `elmer-cli`'s driver.

mod support;

use std::ffi::OsString;

use elmer_cli::driver::{run, Action, Environment};
use elmer_cli::strangler::Variant;
use elmer_cli::ExitCode;
use elmer_core::{SerialComm, Value};
use elmer_run::{run_case, Case, NativeEngine};
use elmer_solver::SolverError;
use support::case::{fe, minimal_example, scratch_case, HandBuilt};
use support::logger;
use support::mesh::FileMesh;

/// `Solver 1 :: Reference Norm` in the pristine case.sif.
const REFERENCE_NORM: f64 = 2.88675135e-01;
/// `Reference Norm Tolerance` default (ElmerSolver.F90 CompareToReferenceSolution).
const TOLERANCE: f64 = 1.0e-5;

#[test]
fn minimal_example_norm_test_passed_and_vtu() {
    let _ = elmer_physics_heat::install_fe(fe());
    let dir = scratch_case("minimal_example", "direct");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let (log, buf) = logger();
    let case = Case {
        model: minimal_example(),
        mesh: &mesh,
        case_dir: dir.clone(),
        mesh_dir: dir.join("triangle"),
    };
    let out = run_case(case, &log, &SerialComm).unwrap();

    let [(1, norm)] = out.norms[..] else {
        panic!("{:?}", out.norms)
    };
    let rel = (norm - REFERENCE_NORM).abs() / REFERENCE_NORM;
    println!("MinimalExample native norm {norm:.17e}, reference {REFERENCE_NORM:e}, rel {rel:e}");
    assert!(rel <= TOLERANCE, "norm {norm}, rel {rel:e}");

    assert_eq!(out.passed, Some(true));
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
    assert!(buf.text().contains("PASSED all 1 tests!"), "{}", buf.text());

    assert_eq!(out.written, vec![dir.join("triangle/case_t0001.vtu")]);
    let vtu = std::fs::read(&out.written[0]).unwrap();
    let head = String::from_utf8_lossy(&vtu[..vtu.len().min(4096)]);
    assert!(head.contains("<VTKFile type=\"UnstructuredGrid\""), "{head}");
    assert!(head.contains("Name=\"temperature\""), "{head}");
    assert!(head.contains("Name=\"GeometryIds\""), "{head}");
}

fn drive(engine: &NativeEngine, dir: &std::path::Path) -> (Action, String, String) {
    let env = Environment::new(dir.to_path_buf()).with_var("ELMER_NATIVE", "1");
    let (log, buf) = logger();
    let mut diag = Vec::new();
    let action = run(
        Variant::Serial,
        &[OsString::from("case.sif")],
        &env,
        engine,
        &log,
        &mut diag,
    );
    (action, buf.text(), String::from_utf8(diag).unwrap())
}

#[test]
fn cli_driver_runs_minimal_example_natively() {
    let dir = scratch_case("minimal_example", "cli");
    let engine = NativeEngine::new(Box::new(HandBuilt(|_| {})));
    let (action, out, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Pass), "{out}\n{diag}");
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
    assert!(dir.join("triangle/case_t0001.vtu").is_file());
}

#[test]
fn cli_driver_refuses_without_side_effects() {
    // No front end in this build.
    let dir = scratch_case("minimal_example", "nofront");
    let (action, _, diag) = drive(&NativeEngine::default(), &dir);
    assert_eq!(action, Action::Exit(ExitCode::Refuse));
    assert!(diag.contains("no .sif and mesh front end"), "{diag}");

    // An iterative solve is not wired until elmer-krylov merges.
    let dir = scratch_case("minimal_example", "iterative");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].0.remove("linear system solver");
    })));
    let (action, _, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Refuse));
    assert!(diag.contains("iterative method bicgstab"), "{diag}");
    assert!(!dir.join("TEST.PASSED").exists());
    assert!(!dir.join("triangle/case_t0001.vtu").exists());

    // A procedure no crate provides is refused before setup, not skipped.
    let dir = scratch_case("minimal_example", "procedure");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].insert("Procedure", Value::Str("\"StatCurrentSolve\" \"StatCurrentSolver\"".into()));
    })));
    let (action, _, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Refuse));
    assert!(diag.contains("StatCurrentSolve"), "{diag}");
}

#[test]
fn transient_post_file_is_refused() {
    let _ = elmer_physics_heat::install_fe(fe());
    let dir = scratch_case("minimal_example", "transient");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut model = minimal_example();
    model.simulation.insert("Simulation Type", Value::Str("transient".into()));
    let (log, _) = logger();
    let case = Case { model, mesh: &mesh, case_dir: dir.clone(), mesh_dir: dir.join("triangle") };
    let err = run_case(case, &log, &SerialComm).unwrap_err();
    assert!(matches!(err, SolverError::Unsupported(_)), "{err}");
}
