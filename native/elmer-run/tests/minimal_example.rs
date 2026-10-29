//! MinimalExample end to end: hand-built case and mesh -> `elmer-assembly`
//! -> `elmer-direct` band LU -> `ComputeNorm` -> `TEST.PASSED` and the
//! `Post File` VTU, directly and through `elmer-cli`'s driver (hand-built
//! case, and the pristine `case.sif` through `elmer-sif`).

mod support;

use std::ffi::OsString;

use elmer_cli::driver::{run, Action, Environment};
use elmer_cli::strangler::Variant;
use elmer_cli::ExitCode;
use elmer_core::{SerialComm, Value};
use elmer_run::{run_case, Case, NativeEngine};
use support::case::{minimal_example, scratch_case, sif_frontend, HandBuilt};
use support::logger;
use support::mesh::FileMesh;

/// `Solver 1 :: Reference Norm` in the pristine case.sif.
const REFERENCE_NORM: f64 = 2.88675135e-01;
/// `Reference Norm Tolerance` default (ElmerSolver.F90 CompareToReferenceSolution).
const TOLERANCE: f64 = 1.0e-5;

#[test]
fn minimal_example_norm_test_passed_and_vtu() {
    let dir = scratch_case("minimal_example", "direct");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let (log, buf) = logger();
    let case = Case {
        model: minimal_example(),
        reloads: vec![],
        mesh: &mesh,
        case_dir: dir.clone(),
        mesh_dir: dir.join("triangle"),
    };
    let out = run_case(case, &log, &SerialComm).unwrap();

    let Some(&(_, norm)) = out.norms.iter().find(|(id, _)| *id == 1) else {
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
fn cli_driver_runs_the_pristine_sif_through_elmer_sif() {
    let dir = scratch_case("minimal_example", "cli-sif");
    let engine = NativeEngine::new(Box::new(sif_frontend()));
    let (action, out, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Pass), "{out}\n{diag}");
    assert!(out.contains("PASSED all 1 tests!"), "{out}");
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

    // A direct method elmer-direct does not provide.
    let dir = scratch_case("minimal_example", "spqr");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].insert("Linear System Direct Method", Value::Str("spqr".into()));
    })));
    let (action, _, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Refuse));
    assert!(diag.contains("spqr"), "{diag}");
    assert!(!dir.join("TEST.PASSED").exists());
    assert!(!dir.join("triangle/case_t0001.vtu").exists());

    // A procedure no crate provides is refused before setup, not skipped.
    let dir = scratch_case("minimal_example", "procedure");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].insert("Procedure", Value::Str("\"MagneticSolve\" \"MagneticSolver\"".into()));
    })));
    let (action, _, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Refuse));
    assert!(diag.contains("MagneticSolve"), "{diag}");
}

#[test]
fn missing_linear_system_solver_is_fatal_without_output() {
    // SolveLinearSystem's Fatal, as in Elmer: exit 1, no TEST.PASSED, no VTU.
    let dir = scratch_case("minimal_example", "nosolver");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].0.remove("linear system solver");
    })));
    let (action, out, _) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Fail));
    assert!(out.contains("Linear System Solver"), "{out}");
    assert!(!dir.join("TEST.PASSED").exists());
    assert!(!dir.join("triangle/case_t0001.vtu").exists());
}

#[test]
fn implicit_legacy_heat_equation_solver() {
    // `Equation = Heat Equation` with no Procedure or Variable: AddEquationBasics
    // supplies HeatSolve and Temperature.
    let dir = scratch_case("minimal_example", "legacy");
    let engine = NativeEngine::new(Box::new(HandBuilt(|m| {
        m.solvers[0].0.remove("procedure");
        m.solvers[0].0.remove("variable");
    })));
    let (action, out, diag) = drive(&engine, &dir);
    assert_eq!(action, Action::Exit(ExitCode::Pass), "{out}\n{diag}");
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
}

#[test]
fn run_blocks_reload_the_simulation() {
    // A second RUN block that only restates `Reference Norm` runs again
    // (Simulation::reload) and is checked again: two tests, both passed.
    let dir = scratch_case("minimal_example", "reload");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut block = elmer_solver::Model::default();
    let mut s1 = elmer_core::ParamMap::new();
    s1.insert("Reference Norm", Value::real(REFERENCE_NORM));
    block.solvers = vec![s1];
    let (log, buf) = logger();
    let case = Case {
        model: minimal_example(),
        reloads: vec![block],
        mesh: &mesh,
        case_dir: dir.clone(),
        mesh_dir: dir.join("triangle"),
    };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{e}\n{}", buf.text()));
    assert_eq!(out.passed, Some(true));
    assert!(buf.text().contains("PASSED all 2 tests!"), "{}", buf.text());
}

#[test]
fn run_block_changing_physics_is_refused_before_running() {
    let dir = scratch_case("minimal_example", "reload-refused");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut block = elmer_solver::Model::default();
    let mut mat = elmer_core::ParamMap::new();
    mat.insert("Heat Conductivity", Value::real(2.0));
    block.materials = vec![mat];
    let (log, _) = logger();
    let case = Case {
        model: minimal_example(),
        reloads: vec![block],
        mesh: &mesh,
        case_dir: dir.clone(),
        mesh_dir: dir.join("triangle"),
    };
    let err = run_case(case, &log, &SerialComm).unwrap_err().to_string();
    assert!(err.contains("RUN block 2 changes Material 1"), "{err}");
    assert!(!dir.join("TEST.PASSED").exists());
}

/// Runs MinimalExample with `edit` applied and returns its VTU.
fn run_minimal_with(tag: &str, edit: impl FnOnce(&mut elmer_solver::Model)) -> String {
    let dir = scratch_case("minimal_example", tag);
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut model = minimal_example();
    model.solvers[0].0.remove("reference norm");
    model.simulation.insert("Ascii Output", Value::Logical(true));
    edit(&mut model);
    let (log, buf) = logger();
    let case = Case { model, reloads: vec![], mesh: &mesh, case_dir: dir.clone(), mesh_dir: dir.join("triangle") };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{e}\n{}", buf.text()));
    std::fs::read_to_string(&out.written[0]).unwrap()
}

#[test]
fn mesh_levels_split_the_mesh_before_solving() {
    // One SplitMeshEqual of the triangle: 6 nodes, 4 triangles.
    let vtu = run_minimal_with("levels", |m| {
        m.simulation.insert("Mesh Levels", Value::Integer(vec![2]));
    });
    assert!(vtu.contains("NumberOfPoints=\"6\""), "{vtu}");
}

#[test]
fn transient_post_file_writes_each_output_step() {
    // The Post File solver runs `after saving`, i.e. at every output step.
    let dir = scratch_case("minimal_example", "transient");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut model = minimal_example();
    model.simulation.insert("Simulation Type", Value::Str("transient".into()));
    model.simulation.insert("Timestep Intervals", Value::Integer(vec![3]));
    model.simulation.insert("Timestep Sizes", Value::real(0.1));
    model.solvers[0].0.remove("reference norm");
    let (log, buf) = logger();
    let case = Case { model, reloads: vec![], mesh: &mesh, case_dir: dir.clone(), mesh_dir: dir.join("triangle") };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{e}\n{}", buf.text()));
    let want: Vec<_> = (1..=3).map(|i| dir.join(format!("triangle/case_t{i:04}.vtu"))).collect();
    assert_eq!(out.written, want);
    assert!(want.iter().all(|p| p.is_file()));
}
