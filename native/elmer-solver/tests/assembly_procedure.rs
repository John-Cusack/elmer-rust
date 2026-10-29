//! `Assembly Solver = True` with a `Bulk Assembly Procedure`, checked
//! against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/CoupledPoisson1`: mesh from `ElmerGrid 1 2 angle`
//! (`tests/data/angle/`), the run's messages in
//! `tests/data/coupled_poisson1_stdout.txt`. The `PoissonAssembly.f90` kernel
//! is `elmer-fixtures`' port; basis and quadrature are the test stand-ins.

mod support;

use std::path::Path;
use std::sync::Arc;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::fixtures::Bulk;
use support::{i, r, s, sec, split_change, Capture};

const STDOUT: &str = include_str!("data/coupled_poisson1_stdout.txt");

fn model() -> Model {
    let solver = sec(&[
        ("Assembly Solver", Value::Logical(true)),
        ("Equation", s("poisson")),
        ("Variable", s("TempA")),
        ("Steady State Convergence Tolerance", r(1.0e-5)),
        ("Nonlinear System Relaxation Factor", r(1.0)),
        ("Nonlinear System Max Iterations", i(&[100])),
        ("Nonlinear System Max Stepsize Tests", i(&[30])),
        ("Nonlinear System Convergence Tolerance", r(1.0e-6)),
        ("Linear System Solver", s("iterative")),
        ("Linear System Iterative Method", s("bicgstab")),
        ("Linear System Convergence Tolerance", r(1.0e-8)),
        ("Linear System Max Iterations", i(&[500])),
        ("Linear System Preconditioning", s("ilu1")),
        ("Bulk Assembly Procedure", Value::File("\"./PoissonAssembly\" \"BulkAssembly\"".into())),
        ("Reference Norm", r(0.77824031e-01)),
    ]);
    let body = sec(&[("Equation", i(&[1])), ("Material", i(&[1]))]);
    Model {
        simulation: sec(&[
            ("Coordinate System", s("cartesian")),
            ("Simulation Type", s("steady state")),
            ("Steady State Max Iterations", i(&[1])),
            ("Steady State Min Iterations", i(&[1])),
            ("Output Intervals", i(&[1])),
        ]),
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![body.clone(), body],
        materials: vec![sec(&[("Conductivity", r(1.0)), ("Reaction", r(0.0)), ("Source", r(1.0))])],
        boundary_conditions: vec![sec(&[("Target Boundaries", i(&[1, 2])), ("TempA", r(0.0))])],
        ..Model::default()
    }
}

#[test]
fn coupled_poisson1_matches_elmer() {
    let mesh = support::mesh::FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/angle"), 2);
    let mut reg = SolverRegistry::new();
    reg.register_assembly("PoissonAssembly", "BulkAssembly", Arc::new(|_| Ok(Box::new(Bulk::default()))));
    let mut sim = Simulation::new(model());
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));

    let prefixes = ["CoupledSolver", "BlockSystemAssembly", "PoissonBulkAssembly", "ComputeChange"];
    let got = log.lines(&prefixes);
    let want: Vec<&str> = STDOUT.lines().filter(|l| prefixes.iter().any(|p| l.starts_with(p))).collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if !w.starts_with("ComputeChange") {
            assert_eq!(g, w);
            continue;
        }
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        // The second iterate's change is BiCGStab round-off (7e-11).
        let tol = if wr < 1e-6 { 1e-8 } else { 1e-7 * wr };
        assert!((gr - wr).abs() <= tol, "{g}\n{w}");
    }
}

#[test]
fn missing_and_unregistered_routines_are_fatal() {
    let mut m = model();
    m.solvers[0].0.remove("bulk assembly procedure");
    let mesh = support::mesh::FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/angle"), 2);
    let log = Capture::default();
    let err = Simulation::new(m).setup(&mesh, &SolverRegistry::new(), &log, &SerialComm).unwrap_err();
    assert_eq!(err.to_string(), "BlockSystemAssembly: Bulk Assembly Procedure not given!");
    let err = Simulation::new(model()).setup(&mesh, &SolverRegistry::new(), &log, &SerialComm).unwrap_err();
    assert_eq!(
        err.to_string(),
        "BlockSystemAssembly: Assembly routine not found: ./PoissonAssembly BulkAssembly"
    );
    let mut m = model();
    m.solvers[0].insert("Block Solver", Value::Logical(true));
    m.solvers[0].insert("Assembly Solver", Value::Logical(false));
    // A Block Solver counts its blocks from `Variable 1..n`; none here.
    let err = Simulation::new(m).setup(&mesh, &SolverRegistry::new(), &log, &SerialComm).unwrap_err();
    assert_eq!(err.to_string(), "BlockSolver: No variables, nothing to do!");
}
