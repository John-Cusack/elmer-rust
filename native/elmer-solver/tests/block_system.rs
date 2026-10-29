//! `Linear System Solver = Block`, dispatched by `solve_linear` to
//! `elmer-block` after the scaling (`BlockSolveExt`, SolveCore.F90:1663), on
//! the interleaved real Helmholtz system of `tests/data/complex/direct.sif`
//! split into its two components (UMFPACK blocks inside an outer GCR).
//! HelmholtzSolve is the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Simulation, SolverError, SolverRegistry};
use support::poisson::helmholtz;
use support::{sif_model, Capture};

fn run(extra: &[(&str, Value)]) -> Result<(Simulation, Capture), SolverError> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let mut model = sif_model(&dir.join("complex/direct.sif"));
    for (k, v) in extra {
        model.solvers[0].insert(k, v.clone());
    }
    let mut reg = SolverRegistry::new();
    reg.register("HelmholtzSolve", "HelmholtzSolver", helmholtz);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm)?;
    sim.run(&mesh, &log, &SerialComm)?;
    Ok((sim, log))
}

fn block_keywords() -> Vec<(&'static str, Value)> {
    let s = |v: &str| Value::Str(v.into());
    vec![
        ("Linear System Solver", s("block")),
        ("Block Preconditioner", Value::Logical(true)),
        ("Outer: Linear System Iterative Method", s("gcr")),
        ("Outer: Linear System Max Iterations", Value::Integer(vec![100])),
        ("Outer: Linear System Convergence Tolerance", Value::real(1e-12)),
        ("block 11: Linear System Solver", s("direct")),
        ("block 11: Linear System Direct Method", s("umfpack")),
        ("block 22: Linear System Solver", s("direct")),
        ("block 22: Linear System Direct Method", s("umfpack")),
    ]
}

#[test]
fn block_systems_go_to_elmer_block() {
    // With and without the monolithic scaling, which the outer iteration
    // then gets as well.
    for scaling in [true, false] {
        let mut kw = block_keywords();
        kw.push(("Linear System Scaling", Value::Logical(scaling)));
        let (sim, log) = run(&kw).unwrap();
        assert!(sim.vars.get("pressure").unwrap().values.iter().all(|v| v.is_finite()));
        assert_eq!(log.lines(&["BlockInitMatrix: Block matrix will be of size 2"]).len(), 1);
        assert_eq!(log.lines(&["BlockKrylovIter: Solution norm:"]).len(), 1);
    }
}

#[test]
fn a_missing_schur_matrix_is_fatal() {
    let mut kw = block_keywords();
    kw.push(("Block Matrix Schur Variable", Value::Str("nosuch".into())));
    let e = run(&kw).err().unwrap();
    assert!(matches!(&e, SolverError::Fatal { caller, .. } if caller == "BlockPrecMatrix"), "{e}");
}

#[test]
fn block_with_complex_arithmetic_is_refused() {
    let mut kw = block_keywords();
    kw.push(("Linear System Complex", Value::Logical(true)));
    let e = run(&kw).err().unwrap();
    assert!(matches!(e, SolverError::Unsupported(_)), "{e}");
}
