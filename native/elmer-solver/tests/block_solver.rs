//! `Block Solver = True` with a `Procedure` (`BlockSolver`,
//! MainUtils.F90:4323).
//!
//! - HeatSolve with `Block Solver` on the LinearSystemNamespace mesh
//!   (`tests/data/blockdrv/`) follows the pristine ElmerSolver (26.2-devel,
//!   Rev 9f6af2f85) run up to its fatal: one block cannot be picked.
//! - Without a Procedure and without `Variable 1`, BlockSolver's fatal.
//! - The two-component Helmholtz system of `tests/data/complex/direct.sif`:
//!   one coupled iteration is the `Linear System Solver = Block` solve of
//!   the same scaled system through `elmer-block`, bit for bit, and the
//!   keywords the driver overrides are put back.
//!
//! HeatSolve and HelmholtzSolve are the test-local stand-ins.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm, Value};
use elmer_solver::{Simulation, SolverError, SolverRegistry};
use support::poisson::{heat_nl, helmholtz};
use support::{sif_model, Capture};

/// The messages up to `Max Output Level = 10`, as the pristine run.
struct Level10(Capture);

impl elmer_core::Logger for Level10 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 10 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

#[test]
fn one_block_stops_where_elmer_does() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("namespace"), 2);
    let model = sif_model(&dir.join("blockdrv/case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Level10(Capture::default());
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let e = sim.run(&mesh, &log, &SerialComm).err().unwrap();
    let mut got = log.0.lines(&["BlockSolver", "DefaultSolve: Solver disabled"]);
    got.push(format!("ERROR:: {e}"));
    let stdout = std::fs::read_to_string(dir.join("blockdrv/stdout.txt")).unwrap();
    assert_eq!(got, stdout.lines().collect::<Vec<_>>());
}

fn helmholtz_run(extra: &[(&str, Value)]) -> Result<Simulation, SolverError> {
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
    Ok(sim)
}

fn block_keywords() -> Vec<(&'static str, Value)> {
    let s = |v: &str| Value::Str(v.into());
    vec![
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
fn a_coupled_iteration_is_the_block_solve() {
    let mut kw = block_keywords();
    kw.push(("Linear System Solver", Value::Str("block".into())));
    let want = helmholtz_run(&kw).unwrap();

    let mut kw = block_keywords();
    kw.push(("Block Solver", Value::Logical(true)));
    let sim = helmholtz_run(&kw).unwrap();
    assert_eq!(
        sim.vars.get("pressure").unwrap().values,
        want.vars.get("pressure").unwrap().values
    );
    let p = &sim.slots[0].params;
    assert_eq!(p.get("linear system solver disabled"), None);
    assert_eq!(p.integer("nonlinear system max iterations"), Some(1));
}

#[test]
fn block_solver_without_procedure_needs_variables() {
    // Without a Procedure the blocks are `Variable 1..n`; none is BlockSolver's
    // fatal.
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let mut model = sif_model(&dir.join("complex/direct.sif"));
    model.solvers[0].0.remove("procedure");
    model.solvers[0].insert("Block Solver", Value::Logical(true));
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    let r = sim.setup(&mesh, &SolverRegistry::new(), &log, &SerialComm);
    assert!(
        matches!(&r, Err(SolverError::Fatal { caller, msg }) if caller == "BlockSolver" && msg == "No variables, nothing to do!"),
        "{r:?}"
    );
}
