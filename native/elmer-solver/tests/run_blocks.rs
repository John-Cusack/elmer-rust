//! Several `RUN` blocks in one simulation (`Simulation::reload`), checked
//! against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/linearsolvers`: 14 blocks, each changing the linear solver
//! and the boundary temperature (`T = k` in block `k`), on the case's own
//! mesh (`tests/data/linearsolvers/`). The heat solver is the test-local
//! steady diffusion stand-in.
//!
//! Elmer prints RELC = 2/(2k-1) in block `k`: the field and its norm carry
//! over from the previous block, which only holds if nothing is reset but
//! the time globals.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{finalize_reference, Simulation, SolverRegistry};
use support::poisson::heat;
use support::{sif_runs, split_change, Capture};

#[test]
fn linearsolvers_fourteen_run_blocks() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/linearsolvers");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut runs = sif_runs(&dir.join("TempDist.sif"));
    assert_eq!(runs.len(), 14);
    // HeatSolver_Init's default variable name.
    runs[0].solvers[0].insert("Variable", elmer_core::Value::Str("temperature".into()));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let log = Capture::default();
    let mut sim = Simulation::new(runs[0].clone());
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    for (k, block) in runs.iter().enumerate() {
        if k > 0 {
            sim.reload(block, &mesh, &reg, &log, &SerialComm).unwrap();
        }
        sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("block {}: {e}", k + 1));
    }
    finalize_reference(&sim.reference, Path::new(env!("CARGO_TARGET_TMPDIR")), &log, &SerialComm).unwrap();

    let got = log.lines(&["ComputeChange", "CompareToReferenceSolution"]);
    let want: Vec<&str> = include_str!("data/linearsolvers/stdout.txt").lines().collect();
    let skip = |l: &&str| !l.contains("Relative Error");
    let got: Vec<&str> = got.iter().map(String::as_str).filter(skip).collect();
    let want: Vec<&str> = want.into_iter().filter(skip).collect();
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            assert!((gr - wr).abs() <= 1e-7 * wr, "{g}\n{w}");
        } else {
            assert_eq!(g, w);
        }
    }
    assert_eq!(got.len(), want.len(), "{got:#?}");
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (14, 14));
}
