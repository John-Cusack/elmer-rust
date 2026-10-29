//! `Set Dirichlet BCs by BC Numbering` and `Dirichlet BC Order`, checked
//! against pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) runs of
//! `tests/data/bcorder/{element,bc,order}.sif`: HeatSolve on the
//! IntegralConstraintPoisson angle mesh (`tests/data/intbc/poisson/`) with
//! three conditions numbered against the element order, so that the value
//! at the nodes boundary 1 shares with 2 and 3 depends on which is set last.
//! HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

fn check(case: &str) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("intbc/poisson"), 2);
    let model = sif_model(&data.join(format!("bcorder/{case}.sif")));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    let stdout = std::fs::read_to_string(data.join(format!("bcorder/{case}.txt"))).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange"]);
    assert_eq!(got.len(), want.len(), "{case}: {got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt), "{case}");
        assert!((gr - wr).abs() <= 1e-9 * wr, "{case}\n{g}\n{w}");
    }
}

#[test]
fn element_order() {
    check("element");
}

#[test]
fn bc_numbering_order() {
    check("bc");
}

#[test]
fn given_bc_order() {
    check("order");
}
