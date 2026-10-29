//! `Linear System Namespace Number`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `fem/tests/LinearSystemNamespace` from
//! its `case.sif` on its ElmerGrid mesh (`tests/data/namespace/`): the
//! number follows `nonlin iter`, so the first nonlinear iteration solves
//! with the `linsys1:` keywords (umfpack) and the second with `linsys2:`
//! (BiCGStab). HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn namespace_follows_the_nonlinear_iteration() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/namespace");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["DefaultSolve: Linear system namespace", "ComputeChange", "CompareToReferenceSolution: Solver"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            // ITER=2 changes by round-off only.
            let tol = if wr < 1e-6 { 1e-12 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    // The namespaced keywords stand in only during the solve.
    let p = &sim.model.solvers[0];
    assert_eq!(p.string("linear system solver"), None);
    assert_eq!(sim.slots[0].params.string("linear system solver"), None);
    assert_eq!(sim.slots[0].params.string("linsys2: linear system solver"), Some("iterative"));
}
