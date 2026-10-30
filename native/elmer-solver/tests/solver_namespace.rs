//! The solver keyword `Namespace` (`ListPushNamespace` in SolverActivate,
//! MainUtils.F90:5816), checked against the pristine ElmerSolver (26.2-devel,
//! Rev 9f6af2f85) run of `tests/data/solver_namespace/case.sif`:
//! `fem/tests/heateq` on its ElmerGrid mesh with `heateq 1: Heat
//! Conductivity = 2` (instead of `equals nonlin iter`) and one nonlinear
//! iteration. Each solver sees its own `heateq k: Heat Conductivity`; the
//! material is restored after each activation. HeatSolve is the test-local
//! stand-in.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn each_solver_reads_its_namespaced_keywords() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/solver_namespace");
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
    let got: Vec<String> = log
        .lines(&["ComputeChange", "CompareToReferenceSolution: Solver"])
        .into_iter()
        .filter(|l| !l.ends_with("did not converge to tolerance"))
        .collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
        } else {
            assert_eq!(g, w);
        }
    }
    let m = &sim.model.materials[0];
    assert_eq!(m.get("heat conductivity"), None);
    assert_eq!(m.real("heateq 1: heat conductivity"), Some(2.0));
}
