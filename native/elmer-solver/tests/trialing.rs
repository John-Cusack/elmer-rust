//! `Linear System Trialing`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `fem/tests/LinearSystemNamespace2`
//! (`tests/data/trialing/case.sif` on the ElmerGrid mesh in
//! `tests/data/namespace/`): two scanning steps whose solve tries BiCGStab
//! without preconditioning (`linsys1:`), then with ILU0 (`linsys2:`), both
//! cut at five iterations without aborting, and converges with umfpack
//! (`linsys3:`); each try computes its change. With `Linear System Trialing
//! Conserve` (`stdout_conserve.txt`, the same run with the sif's commented
//! line on) the second step starts from `linsys3:`. HeatSolve is the
//! test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm, Value};
use elmer_solver::{FeBackend, Simulation, SolverError, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

fn run(edit: impl FnOnce(&mut elmer_solver::Model)) -> (Result<(), SolverError>, Simulation, Capture) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("namespace"), 2);
    let mut model = sif_model(&data.join("trialing/case.sif"));
    edit(&mut model);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let r = sim.run(&mesh, &log, &SerialComm).map(|_| ());
    (r, sim, log)
}

#[test]
fn trialing_moves_on_until_the_linear_solve_converges() {
    let (r, sim, log) = run(|_| {});
    r.unwrap();
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let stdout = std::fs::read_to_string(data.join("trialing/stdout.txt")).unwrap();
    let keep = ["ComputeChange", "CompareToReferenceSolution: Solver", "CompareToReferenceSolution: Relative"];
    let want: Vec<&str> = stdout.lines().filter(|l| keep.iter().any(|k| l.starts_with(k))).collect();
    // The not-converged notice is level 6, above the sif's Max Output Level 5.
    let got: Vec<String> = log
        .lines(&keep)
        .into_iter()
        .filter(|l| !l.ends_with("did not converge to tolerance"))
        .collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            assert!((gr - wr).abs() <= 1e-6 * wr, "line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    let tries: Vec<String> = log.lines(&["DefaultSolve: Linear system namespace"]);
    let want_tries: Vec<String> = [1, 2, 3, 1, 2, 3]
        .iter()
        .map(|n| format!("DefaultSolve: Linear system namespace number: {n}"))
        .collect();
    assert_eq!(tries, want_tries);
    // Without Conserve the next solve starts from the first namespace again.
    assert_eq!(sim.slots[0].params.get("linear system namespace number"), None);
}

#[test]
fn conserve_keeps_the_converged_namespace() {
    let (r, sim, log) = run(|m| {
        m.solvers[0].insert("Linear System Trialing Conserve", Value::Logical(true));
    });
    r.unwrap();
    let tries = log.lines(&["DefaultSolve: Linear system namespace"]);
    let n: Vec<&str> = tries.iter().map(|l| l.rsplit(' ').next().unwrap()).collect();
    assert_eq!(n, ["1", "2", "3", "3"]);
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let stdout = std::fs::read_to_string(data.join("trialing/stdout_conserve.txt")).unwrap();
    let got = log.lines(&["ComputeChange: NS"]);
    assert_eq!(got.len(), stdout.lines().count(), "{got:#?}");
    for (g, w) in got.iter().zip(stdout.lines()) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-6 * wr, "{g}\n{w}");
    }
    assert_eq!(sim.slots[0].params.real("linear system namespace number"), Some(3.0));
}

#[test]
fn running_out_of_namespaces_is_fatal() {
    let (r, _, _) = run(|m| {
        let p = &mut m.solvers[0];
        p.0.retain(|k, _| !k.starts_with("linsys3"));
    });
    match r {
        Err(SolverError::Fatal { msg, .. }) => assert_eq!(msg, "Exhausted all linear system strategies!"),
        other => panic!("expected the strategies exhausted, got {other:?}"),
    }
}
