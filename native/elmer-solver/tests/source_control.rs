//! `Apply Source Control`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) runs of `fem/tests/LinControlCoordHeat`
//! (target at a node by coordinates), `LinControlMaxHeat` (target at the
//! extremum) and `LinControlFourHeaters` (four sources, four targets) from
//! their `case.sif` files on their ElmerGrid meshes (`tests/data/control/`),
//! and of `lincontrolloadheat` (LinControlCoordHeat's mesh, `Control Use
//! Loads`: the target is the nodal load at a Dirichlet node, direct solves):
//! the `ComputeChange` lines, every `Applied Control Amplitude` and the
//! reference-norm check. HeatSolve is the test-local nonlinear stand-in.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

fn check(case: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/control").join(case);
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
    sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("{case}: {e}"));

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let keep = |l: &&str| {
        (l.starts_with("ComputeChange") && l.contains("(NRM,RELC)"))
            || l.starts_with("ControlLinearSystem: Applied")
            || (l.starts_with("CompareTo") && !l.contains("Relative Error") && !l.contains("all 1 tests"))
    };
    let got = log.lines(&["ComputeChange", "ControlLinearSystem", "CompareTo"]);
    let got: Vec<&str> = got.iter().map(String::as_str).filter(keep).collect();
    let want: Vec<&str> = stdout.lines().filter(keep).collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "{case} line {i}");
            let tol = if wr < 1e-6 { 1e-8 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "{case} line {i}\n{g}\n{w}");
        } else if w.starts_with("CompareTo") {
            // The 9-digit norm may differ in its last digit after an iterative solve.
            let cut = w.find("Norm =").unwrap() + 14;
            assert_eq!(g[..cut], w[..cut], "{case} line {i}");
        } else {
            assert_eq!(g, w, "{case} line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "{case}: {got:#?}");
    if case != "lincontrolloadheat" {
        assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "{case}");
    }
}

#[test]
fn control_at_a_node_given_by_coordinates() {
    check("lincontrolcoordheat");
}

#[test]
fn control_at_the_extremum() {
    check("lincontrolmaxheat");
}

#[test]
fn four_controls() {
    check("lincontrolfourheaters");
}

#[test]
fn control_on_the_nodal_loads() {
    check("lincontrolloadheat");
}
