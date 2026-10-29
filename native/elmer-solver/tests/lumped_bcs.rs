//! The lumped Dirichlet conditions of `DefaultDirichletBCs`, run from the
//! `case.sif` files of `fem/tests/FixPotentialLevels` (`Temperature Single
//! Node` on three body forces), `ConstantBCTemperature` (`Temperature
//! Constant`), `ProfileBCTemperature` and `ProfileBCTemperatureRobin`
//! (`Temperature Profile`) on their ElmerGrid meshes (`tests/data/lumped/`).
//! FixPotentialLevels is checked against the pristine ElmerSolver (26.2-devel,
//! Rev 9f6af2f85) run: `ComputeChange` lines and the reference norm. The
//! pristine build segfaults in `List_MoveRow` (ListMatrix.F90:1089 tests
//! `p2 % Index` of a null `p2`) on the other three, so they are checked
//! against their cases' own reference norms and the lumped structure of the
//! solution. HeatSolve is the test-local stand-in; the SaveScalars solver of
//! the cases is left out.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::{MeshView, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

fn run(case: &str) -> (Simulation, support::mesh::FileMesh, Capture) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/lumped").join(case);
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.solvers.truncate(1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("{case}: {e}"));
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "{case}");
    // All nine digits of the reference norm.
    let l = log.lines(&["CompareToReferenceSolution: Solver 1"]).pop().unwrap();
    let (norm, refnorm) = l.split_once("Norm = ").unwrap().1.split_once("  RefNorm = ").unwrap();
    assert_eq!(norm, refnorm, "{case}");
    (sim, mesh, log)
}

/// The nodes of boundary `tag`, sorted.
fn boundary_nodes(mesh: &support::mesh::FileMesh, tag: i32) -> Vec<u32> {
    let nb = mesh.n_bulk() as u32;
    let mut v: Vec<u32> = (nb..nb + mesh.n_boundary() as u32)
        .map(elmer_core::ElemId)
        .filter(|&e| mesh.element(e).tag == tag)
        .flat_map(|e| mesh.element(e).nodes.to_vec())
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

fn temperature(sim: &Simulation, node: u32) -> f64 {
    let v = sim.vars.get("temperature").unwrap();
    v.values[v.perm.as_deref().unwrap()[node as usize] as usize]
}

#[test]
fn single_nodes_fix_three_floating_bodies() {
    let (_, _, log) = run("singlenode");
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/lumped/singlenode");
    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange", "CompareToReferenceSolution: Solver"]);
    let got: Vec<&str> = got.iter().map(String::as_str).filter(|l| !l.contains("Relative Error")).collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            // ITER=2 changes by round-off only.
            let tol = if wr < 1e-6 { 1e-9 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
}

#[test]
fn constant_value_on_a_boundary() {
    let (sim, mesh, log) = run("constant");
    let nodes = boundary_nodes(&mesh, 2);
    let t0 = temperature(&sim, nodes[0]);
    assert!(t0 > 0.0);
    for &k in &nodes {
        assert!((temperature(&sim, k) - t0).abs() <= 1e-9 * t0, "node {k}");
    }
    assert!(log.lines(&["SetDirichletBoundaries: Previous bc 2 lumped flux:"]).len() == 1);
}

fn profile(case: &str) {
    let (sim, mesh, _) = run(case);
    // Temperature Profile = 0.75 - 0.25 x, largest at the smallest x.
    let nodes = boundary_nodes(&mesh, 2);
    let p = |k: u32| 0.75 - 0.25 * mesh.coord(k)[0];
    let k0 = *nodes.iter().max_by(|&&a, &&b| p(a).abs().total_cmp(&p(b).abs())).unwrap();
    let t0 = temperature(&sim, k0);
    for &k in &nodes {
        let want = t0 * p(k) / p(k0);
        assert!((temperature(&sim, k) - want).abs() <= 1e-9 * t0.abs(), "node {k}");
    }
}

#[test]
fn profile_on_a_boundary() {
    profile("profile");
}

#[test]
fn profile_on_a_robin_boundary() {
    profile("profilerobin");
}
