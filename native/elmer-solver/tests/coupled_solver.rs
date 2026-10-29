//! `Coupled Solver = True` (MainUtils.F90 `CoupledSolver`), checked against
//! the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) runs of
//! `fem/tests/CoupledPoisson3`..`9`: the real `case.sif` files
//! (`tests/data/coupled/cpN.sif`, read with `elmer-sif`), meshes from
//! `ElmerGrid 1 2 angle|squares`, and each run's `CoupledSolver`,
//! `BlockSystemAssembly`, `PoissonBulkAssembly`, `ComputeChange` and
//! `CompareToReferenceSolution` lines (`cpN_stdout.txt`).
//!
//! - 3: two block variables, no constraints;
//! - 4, 5: an integral constraint per variable over the body forces;
//! - 6: an integral constraint over a boundary;
//! - 7: an equality constraint between two variables on a boundary;
//! - 8, 9: a floating constraint on a boundary, without and with a value.

mod support;

use std::path::Path;
use std::sync::Arc;

use elmer_core::SerialComm;
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::fixtures::Bulk;
use support::{sif_model, split_change, Capture};

fn run_case(k: usize, mesh_dir: &str, stdout: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join(mesh_dir), 2);
    let model = sif_model(&dir.join(format!("coupled/cp{k}.sif")));
    let mut reg = SolverRegistry::new();
    reg.register_assembly("PoissonAssembly", "BulkAssembly", Arc::new(|_| Ok(Box::new(Bulk::default()))));
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("CoupledPoisson{k}: {e}"));

    let prefixes = ["CoupledSolver", "BlockSystemAssembly", "PoissonBulkAssembly", "ComputeChange", "CompareTo"];
    let got = log.lines(&prefixes);
    // `PASSED all 1 tests!` is `finalize_reference`'s, after the run.
    let want: Vec<&str> = stdout.lines().filter(|l| !l.contains("Relative Error") && !l.contains("all 1 tests")).collect();
    let got: Vec<&String> = got.iter().filter(|l| !l.contains("Relative Error")).collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "CoupledPoisson{k} line {i}");
            let tol = if wr < 1e-6 { 1e-8 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "CoupledPoisson{k} line {i}\n{g}\n{w}");
        } else if w.starts_with("CompareTo") && w.contains("Norm =") {
            // The 9-digit norm may differ in the last digit after an iterative solve.
            assert_eq!(g[..w.find("Norm =").unwrap()], w[..w.find("Norm =").unwrap()], "CoupledPoisson{k}");
        } else {
            assert_eq!(*g, w, "CoupledPoisson{k} line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "CoupledPoisson{k}: {got:#?}");
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "CoupledPoisson{k}");
}

#[test]
fn coupled_poisson3_two_variables() {
    run_case(3, "angle", include_str!("data/coupled/cp3_stdout.txt"));
}

#[test]
fn coupled_poisson4_integral_constraints() {
    run_case(4, "angle", include_str!("data/coupled/cp4_stdout.txt"));
}

#[test]
fn coupled_poisson5_integral_constraints_direct() {
    run_case(5, "angle", include_str!("data/coupled/cp5_stdout.txt"));
}

#[test]
fn coupled_poisson6_boundary_integral_constraint() {
    run_case(6, "angle", include_str!("data/coupled/cp6_stdout.txt"));
}

#[test]
fn coupled_poisson7_equality_constraint() {
    run_case(7, "angle", include_str!("data/coupled/cp7_stdout.txt"));
}

#[test]
fn coupled_poisson8_floating_constraint() {
    run_case(8, "squares", include_str!("data/coupled/cp8_stdout.txt"));
}

#[test]
fn coupled_poisson9_floating_constraint_with_value() {
    run_case(9, "squares", include_str!("data/coupled/cp9_stdout.txt"));
}
