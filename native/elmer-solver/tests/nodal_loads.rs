//! `Calculate Loads` (`<var> Loads` from the bulk system after the solve)
//! and `<var> Load` nodal loads, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `fem/tests/DirichletNeumann` with its
//! real `case.sif` (`Post File` commented out) on `ElmerGrid 1 2
//! internal_slope` (`tests/data/dn_loads/`).
//!
//! Solver 1 (`TempB`, upper body) computes `TempB Loads`; solver 2 (`TempA`,
//! lower body) takes `TempA Load = Opposes "TempB Loads"` on the interface
//! while `TempB = Equals TempA` there. Without the loads `TempA` stays zero.
//! HeatSolve is the test-local diffusion stand-in (unit conductivity).

mod support;

use std::path::Path;

use elmer_core::{SerialComm, Value};
use elmer_solver::{finalize_reference, Simulation, SolverRegistry};
use support::poisson::heat;
use support::{sif_model, split_change, Capture};

#[test]
fn dirichlet_neumann_through_nodal_loads() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dn_loads");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    // AddEquationBasics created the loads field with TempB's perm.
    let loads = sim.vars.get("tempb loads").expect("TempB Loads");
    assert_eq!(loads.perm, sim.vars.get("tempb").unwrap().perm);
    sim.run(&mesh, &log, &SerialComm).unwrap();
    finalize_reference(&sim.reference, Path::new(env!("CARGO_TARGET_TMPDIR")), &log, &SerialComm).unwrap();

    let got = log.lines(&["ComputeChange", "CompareToReferenceSolution"]);
    // Elmer printed at Max Output Level 5; the capture keeps every level.
    let skip = |l: &&str| !l.contains("Relative Error") && (!l.starts_with("ComputeChange") || l.contains("(NRM,RELC)"));
    let got: Vec<&str> = got.iter().map(String::as_str).filter(skip).collect();
    let want: Vec<&str> = include_str!("data/dn_loads/stdout.txt").lines().filter(skip).collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            let tol = if wr < 1e-6 { 1e-10 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "{got:#?}");
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
}

#[test]
fn calculate_loads_slot_other_than_bulk_is_refused() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dn_loads");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.solvers[0].insert("Calculate Loads Slot", Value::Str("boundary assembly".into()));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &Capture::default(), &SerialComm).unwrap();
    let err = sim.run(&mesh, &Capture::default(), &SerialComm).unwrap_err();
    assert!(err.to_string().contains("Calculate Loads Slot"), "{err}");
}
