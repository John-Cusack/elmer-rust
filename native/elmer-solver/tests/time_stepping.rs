//! Variable step sizes (`Timestep Size` as a function of `Time`) and
//! `Adaptive Timestepping`, checked against pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) runs of HeatSolve on the 1dtests mesh: BDF(1),
//! `T = 0` initially, `T = 1` and `T = 0` at the ends
//! (`tests/data/timestep/*.sif`, messages in `*_stdout.txt`). HeatSolve is
//! the test-local diffusion stand-in.
//!
//! The adaptive run covers 3 steps of 0.2 in 66 substep decisions
//! (`Adaptive(cum,ddt,err)`): a full step and two half steps each, the
//! step size halving on rejection and doubling after `Adaptive Keep
//! Smallest` accepted substeps.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat;
use support::{sif_model, split_change, Capture};

fn numbers(l: &str) -> Vec<f64> {
    l.rsplit_once("Adaptive(cum,ddt,err):")
        .unwrap()
        .1
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

fn check(sif: &str, stdout: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let model = sif_model(&dir.join("timestep").join(sif));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let keep = |l: &&str| l.contains("(NRM,RELC)") || l.contains("Adaptive(") || l.starts_with("MAIN: Time");
    let got = log.lines(&["ComputeChange", "ExecSimulation", "MAIN"]);
    let got: Vec<&str> = got.iter().map(String::as_str).filter(keep).collect();
    let want: Vec<&str> = stdout.lines().filter(keep).collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            assert!((gr - wr).abs() <= 1e-7 * wr, "line {i}\n{g}\n{w}");
        } else if w.contains("Adaptive(") {
            for (a, b) in numbers(g).iter().zip(numbers(w)) {
                assert!((a - b).abs() <= 1e-9 * b.abs().max(1e-12), "line {i}\n{g}\n{w}");
            }
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "{got:#?}");
}

#[test]
fn timestep_size_follows_time() {
    check("varstep.sif", include_str!("data/timestep/varstep_stdout.txt"));
}

#[test]
fn adaptive_timestepping() {
    check("adapt.sif", include_str!("data/timestep/adapt_stdout.txt"));
}
