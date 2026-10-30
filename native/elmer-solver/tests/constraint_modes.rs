//! `Constraint Modes Analysis` with `Constraint Modes Temperature = True` on
//! two boundaries, checked against the pristine ElmerSolver (26.2-devel, Rev
//! 9f6af2f85) run of `tests/data/cmodes/case.sif` (HeatSolve without
//! sources on LinControlCoordHeat's square mesh, the third boundary at
//! zero, `Max Output Level = 7`): the mode count, one constrained solve per
//! boundary dof (40 `NS` lines, the iteration count running through the
//! modes) and the steady-state change of the last mode; and of
//! `fluxes.sif` (the same with `Constraint Modes Fluxes`, `... Symmetric`
//! and `... Norm`): the lumped matrix, its norm and the `res:` entry.
//! HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn one_solve_per_mode_dof() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cmodes");
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
    let got = log.lines(&[
        "ComputeChange: NS",
        "ComputeChange: SS",
        "SetConstraintModesBoundaries: Number",
        "ConstraintModesDriver: Setting up",
    ]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            // RELC between symmetric modes is zero up to round-off.
            assert!((gr - wr).abs() <= 1e-9 * wr.max(1e-6), "{g}\n{w}");
        } else {
            assert_eq!(g, w);
        }
    }
}

/// The messages up to the sif's `Max Output Level = 7`.
struct Level7(Capture);

impl elmer_core::Logger for Level7 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 7 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

#[test]
fn lumped_fluxes_matrix() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cmodes");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("fluxes.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Level7(Capture::default());
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(dir.join("fluxes_stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.0.lines(&["StoreLumpedFluxes: Allocating", "FinalizeLumpedMatrix", "ComputeChange: NS", "ComputeChange: SS"]);
    let got: Vec<&String> = got.iter().filter(|l| !l.ends_with("did not converge to tolerance")).collect();
    assert_eq!(got.len(), want.len(), "{:#?}", &got[..got.len().min(8)]);
    let entry = |l: &str| -> f64 { l.split_whitespace().last().unwrap().parse().unwrap() };
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            assert!((gr - wr).abs() <= 1e-9 * wr.max(1e-6), "{g}\n{w}");
        } else if w.starts_with("FinalizeLumpedMatrix:  ") {
            assert_eq!(g[..28], w[..28]);
            // Entries against the matrix scale: the diagonal is round-off.
            assert!((entry(g) - entry(w)).abs() <= 1e-9, "{g}\n{w}");
        } else if w.contains("lumped matrix norm") {
            assert_eq!(g[..g.rfind(' ').unwrap()], w[..w.rfind(' ').unwrap()]);
            assert!((entry(g) - entry(w)).abs() <= 1e-9 * entry(w), "{g}\n{w}");
        } else {
            assert_eq!(g.as_str(), *w);
        }
    }
    let nrm = sim.model.simulation.real("res: temperature lumped matrix norm").unwrap();
    assert!((nrm - 43.9459).abs() < 1e-3, "{nrm}");
}
