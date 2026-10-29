//! `Default Internal BC` / `Default External BC`, checked against the
//! pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/DefaultIntExtBCs` from its `case.sif` on its ElmerGrid mesh
//! (`tests/data/defaultbcs/`): the boundaries no condition targets go to
//! the internal (two parents, 20 elements) and external (one parent, 60
//! elements) sections, with the run's MapInternalExternalBCs messages,
//! `ComputeChange` lines and reference norm. HeatSolve is the test-local
//! stand-in.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

/// The messages up to `Max Output Level = 8`, as the pristine run.
struct Level8(Capture);

impl elmer_core::Logger for Level8 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 8 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

#[test]
fn untargeted_boundaries_get_the_default_sections() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/defaultbcs");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Level8(Capture::default());
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.0.lines(&["MapInternal", "ComputeChange", "CompareToReferenceSolution: Solver"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            assert!((gr - wr).abs() <= 1e-6 * wr, "line {i}\n{g}\n{w}");
        } else if w.starts_with("CompareTo") {
            let cut = w.find("Norm =").unwrap() + 14;
            assert_eq!(g[..cut], w[..cut], "line {i}");
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
    // The mapping is the model's, for every reader of bc_of_tag.
    assert_eq!((sim.model.bc_of_tag(1), sim.model.bc_of_tag(2), sim.model.bc_of_tag(3)), (Some(3), Some(1), Some(2)));
}
