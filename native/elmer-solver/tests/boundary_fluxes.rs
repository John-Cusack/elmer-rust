//! `Calculate Boundary Fluxes`: per boundary condition sums of the nodal
//! loads, stored as `res: <var> Flux over BC i` for SaveScalars. Checked
//! against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `fem/tests/BoundaryFluxes` (its `case.sif`, `ElmerGrid 1 2 angle`;
//! `tests/data/bflux/`): the `CalculateLoads` messages and the two values
//! SaveScalars reports. HeatSolve is the test-local stand-in; SaveScalars
//! (solver 2) is not registered here, the values are read from the
//! Simulation section where it finds them.

mod support;

use std::path::Path;

use elmer_core::{Params, SerialComm};
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn boundary_fluxes_match_elmer() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bflux");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    assert!(sim.vars.get("temperature loads").is_some());
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = include_str!("data/bflux/stdout.txt");
    let want: Vec<&str> = stdout.lines().filter(|l| !l.starts_with("SaveScalars")).collect();
    let got = log.lines(&["CalculateLoads", "ComputeChange: NS"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            assert!((gr - wr).abs() <= 1e-9 * wr);
        } else if w.contains("over all BCs") {
            // The sum of the two is iterative-solver round-off (1e-8).
            let v: f64 = g.rsplit_once(": ").unwrap().1.trim().parse().unwrap();
            assert!(v.abs() < 1e-6, "{g}");
        } else {
            assert_eq!(g, w);
        }
    }
    // SaveScalars' 13-digit values.
    for l in stdout.lines().filter(|l| l.starts_with("SaveScalars")) {
        let (key, val) = l.split_once("res: ").unwrap().1.trim().rsplit_once(' ').unwrap();
        let want: f64 = val.parse().unwrap();
        let got = sim.model.simulation.real(&format!("res: {}", key.trim())).unwrap();
        assert!((got - want).abs() <= 1e-8 * want.abs(), "{key}: {got} vs {want}");
    }
}
