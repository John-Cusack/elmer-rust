//! Point Dirichlet conditions (`Target Nodes`, `Target Coordinates`),
//! checked against a pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run of
//! `tests/data/pointbcs/case.sif` (HeatSolve on the LinearSystemNamespace
//! angle mesh, `tests/data/namespace/`, with two target nodes at 0.05 and two
//! target coordinates at `0.1 x`): the `ComputeChange` lines, and the values
//! at the points. HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::{MeshView, SerialComm};
use elmer_solver::{Simulation, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn target_nodes_and_coordinates() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("namespace"), 2);
    let model = sif_model(&data.join("pointbcs/case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(data.join("pointbcs/stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
    }

    let t = sim.vars.get("temperature").unwrap();
    let at = |node: u32| t.values[t.perm.as_deref().unwrap()[node as usize] as usize];
    for node in [39, 76] {
        assert!((at(node) - 0.05).abs() < 1e-14, "node {}", node + 1);
    }
    for c in [[0.3, 0.4], [0.6, 0.25]] {
        let node = (0..mesh.n_nodes() as u32)
            .min_by(|&a, &b| {
                let d = |n: u32| (mesh.coord(n)[0] - c[0]).powi(2) + (mesh.coord(n)[1] - c[1]).powi(2);
                d(a).total_cmp(&d(b))
            })
            .unwrap();
        assert!((at(node) - 0.1 * mesh.coord(node)[0]).abs() < 1e-12, "near {c:?}");
    }
}
