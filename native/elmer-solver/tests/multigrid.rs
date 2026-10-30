//! Geometric multigrid (`Linear System Solver = Multigrid`, `MG Equal
//! Split`, `MG Levels = 4`) on `fem/tests/gmultg` (`tests/data/gmultg/`, its
//! ElmerGrid mesh; `Active Solvers` instead of the legacy `Heat Equation =
//! True`), against the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) run:
//! the solve runs on the mesh split three times with `SplitMeshEqual`
//! (`elmer-mesh-levels`), its parents being the coarse levels. HeatSolve is
//! the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{MeshView, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, split_change, Capture};

#[test]
fn gmultg_top_level_residuals_and_norm() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gmultg");
    let m1 = support::mesh::FileMesh::read(&dir, 2);
    let m2 = elmer_mesh_levels::split_mesh_equal(&m1).unwrap();
    let m3 = elmer_mesh_levels::split_mesh_equal(&m2).unwrap();
    let m4 = elmer_mesh_levels::split_mesh_equal(&m3).unwrap();
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&m4, &reg, &log, &SerialComm).unwrap();
    let coarse: Vec<Box<dyn MeshView>> = vec![Box::new(m1), Box::new(m2), Box::new(m3)];
    sim.set_multigrid_meshes(1, coarse).unwrap();
    sim.run(&m4, &log, &SerialComm).unwrap();

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    // The not-converged notice is level 6, above the sif's Max Output Level 5.
    let got: Vec<String> = log
        .lines(&["GMGSolve: MG Residual", "ComputeChange", "CompareToReferenceSolution: Solver"])
        .into_iter()
        .filter(|l| !l.ends_with("did not converge to tolerance"))
        .collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("ComputeChange") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt));
            assert!((gr - wr).abs() <= 1e-6 * wr, "{g}\n{w}");
        } else if w.starts_with("GMGSolve") {
            assert_eq!(g[..g.find("is:").unwrap()], w[..w.find("is:").unwrap()]);
            // Levels 3 and 4 to 1e-6 (the Galerkin products are summed in
            // another order); level 2 sits at round-off after the direct
            // solve below it.
            if !w.contains("level: 2 ") {
                let num = |l: &str| -> f64 { l.split_whitespace().last().unwrap().parse().unwrap() };
                assert!((num(g) - num(w)).abs() <= 1e-6 * num(w), "{g}\n{w}");
            }
        } else {
            assert!(g.contains("PASSED"), "{g}");
        }
    }
}
