//! `Block Solver = True` without a `Procedure`, checked against the
//! pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) runs of
//! `fem/tests/BlockPoisson1` (one block, solved in the standard manner) and
//! `BlockPoisson2` (two blocks, `Block Preconditioner = False`: block
//! Jacobi) from their `case.sif` files on their ElmerGrid meshes, at
//! `Max Output Level = 10`: the BlockSolver and BlockStandardIter messages,
//! `ComputeChange` and the reference norm (the block norms to their four
//! printed digits). The `PoissonAssembly.f90` kernel is
//! `elmer-fixtures`' port; basis and quadrature are the test stand-ins.

mod support;

use std::path::Path;
use std::sync::Arc;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, SolverRegistry};
use support::fixtures::Bulk;
use support::{sif_model, split_change, Capture};

/// The messages up to `Max Output Level = 10`, as the pristine runs.
struct Level10(Capture);

impl elmer_core::Logger for Level10 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 10 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

fn check(mesh_dir: &str, sif: &str, stdout: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join(mesh_dir), 2);
    let model = sif_model(&dir.join(sif));
    let mut reg = SolverRegistry::new();
    reg.register_assembly("PoissonAssembly", "BulkAssembly", Arc::new(|_| Ok(Box::new(Bulk::default()))));
    let log = Level10(Capture::default());
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();

    let want = std::fs::read_to_string(dir.join(stdout)).unwrap();
    let want: Vec<&str> = want.lines().collect();
    let got = log
        .0
        .lines(&["BlockSolver", "ComputeChange", "CompareToReferenceSolution: Solver", "BlockStandardIter"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "line {i}");
            let tol = if wr < 1e-6 { 1e-8 } else { 1e-6 * wr };
            assert!((gr - wr).abs() <= tol, "line {i}\n{g}\n{w}");
        } else if w.contains(" norms:") {
            // elmer-block prints these with Rust's exponent form (3.456E-2,
            // not ES 3.456E-02): compare the four digits, and round-off
            // sized changes only as such.
            assert_eq!(g.split("norms:").next(), w.split("norms:").next(), "line {i}");
            let nums = |l: &str| -> Vec<f64> {
                l.split("norms:").nth(1).unwrap().split_whitespace().map(|t| t.parse().unwrap()).collect()
            };
            for (a, b) in nums(g).into_iter().zip(nums(w)) {
                assert!((a < 1e-10 && b < 1e-10) || (a - b).abs() <= 1e-3 * b, "line {i}\n{g}\n{w}");
            }
        } else {
            assert_eq!(g, w, "line {i}");
        }
    }
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
}

#[test]
fn one_block_in_the_standard_manner() {
    check("angle", "blockpoisson/bp1.sif", "blockpoisson/bp1_stdout.txt");
}

#[test]
fn two_blocks_by_block_jacobi() {
    check("blockpoisson/bp2", "blockpoisson/bp2/case.sif", "blockpoisson/bp2/stdout.txt");
}
