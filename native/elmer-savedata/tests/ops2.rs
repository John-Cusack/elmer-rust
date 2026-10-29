//! The remaining operators, file options and `SaveLine` options, on the same
//! two-cell mesh with a steady temperature, against the pristine ElmerSolver
//! (`tests/data/ops2/case.sif`, outputs in `tests/data/ops2/golden`).

mod support;

use std::path::{Path, PathBuf};

use elmer_core::SerialComm;
use elmer_solver::{Model, Solver, SolverCtx, SolverError, SolverRegistry};
use support::backends::SimMesh;
use support::harness::*;
use support::meshio::read_elmer_mesh;

/// The steady temperature of the case: 0 at the left, 10 at the right, 6 in
/// the middle column.
struct Replay;

impl Solver for Replay {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let f = ctx.primary_mut().unwrap();
        for (n, w) in [0.0, 6.0, 10.0, 0.0, 6.0, 10.0].iter().enumerate() {
            f.values[n] = *w;
        }
        ctx.nonlinear_change()?;
        Ok(())
    }
}

fn registry() -> SolverRegistry {
    let mut reg = SolverRegistry::new();
    elmer_savedata::register(&mut reg);
    reg.register("Test", "Replay", |_| Ok(Box::new(Replay)));
    reg
}

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

/// `tests/data/ops2/case.sif` as parsed by `elmer-sif`, with the heat solver
/// replaced by [`Replay`] and the files sent to `out`.
fn model(out: &Path) -> Model {
    let mut m = sif_model(&data("ops2/case.sif"));
    m.solvers[0].insert("Procedure", st("\"Test\" \"Replay\""));
    for s in &mut m.solvers[1..] {
        s.insert("Output Directory", dir_str(out));
    }
    m
}

fn run(name: &str) -> (PathBuf, RecLogger) {
    let out = tmp(name);
    let mesh = SimMesh::new(read_elmer_mesh(&data("two_cell"), 2));
    let log = RecLogger::default();
    let mut sim = heat_simulation(model(&out));
    sim.setup(&mesh, &registry(), &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    (out, log)
}

fn numbers(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .map(|l| l.split(|c: char| c.is_whitespace() || c == ',').filter(|t| !t.is_empty()).map(|t| t.parse().unwrap()).collect())
        .collect()
}

fn assert_close(got: &[Vec<f64>], want: &[Vec<f64>], what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: rows");
    for (r, (g, w)) in got.iter().zip(want).enumerate() {
        assert_eq!(g.len(), w.len(), "{what}: columns of row {r}");
        for (c, (a, b)) in g.iter().zip(w).enumerate() {
            let tol = 1e-9 * b.abs().max(1.0);
            assert!((a - b).abs() <= tol, "{what}: row {r} column {}: {a} vs {b}", c + 1);
        }
    }
}

fn golden(name: &str) -> String {
    std::fs::read_to_string(data("ops2/golden").join(name)).unwrap()
}

fn without_path(text: &str) -> String {
    text
        .lines()
        .filter(|l| !l.starts_with("Metadata for") && !l.starts_with("File started at"))
        .collect::<Vec<_>>().join("\n")
}

#[test]
fn every_operator_matches_elmer() {
    let (out, _) = run("ops2-scalars");
    let got = std::fs::read_to_string(out.join("scalars.dat")).unwrap();
    assert_close(&numbers(&got), &numbers(&golden("scalars.dat")), "scalars.dat");
    let names = std::fs::read_to_string(out.join("scalars.dat.names")).unwrap();
    assert_eq!(without_path(&names), without_path(&golden("scalars.dat.names")));
}

#[test]
fn line_marker_append_csv_and_partition_numbering() {
    let (out, _) = run("ops2-second");
    for f in ["second_np1.dat", "second_np1.dat.csv"] {
        let got = std::fs::read_to_string(out.join(f)).unwrap();
        assert_close(&numbers(&got), &numbers(&golden(f)), f);
    }
    // Integer line marker in I10, then the values.
    let dat = std::fs::read_to_string(out.join("second_np1.dat")).unwrap();
    assert_eq!(dat.lines().next().unwrap().len(), 10 + 3 * 22);
    assert_eq!(
        std::fs::read_to_string(out.join("second_np1.dat.marker")).unwrap(),
        golden("second_np1.dat.marker")
    );
    let names = std::fs::read_to_string(out.join("second_np1.dat.names")).unwrap();
    assert_eq!(without_path(&names), without_path(&golden("second_np1.dat.names")));
}

#[test]
fn saveline_without_optimisation_is_first_touch_order() {
    let (out, _) = run("ops2-lines");
    let got = std::fs::read_to_string(out.join("lines.dat")).unwrap();
    assert_eq!(got.lines().count(), 6);
    assert_close(&numbers(&got), &numbers(&golden("lines.dat")), "lines.dat");
    let names = std::fs::read_to_string(out.join("lines.dat.names")).unwrap();
    assert_eq!(without_path(&names), without_path(&golden("lines.dat.names")));
}

#[test]
fn saveline_optimised_order_follows_the_curve() {
    let (out, _) = run("ops2-other");
    let got = std::fs::read_to_string(out.join("other.dat")).unwrap();
    // Left edge and top edge as one polyline: (3,1) (1,1) (0,1) (0,0).
    assert_close(&numbers(&got), &numbers(&golden("other.dat")), "other.dat");
    let names = std::fs::read_to_string(out.join("other.dat.names")).unwrap();
    assert_eq!(without_path(&names), without_path(&golden("other.dat.names")));
}

#[test]
fn reference_values_and_show_norm_name_set_the_norm() {
    let (out, log) = run("ops2-norm");
    let text = log.text();
    assert!(text.contains("SS (ITER=1) (NRM,RELC): (  5.3333333      0.0000000     ) :: second"), "{text}");
    assert!(text.contains("SS (ITER=1) (NRM,RELC): (  1.1250000      0.0000000     ) :: \n"), "{text}");
    let got = std::fs::read_to_string(out.join("ref.dat")).unwrap();
    assert_close(&numbers(&got), &numbers(&golden("ref.dat")), "ref.dat");
}

#[test]
fn cost_function_lands_in_the_simulation_section() {
    let out = tmp("ops2-cost");
    let mesh = SimMesh::new(read_elmer_mesh(&data("two_cell"), 2));
    let log = RecLogger::default();
    let mut sim = heat_simulation(model(&out));
    sim.setup(&mesh, &registry(), &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    use elmer_core::Params;
    assert_eq!(sim.model.simulation.real("cost function"), Some(0.0)); // min: temperature
    assert!((sim.model.simulation.real("max: temperature").unwrap() - 10.0).abs() < 1e-12);
}
