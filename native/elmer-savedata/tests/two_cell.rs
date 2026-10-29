//! Heat -> SaveScalars -> SaveScalars -> SaveLine through `Simulation::run`
//! on a two-cell transient case, against the files the pristine ElmerSolver
//! (26.2, Rev 9f6af2f85) wrote for the same `tests/data/two_cell/case.sif`.
//!
//! The reference run is serial; its numbers differ from ours only by the
//! linear solver's rounding, so values are compared to 1e-9 relative and
//! everything else (names, labels, row order) exactly.

mod support;

use std::path::{Path, PathBuf};

use elmer_core::{ParamMap, Params, SerialComm, Value};
use elmer_solver::{Model, Solver, SolverCtx, SolverError, SolverRegistry};
use support::backends::SimMesh;
use support::harness::*;
use support::meshio::read_elmer_mesh;

/// Stands in for HeatSolve: writes the temperatures the pristine run
/// printed (`golden/line.dat`) into the primary variable at each step, then
/// runs the nonlinear `ComputeChange` so that `Var % Norm` is Elmer's.
///
/// Elmer's first step of this case is reproduced by a backward-Euler step
/// that starts with the Dirichlet values already in the initial field, which
/// `elmer-solver` does not do; the heat solution of this transient case
/// therefore differs (see `tests/data/README.md`). The SaveData hooks are
/// compared against Elmer's own fields so that this does not matter here.
struct Replay;

impl Solver for Replay {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let step = ctx.scalar("timestep").unwrap() as usize;
        let t = [18.0 / 11.0, 342.0 / 121.0][step - 1];
        let want = [0.0, t, 10.0, 0.0, t, 10.0];
        let f = ctx.primary_mut().unwrap();
        for (n, w) in want.iter().enumerate() {
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
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/two_cell").join(name)
}

/// `tests/data/two_cell/case.sif` as parsed by `elmer-sif`, with the heat
/// solver replaced by [`Replay`] and the files sent to `out`.
fn model(out: &Path) -> Model {
    let mut m = sif_model(&data("case.sif"));
    // `Save Coordinates(1,2)` comes out of elmer-sif as 2x2 with two values
    // (it uses n1 * n2 for the row count); give the intended 1x2 shape.
    m.solvers[1].insert("Save Coordinates", array(1, 2, &[2.9, 0.9]));
    m.solvers[0].insert("Procedure", st("\"Test\" \"Replay\""));
    for s in &mut m.solvers[1..] {
        s.insert("Output Directory", dir_str(out));
    }
    m
}

/// Runs the case; returns the output directory and the recorded messages.
fn run(name: &str) -> (PathBuf, RecLogger, elmer_solver::Simulation) {
    let out = tmp(name);
    let mesh = SimMesh::new(read_elmer_mesh(&data(""), 2));
    let log = RecLogger::default();
    let mut sim = heat_simulation(model(&out));
    sim.setup(&mesh, &registry(), &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    (out, log, sim)
}

fn numbers(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
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

/// A `.names` file without the path line, whose directory differs per run.
fn names(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("Metadata for"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn scalars_file_matches_elmer() {
    let (out, _, _) = run("scalars");
    let got = std::fs::read_to_string(out.join("scalars.dat")).unwrap();
    let want = std::fs::read_to_string(data("golden/scalars.dat")).unwrap();
    // Values with rounding of the linear solver.
    assert_close(&numbers(&got), &numbers(&want), "scalars.dat");
    // Same column widths: two steps, integer columns formatted as integers.
    assert_eq!(got.lines().count(), 2);
    assert_eq!(got.lines().next().unwrap().len(), want.lines().next().unwrap().len());
    assert_eq!(
        names(&out.join("scalars.dat.names")),
        names(&data("golden/scalars.dat.names"))
    );
}

#[test]
fn names_file_header_names_the_file() {
    let (out, _, _) = run("scalars-names");
    let text = std::fs::read_to_string(out.join("scalars.dat.names")).unwrap();
    assert!(text.contains(&format!("Metadata for SaveScalars file: {}/scalars.dat", out.display())));
}

#[test]
fn target_variable_is_visible_to_the_next_hook() {
    let (out, _, sim) = run("reader");
    let got = std::fs::read_to_string(out.join("reader.dat")).unwrap();
    let want = std::fs::read_to_string(data("golden/reader.dat")).unwrap();
    assert_close(&numbers(&got), &numbers(&want), "reader.dat");
    assert_eq!(
        names(&out.join("reader.dat.names")),
        names(&data("golden/reader.dat.names"))
    );
    // The exported variable holds the integral of the last step, and the
    // reader's own norm variable carries it as `Show Norm Index` says.
    let v = sim.vars.get("temperature integral").unwrap();
    assert!((v.values[0] - 14.23966942149).abs() < 1e-8);
    let n = sim.vars.get("savereader_var").unwrap();
    assert!((n.conv.norm - 14.23966942149).abs() < 1e-8);
    assert!(!n.output);
}

#[test]
fn pseudo_norm_line_has_elmers_format() {
    let (_, log, _) = run("pseudo-norm");
    let text = log.text();
    let lines: Vec<&str> = text.lines().filter(|l| l.contains(":: savereader")).collect();
    assert_eq!(lines.len(), 2, "{text}");
    assert_eq!(
        lines[0],
        "ComputeChange: SS (ITER=1) (NRM,RELC): (  12.454545      0.0000000     ) :: savereader"
    );
    assert_eq!(
        lines[1],
        "ComputeChange: SS (ITER=1) (NRM,RELC): (  14.239669      0.0000000     ) :: savereader"
    );
}

#[test]
fn line_file_matches_elmer_row_for_row() {
    let (out, _, _) = run("line");
    let got = std::fs::read_to_string(out.join("line.dat")).unwrap();
    let want = std::fs::read_to_string(data("golden/line.dat")).unwrap();
    // Order, labels (timestep, call count, boundary condition, node) and
    // values: the RCM-ordered boundary nodes come out as in Elmer.
    assert_close(&numbers(&got), &numbers(&want), "line.dat");
    let labels = |t: &str| -> Vec<String> {
        t.lines()
            .map(|l| l.split_whitespace().take(4).collect::<Vec<_>>().join(" "))
            .collect()
    };
    assert_eq!(labels(&got), labels(&want));
    assert_eq!(
        names(&out.join("line.dat.names")),
        names(&data("golden/line.dat.names"))
    );
}

#[test]
fn simulation_section_keeps_the_values_for_later_solvers() {
    let (_, _, sim) = run("simulation");
    let max = sim.model.simulation.get("max: temperature").expect("max: temperature stored");
    let Value::Real { data, .. } = max else { panic!("{max:?}") };
    assert!((data[0] - 10.0).abs() < 1e-9);
    let _: &ParamMap = &sim.model.simulation;
}
