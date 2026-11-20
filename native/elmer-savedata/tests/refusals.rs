//! Everything SaveData can do that this crate does not is refused with a
//! message naming it: when the solver is built for keywords and operators,
//! at the visit for the run-time conditions (no element library, more than
//! one partition).

mod support;

use std::path::PathBuf;

use elmer_core::{Communicator, ParamMap, Params, SerialComm, Value};
use elmer_solver::{Model, SolverError, SolverRegistry};
use support::harness::*;
use support::meshio::{read_elmer_mesh, Mesh};

fn mesh() -> Mesh {
    read_elmer_mesh(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/two_cell"), 2)
}

fn build(entry: &str, section: ParamMap) -> Result<(), SolverError> {
    let mut reg = SolverRegistry::new();
    elmer_savedata::register(&mut reg);
    reg.get("SaveData", entry).unwrap()(&section).map(|_| ())
}

fn refused(entry: &str, section: ParamMap) -> String {
    match build(entry, section) {
        Err(SolverError::Unsupported(m)) => m,
        Err(e) => panic!("expected Unsupported, got {e}"),
        Ok(()) => panic!("accepted"),
    }
}

fn with_oper(op: &str) -> ParamMap {
    sec(&[("Variable 1", st("temperature")), ("Operator 1", st(op))])
}

#[test]
fn operators_that_are_not_implemented_are_refused_by_name() {
    for (op, why) in [
        ("convective energy", "convective energy"),
        ("cpu time", "cpu time"),
        ("wall time", "wall time"),
        ("cpu memory", "cpu memory"),
        ("threads", "threads"),
        ("positive max", "positive/negative"),
        ("negative sum", "positive/negative"),
        ("frobnicate", "frobnicate"),
        ("body boundary sum", "mask prefix"),
        ("material nodes", "mask prefix"),
    ] {
        let m = refused("SaveScalars", with_oper(op));
        assert!(m.contains(why) && m.contains("SaveScalars"), "{op}: {m}");
    }
}

#[test]
fn keywords_that_are_not_implemented_are_refused_by_name() {
    for (key, value) in [
        ("Polyline Coordinates", Value::Real { rows: 2, cols: 2, data: vec![0.0, 0.0, 1.0, 1.0] }),
        ("Exact Coordinates", yes()),
        ("Simulation Timing", yes()),
        ("Json Format", yes()),
        ("Filename Numbering", yes()),
        ("Run Control Numbering", yes()),
        ("Mark Failed Strategy", yes()),
        ("Avoid Failed Strategy", yes()),
        ("Expression 1", st("1+1")),
    ] {
        let mut s = with_oper("max");
        s.insert(key, value);
        let m = refused("SaveScalars", s);
        assert!(m.to_lowercase().contains(&key.to_lowercase()), "{key}: {m}");
    }
    // A keyword switched off is not a request.
    let mut s = with_oper("max");
    s.insert("Json Format", Value::Logical(false));
    build("SaveScalars", s).unwrap();
    // `Coefficient n` weights the integrals; a statistic ignores it.
    let mut s = with_oper("max");
    s.insert("Coefficient 1", st("Density"));
    build("SaveScalars", s).unwrap();
}

#[test]
fn saveline_modes_other_than_existing_boundaries_are_refused() {
    for (key, value) in [
        ("Polyline Coordinates", Value::Real { rows: 2, cols: 2, data: vec![0.0, 0.0, 1.0, 1.0] }),
        ("Circle Coordinates", Value::Real { rows: 2, cols: 2, data: vec![0.0, 0.0, 1.0, 1.0] }),
        ("Save Isocurves", yes()),
        ("Save Axis", yes()),
        ("Save Heat Flux", yes()),
        ("Save Flux", yes()),
        ("Intersect Edge", yes()),
        ("Break Line Loop", yes()),
        ("Filename Numbering", yes()),
    ] {
        let mut s = sec(&[("Variable 1", st("temperature"))]);
        s.insert(key, value);
        let m = refused("SaveLine", s);
        assert!(m.contains("SaveLine"), "{key}: {m}");
    }
}

fn mesh_model() -> Model {
    Model {
        bodies: vec![sec(&[("Target Bodies", ints(&[1]))]), sec(&[("Target Bodies", ints(&[2]))])],
        boundary_conditions: vec![sec(&[("Target Boundaries", ints(&[1])), ("Save Line", yes())])],
        ..Model::default()
    }
}

#[test]
fn integrals_without_an_element_library_are_not_wired() {
    // This binary never calls `install_fe`.
    let mut c = Case::bare(mesh(), mesh_model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c.hook("SaveScalars", with_oper("int")).unwrap();
    let err = c.visit(ix, false).unwrap_err();
    assert!(matches!(err, SolverError::NotWired(_)), "{err}");
    // Statistics need none.
    let ix = c.hook("SaveScalars", with_oper("max")).unwrap();
    c.visit(ix, false).unwrap();
}

/// Rank 0 of two; the other rank holds 5 for every `max_f64` (so `-5` for
/// a minimum, which reduces as the maximum of the negated values).
struct TwoPartitions;

impl Communicator for TwoPartitions {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        2
    }
    fn sum_f64(&self, _: &mut [f64]) {}
    fn max_f64(&self, buf: &mut [f64]) {
        buf.iter_mut().for_each(|x| *x = x.max(5.0));
    }
    fn barrier(&self) {}
}

#[test]
fn mpi_runs_need_parallel_reduce_and_reducible_operators() {
    let mut c = Case::bare(mesh(), mesh_model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    // Per-partition files (a Filename without Parallel Reduce) and SaveLine
    // stay refused.
    let mut f = with_oper("max");
    f.insert("Filename", st("scalars.dat"));
    let a = c.hook("SaveScalars", f).unwrap();
    let b = c.hook("SaveLine", sec(&[("Variable 1", st("temperature"))])).unwrap();
    for (ix, why) in [(a, "Parallel Reduce"), (b, "MPI")] {
        match c.visit_with(ix, &TwoPartitions) {
            Err(SolverError::Unsupported(m)) => assert!(m.contains(why), "{m}"),
            other => panic!("{other:?}"),
        }
    }
    // An operator whose reduction is not implemented.
    let mut p = with_oper("int");
    p.insert("Parallel Reduce", yes());
    let ix = c.hook("SaveScalars", p).unwrap();
    match c.visit_with(ix, &TwoPartitions) {
        Err(SolverError::Unsupported(m)) => assert!(m.contains("in parallel"), "{m}"),
        other => panic!("{other:?}"),
    }
    let _ = SerialComm;
}

/// Without `Parallel Reduce` (and without a file) each rank keeps its own
/// value, unreduced and under the plain name (SaveScalars.F90:283-325).
#[test]
fn without_parallel_reduce_the_values_stay_local() {
    let mut c = Case::bare(mesh(), mesh_model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c.hook("SaveScalars", with_oper("max")).unwrap();
    c.visit_with(ix, &TwoPartitions).unwrap();
    let s = &c.sim.model.simulation;
    let local = s.real("max: temperature").unwrap();
    assert!(local < 5.0, "{local}");
    assert_eq!(s.real("max: temperature : mpi_max"), None);
}

/// AddToSaveList: with `Parallel Reduce` the extrema are reduced over the
/// ranks and their names get ` : mpi_<op>`.
#[test]
fn parallel_reduce_takes_the_extrema_over_the_ranks() {
    let mut c = Case::bare(mesh(), mesh_model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c
        .hook(
            "SaveScalars",
            sec(&[
                ("Variable 1", st("temperature")),
                ("Operator 1", st("max")),
                ("Variable 2", st("temperature")),
                ("Operator 2", st("min")),
                ("Parallel Reduce", yes()),
            ]),
        )
        .unwrap();
    c.visit_with(ix, &TwoPartitions).unwrap();
    let s = &c.sim.model.simulation;
    assert_eq!(s.real("max: temperature : mpi_max"), Some(5.0));
    assert_eq!(s.real("min: temperature : mpi_min"), Some(-5.0));
}

/// Polylines are kept only with `Polyline Divisions`; with them the points
/// along the line are located and interpolated, and `Show Norm Index` sums
/// a column into the solver's norm.
#[test]
fn polyline_points_interpolate_and_feed_the_norm() {
    let mut s = sec(&[("Variable 1", st("temperature"))]);
    s.insert("Polyline Coordinates", Value::Real { rows: 2, cols: 2, data: vec![0.0, 0.0, 1.0, 1.0] });
    assert!(refused("SaveLine", s.clone()).contains("Polyline Divisions"));

    let mut c = Case::bare(mesh(), mesh_model());
    // Bilinear, so the quads interpolate it exactly.
    c.field("Temperature", 1, |x| vec![x[0] + 2.0 * x[1]]);
    // SaveLine writes below the working directory.
    let out = PathBuf::from(format!("polyline-test-{}.dat", std::process::id()));
    let ix = c
        .hook(
            "SaveLine",
            sec(&[
                ("Equation", st("line")),
                ("Filename", st(out.to_str().unwrap())),
                ("Variable 1", st("temperature")),
                ("Polyline Coordinates", Value::Real { rows: 2, cols: 2, data: vec![0.0, 0.5, 3.0, 0.5] }),
                ("Polyline Divisions", ints(&[6])),
                // Columns: call count, line, point, x, y, z, temperature.
                ("Show Norm Index", ints(&[7])),
            ]),
        )
        .unwrap();
    c.visit(ix, false).unwrap();
    let text = std::fs::read_to_string(&out).unwrap();
    let _ = std::fs::remove_file(&out);
    let _ = std::fs::remove_file(format!("{}.names", out.display()));
    // The model's `Save Line` boundary (y = 0, x = 3, 1, 0) comes first;
    // then the seven points x = 0, 0.5, .., 3 at y = 0.5, in order. The norm
    // sums the temperature column of every row: 4 + sum of (x + 1) = 21.5.
    let rows: Vec<&str> = text.lines().collect();
    assert_eq!(rows.len(), 10, "{text}");
    assert!(rows[3].starts_with("1 1 0 "), "{}", rows[3]);
    assert!(rows[9].starts_with("1 1 6 "), "{}", rows[9]);
    let norm = c.sim.vars.get("line_var").unwrap().conv.norm;
    assert!((norm - 21.5).abs() < 1e-12, "{norm}");
}
