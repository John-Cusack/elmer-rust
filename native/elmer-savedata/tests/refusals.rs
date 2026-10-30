//! Everything SaveData can do that this crate does not is refused with a
//! message naming it: when the solver is built for keywords and operators,
//! at the visit for the run-time conditions (no element library, more than
//! one partition).

mod support;

use std::path::PathBuf;

use elmer_core::{Communicator, ParamMap, SerialComm, Value};
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
        ("convective flux", "flux reconstruction"),
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
        ("Save Eigenvalues", yes()),
        ("Save Eigenfrequencies", yes()),
        ("Simulation Timing", yes()),
        ("Save Component Results", yes()),
        ("Json Format", yes()),
        ("Filename Numbering", yes()),
        ("Run Control Numbering", yes()),
        ("Mark Failed Strategy", yes()),
        ("Avoid Failed Strategy", yes()),
        ("Coefficient 2", st("Density")),
        ("Expression 1", st("1+1")),
    ] {
        let mut s = with_oper("max");
        s.insert(key, value);
        let m = refused("SaveScalars", s);
        assert!(m.to_lowercase().contains(&key.to_lowercase()) || m.contains("coefficient") , "{key}: {m}");
    }
    // A keyword switched off is not a request.
    let mut s = with_oper("max");
    s.insert("Json Format", Value::Logical(false));
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
        ("Calculate Weights", yes()),
        ("Filename Numbering", yes()),
        ("Show Norm Index", ints(&[1])),
        ("Parallel Reduce", yes()),
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

struct TwoPartitions;

impl Communicator for TwoPartitions {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        2
    }
    fn sum_f64(&self, _: &mut [f64]) {}
    fn max_f64(&self, _: &mut [f64]) {}
    fn barrier(&self) {}
}

#[test]
fn mpi_reductions_are_refused_at_the_first_visit() {
    let mut c = Case::bare(mesh(), mesh_model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let a = c.hook("SaveScalars", with_oper("max")).unwrap();
    let b = c.hook("SaveLine", sec(&[("Variable 1", st("temperature"))])).unwrap();
    for ix in [a, b] {
        match c.visit_with(ix, &TwoPartitions) {
            Err(SolverError::Unsupported(m)) => assert!(m.contains("MPI"), "{m}"),
            other => panic!("{other:?}"),
        }
    }
    let _ = SerialComm;
}
