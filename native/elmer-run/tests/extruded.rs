//! MinimalExample extruded into a 3D wedge column (`Extruded Mesh Levels`).
//! Its own test binary: the physics crates keep one element basis per
//! process, and this case needs the 3D one.

mod support;

use elmer_core::{SerialComm, Value};
use elmer_run::{run_case, Case};
use support::case::{minimal_example, scratch_case};
use support::logger;
use support::mesh::FileMesh;

#[test]
fn extruded_mesh_runs_on_the_wedge_column() {
    let dir = scratch_case("minimal_example", "extruded");
    let mesh = FileMesh::read(&dir.join("triangle"), 2);
    let mut model = minimal_example();
    model.solvers[0].0.remove("reference norm");
    model.simulation.insert("Ascii Output", Value::Logical(true));
    model.simulation.insert("Extruded Mesh Levels", Value::Integer(vec![3]));
    model.simulation.insert("Coordinate System", Value::Str("cartesian 3d".into()));
    let (log, buf) = logger();
    let case = Case { model, reloads: vec![], mesh: &mesh, case_dir: dir.clone(), mesh_dir: dir.join("triangle") };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{e}\n{}", buf.text()));
    // Three node layers of the triangle: 9 nodes, 2 wedges (706).
    let vtu = std::fs::read_to_string(&out.written[0]).unwrap();
    assert!(vtu.contains("NumberOfPoints=\"9\""), "{vtu}");
    assert!(vtu.contains("NumberOfCells="), "{vtu}");
}
