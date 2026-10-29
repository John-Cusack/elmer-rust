//! `MeshChecksum` on the mesh of the pristine quick test `ElmerGridConvertUnv`
//! (`box.unv` converted by ElmerGrid, 554 nodes, 2308 tetrahedra, 570
//! boundary triangles), against the pristine run's own printout:
//!
//! ```text
//! ThisResults:  5170  1163232  2308  3059754  172710  94  247056  166390.56857056363
//! PseudoNorm:   1.0000000000000000
//! ```

mod support;

use std::path::PathBuf;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_meshtools::{checksum, pseudo_norm, register};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn mesh() -> SimMesh {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/unv");
    SimMesh::new(read_elmer_mesh(&dir, 3))
}

/// The `.sif` of the test: one body, one boundary condition on tag 1, and the
/// solver's `Reference Values(8)`.
fn model(reference: [f64; 8]) -> Model {
    Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("steady state".into())),
            ("Steady State Max Iterations", Value::Integer(vec![1])),
            ("Output Intervals", Value::Integer(vec![0])),
        ]),
        solvers: vec![sec(&[
            ("Equation", Value::Str("CheckMesh".into())),
            (
                "Procedure",
                Value::Str("\"MeshChecksum\" \"MeshChecksum\"".into()),
            ),
            (
                "Reference Values",
                Value::Real {
                    rows: 8,
                    cols: 1,
                    data: reference.to_vec(),
                },
            ),
        ])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![sec(&[
            ("Equation", Value::Integer(vec![1])),
            ("Target Bodies", Value::Integer(vec![1])),
        ])],
        boundary_conditions: vec![sec(&[("Target Boundaries", Value::Integer(vec![1]))])],
        ..Model::default()
    }
}

const PRISTINE: [f64; 8] = [
    5170.0,
    1163232.0,
    2308.0,
    3059754.0,
    172710.0,
    94.0,
    247056.0,
    166390.56857056363,
];

#[test]
fn sums_match_the_pristine_printout() {
    let m = mesh();
    let cs = checksum(&m, &model(PRISTINE));
    assert_eq!(cs[..7], PRISTINE[..7]);
    assert_eq!(cs[7], PRISTINE[7], "coordinate sum to the last bit");
    assert_eq!(pseudo_norm(&cs, &PRISTINE), 1.0);
}

#[test]
fn solver_leaves_a_unit_norm() {
    let m = mesh();
    let mut sim = Simulation::new(model(PRISTINE));
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    sim.setup(&m, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&m, &NullLogger, &SerialComm).unwrap();
    let id = sim
        .slot(1)
        .unwrap()
        .variable
        .expect("MeshChecksum_init adds meshcheck_var");
    let v = sim.vars.by_id(id);
    assert_eq!(v.name, "meshcheck_var");
    assert_eq!(v.conv.norm, 1.0);
    assert!(!v.output);
    // ListAddNewConstReal(Params, 'Reference Norm', 1).
    assert_eq!(
        sim.slot(1).unwrap().params.0["reference norm"],
        Value::real(1.0)
    );
}

#[test]
fn a_different_mesh_gives_a_norm_above_one() {
    let m = mesh();
    let mut wrong = PRISTINE;
    wrong[1] *= 2.0;
    wrong[5] = 0.0;
    let cs = checksum(&m, &model(PRISTINE));
    // slot 2: ratio 1/2 -> 2; slot 6: zero reference -> 1 + 94.
    let want = (7.0 + 2.0 + 0.0 + 95.0 - 1.0) / 8.0;
    assert!(
        (pseudo_norm(&cs, &wrong) - want).abs() < 1e-12,
        "{}",
        pseudo_norm(&cs, &wrong)
    );
}
