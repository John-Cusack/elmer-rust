//! `StructuredProjectToPlane` (and the `StructuredMeshMapper` with a mid
//! surface in front of it) on the pristine quick test `structmap5`: a layered
//! 3D mesh of two bodies. The mapped nodes, the projected fields and the
//! HeatSolve on the mapped mesh against the pristine run's VTU file, and the
//! norms against the `.sif`'s `Reference Norm`s.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::meshio::read_elmer_mesh;
use support::run::*;
use support::simmesh::SimMesh;

const OPERATORS: [&str; 11] = [
    "depth",
    "height",
    "distance",
    "index",
    "upper depth",
    "upper height",
    "upper distance",
    "upper index",
    "upper thickness",
    "int",
    "upper int",
];

fn solver3() -> elmer_core::ParamMap {
    // `Solver 3` appears twice in the `.sif` (HeatB, then HeightDepth); the
    // pristine run merges them: the later `Equation` and `Procedure` win and
    // the first block's `Variable = TempB` stays.
    let mut v = vec![
        ("Equation", st("HeightDepth")),
        ("Variable", st("TempB")),
        (
            "Procedure",
            st("\"StructuredProjectToPlane\" \"StructuredProjectToPlane\""),
        ),
        ("Active Coordinate", ints(&[3])),
        ("Show Norm Index", ints(&[11])),
        ("Reference Norm", Value::real(0.59397225)),
    ];
    let ops: Vec<(String, Value)> = OPERATORS
        .iter()
        .enumerate()
        .map(|(i, o)| (format!("Operator {}", i + 1), st(o)))
        .collect();
    let mut s = sec(&v);
    for (k, val) in ops {
        s.insert(&k, val);
    }
    s.insert("Variable 10", st("Coordinate 1"));
    s.insert("Variable 11", st("Coordinate 1"));
    v.clear();
    s
}

fn model(with_heat: bool) -> Model {
    let mapper = sec(&[
        ("Equation", st("MapCoordinates")),
        (
            "Procedure",
            st("\"StructuredMeshMapper\" \"StructuredMeshMapper\""),
        ),
        ("Active Coordinate", ints(&[3])),
    ]);
    let heat_a = sec(&[
        ("Equation", st("HeatA")),
        ("Variable", st("TempA")),
        ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
        ("Linear System Solver", st("Iterative")),
        ("Linear System Iterative Method", st("BiCGStab")),
        ("Linear System Max Iterations", ints(&[1000])),
        ("Linear System Convergence Tolerance", Value::real(1.0e-12)),
        ("Linear System Preconditioning", st("ILU0")),
        ("Linear System Residual Output", ints(&[0])),
        ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
        ("Stabilize", lg(true)),
        (
            "Nonlinear System Convergence Tolerance",
            Value::real(1.0e-5),
        ),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        (
            "Nonlinear System Newton After Tolerance",
            Value::real(1.0e-2),
        ),
        ("Reference Norm", Value::real(1.5335462336552876)),
    ]);
    let (solvers, active) = if with_heat {
        (vec![mapper, heat_a, solver3()], (ints(&[2]), ints(&[3])))
    } else {
        (
            vec![mapper, sec(&[("Equation", st("none"))]), solver3()],
            (ints(&[1]), ints(&[3])),
        )
    };
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ]),
        solvers,
        equations: vec![
            sec(&[("Active Solvers", active.0)]),
            sec(&[("Active Solvers", active.1)]),
        ],
        bodies: vec![
            sec(&[
                ("Equation", ints(&[1])),
                ("Material", ints(&[1])),
                ("Body Force", ints(&[1])),
            ]),
            sec(&[
                ("Equation", ints(&[2])),
                ("Material", ints(&[2])),
                ("Body Force", ints(&[1])),
            ]),
        ],
        body_forces: vec![sec(&[("Heat Source", Value::real(0.0))])],
        materials: vec![
            sec(&[
                ("Density", Value::real(1.0)),
                ("Heat Conductivity", Value::real(1.0)),
            ]),
            sec(&[
                ("Density", Value::real(1.0)),
                ("Heat Conductivity", Value::real(1.0)),
            ]),
        ],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[5])),
                ("TempA", Value::real(1.0)),
                ("TempB", Value::real(1.0)),
                ("Bottom Surface", Value::real(0.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[6])),
                ("TempA", Value::real(2.0)),
                ("TempB", Value::real(2.0)),
                (
                    "Mid Surface",
                    matc(&["coordinate"], "1-0.1*cos(pi*tx(1))*sin(pi*tx(0))"),
                ),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[7])),
                ("TempA", Value::real(3.0)),
                ("TempB", Value::real(3.0)),
                (
                    "Top Surface",
                    matc(&["coordinate"], "2+0.1*cos(pi*tx(1))*sin(pi*tx(0))"),
                ),
            ]),
        ],
        ..Model::default()
    }
}

fn mesh() -> SimMesh {
    SimMesh::new(read_elmer_mesh(&data("structmap5"), 3))
}

#[test]
fn mapped_nodes_and_projections_match_pristine() {
    let mesh = mesh();
    let n = mesh.mesh.coords.len();
    let sim = run(model(false), &mesh).unwrap();
    let r = Reference::read("structmap5.txt");
    for (k, c) in ["x", "y", "z"].into_iter().enumerate() {
        let d = max_diff(&coords_of(&mesh, k), &r.col(c));
        assert!(d < 1e-12, "{c}: max |diff| {d:e}");
    }
    for col in [
        "depth",
        "height",
        "distance",
        "index",
        "upper depth",
        "upper height",
        "upper distance",
        "upper index",
        "upper thickness",
        "int coordinate 1",
        "upper int coordinate 1",
    ] {
        let d = max_diff(&nodal(&sim, col, n), &r.col(col));
        assert!(d < 1e-12, "{col}: max |diff| {d:e}");
    }
}

#[test]
fn heat_and_the_shown_norm_reach_the_reference_norms() {
    let mesh = mesh();
    let n = mesh.mesh.coords.len();
    let sim = run(model(true), &mesh).unwrap();
    let r = Reference::read("structmap5.txt");
    let d = max_diff(&nodal(&sim, "tempa", n), &r.col("tempa"));
    assert!(d < 1e-8, "tempa max |diff| {d:e}");
    let norm = sim.vars.get("tempa").unwrap().conv.norm;
    assert!(((norm - 1.5335462336552876) / norm).abs() < 1e-8, "{norm}");
    // Solver 3's variable is filled with the norm of the 11th operator's field.
    let shown = sim.vars.get("tempb").unwrap().conv.norm;
    assert!(((shown - 0.59397225) / shown).abs() < 1e-7, "{shown}");
}
