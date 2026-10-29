//! `HeatSolver` end to end through `elmer-solver`'s `Simulation` (see
//! `support::cases`), with the pristine ElmerSolver's solution as reference.

mod support;

use elmer_core::{NullLogger, SerialComm, Value};
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{data, direct, heat_solver, ints, max_diff, nodal, run, sec, simulation, st};
use support::glue::{read_perm, read_vector};
use support::meshio::read_elmer_mesh;

fn solver_section(extra: &[(&str, Value)]) -> elmer_core::ParamMap {
    let mut v = direct();
    v.extend(extra.iter().cloned());
    heat_solver("Temperature", &v)
}

/// Temperature by node.
fn temperature(sim: &elmer_solver::Simulation, n_nodes: usize) -> Vec<f64> {
    nodal(sim, "temperature", n_nodes)
}

/// Fortran's solution dump, by node.
fn reference(case: &str) -> Vec<f64> {
    let perm = read_perm(&data(case).join("linsys_perm.dat"));
    let sol = read_vector(&data(case).join("linsys_sol.dat"));
    perm.iter().map(|&p| sol[p]).collect()
}

fn minimal_model() -> Model {
    Model {
        simulation: sec(&[("Simulation Type", st("steady state"))]),
        solvers: vec![solver_section(&[])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
        ])],
        materials: vec![sec(&[
            ("Heat Conductivity", Value::real(1.0)),
            ("Density", Value::real(1.0)),
        ])],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Temperature", Value::real(0.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[3])),
                ("Heat Flux", Value::real(1.0)),
            ]),
        ],
        ..Model::default()
    }
}

#[test]
fn minimal_example_solution_matches_elmer() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("minimal_example"), 2));
    let sim = run(minimal_model(), &mesh).unwrap();
    let t = temperature(&sim, 3);
    let want = reference("minimal_example");
    assert!(max_diff(&t, &want) < 1e-5, "{t:?} vs {want:?}");
    assert!((t[0] - 0.5).abs() < 1e-13);
    // The sif's own reference norm.
    let norm = (t.iter().map(|x| x * x).sum::<f64>() / 3.0).sqrt();
    assert!((norm - 2.88675135e-1).abs() < 1e-8, "{norm}");
}

#[test]
fn thermal_pipe_2d_solution_matches_elmer() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("thermal_pipe2d"), 2));
    let n = mesh.mesh.coords.len();
    let model = Model {
        simulation: sec(&[("Simulation Type", st("steady state"))]),
        solvers: vec![solver_section(&[])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: (1..=3)
            .map(|b| {
                sec(&[
                    ("Target Bodies", ints(&[b])),
                    ("Equation", ints(&[1])),
                    ("Material", ints(&[1])),
                ])
            })
            .collect(),
        materials: vec![sec(&[
            ("Heat Conductivity", Value::real(1.0)),
            ("Density", Value::real(1000.0)),
        ])],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[9, 10, 11, 12])),
                ("Temperature", Value::real(10.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[5, 6, 7, 8]))]),
            sec(&[
                ("Target Boundaries", ints(&[1, 2, 3, 4])),
                ("Temperature", Value::real(0.0)),
            ]),
        ],
        ..Model::default()
    };
    let sim = run(model, &mesh).unwrap();
    let t = temperature(&sim, n);
    // The dump has 6 significant digits (values up to 10).
    let want = reference("thermal_pipe2d");
    assert!(
        max_diff(&t, &want) < 1e-5,
        "max diff {:e}",
        max_diff(&t, &want)
    );
}

#[test]
fn transient_step_applies_dt_once() {
    // Unit triangle: nodes 2 and 3 held at 0, node 1 starts at 1. With rho cp
    // = 6, k = 1, dt = 1 one backward-Euler step gives T1 = 1/3: `Initialize
    // Dirichlet Conditions` (default on) sets nodes 2 and 3 to 0 before the
    // step, so the previous values are (1, 0, 0) and
    // T1 = (rho cp / (12 dt)) / (K11 + 2 rho cp / (24 dt)) = 0.5 / 1.5.
    // The pristine ElmerSolver gives the same (max T = 3.333333333333E-001).
    let mesh = SimMesh::new(read_elmer_mesh(&data("minimal_example"), 2));
    let mut model = minimal_model();
    model.simulation = sec(&[
        ("Simulation Type", st("transient")),
        ("Timestep Intervals", ints(&[1])),
        ("Timestep Sizes", Value::real(1.0)),
        ("Output Intervals", ints(&[0])),
    ]);
    model.materials[0].insert("Density", Value::real(2.0));
    model.materials[0].insert("Heat Capacity", Value::real(3.0));
    model.boundary_conditions.pop();
    model.bodies[0].insert("Initial Condition", ints(&[1]));
    model.initial_conditions = vec![sec(&[("Temperature", Value::real(1.0))])];
    let sim = run(model, &mesh).unwrap();
    let t = temperature(&sim, 3);
    assert!((t[0] - 1.0 / 3.0).abs() < 1e-13, "{t:?}");
    assert_eq!(t[1], 0.0);
    assert_eq!(t[2], 0.0);
}

#[test]
fn unsupported_cases_are_refused_at_setup() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("minimal_example"), 2));
    let try_model = |m: Model| {
        let (mut sim, reg) = simulation(m);
        sim.setup(&mesh, &reg, &NullLogger, &SerialComm)
            .unwrap_err()
            .to_string()
    };

    let mut m = minimal_model();
    m.boundary_conditions[1].insert("Radiation", st("Idealized"));
    assert!(try_model(m).contains("Radiation"));

    // Reading the flow solution is refused when the solver runs, not at
    // setup: Elmer's own init accepts it, and `Exec Solver = Never` skips it.
    let mut m = minimal_model();
    m.equations[0].insert("Convection", st("Computed"));
    m.solvers[0].insert("Stabilize", Value::Logical(true));
    let err = run(m, &mesh).err().expect("refused").to_string();
    assert!(err.contains("Convection = Computed"), "{err}");

    let mut m = minimal_model();
    m.materials[0].insert("Compressibility Model", st("Thermal"));
    let err = run(m, &mesh).err().expect("refused").to_string();
    assert!(err.contains("Compressibility Model"), "{err}");

    // Convection without SUPG would need Elmer's residual-free bubble.
    let mut m = minimal_model();
    m.equations[0].insert("Convection", st("Constant"));
    assert!(try_model(m).contains("bubble"));

    let mut m = minimal_model();
    m.materials[0].insert("Heat Conductivity At IP", Value::Logical(true));
    assert!(try_model(m).contains("At IP"));
}
