//! Stage S2 heat features end to end against the pristine ElmerSolver
//! (`tests/data/README.md` records the runs): BDF2 transient through
//! `elmer-timeint` (`fem/tests/heateq_bdf2`), a temperature-dependent cubic
//! conductivity table iterated to nonlinear convergence
//! (`fem/tests/heateq_spline_table`), and Robin and axisymmetric variants of
//! `fem/tests/ThermalPipe2D`.

mod support;

use elmer_core::{Dependent, DependentKind, Value};
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{
    bicgstab_ilu0, data, direct, heat_solver, ints, max_diff, nodal, norm, read_solution, reals,
    run, sec, st,
};
use support::meshio::read_elmer_mesh;

fn assert_rel(got: f64, want: f64, tol: f64) {
    let rel = (got - want).abs() / want.abs();
    assert!(rel < tol, "{got:e} vs {want:e}: relative error {rel:e}");
}

/// The single-body `angle` mesh of heateq_bdf2 and heateq_spline_table:
/// heat source 1 per unit mass, `Temp = 0` on boundary 1.
fn angle_model(simulation: elmer_core::ParamMap, solver: elmer_core::ParamMap) -> Model {
    Model {
        simulation,
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Body Force", ints(&[1])),
        ])],
        body_forces: vec![sec(&[("Heat Source", Value::real(1.0))])],
        boundary_conditions: vec![sec(&[
            ("Target Boundaries", ints(&[1])),
            ("Temp", Value::real(0.0)),
        ])],
        ..Model::default()
    }
}

#[test]
fn heateq_bdf2_transient_matches_elmer() {
    // BDF order 2 over 5 steps of 0.1 then 3 of 0.2: the first step runs
    // BDF1, later ones BDF2 with variable step sizes.
    let mesh = SimMesh::new(read_elmer_mesh(&data("heateq_bdf2"), 2));
    let simulation = sec(&[
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[2])),
        ("Timestep Sizes", reals(&[0.1, 0.2])),
        ("Timestep Intervals", ints(&[5, 3])),
        ("Output Intervals", ints(&[1, 1])),
        ("Steady State Max Iterations", ints(&[1])),
    ]);
    let mut solver = direct();
    solver.extend([
        ("Stabilize", Value::Logical(true)),
        ("Nonlinear System Convergence Tolerance", Value::real(1e-5)),
        ("Nonlinear System Max Iterations", ints(&[10])),
        ("Nonlinear System Relaxation Factor", Value::real(1.0)),
    ]);
    let mut model = angle_model(simulation, heat_solver("Temp", &solver));
    model.materials = vec![sec(&[
        ("Density", Value::real(1.0)),
        ("Heat Conductivity", Value::real(1.0)),
        ("Heat Capacity", Value::real(10.0)),
    ])];
    let sim = run(model, &mesh).unwrap();

    let n = mesh.mesh.coords.len();
    let t = nodal(&sim, "temp", n);
    let want = read_solution(&data("heateq_bdf2").join("solution.dat"));
    let d = max_diff(&t, &want);
    assert!(d < 1e-13, "max diff {d:e}");
    // The sif's Reference Norm is 5.214225320136149E-002.
    assert_rel(norm(&sim, "temp"), 0.05214225320136147, 1e-13);
}

/// `cond.dat`: `Heat Conductivity = Variable Temp; Real Cubic`.
fn conductivity_table() -> Value {
    let text = std::fs::read_to_string(data("heateq_spline_table").join("cond.dat")).unwrap();
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for l in text.lines() {
        let t: Vec<f64> = l.split_whitespace().map(|s| s.parse().unwrap()).collect();
        x.push(t[0]);
        y.push(t[1]);
    }
    Value::Dependent(Dependent {
        variables: vec!["Temp".into()],
        kind: DependentKind::Table { x, y, cubic: true },
        rows: 1,
        cols: 1,
    })
}

#[test]
fn heateq_spline_table_nonlinear_matches_elmer() {
    // k(T) = 1/(1+T) as a cubic spline: the Picard loop re-evaluates the
    // table at the current temperature each iteration. Elmer converges in
    // 10 iterations (relative change 1.4e-6 < 1e-5).
    let mesh = SimMesh::new(read_elmer_mesh(&data("heateq_spline_table"), 2));
    let simulation = sec(&[
        ("Simulation Type", st("Steady State")),
        ("Steady State Max Iterations", ints(&[1])),
    ]);
    let mut solver = bicgstab_ilu0();
    solver.extend([
        ("Stabilize", Value::Logical(true)),
        ("Nonlinear System Convergence Tolerance", Value::real(1e-5)),
        ("Nonlinear System Max Iterations", ints(&[100])),
        ("Nonlinear System Relaxation Factor", Value::real(1.0)),
    ]);
    let mut model = angle_model(simulation, heat_solver("Temp", &solver));
    model.body_forces = vec![sec(&[("Heat Source", Value::real(10.0))])];
    model.materials = vec![sec(&[
        ("Density", Value::real(1.0)),
        ("Heat Conductivity", conductivity_table()),
        ("Heat Capacity", Value::real(1.0)),
    ])];
    let sim = run(model, &mesh).unwrap();

    let v = sim.vars.get("temp").unwrap();
    assert_eq!(v.conv.nonlin_iter, 10);
    assert_eq!(v.conv.nonlin_converged, 1);
    // Full-precision norm of the Fortran run's final temperature; the sif's
    // Reference Norm is 1.5131000.
    assert_rel(norm(&sim, "temp"), 1.513099967071337, 1e-10);
}

/// ThermalPipe2D (three bodies, one material, `Temperature = 10` on the
/// exterior) with the pipe wall turned into a Robin condition.
fn pipe_model(coordinates: &str) -> Model {
    Model {
        simulation: sec(&[
            ("Simulation Type", st("Steady state")),
            ("Coordinate System", st(coordinates)),
        ]),
        solvers: vec![heat_solver("Temperature", &direct())],
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
                ("Heat Transfer Coefficient", Value::real(20.0)),
                ("External Temperature", Value::real(2.0)),
            ]),
        ],
        ..Model::default()
    }
}

#[test]
fn thermal_pipe_robin_matches_elmer() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("thermal_pipe2d"), 2));
    let sim = run(pipe_model("Cartesian"), &mesh).unwrap();
    let t = nodal(&sim, "temperature", mesh.mesh.coords.len());
    let want = read_solution(&data("thermal_pipe2d_robin").join("solution.dat"));
    let d = max_diff(&t, &want);
    assert!(d < 1e-11, "max diff {d:e}");
    assert_rel(norm(&sim, "temperature"), 6.6452169458174035, 1e-13);
}

#[test]
fn thermal_pipe_axisymmetric_matches_elmer() {
    // The ThermalPipe2D mesh moved to x in [1, 11] and revolved about the
    // y axis; the Robin pipe wall and a source in the ground (body 1) check
    // the r-weighting of both bulk and boundary integrals.
    let mut m = read_elmer_mesh(&data("thermal_pipe2d"), 2);
    m.coords.iter_mut().for_each(|c| c[0] += 6.0);
    let mesh = SimMesh::new(m);
    let mut model = pipe_model("Axi Symmetric");
    model.bodies[0].insert("Body Force", ints(&[1]));
    model.body_forces = vec![sec(&[("Heat Source", Value::real(0.01))])];
    let sim = run(model, &mesh).unwrap();
    let t = nodal(&sim, "temperature", mesh.mesh.coords.len());
    let want = read_solution(&data("thermal_pipe2d_axi").join("solution.dat"));
    let d = max_diff(&t, &want);
    assert!(d < 1e-10, "max diff {d:e}");
    assert_rel(norm(&sim, "temperature"), 35.21274966573315, 1e-13);
}
