//! StatCurrentSolve end to end against the pristine ElmerSolver (`data/README.md`
//! records the runs): `fem/tests/current` (conductivity table in `Time`,
//! current density), `fem/tests/current_heat_control` (`Current Control`), and
//! two derived cases for tensor conductivity with Robin conditions and the
//! field recovery, and for the displacement-current mass term.

mod support;

use elmer_core::{Dependent, DependentKind, ParamMap, Params, Value};
use elmer_solver::Model;
use support::cases::*;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn stat_current_solver(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("Static Current")),
        ("Variable", st("Potential")),
        ("Variable DOFs", ints(&[1])),
        (
            "Procedure",
            st("\"StatCurrentSolve\" \"StatCurrentSolver\""),
        ),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
        ("Nonlinear System Max Iterations", ints(&[1])),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

fn step_model(
    simulation: ParamMap,
    solver: ParamMap,
    material: ParamMap,
    bcs: Vec<ParamMap>,
) -> Model {
    Model {
        simulation,
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Body Force", ints(&[1])),
        ])],
        materials: vec![material],
        body_forces: vec![ParamMap::new()],
        boundary_conditions: bcs,
        ..Model::default()
    }
}

fn steady(iterations: i32) -> ParamMap {
    sec(&[
        ("Coordinate System", st("Cartesian 2D")),
        ("Simulation Type", st("Steady")),
        ("Steady State Max Iterations", ints(&[iterations])),
    ])
}

fn step_mesh(case: &str) -> SimMesh {
    SimMesh::new(read_elmer_mesh(&data(case), 2))
}

#[test]
fn current_time_table_and_current_density_match_elmer() {
    // Electric Conductivity = Variable Time: 1e-3 at time 1, 1e-4 at 2 and 3;
    // two steady-state iterations, so the last assembly sees time 2.
    let mesh = step_mesh("current");
    let table = Value::Dependent(Dependent {
        variables: vec!["Time".into()],
        kind: DependentKind::Table {
            x: vec![1.0, 2.0, 3.0],
            y: vec![1.0e-3, 1.0e-4, 1.0e-4],
            cubic: false,
        },
        rows: 1,
        cols: 1,
    });
    // elmer-solver judges the steady state reached after the first step
    // (RunReport::steps == 1) where the pristine run takes two, so the
    // minimum is stated here.
    let mut simulation = steady(2);
    simulation.insert("Steady State Min Iterations", ints(&[2]));
    let m = step_model(
        simulation,
        stat_current_solver(&[]),
        sec(&[("Electric Conductivity", table)]),
        vec![
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("Current Density", real(1.0e-3)),
            ]),
            sec(&[("Target Boundaries", ints(&[2])), ("Potential", real(0.0))]),
        ],
    );
    let sim = run(m, &mesh).unwrap();
    let want = &read_columns(&data("current").join("solution.dat"))[0];
    let d = rel_max_diff(&nodal(&sim, "potential", 0, mesh.mesh.coords.len()), want);
    assert!(d < 1e-10, "{d:e}");
    // The sif's Reference Norm is 7.68119284219.
    assert_rel(norm(&sim, "potential"), 7.6811928421903835, 1e-11);
}

#[test]
fn current_control_scales_solution_and_reports_totals() {
    let mesh = step_mesh("current_control");
    let m = step_model(
        steady(1),
        stat_current_solver(&[("Current Control", real(0.1))]),
        sec(&[("Electric Conductivity", real(1.0e-3))]),
        vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Potential", real(1.0))]),
            sec(&[("Target Boundaries", ints(&[2])), ("Potential", real(0.0))]),
        ],
    );
    let sim = run(m, &mesh).unwrap();
    let want = &read_columns(&data("current_control").join("solution.dat"))[0];
    let n = mesh.mesh.coords.len();
    let pot = nodal(&sim, "potential", 0, n);
    assert!(rel_max_diff(&pot, want) < 1e-10);
    let p = &sim.model.simulation;
    assert_rel(
        p.real("res: total joule heating").unwrap(),
        1.3140158060369138e-4,
        1e-10,
    );
    assert_rel(
        p.real("res: effective resistance").unwrap(),
        7610.258532703735,
        1e-10,
    );
    assert_rel(
        p.real("res: currentsolver scaling").unwrap(),
        761.0258532703736,
        1e-10,
    );
}

fn tensor_model() -> Model {
    // Symmetric full conductivity tensor, a Robin condition on boundary 1.
    let tensor = Value::Real {
        rows: 2,
        cols: 2,
        data: vec![1.0e-3, 2.0e-4, 2.0e-4, 5.0e-4],
    };
    step_model(
        steady(1),
        stat_current_solver(&[
            ("Calculate Joule Heating", Value::Logical(true)),
            ("Calculate Volume Current", Value::Logical(true)),
            ("Calculate Electric Field", Value::Logical(true)),
            ("Calculate Nodal Heating", Value::Logical(true)),
            ("Calculate Nodal Current", Value::Logical(true)),
        ]),
        sec(&[
            ("Electric Conductivity", tensor),
            ("Relative Permittivity", real(4.0)),
        ]),
        vec![
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("External Conductivity", real(3.0e-3)),
                ("External Potential", real(2.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[2])), ("Potential", real(0.0))]),
        ],
    )
}

#[test]
fn tensor_conductivity_robin_and_fields_match_elmer() {
    let mesh = step_mesh("current");
    let sim = run(tensor_model(), &mesh).unwrap();
    let n = mesh.mesh.coords.len();
    let cols = read_columns(&data("current_tensor").join("solution.dat"));
    // potential, joule heating, nodal joule heating, volume current 1-2,
    // nodal current 1-2, electric field 1-2.
    let checks: [(&str, usize, usize); 9] = [
        ("potential", 0, 0),
        ("joule heating", 0, 1),
        ("nodal joule heating", 0, 2),
        ("volume current", 0, 3),
        ("volume current", 1, 4),
        ("nodal current", 0, 5),
        ("nodal current", 1, 6),
        ("electric field", 0, 7),
        ("electric field", 1, 8),
    ];
    let bad: Vec<String> = checks
        .iter()
        .filter_map(|&(var, k, col)| {
            let d = rel_max_diff(&nodal(&sim, var, k, n), &cols[col]);
            (d >= 1e-9).then(|| format!("{var} {k}: {d:e}"))
        })
        .collect();
    assert!(bad.is_empty(), "{bad:?}");
    assert_rel(
        sim.model
            .simulation
            .real("res: total joule heating")
            .unwrap(),
        3.306_196_933_777_63e-4,
        1e-10,
    );
}

#[test]
fn transient_displacement_current_matches_elmer() {
    // BDF1, six steps of 5 ns with the mass term eps0 eps_r grad.grad, from a
    // zero potential with 1 V on boundary 1.
    let mesh = step_mesh("current");
    let mut m = step_model(
        sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Transient")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Timestepping Method", st("BDF")),
            ("BDF Order", ints(&[1])),
            ("Timestep Sizes", real(5.0e-9)),
            ("Timestep Intervals", ints(&[6])),
        ]),
        stat_current_solver(&[]),
        sec(&[
            ("Electric Conductivity", real(1.0e-3)),
            ("Relative Permittivity", real(4.0)),
        ]),
        vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Potential", real(1.0))]),
            sec(&[("Target Boundaries", ints(&[2])), ("Potential", real(0.0))]),
        ],
    );
    m.simulation.insert("Output Intervals", ints(&[6]));
    let sim = run(m, &mesh).unwrap();
    let want = &read_columns(&data("current_transient").join("solution.dat"))[0];
    let d = rel_max_diff(&nodal(&sim, "potential", 0, mesh.mesh.coords.len()), want);
    // Six BDF1 solves of an ill-scaled system (1e-11 mass over 5e-9 s against
    // 1e-3 conductivity) agree to the solvers' conditioning, not to round-off.
    assert!(d < 1e-7, "{d:e}");
}
