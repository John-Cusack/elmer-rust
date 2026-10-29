//! A `Real Procedure` keyword value served by the simulation's function
//! registry, with the user function reading the run through the live
//! model view. Checked against the pristine ElmerSolver (26.2-devel, Rev
//! 9f6af2f85) run of `fem/tests/1dtests`: mesh from `ElmerGrid 1 2 1d
//! -decimals 20` (`tests/data/1d/`), messages in
//! `tests/data/1dtests_stdout.txt`.
//!
//! The body force is `elmer-fixtures`' port of `1dtests/Source.f90`: from
//! step 2 on it checks the previous step's `Potential` at every node
//! against its analytic form and panics with the Fortran's `Fatal` text if
//! the view showed it stale values. The Poisson solver is a test-local port
//! of `1dtests/Poisson.f90` (`support::poisson`).

mod support;

use std::path::Path;
use std::sync::Arc;

use elmer_core::{Dependent, DependentKind, SerialComm, Value};
use elmer_expr::FieldLookup;
use elmer_fixtures::ModelView;
use elmer_solver::{LiveModel, Model, Simulation, SolverRegistry};
use support::poisson::poisson;
use support::{i, r, s, sec, split_change, Capture};

const STDOUT: &str = include_str!("data/1dtests_stdout.txt");

/// The fixture crate's model trait over the solver's live view.
struct View(Arc<LiveModel>);

impl FieldLookup for View {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        self.0.value(var, node)
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.0.coord(node)
    }
    fn time(&self) -> f64 {
        self.0.time()
    }
}

impl ModelView for View {
    fn n_nodes(&self) -> usize {
        self.0.n_nodes()
    }
}

fn model() -> Model {
    let source = Value::Dependent(Dependent {
        variables: vec!["Potential".into()],
        kind: DependentKind::Procedure {
            library: "Source".into(),
            entry: "Source".into(),
        },
        rows: 1,
        cols: 1,
    });
    Model {
        simulation: sec(&[
            ("Coordinate System", s("cartesian 1d")),
            ("Simulation Type", s("transient")),
            ("Timestepping Method", s("bdf")),
            ("BDF Order", i(&[1])),
            ("Timestep Intervals", i(&[5])),
            ("Timestep Sizes", r(1.0)),
            ("Output Intervals", i(&[1])),
            ("Steady State Max Iterations", i(&[50])),
        ]),
        solvers: vec![sec(&[
            ("Equation", s("poisson")),
            ("Variable", s("Potential")),
            ("Variable DOFs", i(&[1])),
            ("Procedure", s("\"Poisson\" \"PoissonSolver\"")),
            ("Linear System Solver", s("direct")),
            ("Steady State Convergence Tolerance", r(1e-9)),
            ("Reference Norm", r(1.40257362810)),
            ("Reference Norm Tolerance", r(1e-4)),
        ])],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![sec(&[("Equation", i(&[1])), ("Body Force", i(&[1]))])],
        body_forces: vec![sec(&[("Source", source)])],
        boundary_conditions: vec![sec(&[("Target Boundaries", i(&[1, 2])), ("Potential", r(0.0))])],
        ..Model::default()
    }
}

#[test]
fn one_d_tests_source_reads_the_live_model() {
    let mesh = support::mesh::FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/1d"), 1);
    let mut reg = SolverRegistry::new();
    reg.register("Poisson", "PoissonSolver", poisson);
    let mut sim = Simulation::new(model());
    let view = Arc::new(View(sim.live_model(&mesh)));
    assert_eq!(elmer_fixtures::registry::install("1dtests", view, &sim.backends.udfs), 1);
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let report = sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!(report.steps, 5);
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));

    let got = log.lines(&["ComputeChange"]);
    let want: Vec<&str> = STDOUT.lines().filter(|l| l.starts_with("ComputeChange")).collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-6 * wr.max(1e-3), "{g}\n{w}");
    }
}

#[test]
fn live_model_follows_published_values() {
    let mesh = support::mesh::FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/1d"), 1);
    let mut reg = SolverRegistry::new();
    reg.register("Poisson", "PoissonSolver", poisson);
    let mut sim = Simulation::new(model());
    let live = sim.live_model(&mesh);
    assert_eq!(live.n_nodes(), 61);
    // Nothing set up yet: no Potential, the time global at zero.
    assert_eq!(live.value("potential", 0), None);
    assert_eq!(live.time(), 0.0);
    elmer_fixtures::registry::install("1dtests", Arc::new(View(Arc::clone(&live))), &sim.backends.udfs);
    sim.setup(&mesh, &reg, &Capture::default(), &SerialComm).unwrap();
    sim.run(&mesh, &Capture::default(), &SerialComm).unwrap();
    assert_eq!(live.time(), 5.0);
    // Boundary nodes are fixed to zero; the published values are the solution.
    let v = sim.vars.get("potential").unwrap();
    for node in 0..61 {
        let slot = v.slot(node).unwrap();
        assert_eq!(live.value("potential", node), Some(v.values[slot]));
    }
}
