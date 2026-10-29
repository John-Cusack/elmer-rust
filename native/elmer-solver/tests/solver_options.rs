//! Keyword-driven behaviour without a Fortran capture: MATC keyword values
//! through a shared engine, `Linear System Residual Mode`, the stationary
//! iterative methods, and second-order time stepping. Each test checks an
//! exact property (identical results, a recurrence written out from the
//! Fortran formulas, a solved system), not Elmer parity.

mod support;

use std::cell::RefCell;
use std::path::Path;
use std::rc::Rc;
use std::sync::Arc;

use elmer_core::{
    CsrMatrix, Dependent, DependentKind, ElemId, ElementRef, ElementType, LinearOperator,
    LocalSystem, MeshView, NullLogger, ParamMap, SerialComm, Value,
};
use elmer_expr::MatcEngine;
use elmer_solver::backend::LinearSolverChoice;
use elmer_solver::linsys::elmer_registry;
use elmer_solver::{ElementAssembly, Model, Simulation, SolverCtx, SolverError, SolverRegistry};
use support::mesh::FileMesh;
use support::poisson::poisson;
use support::{i, r, s, sec, Capture};

fn line_mesh() -> FileMesh {
    FileMesh::read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/1d"), 1)
}

/// `-u'' = 1` on the 1dtests mesh with `u = 0` at both ends.
fn poisson_model(transient: bool, solver_extra: &[(&str, Value)]) -> Model {
    let mut simulation = sec(&[("Steady State Max Iterations", i(&[2]))]);
    if transient {
        simulation.insert("Simulation Type", s("transient"));
        simulation.insert("Timestep Intervals", i(&[4]));
        simulation.insert("Timestep Sizes", r(1.0));
        simulation.insert("Steady State Max Iterations", i(&[1]));
    } else {
        simulation.insert("Simulation Type", s("steady state"));
    }
    let mut solver = sec(&[
        ("Equation", s("poisson")),
        ("Variable", s("Potential")),
        ("Procedure", s("\"Poisson\" \"PoissonSolver\"")),
        ("Linear System Solver", s("direct")),
    ]);
    for (k, v) in solver_extra {
        solver.insert(k, v.clone());
    }
    Model {
        simulation,
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![sec(&[("Equation", i(&[1])), ("Body Force", i(&[1]))])],
        body_forces: vec![sec(&[("Source", r(1.0))])],
        boundary_conditions: vec![sec(&[("Target Boundaries", i(&[1, 2])), ("Potential", r(0.0))])],
        ..Model::default()
    }
}

fn run(model: Model, before_setup: impl FnOnce(&mut Simulation)) -> Simulation {
    let mesh = line_mesh();
    let mut reg = SolverRegistry::new();
    reg.register("Poisson", "PoissonSolver", poisson);
    let mut sim = Simulation::new(model);
    before_setup(&mut sim);
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    sim
}

fn potential(sim: &Simulation) -> Vec<f64> {
    sim.vars.get("potential").unwrap().values.clone()
}

/// A MATC stand-in with state, as the real binding has: it is `!Send`
/// (`Rc`/`RefCell`) and records every call, so the test can see the
/// solver used the same engine it holds.
#[derive(Default)]
struct StubMatc {
    calls: RefCell<Vec<(String, String, Vec<f64>)>>,
}

impl MatcEngine for StubMatc {
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.eval_with("tx", expr, args)
    }

    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.calls.borrow_mut().push((var.into(), expr.into(), args.to_vec()));
        match expr {
            "0.5" => Ok(vec![0.5]),
            "st - 2.5" => Ok(vec![args[0] - 2.5]),
            _ => Err(format!("stub cannot evaluate {expr}")),
        }
    }
}

fn matc(expr: &str) -> Value {
    Value::Dependent(Dependent {
        variables: Vec::new(),
        kind: DependentKind::Matc(expr.into()),
        rows: 1,
        cols: 1,
    })
}

#[test]
fn matc_relaxation_factor_reaches_compute_change() {
    let plain = run(poisson_model(false, &[("Nonlinear System Relaxation Factor", r(0.5))]), |_| {});
    let engine = Rc::new(StubMatc::default());
    let shared = Rc::clone(&engine);
    let via_matc = run(
        poisson_model(false, &[("Nonlinear System Relaxation Factor", matc("0.5"))]),
        move |sim| sim.set_matc(shared),
    );
    assert_eq!(potential(&plain), potential(&via_matc));
    // Relaxed: the first steady iteration keeps half of the solution.
    let full = run(poisson_model(false, &[]), |_| {});
    assert_ne!(potential(&plain), potential(&full));
    // ListGetCReal binds `time` as `st` for a MATC value with no Variable.
    let calls = engine.calls.borrow();
    assert!(!calls.is_empty());
    assert!(calls.iter().all(|(var, e, _)| var == "st" && e == "0.5"), "{calls:?}");
}

#[test]
fn matc_exec_condition_follows_time() {
    let engine = Rc::new(StubMatc::default());
    let shared = Rc::clone(&engine);
    let sim = run(poisson_model(true, &[("Exec Condition", matc("st - 2.5"))]), move |sim| {
        sim.set_matc(shared)
    });
    // Steps at t = 1..4; the solver runs where t - 2.5 >= 0.
    assert_eq!(sim.slots[0].times_visited, 2);
    let times: Vec<f64> = engine.calls.borrow().iter().map(|c| c.2[0]).collect();
    assert_eq!(times, vec![1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn matc_without_engine_is_an_error() {
    let mesh = line_mesh();
    let mut reg = SolverRegistry::new();
    reg.register("Poisson", "PoissonSolver", poisson);
    let mut sim = Simulation::new(poisson_model(false, &[("Nonlinear System Relaxation Factor", matc("0.5"))]));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    let err = sim.run(&mesh, &NullLogger, &SerialComm).unwrap_err();
    assert!(err.to_string().contains("no MATC engine installed"), "{err}");
}

#[test]
fn residual_mode_gives_the_plain_solution() {
    // The second steady iteration starts from a nonzero iterate, so the
    // residual system b - A x0 is not b.
    let plain = run(poisson_model(false, &[]), |_| {});
    let residual = run(poisson_model(false, &[("Linear System Residual Mode", Value::Logical(true))]), |_| {});
    for (a, b) in potential(&plain).iter().zip(potential(&residual)) {
        assert!((a - b).abs() <= 1e-13 * a.abs().max(1.0), "{a} vs {b}");
    }
    let relaxed = run(
        poisson_model(
            false,
            &[
                ("Linear System Residual Mode", Value::Logical(true)),
                ("Nonlinear System Relaxation Factor", r(0.5)),
            ],
        ),
        |_| {},
    );
    let want = run(poisson_model(false, &[("Nonlinear System Relaxation Factor", r(0.5))]), |_| {});
    for (a, b) in potential(&want).iter().zip(potential(&relaxed)) {
        assert!((a - b).abs() <= 1e-13 * a.abs().max(1.0), "{a} vs {b}");
    }
}

/// `n x n` tridiagonal `[-0.5, 4, -0.5]`: diagonally dominant, and its
/// row sums keep Richardson's lumped iteration contractive.
fn tridiag(n: usize) -> CsrMatrix<f64> {
    let (mut row_ptr, mut col_idx, mut values, mut diag) = (vec![0], Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        for j in i.saturating_sub(1)..(i + 2).min(n) {
            if i == j {
                diag.push(values.len());
            }
            col_idx.push(j as u32);
            values.push(if i == j { 4.0 } else { -0.5 });
        }
        row_ptr.push(values.len());
    }
    CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr,
        col_idx,
        values,
        diag: Some(diag),
    }
}

#[test]
fn stationary_methods_are_built_from_the_matrix() {
    let reg = elmer_registry();
    let a = tridiag(12);
    let b = vec![1.0; 12];
    for m in ["jacobi", "sgs", "richardson"] {
        let p = sec(&[
            ("Linear System Solver", s("iterative")),
            ("Linear System Iterative Method", s(m)),
            ("Linear System Convergence Tolerance", r(1e-12)),
            ("Linear System Max Iterations", i(&[500])),
        ]);
        assert!(matches!(reg.select(&p, &NullLogger), Err(SolverError::NotWired(_))), "{m}");
        let LinearSolverChoice::Iterative { mut solver, precond } = reg.select_for(&p, &a, &NullLogger).unwrap() else {
            panic!("{m}: not iterative");
        };
        assert!(precond.is_none());
        let mut x = vec![0.0; 12];
        let stats = solver.solve(&a, &elmer_core::IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
        let mut ax = vec![0.0; 12];
        a.apply(&x, &mut ax);
        let res = ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt();
        assert!(res < 1e-11, "{m}: residual {res} after {} iterations", stats.iterations);
    }
}

/// One point carrying a mass-spring oscillator `m u'' + k u = 0`.
struct Point;

impl MeshView for Point {
    fn dim(&self) -> usize {
        1
    }
    fn n_nodes(&self) -> usize {
        1
    }
    fn coord(&self, _: u32) -> [f64; 3] {
        [0.0; 3]
    }
    fn n_bulk(&self) -> usize {
        1
    }
    fn n_boundary(&self) -> usize {
        0
    }
    fn element(&self, _: ElemId) -> ElementRef<'_> {
        ElementRef {
            etype: ElementType(101),
            nodes: &[0],
            tag: 1,
            parents: [None, None],
        }
    }
}

const M: f64 = 2.0;
const K: f64 = 50.0;

struct Spring;

impl ElementAssembly for Spring {
    fn assemble(&mut self, _: &SolverCtx<'_>, _: ElemId, _: f64, _: bool, local: &mut LocalSystem) -> Result<(), SolverError> {
        local.stiff[0] = K;
        local.mass.as_mut().unwrap()[0] = M;
        Ok(())
    }
}

/// The scalar form of `InitializeTimestep` (second order),
/// `Time2ndOrder` and `Update2ndOrder` (TimeIntegrate.F90:431, 568) with
/// the zero initial history Elmer starts from.
fn second_order_reference(param_a: f64, param_b: f64, u0: f64, dt: f64, steps: usize) -> Vec<f64> {
    let (al, be, ga, de) = if param_b >= 0.0 {
        let rb = param_b;
        ((2.0 * rb - 1.0) / (rb + 1.0), 1.0 / ((rb + 1.0) * (rb + 1.0)), (3.0 - rb) / (2.0 * (rb + 1.0)), rb / (rb + 1.0))
    } else {
        (param_a, (1.0 - param_a).powi(2) / 4.0, 0.5 - param_a, 0.0)
    };
    let mut c = [0.0_f64; 8]; // 1-based PrevValues columns
    let mut u = u0;
    let mut out = Vec::new();
    for _ in 0..steps {
        c[3] = u;
        c[4] = c[1];
        c[5] = c[2];
        c[7] = c[6];
        let a2b = if param_b >= 0.0 { c[7] } else { c[5] };
        let f = (1.0 - al) / (1.0 - de) / (be * dt * dt) * M * c[3] + (1.0 - al) / (1.0 - de) / (be * dt) * M * c[4]
            + de / (1.0 - de) * M * c[5]
            - ((1.0 - al) * (1.0 - 1.0 / (2.0 * be)) + al) / (1.0 - de) * M * a2b;
        let keff = K + (1.0 - al) / (1.0 - de) / (be * dt * dt) * M;
        u = f / keff;
        let (dx, v, a, a2) = (u - c[3], c[4], c[5], c[7]);
        if de == 0.0 {
            c[1] = ga / (be * dt) * dx - (ga - be) / be * v - dt * (ga - 2.0 * be) / (2.0 * be) * a;
            c[2] = 1.0 / (be * dt * dt) * dx - 1.0 / (be * dt) * v - (1.0 - 2.0 * be) / (2.0 * be) * a;
        } else {
            c[1] = ga / (be * dt) * dx - (ga - be) / be * v - dt * (ga - 2.0 * be) / (2.0 * be) * a2;
            c[2] = (1.0 - al) / (1.0 - de) / (be * dt * dt) * dx - (1.0 - al) / (1.0 - de) / (be * dt) * v
                - de / (1.0 - de) * a
                - (1.0 - al - 2.0 * be) / (1.0 - de) / (2.0 * be) * a2;
            c[6] = (1.0 - de) / (1.0 - al) * c[2] + de / (1.0 - al) * a - al / (1.0 - al) * a2;
        }
        out.push(u);
    }
    out
}

/// Runs the oscillator through `Assembly Solver` with `Time Derivative
/// Order = 2`, recording `u` after every step.
fn oscillator(method: &[(&str, Value)], dt: f64, steps: i32) -> Vec<f64> {
    let mut solver = sec(&[
        ("Assembly Solver", Value::Logical(true)),
        ("Equation", s("spring")),
        ("Variable", s("U")),
        ("Time Derivative Order", i(&[2])),
        ("Linear System Solver", s("direct")),
        ("Bulk Assembly Procedure", s("\"Spring\" \"Spring\"")),
    ]);
    for (k, v) in method {
        solver.insert(k, v.clone());
    }
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", s("transient")),
            ("Timestep Intervals", i(&[steps])),
            ("Timestep Sizes", r(dt)),
        ]),
        solvers: vec![solver, sec(&[("Exec Solver", s("after timestep")), ("Procedure", s("\"Rec\" \"Rec\""))])],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![sec(&[("Equation", i(&[1])), ("Initial Condition", i(&[1]))])],
        initial_conditions: vec![sec(&[("U", r(1.0))])],
        ..Model::default()
    };
    thread_local!(static SEEN: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) });
    struct Rec;
    impl elmer_solver::Solver for Rec {
        fn solve(&mut self, ctx: &mut SolverCtx<'_>, _: f64, _: bool) -> Result<(), SolverError> {
            let u = ctx.variable("u").unwrap().values[0];
            SEEN.with(|s| s.borrow_mut().push(u));
            Ok(())
        }
    }
    fn rec(_: &ParamMap) -> Result<Box<dyn elmer_solver::Solver>, SolverError> {
        Ok(Box::new(Rec))
    }
    SEEN.with(|s| s.borrow_mut().clear());
    let mut reg = SolverRegistry::new();
    reg.register_assembly("Spring", "Spring", Arc::new(|_| Ok(Box::new(Spring))));
    reg.register("Rec", "Rec", rec);
    let mut sim = Simulation::new(model);
    sim.setup(&Point, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&Point, &NullLogger, &SerialComm).unwrap();
    assert_eq!(sim.vars.get("u").unwrap().prev_values.len(), 7);
    SEEN.with(|s| s.borrow().clone())
}

#[test]
fn bossak_and_generalized_alpha_follow_the_fortran_recurrence() {
    let (dt, steps) = (0.05, 40);
    for (method, a, b) in [
        (vec![], -0.05, -1.0),
        (vec![("Bossak Alpha", r(-0.1))], -0.1, -1.0),
        (vec![("Timestepping Method", s("newmark")), ("Bossak Alpha", r(0.0))], 0.0, -1.0),
        (vec![("Timestepping Method", s("generalized-alpha"))], 0.0, 0.75),
        (
            vec![("Timestepping Method", s("generalized-alpha")), ("Generalized-alpha rinf", r(0.3))],
            0.0,
            0.3,
        ),
    ] {
        let got = oscillator(&method, dt, steps);
        let want = second_order_reference(a, b, 1.0, dt, steps as usize);
        assert_eq!(got.len(), want.len());
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!((g - w).abs() <= 1e-12, "{method:?} step {}: {g} vs {w}", k + 1);
        }
        // The oscillator actually oscillates: period 2 pi sqrt(m/k) = 1.26.
        assert!(got.iter().any(|&u| u < -0.5), "{method:?}: {got:?}");
    }
}
