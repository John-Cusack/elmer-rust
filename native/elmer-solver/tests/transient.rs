//! Transient stepping with the default time backend: history commit
//! (InitializeTimestep) plus `elmer-timeint` local schemes, on du/dt = -k u
//! with one dof. Discrete recurrences, including start-up, are the
//! reference (local correctness, not an Elmer capture).

use std::cell::RefCell;

use elmer_core::{ElemId, ElementRef, ElementType, LocalSystem, MeshView, NullLogger, ParamMap, SerialComm, Value};
use elmer_solver::{Model, Simulation, Solver, SolverCtx, SolverError, SolverRegistry};

const K: f64 = 2.0;
const DT: f64 = 0.1;
const STEPS: i32 = 6;

struct Point;

impl MeshView for Point {
    fn dim(&self) -> usize {
        1
    }
    fn n_nodes(&self) -> usize {
        1
    }
    fn coord(&self, _node: u32) -> [f64; 3] {
        [0.0; 3]
    }
    fn n_bulk(&self) -> usize {
        1
    }
    fn n_boundary(&self) -> usize {
        0
    }
    fn element(&self, _e: ElemId) -> ElementRef<'_> {
        ElementRef {
            etype: ElementType(101),
            nodes: &[0],
            tag: 1,
            parents: [None, None],
        }
    }
}

thread_local! {
    static TRACE: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) };
}

/// M = 1, K = k, F = 0, through `default_first_order_time`; the 1x1 system
/// is solved in place.
struct Decay;

impl Solver for Decay {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        assert!(transient);
        let mut l = LocalSystem::zeros(1);
        l.stiff[0] = K;
        l.mass = Some(vec![1.0]);
        ctx.default_first_order_time(&mut l, ElemId(0), dt)?;
        let u = l.force[0] / l.stiff[0];
        ctx.primary_mut().unwrap().values[0] = u;
        TRACE.with(|t| t.borrow_mut().push(u));
        Ok(())
    }
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn run(method: &[(&str, Value)]) -> Vec<f64> {
    TRACE.with(|t| t.borrow_mut().clear());
    let mut solver = vec![
        ("Equation", Value::Str("decay".into())),
        ("Variable", Value::Str("u".into())),
        ("Procedure", Value::Str("\"Test\" \"Decay\"".into())),
    ];
    solver.extend(method.iter().cloned());
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("transient".into())),
            ("Timestep Intervals", Value::Integer(vec![STEPS])),
            ("Timestep Sizes", Value::real(DT)),
        ]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1])), ("Initial Condition", Value::Integer(vec![1]))])],
        initial_conditions: vec![sec(&[("u", Value::real(1.0))])],
        ..Model::default()
    };
    let mut reg = SolverRegistry::new();
    reg.register("Test", "Decay", |_| Ok(Box::new(Decay)));
    let mut sim = Simulation::new(model);
    sim.setup(&Point, &reg, &NullLogger, &SerialComm).unwrap();
    let report = sim.run(&Point, &NullLogger, &SerialComm).unwrap();
    assert_eq!(report.steps, STEPS as usize);
    let out = TRACE.with(|t| t.borrow().clone());
    assert_eq!(out.len(), STEPS as usize);
    out
}

fn close(got: &[f64], want: &[f64]) {
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!((g - w).abs() <= 1e-14 * w.abs(), "step {}: {g} vs {w}", i + 1);
    }
}

#[test]
fn backward_euler() {
    let got = run(&[("Timestepping Method", Value::Str("implicit euler".into()))]);
    let mut u = 1.0;
    let want: Vec<f64> = (0..STEPS)
        .map(|_| {
            u /= 1.0 + K * DT;
            u
        })
        .collect();
    close(&got, &want);
}

#[test]
fn crank_nicolson_starts_with_backward_euler() {
    let got = run(&[("Timestepping Method", Value::Str("crank-nicolson".into()))]);
    let mut u = 1.0 / (1.0 + K * DT);
    let mut want = vec![u];
    for _ in 1..STEPS {
        u *= (1.0 - K * DT / 2.0) / (1.0 + K * DT / 2.0);
        want.push(u);
    }
    close(&got, &want);
}

#[test]
fn bdf2_starts_with_bdf1() {
    let got = run(&[
        ("Timestepping Method", Value::Str("bdf".into())),
        ("BDF Order", Value::Integer(vec![2])),
    ]);
    let (mut um1, mut u) = (1.0, 1.0 / (1.0 + K * DT));
    let mut want = vec![u];
    for _ in 1..STEPS {
        let next = (4.0 * u - um1) / (3.0 + 2.0 * K * DT);
        um1 = u;
        u = next;
        want.push(u);
    }
    close(&got, &want);
}
