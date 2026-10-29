//! Running `FlowSolver` through `elmer-solver`'s `Simulation` with its
//! default backends (assembly, faer LU for `Direct`), the test-local basis
//! standing in for `elmer-fe` and hand-built model sections standing in for
//! `elmer-sif`.

use std::path::{Path, PathBuf};

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_flow::{install_fe, register};
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry};

use super::mesh::TestMesh;
use super::{TestBasis, TestQuadrature};

pub fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case)
}

pub fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

pub fn ints(v: &[i32]) -> Value {
    Value::Integer(v.to_vec())
}

pub fn real(x: f64) -> Value {
    Value::real(x)
}

pub fn st(s: &str) -> Value {
    Value::Str(s.into())
}

pub fn logical(b: bool) -> Value {
    Value::Logical(b)
}

pub fn install() {
    // Idempotent across tests in this binary.
    let _ = install_fe(|| (Box::new(TestBasis), Box::new(TestQuadrature)));
}

pub fn run(model: Model, mesh: &TestMesh) -> Result<Simulation, SolverError> {
    install();
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// A `FlowSolve` solver section with `extra` keywords and a direct solve.
pub fn flow_solver(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("Navier-Stokes")),
        ("Procedure", st("\"FlowSolve\" \"FlowSolver\"")),
        ("Variable", st("Flow Solution[Velocity:2 Pressure:1]")),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

/// A steady single-body model around `solver`.
pub fn steady_model(solver: ParamMap, material: ParamMap, body_force: Option<ParamMap>, bcs: Vec<ParamMap>) -> Model {
    let mut body = vec![
        ("Target Bodies", ints(&[1])),
        ("Equation", ints(&[1])),
        ("Material", ints(&[1])),
    ];
    if body_force.is_some() {
        body.push(("Body Force", ints(&[1])));
    }
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        constants: sec(&[(
            "Gravity",
            Value::Real {
                rows: 4,
                cols: 1,
                data: vec![0.0, -1.0, 0.0, 9.82],
            },
        )]),
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&body)],
        materials: vec![material],
        body_forces: body_force.into_iter().collect(),
        boundary_conditions: bcs,
        ..Model::default()
    }
}

/// Node-ordered `(u, v, p)` from a pristine ElmerSolver `.result`, as
/// `node u v p` lines.
pub fn read_solution(path: &Path) -> Vec<[f64; 3]> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .enumerate()
        .map(|(i, l)| {
            let t: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(t[0].parse::<usize>().unwrap(), i + 1);
            [t[1].parse().unwrap(), t[2].parse().unwrap(), t[3].parse().unwrap()]
        })
        .collect()
}

/// The flow variable by node, as `(u, v, p)`.
pub fn nodal(sim: &Simulation, n_nodes: usize) -> Vec<[f64; 3]> {
    let v = sim.vars.get("flow solution").unwrap();
    (0..n_nodes)
        .map(|n| {
            let s = v.slot(n as u32).unwrap() * 3;
            [v.values[s], v.values[s + 1], v.values[s + 2]]
        })
        .collect()
}

/// Elmer's norm of the flow variable after the last `ComputeChange`.
pub fn norm(sim: &Simulation) -> f64 {
    sim.vars.get("flow solution").unwrap().conv.norm
}
