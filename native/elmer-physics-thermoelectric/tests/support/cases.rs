//! Running `ThermoElectricSolver` through `elmer-solver`'s `Simulation` with
//! its default backends (assembly, faer LU for `Direct`), the `elmer-fe`
//! basis and Gauss rules, and hand-built model sections standing in for
//! `elmer-sif`.

use std::path::{Path, PathBuf};

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_thermoelectric::register;
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry};

use super::mesh::TestMesh;

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

/// A row-major `rows x cols` real array.
pub fn array(rows: usize, cols: usize, data: &[f64]) -> Value {
    Value::Real {
        rows,
        cols,
        data: data.to_vec(),
    }
}

pub fn st(s: &str) -> Value {
    Value::Str(s.into())
}

pub fn logical(b: bool) -> Value {
    Value::Logical(b)
}

pub fn run(model: Model, mesh: &TestMesh) -> Result<Simulation, SolverError> {
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// A `ThermoElectricSolver` section with `extra` keywords and a direct solve.
pub fn te_solver(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("ThermoElectric")),
        ("Variable", st("POT[Temperature:1 Potential:1]")),
        ("Procedure", st("\"ThermoElectricSolver\" \"ThermoElectricSolver\"")),
        ("Nonlinear System Convergence Tolerance", real(1.0e-6)),
        ("Nonlinear System Max Iterations", ints(&[100])),
        ("Nonlinear System Newton After Iterations", ints(&[1])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-9)),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

/// The `Material 1` of `ThermoElectric*/ptest.sif`.
pub fn material() -> ParamMap {
    sec(&[
        ("Electric Conductivity", real(1.0)),
        ("Density", real(1.0)),
        ("Heat Capacity", real(1.0)),
        ("Heat Conductivity", real(1.0)),
        ("Seebeck Coefficient", real(-0.005)),
    ])
}

pub fn simulation(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![("Coordinate System", st("Cartesian")), ("Steady State Max Iterations", ints(&[1]))];
    v.extend(extra.iter().cloned());
    sec(&v)
}

pub fn steady() -> ParamMap {
    simulation(&[("Simulation Type", st("Steady State"))])
}

/// A single-body model around `solver`.
pub fn model(simulation: ParamMap, solver: ParamMap, material: ParamMap, body_force: ParamMap, bcs: Vec<ParamMap>) -> Model {
    let body = vec![
        ("Target Bodies", ints(&[1])),
        ("Equation", ints(&[1])),
        ("Material", ints(&[1])),
        ("Body Force", ints(&[1])),
    ];
    Model {
        simulation,
        constants: ParamMap::new(),
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&body)],
        materials: vec![material],
        body_forces: vec![body_force],
        boundary_conditions: bcs,
        ..Model::default()
    }
}

/// Node-ordered `temperature potential` lines from a pristine ElmerSolver
/// `.result`, as `node t p` lines.
pub fn read_solution(path: &Path) -> Vec<Vec<f64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .enumerate()
        .map(|(i, l)| {
            let t: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(t[0].parse::<usize>().unwrap(), i + 1);
            t[1..].iter().map(|x| x.parse().unwrap()).collect()
        })
        .collect()
}

/// The coupled variable by node: `[temperature, potential]`.
pub fn nodal(sim: &Simulation, n_nodes: usize) -> Vec<Vec<f64>> {
    let v = sim.vars.get("pot").unwrap();
    (0..n_nodes)
        .map(|n| {
            let s = v.slot(n as u32).unwrap() * v.dofs;
            v.values[s..s + v.dofs].to_vec()
        })
        .collect()
}

/// Elmer's norm of the coupled variable after the last `ComputeChange`.
pub fn norm(sim: &Simulation) -> f64 {
    sim.vars.get("pot").unwrap().conv.norm
}
