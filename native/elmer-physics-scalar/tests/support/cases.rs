//! Running the solvers through `elmer-solver`'s `Simulation` with its default
//! backends (assembly, time, HUTI, ILU, direct) and the test-local basis
//! standing in for `elmer-fe`. Model sections are built by hand until
//! `elmer-sif` is wired in.

use std::path::{Path, PathBuf};

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_scalar::{install_fe, register};
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry};

use super::simmesh::SimMesh;
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

pub fn st(s: &str) -> Value {
    Value::Str(s.into())
}

pub fn real(x: f64) -> Value {
    Value::real(x)
}

pub fn install() {
    // Idempotent across tests in this binary.
    let _ = install_fe(|| (Box::new(TestBasis), Box::new(TestQuadrature)));
}

pub fn simulation(model: Model) -> (Simulation, SolverRegistry) {
    install();
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    (Simulation::new(model), reg)
}

pub fn run(model: Model, mesh: &SimMesh) -> Result<Simulation, SolverError> {
    let (mut sim, reg) = simulation(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// Columns of a `node v1 v2 ...` file written from the pristine ElmerSolver's
/// ASCII `.result` (node order, full precision): `columns[k][node - 1]`.
pub fn read_columns(path: &Path) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut cols: Vec<Vec<f64>> = Vec::new();
    for (i, l) in text.lines().enumerate() {
        let t: Vec<&str> = l.split_whitespace().collect();
        assert_eq!(t[0].parse::<usize>().unwrap(), i + 1);
        cols.resize(t.len() - 1, Vec::new());
        for (c, s) in cols.iter_mut().zip(&t[1..]) {
            c.push(s.parse().unwrap());
        }
    }
    cols
}

/// A solver variable component by node (`dofs` interleaved, component `k`).
pub fn nodal(sim: &Simulation, var: &str, k: usize, n_nodes: usize) -> Vec<f64> {
    let v = sim.vars.get(var).unwrap_or_else(|| panic!("no variable {var}"));
    (0..n_nodes)
        .map(|n| v.values[v.slot(n as u32).unwrap() * v.dofs + k])
        .collect()
}

pub fn norm(sim: &Simulation, var: &str) -> f64 {
    sim.vars.get(var).unwrap().conv.norm
}

/// Largest `|a - b|` relative to the largest `|b|`.
pub fn rel_max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    a.iter().zip(b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs())) / scale
}

pub fn assert_rel(got: f64, want: f64, tol: f64) {
    let rel = (got - want).abs() / want.abs();
    assert!(rel < tol, "{got:e} vs {want:e}: relative error {rel:e}");
}

/// Applies `Coordinate Scaling`.
pub fn scale_coords(mesh: &mut super::meshio::Mesh, s: f64) {
    for c in &mut mesh.coords {
        for x in c.iter_mut() {
            *x *= s;
        }
    }
}
