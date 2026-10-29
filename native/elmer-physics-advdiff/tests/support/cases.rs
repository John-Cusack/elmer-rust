//! Running the solvers through `elmer-solver`'s `Simulation` with its default
//! assembly (`elmer-assembly`) and time (`elmer-timeint`) backends, HUTI
//! BiCGStab with ILU(n) from `elmer-krylov`/`elmer-sparse`, and a dense LU
//! standing in for UMFPACK. Model sections are built by hand.

use std::path::{Path, PathBuf};

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_krylov::{Huti, HutiConfig};
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry};
use elmer_sparse::Ilu;

use super::backends::{DenseLu, SimMesh};

pub fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
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

pub fn simulation(model: Model) -> (Simulation, SolverRegistry) {
    let mut reg = SolverRegistry::new();
    elmer_physics_advdiff::register(&mut reg);
    elmer_physics_heat::register(&mut reg);
    // Idempotent across tests in this binary.
    let _ = elmer_physics_heat::install_fe(|| {
        (Box::new(LagrangeBasis::new(2)), Box::new(ElmerQuadrature))
    });
    let mut sim = Simulation::new(model);
    let linear = &mut sim.backends.linear;
    linear.register_direct("umfpack", |_| Ok(Box::new(DenseLu::default())));
    linear.register_iterative("bicgstab", |p| {
        Ok(Box::new(Huti::<f64>::new(HutiConfig::from_params(p)?)))
    });
    linear.register_precond("ilu0", |_, a| Ok(Box::new(Ilu::ilu_n(a, 0)?)));
    linear.register_precond("ilu1", |_, a| Ok(Box::new(Ilu::ilu_n(a, 1)?)));
    (sim, reg)
}

pub fn run(model: Model, mesh: &SimMesh) -> Result<Simulation, SolverError> {
    let (mut sim, reg) = simulation(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// The solver variable by node.
pub fn nodal(sim: &Simulation, var: &str, n_nodes: usize) -> Vec<f64> {
    let v = sim.vars.get(var).unwrap();
    (0..n_nodes)
        .map(|n| v.values[v.slot(n as u32).unwrap()])
        .collect()
}

/// A node-ordered solution written from the pristine ElmerSolver's ASCII
/// `.result` file: one `node value` line per node, full precision.
pub fn read_solution(path: &Path) -> Vec<f64> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .enumerate()
        .map(|(i, l)| {
            let t: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(t[0].parse::<usize>().unwrap(), i + 1);
            t[1].parse().unwrap()
        })
        .collect()
}

pub fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .fold(0.0, |m, (x, y)| f64::max(m, (x - y).abs()))
}

/// Elmer's norm: `sqrt(sum(x^2) / n)`.
pub fn rms(x: &[f64]) -> f64 {
    (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt()
}
