//! Running `HeatSolver` through `elmer-solver`'s `Simulation` with its
//! default assembly (`elmer-assembly`) and time (`elmer-timeint`) backends,
//! HUTI BiCGStab with ILU(0) from `elmer-krylov`/`elmer-sparse`, a dense LU
//! standing in for UMFPACK, and `elmer-fe`'s Lagrange basis and Gauss rules.
//! Model sections are built by hand until `elmer-sif` merges.

use std::path::{Path, PathBuf};

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_krylov::{Huti, HutiConfig};
use elmer_physics_heat::{install_fe, register};
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

pub fn reals(v: &[f64]) -> Value {
    Value::Real {
        rows: v.len(),
        cols: 1,
        data: v.to_vec(),
    }
}

pub fn st(s: &str) -> Value {
    Value::Str(s.into())
}

pub fn install() {
    // Idempotent across tests in this binary.
    // cdim 3: the test meshes are flat (z = 0), which gives the same bits
    // as Elmer's model dimension 2.
    let _ = install_fe(|| {
        (
            Box::new(elmer_fe::LagrangeBasis::new(3)),
            Box::new(elmer_fe::ElmerQuadrature),
        )
    });
}

pub fn simulation(model: Model) -> (Simulation, SolverRegistry) {
    install();
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    let linear = &mut sim.backends.linear;
    linear.register_direct("umfpack", |_| Ok(Box::new(DenseLu::default())));
    linear.register_iterative("bicgstab", |p| {
        Ok(Box::new(Huti::<f64>::new(HutiConfig::from_params(p)?)))
    });
    linear.register_precond("ilu0", |_, a| Ok(Box::new(Ilu::ilu_n(a, 0)?)));
    (sim, reg)
}

pub fn run(model: Model, mesh: &SimMesh) -> Result<Simulation, SolverError> {
    let (mut sim, reg) = simulation(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// A `HeatSolve` solver section for `variable`, with `extra` keywords.
pub fn heat_solver(variable: &str, extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("heat equation")),
        ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
        ("Variable", st(variable)),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

/// `Linear System Solver = Direct` (UMFPACK in the Fortran runs).
pub fn direct() -> Vec<(&'static str, Value)> {
    vec![("Linear System Solver", st("Direct"))]
}

/// BiCGStab + ILU(0) to a tolerance well below the Fortran direct solve's
/// round-off in the compared norms.
pub fn bicgstab_ilu0() -> Vec<(&'static str, Value)> {
    vec![
        ("Linear System Solver", st("Iterative")),
        ("Linear System Iterative Method", st("BiCGStab")),
        ("Linear System Preconditioning", st("ILU0")),
        ("Linear System Max Iterations", ints(&[2000])),
        ("Linear System Convergence Tolerance", Value::real(1e-14)),
        ("Linear System Residual Output", ints(&[0])),
    ]
}

/// The solver variable by node.
pub fn nodal(sim: &Simulation, var: &str, n_nodes: usize) -> Vec<f64> {
    let v = sim.vars.get(var).unwrap();
    (0..n_nodes)
        .map(|n| v.values[v.slot(n as u32).unwrap()])
        .collect()
}

/// Elmer's norm of the solver variable after the last `ComputeChange`.
pub fn norm(sim: &Simulation, var: &str) -> f64 {
    sim.vars.get(var).unwrap().conv.norm
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
