//! Running `IncompressibleNsSolver` through `elmer-solver`'s `Simulation`
//! with its default backends (assembly, faer LU for `Direct`), the `elmer-fe`
//! basis and Gauss rules, and hand-built model sections standing in for
//! `elmer-sif`.

use std::path::{Path, PathBuf};

use elmer_core::{Dependent, DependentKind, NullLogger, ParamMap, Params, SerialComm, Value};
use elmer_physics_incompressible::{register, solver_defaults};
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

pub fn st(s: &str) -> Value {
    Value::Str(s.into())
}

pub fn logical(b: bool) -> Value {
    Value::Logical(b)
}

/// `Real Variable Coordinate <axis>` as a linear table through `(x, y)`.
pub fn table(axis: usize, pts: &[(f64, f64)]) -> Value {
    Value::Dependent(Dependent {
        variables: vec![format!("coordinate {axis}")],
        kind: DependentKind::Table {
            x: pts.iter().map(|p| p.0).collect(),
            y: pts.iter().map(|p| p.1).collect(),
            cubic: false,
        },
        rows: 1,
        cols: 1,
    })
}

/// `f` sampled at the coordinate `axis` of the nodes of boundary `tag`, as a
/// table that reproduces it exactly at those nodes (what MATC gives there).
pub fn profile(mesh: &TestMesh, tag: i32, axis: usize, f: impl Fn(f64) -> f64) -> Value {
    let mut xs: Vec<f64> = mesh
        .boundary
        .iter()
        .filter(|b| b.tag == tag)
        .flat_map(|b| b.nodes.iter().map(|&n| mesh.coords[n as usize][axis - 1]))
        .collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup();
    table(axis, &xs.iter().map(|&x| (x, f(x))).collect::<Vec<_>>())
}

pub fn run(mut model: Model, mesh: &TestMesh) -> Result<Simulation, SolverError> {
    let transient = model.simulation.string("simulation type").is_some_and(|s| s.eq_ignore_ascii_case("transient"));
    for s in model.solvers.iter_mut() {
        solver_defaults(s, mesh.dim, transient);
    }
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

/// An `IncompressibleNS` solver section with `extra` keywords and a direct
/// solve.
pub fn ns_solver(extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("NS-Vec")),
        ("Procedure", st("\"IncompressibleNS\" \"IncompressibleNSSolver\"")),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

/// A steady `Simulation` section.
pub fn steady(csys: &str) -> ParamMap {
    sec(&[
        ("Coordinate System", st(csys)),
        ("Simulation Type", st("Steady State")),
        ("Steady State Max Iterations", ints(&[1])),
    ])
}

/// A single-body model around `solver`.
pub fn model(simulation: ParamMap, solver: ParamMap, material: ParamMap, body_force: Option<ParamMap>, bcs: Vec<ParamMap>) -> Model {
    let mut body = vec![
        ("Target Bodies", ints(&[1])),
        ("Equation", ints(&[1])),
        ("Material", ints(&[1])),
    ];
    if body_force.is_some() {
        body.push(("Body Force", ints(&[1])));
    }
    Model {
        simulation,
        constants: ParamMap::new(),
        solvers: vec![solver],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&body)],
        materials: vec![material],
        body_forces: body_force.into_iter().collect(),
        boundary_conditions: bcs,
        ..Model::default()
    }
}

/// Node-ordered components (velocity, then pressure) from a pristine
/// ElmerSolver `.result`, as `node c1 c2 ...` lines.
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

/// The flow variable by node.
pub fn nodal(sim: &Simulation, n_nodes: usize) -> Vec<Vec<f64>> {
    let v = sim.vars.get("flow solution").unwrap();
    (0..n_nodes)
        .map(|n| {
            let s = v.slot(n as u32).unwrap() * v.dofs;
            v.values[s..s + v.dofs].to_vec()
        })
        .collect()
}

/// Elmer's norm of the flow variable after the last `ComputeChange`.
pub fn norm(sim: &Simulation) -> f64 {
    sim.vars.get("flow solution").unwrap().conv.norm
}
