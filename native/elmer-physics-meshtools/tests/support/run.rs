//! Running mappers (and the HeatSolve behind them) through `elmer-solver`.

use std::path::PathBuf;
use std::rc::Rc;

use elmer_core::{Dependent, DependentKind, FieldKind, NullLogger, ParamMap, SerialComm, Value};
use elmer_krylov::{Huti, HutiConfig};
use elmer_solver::{FeBackend, Model, Simulation, SolverError, SolverRegistry, TextMatc};
use elmer_sparse::Ilu;

use super::simmesh::SimMesh;

pub fn data(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(rel)
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

pub fn lg(b: bool) -> Value {
    Value::Logical(b)
}

pub fn reals(v: &[f64]) -> Value {
    Value::Real { rows: v.len(), cols: 1, data: v.to_vec() }
}

/// `Real MATC "expr"` depending on `vars`.
pub fn matc(vars: &[&str], expr: &str) -> Value {
    Value::Dependent(Dependent {
        variables: vars.iter().map(|s| s.to_string()).collect(),
        kind: DependentKind::Matc(expr.into()),
        rows: 1,
        cols: 1,
    })
}

/// A simulation with the real element layer, MATC, and BiCGStab with ILU(0..2).
pub fn simulation(model: Model) -> (Simulation, SolverRegistry) {
    let _ = elmer_physics_heat::install_fe(|| {
        (Box::new(elmer_fe::LagrangeBasis::default()), Box::new(elmer_fe::ElmerQuadrature))
    });
    let mut reg = SolverRegistry::new();
    elmer_physics_heat::register(&mut reg);
    elmer_physics_meshtools::register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(elmer_fe::LagrangeBasis::default()),
        quadrature: Box::new(elmer_fe::ElmerQuadrature),
    });
    let linear = &mut sim.backends.linear;
    linear.register_iterative("bicgstab", |p| Ok(Box::new(Huti::<f64>::new(HutiConfig::from_params(p)?))));
    linear.register_precond("ilu0", |_, a| Ok(Box::new(Ilu::ilu_n(a, 0)?)));
    linear.register_precond("ilu1", |_, a| Ok(Box::new(Ilu::ilu_n(a, 1)?)));
    (sim, reg)
}

pub fn run(model: Model, mesh: &SimMesh) -> Result<Simulation, SolverError> {
    let (mut sim, reg) = simulation(model);
    let prev = elmer_physics_meshtools::motion::install(Some(mesh.coords.clone()));
    let r = (|| {
        sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
        sim.run(mesh, &NullLogger, &SerialComm)?;
        Ok(())
    })();
    elmer_physics_meshtools::motion::install(prev);
    r.map(|()| sim)
}

/// A pristine VTU dump (`# columns: ...` header, one row per node).
pub struct Reference {
    cols: Vec<String>,
    rows: Vec<Vec<f64>>,
}

impl Reference {
    pub fn read(file: &str) -> Self {
        let text = std::fs::read_to_string(data(file)).unwrap();
        let mut cols = Vec::new();
        let mut rows = Vec::new();
        for line in text.lines() {
            if let Some(c) = line.strip_prefix("# columns: ") {
                cols = c.split(" | ").map(|s| s.trim().to_string()).collect();
            } else if !line.trim().is_empty() {
                rows.push(line.split_whitespace().map(|x| x.parse().unwrap()).collect());
            }
        }
        Self { cols, rows }
    }

    pub fn col(&self, name: &str) -> Vec<f64> {
        let i = self.cols.iter().position(|c| c == name).unwrap_or_else(|| panic!("no column {name}: {:?}", self.cols));
        self.rows.iter().map(|r| r[i]).collect()
    }
}

pub fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).fold(0.0, |m, (x, y)| f64::max(m, (x - y).abs()))
}

/// A nodal variable by node.
pub fn nodal(sim: &Simulation, var: &str, n_nodes: usize) -> Vec<f64> {
    let v = sim.vars.get(var).unwrap_or_else(|| panic!("no variable {var}"));
    assert_eq!(v.kind, FieldKind::Nodal);
    // Nodes outside the solver's bodies have no dof; the VTU writer shows 0 there.
    (0..n_nodes).map(|n| v.slot(n as u32).map_or(0.0, |s| v.values[s])).collect()
}

pub fn coords_of(mesh: &SimMesh, k: usize) -> Vec<f64> {
    mesh.coords.snapshot().iter().map(|c| c[k]).collect()
}
