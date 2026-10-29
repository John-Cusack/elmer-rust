//! Running the hooks: a recording logger, section builders, a case that holds
//! a hand-filled `Simulation`, and the full heat + SaveData registry.

use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use elmer_core::{FieldKind, Logger, NullLogger, ParamMap, SerialComm, Value};
use elmer_saveio::savedata::NamesMeta;
use elmer_solver::{Model, Simulation, SolverCtx, SolverError, SolverRegistry, SolverSlot, Variable};

use super::backends::{BackwardEuler, DenseLu, SimMesh, TestAssembly};
use super::meshio::Mesh;
use super::{TestBasis, TestQuadrature};

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

pub fn yes() -> Value {
    Value::Logical(true)
}

/// A real array `rows x cols`, row-major.
pub fn array(rows: usize, cols: usize, data: &[f64]) -> Value {
    Value::Real {
        rows,
        cols,
        data: data.to_vec(),
    }
}

/// Collects every message.
#[derive(Default)]
pub struct RecLogger {
    pub lines: Mutex<Vec<String>>,
}

impl Logger for RecLogger {
    fn info(&self, caller: &str, msg: &str, _level: u8) {
        self.lines.lock().unwrap().push(format!("{caller}: {msg}"));
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.lines.lock().unwrap().push(format!("WARNING:: {caller}: {msg}"));
    }
}

impl RecLogger {
    pub fn text(&self) -> String {
        self.lines.lock().unwrap().join("\n")
    }
}

pub fn names_meta() -> NamesMeta {
    NamesMeta {
        version: "26.2".into(),
        revision: Some("9f6af2f85".into()),
        compilation_date: Some("2026-09-23".into()),
        sif_name: Some("case.sif".into()),
        started_at: Some("2026/10/29 11:11:35".into()),
        ..NamesMeta::default()
    }
}

/// Installs the stand-in element library and the file header, once per test
/// binary.
pub fn install() {
    let _ = elmer_savedata::install_fe(|| (Box::new(TestBasis), Box::new(TestQuadrature)));
    let _ = elmer_physics_heat::install_fe(|| (Box::new(TestBasis), Box::new(TestQuadrature)));
    let _ = elmer_savedata::install_names_meta(names_meta);
}

/// A fresh directory under cargo's test tmp.
pub fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("elmer-savedata").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

pub fn dir_str(d: &std::path::Path) -> Value {
    Value::Str(d.display().to_string())
}

/// A hand-filled simulation: fields set from functions of the node, hooks
/// built from sections and visited one by one.
pub struct Case {
    pub sim: Simulation,
    pub mesh: SimMesh,
    pub log: RecLogger,
}

impl Case {
    pub fn new(mesh: Mesh, model: Model) -> Case {
        install();
        Case::bare(mesh, model)
    }

    /// Without installing the element library or the file header.
    pub fn bare(mesh: Mesh, model: Model) -> Case {
        Case {
            sim: Simulation::new(model),
            mesh: SimMesh::new(mesh),
            log: RecLogger::default(),
        }
    }

    /// A nodal field with every node in the system.
    pub fn field(&mut self, name: &str, dofs: usize, f: impl Fn([f64; 3]) -> Vec<f64>) {
        let n = self.mesh.mesh.coords.len();
        let perm: Arc<[u32]> = (0..n as u32).collect::<Vec<_>>().into();
        let values: Vec<f64> = self.mesh.mesh.coords.iter().flat_map(|&x| f(x)).collect();
        self.sim
            .vars
            .add(Variable::new(name, FieldKind::Nodal, dofs, Some(perm), values))
            .unwrap();
    }

    /// Adds a `Solver n` slot running the hook `entry` ("SaveScalars" or
    /// "SaveLine"); returns its index.
    pub fn hook(&mut self, entry: &str, section: ParamMap) -> Result<usize, SolverError> {
        let ix = self.sim.slots.len();
        let mut params = section;
        params.insert("Procedure", st(&format!("\"SaveData\" \"{entry}\"")));
        let mut reg = SolverRegistry::new();
        elmer_savedata::register(&mut reg);
        let factory = reg.get("SaveData", entry).expect("registered");
        let mut inst = factory(&params)?;
        self.sim.slots.push(SolverSlot::new(ix + 1, params)?);
        inst.init(&mut SolverCtx::new(&mut self.sim, ix, &self.mesh, &self.log, &SerialComm))?;
        self.sim.slots[ix].set_instance(inst);
        Ok(ix)
    }

    /// One visit, as the scheduler would make it.
    pub fn visit(&mut self, ix: usize, transient: bool) -> Result<(), SolverError> {
        self.sim
            .solver_activate(ix, 1.0, transient, &self.mesh, &self.log, &SerialComm)
    }

    /// One visit on a communicator of another size.
    pub fn visit_with(&mut self, ix: usize, comm: &dyn elmer_core::Communicator) -> Result<(), SolverError> {
        self.sim.solver_activate(ix, 1.0, false, &self.mesh, &self.log, comm)
    }
}

/// The full registry: heat plus both hooks, with the test backends.
pub fn heat_simulation(model: Model) -> Simulation {
    install();
    let mut sim = Simulation::new(model);
    sim.backends.assembly = Some(Box::new(TestAssembly));
    sim.backends.time = Some(Box::new(BackwardEuler));
    sim.backends
        .linear
        .register_direct("direct", |_| Ok(Box::new(DenseLu::default())));
    sim.backends
        .linear
        .register_direct("umfpack", |_| Ok(Box::new(DenseLu::default())));
    sim
}

pub fn full_registry() -> SolverRegistry {
    let mut reg = SolverRegistry::new();
    elmer_physics_heat::register(&mut reg);
    elmer_savedata::register(&mut reg);
    reg
}

pub fn quiet() -> NullLogger {
    NullLogger
}

/// The sections of a `.sif` as the solver framework reads them, through
/// `elmer-sif`'s parser.
pub fn sif_model(path: &std::path::Path) -> Model {
    let parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
    let m = parsed.sif.merged();
    // Elmer's `ListAddString` lowercases every string that is not a file
    // name, quoted or not; `elmer-sif` keeps the case of quoted text.
    let convert = |s: &elmer_sif::Section| -> ParamMap {
        let mut p = elmer_sif::section_to_param_map(s);
        for v in p.0.values_mut() {
            if let Value::Str(t) = v {
                *t = t.to_lowercase();
            }
        }
        p
    };
    let numbered = |v: &[elmer_sif::Section]| -> Vec<ParamMap> {
        let n = v.iter().map(|s| s.id as usize).max().unwrap_or(0);
        let mut out = vec![ParamMap::new(); n];
        for s in v {
            out[s.id as usize - 1] = convert(s);
        }
        out
    };
    let one = |s: &Option<elmer_sif::Section>| s.as_ref().map(convert).unwrap_or_default();
    Model {
        header: ParamMap::new(),
        simulation: one(&m.simulation),
        constants: one(&m.constants),
        solvers: numbered(&m.solvers),
        equations: numbered(&m.equations),
        bodies: numbered(&m.bodies),
        materials: numbered(&m.materials),
        body_forces: numbered(&m.body_forces),
        initial_conditions: numbered(&m.ics),
        boundary_conditions: numbered(&m.bcs),
    }
}
