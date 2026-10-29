//! A case run end to end with the fixtures installed: the real `.sif`
//! through `elmer-run`'s `SifFrontend` (`elmer-sif`), the mesh DB through
//! `elmer-mesh`, Elmer's mesh-load transforms (`elmer-mesh-levels`), then
//! `elmer-solver`'s `Simulation` with `elmer-run`'s solver registry plus
//! [`elmer_fixtures::registry::register`] and the case's functions
//! ([`elmer_fixtures::registry::install_live`]).
//!
//! This is `elmer_run::run_case` with the fixtures added; `elmer-run` does
//! not take them yet. `install_fe` keeps one element basis per process, so
//! each test binary runs cases of one coordinate dimension only.
#![allow(dead_code)]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_cli::logger::OutputConfig;
use elmer_cli::ElmerLogger;
use elmer_core::{ElemId, ElementRef, MeshView, ParamMap, Params, SerialComm, Value};
use elmer_mesh_levels::{extrude, mesh_levels, ExtrudeOptions, MeshLevelsError};
use elmer_physics_meshtools::motion::{self, MeshMotion, SharedCoords};
use elmer_run::output::is_result_output;
use elmer_run::{Frontend, LoadedCase, SifFrontend, SifMatc};
use elmer_solver::solver::procedure_of;
use elmer_solver::{finalize_reference, Model, Simulation, SolverError, SolverRegistry};

/// The case's files under `tests/data`, with any shared mesh, copied into a
/// fresh directory laid out as the `.sif`'s `Mesh DB` expects.
pub fn stage(case: &str, shared_mesh: Option<&str>) -> PathBuf {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let dst = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("elmer-fixtures-{case}"));
    let _ = std::fs::remove_dir_all(&dst);
    copy_dir(&data.join(case), &dst);
    if let Some(m) = shared_mesh {
        copy_dir(&data.join("meshes").join(m), &dst.join(m));
    }
    dst
}

pub fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for e in std::fs::read_dir(src).unwrap() {
        let e = e.unwrap();
        let to = dst.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &to);
        } else {
            std::fs::copy(e.path(), to).unwrap();
        }
    }
}

fn load_mesh(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    if simulation.get("coordinate scaling").is_some() {
        return Err("Coordinate Scaling (not needed by these cases)".into());
    }
    let dir = dir.to_str().ok_or("non-UTF-8 mesh path")?;
    Ok(Box::new(elmer_mesh::load(dir).map_err(|e| e.to_string())?))
}

/// The log, readable after the run.
#[derive(Clone, Default)]
pub struct Buf(Arc<Mutex<Vec<u8>>>);

impl Write for Buf {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Buf {
    pub fn text(&self) -> String {
        String::from_utf8_lossy(&self.0.lock().unwrap()).into_owned()
    }
}

/// What a run produced.
pub struct Outcome {
    /// `(solver id, ComputeNorm)` of every solver with a primary variable.
    pub norms: Vec<(usize, f64)>,
    /// The `TEST.PASSED` verdict.
    pub passed: Option<bool>,
    pub log: String,
}

impl Outcome {
    pub fn norm(&self, solver: usize) -> f64 {
        self.norms.iter().find(|n| n.0 == solver).unwrap_or_else(|| panic!("no norm for solver {solver}")).1
    }
}

/// `MapBodiesAndBCs` tags, which extrusion works on (`elmer-run`'s
/// `MappedTags`).
struct MappedTags<'a> {
    mesh: &'a dyn MeshView,
    model: &'a Model,
}

impl MeshView for MappedTags<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let mut el = self.mesh.element(e);
        el.tag = if (e.0 as usize) < self.mesh.n_bulk() {
            self.model.body_of_tag(el.tag).map_or(el.tag, |b| b as i32)
        } else {
            self.model.bc_of_tag(el.tag).map_or(0, |b| b as i32)
        };
        el
    }
}

/// The mesh with its coordinates in a [`SharedCoords`] the mesh mappers
/// move (`elmer-run`'s `MovableMesh`).
struct MovableMesh<'a> {
    mesh: &'a dyn MeshView,
    coords: Arc<SharedCoords>,
}

impl MeshView for MovableMesh<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords.get(node)
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        self.mesh.element(e)
    }
    fn element_edges(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_edges(e)
    }
    fn element_faces(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_faces(e)
    }
}

/// Restores the previously installed mesh motion on drop.
struct MotionGuard(Option<Arc<dyn MeshMotion>>);

impl Drop for MotionGuard {
    fn drop(&mut self) {
        motion::install(self.0.take());
    }
}

fn levels_err(e: MeshLevelsError) -> SolverError {
    match e {
        MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshLevelsError::Invalid(m) => SolverError::fatal("LoadModel", m),
    }
}

/// `elmer-run`'s procedure preflight, which also accepts the solver mode
/// sections (`Assembly Solver`, `Coupled Solver`) `elmer-solver` runs
/// without a `Procedure`.
fn preflight(model: &Model, registry: &SolverRegistry) -> Result<(), SolverError> {
    for (i, p) in model.solvers.iter().enumerate() {
        if p.0.is_empty() || p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never")) {
            continue;
        }
        if is_result_output(p) {
            return Err(SolverError::Unsupported(format!("Solver {}: VTU output (not composed here)", i + 1)));
        }
        let mode = ["coupled solver", "block solver", "assembly solver"].iter().any(|k| p.logical(k) == Some(true));
        match procedure_of(p) {
            None if mode => {}
            None => return Err(SolverError::Unsupported(format!("Solver {}: no Procedure", i + 1))),
            Some((lib, entry)) if registry.get(&lib, &entry).is_none() => {
                return Err(SolverError::Unsupported(format!("Solver {}: Procedure \"{lib}\" \"{entry}\"", i + 1)))
            }
            Some(_) => {}
        }
    }
    elmer_run::unported::check(model)
}

/// `AddEquationBasics` for `Equation = Heat Equation` without a
/// `Procedure`, as `elmer-run` completes it.
fn legacy_heat(model: &mut Model) {
    for p in &mut model.solvers {
        let heat = p.string("equation").is_some_and(|e| e.trim().eq_ignore_ascii_case("heat equation"));
        if heat && p.get("procedure").is_none() {
            if p.get("radiation solver").is_none() {
                p.insert("Radiation Solver", Value::Logical(true));
            }
            if p.get("variable").is_none() {
                p.insert("Variable", Value::Str("Temperature".into()));
            }
            p.insert("Procedure", Value::Str("HeatSolve HeatSolver".into()));
        }
    }
}

/// `StressSolver_Init0` (StressSolve.F90:292) names the variable of a
/// section that has none: `Variable = -dofs <dim> Displacement`. Neither
/// `elmer-run` nor `elmer-physics-elastic` does yet (the `tresca` case needs
/// it), so this stands in for that step.
fn stress_solver_variable(model: &mut Model, dim: usize) {
    for p in &mut model.solvers {
        let stress = procedure_of(p).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("stresssolve") && entry.eq_ignore_ascii_case("stresssolver")
        });
        if stress && p.get("variable").is_none() {
            p.insert("Variable", Value::Str("Displacement".into()));
            p.insert("Variable Dofs", Value::Integer(vec![dim as i32]));
        }
    }
}

/// A later `RUN` block may only change reference norms, as `elmer-run`
/// requires.
fn check_reloads(reloads: &[Model]) -> Result<(), SolverError> {
    const ALLOWED: [&str; 2] = ["reference norm", "reference norm tolerance"];
    for block in reloads {
        let others = [&block.equations, &block.bodies, &block.materials, &block.body_forces, &block.boundary_conditions];
        let changed = block.solvers.iter().any(|p| p.0.keys().any(|k| !ALLOWED.contains(&k.as_str())))
            || others.iter().any(|v| v.iter().any(|p| !p.0.is_empty()))
            || block.simulation.0.keys().any(|k| k != "solver input file");
        if changed {
            return Err(SolverError::Unsupported("RUN block changes more than reference norms".into()));
        }
    }
    Ok(())
}

/// Runs `sif` in `dir` as case `case` (the `fem/tests` directory name that
/// selects the fixtures). A `Post File` is dropped: VTU output is
/// `elmer-run`'s, not what these tests check.
pub fn run(case: &str, dir: &Path, sif: &str) -> Result<Outcome, SolverError> {
    let LoadedCase { mut model, reloads, mesh, case_dir, .. } =
        SifFrontend { mesh: load_mesh }.load(&dir.join(sif)).map_err(SolverError::Unsupported)?;
    model.solvers.retain(|p| !is_result_output(p));
    legacy_heat(&mut model);

    let split = mesh_levels(mesh.as_ref(), &model.simulation, &model.boundary_conditions).map_err(levels_err)?;
    let current: &dyn MeshView = split.as_ref().map_or(mesh.as_ref(), |m| m as &dyn MeshView);
    let extruded = match ExtrudeOptions::from_params(&model.simulation, current.dim(), None).map_err(levels_err)? {
        Some(opts) => {
            let m = extrude(&MappedTags { mesh: current, model: &model }, &opts).map_err(levels_err)?;
            for (i, body) in model.bodies.iter_mut().enumerate() {
                body.insert("Target Bodies", Value::Integer(vec![i as i32 + 1]));
            }
            for (i, bc) in model.boundary_conditions.iter_mut().enumerate() {
                bc.insert("Target Boundaries", Value::Integer(vec![i as i32 + 1]));
            }
            Some(m)
        }
        None => None,
    };
    let mesh: &dyn MeshView = match (&extruded, &split) {
        (Some(m), _) => m,
        (None, Some(m)) => m,
        (None, None) => mesh.as_ref(),
    };
    let coords = Arc::new(SharedCoords::new((0..mesh.n_nodes() as u32).map(|n| mesh.coord(n)).collect()));
    let shared: Arc<dyn MeshMotion> = coords.clone();
    let _motion = MotionGuard(motion::install(Some(shared)));
    let movable = MovableMesh { mesh, coords };
    let mesh: &dyn MeshView = &movable;

    let mut registry = elmer_run::solver_registry();
    elmer_fixtures::registry::register(case, &mut registry);
    elmer_fixtures::registry::init_sections(case, &mut model);
    preflight(&model, &registry)?;
    check_reloads(&reloads)?;
    let coords = elmer_run::coordinate_system(&model, mesh)?;
    stress_solver_variable(&mut model, coords.dim);
    elmer_run::install_fe(coords.dim)?;

    let buf = Buf::default();
    let log = ElmerLogger::new(OutputConfig::default(), Box::new(buf.clone()));
    let comm = SerialComm;
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(elmer_run::fe_backend(coords.dim));
    sim.set_matc(Rc::new(SifMatc));
    elmer_fixtures::registry::install_live(case, &mut sim, mesh);
    sim.setup(mesh, &registry, &log, &comm)?;
    sim.run(mesh, &log, &comm)?;
    for block in &reloads {
        sim.reload(block, mesh, &registry, &log, &comm)?;
        sim.run(mesh, &log, &comm)?;
    }
    let passed = finalize_reference(&sim.reference, &case_dir, &log, &comm)?;
    let norms = sim.slots.iter().filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm))).collect();
    Ok(Outcome { norms, passed, log: buf.text() })
}

/// `CompareToReferenceSolution`'s `Norm = ` value (`ES16.8`) of a pristine
/// run, against the native norm printed the same way.
pub fn assert_norm_prints(case: &str, got: f64, pristine: &str) {
    let printed = format!("{got:.8E}");
    let (m, e) = printed.split_once('E').unwrap();
    let e: i32 = e.parse().unwrap();
    let printed = format!("{m}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs());
    assert_eq!(printed, pristine, "{case}: native norm {got:e} prints {printed}, pristine printed {pristine}");
}
