//! One case through the native path: preflight, `Simulation` setup and run,
//! `Post File` output, then `CompareToReferenceSolution`'s `TEST.PASSED`.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_core::{normalize_key, Communicator, DependentKind, Logger, MeshView, ParamMap, Params, SerialComm, Value};
use elmer_fe::{CoordinateKind, CoordinateSystem};
use elmer_solver::{
    coupled::solver_for_mode, finalize_reference, solver::procedure_of, Model, RunReport, Simulation, Solver, SolverCtx,
    SolverError, SolverRegistry,
};

use crate::compose::{coordinate_system, fe_backend, install_fe, solver_registry_for};
use crate::matc::SifMatc;
use elmer_expr::MatcEngine;
use elmer_pdofs::PDofsBackend;
use crate::motion::MovableMesh;
use crate::output::{is_result_output, prepare_vtu_output, MappedTags};
use elmer_mortar::{generate_periodic_projectors, PeriodicAssembly, ProjectorEnv};
use elmer_mesh_levels::{extrude, mesh_levels, ExtrudeOptions, Mesh as LevelsMesh, MeshLevelsError};
use elmer_mesh_ops::{create_discont_mesh, increase_element_order, read_target_names, MeshOpsError};

/// A parsed case and its mesh.
pub struct Case<'a> {
    pub model: Model,
    /// Further `RUN` blocks (`Simulation::reload`), each with only its own
    /// entries; empty for a single simulation.
    pub reloads: Vec<Model>,
    pub mesh: &'a dyn MeshView,
    /// Where `TEST.PASSED` goes and a relative `Output Directory` starts.
    pub case_dir: PathBuf,
    /// The `Mesh DB` directory, the default VTU output directory.
    pub mesh_dir: PathBuf,
}

/// What a finished run produced.
#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    pub report: RunReport,
    /// `(solver id, ComputeNorm)` of every solver with a primary variable.
    pub norms: Vec<(usize, f64)>,
    /// The `TEST.PASSED` verdict; `None` when no solver has a reference.
    pub passed: Option<bool>,
    /// Output files written.
    pub written: Vec<PathBuf>,
}

/// The procedure `AddEquationBasics` gives a solver section whose
/// `Equation` is the historical solver `equation`, for the ones this build
/// completes (see [`legacy_solvers`]).
pub fn legacy_procedure(equation: &str) -> Option<(&'static str, &'static str)> {
    match equation.trim().to_ascii_lowercase().as_str() {
        "heat equation" => Some(("HeatSolve", "HeatSolver")),
        "stress analysis" => Some(("StressSolve", "StressSolver")),
        "navier-stokes" => Some(("FlowSolve", "FlowSolver")),
        _ => None,
    }
}

/// Ids of the solvers other solvers call as slaves (`Pre Solvers`,
/// `Post Solvers`, their `Nonlinear` forms, `Slave Solvers`): they run
/// whatever their own `Exec Solver` says.
fn slave_solvers(model: &Model) -> Vec<usize> {
    const KEYS: [&str; 5] = ["pre solvers", "post solvers", "nonlinear pre solvers", "nonlinear post solvers", "slave solvers"];
    let mut ids = Vec::new();
    for p in &model.solvers {
        for key in KEYS {
            if let Some(Value::Integer(v)) = p.get(key) {
                ids.extend(v.iter().filter(|&&i| i > 0).map(|&i| i as usize));
            }
        }
    }
    ids
}

/// A SaveData section with `Exec Solver = Never` never runs and creates no
/// field another solver reads, but `elmer-solver` builds every slot and
/// `elmer-savedata` refuses some operators (`cpu time`) when built. Such a
/// section is blanked, as an unset solver id is (MortarPoisson2D's timing
/// SaveScalars).
fn drop_never_savedata(model: &mut Model) {
    let slaves = slave_solvers(model);
    for (i, p) in model.solvers.iter_mut().enumerate() {
        let never = p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"));
        let savedata = procedure_of(p).is_some_and(|(lib, _)| lib.eq_ignore_ascii_case("savedata"));
        if never && savedata && !slaves.contains(&(i + 1)) {
            *p = ParamMap::new();
        }
    }
}

/// The historical solvers `AddEquationBasics` (MainUtils.F90:1221-1305)
/// completes from the solver's `Equation` name when it has no `Procedure`,
/// for the ones whose procedure this build has. Others stay without a
/// procedure and are refused by [`preflight`].
fn legacy_solvers(model: &mut Model, coords: CoordinateSystem) {
    let cdim = coords.dim;
    for params in &mut model.solvers {
        // `ElasticSolver_Init` (StressSolve delegates to it): without a
        // `Variable`, the displacement of the coordinate dimension.
        let elastic = procedure_of(params).is_some_and(|(lib, entry)| {
            [("stresssolve", "stresssolver"), ("elasticsolve", "elasticsolver")]
                .iter()
                .any(|(l, e)| lib.eq_ignore_ascii_case(l) && entry.eq_ignore_ascii_case(e))
        });
        let mixed = params.logical("mixed formulation") == Some(true) && params.logical("neo-hookean material") == Some(true);
        if elastic && !mixed && params.get("variable").is_none() {
            params.insert("Variable", Value::Str("Displacement".into()));
            params.insert("Variable Dofs", Value::Integer(vec![cdim as i32]));
        }
        // `HeatSolver_init` (HeatSolve.F90:167): `Variable = Temperature`.
        let heat = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("heatsolve") && entry.eq_ignore_ascii_case("heatsolver")
        });
        if heat && params.get("variable").is_none() {
            params.insert("Variable", Value::Str("Temperature".into()));
        }
        // `IncompressibleNSSolver_Init` (IncompressibleNS.F90:2390).
        let ns = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("incompressiblens") && entry.eq_ignore_ascii_case("incompressiblenssolver")
        });
        if ns && params.get("variable").is_none() {
            let v = if cdim == 2 { "Flow Solution[Velocity:2 Pressure:1]" } else { "Flow Solution[Velocity:3 Pressure:1]" };
            params.insert("Variable", Value::Str(v.into()));
        }
        if params.get("procedure").is_some() {
            continue;
        }
        let name = params.string("equation").map(|s| s.trim().to_ascii_lowercase());
        let (var, dofs, procedure) = match name.as_deref() {
            Some("heat equation") => {
                if params.get("radiation solver").is_none() {
                    params.insert("Radiation Solver", Value::Logical(true));
                }
                ("Temperature", 1, "HeatSolve HeatSolver")
            }
            Some("stress analysis") => ("Displacement", cdim, "StressSolve StressSolver"),
            // `dofs` is reset to 1: the components are in the name. The 1e-6
            // start value is applied by elmer-physics-flow.
            Some("navier-stokes") => {
                let velocity = cdim + usize::from(coords.kind == CoordinateKind::CylindricSymmetric);
                if velocity == 3 {
                    ("Flow Solution[Velocity:3 Pressure:1]", 1, "FlowSolve FlowSolver")
                } else {
                    ("Flow Solution[Velocity:2 Pressure:1]", 1, "FlowSolve FlowSolver")
                }
            }
            _ => continue,
        };
        if params.get("variable").is_none() {
            params.insert("Variable", Value::Str(var.into()));
            if dofs > 1 {
                params.insert("Variable Dofs", Value::Integer(vec![dofs as i32]));
            }
        }
        params.insert("Procedure", Value::Str(procedure.into()));
    }
}

/// Refuses solver sections this build cannot run: `Simulation` skips a
/// procedure missing from the registry, which would drop a solver silently.
/// `Exec Solver = Never` sections are exempt. Then refuses the model
/// features no provider implements ([`crate::unported`]).
fn preflight(model: &Model, registry: &SolverRegistry) -> Result<(), SolverError> {
    let slaves = slave_solvers(model);
    for (i, params) in model.solvers.iter().enumerate() {
        // A section that never runs needs no implementation (heateq_bdf2's
        // SaveScalars), unless another solver calls it as a slave.
        let never = params.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"));
        if params.0.is_empty() || (never && !slaves.contains(&(i + 1))) {
            continue;
        }
        // Checked and served by [`crate::output`].
        if is_result_output(params) {
            continue;
        }
        let Some((lib, entry)) = procedure_of(params) else {
            // `Coupled Solver` / `Assembly Solver` sections run without one
            // (`elmer-solver`'s built-in coupled and assembly solvers).
            if solver_for_mode(params, registry)?.is_some() {
                continue;
            }
            return Err(SolverError::Unsupported(format!("Solver {}: no Procedure", i + 1)));
        };
        if registry.get(&lib, &entry).is_none() {
            return Err(SolverError::Unsupported(format!("Solver {}: Procedure \"{lib}\" \"{entry}\"", i + 1)));
        }
    }
    check_calcfields_scalars(model)?;
    check_enclosed_flow(model)?;
    crate::unported::check(model)
}

/// FlowSolve with every boundary condition fixing all velocity components
/// and none setting the pressure: the pressure is defined only up to a
/// constant, which Elmer's solve tolerates. The composition diverges on it
/// (CavityLid: norm 5.5e20 against 0.2148, with or without BC-numbering
/// order, stabilization or Newton).
fn check_enclosed_flow(model: &Model) -> Result<(), SolverError> {
    let flow = model.solvers.iter().any(|p| {
        procedure_of(p).is_some_and(|(lib, _)| lib.eq_ignore_ascii_case("flowsolve"))
            && !p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
    });
    if !flow || model.boundary_conditions.is_empty() {
        return Ok(());
    }
    let dim = if model.boundary_conditions.iter().any(|bc| bc.get("velocity 3").is_some()) { 3 } else { 2 };
    let enclosed = model.boundary_conditions.iter().all(|bc| {
        (1..=dim).all(|i| bc.get(&format!("velocity {i}")).is_some())
            && !bc.0.keys().any(|k| k.contains("pressure") || k.starts_with("normal-tangential"))
    });
    if enclosed {
        return Err(SolverError::Unsupported(
            "FlowSolve on an enclosed domain (velocity fixed on every boundary, pressure only up to a constant)".into(),
        ));
    }
    Ok(())
}

/// `MagnetoDynamicsCalcFields` also stores its integral results (field
/// energies, losses) as `res:` Simulation values, which a SaveScalars
/// section without its own `Variable 1` saves. `elmer-physics-magneto`
/// computes the nodal fields only, so such a SaveScalars would report zeros
/// (mgdyn_airgap's Solver 3).
fn check_calcfields_scalars(model: &Model) -> Result<(), SolverError> {
    let runs = |p: &&ParamMap| !p.0.is_empty() && !p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"));
    let entry = |p: &ParamMap, name: &str| procedure_of(p).is_some_and(|(_, e)| e.eq_ignore_ascii_case(name));
    let calc = model.solvers.iter().filter(runs).any(|p| entry(p, "magnetodynamicscalcfields"));
    let default_scalars =
        model.solvers.iter().filter(runs).any(|p| entry(p, "savescalars") && p.get("variable 1").is_none());
    if calc && default_scalars {
        return Err(SolverError::Unsupported(
            "SaveScalars of MagnetoDynamicsCalcFields' integral results (res: energies) is not implemented".into(),
        ));
    }
    Ok(())
}

/// HeatSolve convects relative to the mesh (`U - MU` with the `Mesh
/// Velocity` variable, HeatSolve.F90); `elmer-physics-heat` uses the
/// material velocity alone, so a convected heat equation on a moving mesh
/// is refused once the variables exist.
fn check_mesh_velocity(sim: &Simulation) -> Result<(), SolverError> {
    if sim.vars.id("mesh velocity").is_none() {
        return Ok(());
    }
    let heat = sim.model.solvers.iter().any(|p| {
        procedure_of(p).is_some_and(|(lib, _)| lib.eq_ignore_ascii_case("heatsolve"))
            && !p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
    });
    let convected =
        sim.model.equations.iter().any(|p| p.string("convection").is_some_and(|c| !c.eq_ignore_ascii_case("none")));
    if heat && convected {
        return Err(SolverError::Unsupported(
            "HeatSolve: convection relative to the Mesh Velocity (moving mesh) is not implemented".into(),
        ));
    }
    Ok(())
}

/// Section-wise merge of `block` into `model`, as `Simulation::reload` does
/// (`ReloadInputFile`): entries replace same-name ones, new sections append.
fn merged_with(model: &Model, block: &Model) -> Model {
    fn merge(into: &mut ParamMap, from: &ParamMap) {
        for (k, v) in &from.0 {
            into.0.insert(k.clone(), v.clone());
        }
    }
    fn merge_vec(into: &mut Vec<ParamMap>, from: &[ParamMap]) {
        if into.len() < from.len() {
            into.resize_with(from.len(), ParamMap::new);
        }
        for (a, b) in into.iter_mut().zip(from) {
            merge(a, b);
        }
    }
    let mut m = model.clone();
    merge(&mut m.simulation, &block.simulation);
    merge(&mut m.constants, &block.constants);
    merge_vec(&mut m.solvers, &block.solvers);
    merge_vec(&mut m.equations, &block.equations);
    merge_vec(&mut m.bodies, &block.bodies);
    merge_vec(&mut m.materials, &block.materials);
    merge_vec(&mut m.body_forces, &block.body_forces);
    merge_vec(&mut m.initial_conditions, &block.initial_conditions);
    merge_vec(&mut m.boundary_conditions, &block.boundary_conditions);
    m
}

/// Refuses, before anything runs, what a later `RUN` block would make
/// unsupported. Fortran solvers read their keywords on every visit; native
/// solvers are rebuilt after the reload ([`rebuild_changed`]). Not
/// composed: Simulation and Constants changes, a new solver, and changes to
/// the VTU output and mortar sections `run_case` installed its own
/// instances for.
fn check_reloads(model: &Model, reloads: &[Model], registry: &SolverRegistry) -> Result<(), SolverError> {
    const ALLOWED: [&str; 2] = ["reference norm", "reference norm tolerance"];
    let mut merged = model.clone();
    for (k, block) in reloads.iter().enumerate() {
        let run = k + 2;
        let refuse = |what: String| {
            Err(SolverError::Unsupported(format!(
                "RUN block {run} changes {what} (re-reading changed keywords after a reload is not composed)"
            )))
        };
        // HelmholtzFEM's second block (`Simulation :: Angular Frequency`)
        // runs at the new frequency but misses its reference (1.28e-2 vs
        // 9.90e-3); Simulation and Constants changes stay refused.
        if let Some(key) = block.simulation.0.keys().find(|k| *k != "solver input file") {
            return refuse(format!("Simulation {key}"));
        }
        if let Some(key) = block.constants.0.keys().next() {
            return refuse(format!("Constants {key}"));
        }
        for (i, params) in block.solvers.iter().enumerate() {
            if params.0.keys().all(|k| ALLOWED.contains(&k.as_str())) {
                continue;
            }
            if i >= model.solvers.len() {
                return refuse(format!("Solver {} (a new solver)", i + 1));
            }
            let own = &merged.solvers[i];
            if is_result_output(own) || crate::mortar::applies(own) {
                let key = params.0.keys().find(|k| !ALLOWED.contains(&k.as_str())).cloned().unwrap_or_default();
                return refuse(format!("Solver {} {key}", i + 1));
            }
        }
        merged = merged_with(&merged, block);
        preflight(&merged, registry)?;
    }
    Ok(())
}

/// After `Simulation::reload(block)`: a solver whose own section `block`
/// changes, or every solver when `block` changes any other section, gets a
/// new instance from the registry, so a solver that reads its sections when
/// built (SaveScalars, HelmholtzSolve's frequency, ...) sees the new
/// keywords as the Fortran does on every visit. Its `init` runs on the first
/// solve. The VTU output and mortar instances `run_case` installed stay.
fn rebuild_changed(sim: &mut Simulation, block: &Model, registry: &SolverRegistry) -> Result<(), SolverError> {
    const ALLOWED: [&str; 2] = ["reference norm", "reference norm tolerance"];
    let others = [&block.simulation, &block.constants]
        .into_iter()
        .chain(&block.equations)
        .chain(&block.bodies)
        .chain(&block.materials)
        .chain(&block.body_forces)
        .chain(&block.initial_conditions)
        .chain(&block.boundary_conditions)
        .any(|p| p.0.keys().any(|k| k != "solver input file"));
    for slot in &mut sim.slots {
        let own = block.solvers.get(slot.id - 1).is_some_and(|p| p.0.keys().any(|k| !ALLOWED.contains(&k.as_str())));
        let params = &sim.model.solvers[slot.id - 1];
        if !(own || others) || is_result_output(params) || crate::mortar::applies(params) {
            continue;
        }
        let Some(factory) = procedure_of(params).and_then(|(lib, entry)| registry.get(&lib, &entry)) else {
            continue; // Built-in coupled / assembly solvers read their section live.
        };
        slot.set_instance(Box::new(Reinit { inner: factory(params)?, ready: false }));
    }
    Ok(())
}

/// A solver built after setup: `init` on its first visit.
struct Reinit {
    inner: Box<dyn Solver>,
    ready: bool,
}

impl Solver for Reinit {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if !self.ready {
            self.inner.init(ctx)?;
            self.ready = true;
        }
        self.inner.solve(ctx, dt, transient)
    }
}

/// Inputs no provider applies yet: restart files, mesh-load transforms
/// beyond `Coordinate Scaling`, `Mesh Levels` and extrusion, and element
/// definitions other than the Lagrange (`n:1`) basis `elmer-fe` evaluates.
fn check_mesh_and_elements(model: &Model) -> Result<(), SolverError> {
    let sim = &model.simulation;
    // `Restart File` sets initial fields from a .result file; nothing reads
    // it yet, and running without it is a different case.
    if sim.get("restart file").is_some() {
        return Err(SolverError::Unsupported("Simulation: Restart File (restart from a .result file)".into()));
    }
    if sim.get("coordinate transformation").is_some() {
        return Err(SolverError::Unsupported("Simulation: Coordinate Transformation".into()));
    }
    match sim.get("coordinate mapping") {
        None => {}
        Some(Value::Integer(v)) if v[..] == [1, 2, 3] => {}
        Some(v) => return Err(SolverError::Unsupported(format!("Simulation: Coordinate Mapping {v:?}"))),
    }
    // Element definitions and the solvers they reach (`NonNodalElements`,
    // MeshLoad.F90): the solver's own `Element`; an equation's `Element`
    // (every solver) or `Element{s}`; a body's `Solver s: Element`. Each is
    // judged as if it were the solver's own.
    let runs = |p: &ParamMap| {
        !p.0.is_empty()
            && !p.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
            && !is_result_output(p)
    };
    let solver_ix = |s: &str| s.trim().parse::<usize>().ok().filter(|&s| s >= 1 && s <= model.solvers.len()).map(|s| s - 1);
    let mut defs: Vec<(String, &str, Vec<usize>)> = Vec::new();
    for (i, p) in model.solvers.iter().enumerate() {
        if let Some(e) = p.string("element") {
            defs.push((format!("Solver section {}", i + 1), e, vec![i]));
        }
    }
    for (kind, sections) in [("Equation", &model.equations), ("Body", &model.bodies)] {
        for (i, p) in sections.iter().enumerate() {
            for (key, v) in &p.0 {
                let Value::Str(e) = v else { continue };
                let targets = match (kind, key.as_str()) {
                    ("Equation", "element") => (0..model.solvers.len()).collect(),
                    ("Equation", k) => match k.strip_prefix("element{").and_then(|r| r.strip_suffix('}')) {
                        Some(s) => solver_ix(s).into_iter().collect(),
                        None => continue,
                    },
                    (_, k) => match k.strip_prefix("solver ").and_then(|r| r.strip_suffix(": element")) {
                        Some(s) => solver_ix(s).into_iter().collect(),
                        None => continue,
                    },
                };
                defs.push((format!("{kind} section {}: {key}", i + 1), e, targets));
            }
        }
    }
    for (place, e, targets) in defs {
        for s in targets.into_iter().filter(|&s| runs(&model.solvers[s])) {
            let mut solver = model.solvers[s].clone();
            solver.insert("Element", Value::Str(e.to_string()));
            if !native_element(&solver) {
                return Err(SolverError::Unsupported(format!(
                    "{place}: Element = {e} for Solver {} (p-, bubble, edge or DG elements)",
                    s + 1
                )));
            }
        }
    }
    Ok(())
}

/// Whether the section's `Element` (if any) is one a provider evaluates:
/// the Lagrange `n:1`; `p:1` for ThermoElectricSolver, whose `b:N`
/// bubbles stay refused (the solver condenses them, while Elmer keeps them
/// in the global system and its solution norm, the `Reference Norm`); and
/// any definition for the solvers that read their DOFs from `elmer-pdofs`
/// (`ctx.solver_dofs()`) and refuse what they do not port: StressSolve
/// (#206), StatElecSolve (#208), ModelPDE's AdvDiffSolver and the test-local `Poisson` / `PoissonML`
/// `PoissonSolver`, `ModelPDEsplitted`, `SingleElements`, `BenchmarkShell`
/// and `Stokes1` ports.
pub fn native_element(params: &ParamMap) -> bool {
    let Some(e) = params.string("element") else { return true };
    if normalize_key(e) == "n:1" {
        return true;
    }
    let Some((lib, entry)) = procedure_of(params) else { return false };
    let lib = lib.trim_start_matches("./");
    const P_DOFS: [(&str, &str); 10] = [
        ("benchmarkshell", "shellsolver"),
        ("statelecsolve", "statelecsolver"),
        ("statelecsolvevec", "statelecsolver"),
        ("stresssolve", "stresssolver"),
        ("stokes1", "stokessolver"),
        ("poisson", "poissonsolver"),
        ("poissonml", "poissonsolver"),
        ("modelpde", "advdiffsolver"),
        ("modelpdesplitted", "advdiffsolver"),
        ("singleelements", "singleelements"),
    ];
    if P_DOFS.iter().any(|(l, n)| lib.eq_ignore_ascii_case(l) && entry.eq_ignore_ascii_case(n)) {
        return true;
    }
    lib.eq_ignore_ascii_case("thermoelectricsolver")
        && elmer_physics_thermoelectric::parse_element(e).is_ok_and(|nb| nb == 0)
}

/// `ListGetFun(Simulation, 'Extruded Mesh Density' | '1D Mesh Density', t)`:
/// a table (linear, or cubic for `Real Cubic`) or MATC function of the
/// position `t` on the unit segment.
fn mesh_density(sim: &ParamMap, t: f64) -> Result<f64, String> {
    let v = sim.get("extruded mesh density").or_else(|| sim.get("1d mesh density"));
    let Some(Value::Dependent(d)) = v else {
        return Err("Extruded Mesh Density is not a function".into());
    };
    match &d.kind {
        DependentKind::Table { x, y, cubic } => {
            let spline = cubic.then(|| elmer_expr::curve::cubic_spline(x, y, false));
            Ok(elmer_expr::curve::interpolate_curve(x, y, t, spline.as_deref()))
        }
        DependentKind::Matc(expr) => SifMatc.eval(expr, &[t])?.first().copied().ok_or_else(|| format!("MATC {expr}: no value")),
        _ => Err("Extruded Mesh Density: only tables and MATC functions are evaluated".into()),
    }
}

/// Elmer's mesh-load transforms, in Elmer's order: `Mesh Levels` splits
/// (`LoadModel`, ModelDescription.F90:3049), then `Extruded Mesh *`
/// extrusion (`CreateExtrudedMesh`, ElmerSolver.F90:897). Extrusion works on
/// the `MapBodiesAndBCs` tags and leaves body and BC indices as the new
/// mesh's tags, so the model's `Target Bodies` / `Target Boundaries` become
/// the identity. `None` when the mesh is used as loaded.
fn transform_mesh(model: &mut Model, mesh: &dyn MeshView) -> Result<Option<LevelsMesh>, SolverError> {
    let err = |e: MeshLevelsError| match e {
        MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshLevelsError::Invalid(m) => SolverError::fatal("LoadModel", m),
    };
    let ops_err = |e: MeshOpsError| match e {
        MeshOpsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshOpsError::Invalid(m) => SolverError::fatal("MeshOps", m),
    };
    let identity_targets = |model: &mut Model| {
        for (i, body) in model.bodies.iter_mut().enumerate() {
            body.insert("Target Bodies", Value::Integer(vec![i as i32 + 1]));
        }
        for (i, bc) in model.boundary_conditions.iter_mut().enumerate() {
            bc.insert("Target Boundaries", Value::Integer(vec![i as i32 + 1]));
        }
    };
    // `CreateDiscontMesh` (LoadMesh2) comes first, on the mapped mesh.
    let mut base: Option<LevelsMesh> = None;
    if model.boundary_conditions.iter().any(|bc| bc.logical("discontinuous boundary") == Some(true)) {
        let d = create_discont_mesh(&MappedTags::new(mesh, model), &model.simulation, &model.boundary_conditions)
            .map_err(ops_err)?;
        if let Some(d) = d {
            base = Some(d.mesh);
            identity_targets(model);
        }
    }
    let discont: &dyn MeshView = match &base {
        Some(m) => m,
        None => mesh,
    };
    let split = mesh_levels(discont, &model.simulation, &model.boundary_conditions).map_err(err)?;
    let split = split.or(base);
    let current: &dyn MeshView = match &split {
        Some(m) => m,
        None => mesh,
    };
    let mut density = |t: f64| mesh_density(&model.simulation, t);
    let Some(opts) = ExtrudeOptions::from_params(&model.simulation, current.dim(), Some(&mut density)).map_err(err)? else {
        if model.simulation.logical("increase element order") == Some(true) {
            return Ok(Some(increase_element_order(current).map_err(ops_err)?));
        }
        return Ok(split);
    };
    let extruded = extrude(&MappedTags::new(current, model), &opts).map_err(err)?;
    identity_targets(model);
    Ok(Some(extruded))
}

/// `GeneratePeriodicProjectors` (ElmerSolver.F90) when a boundary
/// condition has `Periodic BC`: `elmer-mortar`'s projectors, applied by its
/// [`PeriodicAssembly`] in place of the default assembly. `Save Projector`
/// files go to the case directory.
fn periodic_assembly(
    model: &Model,
    mesh: &dyn MeshView,
    coord_dim: usize,
    case_dir: &Path,
    log: &dyn Logger,
) -> Result<Option<PeriodicAssembly>, SolverError> {
    if !model.boundary_conditions.iter().any(|bc| bc.get("periodic bc").is_some()) {
        return Ok(None);
    }
    let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
    let bc_of = |tag: i32| model.bc_of_tag(tag);
    let env = ProjectorEnv {
        mesh,
        bcs: &bcs,
        bc_of: &bc_of,
        simulation: &model.simulation,
        solver: model.solvers.first().map(|s| s as &dyn Params),
        coord_dim,
        log,
    };
    let (projectors, stop) = generate_periodic_projectors(&env, case_dir, 0)?;
    if stop {
        return Err(SolverError::Unsupported("Save Projector And Stop".into()));
    }
    Ok(Some(PeriodicAssembly::new(model.boundary_conditions.clone(), projectors, mesh.n_nodes())))
}

/// The `fem/tests` case name `elmer-fixtures` keys its ports by: the case
/// directory's name (the `.sif`'s directory; empty means the working one).
fn fixture_case(case_dir: &Path) -> String {
    let dir = if case_dir.as_os_str().is_empty() { std::env::current_dir().ok() } else { case_dir.canonicalize().ok() };
    dir.and_then(|d| d.file_name().and_then(|n| n.to_str()).map(str::to_string)).unwrap_or_default()
}

/// Runs `case` serially. Unsupported input comes back as
/// [`SolverError::Unsupported`] or [`SolverError::NotWired`] before any
/// file is written.
pub fn run_case(case: Case<'_>, log: &dyn Logger, comm: &dyn Communicator) -> Result<Outcome, SolverError> {
    if comm.size() > 1 {
        return Err(SolverError::Unsupported("more than one MPI rank".into()));
    }
    run_case_with(case, log, &RunEnv { job: comm, globals: Vec::new() })
}

/// The parallel job a run belongs to when every rank solves its own copy of
/// the whole mesh (`Single Mesh = True` under `ElmerSolver_mpi`): the
/// `Simulation` itself stays serial, as Elmer's `Mesh % SingleMesh` makes
/// it. From `job`'s rank and size, rank 0 writes `TEST.PASSED_<np>`, VTU
/// output is named per partition with a `.pvtu` from rank 0 (`VtuFileNaming`
/// with `Part`), and the `Partition` and `Parallel Slices` globals are set.
/// `globals` are further global variables, set after those.
pub struct RunEnv<'a> {
    pub job: &'a dyn Communicator,
    pub globals: Vec<(String, f64)>,
}

/// `ExecSimulation`'s run-time globals: `Partition` (ElmerSolver.F90:1774;
/// only in a parallel job here, so serial cases keep their variable set)
/// and, when `Parallel Slices` is given, one slice per rank: `slice`,
/// `slice ratio` and `slice weight` (:2785-2875), with `Number Of Slices`
/// and `Number Of Times` added to the Simulation section.
fn job_globals(simulation: &mut ParamMap, rank: usize, size: usize) -> Result<Vec<(String, f64)>, SolverError> {
    let mut out = Vec::new();
    if size > 1 {
        out.push(("partition".to_string(), rank as f64));
    }
    if simulation.get("parallel timestepping").is_some() {
        return Err(SolverError::Unsupported("Simulation: Parallel Timestepping".into()));
    }
    if simulation.get("parallel slices").is_none() {
        return Ok(out);
    }
    let slices = simulation.logical("parallel slices") == Some(true);
    if slices && size > 1 && simulation.logical("single mesh") != Some(true) {
        return Err(SolverError::fatal("ExecSimulation", "Parallel time and slices only available with \"Single Mesh\""));
    }
    let n = if slices { size } else { 1 };
    if slices {
        simulation.insert("Number Of Times", Value::Integer(vec![1]));
        simulation.insert("Number Of Slices", Value::Integer(vec![n as i32]));
    }
    if n > 1 && simulation.get("extruded max coordinate").is_some() {
        return Err(SolverError::Unsupported("Parallel Slices moved by Extruded Min/Max Coordinate".into()));
    }
    let (slice, ratio, weight) = if n <= 1 {
        (0.0, 0.0, 1.0)
    } else {
        (rank as f64, (rank as f64 + 0.5) / n as f64 - 0.5, 1.0 / n as f64)
    };
    out.extend([("slice".to_string(), slice), ("slice ratio".to_string(), ratio), ("slice weight".to_string(), weight)]);
    Ok(out)
}

/// [`run_case`] as one rank of `env.job`; see [`RunEnv`].
pub fn run_case_with(case: Case<'_>, log: &dyn Logger, env: &RunEnv<'_>) -> Result<Outcome, SolverError> {
    let (rank, size) = (env.job.rank(), env.job.size());
    let comm: &dyn Communicator = if size == 1 { env.job } else { &SerialComm };
    let Case { mut model, reloads, mesh, case_dir, mesh_dir } = case;
    let globals = job_globals(&mut model.simulation, rank, size)?;
    let fixture_case = fixture_case(&case_dir);
    let registry = solver_registry_for(&fixture_case);
    read_target_names(&mesh_dir.join("mesh.names"), &model.simulation, &mut model.bodies, &mut model.boundary_conditions)
        .map_err(|e| SolverError::fatal("ReadTargetNames", e.to_string()))?;
    let transformed = transform_mesh(&mut model, mesh)?;
    let mesh: &dyn MeshView = match &transformed {
        Some(m) => m,
        None => mesh,
    };
    let mesh = MovableMesh::new(mesh);
    let _motion = mesh.install();
    let mesh: &dyn MeshView = &mesh;
    let mut outputs = prepare_vtu_output(&mut model, &case_dir, &mesh_dir, size > 1)?;
    let coords = coordinate_system(&model, mesh)?;
    legacy_solvers(&mut model, coords);
    drop_never_savedata(&mut model);
    preflight(&model, &registry)?;
    check_mesh_and_elements(&model)?;
    check_reloads(&model, &reloads, &registry)?;
    install_fe(coords.dim)?;

    let periodic = periodic_assembly(&model, mesh, coords.dim, &case_dir, log)?;

    // The test-local plugins' `<Proc>_init` section edits (before Elmer
    // creates their variables).
    elmer_fixtures::registry::init_sections(&fixture_case, &mut model);
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_datatofield::is_section(p)) {
        elmer_physics_datatofield::init_section(p);
    }
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(fe_backend(coords.dim));
    // Solver sections with a non-nodal `Element`; `p:%f` degrees are
    // `GetMatcReal(f, 4, [body, x, y, z], 'cx')` (MeshBasics.F90).
    let degree: Arc<elmer_pdofs::solver::MatcFn> = Arc::new(|f: &str, cx: [f64; 4]| {
        SifMatc.eval_with("cx", f, &cx)?.first().copied().ok_or_else(|| format!("MATC {f}: no value"))
    });
    sim.backends.dofs = Some(Box::new(PDofsBackend::new().with_matc(degree)));
    if let Some(periodic) = periodic {
        sim.backends.assembly = Some(Box::new(periodic));
    }
    sim.set_matc(Rc::new(SifMatc));
    for (name, value) in globals.iter().chain(&env.globals) {
        sim.vars.set_scalar(name, *value);
    }
    elmer_fixtures::registry::install_live(&fixture_case, &mut sim, mesh);
    sim.setup(mesh, &registry, log, comm)?;
    check_mesh_velocity(&sim)?;
    // Field lists are fixed once every variable exists, before the run, so a
    // refusal still comes before any file.
    for out in &mut outputs {
        out.select_fields(&sim, mesh.dim())?;
    }
    let written = Arc::new(Mutex::new(Vec::new()));
    for out in outputs {
        let id = out.index + 1;
        let slot = sim
            .slots
            .iter_mut()
            .find(|s| s.id == id)
            .ok_or_else(|| SolverError::fatal("ResultOutputSolver", format!("no slot for Solver {id}")))?;
        slot.set_instance(Box::new(out.into_solver(Arc::clone(&written)).in_job(rank, size)));
    }
    // `Apply Mortar BCs`: the slot's solver behind the mortar constraints.
    for slot in &mut sim.slots {
        let params = &sim.model.solvers[slot.id - 1];
        if !crate::mortar::applies(params) {
            continue;
        }
        let factory = procedure_of(params)
            .and_then(|(lib, entry)| registry.get(&lib, &entry))
            .ok_or_else(|| SolverError::Unsupported(format!("Solver {}: Apply Mortar BCs without a registered Procedure", slot.id)))?;
        let inner = factory(params)?;
        slot.set_instance(Box::new(crate::mortar::MortarSolver {
            inner,
            coord_dim: coords.dim,
            case_dir: case_dir.clone(),
            ready: false,
        }));
    }
    let mut report = sim.run(mesh, log, comm)?;
    for block in &reloads {
        sim.reload(block, mesh, &registry, log, comm)?;
        rebuild_changed(&mut sim, block, &registry)?;
        report = sim.run(mesh, log, comm)?;
    }
    let written = std::mem::take(&mut *written.lock().unwrap());
    let passed = finalize_reference(&sim.reference, &case_dir, log, env.job)?;
    let norms = sim
        .slots
        .iter()
        .filter_map(|s| Some((s.id, sim.vars.by_id(s.variable?).conv.norm)))
        .collect();
    Ok(Outcome { report, norms, passed, written })
}

#[cfg(test)]
mod job_globals_tests {
    use super::*;

    fn sim(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    #[test]
    fn serial_case_without_slices_gets_no_globals() {
        let mut p = sim(&[]);
        assert!(job_globals(&mut p, 0, 1).unwrap().is_empty());
    }

    #[test]
    fn one_slice_per_rank() {
        let mut p = sim(&[("Parallel Slices", Value::Logical(true)), ("Single Mesh", Value::Logical(true))]);
        let g = job_globals(&mut p, 1, 4).unwrap();
        let get = |n: &str| g.iter().find(|(k, _)| k == n).map(|(_, v)| *v);
        assert_eq!(get("partition"), Some(1.0));
        assert_eq!(get("slice"), Some(1.0));
        assert_eq!(get("slice ratio"), Some(1.5 / 4.0 - 0.5));
        assert_eq!(get("slice weight"), Some(0.25));
        assert_eq!(p.integer("number of slices"), Some(4));
    }

    #[test]
    fn serial_slices_are_slice_zero() {
        let mut p = sim(&[("Parallel Slices", Value::Logical(true))]);
        let g = job_globals(&mut p, 0, 1).unwrap();
        assert_eq!(g, vec![("slice".into(), 0.0), ("slice ratio".into(), 0.0), ("slice weight".into(), 1.0)]);
    }

    #[test]
    fn parallel_slices_need_single_mesh() {
        let mut p = sim(&[("Parallel Slices", Value::Logical(true))]);
        assert!(matches!(job_globals(&mut p, 0, 2), Err(SolverError::Fatal { .. })));
    }
}
