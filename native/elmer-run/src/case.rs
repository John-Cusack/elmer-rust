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
use elmer_mortar::{generate_conforming_map, generate_periodic_projectors, PeriodicAssembly, ProjectorEnv};
use elmer_mesh_levels::{extrude, mesh_levels, split_mesh_equal, ExtrudeOptions, Mesh as LevelsMesh, MeshLevelsError};
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
        "mesh update" => Some(("MeshSolve", "MeshSolver")),
        _ => None,
    }
}

/// Ids of the solvers other solvers call as slaves (`Pre Solvers`,
/// `Post Solvers`, their `Nonlinear` forms, `Slave Solvers`): they run
/// whatever their own `Exec Solver` says.
fn slave_solvers(model: &Model) -> Vec<usize> {
    const KEYS: [&str; 6] =
        ["pre solvers", "post solvers", "nonlinear pre solvers", "nonlinear post solvers", "slave solvers", "assembly solvers"];
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
/// `ReynoldsSolver_init` and `ReynoldsPostprocess_Init`
/// (ReynoldsSolver.F90:983, :1639): the film pressure variable, and the
/// postprocessor's `<p> Heating` variable, `-dofs d <p> Force / Flux / Mean
/// Velocity` exports (`d` = 2 in 1D or for a `Dim` of 2, else 3) and its
/// linear system defaults.
fn reynolds_init(params: &mut ParamMap, cdim: usize) {
    let Some((lib, entry)) = procedure_of(params) else { return };
    if !lib.eq_ignore_ascii_case("reynoldssolver") {
        return;
    }
    if entry.eq_ignore_ascii_case("reynoldssolver") {
        if params.get("variable").is_none() {
            params.insert("Variable", Value::Str("FilmPressure".into()));
        }
        return;
    }
    if !entry.eq_ignore_ascii_case("reynoldspostprocess") {
        return;
    }
    let pname = params.string("reynolds pressure variable name").unwrap_or("FilmPressure").to_string();
    if params.logical("calculate heating") == Some(true) {
        params.insert("Variable", Value::Str(format!("{pname} Heating")));
    } else if params.get("variable").is_none() {
        params.insert("Variable", Value::Str("-nooutput ReynoldsPost".into()));
    }
    for (flag, dim_key, field) in [
        ("calculate force", "calculate force dim", "Force"),
        ("calculate flux", "calculate flux dim", "Flux"),
        ("calculate mean velocity", "calculate mean velocity dim", "Mean Velocity"),
    ] {
        if params.logical(flag) != Some(true) {
            continue;
        }
        let dofs = if cdim == 1 || params.integer(dim_key) == Some(2) { 2 } else { 3 };
        let k = (1..).find(|k| params.get(&format!("exported variable {k}")).is_none()).expect("a free index");
        params.insert(&format!("Exported Variable {k}"), Value::Str(format!("-dofs {dofs} {pname} {field}")));
    }
    params.insert("Time derivative order", Value::Integer(vec![0]));
    for (key, value) in [
        ("linear system solver", Value::Str("Iterative".into())),
        ("linear system iterative method", Value::Str("cg".into())),
        ("linear system preconditioning", Value::Str("ILU0".into())),
        ("linear system max iterations", Value::Integer(vec![500])),
        ("linear system residual output", Value::Integer(vec![10])),
        ("linear system convergence tolerance", Value::real(1.0e-10)),
    ] {
        if params.get(key).is_none() {
            params.insert(key, value);
        }
    }
}

/// The solver sections' `_init` defaults (variables, exported fields and
/// keywords a procedure's `_init` adds), applied before setup; the
/// partitioned runner applies them too.
pub fn legacy_solvers(model: &mut Model, coords: CoordinateSystem) {
    let cdim = coords.dim;
    let transient = model.simulation.string("simulation type").is_some_and(|t| t.eq_ignore_ascii_case("transient"));
    let serendipity = model.simulation.logical("serendipity p elements").unwrap_or(true);
    let equations = model.equations.clone();
    let materials = model.materials.clone();
    for params in &mut model.solvers {
        reynolds_init(params, cdim);
        if crate::saveprojection::is_section(params) {
            crate::saveprojection::init_section(params);
        }
        if crate::gmsh_reader::is_section(params) {
            crate::gmsh_reader::init_section(params);
        }
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
        // ElasticSolver_Init0 / ElasticSolver_Init (ElasticSolve.F90:56-60,
        // 101-115): the neo-Hookean mixed formulation defaults to `p:2`
        // elements and the variable `MixedSol[Disp:dim Pres:1]`.
        let elastic_only = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("elasticsolve") && entry.eq_ignore_ascii_case("elasticsolver")
        });
        if elastic_only && mixed {
            if params.get("element").is_none() {
                params.insert("Element", Value::Str("p:2".into()));
            }
            if params.get("variable").is_none() {
                params.insert("Variable", Value::Str(format!("MixedSol[Disp:{cdim} Pres:1]")));
                params.insert("Variable Dofs", Value::Integer(vec![cdim as i32 + 1]));
            }
        }
        // ElasticSolve.F90:180: `ListAddNewLogical('Bubbles in Global
        // System', .TRUE.)`, which also sizes a nodal instance's permutation
        // for the mesh's bubbles (AddEquationBasics).
        let elasticsolve = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("elasticsolve") && entry.eq_ignore_ascii_case("elasticsolver")
        });
        if elasticsolve && params.get("bubbles in global system").is_none() {
            params.insert("Bubbles in Global System", Value::Logical(true));
        }
        // ElasticSolver_Init: a UMAT material's integration-point fields.
        if elasticsolve && materials.iter().any(|m| m.get("umat subroutine").is_some()) {
            let out = params.logical("output state variables") == Some(true);
            let no = if out { "" } else { "-nooutput " };
            let stress = if cdim == 3 {
                "UmatStress[UmatStress_xx:1 UmatStress_yy:1 UmatStress_zz:1 UmatStress_xy:1 UmatStress_yz:1 UmatStress_xz:1]"
            } else {
                "UmatStress[UmatStress_xx:1 UmatStress_zz:1 UmatStress_yy:1 UmatStress_xy:1]"
            };
            let mut defs = vec![format!("{no}-ip {stress}"), format!("{no}-dofs 3 -ip UmatEnergy")];
            let nstate = materials
                .iter()
                .filter(|m| m.get("umat subroutine").is_some())
                .map(|m| m.integer("number of state variables").unwrap_or(0))
                .max()
                .unwrap_or(0);
            if nstate > 0 {
                defs.push(format!("{no}-dofs {nstate} -ip UmatState"));
            }
            for def in defs {
                let k = (1..).find(|k| params.get(&format!("exported variable {k}")).is_none()).unwrap_or(1);
                params.insert(&format!("Exported Variable {k}"), Value::Str(def));
            }
        }
        // ElasticSolver_Init / StressSolver_Init (ElasticSolve.F90:145-168,
        // StressSolve.F90:294-336): a Maxwell material, in any material or
        // the solver section, integrates first order (BDF1) and keeps its lag
        // stress in `ve_stress` at the integration points; at least two
        // nonlinear iterates.
        let maxwell = materials.iter().any(|m| m.logical("maxwell material") == Some(true))
            || params.logical("maxwell material") == Some(true);
        if elastic && maxwell {
            let axi = matches!(
                coords.kind,
                CoordinateKind::AxisSymmetric | CoordinateKind::CylindricSymmetric
            );
            let nve = match (cdim, axi) {
                (2, true) => 4,
                (2, false) => 3,
                (3, _) => 6,
                _ => 1,
            };
            params.insert("Timestepping Method", Value::Str("BDF".into()));
            params.insert("BDF Order", Value::Integer(vec![1]));
            params.insert("Time derivative Order", Value::Integer(vec![1]));
            let k = (1..).find(|k| params.get(&format!("exported variable {k}")).is_none()).unwrap_or(1);
            params.insert(&format!("Exported Variable {k}"), Value::Str(format!("-dofs {nve} -ip ve_stress")));
            for key in ["Nonlinear System Min Iterations", "Nonlinear System Max Iterations"] {
                let i = params.integer(key).unwrap_or(0);
                params.insert(key, Value::Integer(vec![i.max(2)]));
            }
        }
        // `HeatSolver_init` (HeatSolve.F90:167): `Variable = Temperature`.
        let heat = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case("heatsolve") && entry.eq_ignore_ascii_case("heatsolver")
        });
        if heat && params.get("variable").is_none() {
            params.insert("Variable", Value::Str("Temperature".into()));
        }
        // `HeatSolver_Init0` / `_Init`: the residual-free bubble element of
        // unstabilized convection.
        if heat {
            elmer_physics_heat::solver_defaults(params, &equations, serendipity, transient);
        }
        // `IncompressibleNSSolver_Init0` / `_Init` (IncompressibleNS.F90):
        // the `Variable`, `Element` and discretisation defaults.
        let ns = procedure_of(params).is_some_and(|(lib, entry)| {
            (lib.eq_ignore_ascii_case("incompressiblens") || lib.eq_ignore_ascii_case("incompressiblensvec"))
                && entry.eq_ignore_ascii_case("incompressiblenssolver")
        });
        if ns {
            elmer_physics_incompressible::solver_defaults_with(params, cdim, transient, serendipity);
        }
        // `FilmFlowSolver_Init0` / `_Init`: the MINI element and the `Flow`
        // variable of the `Model Dimension`.
        let film = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case(elmer_physics_filmflow::PROCEDURE.0)
                && entry.eq_ignore_ascii_case(elmer_physics_filmflow::PROCEDURE.1)
        });
        if film {
            elmer_physics_filmflow::solver_defaults(params, serendipity);
        }
        // `CoilSolver_init` and `DivergenceSolver_Init`: their variables.
        let is = |p: &ParamMap, (l, e): (&str, &str)| {
            procedure_of(p).is_some_and(|(lib, entry)| lib.eq_ignore_ascii_case(l) && entry.eq_ignore_ascii_case(e))
        };
        if is(params, elmer_physics_coil::PROCEDURE) {
            elmer_physics_coil::solver_defaults(params, cdim);
        }
        if is(params, elmer_physics_flux::DIVERGENCE_PROCEDURE) {
            elmer_physics_flux::divergence_defaults(params);
        }
        // `MarchingODESolver_init`: no matrix.
        if is(params, ("MarchingODESolver", "MarchingODESolver")) {
            elmer_physics_meshtools::marching_defaults(params);
        }
        // `ShellSolver_Init0`: the deflection variable, large deflection
        // and the Simulation's `Initialize Dirichlet Conditions`.
        if is(params, elmer_physics_shell::PROCEDURE) {
            elmer_physics_shell::solver_defaults(params, &mut model.simulation);
        }
        // `ParticleDynamics_Init`: the global pseudo variable.
        if is(params, elmer_physics_particle::PROCEDURE) {
            elmer_physics_particle::solver_defaults(params);
        }
        // `MixedPoisson_Init0`: the face element, the global bubbles and the
        // exported potential and flux.
        if is(params, elmer_pdofs::mixed_poisson::PROCEDURE) {
            elmer_pdofs::mixed_poisson::solver_defaults(params, model.simulation.string("coordinate system"));
        }
        // `SmitcSolver_Init0`: `-dofs 3 Deflection` and the global mass matrix.
        if is(params, elmer_physics_shell::SMITC_PROCEDURE) {
            elmer_physics_shell::smitc_defaults(params);
        }
        // FacetShellSolve's `ShellSolver_Init0`: six-dof `Deflection`.
        if is(params, elmer_physics_shell::FACET_PROCEDURE) {
            elmer_physics_shell::facet_defaults(params);
        }
        // `MeshSolver_Init`: a transient run's `-dofs dim Mesh Velocity`,
        // unless `Compute Mesh Velocity = False` (or the legacy equation
        // already exported it).
        let mesh_solve = procedure_of(params)
            .is_some_and(|(l, e)| l.eq_ignore_ascii_case("meshsolve") && e.eq_ignore_ascii_case("meshsolver"))
            || params.string("equation").is_some_and(|q| q.trim().eq_ignore_ascii_case("mesh update"));
        if mesh_solve && transient && params.logical("compute mesh velocity") != Some(false) {
            let exported = params.0.iter().any(|(k, v)| {
                k.starts_with("exported variable")
                    && matches!(v, Value::Str(t) if normalize_key(t).ends_with("mesh velocity"))
            });
            if !exported {
                let k = (1..).find(|k| params.get(&format!("exported variable {k}")).is_none()).expect("a free index");
                params.insert(&format!("Exported Variable {k}"), Value::Str(format!("-dofs {cdim} Mesh Velocity")));
            }
        }
        // `TransientPhaseChange_Init`: the interface velocity field.
        if is(params, elmer_physics_meshtools::TRANSIENT_PHASE_CHANGE) {
            elmer_physics_meshtools::transient_defaults(params);
        }
        // `OdeSolver_init`: the global ODE matrix and its time terms.
        if is(params, elmer_physics_scalar::ode_solver::PROCEDURE) {
            elmer_physics_scalar::ode_solver::solver_defaults(params);
        }
        // `NormalSolver_Init`: its variable and the normals' result field.
        if is(params, elmer_physics_scalar::normal_solver::PROCEDURE) {
            elmer_physics_scalar::normal_solver::solver_defaults(params, cdim);
        }
        // `TopoOpt_init0` / `_init`: the elemental topology fields.
        if is(params, elmer_physics_topoopt::PROCEDURE) {
            elmer_physics_topoopt::solver_defaults(params);
        }
        // AddEquationBasics (MainUtils.F90:1231): a legacy `Equation` name
        // gives the `Variable` when there is none, with or without a
        // `Procedure`, and the `Procedure` when there is none.
        let has_procedure = params.get("procedure").is_some();
        let name = params.string("equation").map(|s| s.trim().to_ascii_lowercase());
        let (var, dofs, procedure) = match name.as_deref() {
            Some("heat equation") => {
                if !has_procedure && params.get("radiation solver").is_none() {
                    params.insert("Radiation Solver", Value::Logical(true));
                }
                ("Temperature", 1, "HeatSolve HeatSolver")
            }
            Some("stress analysis") => ("Displacement", cdim, "StressSolve StressSolver"),
            // A transient run's `Mesh Velocity` is exported above, as
            // AddEquationBasics and MeshSolver_Init do (once).
            Some("mesh update") => ("Mesh Update", cdim, "MeshSolve MeshSolver"),
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
        if !has_procedure {
            params.insert("Procedure", Value::Str(procedure.into()));
        }
    }
}

/// The `_init` expansions of the lumped-matrix flags (`Calculate
/// Conductance / Conductivity / Capacitance Matrix` of HeatSolve,
/// StatCurrentSolve, StatElecSolve, and their Vec versions) into
/// constraint-mode lumping, before the variables and permutations exist
/// ([`elmer_solver::lumping_init`]).
pub fn lumped_matrix_inits(model: &mut Model) {
    for i in 0..model.solvers.len() {
        let Some((lib, entry)) = procedure_of(&model.solvers[i]) else { continue };
        let (lib, entry) = (lib.to_ascii_lowercase(), entry.to_ascii_lowercase());
        let l = match (lib.trim_end_matches("vec"), entry.as_str()) {
            ("heatsolve", "heatsolver") => &elmer_solver::CONDUCTANCE,
            ("statcurrentsolve", "statcurrentsolver") => &elmer_solver::CONDUCTIVITY,
            ("statelecsolve", "statelecsolver") => &elmer_solver::CAPACITANCE,
            _ => continue,
        };
        elmer_solver::lumping_init(&mut model.solvers[i], &mut model.boundary_conditions, &mut model.body_forces, l);
    }
}

/// Whether a keyword of `p` is a `Procedure` of library `lib` (any case,
/// with or without `.so`).
fn names_library(p: &ParamMap, lib: &str) -> bool {
    p.0.values().any(|v| match v {
        Value::Dependent(d) => match &d.kind {
            DependentKind::Procedure { library, .. } => elmer_expr::udf::normalize_library(library).eq_ignore_ascii_case(lib),
            _ => false,
        },
        _ => false,
    })
}

/// `StatElecBoundary.F90`'s functions, when a Boundary Condition or a
/// Material names them: they read the section of `Model % CurrentElement`,
/// so they get every element's section, the live fields and MATC. Both
/// the simulation's registry and the process-wide one get them (some
/// solvers resolve keyword procedures in the latter only).
fn install_statelec_boundary(sim: &mut Simulation, mesh: &dyn MeshView) {
    use elmer_physics_scalar::statelec_boundary::{register, ListOf, Lists, LIBRARY};
    let m = &sim.model;
    if !m.boundary_conditions.iter().chain(&m.materials).any(|p| names_library(p, LIBRARY)) {
        return;
    }
    let nb = mesh.n_bulk();
    let list_of = (0..nb + mesh.n_boundary())
        .map(|t| {
            let tag = mesh.element(elmer_core::ElemId(t as u32)).tag;
            if t < nb {
                let body = m.body_of_tag(tag)?;
                let mat = m.bodies.get(body - 1)?.integer("material")?;
                (mat >= 1).then(|| ListOf::Material(mat as usize - 1))
            } else {
                m.bc_of_tag(tag).map(|b| ListOf::Bc(b - 1))
            }
        })
        .collect();
    let lists = Lists {
        bcs: m.boundary_conditions.clone(),
        materials: m.materials.clone(),
        list_of,
        eps0: m.constants.real("permittivity of vacuum").unwrap_or(elmer_physics_scalar::PHYSICAL_EPS0),
    };
    let fields: Arc<dyn elmer_expr::FieldLookup + Send + Sync> = sim.live_model(mesh);
    let udfs = Arc::clone(&sim.backends.udfs);
    for into in [&*sim.backends.udfs, elmer_expr::UdfRegistry::global()] {
        let lists = Lists { bcs: lists.bcs.clone(), materials: lists.materials.clone(), list_of: lists.list_of.clone(), eps0: lists.eps0 };
        register(into, lists, Arc::clone(&fields), Arc::new(SifMatc), Arc::clone(&udfs));
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
    crate::unported::check(model)
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
        // `Angular Frequency` / `Frequency` are read by the solvers, which
        // are rebuilt after the reload (HelmholtzFEM's second block); other
        // Simulation and Constants changes stay refused.
        const SIM_ALLOWED: [&str; 3] = ["solver input file", "angular frequency", "frequency"];
        if let Some(key) = block.simulation.0.keys().find(|k| !SIM_ALLOWED.contains(&k.as_str())) {
            return refuse(format!("Simulation {key}"));
        }
        if let Some(key) = block.constants.0.keys().next() {
            return refuse(format!("Constants {key}"));
        }
        // The periodic projectors are built once, at setup
        // (`periodic_assembly`): a change of the periodic pairing itself
        // (`Periodic BC`, `Periodic BC <var>`, ...) stays refused.
        // `Periodic BC Offset <var>` is handed to the assembly after the
        // reload (`reload_boundary_conditions`).
        for (i, bc) in block.boundary_conditions.iter().enumerate() {
            if let Some(key) = bc.0.keys().find(|k| k.starts_with("periodic bc") && !k.starts_with("periodic bc offset ")) {
                return refuse(format!("Boundary Condition {} {key}", i + 1));
            }
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

/// Inputs no provider applies yet: mesh-load transforms
/// beyond `Coordinate Scaling`, `Mesh Levels` and extrusion, and element
/// definitions other than the Lagrange (`n:1`) basis `elmer-fe` evaluates.
fn check_mesh_and_elements(model: &Model) -> Result<(), SolverError> {
    let sim = &model.simulation;
    // `Coordinate Mapping` and `Coordinate Transformation` are applied at
    // load (`mesh_db::load`), before any extrusion; Elmer transforms after it.
    let extruded = sim.get("extruded mesh levels").is_some() || sim.get("extruded mesh density").is_some();
    if sim.get("coordinate transformation").is_some() && extruded {
        return Err(SolverError::Unsupported("Simulation: Coordinate Transformation of an extruded mesh".into()));
    }
    match sim.get("coordinate mapping") {
        None => {}
        Some(Value::Integer(v)) if v[..] == [1, 2, 3] => {}
        Some(Value::Integer(_)) if !extruded => {}
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
/// the Lagrange `n:1`; `p:1` with or without `b:N` bubbles for
/// ThermoElectricSolver, which takes them from `elmer-pdofs` and keeps the
/// bubbles in the global system and its norm as Elmer does (#228); and
/// any definition for the solvers that read their DOFs from `elmer-pdofs`
/// (`ctx.solver_dofs()`) and refuse what they do not port: StressSolve
/// (#206), StatElecSolve (#208), IncompressibleNS (Taylor-Hood), HeatSolve
/// (its residual-free bubble; nothing beyond `p:1`), ModelPDE's
/// AdvDiffSolver and the test-local `Poisson` / `PoissonML`
/// `PoissonSolver`, `ModelPDEsplitted`, `SingleElements`, `BenchmarkShell`,
/// `Stokes1` and `AdvReactDG_P` (DG on p elements) ports.
pub fn native_element(params: &ParamMap) -> bool {
    let Some(e) = params.string("element") else { return true };
    if normalize_key(e) == "n:1" {
        return true;
    }
    // A Block Solver assembled by `Bulk Assembly Procedure ij` routines:
    // BlockSystemAssembly lays the blocks out from elmer-pdofs, and a
    // routine that is not ported stops at "Assembly routine not found".
    if params.logical("block solver") == Some(true) && params.0.keys().any(|k| k.starts_with("bulk assembly procedure")) {
        return true;
    }
    let Some((lib, entry)) = procedure_of(params) else { return false };
    let lib = lib.trim_start_matches("./");
    const P_DOFS: [(&str, &str); 27] = [
        ("advreactdg_p", "advectionreactionsolver"),
        ("benchmarkshell", "shellsolver"),
        ("elasticsolve", "elasticsolver"),
        ("filmflowsolver", "filmflowsolver"),
        ("flowsolve", "flowsolver"),
        ("heatsolve", "heatsolver"),
        ("incompressiblens", "incompressiblenssolver"),
        ("incompressiblensvec", "incompressiblenssolver"),
        ("statelecsolve", "statelecsolver"),
        ("statelecsolvevec", "statelecsolver"),
        ("stresssolve", "stresssolver"),
        ("stokes", "stokessolver"),
        ("stokes1", "stokessolver"),
        ("q1q0", "stokessolver"),
        ("poisson", "poissonsolver"),
        ("poissonml", "poissonsolver"),
        ("modelpde", "advdiffsolver"),
        ("modelpdesplitted", "advdiffsolver"),
        ("singleelements", "singleelements"),
        ("modelmixedpoisson", "mixedpoisson"),
        ("facefeinterpolationtest", "bestapproximationsolver"),
        ("edgeelementtest", "edgeelementsolver"),
        ("edgefeinterpolationtest", "bestapproximationsolver"),
        ("magnetodynamics", "whitneyavsolver"),
        ("hrot", "hrotsolver"),
        ("helmholtz", "helmholtzhdivsolver"),
        // RefineMesh's `RTFlux` carrier (`Flux Recovery`).
        ("vectorhelmholtz", "vectorhelmholtz_dummy"),
    ];
    if P_DOFS.iter().any(|(l, n)| lib.eq_ignore_ascii_case(l) && entry.eq_ignore_ascii_case(n)) {
        return true;
    }
    lib.eq_ignore_ascii_case("thermoelectricsolver") && elmer_physics_thermoelectric::parse_element(e).is_ok()
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
        DependentKind::Table { x, y, cubic, monotone } => {
            let spline = cubic.then(|| elmer_expr::curve::cubic_spline(x, y, *monotone));
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
/// `discont_perm` is set to `Mesh % DiscontPerm` when `CreateDiscontMesh`
/// split the mesh and no later transform renumbered its nodes.
fn transform_mesh(model: &mut Model, mesh: &dyn MeshView, discont_perm: &mut Option<Vec<u32>>) -> Result<Option<LevelsMesh>, SolverError> {
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
            if d.discont_mesh {
                *discont_perm = Some(d.perm);
            }
            base = Some(d.mesh);
            identity_targets(model);
        }
    }
    let discont: &dyn MeshView = match &base {
        Some(m) => m,
        None => mesh,
    };
    // `Follow * Boundary`: the loaded mesh and each split follow the
    // curved boundaries (FollowCurvedBoundary).
    let split = if crate::follow::any_follow(&model.boundary_conditions) {
        if model.simulation.integer("mesh levels").unwrap_or(1) > 1 {
            *discont_perm = None;
        }
        Some(crate::follow::follow_levels(discont, model)?)
    } else {
        let split = mesh_levels(discont, &model.simulation, &model.boundary_conditions).map_err(err)?;
        if split.is_some() {
            *discont_perm = None;
        }
        split
    };
    let split = split.or(base);
    let current: &dyn MeshView = match &split {
        Some(m) => m,
        None => mesh,
    };
    let mut density = |t: f64| mesh_density(&model.simulation, t);
    let Some(opts) = ExtrudeOptions::from_params(&model.simulation, current.dim(), Some(&mut density)).map_err(err)? else {
        if model.simulation.logical("increase element order") == Some(true) {
            *discont_perm = None;
            return Ok(Some(increase_element_order(current).map_err(ops_err)?));
        }
        return Ok(split);
    };
    *discont_perm = None;
    let extruded = extrude(&MappedTags::new(current, model), &opts).map_err(err)?;
    identity_targets(model);
    Ok(Some(extruded))
}

/// A geometric multigrid solver's mesh hierarchy: with `Linear System
/// Solver = Multigrid` and `MG Equal Split`, `AddEquationBasics`
/// (MainUtils.F90) splits the solver's mesh `MG Levels - 1` times and moves
/// the solver to the finest. The solver id, its coarse meshes (coarsest
/// first) and the finest mesh; `None` when no solver asks for it, so
/// `elmer-solver` refuses the other multigrid variants itself.
type MultigridMeshes = (usize, Vec<Box<dyn MeshView>>, LevelsMesh);

fn multigrid_meshes(model: &Model, mesh: &dyn MeshView, ranks: usize) -> Result<Option<MultigridMeshes>, SolverError> {
    let is_mg = |p: &ParamMap| {
        p.string("linear system solver").is_some_and(|s| normalize_key(s) == "multigrid")
            && p.logical("mg equal split") == Some(true)
    };
    let Some(ix) = model.solvers.iter().position(is_mg) else {
        return Ok(None);
    };
    // Every other solver would stay on the unsplit mesh.
    if model.solvers.len() > 1 {
        return Err(SolverError::Unsupported("MG Equal Split alongside other solvers".into()));
    }
    if ranks > 1 {
        return Err(SolverError::Unsupported("MG Equal Split in parallel".into()));
    }
    let p = &model.solvers[ix];
    let levels = p.integer("mg levels").or_else(|| p.integer("multigrid levels")).unwrap_or(1).max(1) as usize;
    let err = |e: MeshLevelsError| match e {
        MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        MeshLevelsError::Invalid(m) => SolverError::fatal("SplitMeshEqual", m),
    };
    let mut coarse: Vec<Box<dyn MeshView>> = Vec::with_capacity(levels - 1);
    let mut fine = LevelsMesh::from_view(mesh);
    for _ in 1..levels {
        let next = split_mesh_equal(&fine).map_err(err)?;
        coarse.push(Box::new(std::mem::replace(&mut fine, next)));
    }
    Ok(Some((ix + 1, coarse, fine)))
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
    discont_perm: Option<&[u32]>,
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
        discont_perm,
        edges: None,
    };
    let (projectors, stop) = generate_periodic_projectors(&env, case_dir, 0)?;
    if stop {
        return Err(SolverError::Unsupported("Save Projector And Stop".into()));
    }
    Ok(Some(PeriodicAssembly::new(model.boundary_conditions.clone(), projectors, mesh.n_nodes())))
}

/// `Mesh % PeriodicPerm` / `PeriodicFlip` when a boundary condition has
/// `Conforming BC` (`GeneratePeriodicProjectors`): `elmer-mortar`'s node
/// map, which the solvers with `Apply Conforming BCs` number their dofs by.
pub fn conforming_map(
    model: &Model,
    mesh: &dyn MeshView,
    coord_dim: usize,
    log: &dyn Logger,
) -> Result<Option<elmer_solver::conforming::ConformingMap>, SolverError> {
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
        discont_perm: None,
        edges: None,
    };
    Ok(generate_conforming_map(&env)?)
}

/// `SetNormalizedKeywords` (MainUtils.F90:226) for the `-distribute`
/// keywords: `CalculateEntityWeights` (SolverBasics.F90:10444) sums each
/// body's, body force's, material's (bulk elements) and boundary
/// condition's (boundary elements) measure, point by point over the
/// elements' default Gauss rule, and `ListSetParameters` gives each tagged
/// keyword the coefficient `1 / weight` that every `ListGet*` multiplies
/// its value by (`F = coeff * F`). Refused: a tagged value that is not a
/// constant, and coordinate systems other than Cartesian.
fn apply_distribute(model: &mut Model, mesh: &dyn MeshView, coords: CoordinateSystem) -> Result<(), SolverError> {
    use elmer_solver::model::DistEntity as D;
    if model.distributed.is_empty() {
        return Ok(());
    }
    if coords.kind != CoordinateKind::Cartesian {
        return Err(SolverError::Unsupported("-distribute keywords outside Cartesian coordinates".into()));
    }
    let (nbc, nbody, nmat, nbf) =
        (model.boundary_conditions.len(), model.bodies.len(), model.materials.len(), model.body_forces.len());
    let (mut wbc, mut wbody, mut wmat, mut wbf) = (vec![0.0f64; nbc], vec![0.0f64; nbody], vec![0.0f64; nmat], vec![0.0f64; nbf]);
    use elmer_core::ElementBasis;
    let basis = elmer_fe::LagrangeBasis::new(coords.dim);
    let mut ev = elmer_core::BasisEval::default();
    let nb = mesh.n_bulk();
    for t in 0..nb + mesh.n_boundary() {
        let el = mesh.element(elmer_core::ElemId(t as u32));
        let (mut body, mut bf, mut mat, mut bc) = (0usize, 0usize, 0usize, 0usize);
        if t < nb {
            let Some(b) = model.body_of_tag(el.tag) else { continue };
            body = b;
            let p = &model.bodies[b - 1];
            bf = p.integer("body force").map_or(0, |i| i.max(0) as usize);
            mat = p.integer("material").map_or(0, |i| i.max(0) as usize);
        } else {
            let Some(b) = model.bc_of_tag(el.tag) else { continue };
            bc = b;
        }
        let def = elmer_fe::element_type(el.etype.0).map_err(|e| SolverError::fatal("CalculateEntityWeights", e.to_string()))?;
        let rule = elmer_fe::quadrature::gauss_points(el.etype, Some(def.gauss_points))
            .map_err(|e| SolverError::fatal("CalculateEntityWeights", e.to_string()))?;
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
        for i in 0..rule.len() {
            basis
                .eval(el.etype, &xyz, [rule.u[i], rule.v[i], rule.w[i]], &mut ev)
                .map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let s = 1.0 * ev.det_j * rule.s[i];
            if (1..=nbc).contains(&bc) {
                wbc[bc - 1] += s;
            }
            if (1..=nbody).contains(&body) {
                wbody[body - 1] += s;
            }
            if (1..=nmat).contains(&mat) {
                wmat[mat - 1] += s;
            }
            if (1..=nbf).contains(&bf) {
                wbf[bf - 1] += s;
            }
        }
    }
    for (kind, i, key) in std::mem::take(&mut model.distributed) {
        let (w, sec) = match kind {
            D::Body => (wbody.get(i - 1), model.bodies.get_mut(i - 1)),
            D::BodyForce => (wbf.get(i - 1), model.body_forces.get_mut(i - 1)),
            D::BoundaryCondition => (wbc.get(i - 1), model.boundary_conditions.get_mut(i - 1)),
            D::Material => (wmat.get(i - 1), model.materials.get_mut(i - 1)),
        };
        let (Some(&w), Some(sec)) = (w, sec) else { continue };
        // `Refusing division with zero!`: the keyword keeps its value.
        if w <= f64::MIN_POSITIVE {
            continue;
        }
        let coeff = 1.0 / w;
        match sec.0.get_mut(&key) {
            Some(Value::Real { data, .. }) => data.iter_mut().for_each(|x| *x *= coeff),
            Some(_) => return Err(SolverError::Unsupported(format!("-distribute {key} that is not a constant"))),
            None => {}
        }
    }
    Ok(())
}

/// The `fem/tests` case name `elmer-fixtures` keys its ports by: the case
/// directory's name (the `.sif`'s directory; empty means the working one).
pub fn fixture_case(case_dir: &Path) -> String {
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
/// Whether the run assembles through a periodic backend.
fn periodic_in_use(sim: &Simulation) -> bool {
    sim.model.boundary_conditions.iter().any(|b| b.get("periodic bc").is_some())
}

pub fn run_case_with(case: Case<'_>, log: &dyn Logger, env: &RunEnv<'_>) -> Result<Outcome, SolverError> {
    run_case_meshes(case, Vec::new(), log, env)
}

/// [`run_case`] with the meshes the case's solvers name (`Solver i :: Mesh`),
/// as the front end read them.
pub fn run_case_with_meshes(
    case: Case<'_>,
    solver_meshes: Vec<(String, Box<dyn MeshView>)>,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<Outcome, SolverError> {
    if comm.size() > 1 {
        return Err(SolverError::Unsupported("more than one MPI rank".into()));
    }
    run_case_meshes(case, solver_meshes, log, &RunEnv { job: comm, globals: Vec::new() })
}

fn run_case_meshes(
    case: Case<'_>,
    solver_meshes: Vec<(String, Box<dyn MeshView>)>,
    log: &dyn Logger,
    env: &RunEnv<'_>,
) -> Result<Outcome, SolverError> {
    let (rank, size) = (env.job.rank(), env.job.size());
    let comm: &dyn Communicator = if size == 1 { env.job } else { &SerialComm };
    let Case { mut model, reloads, mesh, case_dir, mesh_dir } = case;
    let globals = job_globals(&mut model.simulation, rank, size)?;
    let fixture_case = fixture_case(&case_dir);
    let registry = solver_registry_for(&fixture_case);
    read_target_names(&mesh_dir.join("mesh.names"), &model.simulation, &mut model.bodies, &mut model.boundary_conditions)
        .map_err(|e| SolverError::fatal("ReadTargetNames", e.to_string()))?;
    let mut discont_perm: Option<Vec<u32>> = None;
    let transformed = transform_mesh(&mut model, mesh, &mut discont_perm)?;
    // ReadSurfaceDirector: ShellSolver's `Mesh % Name/mesh.director`.
    // ParticleOutputSolver writes next to the mesh (`SolverOutputDirectory(
    // ..., UseMeshDir = .TRUE. )`).
    for params in &mut model.solvers {
        let out = procedure_of(params).is_some_and(|(lib, entry)| {
            lib.eq_ignore_ascii_case(elmer_physics_particle::OUTPUT_PROCEDURE.0)
                && entry.eq_ignore_ascii_case(elmer_physics_particle::OUTPUT_PROCEDURE.1)
        });
        if out {
            params.insert(
                elmer_physics_particle::MESH_DIRECTORY_KEY,
                Value::Str(mesh_dir.to_string_lossy().into_owned()),
            );
        }
    }
    for params in &mut model.solvers {
        if params.logical("adaptive mesh refinement") == Some(true) {
            params.insert(
                elmer_adaptive::OUTPUT_DIRECTORY_KEY,
                Value::Str(case_dir.to_string_lossy().into_owned()),
            );
        }
    }
    let director = mesh_dir.join("mesh.director");
    if director.exists() {
        for params in &mut model.solvers {
            let shell = procedure_of(params).is_some_and(|(lib, entry)| {
                lib.eq_ignore_ascii_case(elmer_physics_shell::PROCEDURE.0) && entry.eq_ignore_ascii_case(elmer_physics_shell::PROCEDURE.1)
            });
            if !shell {
                continue;
            }
            if transformed.is_some() {
                return Err(SolverError::Unsupported("ShellSolver: mesh.director on a transformed mesh".into()));
            }
            params.insert("Mesh Director File", Value::File(director.to_string_lossy().into_owned()));
        }
    }
    let mesh: &dyn MeshView = match &transformed {
        Some(m) => m,
        None => mesh,
    };
    let (mg_coarse, mg_fine) = match multigrid_meshes(&model, mesh, size)? {
        Some((id, coarse, fine)) => (Some((id, coarse)), Some(fine)),
        None => (None, None),
    };
    if mg_fine.is_some() {
        discont_perm = None;
    }
    let mesh: &dyn MeshView = match &mg_fine {
        Some(m) => m,
        None => mesh,
    };
    let mesh = MovableMesh::new(mesh);
    let _motion = mesh.install();
    let mesh: &dyn MeshView = &mesh;
    let mut outputs = prepare_vtu_output(&mut model, &case_dir, &mesh_dir, size > 1)?;
    // The Output File directory: `Mesh % Name` under the case (the mesh
    // directory itself when the Mesh DB is "." "<name>").
    let result_dir = match model.header.string("mesh name") {
        Some(n) => case_dir.join(n),
        None => mesh_dir.clone(),
    };
    let coords = coordinate_system(&model, mesh)?;
    apply_distribute(&mut model, mesh, coords)?;
    legacy_solvers(&mut model, coords);
    lumped_matrix_inits(&mut model);
    drop_never_savedata(&mut model);
    preflight(&model, &registry)?;
    check_mesh_and_elements(&model)?;
    check_reloads(&model, &reloads, &registry)?;
    install_fe(coords.dim)?;

    let periodic = periodic_assembly(&model, mesh, coords.dim, &case_dir, log, discont_perm.as_deref())?;
    let conforming = conforming_map(&model, mesh, coords.dim, log)?;

    // The test-local plugins' `<Proc>_init` section edits (before Elmer
    // creates their variables).
    elmer_fixtures::registry::init_sections(&fixture_case, &mut model);
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_datatofield::is_section(p)) {
        elmer_physics_datatofield::init_section(p);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_advdiff::model_pde_cutfem::is_section(p)) {
        elmer_physics_advdiff::model_pde_cutfem::init_section(p);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_isosurface::is_section(p)) {
        elmer_physics_isosurface::init_section(p);
    }
    let transient_sim = model
        .simulation
        .string("simulation type")
        .is_some_and(|t| t.trim().eq_ignore_ascii_case("transient"));
    for p in &mut model.solvers {
        elmer_physics_magneto::solver::whitney_init_section(p, &model.body_forces, &model.boundary_conditions, transient_sim);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_magneto::mgdyn2d_steady::is_section(p)) {
        elmer_physics_magneto::mgdyn2d_steady::init_section(p);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_magneto::statmag::is_section(p)) {
        elmer_physics_magneto::statmag::init_section(p);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_magneto::emport::is_section(p)) {
        elmer_physics_magneto::emport::init_section(p, &mut model.boundary_conditions);
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_magneto::emwave::is_section(p)) {
        elmer_physics_magneto::emwave::init_section(p);
    }
    for i in 0..model.solvers.len() {
        if elmer_physics_magneto::emwave::is_calc_section(&model.solvers[i]) {
            elmer_physics_magneto::emwave::calc_init_section(&mut model.solvers, i);
        }
    }
    for p in model.solvers.iter_mut().filter(|p| elmer_physics_magneto::vhnodal::is_section(p)) {
        elmer_physics_magneto::vhnodal::init_section(p);
    }
    for i in 0..model.solvers.len() {
        if elmer_physics_magneto::calcfields2d::is_section(&model.solvers[i]) {
            elmer_physics_magneto::calcfields2d::init_section(&mut model, i, coords.dim as usize);
        }
    }
    let mut sim = Simulation::new(model);
    let mesh_names: Vec<String> = solver_meshes.iter().map(|(n, _)| n.clone()).collect();
    for (name, m) in solver_meshes {
        sim.add_solver_mesh(&name, Arc::from(m));
    }
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
    sim.backends.conforming = conforming.map(Arc::new);
    // Mortar contact (`Apply Contact BCs`, elmer-mortar's DetermineContact).
    if sim.model.solvers.iter().any(|s| s.logical("apply contact bcs") == Some(true)) {
        sim.backends.contact = Some(elmer_mortar::contact::MortarContact::new(case_dir.clone(), coords.dim));
    }
    sim.set_matc(Rc::new(SifMatc));
    for (name, value) in globals.iter().chain(&env.globals) {
        sim.vars.set_scalar(name, *value);
    }
    elmer_fixtures::registry::install_live(&fixture_case, &mut sim, mesh);
    // UMAT_StVenant_2D's `MyUMAT.F90` is UMATLib's `stvenant_kirchhoff`
    // verbatim (only a `USE Messages` line differs).
    if fixture_case == "UMAT_StVenant_2D" {
        elmer_physics_elasticsolve::umat::register_umat(
            "./MyUMAT",
            "stvenant_kirchhoff",
            elmer_physics_elasticsolve::umat::stvenant_kirchhoff,
        );
    }
    // UMAT_StVenant_axials' `LocalUMAT.F90` `stvenant_kirchhoff_axials` is
    // UMATLib's `stvenant_kirchhoff` written out for one shear component,
    // the axisymmetric `NSHR = 1`.
    if fixture_case == "UMAT_StVenant_axials" {
        elmer_physics_elasticsolve::umat::register_umat(
            "./LocalUMAT",
            "stvenant_kirchhoff_axials",
            elmer_physics_elasticsolve::umat::stvenant_kirchhoff,
        );
    }
    install_statelec_boundary(&mut sim, mesh);
    sim.setup(mesh, &registry, log, comm)?;
    // `SetInitialConditions`: the initial conditions (in `setup`), then the
    // Simulation's `Restart File` under `Mesh % Name`.
    elmer_solver::restart(&mut sim, &result_dir, log)?;
    if let Some(writer) = crate::result_file::writer(&sim.model, &result_dir, rank, size, log) {
        sim.set_save_current(Box::new(writer));
    }
    if sim.model.solvers.iter().any(|p| p.logical("adaptive mesh refinement") == Some(true)) {
        // Adaptive mesh refinement: each refined mesh saves to its own
        // directory (`Mesh % Name // '/' // Output File`).
        // VTU outputs pick their fields and the mesh at each visit, so they
        // follow the swaps.
        if periodic_in_use(&sim) || size > 1 {
            return Err(SolverError::Unsupported(
                "adaptive mesh refinement with periodic assembly or more than one partition".into(),
            ));
        }
        let case_dir2 = case_dir.clone();
        let log_level = log;
        let _ = log_level;
        // `Displace Mesh` and the mesh mappers move the mesh the run is on.
        sim.set_mesh_wrapper(Box::new(|m| {
            let om = crate::motion::OwnedMovableMesh::new(m);
            om.install();
            Arc::new(om) as Arc<dyn MeshView + Send + Sync>
        }));
        // `LoadMesh2` of an external mesher's output (`Adaptive Remesh`).
        let simulation = sim.model.simulation.clone();
        sim.set_mesh_loader(Box::new(move |dir: &Path| {
            crate::mesh_db::load(dir, &simulation)
                .map_err(|e| SolverError::fatal("LoadMesh", format!("{}: {e}", dir.display())))
        }));
        sim.set_mesh_writer_factory(Box::new(move |name: &str, sim: &elmer_solver::Simulation| {
            let dir = case_dir2.join(name);
            let _ = std::fs::create_dir_all(&dir);
            Ok(crate::result_file::writer(&sim.model, &dir, 0, 1, &elmer_core::NullLogger)
                .map(|w| Box::new(w) as Box<dyn elmer_solver::sim::SaveCurrent>))
        }));
    }
    // Each solver mesh saves its own `<mesh name>/<Output File>`.
    for (k, name) in mesh_names.iter().enumerate() {
        if let Some(writer) = crate::result_file::writer(&sim.model, &case_dir.join(name), rank, size, log) {
            sim.set_solver_mesh_writer(k + 1, Box::new(writer));
        }
    }
    if let Some((id, coarse)) = mg_coarse {
        sim.set_multigrid_meshes(id, coarse)?;
    }
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
    // ResultOutputSolver sections with `STL Format` only.
    for slot in &mut sim.slots {
        let params = &sim.model.solvers[slot.id - 1];
        if crate::output::is_stl_output(params) {
            slot.set_instance(Box::new(crate::stl::StlSolver { params: params.clone(), case_dir: case_dir.clone() }));
        }
    }
    // `Apply Mortar BCs`, or a `Periodic BC Use Lagrange Coefficient` of its
    // variable: the slot's solver behind the mortar constraints.
    let discont_perm = discont_perm.map(Arc::new);
    for slot in &mut sim.slots {
        let params = &sim.model.solvers[slot.id - 1];
        if crate::openfoam::is_section(params) {
            slot.set_instance(Box::new(crate::openfoam::Elmer2OpenFoamWrite { case_dir: case_dir.to_path_buf() }));
            continue;
        }
        if crate::gmsh_reader::is_section(params) {
            slot.set_instance(Box::new(crate::gmsh_reader::GmshOutputReader { case_dir: case_dir.to_path_buf() }));
            continue;
        }
        // SaveProjection uses the projectors itself, without constraints.
        if crate::saveprojection::is_section(params) {
            slot.set_instance(Box::new(crate::saveprojection::SaveProjection {
                coord_dim: coords.dim,
                case_dir: case_dir.to_path_buf(),
                discont_perm: discont_perm.clone(),
            }));
            continue;
        }
        if !crate::mortar::applies(params) && !crate::mortar::lagrange_periodic(params, &sim.model.boundary_conditions) {
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
            discont_perm: discont_perm.clone(),
            jump: Default::default(),
        }));
    }
    let mut report = sim.run(mesh, log, comm)?;
    for block in &reloads {
        sim.reload(block, mesh, &registry, log, comm)?;
        rebuild_changed(&mut sim, block, &registry)?;
        // The Fortran reads `Periodic BC Offset` at every assembly, so a later
        // block's offset applies from that block on (periodic_homogenize).
        if let Some(assembly) = sim.backends.assembly.as_mut() {
            assembly.reload_boundary_conditions(&sim.model.boundary_conditions);
        }
        // `ReloadInputFile` sets `Mesh % SavesDone = 0` (ElmerSolver.F90:501):
        // each RUN block starts the Output File afresh.
        if let Some(writer) = crate::result_file::writer(&sim.model, &result_dir, rank, size, log) {
            sim.set_save_current(Box::new(writer));
        }
        report = sim.run(mesh, log, comm)?;
    }
    let written = std::mem::take(&mut *written.lock().unwrap());
    let passed = finalize_reference(&sim.reference, &case_dir, log, env.job)?;
    let norms = (0..sim.slots.len()).filter_map(|ix| Some((sim.slots[ix].id, sim.slot_norm(ix)?))).collect();
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
