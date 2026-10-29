//! Case-qualified fixture selection.
//!
//! Test plugins are named by file, not by purpose: `Source` is the 1dtests
//! right-hand side, a manufactured source in `BDM2D` and `MixedPoisson*`, and
//! more, each with different arguments. A fixture is therefore keyed by
//! `(case, library, entry)` where `case` is the `fem/tests/<dir>` name, and
//! the caller picks the case it is running. Nothing is registered process-wide
//! on its own.
//!
//! Each call builds new fixture state, so two simulations (or two cases in
//! one test process) never share `SAVE`d Fortran variables.
//!
//! Scalar functions go into an [`UdfRegistry`] ([`install`], or
//! [`install_live`] for a [`Simulation`]); solver procedures and element
//! assembly routines into a [`SolverRegistry`] ([`register`]). A composition
//! running case `case` calls both, and [`init_sections`], before
//! `Simulation::setup`.

use std::sync::Arc;

use elmer_core::{MeshView, ParamMap};
use elmer_expr::{Udf, UdfRegistry};
use elmer_solver::solver::procedure_of;
use elmer_solver::{ElementAssembly, Model, Simulation, SolverFactory, SolverRegistry};

use crate::bentonite;
use crate::hooks::{PoissonBoundary, PoissonBulk};
use crate::model_pde_handle;
use crate::poisson_solver;
use crate::single_elements;
use crate::stokes_proj;
use crate::tresca;
use crate::udf::{init_dist_fn, my_source_fn, OneDSource, OneDTolerances};
use crate::ModelView;

type Build = fn(Arc<dyn ModelView>) -> Udf;

/// `(case, library, entry, builder)`. Library and entry are the strings the
/// case's `Procedure` keyword uses.
const FIXTURES: &[(&str, &str, &str, Build)] = &[
    ("ConstantParamFunc", "MySource", "MySource", |_| Udf::Node(my_source_fn())),
    ("AdvDiffFCT", "InitDist", "InitDist", |m| Udf::Node(init_dist_fn(m))),
    ("1dtests", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::ONE_D_TESTS).into_fn())
    }),
    ("1dtests_MeshLevels", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::MESH_LEVELS).into_fn())
    }),
    // `HeatControlExplicit`'s `AlternatingSource` (`udf::AlternatingSource`)
    // is ported but not listed: its argument `cpar` exists only through
    // `Apply Explicit Control`, which `elmer-solver` does not implement, and
    // `elmer-expr` evaluates a missing variable as zero where Elmer stops
    // (`Can't find independent variable`). Listing it would turn the case's
    // refusal into a wrong answer.
];

/// The `(library, entry)` pairs ported for `case`, empty when none are.
pub fn procedures(case: &str) -> Vec<(&'static str, &'static str)> {
    FIXTURES.iter().filter(|f| f.0 == case).map(|f| (f.1, f.2)).collect()
}

/// A new instance of the fixture `case` uses for `Procedure "library" "entry"`,
/// bound to `model`. `None` when that plugin is not ported.
pub fn lookup(case: &str, library: &str, entry: &str, model: Arc<dyn ModelView>) -> Option<Udf> {
    FIXTURES
        .iter()
        .find(|f| f.0 == case && f.1 == library && f.2.eq_ignore_ascii_case(entry))
        .map(|f| f.3(model))
}

/// Registers every fixture of `case` into `registry` under the names the
/// case's SIF uses, as new instances. Returns how many were registered; zero
/// means the case has no ported plugins.
pub fn install(case: &str, model: Arc<dyn ModelView>, registry: &UdfRegistry) -> usize {
    let mut n = 0;
    for (c, lib, entry, build) in FIXTURES {
        if *c == case {
            registry.register_fortran(lib, entry, build(Arc::clone(&model)));
            n += 1;
        }
    }
    n
}

/// [`install`] with the simulation's own live model view and function
/// registry (`Simulation::live_model`, `Backends::udfs`). Call it before
/// `Simulation::setup`, so initial conditions can use the functions too.
pub fn install_live(case: &str, sim: &mut Simulation, mesh: &dyn MeshView) -> usize {
    let model: Arc<dyn ModelView> = sim.live_model(mesh);
    install(case, model, &sim.backends.udfs)
}

/// A solver procedure or an element assembly routine.
#[derive(Clone, Copy)]
enum Hook {
    Solver(SolverFactory),
    Assembly(fn() -> Box<dyn ElementAssembly>),
}

/// The `PoissonAssembly.f90` cases. The file is the same in all of them (up
/// to declaration order and unused `USE` lines).
const POISSON_ASSEMBLY_CASES: [&str; 11] = [
    "CoupledPoisson1",
    "CoupledPoisson2",
    "CoupledPoisson3",
    "CoupledPoisson4",
    "CoupledPoisson5",
    "CoupledPoisson6",
    "CoupledPoisson7",
    "CoupledPoisson8",
    "CoupledPoisson9",
    "BlockPoisson1",
    "BlockPoisson2",
];

/// `(case, library, entry, hook)` for the cases with one source per case.
const HOOKS: &[(&str, &str, &str, Hook)] = &[
    ("1dtests", "Poisson", "PoissonSolver", Hook::Solver(poisson_solver::nodal_factory)),
    ("1dtests_MeshLevels", "Poisson", "PoissonSolver", Hook::Solver(poisson_solver::nodal_factory)),
    ("ExtrudedMesh", "Poisson", "PoissonSolver", Hook::Solver(poisson_solver::condensing_factory)),
    ("ExtrudedMeshBaseline1", "Poisson", "PoissonSolver", Hook::Solver(poisson_solver::condensing_factory)),
    ("tresca", "Tresca", "TrescaSolver", Hook::Solver(tresca::factory)),
    ("ElementUnitTest_h1", "SingleElements", "SingleElements", Hook::Solver(single_elements::factory)),
    ("ModelPDEhandle", "ModelPDEhandle", "AdvDiffSolver", Hook::Solver(model_pde_handle::factory)),
    ("StokesProj", "Pressure", "PressureSolver", Hook::Solver(stokes_proj::pressure_factory)),
    ("StokesProj", "Velocity", "VelocitySolver", Hook::Solver(stokes_proj::velocity_factory)),
    ("bentonite", "BentoniteSolver", "BentoniteSolver", Hook::Solver(bentonite::factory)),
];

/// Edits a solver section.
type SectionInit = fn(&mut ParamMap);

/// `<Proc>_init` routines that only edit their solver section, before its
/// variable exists: `(case, library, entry, edit)`.
const SECTION_INITS: &[(&str, &str, &str, SectionInit)] =
    &[("ElementUnitTest_h1", "SingleElements", "SingleElements", single_elements::init_section)];

fn hooks(case: &str) -> impl Iterator<Item = (&'static str, &'static str, Hook)> + '_ {
    let assembly = POISSON_ASSEMBLY_CASES.contains(&case).then_some([
        ("PoissonAssembly", "BulkAssembly", Hook::Assembly(|| Box::new(PoissonBulk::default()))),
        ("PoissonAssembly", "BoundaryAssembly", Hook::Assembly(|| Box::new(PoissonBoundary::default()))),
    ]);
    assembly
        .into_iter()
        .flatten()
        .chain(HOOKS.iter().filter(move |h| h.0 == case).map(|h| (h.1, h.2, h.3)))
}

/// The `(library, entry)` solver procedures and assembly routines ported for
/// `case`, empty when none are.
pub fn hook_procedures(case: &str) -> Vec<(&'static str, &'static str)> {
    hooks(case).map(|(lib, entry, _)| (lib, entry)).collect()
}

/// Registers the solver procedures and element assembly routines of `case`
/// into `registry` under the names its SIF uses (`registry` ignores a
/// `./` directory part and a `.so` suffix). Every routine instance the
/// solver builds from them starts with fresh `SAVE`d state. Returns how many
/// were registered.
pub fn register(case: &str, registry: &mut SolverRegistry) -> usize {
    let mut n = 0;
    for (lib, entry, hook) in hooks(case) {
        match hook {
            Hook::Solver(f) => registry.register(lib, entry, f),
            Hook::Assembly(build) => registry.register_assembly(lib, entry, Arc::new(move |_| Ok(build()))),
        }
        n += 1;
    }
    n
}

/// Runs the `<Proc>_init` section edits of `case`'s plugins on the solver
/// sections that name them (Elmer calls `_init` before it creates the
/// solver's variable). Call it on the model before `Simulation::new`.
/// Returns how many sections were edited.
pub fn init_sections(case: &str, model: &mut Model) -> usize {
    let mut n = 0;
    for p in &mut model.solvers {
        let Some((lib, entry)) = procedure_of(p) else { continue };
        let lib = lib.rsplit('/').next().unwrap_or(&lib).to_string();
        let found = SECTION_INITS
            .iter()
            .find(|i| i.0 == case && i.1.eq_ignore_ascii_case(&lib) && i.2.eq_ignore_ascii_case(&entry));
        if let Some(init) = found {
            (init.3)(p);
            n += 1;
        }
    }
    n
}
