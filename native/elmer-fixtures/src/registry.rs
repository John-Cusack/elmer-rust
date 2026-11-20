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

use crate::adv_react_dg_p;
use crate::benchmark_shell;
use crate::block_lin_elast;
use crate::bentonite;
use crate::edge_element_test;
use crate::edge_fe_interpolation;
use crate::hooks::{PoissonBoundary, PoissonBulk};
use crate::hrot;
use crate::l2nd_derivatives;
use crate::model_pde_handle;
use crate::model_pde_splitted;
use crate::poisson_ml;
use crate::poisson_solver;
use crate::compute_error;
use crate::face_fe_interpolation;
use crate::count_solver;
use crate::sun_angle;
use crate::mesh_colour;
use crate::curved_bndry_pfem;
use crate::single_edge_elements;
use crate::single_elements;
use crate::q1q0;
use crate::stokes_pfem;
use crate::stokes_proj;
use crate::tresca;
use crate::udf::{init_dist_fn, mixed_source_2d_fn, mixed_source_3d_fn, my_source_fn, stdout_sink, AlternatingSource, OneDSource, OneDTolerances};
use crate::ModelView;

type Build = fn(Arc<dyn ModelView>) -> Udf;

/// The `fem/tests` directory of a ctest name: `ADD_ELMER_TEST(name NPROCS 1
/// 3)` runs one directory as `name` and `name_np3`, so the `_np<N>` suffix
/// is dropped.
fn source_case(case: &str) -> &str {
    match case.rsplit_once("_np") {
        Some((base, n)) if !base.is_empty() && !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => base,
        _ => case,
    }
}

/// `(case, library, entry, builder)`. Library and entry are the strings the
/// case's `Procedure` keyword uses.
const FIXTURES: &[(&str, &str, &str, Build)] = &[
    // `CurvedBndryPFEM2/Poisson.f90`'s `CircX` / `CircY`: the displacement
    // that takes a point radially onto the circle of radius sqrt(2).
    ("CurvedBndryPFEM2", "./Poisson", "CircX", |_| Udf::Const(Arc::new(curved_bndry_pfem::circ_x))),
    ("CurvedBndryPFEM2", "./Poisson", "CircY", |_| Udf::Const(Arc::new(curved_bndry_pfem::circ_y))),
    ("SD_CurvedBndryPFEM2", "./Poisson", "CircX", |_| Udf::Const(Arc::new(curved_bndry_pfem::circ_x))),
    ("SD_CurvedBndryPFEM2", "./Poisson", "CircY", |_| Udf::Const(Arc::new(curved_bndry_pfem::circ_y))),
    ("ConstantParamFunc", "MySource", "MySource", |_| Udf::Node(my_source_fn())),
    ("AdvDiffFCT", "InitDist", "InitDist", |m| Udf::Node(init_dist_fn(m))),
    // `Source.F90` of the 2D mixed Poisson tests: zero, and the source of
    // the manufactured `u = 16 (x - x^2)(y - y^2)`.
    ("MixedPoisson2D", "Source", "SourceFun", |_| Udf::Node(Arc::new(|_, _| 0.0))),
    ("RaviartThomas2D", "Source", "SourceFun", |m| Udf::Node(mixed_source_2d_fn(m))),
    ("BDM2D", "Source", "SourceFun", |m| Udf::Node(mixed_source_2d_fn(m))),
    ("NormalTangentialDisplacement", "TangentForce", "TangentForce", |m| {
        Udf::Node(crate::tangent_force::tangent_force_fn(m))
    }),
    // The 3D ones: zero, and the source of `u = -64 (x^2 - 1/4)(y^2 -
    // 1/4)(z^2 - 1/4)`.
    ("MixedPoisson3D", "Source", "SourceFun", |_| Udf::Node(Arc::new(|_, _| 0.0))),
    ("MixedPoisson3D_Brick", "Source", "SourceFun", |_| Udf::Node(Arc::new(|_, _| 0.0))),
    ("RaviartThomas3D", "Source", "SourceFun", |m| Udf::Node(mixed_source_3d_fn(m))),
    ("1dtests", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::ONE_D_TESTS).into_fn())
    }),
    ("1dtests_MeshLevels", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::MESH_LEVELS).into_fn())
    }),
    // Its argument `cpar` exists only through `Apply Explicit Control`,
    // which `elmer-run` refuses by name until it is implemented; with it
    // missing, `elmer-expr` would evaluate the argument as zero.
    // Serial: rank 0 prints.
    ("HeatControlExplicit", "AlternatingSource", "AlternatingSource", |_| {
        Udf::Node(AlternatingSource::new(true, stdout_sink()).into_fn())
    }),
];

/// The `(library, entry)` pairs ported for `case`, empty when none are.
pub fn procedures(case: &str) -> Vec<(&'static str, &'static str)> {
    let case = source_case(case);
    FIXTURES.iter().filter(|f| f.0 == case).map(|f| (f.1, f.2)).collect()
}

/// A new instance of the fixture `case` uses for `Procedure "library" "entry"`,
/// bound to `model`. `None` when that plugin is not ported.
pub fn lookup(case: &str, library: &str, entry: &str, model: Arc<dyn ModelView>) -> Option<Udf> {
    let case = source_case(case);
    FIXTURES
        .iter()
        .find(|f| f.0 == case && f.1 == library && f.2.eq_ignore_ascii_case(entry))
        .map(|f| f.3(model))
}

/// Registers every fixture of `case` into `registry` under the names the
/// case's SIF uses, as new instances. Returns how many were registered; zero
/// means the case has no ported plugins.
pub fn install(case: &str, model: Arc<dyn ModelView>, registry: &UdfRegistry) -> usize {
    let case = source_case(case);
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
    // The same source as `_h1`, on `Element = "p:1"`.
    ("ElementUnitTest_p1", "SingleElements", "SingleElements", Hook::Solver(single_elements::factory)),
    ("SD_ElementUnitTest_p1", "SingleElements", "SingleElements", Hook::Solver(single_elements::factory)),
    ("ModelPDEhandle", "ModelPDEhandle", "AdvDiffSolver", Hook::Solver(model_pde_handle::factory)),
    // No `AdvDiffSolver` in the library: Elmer runs `_bulk` and `_boundary` in steps.
    ("ModelPDEsplitted", "ModelPDEsplitted", "AdvDiffSolver", Hook::Solver(model_pde_splitted::factory)),
    ("SD_ModelPDEsplitted", "ModelPDEsplitted", "AdvDiffSolver", Hook::Solver(model_pde_splitted::factory)),
    ("PoissonML", "PoissonML", "PoissonSolver", Hook::Solver(poisson_ml::factory)),
    ("StokesPFEM", "Stokes1", "StokesSolver", Hook::Solver(stokes_pfem::factory)),
    ("Q1Q0", "Q1Q0", "StokesSolver", Hook::Solver(q1q0::factory)),
    ("SlaveSolverSlots", "CountSolver", "CountSolver", Hook::Solver(count_solver::factory)),
    // The same source as `CountSolver2.F90`.
    ("SlaveSolverSlotsScan", "CountSolver2", "CountSolver", Hook::Solver(count_solver::factory)),
    ("AdvReactDG_P", "AdvReactDG_P", "AdvectionReactionSolver", Hook::Solver(adv_react_dg_p::factory)),
    ("BlockLinElast1", "BlockLinElast", "BulkAssembly", Hook::Assembly(|| Box::new(block_lin_elast::BlockLinElastBulk))),
    ("SD_BlockLinElast1", "BlockLinElast", "BulkAssembly", Hook::Assembly(|| Box::new(block_lin_elast::BlockLinElastBulk))),
    // One `BenchmarkShell.F90` in all five shell cases.
    ("p-FEM_two_solvers", "BenchmarkShell", "ShellSolver", Hook::Solver(benchmark_shell::factory)),
    ("SD_p-FEM_two_solvers", "BenchmarkShell", "ShellSolver", Hook::Solver(benchmark_shell::factory)),
    ("p-FEM_with_varying_p", "BenchmarkShell", "ShellSolver", Hook::Solver(benchmark_shell::factory)),
    ("SD_p-FEM_with_varying_p", "BenchmarkShell", "ShellSolver", Hook::Solver(benchmark_shell::factory)),
    ("SD_Classical2DShell", "BenchmarkShell", "ShellSolver", Hook::Solver(benchmark_shell::factory)),
    ("L2ndDerivatives", "L2ndDerivativesTest", "TestSolver", Hook::Solver(l2nd_derivatives::factory)),
    ("StokesProj", "Pressure", "PressureSolver", Hook::Solver(stokes_proj::pressure_factory)),
    ("StokesProj", "Velocity", "VelocitySolver", Hook::Solver(stokes_proj::velocity_factory)),
    ("bentonite", "BentoniteSolver", "BentoniteSolver", Hook::Solver(bentonite::factory)),
    ("ElementUnitTest_edge", "SingleEdgeElements", "SingleEdgeElements", Hook::Solver(single_edge_elements::factory)),
    ("ElementUnitTest_edge_piola", "SingleEdgeElements", "SingleEdgeElements", Hook::Solver(single_edge_elements::factory)),
    (
        "ElementUnitTest_edge_quadratic",
        "SingleEdgeElements",
        "SingleEdgeElements",
        Hook::Solver(single_edge_elements::factory),
    ),
    ("RaviartThomas2D", "ComputeError", "IntegrateError", Hook::Solver(compute_error::factory)),
    ("BDM2D", "ComputeError", "IntegrateError", Hook::Solver(compute_error::factory)),
    ("RaviartThomas3D", "ComputeError", "IntegrateError", Hook::Solver(compute_error::factory_3d)),
    ("SecondOrderEdgeElement2D", "EdgeElementTest", "EdgeElementSolver", Hook::Solver(edge_element_test::factory_2d)),
    ("HelmholtzEdge", "Hrot", "HrotSolver", Hook::Solver(hrot::factory)),
    ("HelmholtzFace", "Helmholtz", "HelmholtzHDIVSolver", Hook::Solver(hrot::factory_hdiv)),
    ("SecondOrderEdgeElement3D", "EdgeElementTest", "EdgeElementSolver", Hook::Solver(edge_element_test::factory_3d)),
    (
        "SecondOrderEdgeElement2D_BCs",
        "EdgeFEInterpolationTest",
        "BestApproximationSolver",
        Hook::Solver(edge_fe_interpolation::factory),
    ),
    ("FaceElement2D_BCs", "FaceFEInterpolationTest", "BestApproximationSolver", Hook::Solver(face_fe_interpolation::factory)),
    // fem/src/modules/SunAngle.F90, built into the case as `./SunAngle`.
    ("SunAngle", "./SunAngle", "SunAngleSolver", Hook::Solver(sun_angle::factory)),
    ("MeshColour", "MeshColour", "MeshColour", Hook::Solver(mesh_colour::factory)),
    ("FaceElement3D_BCs", "FaceFEInterpolationTest", "BestApproximationSolver", Hook::Solver(face_fe_interpolation::factory)),
    // `PoissonPFEM/Poisson.f90` on p-elements, ported in `elmer-pdofs`.
    ("PoissonPFEM", "Poisson", "PoissonSolver", Hook::Solver(elmer_pdofs::solver::poisson_pfem)),
    // `CurvedBndryPFEM/Poisson.f90`: the mesh moved onto a circle first.
    ("CurvedBndryPFEM", "./Poisson", "PoissonSolver", Hook::Solver(curved_bndry_pfem::factory)),
    ("SD_CurvedBndryPFEM", "./Poisson", "PoissonSolver", Hook::Solver(curved_bndry_pfem::factory)),
    // `CurvedBndryPFEM2/Poisson.f90`: the plain p Poisson on the mesh a p
    // MeshSolve has curved.
    ("CurvedBndryPFEM2", "./Poisson", "PoissonSolver", Hook::Solver(curved_bndry_pfem::factory2)),
    ("SD_CurvedBndryPFEM2", "./Poisson", "PoissonSolver", Hook::Solver(curved_bndry_pfem::factory2)),
    // The same source up to blank lines.
    ("SD_PoissonPFEM", "Poisson", "PoissonSolver", Hook::Solver(elmer_pdofs::solver::poisson_pfem)),
    // The same `Poisson.f90` as SD_PoissonPFEM's, on solver-specific meshes.
    ("SolverSpecificMeshes", "Poisson", "PoissonSolver", Hook::Solver(elmer_pdofs::solver::poisson_pfem)),
    ("SD_SolverSpecificMeshes", "Poisson", "PoissonSolver", Hook::Solver(elmer_pdofs::solver::poisson_pfem)),
];

/// Edits a solver section.
type SectionInit = fn(&mut ParamMap);

/// `<Proc>_init` routines that only edit their solver section, before its
/// variable exists: `(case, library, entry, edit)`.
const SECTION_INITS: &[(&str, &str, &str, SectionInit)] = &[
    ("ElementUnitTest_h1", "SingleElements", "SingleElements", single_elements::init_section),
    ("SunAngle", "SunAngle", "SunAngleSolver", sun_angle::init_section),
    ("MeshColour", "MeshColour", "MeshColour", mesh_colour::init_section),
    ("ElementUnitTest_p1", "SingleElements", "SingleElements", single_elements::init_section),
    ("SD_ElementUnitTest_p1", "SingleElements", "SingleElements", single_elements::init_section),
    ("ElementUnitTest_edge", "SingleEdgeElements", "SingleEdgeElements", single_edge_elements::init_section),
    ("HelmholtzEdge", "Hrot", "HrotSolver", hrot::init_section),
    ("HelmholtzFace", "Helmholtz", "HelmholtzHDIVSolver", hrot::init_section),
    ("ElementUnitTest_edge_piola", "SingleEdgeElements", "SingleEdgeElements", single_edge_elements::init_section),
    ("ElementUnitTest_edge_quadratic", "SingleEdgeElements", "SingleEdgeElements", single_edge_elements::init_section),
    (
        "SecondOrderEdgeElement2D_BCs",
        "EdgeFEInterpolationTest",
        "BestApproximationSolver",
        edge_fe_interpolation::init_section,
    ),
];

fn hooks(case: &str) -> impl Iterator<Item = (&'static str, &'static str, Hook)> + '_ {
    let case = source_case(case);
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
    let case = source_case(case);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mpi_test_names_use_their_directory() {
        assert_eq!(source_case("HeatControlExplicit_np3"), "HeatControlExplicit");
        assert_eq!(source_case("passive_np8"), "passive");
        for keep in ["HeatControlExplicit", "x_np", "_np3", "a_npx", "1dtests"] {
            assert_eq!(source_case(keep), keep);
        }
        assert_eq!(procedures("HeatControlExplicit_np3"), procedures("HeatControlExplicit"));
        assert!(!procedures("HeatControlExplicit_np3").is_empty());
    }
}
