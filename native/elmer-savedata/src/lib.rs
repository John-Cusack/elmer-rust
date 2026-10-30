//! `SaveScalars` and `SaveLine` of Elmer's `SaveData` module as scheduled
//! [`Solver`] hooks.
//!
//! Both are ordinary solvers of the coupled-system loop in `elmer-solver`
//! (`Procedure = "SaveData" "SaveScalars"` / `"SaveLine"`, placed by `Exec
//! Solver` and `Exec Intervals`); there is no second simulation loop. Files
//! are written through `elmer-saveio` (`.dat`, `.names`), so the layout is the
//! one the Fortran writes.
//!
//! | Hook | Fortran (`fem/src/modules/SaveData/`) | What it does |
//! |---|---|---|
//! | [`SaveScalars`] | `SaveScalars.F90` | one row per visit of nodal statistics, bulk and boundary integrals, counts, point values and Simulation-section scalars; exports `Target Variable n` as a global variable and every value to the Simulation section |
//! | [`SaveLine`] | `SaveLine.F90` | the nodes of the boundary conditions flagged `Save Line`, ordered along the line, with coordinates and field values |
//!
//! # What is implemented
//!
//! `SaveScalars` operators: `sum`, `sum abs`, `sum square`, `mean`, `mean
//! abs`, `mean square`, `rms`, `variance`, `deviation`, `range`, `max`, `max
//! abs`, `min`, `min abs` (nodal, optionally masked with `body `, `body
//! force ` or `material ` and `Nodal Variable n`); `int`, `int mean`, `int
//! square`, `int square mean`, `int rms`, `int abs`, `int abs mean`, `int
//! variance`, `volume`, `potential energy`, `diffusive energy` (element
//! integrals, Cartesian and axisymmetric); `boundary sum`, `boundary dofs`,
//! `boundary max`, `boundary max abs`, `boundary min`, `boundary min abs`,
//! `boundary mean`, `boundary int`, `boundary int mean`, `area`; `nodes`,
//! `elements`, `bounding box`, `partitions`, `norm`, `nonlin change`,
//! `steady state change`, `nonlin iter`, `nonlin converged`, `steady
//! converged`, `dofs`, `nans`. Also `Save Points`, `Save Coordinates`
//! (nearest node), the `res:` scalars other solvers leave in the Simulation
//! section, `Target Variable n`, `Show Norm Index/Name`, `Reference Values`,
//! `Cost Function Index`, `Line Marker`, `File Append`, `Live Graph`
//! (`.csv`) and `Partition Numbering`.
//!
//! `SaveLine`: existing boundaries only (`Save Mask`), `Optimize Node
//! Ordering`, `Line Marker`, `Skip Boundary Info`, `File Append`.
//!
//! # What is refused
//!
//! Anything else in the two sections fails when the solver is built, with
//! [`SolverError::Unsupported`]
//! naming the keyword: flux reconstruction (`diffusive flux`, `convective
//! flux`, `Save Heat Flux`), material `Coefficient n`, `positive `/`negative `
//! operators, polylines, circles, isocurves and axis lines, exact point
//! coordinates, eigenvalue columns, MATC `Expression n`, run-dependent
//! values (`cpu time`, `wall time`, `cpu memory`, `threads`), JSON output,
//! numbered file names, failed-strategy marker files and every MPI
//! reduction (a run with more than one partition stops at the first visit).
//! The sections' cross-field requirements - a nodal field, Cartesian or axi
//! symmetric coordinates for integrals - are checked at the visit.
//!
//! # Where Elmer's behaviour is odd, and what is copied
//!
//! The port is bug-compatible where the pristine ElmerSolver's output
//! depends on it, and says so in the code:
//! - `AddToSaveList` overwrites the host-associated `VariableName`, so the
//!   column names of the second and later boundary conditions of `boundary
//!   int`, `boundary int mean` and `area` carry no variable name;
//! - `variance` is the standard deviation, `boundary min abs` and `boundary
//!   max abs` compare the wrong way round, and a vector `boundary max` sums
//!   the components without squaring them;
//! - `Show Norm Index` makes the hook's own global variable
//!   `<Equation>_var` carry the value.
//!
//! # Wiring
//!
//! [`register`] adds both procedures to a
//! [`SolverRegistry`]. The factories are plain
//! `fn` pointers, so the element basis and quadrature the integrals need
//! come from [`install_fe`], and the version text of the `.names` files from
//! [`install_names_meta`]. Statistics, counts, point values and `SaveLine`
//! need neither.
#![forbid(unsafe_code)]

mod integrals;
mod line;
mod materials;
mod plan;
mod scalars;
mod view;

use std::sync::OnceLock;

use elmer_core::{ElementBasis, ParamMap, Quadrature};
use elmer_saveio::savedata::NamesMeta;
use elmer_solver::{Solver, SolverError, SolverRegistry};

pub use line::SaveLine;
pub use materials::SaveMaterials;
pub use plan::{LinePlan, ScalarPlan};
pub use scalars::SaveScalars;

/// `Procedure = "SaveData" "SaveScalars"`.
pub const SAVE_SCALARS: (&str, &str) = ("SaveData", "SaveScalars");
/// `Procedure = "SaveData" "SaveLine"`.
pub const SAVE_LINE: (&str, &str) = ("SaveData", "SaveLine");
/// `Procedure = "SaveData" "SaveMaterials"`.
pub const SAVE_MATERIALS: (&str, &str) = ("SaveData", "SaveMaterials");

/// Builds the element basis and quadrature the integrals use.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();
static META: OnceLock<fn() -> NamesMeta> = OnceLock::new();

/// Installs the basis and quadrature provider. Returns `false` if one is
/// already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Installs the header the `.names` files carry (`Elmer version`, the `.sif`
/// name, a fixed `File started at:` for reproducible tests). Returns `false`
/// if one is already installed.
pub fn install_names_meta(meta: fn() -> NamesMeta) -> bool {
    META.set(meta).is_ok()
}

fn default_meta() -> NamesMeta {
    NamesMeta {
        version: "9.0".into(),
        ..NamesMeta::default()
    }
}

fn names_meta() -> fn() -> NamesMeta {
    META.get().copied().unwrap_or(default_meta)
}

/// Registers `("SaveData", "SaveScalars")`, `("SaveData", "SaveLine")` and
/// `("SaveData", "SaveMaterials")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(SAVE_SCALARS.0, SAVE_SCALARS.1, scalars_factory);
    registry.register(SAVE_LINE.0, SAVE_LINE.1, line_factory);
    registry.register(SAVE_MATERIALS.0, SAVE_MATERIALS.1, materials_factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`SaveScalars`]: reads
/// the section once and refuses what is not implemented.
pub fn scalars_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let plan = ScalarPlan::from_params(params)?;
    let fe = FE.get().map(|p| p()).map(|(basis, quad)| integrals::Fe { basis, quad });
    Ok(Box::new(SaveScalars::new(plan, fe, names_meta())))
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`SaveLine`].
pub fn line_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SaveLine::new(LinePlan::from_params(params)?, names_meta())))
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`SaveMaterials`].
pub fn materials_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SaveMaterials::from_params(params)?))
}
