//! Rust ports of the Fortran plugins that Elmer's own test suite compiles
//! next to a case (`ADD_ELMERTEST_MODULE` in `fem/tests/*/CMakeLists.txt`).
//!
//! The Rust-native Elmer loads no shared objects, so a case that names
//! `Procedure "Source" "Source"` needs the function registered under that
//! name. Ported fixtures stay excluded from the held-out set (strategy
//! decision D12): they exist so those cases can run, not as product code.
//!
//! * [`udf`]: scalar `Procedure` functions (`MySource`, `InitDist`, the
//!   1dtests `Source`, `AlternatingSource`), registered into an
//!   [`elmer_expr::UdfRegistry`].
//! * [`poisson_assembly`]: the `BulkAssembly` / `BoundaryAssembly` element
//!   kernels shared by the nine `CoupledPoisson*` cases, and [`hooks`], the
//!   same kernels as `elmer-solver` [`ElementAssembly`] routines on
//!   `elmer-fe`.
//! * [`poisson_solver`]: the `Poisson` `PoissonSolver` solver procedure of
//!   `1dtests` and `ExtrudedMesh`; [`tresca`]: the `tresca` case's
//!   `TrescaSolver`; [`single_elements`]: the `ElementUnitTest_h1` element
//!   checks and [`single_edge_elements`] their H(curl) version;
//!   [`model_pde_handle`]: `ModelPDEhandle`'s `AdvDiffSolver`;
//!   [`model_pde_splitted`]: ModelPDE as `_bulk`/`_boundary` steps;
//!   [`poisson_ml`]: the matrix-free p-element Poisson of `PoissonML`;
//!   [`benchmark_shell`]: the p-element cylindrical shell of the shell cases;
//!   [`stokes_pfem`]: `StokesPFEM`'s p-element Stokes solver;
//!   [`stokes_proj`]: `StokesProj`'s velocity and pressure solvers;
//!   [`bentonite`]: the `bentonite` case's coupled 3-dof solver;
//!   [`l2nd_derivatives`]: the basis and second-derivative check of
//!   `L2ndDerivatives`.
//! * [`registry`]: case-qualified keys. Plugin libraries share names across
//!   cases (`Source` is four different functions), so a fixture is selected
//!   by case, library and entry together, and every `install` builds fresh
//!   per-instance state like a fresh process would.
//!
//! [`ElementAssembly`]: elmer_solver::ElementAssembly
//!
//! The README lists every other quick-test plugin with its blockers.
#![forbid(unsafe_code)]

pub mod benchmark_shell;
pub mod bentonite;
pub mod hooks;
pub mod l2nd_derivatives;
pub mod linalg;
pub mod model_pde_handle;
pub mod model_pde_splitted;
pub mod poisson_assembly;
pub mod poisson_ml;
pub mod poisson_solver;
pub mod registry;
pub mod single_edge_elements;
pub mod single_elements;
pub mod stokes_pfem;
pub mod stokes_proj;
pub mod tresca;
pub mod udf;

use elmer_expr::FieldLookup;
use elmer_solver::LiveModel;

/// The part of Elmer's `Model_t` a scalar fixture reads: node coordinates
/// (`Model % Nodes`), the node count, `Time` and nodal variable values
/// (`VariableGet` + `Perm`).
///
/// The view is live: a fixture asks it again on every call, as the Fortran
/// reads `Model` on every call, so the solver side keeps one shared object
/// whose state follows the run.
pub trait ModelView: FieldLookup + Send + Sync {
    /// `Model % NumberOfNodes`.
    fn n_nodes(&self) -> usize;
}

/// The simulation's own view, kept current by `elmer-solver`.
impl ModelView for LiveModel {
    fn n_nodes(&self) -> usize {
        LiveModel::n_nodes(self)
    }
}

/// A fixture stopped the run, or was called in a state its Fortran leaves
/// undefined.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FixtureError {
    /// `CALL Fatal(caller, message)`: Elmer prints this and halts.
    #[error("{caller}: {message}")]
    Fatal { caller: &'static str, message: &'static str },
    /// The Fortran indexes outside an array or returns an unset result here.
    #[error("{caller}: {what}")]
    Undefined { caller: &'static str, what: String },
}
