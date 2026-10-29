//! Rust ports of the Fortran plugins that Elmer's own test suite compiles
//! next to a case (`ADD_ELMERTEST_MODULE` in `fem/tests/*/CMakeLists.txt`).
//!
//! The Rust-native Elmer loads no shared objects, so a case that names
//! `Procedure "Source" "Source"` needs the function registered under that
//! name. Ported fixtures stay excluded from the held-out set (strategy
//! decision D12): they exist so those cases can run, not as product code.
//!
//! * [`udf`]: scalar `Procedure` functions (`MySource`, `InitDist`, the
//!   1dtests `Source`), registered into an [`elmer_expr::UdfRegistry`].
//! * [`poisson_assembly`]: the `BulkAssembly` / `BoundaryAssembly` element
//!   kernels shared by the nine `CoupledPoisson*` cases.
//! * [`registry`]: case-qualified keys. Plugin libraries share names across
//!   cases (`Source` is four different functions), so a fixture is selected
//!   by case, library and entry together, and every `install` builds fresh
//!   per-instance state like a fresh process would.
//!
//! The README lists every other quick-test plugin with its blockers.
#![forbid(unsafe_code)]

pub mod poisson_assembly;
pub mod registry;
pub mod udf;

use elmer_expr::FieldLookup;

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
