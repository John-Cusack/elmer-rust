//! The p-element degree-of-freedom layer: what ElmerSolver adds on top of
//! a nodal mesh before a solver with `Element = "p:..."` can assemble.
//!
//! The pipeline, each step a port of the named Elmer routine:
//! 1. [`defs`]: the `Element` keyword into `Solver % Def_Dofs`
//!    (`GetDefs`, `GetMaxDefs`);
//! 2. [`pmesh`]: mesh edges, p definitions, edge and bubble DOF counts and
//!    the boundary-element orientation (`NonNodalElements`,
//!    `FindMeshEdges2D`, `SetMeshEdgeFaceDOFs`, `AssignLocalNumber`,
//!    `SetMeshMaxDOFs`);
//! 3. [`dofs`]: the global DOF slots of one element (`mGetElementDOFs`,
//!    `mGetBoundaryIndexesFromParent`);
//! 4. [`perm`]: the permutation and matrix pattern over those slots
//!    (`InitialPermutation`, `MakeListMatrix`, `CreateMatrix`), reusing
//!    `elmer-assembly`'s bandwidth optimisation and CSR setup;
//! 5. [`basis`]: the `elmer-fe-hp` basis and Gauss rule of one element as
//!    `ElementInfo` / `GaussPoints` build them for p elements;
//! 6. [`dirichlet`]: the edge-DOF part of `DefaultDirichletBCs`.
//!
//! Element systems are glued with `elmer_assembly::CsrAssembler` on the
//! row indices [`dofs::DofLayout::element_rows`] gives, so no assembly
//! hook is needed.
//!
//! Supported: first-order 1D and 2D meshes (`202`, `303`, `404` bulk
//! elements), global bubbles or none, scalar or multi-component fields.
//! Refused or not covered: 3D meshes (faces), DG, radiation, periodic and
//! mortar couplings, and element definitions in equation sections.
#![forbid(unsafe_code)]

pub mod basis;
pub mod defs;
pub mod dirichlet;
pub mod dofs;
pub mod perm;
pub mod pmesh;

pub use defs::{initial_max_defs, ElementDefs, Kind, MaxDefs};
pub use dofs::DofLayout;
pub use perm::{create_p_system, PSystem};
pub use pmesh::{PMesh, PMeshOptions, SolverElements};

/// Errors of the p-element setup.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    /// An `Element` definition that cannot be parsed or evaluated.
    #[error("bad element definition: {0}")]
    BadDefinition(String),
    /// A case outside what this crate ports.
    #[error("unsupported: {0}")]
    Unsupported(String),
}
