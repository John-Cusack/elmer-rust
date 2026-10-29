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
//!    `ElementInfo` / `GaussPoints` build them for p elements, with the
//!    p-only rules of [`quadrature`];
//! 6. [`dirichlet`]: the edge-DOF part of `DefaultDirichletBCs`;
//! 7. [`setup`]: all of it for one solver section from its keywords
//!    ([`PSetup`]), the object a solver framework adapter holds;
//! 8. [`solver`]: `elmer-solver`'s `SolverDofs` for [`PSetup`], the
//!    [`PDofsBackend`] to install in `Backends::dofs`, and the PoissonPFEM
//!    test plugin.
//!
//! Element systems are glued with `elmer_assembly::CsrAssembler` on the
//! row indices [`dofs::DofLayout::element_rows`] gives; `elmer-solver`
//! reaches all of it through `Backends::dofs`.
//!
//! Supported: first-order meshes; the p basis and rules for lines,
//! triangles, quadrilaterals (also serendipity), tetrahedra
//! (`ConvertToACTetra`), pyramids (full elements only, as in Elmer), wedges
//! and bricks; global bubbles or none; scalar or multi-component fields.
//! Element definitions are read from the solver section, the equations
//! (`Element`, `Element{id}`) and the bodies (`Solver id: Element`) with
//! Elmer's precedence. Refused or not covered: DG, radiation, periodic and
//! mortar couplings, and boundary elements a BC `Body Id` makes active.
#![forbid(unsafe_code)]

pub mod basis;
pub mod defs;
pub mod dirichlet;
pub mod dofs;
pub mod perm;
pub mod pmesh;
pub mod quadrature;
pub mod setup;
pub mod solver;

pub use defs::{initial_max_defs, ElementDefs, Kind, MaxDefs};
pub use dofs::DofLayout;
pub use perm::{create_p_system, PSystem};
pub use pmesh::{MatcDegree, PMesh, PMeshOptions, SolverElements};
pub use setup::{global_bubbles, is_p_section, ModelFacts, OtherSolver, PSetup};
pub use solver::{register_poisson_pfem, PDofsBackend, PoissonPfem};

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
