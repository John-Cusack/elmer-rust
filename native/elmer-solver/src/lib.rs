//! Solver framework for the Rust-native Elmer: the [`Simulation`] that owns
//! the model sections, the field registry and the solver slots; the
//! [`Solver`] trait physics crates implement and the [`SolverCtx`] they get;
//! Elmer's norm and change semantics (`ComputeNorm`, `ComputeChange`); the
//! coupled-system scheduler; and `CompareToReferenceSolution`.
//!
//! Matrix structure/assembly, local time discretisation and linear solvers
//! are pluggable backends ([`backend`]). Assembly defaults to
//! `elmer-assembly` ([`ElmerAssembly`]) and local time discretisation to
//! `elmer-timeint` ([`ElmerTime`]); keyword values that are not constants
//! are evaluated with `elmer-expr` ([`expr`]): `Procedure` values through the
//! simulation's function registry, with [`LiveModel`] as the model view user
//! functions read, and MATC values through one shared engine
//! ([`Simulation::set_matc`]). `Assembly Solver` sections run [`AssemblySolver`]
//! over a registered `Bulk Assembly Procedure` ([`ElementAssembly`]).
//! Ported behaviour cites the Fortran it
//! follows; options not ported yet are refused with
//! [`SolverError::Unsupported`] or [`SolverError::NotWired`], never
//! approximated.
#![forbid(unsafe_code)]

pub mod assembly;
pub mod backend;
mod components;
mod block_alias;
mod block_driver;
mod block_system;
pub mod change;
pub mod conforming;
mod child;
mod nt;
pub use nt::{boundary_normal, tangent_directions};
mod amg;
mod constraint_modes;
mod constraints;
mod eigen;
mod control;
pub use control::lu_solve;
pub mod coupled;
pub mod error;
mod exported;
pub use exported::{dg_edges, dg_interfaces, parse_elemental_gauss_rules, DgIndexes};
mod global_time;
pub mod expr;
pub mod fortran_fmt;
pub mod linsys;
pub mod live;
mod loads;
mod lumped;
mod matvec;
mod monolithic;
mod multigrid;
pub mod model;
pub mod norm;
mod normalize;
pub mod reference;
pub mod remesh;
pub mod restart;
pub mod sim;
pub mod solver;
pub mod time;
pub mod variable;

pub use assembly::ElmerAssembly;
/// Elmer's file formats and Fortran edit descriptors, for the run's writers.
pub use elmer_saveio;
pub use restart::{load_restart, restart};
/// The dense symmetric eigensolver (`DSYEV`) and the eigen analysis, for
/// solver crates that depend on this one only.
pub use elmer_eigen;
pub use backend::{
    AssemblyBackend, Backends, DofsBackend, FeBackend, LinearSolverRegistry, Perm, SolverDofs, TimeBackend, TimeTerms,
};
pub use coupled::{AssemblySolver, CoupledSystem};
pub use constraint_modes::{constraint_modes_lumping, lumping_init, Lumping, CAPACITANCE, CONDUCTANCE, CONDUCTIVITY};
pub use constraints::ConstraintMatrix;
pub use change::{compute_change, ChangeEnv, ChangeReport, SolverState};
pub use error::SolverError;
pub use expr::{SharedMatc, TextMatc};
pub use live::LiveModel;
pub use matvec::crs_matrix_vector_multiply;
pub use model::{BodyParams, Model};
pub use norm::{compute_norm, NormSpec};
pub use reference::{compare_to_reference, finalize_reference, ReferenceTally};
pub use time::ElmerTime;
pub use sim::{ExecWhen, LinearSystem, RunReport, SaveCurrent, Simulation, SolverSlot};
pub use solver::{AssemblyFactory, ElementAssembly, Solver, SolverCtx, SolverFactory, SolverRegistry};
pub use variable::{FieldMut, Strided, StridedMut, VarHandle, VarId, Variable, VariableStore};
