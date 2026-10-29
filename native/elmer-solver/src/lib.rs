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
pub mod change;
pub mod coupled;
pub mod error;
pub mod expr;
pub mod fortran_fmt;
pub mod linsys;
pub mod live;
pub mod model;
pub mod norm;
pub mod reference;
pub mod sim;
pub mod solver;
pub mod time;
pub mod variable;

pub use assembly::ElmerAssembly;
pub use backend::{AssemblyBackend, Backends, FeBackend, LinearSolverRegistry, Perm, TimeBackend, TimeTerms};
pub use coupled::{AssemblySolver, CoupledSystem};
pub use change::{compute_change, ChangeEnv, ChangeReport, SolverState};
pub use error::SolverError;
pub use expr::{SharedMatc, TextMatc};
pub use live::LiveModel;
pub use model::{BodyParams, Model};
pub use norm::{compute_norm, NormSpec};
pub use reference::{compare_to_reference, finalize_reference, ReferenceTally};
pub use time::ElmerTime;
pub use sim::{ExecWhen, LinearSystem, RunReport, Simulation, SolverSlot};
pub use solver::{AssemblyFactory, ElementAssembly, Solver, SolverCtx, SolverFactory, SolverRegistry};
pub use variable::{FieldMut, Strided, StridedMut, VarHandle, VarId, Variable, VariableStore};
