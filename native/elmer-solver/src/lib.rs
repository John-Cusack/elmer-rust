//! Solver framework for the Rust-native Elmer: the [`Simulation`] that owns
//! the model sections, the field registry and the solver slots; the
//! [`Solver`] trait physics crates implement and the [`SolverCtx`] they get;
//! Elmer's norm and change semantics (`ComputeNorm`, `ComputeChange`); the
//! coupled-system scheduler; and `CompareToReferenceSolution`.
//!
//! Matrix structure/assembly, local time discretisation and linear solvers
//! are pluggable backends ([`backend`]) so this crate does not depend on
//! the lanes that implement them. Ported behaviour cites the Fortran it
//! follows; options not ported yet are refused with
//! [`SolverError::Unsupported`] or [`SolverError::NotWired`], never
//! approximated.
#![forbid(unsafe_code)]

pub mod backend;
pub mod change;
pub mod error;
pub mod fortran_fmt;
pub mod model;
pub mod norm;
pub mod reference;
pub mod sim;
pub mod solver;
pub mod variable;

pub use backend::{AssemblyBackend, Backends, LinearSolverRegistry, Perm, TimeBackend, TimeTerms};
pub use change::{compute_change, ChangeEnv, ChangeReport, SolverState};
pub use error::SolverError;
pub use model::{BodyParams, Model};
pub use norm::{compute_norm, NormSpec};
pub use reference::{compare_to_reference, finalize_reference, ReferenceTally};
pub use sim::{ExecWhen, LinearSystem, RunReport, Simulation, SolverSlot};
pub use solver::{Solver, SolverCtx, SolverFactory, SolverRegistry};
pub use variable::{FieldMut, Strided, StridedMut, VarHandle, VarId, Variable, VariableStore};
