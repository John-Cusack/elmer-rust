//! Serial end-to-end composition of the Rust-native Elmer: a parsed case and
//! its mesh go through `elmer-solver`'s `Simulation` with the merged
//! providers installed, and come out as solver norms, the
//! `CompareToReferenceSolution` verdict (`TEST.PASSED`) and the `Post File`
//! VTU.
//!
//! * [`compose`]: the solver registry (HeatSolve, FluxSolver, StressSolver,
//!   ElasticSolve, FlowSolve, IncompressibleNS, StatElecSolve, StatCurrentSolve, HelmholtzSolve,
//!   AdvectionDiffusion, ModelPDE, the mesh tools (MeshChecksum,
//!   RigidMeshMapper, StructuredMeshMapper, StructuredProjectToPlane),
//!   ThermoElectricSolver, WhitneyAVSolver, WhitneyAVHarmonicSolver,
//!   SaveScalars, SaveLine) and `elmer-fe` as
//!   their shared element library;
//!   assembly, time and linear solvers are `elmer-solver`'s default backends.
//! * [`matc`]: [`SifMatc`], run-time MATC on `elmer-sif`'s binding.
//! * [`case`]: [`run_case`], the `ExecSimulation` path from a [`Case`] to an
//!   [`Outcome`], refusing unsupported input before any side effect.
//! * `motion`: node coordinates the mesh mappers move.
//! * [`output`]: the `Post File = *.vtu` auto-solver, written with `elmer-vtu`.
//! * [`sif`]: [`SifFrontend`], `elmer-sif` as the case reader, with the mesh
//!   reader pluggable until `elmer-mesh` merges.
//! * [`unported`]: keywords whose feature no provider implements, refused
//!   before setup so they cannot turn into a wrong answer.
//! * [`engine`]: [`NativeEngine`], the `elmer-cli` [`Engine`] over a
//!   [`Frontend`] that parses the `.sif` and loads the mesh (S1's
//!   `elmer-sif`/`elmer-mesh`).
//!
//! This crate owns no numerics; anything a provider lacks is refused here
//! rather than filled in.
//!
//! [`Engine`]: elmer_cli::driver::Engine
#![forbid(unsafe_code)]

pub mod case;
pub mod compose;
pub mod engine;
pub mod matc;
mod motion;
pub mod output;
pub mod sif;
pub mod unported;

pub use case::{run_case, Case, Outcome};
pub use compose::{coordinate_system, fe_backend, fe_provider, install_fe, solver_registry, solver_registry_for};
pub use matc::SifMatc;
pub use engine::{Frontend, LoadedCase, NativeEngine};
pub use sif::{MeshLoader, SifFrontend};
