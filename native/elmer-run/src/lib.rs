//! Serial end-to-end composition of the Rust-native Elmer: a parsed case and
//! its mesh go through `elmer-solver`'s `Simulation` with the merged
//! providers installed, and come out as solver norms, the
//! `CompareToReferenceSolution` verdict (`TEST.PASSED`) and the `Post File`
//! VTU.
//!
//! * [`compose`]: the solver registry (HeatSolve) and `elmer-direct` for
//!   `Linear System Solver = Direct`, next to `elmer-solver`'s default
//!   `elmer-assembly` and `elmer-timeint` backends.
//! * [`case`]: [`run_case`], the `ExecSimulation` path from a [`Case`] to an
//!   [`Outcome`], refusing unsupported input before any side effect.
//! * [`output`]: the `Post File = *.vtu` auto-solver, written with `elmer-vtu`.
//! * [`engine`]: [`NativeEngine`], the `elmer-cli` [`Engine`] over a
//!   [`Frontend`] that parses the `.sif` and loads the mesh (S1's
//!   `elmer-sif`/`elmer-mesh`/`elmer-fe`).
//!
//! This crate owns no numerics; anything a provider lacks is refused here
//! rather than filled in.
//!
//! [`Engine`]: elmer_cli::driver::Engine
#![forbid(unsafe_code)]

pub mod case;
pub mod compose;
pub mod engine;
pub mod output;

pub use case::{run_case, Case, Outcome};
pub use compose::{install_backends, solver_registry};
pub use engine::{Frontend, LoadedCase, NativeEngine};
