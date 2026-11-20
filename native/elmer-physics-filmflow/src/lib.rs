//! `FilmFlowSolver` (FilmFlowSolver.F90): the reduced-dimensional
//! Navier-Stokes equations of a thin film or channel of `Gap Height`, with
//! laminar (Poiseuille), Darcy or Manning friction, on the MINI element of
//! `_Init0` (`p:1` with velocity bubbles, condensed).
//!
//! Ported: `Model Dimension` 1 and 2, `Axi Symmetric`, `GradP
//! Discretization`, `Lateral Strain`, `Convect` (with Newton), the
//! `Friction Model`s, gravity on the `Bedrock Height` slope, the body forces
//! (`FilmFlow Bodyforce i`, `Normal Velocity`, `Fsi Velocity`, `Flow
//! Admittance` with `External FilmPressure`), the boundary loads
//! (`FilmPressure i`, `Mass Flux`), `Gap Addition Factor`, `Calculate
//! Friction Heating` (the `FilmFriction HeatFlux` and `Heating Energy`
//! fields) and `Use Heating Source`.
//!
//! Refused: artificial compressibility, `Use Gap Average`, `Calculate
//! Pressure Heating`, `Calculate Sensitivity`, `Skip First Solution`, and a
//! non-nodal element beyond the MINI bubbles.
#![forbid(unsafe_code)]

pub mod kernel;
mod solver;

pub use solver::{factory, register, solver_defaults, FilmFlowSolver, PROCEDURE};

#[derive(Debug, thiserror::Error)]
pub enum FilmFlowError {
    #[error("element basis failed: {0}")]
    Basis(String),
    #[error("invalid input: {0}")]
    Invalid(String),
}
