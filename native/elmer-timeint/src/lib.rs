//! Time discretisation and time-step control for the Rust-native Elmer.
//!
//! Ports of TimeIntegrate.F90 (local and CRS forms), the time parts of
//! `Add1stOrderTime`/`Add2ndOrderTime` (SolverBasics.F90) and the step loop of
//! `ExecSimulation` (ElmerSolver.F90) with its adaptive and predictor-corrector
//! controllers. Everything is a pure function or state machine: the solution
//! history lives with the solver (`elmer-solver`) and is passed in as slices.
#![forbid(unsafe_code)]
// The kernels are line-for-line transliterations of TimeIntegrate.F90: indexed
// loops and `s = s + x` keep each statement comparable with the Fortran it
// ports, and the evaluation order is what the bitwise tests pin down.
#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]

pub mod adaptive;
pub mod crs;
pub mod error;
pub mod local;
pub mod scheme;
pub mod second;
pub mod stepper;

pub use adaptive::{
    Decision, PredCorrConfig, PredCorrControl, Substep, SubstepConfig, SubstepController,
    SubstepMode,
};
pub use crs::{
    add_1st_order_time_crs, bdf_crs, fractional_step_crs, newmark_beta_crs, runge_kutta_crs,
    vbdf_crs, CrsTimeData,
};
pub use error::TimeError;
pub use local::{add_1st_order_time, lump_matrix, AdamsState, Aux1st};
pub use scheme::{effective_pred_corr_order, FsCoeffs, Scheme};
pub use second::{
    add_2nd_order_time, time_2nd_order_crs, update_2nd_order, Newmark2, SecondOrderState,
};
pub use stepper::{fortran_es, geometric_timesteps, RunKind, Step, TimeStepper};
