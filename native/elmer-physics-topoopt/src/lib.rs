//! `TopoOpt` (TopoOpt.F90): SIMP topology optimisation of the compliance of
//! a primary solver. The primary solver (`Primary Solver Index`, default 1)
//! keeps its local matrices (`Local Matrix Storage`) and scales them by the
//! elemental `topo mult = emin + (1-emin) rho^p` (`Matrix Multiplier
//! Name`); each visit computes the element compliances `u_e^T K_e u_e`, the
//! sensitivities, filters them (`Filter Method` density / sensitivity /
//! none with a `Filter Type` simple, distance or PDE filter), updates the
//! densities by the optimality-criteria bisection for the `Volume Fraction`
//! and the black-and-white field `topo bw` (`Create BW Topology`).
//!
//! Refused: `Solve Adjoint Problem`, `Use Adjoint Solution`, constraint
//! modes (`No Modes`), `Periodic PhysSolver`, `Periodic Timesteps` cycles,
//! `Create Zero Levelset` and more than one partition.
#![forbid(unsafe_code)]
// The ports index several parallel arrays with the Fortran loop variables.
#![allow(clippy::needless_range_loop)]

mod filter;
mod solver;

pub use solver::{factory, register, solver_defaults, TopoOpt, PROCEDURE};
