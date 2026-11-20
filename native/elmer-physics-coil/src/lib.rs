//! `CoilSolver` (CoilSolver.F90): the current density of a coil, from a
//! dummy potential `CoilPot` that is 1 at the `Coil Start` and 0 at the
//! `Coil End` boundary, scaled to the `Desired Coil Current` by the nodal
//! loads of the ends. The current `CoilCurrent` is the Galerkin projection
//! of the conductivity times the potential gradient, optionally normalised
//! to a uniform density (`Normalize Coil Current`). With `Coil Conductivity
//! Fix` the conductivity is first fixed by the inverse of the gradient
//! magnitude (`CoilFix`).
//!
//! Ported: one coil (the solver section), open (`Coil Start` / `Coil End`)
//! or closed (`Coil Closed`: two potentials cut on either half of the coil,
//! the coil centre and the inertia fit of its axis by `DSYEV`, and
//! `PotSelect` choosing between them),
//! `Desired Coil Current`, `Desired Current Density`, `Coil Cross Section`,
//! `Normalize Coil Current`, `Coil Conductivity Fix`, `Coil Anisotropic`,
//! `CFix Diffusion`, `Electric Conductivity Name`, `Use Unity
//! Conductivity`, `Calculate Coil Current` and `Save Coil Set`.
//!
//! In a partitioned run the coil centre and inertia, the fixing sets
//! (`ParallelSumVectorInt` max and min over the holders of a node), the
//! end currents (owned nodes) and the counts are reduced across ranks.
//!
//! Refused: `Single Coil Cut` / `Test Coil Cut`, elemental
//! current fields (`Calculate Elemental Fields`, `Fix Input Current
//! Density`), `Use Wall Distance`, `Coil Helicity` and `Save Coil Index`.
#![forbid(unsafe_code)]
// The ports index several parallel arrays with the Fortran loop variables.
#![allow(clippy::needless_range_loop)]

mod fit;
mod solver;

pub use solver::{factory, register, solver_defaults, CoilSolver, PROCEDURE};
