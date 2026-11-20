//! `CoilSolver` (CoilSolver.F90): the current density of a coil, from a
//! dummy potential `CoilPot` that is 1 at the `Coil Start` and 0 at the
//! `Coil End` boundary, scaled to the `Desired Coil Current` by the nodal
//! loads of the ends. The current `CoilCurrent` is the Galerkin projection
//! of the conductivity times the potential gradient, optionally normalised
//! to a uniform density (`Normalize Coil Current`). With `Coil Conductivity
//! Fix` the conductivity is first fixed by the inverse of the gradient
//! magnitude (`CoilFix`).
//!
//! Ported: open coils of one coil (the solver section, no coil components),
//! `Desired Coil Current`, `Desired Current Density`, `Coil Cross Section`,
//! `Normalize Coil Current`, `Coil Conductivity Fix`, `Coil Anisotropic`,
//! `CFix Diffusion`, `Electric Conductivity Name`, `Use Unity
//! Conductivity`, `Calculate Coil Current` and `Save Coil Set`.
//!
//! Refused: closed coils (`Coil Closed`, the coil cuts and their inertia
//! fit), elemental
//! current fields (`Calculate Elemental Fields`, `Fix Input Current
//! Density`), `Use Wall Distance`, `Coil Helicity`, `Save Coil Index` and
//! more than one partition.
#![forbid(unsafe_code)]
// The ports index several parallel arrays with the Fortran loop variables.
#![allow(clippy::needless_range_loop)]

mod solver;

pub use solver::{factory, register, solver_defaults, CoilSolver, PROCEDURE};
