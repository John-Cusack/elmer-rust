//! `ParticleDynamics` (ParticleDynamics.F90) on the particle machinery of
//! ParticleUtils.F90: particles initialised in the mesh's bounding box
//! (`box random` from the run's `RANDOM_NUMBER` stream, or `box uniform
//! cubic`), located by marching face to face through a 2-D mesh, moved by
//! gravity, damping and the electric (`Potential Variable Name`) and
//! magnetic (`Magnetic Field Variable Name`) fields, wrapped in a periodic
//! box, with path integrals of body-force sources and the particle data
//! summed to nodal `Particle <field>` variables (`Particle To Field`).
//!
//! Refused: groups, collisions and contacts, walls and the particle
//! boundary conditions, scattering, Runge-Kutta, variable particle steps,
//! velocity fields (drag), release schedules, masks and weights, the other
//! initialisation methods, `Output Interval` (SaveParticleData), 3-D meshes
//! and more than one partition. `ParticleStatistics` and
//! `ParticleInformation` only write Info lines and are not ported.
#![forbid(unsafe_code)]
// The ports keep the Fortran's loops and status tests.
#![allow(clippy::needless_range_loop, clippy::manual_range_contains)]

mod faces;
mod particles;
mod solver;

pub use solver::{factory, register, solver_defaults, ParticleDynamics, PROCEDURE};
