//! ElmerGUI's process-boundary contract as executable checks.
//!
//! ElmerGUI never links the solver. It probes `ElmerSolver -v` (and
//! `ElmerSolver_mpi -v`), writes a case directory plus
//! `ELMERSOLVER_STARTINFO`, starts `ElmerSolver` with no arguments, scrapes
//! `Time:` and `ComputeChange` lines from the log, and loads the files the
//! case wrote. This crate reproduces each step against any solver executable
//! given by absolute path: the pristine Fortran build or the Rust-native
//! `ElmerSolver` binaries of `elmer-cli`.
//!
//! Linux "Compile" is not implemented in this ElmerGUI fork and is not
//! modelled here.
#![forbid(unsafe_code)]

pub mod check;
pub mod edf;
pub mod fixtures;
pub mod logparse;
pub mod mesh;
pub mod run;
pub mod sif;
pub mod typecheck;
pub mod vtu;
