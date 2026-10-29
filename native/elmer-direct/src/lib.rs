//! Direct linear solvers for the Rust-native Elmer rewrite.
//!
//! * [`sparse`]: `FaerLu` and `FaerCholesky` (UMFPACK / MUMPS / SuperLU /
//!   PARDISO / CHOLMOD replacements) on faer's sparse factorisations.
//! * [`band`]: Elmer/LAPACK band storage with band LU (`DGBTRF`/`DGBTRS`) and
//!   band Cholesky (`DPBTRF`/`DPBTRS`), transliterated from the reference
//!   algorithms (`SolveBand.F90`, `SolveSBand.F90`, `BandMatrix.F90`).
//! * [`dense`]: thin faer adapters for the LAPACK dense routines Elmer imports.
//! * [`select_direct`]: `Linear System Direct Method` dispatch
//!   (`DirectSolve.F90:4266`).
#![forbid(unsafe_code)]

pub mod band;
pub mod dense;
pub mod scalar;
pub mod select;
pub mod sparse;

pub use band::{BandCholesky, BandLu, BandMatrix, SymBandMatrix};
pub use scalar::DirectScalar;
pub use select::{kind_for_method, select_direct, try_select_direct, try_select_direct_for, DirectKind};
pub use sparse::{set_single_thread, FaerCholesky, FaerLu};
