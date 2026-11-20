//! A Rust port of the ARPACK routines Elmer's `ArpackEigenSolve` drives
//! (the implicitly restarted Lanczos method `DSAUPD` / `DSEUPD`; the
//! non-symmetric path's LAPACK routines are in [`lapack_ns`]), operation
//! for operation, with the reference BLAS and LAPACK routines they call, so
//! that eigenpairs come out as the pristine build's do to round-off.
//!
//! ARPACK is BSD-licensed (Rice University; Lehoucq, Sorensen, Yang);
//! LAPACK and the reference BLAS are BSD-licensed (University of Tennessee
//! et al.). This is a translation; no Fortran is compiled or linked.

pub mod blas;
pub mod lapack;
pub mod lapack_ns;
pub mod naupd;
pub mod neupd;
pub mod saupd;
pub mod seupd;
