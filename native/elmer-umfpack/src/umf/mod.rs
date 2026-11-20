//! UMFPACK 4.4 (`umfpack_di_*`, real double, `int` indices).
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  All Rights
//! Reserved.  See the crate documentation for the licence.  Translated to
//! Rust; see the modification notice there.

pub mod analyze;
pub mod colamd;
pub mod consts;
pub mod singletons;
pub mod symbolic;
pub mod twobytwo;
pub mod util;

pub use consts::*;
pub use symbolic::{umfpack_qsymbolic, umfpack_symbolic, Symbolic};

/// `UMFPACK_defaults`.
pub fn umfpack_defaults(control: &mut [f64]) {
    for c in control.iter_mut().take(UMFPACK_CONTROL) {
        *c = 0.0;
    }
    control[UMFPACK_PRL] = UMFPACK_DEFAULT_PRL;
    control[UMFPACK_DENSE_ROW] = UMFPACK_DEFAULT_DENSE_ROW;
    control[UMFPACK_DENSE_COL] = UMFPACK_DEFAULT_DENSE_COL;
    control[UMFPACK_AMD_DENSE] = UMFPACK_DEFAULT_AMD_DENSE;
    control[UMFPACK_STRATEGY] = UMFPACK_DEFAULT_STRATEGY as f64;
    control[UMFPACK_2BY2_TOLERANCE] = UMFPACK_DEFAULT_2BY2_TOLERANCE;
    control[UMFPACK_AGGRESSIVE] = UMFPACK_DEFAULT_AGGRESSIVE;
    control[UMFPACK_PIVOT_TOLERANCE] = UMFPACK_DEFAULT_PIVOT_TOLERANCE;
    control[UMFPACK_SYM_PIVOT_TOLERANCE] = UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE;
    control[UMFPACK_BLOCK_SIZE] = UMFPACK_DEFAULT_BLOCK_SIZE;
    control[UMFPACK_ALLOC_INIT] = UMFPACK_DEFAULT_ALLOC_INIT;
    control[UMFPACK_FRONT_ALLOC_INIT] = UMFPACK_DEFAULT_FRONT_ALLOC_INIT;
    control[UMFPACK_SCALE] = UMFPACK_DEFAULT_SCALE as f64;
    control[UMFPACK_IRSTEP] = UMFPACK_DEFAULT_IRSTEP;
    control[UMFPACK_COMPILED_WITH_BLAS] = 0.0; // NBLAS build
    control[UMFPACK_COMPILED_FOR_MATLAB] = 0.0;
    control[UMFPACK_COMPILED_IN_DEBUG_MODE] = 0.0;
}
