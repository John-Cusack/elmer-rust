//! Sparse-matrix layer of the Rust-native Elmer rewrite: everything Elmer does
//! to a CRS matrix except solve it.
//!
//! * [`ops`]: `spmv` (Elmer's accumulator order), sorting, triplet assembly,
//!   element access, `glue_local` (`CRS_GlueLocalMatrix`), [`CsrOp`].
//! * [`scaling`]: `ScaleLinearSystem` / `BackScaleLinearSystem`.
//! * [`ilu`]: ILU(n), ILUT, incomplete Cholesky and the diagonal
//!   preconditioner, real and complex, with `CRS_LUSolve`.
//! * [`save`]: the `Linear System Save` file format, read and write.
//!
//! Every loop is a transliteration of the Fortran (CRSMatrix.F90,
//! MatrixScaling.F90, SolveCore.F90) with the same floating-point association
//! order, targeting bitwise agreement with the serial, non-MKL, non-OpenMP
//! build. Known gaps: complex division follows libgcc's classic `__divdc3`
//! without its GCC 12+ overflow pre-scaling; complex incomplete Cholesky and
//! the block-view (`BCols`) fast paths are not ported; matrices must have
//! strictly ascending columns and a stored diagonal for the factorizations.
#![forbid(unsafe_code)]
// The kernels index several arrays by one loop variable exactly as the Fortran does.
#![allow(clippy::needless_range_loop)]

pub mod ilu;
pub mod ops;
pub mod save;
pub mod scalar;
pub mod scaling;

pub use ilu::Ilu;
pub use ops::{
    add_to_element, contains, find_diag, find_entry, from_triplets, get_element, glue_local,
    glue_local_into, row_sum, row_sum_abs, set_element, sort_rows, spmv, spmv_blocked,
    spmv_transpose, spmv_with_values, zero_matrix, zero_row, CsrOp, NO_DIAG, NO_INDEX,
};
pub use save::{
    read_linear_system_save, write_linear_system_save, Columns, SaveOptions, SavedSystem,
};
pub use scalar::{ElmerScalar, AEPS};
pub use scaling::{scale_system, RowScaling, ScaleOptions, ScalingMethod};
