//! The compressed-row matrix shared by assembly, sparse ops and solvers.
//!
//! Only the struct and a plain serial `y = A x` live here (the orphan rule
//! means only this crate can implement `LinearOperator` for `CsrMatrix`).
//! Real operations (scaling, ILU, I/O) live in `elmer-sparse`.

use crate::linalg::{LinalgError, LinearOperator};
use crate::scalar::Scalar;

/// Elmer CRS (CRSMatrix.F90) in 0-based form. Columns sorted ascending per row
/// (CRS_SortMatrix, CRSMatrix.F90:175). `diag[i]` indexes A[i,i] in `values`.
#[derive(Clone, Debug)]
pub struct CsrMatrix<T: Scalar = f64> {
    pub nrows: usize,
    pub ncols: usize,
    /// Length `nrows + 1`.
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<u32>,
    pub values: Vec<T>,
    pub diag: Option<Vec<usize>>,
}

impl<T: Scalar> CsrMatrix<T> {
    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
}

impl<T: Scalar> LinearOperator<T> for CsrMatrix<T> {
    fn nrows(&self) -> usize {
        self.nrows
    }

    fn ncols(&self) -> usize {
        self.ncols
    }

    fn apply(&self, x: &[T], y: &mut [T]) {
        for (i, yi) in y.iter_mut().enumerate().take(self.nrows) {
            let mut s = T::zero();
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                s = s + self.values[k] * x[self.col_idx[k] as usize];
            }
            *yi = s;
        }
    }

    /// `y = A^T x` (plain transpose, no conjugation).
    fn apply_transpose(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() < self.nrows || y.len() < self.ncols {
            return Err(LinalgError::Dim(format!(
                "transpose of {}x{} needs x >= {}, y >= {}; got {}, {}",
                self.nrows,
                self.ncols,
                self.nrows,
                self.ncols,
                x.len(),
                y.len()
            )));
        }
        y[..self.ncols].fill(T::zero());
        for (i, &xi) in x.iter().enumerate().take(self.nrows) {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k] as usize;
                y[j] = y[j] + self.values[k] * xi;
            }
        }
        Ok(())
    }
}
