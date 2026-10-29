//! Sparse direct solvers on faer: `FaerLu` (the UMFPACK / MUMPS / SuperLU /
//! PARDISO slot) and `FaerCholesky` (the CHOLMOD slot).
//!
//! Both take Elmer's CSR matrix, convert it to the compressed-column form faer
//! wants, and keep the factorisation for repeated solves. `factor` may be
//! called again on a matrix with the same pattern; the symbolic analysis is
//! then reused (Elmer's `Factorize = False` re-solve path is just `solve_in_place`
//! on the kept factors).

use core::marker::PhantomData;

use elmer_core::{CsrMatrix, DirectSolver, LinalgError};
use faer::linalg::solvers::{Solve, SolveCore};
use faer::sparse::linalg::solvers::{Llt, Lu, SymbolicLlt, SymbolicLu};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::linalg::cholesky::llt::factor::LltError as NumericLltError;
use faer::sparse::linalg::{LltError, LuError};
use faer::{MatMut, Side};

use crate::scalar::DirectScalar;

/// Run faer's sparse factorisations and solves on the calling thread only.
/// faer reads a process-global parallelism setting that `RAYON_NUM_THREADS` and
/// `OMP_NUM_THREADS` do not control; call this once at start-up (tests do).
pub fn set_single_thread() {
    faer::set_global_parallelism(faer::Par::Seq);
}

/// Compressed-column copy of a CSR matrix (the CSC form of `A`, built by a
/// counting-sort transpose of the CSR arrays; rows ascend within each column).
struct Csc<T> {
    n: usize,
    col_ptr: Vec<usize>,
    row_idx: Vec<usize>,
    values: Vec<T>,
}

impl<T: DirectScalar> Csc<T> {
    fn from_csr(a: &CsrMatrix<T>) -> Result<Self, LinalgError> {
        if a.nrows != a.ncols {
            return Err(LinalgError::Dim(format!(
                "direct solve needs a square matrix, got {}x{}",
                a.nrows, a.ncols
            )));
        }
        if a.row_ptr.len() != a.nrows + 1
            || a.row_ptr.last().copied() != Some(a.col_idx.len())
            || a.col_idx.len() != a.values.len()
        {
            return Err(LinalgError::Dim("inconsistent CSR arrays".into()));
        }
        let n = a.nrows;
        let mut col_ptr = vec![0usize; n + 1];
        for &c in &a.col_idx {
            let c = c as usize;
            if c >= n {
                return Err(LinalgError::Dim(format!("column index {c} out of range for n = {n}")));
            }
            col_ptr[c + 1] += 1;
        }
        for c in 0..n {
            col_ptr[c + 1] += col_ptr[c];
        }
        let mut next = col_ptr.clone();
        let mut row_idx = vec![0usize; a.col_idx.len()];
        let mut values = vec![T::zero(); a.col_idx.len()];
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let c = a.col_idx[k] as usize;
                let dst = next[c];
                next[c] += 1;
                row_idx[dst] = i;
                values[dst] = a.values[k];
            }
        }
        Ok(Self { n, col_ptr, row_idx, values })
    }

    fn symbolic(&self) -> SymbolicSparseColMatRef<'_, usize> {
        SymbolicSparseColMatRef::new_checked(self.n, self.n, &self.col_ptr, None, &self.row_idx)
    }

    fn view(&self) -> SparseColMatRef<'_, usize, T> {
        SparseColMatRef::new(self.symbolic(), &self.values)
    }

    /// Same sparsity pattern as `other`? (Decides whether a symbolic analysis can be reused.)
    fn same_pattern(&self, other: &Self) -> bool {
        self.n == other.n && self.col_ptr == other.col_ptr && self.row_idx == other.row_idx
    }
}

fn check_rhs(n: usize, b_len: usize) -> Result<(), LinalgError> {
    if b_len == n {
        Ok(())
    } else {
        Err(LinalgError::Dim(format!("rhs length {b_len} does not match system size {n}")))
    }
}

/// Sparse LU with partial pivoting and a fill-reducing column ordering
/// (faer's supernodal-capable sparse LU). Replaces UMFPACK, MUMPS, SuperLU and
/// PARDISO for real and complex systems.
pub struct FaerLu<T: DirectScalar = f64> {
    state: Option<LuState<T>>,
    _t: PhantomData<T>,
}

struct LuState<T: DirectScalar> {
    pattern: Csc<T>,
    symbolic: SymbolicLu<usize>,
    lu: Lu<usize, T>,
}

impl<T: DirectScalar> Default for FaerLu<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: DirectScalar> FaerLu<T> {
    pub fn new() -> Self {
        Self { state: None, _t: PhantomData }
    }

    /// Order of the factored system, if `factor` has succeeded.
    pub fn size(&self) -> Option<usize> {
        self.state.as_ref().map(|s| s.pattern.n)
    }
}

fn map_lu_error(e: LuError) -> LinalgError {
    match e {
        LuError::Generic(e) => LinalgError::Breakdown {
            method: "faer sparse lu",
            detail: format!("{e:?}"),
        },
        LuError::SymbolicSingular { index } => LinalgError::Singular(index),
    }
}

fn map_llt_error(e: LltError) -> LinalgError {
    match e {
        LltError::Numeric(NumericLltError::NonPositivePivot { index }) => LinalgError::Breakdown {
            method: "faer sparse cholesky",
            detail: format!("non-positive pivot at row {index}"),
        },
        LltError::Generic(e) => LinalgError::Breakdown {
            method: "faer sparse cholesky",
            detail: format!("{e:?}"),
        },
    }
}

impl<T: DirectScalar> DirectSolver<T> for FaerLu<T> {
    fn factor(&mut self, a: &CsrMatrix<T>) -> Result<(), LinalgError> {
        let csc = Csc::from_csr(a)?;
        let symbolic = match self.state.take() {
            Some(old) if old.pattern.same_pattern(&csc) => old.symbolic,
            _ => SymbolicLu::try_new(csc.symbolic()).map_err(|e| LinalgError::Breakdown {
                method: "faer sparse lu",
                detail: format!("symbolic analysis failed: {e:?}"),
            })?,
        };
        let lu = Lu::try_new_with_symbolic(symbolic.clone(), csc.view()).map_err(map_lu_error)?;
        self.state = Some(LuState { pattern: csc, symbolic, lu });
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError> {
        let st = self
            .state
            .as_ref()
            .ok_or_else(|| LinalgError::Unsupported("solve before factor".into()))?;
        check_rhs(st.pattern.n, b.len())?;
        let n = b.len();
        st.lu.solve_in_place_with_conj(faer::Conj::No, MatMut::from_column_major_slice_mut(b, n, 1));
        Ok(())
    }
}

/// Sparse Cholesky `A = L L^H` for symmetric (Hermitian) positive definite
/// systems. Replaces CHOLMOD. Only the lower triangle of the CSR input is read.
pub struct FaerCholesky<T: DirectScalar = f64> {
    state: Option<LltState<T>>,
}

struct LltState<T: DirectScalar> {
    pattern: Csc<T>,
    symbolic: SymbolicLlt<usize>,
    llt: Llt<usize, T>,
}

impl<T: DirectScalar> Default for FaerCholesky<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: DirectScalar> FaerCholesky<T> {
    pub fn new() -> Self {
        Self { state: None }
    }
}

impl<T: DirectScalar> DirectSolver<T> for FaerCholesky<T> {
    fn factor(&mut self, a: &CsrMatrix<T>) -> Result<(), LinalgError> {
        let csc = Csc::from_csr(a)?;
        let symbolic = match self.state.take() {
            Some(old) if old.pattern.same_pattern(&csc) => old.symbolic,
            _ => SymbolicLlt::try_new(csc.symbolic(), Side::Lower).map_err(|e| {
                LinalgError::Breakdown {
                    method: "faer sparse cholesky",
                    detail: format!("symbolic analysis failed: {e:?}"),
                }
            })?,
        };
        let llt = Llt::try_new_with_symbolic(symbolic.clone(), csc.view(), Side::Lower)
            .map_err(map_llt_error)?;
        self.state = Some(LltState { pattern: csc, symbolic, llt });
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError> {
        let st = self
            .state
            .as_ref()
            .ok_or_else(|| LinalgError::Unsupported("solve before factor".into()))?;
        check_rhs(st.pattern.n, b.len())?;
        let n = b.len();
        st.llt.solve_in_place(MatMut::from_column_major_slice_mut(b, n, 1));
        Ok(())
    }
}
