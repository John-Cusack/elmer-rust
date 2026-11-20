//! `UmfpackLu`: Elmer's `UMFPack_SolveSystem` (DirectSolve.F90) on the port,
//! as an [`elmer_core::DirectSolver`].
//!
//! Elmer hands UMFPACK its CRS arrays as if they were compressed columns
//! (`umf4sym (n, n, Rows-1, Cols-1, Values, ...)`), so UMFPACK factors `A'`,
//! and then solves with `sys = 2` (`UMFPACK_Aat`, i.e. `A.' x = b` of that
//! factor, which is `A x = b`) without iterative refinement (`umf4sol` sets
//! `Control [UMFPACK_IRSTEP] = 0`), all with `umf4def` default Control.

use std::cell::RefCell;

use elmer_core::{CsrMatrix, DirectSolver, LinalgError};

use crate::umf::{
    self, Numeric, UMFPACK_Aat, UMFPACK_ERROR_invalid_matrix, UMFPACK_CONTROL, UMFPACK_INFO,
};

/// True for the `Linear System Direct Method` values Elmer serves with
/// UMFPACK (`umfpack`, `big umfpack`; the long-index build gives the same
/// factorisation for systems that fit `int`).
pub fn is_umfpack_method(method: &str) -> bool {
    matches!(
        method.trim().to_ascii_lowercase().as_str(),
        "umfpack" | "big umfpack"
    )
}

/// What [`UmfpackLu::factor`] did.
enum Factor {
    None,
    Umfpack(Numeric),
    /// The CRS arrays are not a valid UMFPACK matrix (unsorted or duplicate
    /// column indices), where Elmer's own CRS always is; served by faer's LU
    /// instead of failing.
    Fallback(Box<dyn DirectSolver<f64> + Send>),
}

/// UMFPACK LU of a real CRS matrix, bitwise as Elmer's bundled UMFPACK 4.4.
pub struct UmfpackLu {
    n: usize,
    factor: RefCell<Factor>,
    fallback: fn() -> Box<dyn DirectSolver<f64> + Send>,
}

impl UmfpackLu {
    /// `fallback` builds the solver used for matrices UMFPACK rejects as
    /// invalid (see [`Factor::Fallback`]).
    pub fn new(fallback: fn() -> Box<dyn DirectSolver<f64> + Send>) -> Self {
        UmfpackLu {
            n: 0,
            factor: RefCell::new(Factor::None),
            fallback,
        }
    }
}

fn err(stage: &str, status: i32) -> LinalgError {
    LinalgError::Breakdown {
        method: "UMFPACK",
        detail: format!("Error occurred in {stage}: {status}"),
    }
}

impl DirectSolver<f64> for UmfpackLu {
    fn factor(&mut self, a: &CsrMatrix<f64>) -> Result<(), LinalgError> {
        if a.nrows != a.ncols {
            return Err(LinalgError::Dim(format!(
                "direct solve needs a square matrix, got {}x{}",
                a.nrows, a.ncols
            )));
        }
        let n = a.nrows;
        self.n = n;
        if n == 0 {
            *self.factor.get_mut() = Factor::None;
            return Ok(());
        }
        let nnz = a.col_idx.len();
        if n > i32::MAX as usize || nnz > i32::MAX as usize {
            return Err(LinalgError::Unsupported(
                "UMFPACK: matrix too large for int indices".into(),
            ));
        }
        let ap: Vec<i32> = a.row_ptr.iter().map(|&p| p as i32).collect();
        let ai: Vec<i32> = a.col_idx.iter().map(|&c| c as i32).collect();
        let mut control = [0.0; UMFPACK_CONTROL];
        umf::umfpack_defaults(&mut control);
        let mut info = [0.0; UMFPACK_INFO];
        let ni = n as i32;
        let sym = match umf::umfpack_symbolic(
            ni,
            ni,
            &ap,
            &ai,
            Some(&a.values),
            Some(&control),
            Some(&mut info),
        ) {
            Ok(s) => s,
            Err(st) if st == UMFPACK_ERROR_invalid_matrix => {
                let mut f = (self.fallback)();
                f.factor(a)?;
                *self.factor.get_mut() = Factor::Fallback(f);
                return Ok(());
            }
            Err(st) => return Err(err("umf4sym", st)),
        };
        // a singular-matrix warning (status 1) is not an error in Elmer
        let (num, _status) =
            umf::umfpack_numeric(&ap, &ai, &a.values, &sym, Some(&control), Some(&mut info))
                .map_err(|st| err("umf4num", st))?;
        *self.factor.get_mut() = Factor::Umfpack(num);
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), LinalgError> {
        if b.len() != self.n {
            return Err(LinalgError::Dim(format!(
                "rhs has length {}, matrix order {}",
                b.len(),
                self.n
            )));
        }
        match &mut *self.factor.borrow_mut() {
            Factor::None => {
                if self.n == 0 {
                    Ok(())
                } else {
                    Err(LinalgError::Unsupported(
                        "UMFPACK: solve before factor".into(),
                    ))
                }
            }
            Factor::Fallback(f) => f.solve_in_place(b),
            Factor::Umfpack(num) => {
                let mut x = vec![0.0; self.n];
                let mut info = [0.0; UMFPACK_INFO];
                let st = umf::umfpack_solve(UMFPACK_Aat, &mut x, b, num, Some(&mut info));
                if st < 0 {
                    return Err(err("umf4sol", st));
                }
                b.copy_from_slice(&x);
                Ok(())
            }
        }
    }
}
