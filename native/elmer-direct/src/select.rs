//! `Linear System Direct Method` dispatch (`DirectSolve.F90:4266`).

use elmer_core::{DirectSolver, LinalgError, Params};

use crate::band::{BandCholesky, BandLu};
use crate::scalar::DirectScalar;
use crate::sparse::{FaerCholesky, FaerLu};

/// The keyword Elmer reads.
pub const METHOD_KEY: &str = "Linear System Direct Method";

/// Which implementation a direct-method keyword value maps to.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DirectKind {
    /// `banded`: LAPACK band LU.
    BandLu,
    /// `symmetric banded`: LAPACK band Cholesky (lower triangle).
    BandCholesky,
    /// `umfpack`, `big umfpack`, `mumps` family, `superlu`, `pardiso`, `cpardiso`:
    /// general sparse LU.
    SparseLu,
    /// `cholmod`: sparse Cholesky for SPD systems.
    SparseCholesky,
}

/// Map a `Linear System Direct Method` value to an implementation. Elmer's
/// default (keyword absent) is `banded`. Methods with no native counterpart
/// (`spqr`, `permon`, ...) are refused rather than silently served by another
/// method; an explicit band method never falls back to sparse LU.
pub fn kind_for_method(method: &str) -> Result<DirectKind, LinalgError> {
    match method.trim().to_ascii_lowercase().as_str() {
        "banded" => Ok(DirectKind::BandLu),
        "symmetric banded" => Ok(DirectKind::BandCholesky),
        "umfpack" | "big umfpack" | "smumps" | "cmumps" | "mumps" | "dmumps" | "zmumps"
        | "mumpslocal" | "superlu" | "pardiso" | "cpardiso" => Ok(DirectKind::SparseLu),
        "cholmod" => Ok(DirectKind::SparseCholesky),
        _ => Err(LinalgError::Unsupported(format!("unknown direct solver method '{method}'"))),
    }
}

/// Build the direct solver a solver section asks for, for either scalar field.
pub fn try_select_direct_for<T: DirectScalar>(
    params: &dyn Params,
) -> Result<Box<dyn DirectSolver<T>>, LinalgError> {
    let kind = match params.string(METHOD_KEY) {
        Some(m) => kind_for_method(m)?,
        None => DirectKind::BandLu,
    };
    Ok(match kind {
        DirectKind::BandLu => Box::new(BandLu::<T>::new()),
        DirectKind::BandCholesky => Box::new(BandCholesky::<T>::new()),
        DirectKind::SparseLu => Box::new(FaerLu::<T>::new()),
        DirectKind::SparseCholesky => Box::new(FaerCholesky::<T>::new()),
    })
}

/// Real-valued [`try_select_direct_for`].
pub fn try_select_direct(params: &dyn Params) -> Result<Box<dyn DirectSolver<f64>>, LinalgError> {
    try_select_direct_for::<f64>(params)
}

/// The real direct solver for a solver section. Panics on an unknown method,
/// as Elmer's `Fatal( 'DirectSolver', 'Unknown direct solver method.' )` stops
/// the run; use [`try_select_direct`] to handle that case.
pub fn select_direct(params: &dyn Params) -> Box<dyn DirectSolver<f64>> {
    try_select_direct(params).unwrap_or_else(|e| panic!("DirectSolver: {e}"))
}
