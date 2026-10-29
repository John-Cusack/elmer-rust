/// Why a kernel call or a preflight check was refused.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum AdvDiffError {
    /// A nodal input array does not have one entry per element node.
    #[error("{what}: expected {expected} nodal values, got {got}")]
    NodalLength {
        what: &'static str,
        expected: usize,
        got: usize,
    },
    #[error("element basis failed: {0}")]
    Basis(String),
    #[error("quadrature returned no points for element type {0}")]
    NoQuadrature(u16),
    /// The case uses a feature this crate does not implement.
    #[error("unsupported {section} keyword \"{keyword}\": {reason}")]
    Unsupported {
        section: &'static str,
        keyword: String,
        reason: &'static str,
    },
    /// Elmer `CALL Fatal(...)` of the Fortran solver.
    #[error("{0}")]
    Fatal(String),
    #[error("invalid input: {0}")]
    Invalid(String),
}

pub(crate) fn check_len(
    what: &'static str,
    expected: usize,
    got: usize,
) -> Result<(), AdvDiffError> {
    if expected == got {
        Ok(())
    } else {
        Err(AdvDiffError::NodalLength {
            what,
            expected,
            got,
        })
    }
}

/// `sum_i basis[i] * nodal[i]`, the interpolation `ListGetElementReal` does.
pub(crate) fn interp(basis: &[f64], nodal: &[f64]) -> f64 {
    basis.iter().zip(nodal).map(|(b, v)| b * v).sum()
}
