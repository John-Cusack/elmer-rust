/// Why a kernel call or a preflight check was refused.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum HeatError {
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
    /// The case uses a HeatSolve feature this crate does not implement.
    #[error("unsupported {section} keyword \"{keyword}\": {reason}")]
    Unsupported {
        section: &'static str,
        keyword: String,
        reason: &'static str,
    },
    #[error("invalid input: {0}")]
    Invalid(String),
}
