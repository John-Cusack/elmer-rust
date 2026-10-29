/// Why a kernel call was refused.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum IncompressibleError {
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
    /// The stabilization length scale is zero, so the stabilization would
    /// silently vanish (`LocalBulkMatrix` makes this fatal).
    #[error("Stabilization element size is zero, so the stabilisation would silently vanish")]
    ZeroElementSize,
    /// The element type or option is outside what this crate ports.
    #[error("unsupported: {0}")]
    Unsupported(String),
}
