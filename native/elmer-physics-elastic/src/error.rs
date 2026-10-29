/// Failure of a local kernel.
#[derive(Debug, thiserror::Error)]
pub enum StressError {
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("bad input: {0}")]
    Input(String),
    #[error("basis evaluation failed: {0}")]
    Basis(String),
}
