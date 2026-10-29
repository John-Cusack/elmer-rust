#[derive(Debug, thiserror::Error)]
pub enum MeshLevelsError {
    /// A keyword or element the port refuses rather than approximates.
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// Input the Fortran stops on with `Fatal`.
    #[error("invalid: {0}")]
    Invalid(String),
}
