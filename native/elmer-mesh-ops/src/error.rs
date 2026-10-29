#[derive(Debug, thiserror::Error)]
pub enum MeshOpsError {
    /// A keyword or element the port refuses rather than approximates.
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// Input the Fortran stops on with `Fatal`.
    #[error("invalid: {0}")]
    Invalid(String),
}

impl From<elmer_mesh_levels::MeshLevelsError> for MeshOpsError {
    fn from(e: elmer_mesh_levels::MeshLevelsError) -> Self {
        match e {
            elmer_mesh_levels::MeshLevelsError::Unsupported(m) => MeshOpsError::Unsupported(m),
            elmer_mesh_levels::MeshLevelsError::Invalid(m) => MeshOpsError::Invalid(m),
        }
    }
}
