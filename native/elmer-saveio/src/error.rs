//! Error type shared by every reader and writer in the crate.

#[derive(Debug, thiserror::Error)]
pub enum SaveIoError {
    #[error("i/o error: {0}")]
    Io(#[from] std::io::Error),
    /// The bytes are not a well-formed file of the expected kind.
    #[error("malformed file: {0}")]
    Format(String),
    /// Well-formed but outside what this crate supports (old file versions,
    /// complex fields, ...). Never silently approximated.
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// The inputs handed to a writer are inconsistent.
    #[error("invalid input: {0}")]
    Input(String),
}

pub type Result<T> = std::result::Result<T, SaveIoError>;

pub(crate) fn format_err<T>(msg: impl Into<String>) -> Result<T> {
    Err(SaveIoError::Format(msg.into()))
}
