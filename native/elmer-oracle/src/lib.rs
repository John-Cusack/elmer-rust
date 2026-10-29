//! External oracle for the Rust-native Elmer rewrite.
//!
//! This crate drives pristine Fortran Elmer builds directly (no ctest),
//! captures golden fields, and grades candidate output against them.
//! It must never depend on candidate crates; `tests/isolation.rs` enforces
//! that over the crate manifest and lockfile.

pub mod grade;
pub mod manifest;
pub mod norm;
pub mod result;
pub mod runner;

use std::fmt;

/// Errors produced anywhere in the oracle.
#[derive(Debug)]
pub enum OracleError {
    Io(std::io::Error, String),
    Parse(String),
    Unsupported(String),
    SolverFailed(String),
    GradeFailed(String),
}

impl fmt::Display for OracleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OracleError::Io(e, ctx) => write!(f, "io error in {ctx}: {e}"),
            OracleError::Parse(msg) => write!(f, "parse error: {msg}"),
            OracleError::Unsupported(msg) => write!(f, "unsupported test construct: {msg}"),
            OracleError::SolverFailed(msg) => write!(f, "solver run failed: {msg}"),
            OracleError::GradeFailed(msg) => write!(f, "grade failed: {msg}"),
        }
    }
}

impl std::error::Error for OracleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            OracleError::Io(e, _) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for OracleError {
    fn from(e: std::io::Error) -> Self {
        OracleError::Io(e, String::from("unspecified"))
    }
}

/// Convenience alias used across the crate.
pub type Result<T, E = OracleError> = std::result::Result<T, E>;

/// Attach file context to an I/O error.
pub(crate) fn io_with(path: &std::path::Path, e: std::io::Error) -> OracleError {
    OracleError::Io(e, path.display().to_string())
}
