use elmer_core::LinalgError;

/// Errors from the solver framework and from `Solver` implementations.
#[derive(Debug, thiserror::Error)]
pub enum SolverError {
    /// Elmer `CALL Fatal(caller, msg)`.
    #[error("{caller}: {msg}")]
    Fatal { caller: String, msg: String },
    /// Elmer `CALL NumericalError(caller, msg)` (NaN norm, norm over the max).
    #[error("{caller}: {msg}")]
    Numerical { caller: String, msg: String },
    #[error(transparent)]
    Linalg(#[from] LinalgError),
    #[error("variable not found: {0}")]
    NoVariable(String),
    /// A pluggable backend (assembly, linear solver, ...) the call needs has
    /// not been installed in the `Simulation`.
    #[error("backend not installed: {0}")]
    NotWired(String),
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error(transparent)]
    Expr(#[from] elmer_expr::ExprError),
    #[error(transparent)]
    Time(#[from] elmer_timeint::TimeError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

impl SolverError {
    pub fn fatal(caller: &str, msg: impl Into<String>) -> Self {
        SolverError::Fatal {
            caller: caller.into(),
            msg: msg.into(),
        }
    }

    pub fn numerical(caller: &str, msg: impl Into<String>) -> Self {
        SolverError::Numerical {
            caller: caller.into(),
            msg: msg.into(),
        }
    }
}
