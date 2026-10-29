//! Errors for the time-integration lane. Most are `Fatal` calls in the Fortran.

#[derive(Debug, thiserror::Error)]
pub enum TimeError {
    /// `BDFLocal` / `BDF_CRS` `CASE DEFAULT`, and the `BDF Order` range check (1..5).
    #[error("Invalid order BDF {0}")]
    InvalidBdfOrder(i64),
    #[error("unknown timestepping method: {0}")]
    UnknownMethod(String),
    #[error("{0}: mass matrix is required")]
    MissingMass(&'static str),
    #[error("{0}: lumped mass values are required")]
    MissingLumpedMass(&'static str),
    #[error("{caller}: needs {needed} previous solution vector(s), got {got}")]
    MissingHistory {
        caller: &'static str,
        needed: usize,
        got: usize,
    },
    #[error("{0}: current solution is required")]
    MissingCurrent(&'static str),
    #[error("Adams-Moulton method must be executed after Adams-Bashforth method!")]
    AdamsMoultonBeforeBashforth,
    #[error("{0}: not implemented for this scheme")]
    NotImplemented(&'static str),
    #[error("fractional step: invalid fsstep {0}")]
    InvalidFsStep(f64),
    #[error("dimension mismatch: {0}")]
    Dim(String),
    #[error("{0}")]
    Config(String),
    #[error("Timestep smaller than epsilon: {0:e}")]
    TimestepTooSmall(f64),
    #[error("Could not find stable timestep above given minimum")]
    NoStableTimestep,
}
