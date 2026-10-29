use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum FeError {
    #[error("element type {0} is not defined in elements.def")]
    UnknownType(u16),
    #[error("element type {0} is not a Lagrange (nodal) element supported by elmer-fe")]
    Unsupported(u16),
    #[error("invalid number of integration points {np} for element family {family}")]
    InvalidPoints { family: u16, np: usize },
    #[error("relative integration order {0} outside the -1..=1 rules of a non-p element")]
    InvalidRelOrder(i32),
    #[error("expected {expected} nodal coordinates for element type {etype}, got {got}")]
    NodeCount {
        etype: u16,
        expected: usize,
        got: usize,
    },
    #[error("degenerate {dim}D element (DetG = {detg:e})")]
    Degenerate { dim: usize, detg: f64 },
    #[error("element basis matrix of type {0} is singular")]
    SingularBasis(u16),
}
