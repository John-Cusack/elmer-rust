//! Cross-lane types and traits for the Rust-native Elmer rewrite.
//!
//! Only types that cross a crate (lane) boundary live here; each crate keeps
//! its own local types. Owners and consumers of every item are listed in
//! `native/docs/interfaces.md`, and changes follow the `core-change` rule there.
#![forbid(unsafe_code)]

pub mod comm;
pub mod csr;
pub mod fe;
pub mod field;
pub mod linalg;
pub mod local;
pub mod log;
pub mod mesh;
pub mod params;
pub mod scalar;

pub use comm::{Communicator, SerialComm};
pub use csr::CsrMatrix;
pub use fe::{BasisEval, ElementBasis, QuadPoint, Quadrature};
pub use field::{FieldKind, FieldRef, FieldSource, NO_DOF};
pub use linalg::{
    DirectSolver, IdentityPrecond, IterativeSolver, LinalgError, LinearOperator, Preconditioner,
    SolveStats, SolveStatus,
};
pub use local::{GlobalAssembler, LocalSystem};
pub use log::{Logger, NullLogger};
pub use mesh::{ElemId, ElementRef, ElementType, MeshView};
pub use params::{normalize_key, Dependent, DependentKind, ParamMap, Params, Value};
pub use scalar::Scalar;
