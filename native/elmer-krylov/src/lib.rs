//! HUTI Krylov solvers for the Rust-native Elmer: CG, CGS, BiCGStab,
//! BiCGStab(2), GMRES(m), QMR and TFQMR, real and complex, with HUTI's
//! stopping rules and stdout lines, plus the IterSolve keyword mapping that
//! configures them (fhutiter/src, fem/src/IterSolve.F90).
//!
//! ```
//! use elmer_core::{IterativeSolver, IdentityPrecond, ParamMap, SerialComm, Value};
//! use elmer_krylov::{Huti, HutiConfig};
//!
//! let mut p = ParamMap::new();
//! p.insert("Linear System Iterative Method", Value::Str("CG".into()));
//! p.insert("Linear System Max Iterations", Value::Integer(vec![50]));
//! p.insert("Linear System Convergence Tolerance", Value::real(1e-10));
//! p.insert("Linear System Residual Output", Value::Integer(vec![0]));
//! let mut solver: Huti = Huti::new(HutiConfig::from_params(&p).unwrap());
//!
//! // 1D Laplacian, 3 unknowns.
//! let a = elmer_core::CsrMatrix {
//!     nrows: 3, ncols: 3,
//!     row_ptr: vec![0, 2, 5, 7],
//!     col_idx: vec![0, 1, 0, 1, 2, 1, 2],
//!     values: vec![2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0],
//!     diag: None,
//! };
//! let b = [1.0, 0.0, 1.0];
//! let mut x = [0.0; 3];
//! let stats = solver.solve(&a, &IdentityPrecond, &b, &mut x, &SerialComm).unwrap();
//! assert_eq!(stats.status, elmer_core::SolveStatus::Converged);
//! assert!((x[1] - 1.0).abs() < 1e-12);
//! ```
#![forbid(unsafe_code)]

mod config;
mod huti;
mod methods;
mod output;
mod scalar;

pub use config::{HutiConfig, HutiMethod, HutiStop};
pub use huti::{Huti, HutiInfo, OutputSink};
pub use output::{format_e, format_i, iteration_line};
pub use scalar::HutiScalar;
