//! Test-only stand-ins for `elmer-fe` (basis, quadrature), `elmer-mesh` (mesh
//! file reader, structured meshes, `MeshView`) and `elmer-sif` (hand-built
//! model sections), plus the kernel-level glue loop the matrix comparisons
//! use. Replace with the real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]

pub mod backends;
pub mod basis;
pub mod cases;
pub mod glue;
pub mod meshio;
pub mod quad;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
