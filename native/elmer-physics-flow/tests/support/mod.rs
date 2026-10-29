//! Test-only stand-ins for `elmer-fe` (basis, quadrature), `elmer-mesh`
//! (mesh-file reader, `MeshView` with boundary parents) and `elmer-sif`
//! (hand-built model sections). Replace with the real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop)]

pub mod basis;
pub mod cases;
pub mod mesh;
pub mod quad;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
