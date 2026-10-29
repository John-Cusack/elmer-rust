//! Test-only stand-ins for `elmer-fe` (basis, quadrature), `elmer-mesh` (mesh
//! file reader, `MeshView`) and `elmer-sif` (hand-built model sections).
//! Replace with the real crates once S1 merges.
#![allow(dead_code)]

pub mod basis;
pub mod cases;
pub mod meshio;
pub mod quad;
pub mod simmesh;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
