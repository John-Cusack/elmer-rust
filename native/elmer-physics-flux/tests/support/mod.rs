//! Test-only stand-ins for `elmer-fe` (basis, quadrature) and `elmer-mesh`
//! (mesh file reader, `MeshView`). Replace with the real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]

pub mod basis;
pub mod meshio;
pub mod quad;
pub mod simmesh;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
