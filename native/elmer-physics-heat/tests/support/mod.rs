//! Test-only stand-ins for `elmer-fe` (basis, quadrature), `elmer-mesh` (mesh
//! file reader, structured meshes) and `elmer-assembly` (glue loop, small CG
//! solve). Replace with the real crates once S1 and L06 merge.
#![allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]

pub mod backends;
pub mod basis;
pub mod glue;
pub mod meshio;
pub mod quad;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
