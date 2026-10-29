//! Test-only stand-ins for `elmer-fe` (basis, quadrature), `elmer-mesh` (mesh
//! file reader) and `elmer-assembly` / `elmer-timeint` / a direct solver,
//! copied from `elmer-physics-heat`'s tests. Replace with the real crates
//! once S1 merges.
#![allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]

pub mod backends;
pub mod basis;
pub mod meshio;
pub mod quad;

pub use basis::TestBasis;
pub use quad::TestQuadrature;
pub mod harness;
