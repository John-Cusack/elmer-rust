//! Test-only stand-ins for `elmer-mesh` (mesh file reader, `MeshView`) and
//! `elmer-sif` (hand-built model sections) until those merge, plus a dense LU
//! standing in for UMFPACK. The element basis, Gauss rules and `StabParam` are
//! the real `elmer-fe`.
#![allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]

pub mod backends;
pub mod cases;
pub mod meshio;
