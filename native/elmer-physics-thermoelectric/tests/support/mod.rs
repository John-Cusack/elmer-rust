//! Test-only stand-ins for `elmer-mesh` (mesh-file reader, `MeshView` with
//! boundary parents) and `elmer-sif` (hand-built model sections). The basis
//! and Gauss rules are `elmer-fe`'s, the solver's defaults.
#![allow(dead_code, clippy::needless_range_loop)]

pub mod cases;
pub mod mesh;
