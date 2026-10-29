//! Test-only stand-ins for `elmer-mesh` (mesh file reader, `MeshView`) and
//! `elmer-sif` (hand-built model sections), plus the pristine reference
//! readers. Replace with the real crates once S1 merges.
#![allow(dead_code)]

pub mod meshio;
pub mod run;
pub mod simmesh;
