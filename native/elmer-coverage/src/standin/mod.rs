//! Test-local stand-ins for S1's `elmer-mesh` (mesh reader) and `elmer-fe`
//! (Lagrange bases and Gauss rules, copied from `elmer-run`'s tests). Replace
//! them with the real crates once S1 merges; every refusal they cause is
//! labelled `stand-in` in the report.
#![allow(clippy::needless_range_loop)]

pub mod basis;
pub mod mesh;
pub mod quad;

use elmer_physics_heat::FeProvider;

/// Element types `TestBasis` evaluates.
pub const SUPPORTED_ELEMENTS: [u16; 7] = [202, 203, 303, 306, 404, 408, 504];

pub fn fe() -> FeProvider {
    || (Box::new(basis::TestBasis), Box::new(quad::TestQuadrature))
}
