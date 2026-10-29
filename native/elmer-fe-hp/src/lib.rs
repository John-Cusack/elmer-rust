//! p-hierarchic H1 bases and Whitney/Nédélec edge bases (lane L15).
//!
//! A structural transliteration of Elmer's `H1Basis.F90`, `PElementBase.F90`,
//! `PElementMaps.F90` and the edge-element parts of `ElemInfo.F90` on the
//! reference elements, with no mesh or solver dependencies.
#![forbid(unsafe_code)]

pub mod brick;
pub mod edge;
pub mod legendre;
pub mod line;
pub mod maps;
pub mod pbasis;
pub mod pyramid;
pub mod quad;
pub mod tetra;
pub mod triangle;
pub mod types;
pub mod wedge;

pub use pbasis::{h1_p_basis, h1_p_basis_spec, p_dof_counts, PElementBasis, PElementSpec};
pub use types::{pow_int, EdgeDir, FaceDir, MAX_EDGES, MAX_FACES};
