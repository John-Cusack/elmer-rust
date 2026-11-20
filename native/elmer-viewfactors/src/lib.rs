//! View factors and Gebhart factors of diffuse gray radiation for the
//! Rust-native Elmer: the numerics of the `ViewFactors` program
//! (`fem/src/ViewFactors.F90`, `fem/src/viewaxis`) and of
//! `RadiationFactors.F90`, without the files in between.
//!
//! * [`view_factors_axis`]: factors of axisymmetric boundaries;
//! * [`normalize_factors`]: the normalisation `ViewFactors` applies before it
//!   writes them, and [`saved_rows`] what it keeps;
//! * [`gebhart_factors`]: `RadiationFactors`' Gebhart factors from them.
//!
//! The mesh side (which boundary elements radiate, their orientation, areas
//! and emissivities) belongs to the caller.
#![forbid(unsafe_code)]

mod axisym;
mod gebhart;
mod normalize;
mod orient;
mod radiator;

pub use axisym::view_factors_axis;
pub use gebhart::{gebhart_factors, gebhart_factors_with, Gebhart, GebhartInput, GebhartMatrix, GebhartSolve};
pub use normalize::{normalize_factors, normalize_factors_with, saved_rows, IterSolv};
pub use orient::{line_direction, surface_direction, ParentInfo};
pub use radiator::{radiator_factors, surface_areas, Radiator, RadiatorParams, RadiatorSurface};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    Invalid(String),
    #[error("the radiation matrix is singular at row {0}")]
    Singular(usize),
}
