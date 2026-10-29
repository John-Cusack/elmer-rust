//! Lagrange element library: Elmer's reference elements, Gauss rules,
//! nodal bases, element metric and stabilization parameters.
//!
//! A structural transliteration of Integration.F90, ElementBasis.F90 /
//! PElementBase.F90 (nodal bases), ElemInfo.F90 (`ElementInfo`,
//! `ElementMetric`, `StabParam`, with the reference-LAPACK `DSYGV` path it
//! calls) and CoordinateSystems.F90. Expression order
//! follows the Fortran so values match the pristine binary bit for bit; the
//! reference tests compare against values captured from it. The
//! p-hierarchic and edge bases live in `elmer-fe-hp`.
#![forbid(unsafe_code)]

pub mod basis;
pub mod coords;
mod dd;
pub mod element_type;
pub mod error;
mod gauss1d;
pub mod info;
mod lapack;
mod lu;
pub mod metric;
pub mod quadrature;
pub mod stab;
mod tables;

pub use coords::{CoordinateInfo, CoordinateKind, CoordinateSystem};
pub use element_type::{element_type, element_types, BasisTerm, ElementTypeDef};
pub use error::FeError;
pub use info::{element_info, LagrangeBasis};
pub use metric::MetricStatus;
pub use quadrature::{gauss_points, gauss_points_rel, ElmerQuadrature, GaussRule};
pub use stab::{element_diameter, stab_param, stab_param_mk_only, type_stabilization_mk};
