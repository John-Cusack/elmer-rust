//! Assembly with Elmer semantics: DOF permutation and bandwidth
//! optimisation, CSR pattern, element-to-global glue and Dirichlet
//! conditions.
//!
//! The usual sequence for one nodal solver is
//! [`create_system`] (Perm, bandwidth decision, zero-valued pattern),
//! [`assemble_elements`] running a physics kernel and feeding [`CsrAssembler`] through
//! [`GlobalAssembler`](elmer_core::GlobalAssembler), then
//! [`DirichletSet`] to record constraints and [`apply_dirichlet`] to
//! enforce them. Each step ports a named Elmer routine; see the module
//! docs for the source lines and for what is out of scope.
#![forbid(unsafe_code)]

pub mod bandwidth;
pub mod dirichlet;
pub mod driver;
pub mod glue;
pub mod list;
pub mod perm;

pub use bandwidth::{bandwidth, bandwidth_reordered, rcm_order, BandwidthReport};
pub use dirichlet::{apply_dirichlet, DirichletBc, DirichletMode, DirichletScaling, DirichletSet};
pub use driver::{assemble_elements, boundary_elements, bulk_elements};
pub use glue::CsrAssembler;
pub use list::ListMatrix;
pub use perm::{
    build_perm, build_perm_coupled, build_perm_with, create_matrix, create_system, initial_permutation,
    initialize_matrix, list_matrix, Perm, PermOptions, System,
};
