//! S1 adapter check (swarm plan 4.4 item 6), landed commented out by L00.
//!
//! When S1 merges `elmer-sif`, `elmer-mesh` and `elmer-fe`, firstmate enables
//! this file by uncommenting the body below and adding to
//! `native/elmer-core/Cargo.toml`:
//!
//! ```toml
//! [dev-dependencies]
//! elmer-sif = { path = "../elmer-sif" }
//! elmer-mesh = { path = "../elmer-mesh" }
//! elmer-fe = { path = "../elmer-fe" }
//! ```
//!
//! then runs `cargo check --tests -p elmer-core`. Replace the placeholder type
//! names with S1's real ones. A missing impl is one S1 follow-up PR, not a
//! change in every consumer.

// fn assert_mesh_view<T: elmer_core::MeshView>() {}
// fn assert_quadrature<T: elmer_core::Quadrature>() {}
// fn assert_element_basis<T: elmer_core::ElementBasis>() {}
// fn assert_params<T: elmer_core::Params>() {}
//
// #[test]
// fn s1_crates_implement_core_views() {
//     assert_mesh_view::<elmer_mesh::Mesh>();
//     assert_quadrature::<elmer_fe::GaussQuadrature>();
//     assert_element_basis::<elmer_fe::LagrangeBasis>();
//     // Either the section type implements Params directly ...
//     assert_params::<elmer_sif::Section>();
//     // ... or converts losslessly into a ParamMap:
//     // let _: fn(&elmer_sif::Section) -> elmer_core::ParamMap = |s| s.into();
// }
