//! `StressSolver` on p:2 tetrahedra through `elmer_solver::Simulation` with
//! `elmer-pdofs`' `PDofsBackend` (a lane-made case on the StokesStabilized
//! tet mesh, `tests/data/pelem3d_tet/case.sif`): p tetrahedra take their
//! nodes in the `ConvertToACTetra` order of their dofs. A test binary of its
//! own, as the solver's basis is installed once per process.

#[path = "support/pcase.rs"]
mod pcase;

use elmer_physics_elastic::{install_fe, FeServices};

#[test]
fn elast_p2_tetrahedra_match_pristine() {
    let _ = install_fe(FeServices {
        basis: Box::new(elmer_fe::LagrangeBasis::new(3)),
        quad: Box::new(elmer_fe::ElmerQuadrature),
    });
    // 143 nodes + 660 edges; the pristine run printed 0.74207688E-03.
    pcase::check_case("pelem3d_tet", 3, 803, 0.0007420768781849677, 7.4207688e-4);
}
