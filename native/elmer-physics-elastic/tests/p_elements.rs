//! `StressSolver` on p elements through `elmer_solver::Simulation` with
//! `elmer-pdofs`' `PDofsBackend`: `fem/tests/ElastPelem2dQuad` (`Element =
//! p:4`, full elements) with only its StressSolve solver (Calculate Stresses
//! removed: it does not change the displacement) and `ElastPelem2dTri`
//! (`p:3` triangles), against the pristine runs' Perm, every slot's
//! displacement and the norms 5.99181622E-03 / 7.71519407E-03
//! (`tests/data/pelem2d_quad`, `tests/data/pelem2d_tri`).

#[path = "support/pcase.rs"]
mod pcase;

use elmer_physics_elastic::{install_fe, FeServices};

/// [`pcase::check_case`] in 2D.
fn check_case(case: &str, n_slots: usize, norm_want: f64, reference: f64) {
    let _ = install_fe(FeServices {
        basis: Box::new(elmer_fe::LagrangeBasis::new(2)),
        quad: Box::new(elmer_fe::ElmerQuadrature),
    });
    pcase::check_case(case, 2, n_slots, norm_want, reference);
}

#[test]
fn elast_pelem_2d_quad_matches_pristine() {
    // p:4 on 12 quads, full elements: 225 slots.
    check_case("pelem2d_quad", 225, 0.005991816215121997, 5.99181622e-3);
}

#[test]
fn elast_pelem_2d_tri_matches_pristine() {
    // p:3 on triangles (edge functions of degree 2 and 3, one bubble each).
    check_case("pelem2d_tri", 133, 0.007715194063294897, 7.71519407e-3);
}
