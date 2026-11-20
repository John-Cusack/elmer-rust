//! `LocalBulkMatrix` against hand-computed integrals on a unit line element.

use elmer_core::{ElementType, LocalSystem, Quadrature};
use elmer_physics_filmflow::kernel::{
    local_bulk_matrix, ns_condensate, BulkInput, BulkOptions, Friction,
};

fn opts() -> BulkOptions {
    BulkOptions {
        mdim: 1,
        csymmetry: false,
        friction: Friction::Laminar,
        roughness: 0.0,
        gravity: None,
        convect: false,
        newton: false,
        lateral_strain: false,
        gradp: false,
        transient: false,
        min_speed: 1e-6,
        calc_friction_heating: false,
        melt_heat: None,
        first_round: false,
    }
}

/// With unit gap, density and viscosity, the momentum row is the
/// Poiseuille resistance `12 mu / gap^2` times the consistent mass and the
/// pressure gradient, the continuity row `d(gap rho u)/dx`.
#[test]
fn laminar_line_element() {
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let one = vec![1.0; 2];
    let zero = vec![vec![0.0; 2]; 5];
    let velo = vec![vec![0.0; 2]; 2];
    let inp = BulkInput {
        etype: ElementType(202),
        xyz: &xyz,
        ntot: 2,
        rho: &one,
        mu: &one,
        gap: &one,
        height: None,
        pres: &[0.0, 0.0],
        velo: &velo,
        load: &zero,
    };
    let basis = elmer_fe::LagrangeBasis::new(1);
    let rule = elmer_fe::ElmerQuadrature.gauss_points(ElementType(202), None);
    let mut out = LocalSystem::zeros(0);
    local_bulk_matrix(&basis, &rule, &inp, &opts(), &mut out).unwrap();
    let k = |r: usize, c: usize| out.stiff[r * 4 + c];
    let close = |a: f64, b: f64| assert!((a - b).abs() < 1e-14, "{a} vs {b}");
    // u-u: 12 * [1/3 1/6; 1/6 1/3].
    close(k(0, 0), 4.0);
    close(k(0, 2), 2.0);
    // u-p: -int phi_q dphi_p/dx; p-u: int gap rho dphi_q/dx phi_p.
    close(k(0, 1), 0.5);
    close(k(0, 3), 0.5);
    close(k(1, 0), -0.5);
    close(k(1, 2), 0.5);
    close(k(1, 1), 0.0);
}

/// The condensed system keeps the retained rows' solution of the full one.
#[test]
fn condensation_eliminates_the_bubble_velocity() {
    // Two functions (one retained, one bubble), two dofs each: the bubble's
    // pressure row is already pinned.
    let mut l = LocalSystem {
        n: 4,
        stiff: vec![
            4.0, 1.0, 1.0, 0.0, //
            1.0, 0.0, 0.5, 0.0, //
            1.0, 0.5, 2.0, 0.0, //
            0.0, 0.0, 0.0, 1.0,
        ],
        force: vec![1.0, 0.0, 2.0, 0.0],
        mass: None,
        damp: None,
    };
    ns_condensate(&mut l, 1, 1, 1).unwrap();
    assert_eq!(l.n, 2);
    // K_cc - K_cb K_bb^-1 K_bc with K_bb = 2.
    assert_eq!(
        l.stiff,
        vec![4.0 - 0.5, 1.0 - 0.25, 1.0 - 0.25, 0.0 - 0.125]
    );
    assert_eq!(l.force, vec![1.0 - 1.0, 0.0 - 0.5]);
}
