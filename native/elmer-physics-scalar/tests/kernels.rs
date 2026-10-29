//! Element kernels on single elements, `DenseLu`, the complex-to-real layout
//! and the registry names, without a mesh or a simulation.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_scalar::*;
use num_complex::Complex64;
use support::{TestBasis, TestQuadrature};

const TRI: ElementType = ElementType(303);
const QUAD: ElementType = ElementType(404);

fn tri() -> Vec<[f64; 3]> {
    vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
}

fn unit_quad() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
}

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-300), "{a:e} vs {b:e}");
}

#[test]
fn stat_elec_bulk_is_eps_stiffness_and_charge_load() {
    let xyz = tri();
    let mat = StatElecAtNodes {
        relative_permittivity: vec![3.0; 3],
        charge_density: Some(vec![2.0; 3]),
    };
    let mut cfg = StatElecConfig::new(2);
    cfg.eps0 = 0.5;
    let mut out = LocalSystem::zeros(0);
    stat_elec_bulk_local(&TestBasis, &TestQuadrature, TRI, &xyz, &mat, &cfg, &mut out).unwrap();
    // Linear shape functions on this triangle (area 1): grad phi = (-1/2, -1),
    // (1/2, 0), (0, 1). K = eps0 eps_r * area * G G^T.
    let g = [[-0.5, -1.0], [0.5, 0.0], [0.0, 1.0]];
    for p in 0..3 {
        for q in 0..3 {
            let want = 0.5 * 3.0 * 1.0 * (g[p][0] * g[q][0] + g[p][1] * g[q][1]);
            close(out.stiff[p * 3 + q], want, 1e-13);
        }
        close(out.force[p], 2.0 / 3.0, 1e-13);
    }
    assert!(out.mass.is_none());
}

#[test]
fn stat_elec_axisymmetric_weights_with_two_pi_r() {
    let xyz = vec![
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let mat = StatElecAtNodes {
        relative_permittivity: vec![1.0; 4],
        charge_density: Some(vec![1.0; 4]),
    };
    let mut cfg = StatElecConfig::new(2);
    cfg.eps0 = 1.0;
    let mut cart = LocalSystem::zeros(0);
    stat_elec_bulk_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &mat,
        &cfg,
        &mut cart,
    )
    .unwrap();
    cfg.axisymmetric = true;
    let mut axi = LocalSystem::zeros(0);
    stat_elec_bulk_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &mat,
        &cfg,
        &mut axi,
    )
    .unwrap();
    // The total load is int 2 pi r dr dz = 2 pi * 1.5.
    close(
        axi.force.iter().sum::<f64>(),
        2.0 * std::f64::consts::PI * 1.5,
        1e-13,
    );
    close(cart.force.iter().sum::<f64>(), 1.0, 1e-13);
}

#[test]
fn stat_elec_layer_boundary_gives_robin_terms() {
    // A unit-length edge with alpha = eps_l / h = 4 and beta = alpha V + rho h / (2 eps0).
    let xyz = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let bc = StatElecBcAtNodes {
        electric_flux: Some(vec![1.0, 1.0]),
        layer_relative_permittivity: Some(vec![2.0; 2]),
        layer_thickness: Some(vec![0.5; 2]),
        electrode_potential: Some(vec![3.0; 2]),
        layer_charge_density: Some(vec![8.0; 2]),
    };
    let mut cfg = StatElecConfig::new(2);
    cfg.eps0 = 2.0;
    let mut out = LocalSystem::zeros(0);
    stat_elec_boundary_local(
        &TestBasis,
        &TestQuadrature,
        ElementType(202),
        &xyz,
        &bc,
        &cfg,
        &mut out,
    )
    .unwrap();
    // Mass matrix of a unit edge: [[1/3, 1/6], [1/6, 1/3]].
    close(out.stiff[0], 4.0 / 3.0, 1e-13);
    close(out.stiff[1], 4.0 / 6.0, 1e-13);
    let beta = 1.0 + 4.0 * 3.0 + 0.5 * 8.0 * 0.5 / 2.0;
    close(out.force[0], beta / 2.0, 1e-13);
    // A layer without a thickness is refused.
    let bad = StatElecBcAtNodes {
        layer_thickness: None,
        ..bc
    };
    assert!(stat_elec_boundary_local(
        &TestBasis,
        &TestQuadrature,
        ElementType(202),
        &xyz,
        &bad,
        &cfg,
        &mut out
    )
    .is_err());
}

#[test]
fn potential_post_recovers_uniform_field_and_energy() {
    // V = 3 x on a unit square: E = (-3, 0), D = eps0 eps_r E, E.D = 9 eps.
    let xyz = unit_quad();
    let pot: Vec<f64> = xyz.iter().map(|p| 3.0 * p[0]).collect();
    let cfg = PostConfig {
        dim: 2,
        axisymmetric: false,
        quad_points: None,
    };
    let mut out = PotentialPostLocal::default();
    potential_post_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &pot,
        &Coefficient::Scalar(vec![2.0; 4]),
        0.5,
        &cfg,
        &mut out,
    )
    .unwrap();
    close(out.volume, 1.0, 1e-13);
    close(out.heat_integral, 9.0, 1e-13);
    close(out.potential_integral, 1.5, 1e-13);
    // The projection of the constant field: M^{-1} b = -grad V at every node.
    let lu = DenseLu::new(out.mass.clone(), 4).unwrap();
    let mut ex = out.field[0].clone();
    lu.solve(&mut ex);
    let mut dx = out.flux[0].clone();
    lu.solve(&mut dx);
    let mut h = out.heat.clone();
    lu.solve(&mut h);
    for i in 0..4 {
        close(ex[i], -3.0, 1e-12);
        close(dx[i], -3.0, 1e-12);
        close(h[i], 9.0, 1e-12);
    }
    assert!(out.field[1].iter().all(|v| v.abs() < 1e-14));
}

#[test]
fn tensor_coefficient_applies_to_the_gradient() {
    let xyz = unit_quad();
    let pot: Vec<f64> = xyz.iter().map(|p| p[0] + 2.0 * p[1]).collect();
    let cfg = PostConfig {
        dim: 2,
        axisymmetric: false,
        quad_points: None,
    };
    let mut out = PotentialPostLocal::default();
    let c = Coefficient::Full(vec![[[2.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 1.0]]; 4]);
    potential_post_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &pot,
        &c,
        1.0,
        &cfg,
        &mut out,
    )
    .unwrap();
    // grad V = (1, 2); C grad V = (4, 7); E.J = 1*4 + 2*7 = 18.
    close(out.heat_integral, 18.0, 1e-13);
    close(out.flux[0].iter().sum::<f64>(), -4.0, 1e-13);
    close(out.flux[1].iter().sum::<f64>(), -7.0, 1e-13);
}

#[test]
fn dense_lu_pivots_like_lapack() {
    // A zero leading entry and several later interchanges.
    let n = 4;
    let a = vec![
        0.0, 2.0, 1.0, -1.0, //
        3.0, 0.5, 4.0, 2.0, //
        1.0, -3.0, 0.0, 5.0, //
        2.0, 1.0, -2.0, 0.25,
    ];
    let x_true = [1.0, -2.0, 0.5, 3.0];
    let mut b: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * x_true[j]).sum())
        .collect();
    DenseLu::new(a, n).unwrap().solve(&mut b);
    for (x, t) in b.iter().zip(x_true) {
        close(*x, t, 1e-12);
    }
    assert!(DenseLu::new(vec![1.0, 2.0, 2.0, 4.0], 2).is_err());
}

fn helmholtz_tri(omega: f64, c: f64, rho: f64) -> HelmholtzLocal {
    let xyz = tri();
    let mat = HelmholtzBulkAtNodes {
        sound_speed: vec![c; 3],
        density: vec![rho; 3],
        damping: Some(vec![0.1; 3]),
        load: Some(vec![Complex64::new(1.0, 2.0); 3]),
    };
    let mut out = HelmholtzLocal::default();
    helmholtz_bulk_local(
        &TestBasis,
        &TestQuadrature,
        TRI,
        &xyz,
        &mat,
        &HelmholtzConfig::new(2, omega),
        &mut out,
    )
    .unwrap();
    out
}

#[test]
fn helmholtz_bulk_is_stiffness_minus_k2_mass_over_density() {
    let (omega, c, rho) = (6.0, 3.0, 2.0);
    let out = helmholtz_tri(omega, c, rho);
    let k = omega / c;
    let g = [[-0.5, -1.0], [0.5, 0.0], [0.0, 1.0]];
    // Mass matrix of a linear triangle of area 1: (1 + delta_pq) / 12.
    for p in 0..3 {
        for q in 0..3 {
            let mass = if p == q { 2.0 } else { 1.0 } / 12.0;
            let grad = g[p][0] * g[q][0] + g[p][1] * g[q][1];
            let want = Complex64::new(grad - k * k * mass, k * 0.1 * mass) / rho;
            let got = out.stiff[p * 3 + q];
            close(got.re, want.re, 1e-12);
            close(got.im, want.im, 1e-12);
        }
        let f = out.force[p];
        close(f.re, 1.0 / 3.0 / rho, 1e-13);
        close(f.im, 2.0 / 3.0 / rho, 1e-13);
    }
}

#[test]
fn complex_local_system_is_interleaved_as_two_by_two_blocks() {
    let out = helmholtz_tri(6.0, 3.0, 2.0);
    let mut real = LocalSystem::zeros(0);
    out.to_real(&mut real);
    assert_eq!(real.n, 6);
    for p in 0..3 {
        for q in 0..3 {
            let a = out.stiff[p * 3 + q];
            let at = |r: usize, c: usize| real.stiff[(2 * p + r) * 6 + 2 * q + c];
            assert_eq!(
                (at(0, 0), at(0, 1), at(1, 0), at(1, 1)),
                (a.re, -a.im, a.im, a.re)
            );
        }
        assert_eq!(
            (real.force[2 * p], real.force[2 * p + 1]),
            (out.force[p].re, out.force[p].im)
        );
    }
}

#[test]
fn helmholtz_boundary_admittance() {
    let xyz = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let omega = 5.0;
    let bc = HelmholtzBcAtNodes {
        density: vec![2.0; 2],
        impedance: Some(vec![Complex64::new(3.0, 4.0); 2]),
        wave_flux: Some(vec![Complex64::new(1.0, -1.0); 2]),
    };
    let mut out = HelmholtzLocal::default();
    helmholtz_boundary_local(
        &TestBasis,
        &TestQuadrature,
        ElementType(202),
        &xyz,
        &bc,
        &HelmholtzConfig::new(2, omega),
        &mut out,
    )
    .unwrap();
    let y = Complex64::new(0.0, omega) / Complex64::new(3.0, 4.0);
    let want = y / 2.0 / 3.0;
    close(out.stiff[0].re, want.re, 1e-13);
    close(out.stiff[0].im, want.im, 1e-13);
    close(out.force[0].re, 0.5 / 2.0, 1e-13);
    close(out.force[0].im, -0.5 / 2.0, 1e-13);

    // No impedance: no admittance, only the flux.
    let bc = HelmholtzBcAtNodes {
        impedance: None,
        ..bc
    };
    helmholtz_boundary_local(
        &TestBasis,
        &TestQuadrature,
        ElementType(202),
        &xyz,
        &bc,
        &HelmholtzConfig::new(2, omega),
        &mut out,
    )
    .unwrap();
    assert!(out.stiff.iter().all(|v| v.norm() == 0.0));
}

#[test]
fn wrong_nodal_lengths_are_refused() {
    let mat = StatElecAtNodes {
        relative_permittivity: vec![1.0; 2],
        charge_density: None,
    };
    let mut out = LocalSystem::zeros(0);
    let e = stat_elec_bulk_local(
        &TestBasis,
        &TestQuadrature,
        TRI,
        &tri(),
        &mat,
        &StatElecConfig::new(2),
        &mut out,
    );
    assert!(matches!(e, Err(ScalarError::NodalLength { .. })));
}

#[test]
fn stat_current_kernel_stiffness_mass_and_robin() {
    let xyz = tri();
    let mat = StatCurrentAtNodes {
        conductivity: Coefficient::Diagonal(vec![[2.0, 3.0, 0.0]; 3]),
        relative_permittivity: Some(vec![4.0; 3]),
        current_source: Some(vec![1.0; 3]),
    };
    let mut cfg = StatCurrentConfig::new(2);
    cfg.eps0 = 0.25;
    cfg.transient = true;
    let mut out = LocalSystem::zeros(0);
    stat_current_bulk_local(&TestBasis, &TestQuadrature, TRI, &xyz, &mat, &cfg, &mut out).unwrap();
    let g = [[-0.5, -1.0], [0.5, 0.0], [0.0, 1.0]];
    let mass = out.mass.as_ref().unwrap();
    for p in 0..3 {
        for q in 0..3 {
            close(
                out.stiff[p * 3 + q],
                2.0 * g[p][0] * g[q][0] + 3.0 * g[p][1] * g[q][1],
                1e-13,
            );
            close(
                mass[p * 3 + q],
                0.25 * 4.0 * (g[p][0] * g[q][0] + g[p][1] * g[q][1]),
                1e-13,
            );
        }
    }
    cfg.transient = false;
    stat_current_bulk_local(&TestBasis, &TestQuadrature, TRI, &xyz, &mat, &cfg, &mut out).unwrap();
    assert!(out.mass.is_none());

    let bc = StatCurrentBcAtNodes {
        current_density: Some(vec![1.0; 2]),
        external_conductivity: Some(vec![6.0; 2]),
        external_potential: Some(vec![2.0; 2]),
    };
    stat_current_boundary_local(
        &TestBasis,
        &TestQuadrature,
        ElementType(202),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        &bc,
        &cfg,
        &mut out,
    )
    .unwrap();
    close(out.stiff[0], 6.0 / 3.0, 1e-13);
    close(out.force[0], (1.0 + 12.0) / 2.0, 1e-13);
}

#[test]
fn procedures_are_registered_under_their_elmer_names() {
    let mut reg = elmer_solver::SolverRegistry::new();
    register(&mut reg);
    for (lib, entry) in [
        ("StatElecSolve", "StatElecSolver"),
        ("StatElecSolveVec", "StatElecSolver"),
        ("StatCurrentSolve", "StatCurrentSolver"),
        ("HelmholtzSolve", "HelmholtzSolver"),
        ("helmholtzsolve.so", "helmholtzsolver"),
    ] {
        assert!(reg.get(lib, entry).is_some(), "{lib} {entry}");
    }
}
