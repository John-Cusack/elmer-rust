//! The element kernel on single elements, against closed forms.

mod support;

use elmer_core::ElementType;
use elmer_physics_flux::{flux_bulk_local, FluxCoefficient, FluxConfig, FluxError, FluxLocal};
use support::{TestBasis, TestQuadrature};

const QUAD: ElementType = ElementType(404);
const TRI: ElementType = ElementType(303);

/// `[x0, x0+w] x [y0, y0+h]`, counter-clockwise from the lower left.
fn rect(x0: f64, y0: f64, w: f64, h: f64) -> Vec<[f64; 3]> {
    vec![
        [x0, y0, 0.0],
        [x0 + w, y0, 0.0],
        [x0 + w, y0 + h, 0.0],
        [x0, y0 + h, 0.0],
    ]
}

fn run(
    etype: ElementType,
    xyz: &[[f64; 3]],
    u: &[f64],
    coeff: &FluxCoefficient,
    cfg: &FluxConfig,
) -> FluxLocal {
    let mut out = FluxLocal::default();
    flux_bulk_local(
        &TestBasis,
        &TestQuadrature,
        etype,
        xyz,
        u,
        coeff,
        cfg,
        &mut out,
    )
    .unwrap();
    out
}

/// Solves `M x = f` by Gaussian elimination.
fn solve(m: &[f64], f: &[f64]) -> Vec<f64> {
    let n = f.len();
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut r = m[i * n..(i + 1) * n].to_vec();
            r.push(f[i]);
            r
        })
        .collect();
    for c in 0..n {
        let p = (c..n)
            .max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))
            .unwrap();
        a.swap(c, p);
        for i in c + 1..n {
            let k = a[i][c] / a[c][c];
            for j in c..=n {
                a[i][j] -= k * a[c][j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = (a[i][n] - s) / a[i][i];
    }
    x
}

fn cfg(flux: bool, grad: bool, flux_mag: bool, grad_mag: bool) -> FluxConfig {
    FluxConfig {
        calculate_flux: flux,
        calculate_grad: grad,
        calculate_flux_magnitude: flux_mag,
        calculate_grad_magnitude: grad_mag,
        ..FluxConfig::new(2)
    }
}

#[test]
fn dofs_and_right_hand_side_order() {
    let c = cfg(true, true, true, true);
    assert_eq!(c.dofs(), 6);
    assert_eq!(c.first_magnitude(), 4);
    assert_eq!(cfg(false, true, false, true).dofs(), 3);
    assert_eq!(cfg(false, false, false, false).dofs(), 0);
    let c3 = FluxConfig {
        calculate_flux: true,
        ..FluxConfig::new(3)
    };
    assert_eq!(c3.dofs(), 3);
}

#[test]
fn mass_matrix_of_a_rectangle() {
    // Q1 consistent mass: (w h / 36) [4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4].
    let xyz = rect(1.0, 0.0, 2.0, 3.0);
    let out = run(
        QUAD,
        &xyz,
        &[0.0; 4],
        &FluxCoefficient::Identity,
        &cfg(false, true, false, false),
    );
    let s = 2.0 * 3.0 / 36.0;
    let want = [
        4., 2., 1., 2., 2., 4., 2., 1., 1., 2., 4., 2., 2., 1., 2., 4.,
    ];
    for (g, w) in out.mass.iter().zip(want) {
        assert!((g - s * w).abs() < 1e-14, "{g} vs {}", s * w);
    }
}

#[test]
fn linear_potential_projects_exactly() {
    // u = 2 x - 3 y + 1: the gradient (2, -3) is in the element space, so the
    // local projection M g = f returns it at every node.
    let xyz = rect(0.5, -1.0, 1.5, 0.75);
    let u: Vec<f64> = xyz.iter().map(|p| 2.0 * p[0] - 3.0 * p[1] + 1.0).collect();
    let out = run(
        QUAD,
        &xyz,
        &u,
        &FluxCoefficient::Identity,
        &cfg(false, true, false, true),
    );
    assert_eq!(out.force.len(), 3);
    let gx = solve(&out.mass, &out.force[0]);
    let gy = solve(&out.mass, &out.force[1]);
    let gm = solve(&out.mass, &out.force[2]);
    for i in 0..4 {
        assert!((gx[i] - 2.0).abs() < 1e-13);
        assert!((gy[i] + 3.0).abs() < 1e-13);
        assert!((gm[i] - 13f64.sqrt()).abs() < 1e-13);
    }
}

#[test]
fn triangle_with_scalar_coefficient_gives_minus_k_grad() {
    let xyz = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let u: Vec<f64> = xyz.iter().map(|p| p[0] + 4.0 * p[1]).collect();
    let k = FluxCoefficient::Scalar(vec![3.0; 3]);
    let out = run(TRI, &xyz, &u, &k, &cfg(true, false, true, false));
    let fx = solve(&out.mass, &out.force[0]);
    let fy = solve(&out.mass, &out.force[1]);
    let fm = solve(&out.mass, &out.force[2]);
    for i in 0..3 {
        assert!((fx[i] + 3.0).abs() < 1e-13, "{fx:?}");
        assert!((fy[i] + 12.0).abs() < 1e-13);
        // Flux_mag is the magnitude of C grad, positive.
        assert!((fm[i] - 3.0 * 17f64.sqrt()).abs() < 1e-12);
    }
}

#[test]
fn diagonal_and_full_tensors() {
    let xyz = rect(0.0, 0.0, 1.0, 1.0);
    let u: Vec<f64> = xyz.iter().map(|p| p[0] + 2.0 * p[1]).collect(); // grad (1, 2)
    let c = cfg(true, false, false, false);

    let diag = FluxCoefficient::Diagonal(vec![[2.0, 5.0, 9.0]; 4]);
    let out = run(QUAD, &xyz, &u, &diag, &c);
    let (fx, fy) = (
        solve(&out.mass, &out.force[0]),
        solve(&out.mass, &out.force[1]),
    );
    assert!((fx[0] + 2.0).abs() < 1e-13 && (fy[0] + 10.0).abs() < 1e-13);

    let k = [[2.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 1.0]];
    let full = FluxCoefficient::Full(vec![k; 4]);
    let out = run(QUAD, &xyz, &u, &full, &c);
    let (fx, fy) = (
        solve(&out.mass, &out.force[0]),
        solve(&out.mass, &out.force[1]),
    );
    // -(2*1 + 1*2, 1*1 + 3*2) = (-4, -7)
    assert!((fx[2] + 4.0).abs() < 1e-13 && (fy[2] + 7.0).abs() < 1e-13);
}

#[test]
fn variable_coefficient_is_interpolated_at_integration_points() {
    // k = 1 + x on [0,1]^2, u = y: flux_y = -(1 + x). The right-hand side at
    // the Gauss points sees the interpolated k, so the projection is the
    // linear function -(1 + x) exactly.
    let xyz = rect(0.0, 0.0, 1.0, 1.0);
    let u: Vec<f64> = xyz.iter().map(|p| p[1]).collect();
    let k = FluxCoefficient::Scalar(xyz.iter().map(|p| 1.0 + p[0]).collect());
    let out = run(QUAD, &xyz, &u, &k, &cfg(true, false, false, false));
    let fy = solve(&out.mass, &out.force[1]);
    for (v, p) in fy.iter().zip(&xyz) {
        assert!((v + 1.0 + p[0]).abs() < 1e-13, "{fy:?}");
    }
}

#[test]
fn flux_mag_is_not_the_magnitude_of_projected_flux() {
    // For a non-constant gradient the projection of |C grad| differs from
    // |projection of C grad|: FluxSolver projects the pointwise magnitude.
    let xyz = rect(0.0, 0.0, 1.0, 1.0);
    let u: Vec<f64> = xyz.iter().map(|p| p[0] * p[1]).collect(); // grad (y, x)
    let out = run(
        QUAD,
        &xyz,
        &u,
        &FluxCoefficient::Identity,
        &cfg(false, true, false, true),
    );
    let (gx, gy) = (
        solve(&out.mass, &out.force[0]),
        solve(&out.mass, &out.force[1]),
    );
    let gm = solve(&out.mass, &out.force[2]);
    let i = 2; // the corner (1, 1)
    assert!((gx[i] - 1.0).abs() < 1e-13 && (gy[i] - 1.0).abs() < 1e-13);
    assert!(gm[i] < 2f64.sqrt() - 1e-3, "{gm:?}");
}

#[test]
fn axisymmetric_weight_is_the_radius() {
    // Sum of the mass matrix is the integral of r over [1,3] x [0,2] = 8.
    let xyz = rect(1.0, 0.0, 2.0, 2.0);
    let mut c = cfg(false, true, false, false);
    c.axisymmetric = true;
    let out = run(QUAD, &xyz, &[0.0; 4], &FluxCoefficient::Identity, &c);
    let total: f64 = out.mass.iter().sum();
    assert!((total - 8.0).abs() < 1e-13, "{total}");
    c.axisymmetric = false;
    let out = run(QUAD, &xyz, &[0.0; 4], &FluxCoefficient::Identity, &c);
    assert!((out.mass.iter().sum::<f64>() - 4.0).abs() < 1e-13);
}

#[test]
fn bad_inputs_are_refused() {
    let xyz = rect(0.0, 0.0, 1.0, 1.0);
    let mut out = FluxLocal::default();
    let c = cfg(true, false, false, false);
    let e = flux_bulk_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &[0.0; 3],
        &FluxCoefficient::Identity,
        &c,
        &mut out,
    )
    .unwrap_err();
    assert!(matches!(
        e,
        FluxError::NodalLength {
            what: "potential",
            ..
        }
    ));
    let e = flux_bulk_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &[0.0; 4],
        &FluxCoefficient::Scalar(vec![1.0; 2]),
        &c,
        &mut out,
    )
    .unwrap_err();
    assert!(matches!(
        e,
        FluxError::NodalLength {
            what: "flux coefficient",
            ..
        }
    ));
    let c1 = FluxConfig::new(1);
    let e = flux_bulk_local(
        &TestBasis,
        &TestQuadrature,
        QUAD,
        &xyz,
        &[0.0; 4],
        &FluxCoefficient::Identity,
        &c1,
        &mut out,
    )
    .unwrap_err();
    assert!(matches!(e, FluxError::Invalid(_)));
}
