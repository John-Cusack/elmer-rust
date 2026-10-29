//! Kernel-level checks that do not need a pristine run: the closed forms of
//! `ElementDiameter`/`StabParam`, normal orientation, the pressure
//! stabilization block, and the Newton linearization's consistency with
//! Picard at the linearization point.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_flow::{
    element_diameter, flow_bulk_local, outward_normal, stab_params, Discretization, FlowBulkInput, FlowError,
    StabParams,
};
use support::{TestBasis, TestQuadrature};

#[test]
fn element_diameter_closed_forms() {
    let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    // 16 A^2 / (3 sum |x_i - c|^2) = 16 * 0.25 / (3 * 4/3) = 1.
    assert!((element_diameter(ElementType(303), &tri, false).unwrap() - 1.0).abs() < 1e-15);
    let quad = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
    // 2 |e1|^2 |e4|^2 / (|e1|^2 + |e4|^2) = 2 * 4 * 1 / 5.
    assert!((element_diameter(ElementType(404), &quad, false).unwrap() - (8.0f64 / 5.0).sqrt()).abs() < 1e-15);
    let tet = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
    assert_eq!(element_diameter(ElementType(504), &tet, false).unwrap(), 1.0);
    assert_eq!(element_diameter(ElementType(504), &tet, true).unwrap(), 13f64.sqrt());
}

#[test]
fn stabilization_constant_by_type() {
    let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert_eq!(stab_params(ElementType(303), &tri, false).unwrap().mk, 1.0 / 3.0);
    let brick: Vec<[f64; 3]> = (0..8)
        .map(|i| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64])
        .collect();
    assert_eq!(stab_params(ElementType(808), &brick, false).unwrap().mk, 1.0 / 6.0);
    let tri6 = [[0.0; 3]; 6];
    assert!(matches!(stab_params(ElementType(306), &tri6, false), Err(FlowError::Unsupported(_))));
}

#[test]
fn normals_point_away_from_the_parent() {
    let bottom = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let n = outward_normal(ElementType(202), &bottom, [0.0; 3], [0.5, 0.5, 0.0]).unwrap();
    assert!((n[0]).abs() < 1e-15 && (n[1] + 1.0).abs() < 1e-15);
    let face = [[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 1.0, 1.0], [2.0, 0.0, 1.0]];
    let n = outward_normal(ElementType(404), &face, [0.3, -0.2, 0.0], [1.5, 0.5, 0.5]).unwrap();
    assert!((n[0] - 1.0).abs() < 1e-15 && n[1].abs() < 1e-15 && n[2].abs() < 1e-15);
}

fn triangle() -> Vec<[f64; 3]> {
    vec![[0.1, 0.0, 0.0], [1.0, 0.2, 0.0], [0.3, 0.9, 0.0]]
}

fn quad() -> Vec<[f64; 3]> {
    vec![[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [1.0, 1.1, 0.0], [-0.1, 0.9, 0.0]]
}

#[allow(clippy::too_many_arguments)]
fn input<'a>(
    etype: u16,
    xyz: &'a [[f64; 3]],
    velocity: &'a [[f64; 3]],
    load: &'a [[f64; 4]],
    props: &'a [f64],
    convect: bool,
    newton: bool,
    mk: f64,
) -> FlowBulkInput<'a> {
    let et = ElementType(etype);
    let mut stab = stab_params(et, xyz, false).unwrap();
    stab.mk = mk;
    FlowBulkInput {
        etype: et,
        xyz,
        dim: 2,
        density: props,
        viscosity: props,
        velocity,
        mesh_velocity: None,
        load,
        convect,
        newton,
        stab,
        discretization: Discretization::default(),
        want_mass: false,
        np: None,
    }
}

/// For P1 triangles without convection and with constant viscosity, the
/// pressure-pressure block is exactly the stabilization `Tau (grad q,
/// grad p)` with `Tau = mK hK^2 / (8 mu)`; with `mK = 0` (no `Mesh %
/// Stabilize`) it vanishes.
#[test]
fn pressure_block_is_tau_laplacian() {
    let xyz = triangle();
    let velocity = vec![[0.0; 3]; 3];
    let load = vec![[0.0; 4]; 3];
    let mu = vec![0.7; 3];
    let mut out = LocalSystem::zeros(0);
    let inp = input(303, &xyz, &velocity, &load, &mu, false, false, 1.0 / 3.0);
    flow_bulk_local(&TestBasis, &TestQuadrature, &inp, &mut out).unwrap();

    let StabParams { hk, mk } = inp.stab;
    let tau = mk * hk * hk / (8.0 * 0.7);
    let (a, b, c) = (xyz[0], xyz[1], xyz[2]);
    let det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
    let grads = [
        [(b[1] - c[1]) / det, (c[0] - b[0]) / det],
        [(c[1] - a[1]) / det, (a[0] - c[0]) / det],
        [(a[1] - b[1]) / det, (b[0] - a[0]) / det],
    ];
    let area = det.abs() / 2.0;
    for p in 0..3 {
        for q in 0..3 {
            let want = tau * area * (grads[p][0] * grads[q][0] + grads[p][1] * grads[q][1]);
            let got = out.stiff[(p * 3 + 2) * 9 + q * 3 + 2];
            assert!((got - want).abs() < 1e-14 * want.abs().max(1.0), "({p},{q}): {got} vs {want}");
        }
    }

    let inp = input(303, &xyz, &velocity, &load, &mu, false, false, 0.0);
    flow_bulk_local(&TestBasis, &TestQuadrature, &inp, &mut out).unwrap();
    for p in 0..3 {
        for q in 0..3 {
            assert_eq!(out.stiff[(p * 3 + 2) * 9 + q * 3 + 2], 0.0);
        }
    }
}

/// Newton adds `rho (grad u_k) delta u` to the matrix and `rho (grad u_k) u_k`
/// to the force, in the Galerkin and in the stabilization terms, so at the
/// linearization point `A_newton u_k - f_newton = A_picard u_k - f_picard`.
#[test]
fn newton_residual_equals_picard_residual_at_the_iterate() {
    for (etype, xyz) in [(303u16, triangle()), (404, quad())] {
        let n = xyz.len();
        let velocity: Vec<[f64; 3]> = (0..n)
            .map(|k| [0.3 + 0.7 * k as f64, -0.4 + 0.5 * (k * k) as f64, 0.0])
            .collect();
        let pressure: Vec<f64> = (0..n).map(|k| 1.0 - 0.3 * k as f64).collect();
        let load: Vec<[f64; 4]> = (0..n).map(|k| [0.5, -1.0 + k as f64, 0.0, 0.0]).collect();
        let props: Vec<f64> = (0..n).map(|k| 0.05 + 0.01 * k as f64).collect();
        let x: Vec<f64> = (0..n)
            .flat_map(|k| [velocity[k][0], velocity[k][1], pressure[k]])
            .collect();
        let residual = |newton: bool| {
            let mut out = LocalSystem::zeros(0);
            let inp = input(etype, &xyz, &velocity, &load, &props, true, newton, 1.0 / 3.0);
            flow_bulk_local(&TestBasis, &TestQuadrature, &inp, &mut out).unwrap();
            let nd = out.n;
            (0..nd)
                .map(|r| (0..nd).map(|c| out.stiff[r * nd + c] * x[c]).sum::<f64>() - out.force[r])
                .collect::<Vec<f64>>()
        };
        let (rn, rp) = (residual(true), residual(false));
        let scale = rp.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for (a, b) in rn.iter().zip(&rp) {
            assert!((a - b).abs() < 1e-12 * scale, "{etype}: {a} vs {b}");
        }
    }
}
