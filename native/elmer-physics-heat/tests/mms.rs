//! Manufactured-solution convergence: `-div(k grad u) = f` on the unit square
//! with `u = sin(pi x) sin(pi y)` and homogeneous Dirichlet conditions.

mod support;

use std::f64::consts::PI;

use elmer_core::{BasisEval, ElementBasis, Quadrature};
use elmer_physics_heat::{Conductivity, HeatBoundaryConfig, HeatBulkConfig, HeatMaterialAtNodes};
use support::glue::{assemble, cg};
use support::meshio::{unit_square, Elem};
use support::{TestBasis, TestQuadrature};

fn exact(x: f64, y: f64) -> f64 {
    (PI * x).sin() * (PI * y).sin()
}

fn l2_error(n: usize, quad: bool, second_order: bool) -> f64 {
    let mesh = unit_square(n, quad, second_order);
    let nn = mesh.coords.len();
    let perm: Vec<usize> = (0..nn).collect();
    let k = 2.0;
    let sys = assemble(
        &mesh,
        &perm,
        &TestBasis,
        &TestQuadrature,
        &HeatBulkConfig::new(2),
        &|e: &Elem| HeatMaterialAtNodes {
            conductivity: Some(Conductivity::Scalar(vec![k; e.nodes.len()])),
            volumetric_heat_source: Some(
                e.nodes
                    .iter()
                    .map(|&i| {
                        let c = mesh.coords[i];
                        2.0 * k * PI * PI * exact(c[0], c[1])
                    })
                    .collect(),
            ),
            ..Default::default()
        },
        &HeatBoundaryConfig::default(),
        &|_| None,
    );

    // Symmetric elimination of the homogeneous Dirichlet boundary.
    let mut fixed = vec![false; nn];
    for e in &mesh.boundary {
        for &i in &e.nodes {
            fixed[i] = true;
        }
    }
    let mut a = std::collections::BTreeMap::new();
    for (&(i, j), &v) in &sys.a {
        if !fixed[i] && !fixed[j] {
            a.insert((i, j), v);
        }
    }
    let mut b = sys.b.clone();
    for i in 0..nn {
        if fixed[i] {
            a.insert((i, i), 1.0);
            b[i] = 0.0;
        }
    }
    let u = cg(nn, &a, &b);

    let mut err2 = 0.0;
    let mut be = BasisEval::default();
    for e in &mesh.bulk {
        let xyz = mesh.xyz(e);
        let np = if quad { 16 } else { 6 };
        for qp in TestQuadrature.gauss_points(e.etype, Some(np)) {
            TestBasis.eval(e.etype, &xyz, qp.uvw, &mut be).unwrap();
            let (mut x, mut y, mut uh) = (0.0, 0.0, 0.0);
            for (m, &i) in e.nodes.iter().enumerate() {
                x += be.basis[m] * mesh.coords[i][0];
                y += be.basis[m] * mesh.coords[i][1];
                uh += be.basis[m] * u[i];
            }
            err2 += qp.weight * be.det_j * (uh - exact(x, y)).powi(2);
        }
    }
    err2.sqrt()
}

fn rate(quad: bool, second_order: bool, expected: f64) {
    let e: Vec<f64> = [8, 16, 32]
        .iter()
        .map(|&n| l2_error(n, quad, second_order))
        .collect();
    let r1 = (e[0] / e[1]).log2();
    let r2 = (e[1] / e[2]).log2();
    assert!(
        (r2 - expected).abs() <= 0.1,
        "quad={quad} p2={second_order}: errors {e:?}, rates {r1:.3} {r2:.3}, expected {expected}"
    );
}

#[test]
fn p1_triangle_l2_rate_is_two() {
    rate(false, false, 2.0);
}

#[test]
fn p2_triangle_l2_rate_is_three() {
    rate(false, true, 3.0);
}

#[test]
fn p1_quad_l2_rate_is_two() {
    rate(true, false, 2.0);
}

#[test]
fn p2_quad_l2_rate_is_three() {
    rate(true, true, 3.0);
}
