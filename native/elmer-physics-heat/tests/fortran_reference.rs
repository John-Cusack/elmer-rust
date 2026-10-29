//! Globally assembled matrix and right-hand side against the pristine
//! ElmerSolver's `Linear System Save` dump (see `data/README.md`), for
//! `fem/tests/MinimalExample` and `fem/tests/ThermalPipe2D`.
//!
//! The assembly here is a test-local glue loop over the heat kernels, so a
//! pass shows the kernels reproduce HeatSolve's element matrices; the
//! bitwise-identity goal needs the real `elmer-fe` quadrature and belongs to
//! the oracle lane.

mod support;

use std::collections::HashMap;
use std::path::PathBuf;

use elmer_physics_heat::{
    Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, HeatMaterialAtNodes,
};
use support::glue::{assemble, matrix_diff, read_matrix, read_perm, read_vector};
use support::meshio::{read_elmer_mesh, Elem};
use support::{TestBasis, TestQuadrature};

fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
}

struct Case {
    dir: &'static str,
    conductivity: f64,
    density: f64,
    /// boundary tag -> flux
    flux: &'static [(i32, f64)],
    /// boundary tag -> Dirichlet temperature
    dirichlet: &'static [(i32, f64)],
}

fn check(case: Case) {
    let dir = data(case.dir);
    let mesh = read_elmer_mesh(&dir, 2);
    let perm = read_perm(&dir.join("linsys_perm.dat"));
    assert_eq!(perm.len(), mesh.coords.len());

    let flux: HashMap<i32, f64> = case.flux.iter().copied().collect();
    let mut sys = assemble(
        &mesh,
        &perm,
        &TestBasis,
        &TestQuadrature,
        &HeatBulkConfig::new(2),
        &|e: &Elem| HeatMaterialAtNodes {
            conductivity: Some(Conductivity::Scalar(vec![case.conductivity; e.nodes.len()])),
            density: Some(vec![case.density; e.nodes.len()]),
            ..Default::default()
        },
        &HeatBoundaryConfig::default(),
        &|e: &Elem| {
            flux.get(&e.tag).map(|&q| HeatBcAtNodes {
                heat_flux: Some(vec![q; e.nodes.len()]),
                ..Default::default()
            })
        },
    );

    let dirichlet: HashMap<i32, f64> = case.dirichlet.iter().copied().collect();
    let mut fixed = 0;
    for e in &mesh.boundary {
        if let Some(&v) = dirichlet.get(&e.tag) {
            for &node in &e.nodes {
                sys.dirichlet_row_scaled(perm[node], v);
                fixed += 1;
            }
        }
    }
    assert!(fixed > 0);

    let a_ref = read_matrix(&dir.join("linsys_a.dat"));
    let b_ref = read_vector(&dir.join("linsys_b.dat"));
    assert_eq!(b_ref.len(), sys.b.len());

    let (diff, scale) = matrix_diff(&sys.a, &a_ref);
    assert!(
        diff <= 1e-13 * scale,
        "matrix: max abs diff {diff:e}, scale {scale:e}"
    );
    let b_scale = b_ref.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(scale);
    let b_diff = sys
        .b
        .iter()
        .zip(&b_ref)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(
        b_diff <= 1e-13 * b_scale,
        "rhs: max abs diff {b_diff:e}, scale {b_scale:e}"
    );
}

#[test]
fn minimal_example_matrix_and_rhs_match_elmer() {
    check(Case {
        dir: "minimal_example",
        conductivity: 1.0,
        density: 1.0,
        flux: &[(3, 1.0)],
        dirichlet: &[(2, 0.0)],
    });
}

#[test]
fn thermal_pipe_2d_matrix_and_rhs_match_elmer() {
    check(Case {
        dir: "thermal_pipe2d",
        conductivity: 1.0,
        density: 1000.0,
        flux: &[],
        dirichlet: &[
            (9, 10.0),
            (10, 10.0),
            (11, 10.0),
            (12, 10.0),
            (1, 0.0),
            (2, 0.0),
            (3, 0.0),
            (4, 0.0),
        ],
    });
}
