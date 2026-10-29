mod support;

use approx::assert_abs_diff_eq;
use elmer_core::{LocalSystem, NullLogger};
use elmer_fixtures::poisson_assembly::{
    AssemblyError, BoundaryProps, BulkProps, Element, PoissonAssembly,
};
use support::{RecordingLogger, TestBasis, TestQuadrature, LINE, TRI};

const UNIT_TRI: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
const AREA: f64 = 0.5;

fn tri() -> Element<'static> {
    Element { etype: TRI, nodes: &UNIT_TRI }
}

fn bulk(props: BulkProps<'_>) -> LocalSystem {
    let mut local = LocalSystem::zeros(3);
    PoissonAssembly::new()
        .bulk(tri(), props, &TestBasis, &TestQuadrature, &NullLogger, &mut local)
        .unwrap();
    local
}

fn assert_all_close(got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert_abs_diff_eq!(g, w, epsilon = 1e-14);
    }
}

#[test]
fn bulk_conduction_matches_the_p1_stiffness() {
    let local = bulk(BulkProps { conductivity: Some(&[2.0; 3]), ..Default::default() });
    // grad phi = (-1,-1), (1,0), (0,1); K = cond * area * g_p . g_q
    let want = [2.0, -1.0, -1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 1.0].map(|v| v * 2.0 * AREA);
    assert_all_close(&local.stiff, &want);
    assert_eq!(local.force, vec![0.0; 3]);
}

#[test]
fn bulk_reaction_and_source_use_the_consistent_mass() {
    let local = bulk(BulkProps {
        reaction: Some(&[3.0; 3]),
        source: Some(&[1.0, 2.0, 3.0]),
        ..Default::default()
    });
    // M = area/12 (1 + delta_pq); F_p = area/12 (sum s + s_p)
    let m = AREA / 12.0;
    let want_k: Vec<f64> = (0..9).map(|i| 3.0 * m * if i % 4 == 0 { 2.0 } else { 1.0 }).collect();
    assert_all_close(&local.stiff, &want_k);
    assert_all_close(&local.force, &[m * (6.0 + 1.0), m * (6.0 + 2.0), m * (6.0 + 3.0)]);
}

#[test]
fn bulk_interpolates_nodal_conductivity_to_the_points() {
    // A conductivity linear in x: cond(x) = 1 + x. On this triangle
    // K00 = integral of (1 + x) * 2 = 2 * (area + int x) = 2 * (1/2 + 1/6).
    let local = bulk(BulkProps { conductivity: Some(&[1.0, 2.0, 1.0]), ..Default::default() });
    assert_abs_diff_eq!(local.stiff[0], 2.0 * (0.5 + 1.0 / 6.0), epsilon = 1e-14);
}

#[test]
fn bulk_absent_keywords_contribute_zero() {
    let local = bulk(BulkProps::default());
    assert_eq!(local.stiff, vec![0.0; 9]);
    assert_eq!(local.force, vec![0.0; 3]);
}

#[test]
fn bulk_accumulates_into_the_given_system() {
    let mut local = LocalSystem::zeros(3);
    local.force = vec![1.0, 1.0, 1.0];
    let mut k = PoissonAssembly::new();
    let p = BulkProps { source: Some(&[6.0; 3]), ..Default::default() };
    k.bulk(tri(), p, &TestBasis, &TestQuadrature, &NullLogger, &mut local).unwrap();
    k.bulk(tri(), p, &TestBasis, &TestQuadrature, &NullLogger, &mut local).unwrap();
    // each call adds source * area / 3 = 1 per node
    assert_all_close(&local.force, &[3.0; 3]);
}

#[test]
fn boundary_adds_the_robin_terms_on_a_line() {
    let nodes = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let el = Element { etype: LINE, nodes: &nodes };
    let mut local = LocalSystem::zeros(2);
    let p = BoundaryProps {
        flux: Some(&[5.0; 2]),
        coefficient: Some(&[3.0; 2]),
        external: Some(&[7.0; 2]),
    };
    let added = PoissonAssembly::new()
        .boundary(el, p, &TestBasis, &TestQuadrature, &NullLogger, &mut local)
        .unwrap();
    assert!(added);
    // length 2: K = coeff * L/6 [[2,1],[1,2]], F_p = L/2 (flux + coeff ext)
    assert_all_close(&local.stiff, &[2.0, 1.0, 1.0, 2.0].map(|v| v * 3.0 * 2.0 / 6.0));
    assert_all_close(&local.force, &[26.0, 26.0]);
}

#[test]
fn boundary_with_only_flux_has_no_stiffness() {
    let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let el = Element { etype: LINE, nodes: &nodes };
    let mut local = LocalSystem::zeros(2);
    let p = BoundaryProps { flux: Some(&[4.0; 2]), ..Default::default() };
    assert!(PoissonAssembly::new()
        .boundary(el, p, &TestBasis, &TestQuadrature, &NullLogger, &mut local)
        .unwrap());
    assert_eq!(local.stiff, vec![0.0; 4]);
    assert_all_close(&local.force, &[2.0, 2.0]);
}

#[test]
fn boundary_returns_early_without_flux_or_coefficient() {
    let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let el = Element { etype: LINE, nodes: &nodes };
    let mut local = LocalSystem::zeros(2);
    // External alone does not count as a boundary load.
    let p = BoundaryProps { external: Some(&[9.0; 2]), ..Default::default() };
    assert!(!PoissonAssembly::new()
        .boundary(el, p, &TestBasis, &TestQuadrature, &NullLogger, &mut local)
        .unwrap());
    assert_eq!(local.stiff, vec![0.0; 4]);
    assert_eq!(local.force, vec![0.0; 2]);
}

#[test]
fn first_time_message_is_once_per_routine_and_instance() {
    let log = RecordingLogger::default();
    let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let line = Element { etype: LINE, nodes: &nodes };
    let flux = BoundaryProps { flux: Some(&[1.0; 2]), ..Default::default() };
    let mut k = PoissonAssembly::new();
    for _ in 0..3 {
        k.bulk(tri(), BulkProps::default(), &TestBasis, &TestQuadrature, &log, &mut LocalSystem::zeros(3))
            .unwrap();
        k.boundary(line, flux, &TestBasis, &TestQuadrature, &log, &mut LocalSystem::zeros(2)).unwrap();
    }
    let seen = log.0.lock().unwrap().clone();
    assert_eq!(
        seen,
        vec![
            ("PoissonBulkAssembly".to_string(), "1st time".to_string(), 4),
            ("PoissonBoundaryAssembly".to_string(), "1st time".to_string(), 4),
        ]
    );
    // A second instance reports again, like a second process.
    PoissonAssembly::new()
        .bulk(tri(), BulkProps::default(), &TestBasis, &TestQuadrature, &log, &mut LocalSystem::zeros(3))
        .unwrap();
    assert_eq!(log.0.lock().unwrap().len(), 3);
}

#[test]
fn mismatched_shapes_are_rejected() {
    let mut k = PoissonAssembly::new();
    let mut small = LocalSystem::zeros(2);
    let err = k
        .bulk(tri(), BulkProps::default(), &TestBasis, &TestQuadrature, &NullLogger, &mut small)
        .unwrap_err();
    assert!(matches!(err, AssemblyError::Shape { what: "LocalSystem", .. }));
    let mut local = LocalSystem::zeros(3);
    let err = k
        .bulk(
            tri(),
            BulkProps { source: Some(&[1.0, 2.0]), ..Default::default() },
            &TestBasis,
            &TestQuadrature,
            &NullLogger,
            &mut local,
        )
        .unwrap_err();
    assert!(matches!(err, AssemblyError::Shape { what: "Source", .. }));
}

#[test]
fn basis_failure_names_the_element_type() {
    let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let err = PoissonAssembly::new()
        .bulk(
            Element { etype: TRI, nodes: &flat },
            BulkProps::default(),
            &TestBasis,
            &TestQuadrature,
            &NullLogger,
            &mut LocalSystem::zeros(3),
        )
        .unwrap_err();
    assert!(matches!(err, AssemblyError::Basis { etype, .. } if etype == TRI));
}
