//! Trait wiring and the behaviours the capture does not exercise: default
//! and relative rules, the coordinate-system table, degenerate and
//! ill-conditioned elements.

use elmer_core::{BasisEval, ElementBasis, ElementType, Quadrature};
use elmer_fe::{
    element_info, element_type, element_types, gauss_points, gauss_points_rel, stab_param,
    CoordinateKind, CoordinateSystem, ElmerQuadrature, FeError, LagrangeBasis, MetricStatus,
};

const LAGRANGE: [u16; 21] = [
    202, 203, 204, 303, 306, 310, 404, 408, 409, 412, 416, 504, 510, 605, 613, 706, 715, 718, 808,
    820, 827,
];

fn reference_nodes(code: u16) -> Vec<[f64; 3]> {
    let d = element_type(code).unwrap();
    (0..d.n_nodes)
        .map(|k| [d.node_u[k], d.node_v[k], d.node_w[k]])
        .collect()
}

#[test]
fn default_rules_integrate_the_reference_measure() {
    let measure = |family: u16| match family {
        2 => 2.0,
        3 => 0.5,
        4 => 4.0,
        5 => 1.0 / 6.0,
        6 => 4.0 / 3.0,
        7 => 1.0,
        _ => 8.0,
    };
    for code in LAGRANGE {
        let et = ElementType(code);
        let pts = ElmerQuadrature.gauss_points(et, None);
        let def = element_type(code).unwrap();
        assert_eq!(pts.len(), gauss_points(et, None).unwrap().len());
        let sum: f64 = pts.iter().map(|p| p.weight).sum();
        assert!(
            (sum - measure(def.family())).abs() < 1e-13,
            "{code}: weights sum to {sum}"
        );
    }
    assert_eq!(gauss_points(ElementType(808), None).unwrap().len(), 8);
    assert_eq!(gauss_points(ElementType(820), None).unwrap().len(), 27);
    assert_eq!(gauss_points(ElementType(303), Some(7)).unwrap().len(), 7);
}

#[test]
fn relative_rules_follow_elements_def() {
    let n = |code, rel| gauss_points_rel(ElementType(code), rel).unwrap().len();
    assert_eq!((n(303, -1), n(303, 0), n(303, 1)), (1, 3, 11));
    assert_eq!((n(504, -1), n(504, 0), n(504, 1)), (1, 4, 11));
    assert_eq!((n(404, 0), n(404, 1)), (4, 16));
    // InitializeElementDescriptions keeps the record of the previous
    // ELEMENT block: 306 gives two counts, so it inherits 303's
    // GaussPoints0 = 1 rather than defaulting to its own 6.
    assert_eq!(element_type(306).unwrap().gauss_points0, 1);
    assert_eq!(n(306, -1), 1);
    assert_eq!(
        gauss_points_rel(ElementType(303), 2),
        Err(FeError::InvalidRelOrder(2))
    );
}

#[test]
fn basis_trait_matches_element_info_and_partitions_unity() {
    for code in LAGRANGE {
        let def = element_type(code).unwrap();
        // A mildly sheared copy of the reference element in 3D.
        let xyz: Vec<[f64; 3]> = reference_nodes(code)
            .iter()
            .map(|p| [p[0] + 0.1 * p[1], p[1] + 0.05 * p[2], p[2] + 0.02 * p[0]])
            .collect();
        let cdim = def.dim.max(2);
        let rule = gauss_points(ElementType(code), None).unwrap();
        let basis = LagrangeBasis::new(cdim);
        for t in 0..rule.len() {
            let uvw = [rule.u[t], rule.v[t], rule.w[t]];
            let mut a = BasisEval::default();
            basis.eval(ElementType(code), &xyz, uvw, &mut a).unwrap();
            let mut b = BasisEval::default();
            let mut st = MetricStatus::Unchecked;
            element_info(def, &xyz, uvw, cdim, &mut st, &mut b).unwrap();
            assert_eq!(a.det_j.to_bits(), b.det_j.to_bits());
            assert_eq!(a.basis, b.basis);
            assert_eq!(a.dbasis_dx, b.dbasis_dx);
            let s: f64 = a.basis.iter().sum();
            assert!((s - 1.0).abs() < 1e-12, "{code}: sum of basis {s}");
            for j in 0..3 {
                let g: f64 = a.dbasis_dx.iter().map(|d| d[j]).sum();
                assert!(g.abs() < 1e-11, "{code}: sum of d/dx{j} = {g}");
            }
        }
    }
}

#[test]
fn errors_instead_of_fatal() {
    let mut ev = BasisEval::default();
    let basis = LagrangeBasis::default();
    assert!(basis
        .eval(ElementType(999), &[[0.0; 3]], [0.0; 3], &mut ev)
        .is_err());
    // Wrong node count.
    assert!(basis
        .eval(ElementType(303), &[[0.0; 3]; 2], [0.2, 0.2, 0.0], &mut ev)
        .is_err());
    // Collapsed triangle: degenerate in double and quad precision.
    let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let def = element_type(303).unwrap();
    let mut st = MetricStatus::Unchecked;
    let r = element_info(def, &flat, [0.2, 0.2, 0.0], 2, &mut st, &mut ev);
    assert!(
        matches!(r, Err(FeError::Degenerate { dim: 2, .. })),
        "{r:?}"
    );
    assert!(gauss_points(ElementType(404), Some(14 * 14)).is_err());
}

#[test]
fn ill_conditioned_element_takes_the_quad_precision_metric() {
    // A needle triangle far from the origin: det(J^T J) cancels badly in
    // double precision, so CheckMetric selects the quad-precision metric.
    let (o, h) = (1.0e7, 67.0 * 2f64.powi(-29));
    let xyz = [[o, o, 0.0], [o + 1.0, o, 0.0], [o + 1.0, o + h, 0.0]];
    let def = element_type(303).unwrap();
    let mut st = MetricStatus::Unchecked;
    let mut ev = BasisEval::default();
    element_info(def, &xyz, [0.2, 0.3, 0.0], 3, &mut st, &mut ev).unwrap();
    assert_eq!(st, MetricStatus::Quad);
    assert!((ev.det_j / h - 1.0).abs() < 1e-9, "detJ {}", ev.det_j);
    // dphi/dy of node 3 is 1/h on this triangle.
    assert!(
        (ev.dbasis_dx[2][1] * h - 1.0).abs() < 1e-9,
        "{:?}",
        ev.dbasis_dx[2]
    );
}

#[test]
fn every_defined_type_parses() {
    let codes: Vec<u16> = element_types().map(|t| t.code).collect();
    for c in LAGRANGE {
        assert!(codes.contains(&c), "{c}");
    }
    for t in element_types().filter(|t| t.code >= 200) {
        assert_eq!(t.basis.len(), t.n_nodes, "{}", t.code);
    }
}

#[test]
fn coordinate_systems() {
    let c = |s, d| CoordinateSystem::from_keyword(s, d).unwrap();
    assert_eq!(c(None, 2).dim, 2);
    assert_eq!(c(Some("Cartesian 3D"), 2).dim, 3);
    let axi = c(Some("Axi Symmetric"), 3);
    assert_eq!((axi.kind, axi.dim), (CoordinateKind::AxisSymmetric, 2));
    assert!(axi.is_axisymmetric());
    assert_eq!(axi.sqrt_metric(0.25, 1.0, 0.0), 0.25);
    assert_eq!(
        c(Some("cylindric symmetric"), 3).sqrt_metric(2.0, 0.0, 0.0),
        2.0
    );
    assert_eq!(c(Some("cartesian"), 3).sqrt_metric(2.0, 0.0, 0.0), 1.0);
    assert!(CoordinateSystem::from_keyword(Some("spherical"), 3).is_none());
    let info = axi.info(2.0, 0.0, 0.0);
    assert_eq!(info.metric[2][2], 0.25);
    assert_eq!(info.symbols[2][2][0], -2.0);
}

#[test]
fn stab_param_linear_closed_forms() {
    // Unit right triangle: hK^2 = 16 A^2 / (3 s) with A = 1/2, s = 4/3.
    let (mk, hk) = stab_param(ElementType(303), &reference_nodes(303), 2, false).unwrap();
    assert_eq!(mk, 1.0 / 3.0);
    assert!((hk - 1.0).abs() < 1e-15, "{hk}");
    // Brick: shortest / longest edge.
    let xyz: Vec<[f64; 3]> = reference_nodes(808)
        .iter()
        .map(|p| [p[0], 0.5 * p[1], 2.0 * p[2]])
        .collect();
    assert_eq!(
        stab_param(ElementType(808), &xyz, 3, false).unwrap(),
        (1.0 / 6.0, 1.0)
    );
    assert_eq!(
        stab_param(ElementType(808), &xyz, 3, true).unwrap(),
        (1.0 / 6.0, 4.0)
    );
}
