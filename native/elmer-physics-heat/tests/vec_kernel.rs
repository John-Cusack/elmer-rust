//! `heat_bulk_local_vec` (LocalMatrixVec's summation order) against the
//! point-by-point `heat_bulk_local`: the same matrices up to rounding.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_heat::{
    heat_bulk_local, heat_bulk_local_vec, CompressModel, Compressibility, Conductivity,
    HeatBulkConfig, HeatMaterialAtNodes, JouleHeat, Stabilization,
};
use support::{TestBasis, TestQuadrature};

const QUAD: ElementType = ElementType(404);

/// A distorted quadrilateral.
fn quad() -> Vec<[f64; 3]> {
    vec![
        [0.1, 0.0, 0.0],
        [1.2, 0.1, 0.0],
        [1.0, 0.9, 0.0],
        [0.0, 1.1, 0.0],
    ]
}

fn varying(base: f64) -> Vec<f64> {
    (0..4).map(|i| base * (1.0 + 0.1 * i as f64)).collect()
}

fn both(mat: &HeatMaterialAtNodes, cfg: &HeatBulkConfig) -> (LocalSystem, LocalSystem) {
    let run = |vec: bool| {
        let mut out = LocalSystem::zeros(0);
        let f = if vec {
            heat_bulk_local_vec
        } else {
            heat_bulk_local
        };
        f(
            &TestBasis,
            &TestQuadrature,
            QUAD,
            &quad(),
            mat,
            cfg,
            &mut out,
        )
        .unwrap();
        out
    };
    (run(true), run(false))
}

fn agree(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() <= 1e-13 * scale, "{x} vs {y}");
    }
}

fn check(mat: &HeatMaterialAtNodes, cfg: &HeatBulkConfig) {
    let (v, p) = both(mat, cfg);
    agree(&v.stiff, &p.stiff);
    agree(&v.force, &p.force);
    assert_eq!(v.mass.is_some(), p.mass.is_some());
    if let (Some(a), Some(b)) = (&v.mass, &p.mass) {
        agree(a, b);
    }
}

fn convected(computed: bool, compressibility: Option<Compressibility>) -> HeatMaterialAtNodes {
    HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(varying(0.5))),
        density: Some(varying(2.0)),
        heat_capacity: Some(varying(4.0)),
        heat_source: Some(varying(1.5)),
        convection_velocity: Some(
            (0..4)
                .map(|i| [2.0 + 0.3 * i as f64, 1.0 - 0.2 * i as f64, 0.0])
                .collect(),
        ),
        convection_computed: computed,
        mesh_velocity: Some(vec![[0.1, -0.05, 0.0]; 4]),
        heat_transfer_multiplier: Some(vec![1.5; 4]),
        compressibility,
        ..Default::default()
    }
}

/// Every configuration on both `LinearForms` paths: 4 points (the plain
/// loops) and 16 (the BLAS order).
fn configs() -> Vec<HeatBulkConfig> {
    let mut out = Vec::new();
    for np in [4, 16] {
        for transient in [false, true] {
            for supg in [false, true] {
                for axisymmetric in [false, true] {
                    let mut cfg = HeatBulkConfig::new(2);
                    cfg.quad_points = Some(np);
                    cfg.transient = transient;
                    cfg.axisymmetric = axisymmetric;
                    if supg {
                        cfg.stabilization = Stabilization::Supg {
                            h_k: 0.7,
                            m_k: 1.0 / 3.0,
                        };
                    }
                    out.push(cfg);
                }
            }
        }
    }
    out
}

#[test]
fn diffusion_and_sources() {
    for cfg in configs() {
        let mut mat = HeatMaterialAtNodes {
            conductivity: Some(Conductivity::Scalar(varying(3.0))),
            density: Some(varying(2.0)),
            heat_capacity: Some(varying(5.0)),
            heat_source: Some(varying(0.7)),
            ..Default::default()
        };
        check(&mat, &cfg);
        mat.volumetric_heat_source = Some(varying(-1.1));
        check(&mat, &cfg);
        mat.conductivity = Some(Conductivity::Diagonal(vec![[2.0, 0.5, 0.0]; 4]));
        check(&mat, &cfg);
        mat.conductivity = Some(Conductivity::Full(vec![
            [
                [2.0, 0.3, 0.0],
                [0.1, 0.5, 0.0],
                [0.0; 3]
            ];
            4
        ]));
        check(&mat, &cfg);
    }
}

#[test]
fn convection_constant_and_computed() {
    for cfg in configs() {
        check(&convected(false, None), &cfg);
        check(&convected(true, None), &cfg);
    }
}

#[test]
fn compressible_with_pressure_work() {
    let pressure = vec![1.0, 2.0, 3.5, 2.5];
    let thermal = |prev: Option<(Vec<f64>, f64)>| Compressibility {
        model: CompressModel::Thermal {
            heat_expansion: Some(varying(0.01)),
            reference_temperature: Some(vec![1.0; 4]),
        },
        temperature: varying(3.0),
        pressure: Some(pressure.clone()),
        prev_pressure: prev,
    };
    let gas = |prev: Option<(Vec<f64>, f64)>| Compressibility {
        model: CompressModel::PerfectGas {
            specific_heat_ratio: Some(vec![1.4; 4]),
            reference_pressure: Some(vec![3.0; 4]),
        },
        temperature: varying(2.0),
        pressure: Some(pressure.clone()),
        prev_pressure: prev,
    };
    for cfg in configs() {
        let prev = cfg.transient.then(|| (vec![0.9, 2.1, 3.0, 2.0], 0.25));
        check(&convected(true, Some(thermal(prev.clone()))), &cfg);
        check(&convected(true, Some(gas(prev))), &cfg);
    }
}

/// A keyword whose nodal values agree is used as that constant, as
/// `ListGetElementRealVec` returns a constant keyword, so `rho * source` is
/// exactly 1 at every point and the force is `sum_t phi_p(t) w(t)` with no
/// interpolation rounding. A solution field is always interpolated.
#[test]
fn a_constant_keyword_is_used_exactly() {
    let mut cfg = HeatBulkConfig::new(2);
    cfg.quad_points = Some(4);
    let mat = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(vec![1.0; 4])),
        density: Some(vec![1.0; 4]),
        heat_source: Some(vec![1.0; 4]),
        ..Default::default()
    };
    let (v, _) = both(&mat, &cfg);
    let unit = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(vec![1.0; 4])),
        volumetric_heat_source: Some(vec![1.0; 4]),
        ..Default::default()
    };
    let (u, _) = both(&unit, &cfg);
    assert_eq!(v.force, u.force);
    assert!(v.force.iter().all(|f| *f > 0.0));
}

/// Joule heat from a potential: `sigma |grad V|^2` at every point, added
/// against each basis function. On the unit triangle with `V = 2x + 3y` and
/// `sigma = 0.5` that is 6.5 everywhere, so each node gets 6.5 * area / 3.
#[test]
fn joule_heat_from_the_potential_gradient() {
    let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mat = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(vec![1.0; 3])),
        joule_heat: Some(JouleHeat::Potential {
            electric_conductivity: vec![0.5; 3],
            potential: vec![0.0, 2.0, 3.0],
        }),
        ..Default::default()
    };
    let mut out = LocalSystem::zeros(0);
    heat_bulk_local_vec(
        &TestBasis,
        &TestQuadrature,
        ElementType(303),
        &tri,
        &mat,
        &HeatBulkConfig::new(2),
        &mut out,
    )
    .unwrap();
    for f in &out.force {
        assert!((f - 6.5 / 6.0).abs() < 1e-14, "{f}");
    }
    // A field instead: sigma * field = 0.5 * 4.
    let mat = HeatMaterialAtNodes {
        joule_heat: Some(JouleHeat::Field {
            electric_conductivity: vec![0.5; 3],
            field: vec![4.0; 3],
        }),
        ..mat
    };
    heat_bulk_local_vec(
        &TestBasis,
        &TestQuadrature,
        ElementType(303),
        &tri,
        &mat,
        &HeatBulkConfig::new(2),
        &mut out,
    )
    .unwrap();
    for f in &out.force {
        assert!((f - 2.0 / 6.0).abs() < 1e-14, "{f}");
    }
}
