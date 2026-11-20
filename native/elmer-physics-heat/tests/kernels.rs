//! Element matrices against independently known integrals.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_heat::{
    heat_boundary_local, heat_bulk_local, CompressModel, Compressibility, Conductivity,
    HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, HeatError, HeatMaterialAtNodes,
    Stabilization,
};
use support::{TestBasis, TestQuadrature};

const TRI: ElementType = ElementType(303);
const LINE: ElementType = ElementType(202);

fn unit_tri() -> Vec<[f64; 3]> {
    vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
}

fn bulk(
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &HeatMaterialAtNodes,
    cfg: &HeatBulkConfig,
) -> LocalSystem {
    let mut out = LocalSystem::zeros(0);
    heat_bulk_local(&TestBasis, &TestQuadrature, etype, xyz, mat, cfg, &mut out).unwrap();
    out
}

fn scalar_k(k: f64, n: usize) -> HeatMaterialAtNodes {
    HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(vec![k; n])),
        ..Default::default()
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-13 * (1.0 + b.abs()), "{a} != {b}");
}

#[test]
fn p1_triangle_stiffness_is_k_times_the_textbook_matrix() {
    let out = bulk(TRI, &unit_tri(), &scalar_k(3.0, 3), &HeatBulkConfig::new(2));
    let expect = [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]];
    for p in 0..3 {
        for q in 0..3 {
            close(out.stiff[p * 3 + q], 3.0 * expect[p][q]);
        }
    }
    assert!(out.mass.is_none());
    assert_eq!(out.force, vec![0.0; 3]);
}

#[test]
fn transient_mass_is_rho_cp_consistent_mass_and_dt_is_not_applied() {
    let mut mat = scalar_k(1.0, 3);
    mat.density = Some(vec![2.0; 3]);
    mat.heat_capacity = Some(vec![5.0; 3]);
    let mut cfg = HeatBulkConfig::new(2);
    cfg.transient = true;
    let out = bulk(TRI, &unit_tri(), &mat, &cfg);
    let mass = out.mass.as_ref().unwrap();
    for p in 0..3 {
        for q in 0..3 {
            let m = if p == q { 2.0 } else { 1.0 };
            close(mass[p * 3 + q], 10.0 * m / 24.0);
        }
    }
    // Stiffness is untouched by the transient switch.
    let steady = bulk(TRI, &unit_tri(), &mat, &HeatBulkConfig::new(2));
    assert_eq!(steady.stiff, out.stiff);
}

#[test]
fn sources_volumetric_wins_and_heat_source_scales_with_density() {
    let mut mat = scalar_k(1.0, 3);
    mat.density = Some(vec![4.0; 3]);
    mat.heat_source = Some(vec![3.0; 3]);
    let out = bulk(TRI, &unit_tri(), &mat, &HeatBulkConfig::new(2));
    // rho * q * area / 3
    for f in &out.force {
        close(*f, 4.0 * 3.0 * 0.5 / 3.0);
    }
    mat.volumetric_heat_source = Some(vec![7.0; 3]);
    let out = bulk(TRI, &unit_tri(), &mat, &HeatBulkConfig::new(2));
    for f in &out.force {
        close(*f, 7.0 * 0.5 / 3.0);
    }
}

#[test]
fn full_tensor_conductivity_uses_off_diagonals_and_diagonal_matches_full() {
    let k = [[2.0, 0.7, 0.0], [0.7, 3.0, 0.0], [0.0, 0.0, 1.0]];
    let mat = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Full(vec![k; 3])),
        ..Default::default()
    };
    let out = bulk(TRI, &unit_tri(), &mat, &HeatBulkConfig::new(2));
    let g = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
    for p in 0..3 {
        for q in 0..3 {
            let mut a = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    a += k[i][j] * g[p][i] * g[q][j];
                }
            }
            close(out.stiff[p * 3 + q], 0.5 * a);
        }
    }
    // An off-diagonal really contributes.
    let diag_only = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Diagonal(vec![[2.0, 3.0, 1.0]; 3])),
        ..Default::default()
    };
    let d = bulk(TRI, &unit_tri(), &diag_only, &HeatBulkConfig::new(2));
    assert!((d.stiff[1] - out.stiff[1]).abs() > 0.1);
    // Diagonal(kxx,kyy) equals a Full tensor with zero off-diagonals.
    let kd = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]];
    let full_d = HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Full(vec![kd; 3])),
        ..Default::default()
    };
    let f = bulk(TRI, &unit_tri(), &full_d, &HeatBulkConfig::new(2));
    for (a, b) in d.stiff.iter().zip(&f.stiff) {
        close(*a, *b);
    }
}

#[test]
fn p1_tetrahedron_stiffness() {
    let xyz = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let out = bulk(
        ElementType(504),
        &xyz,
        &scalar_k(2.0, 4),
        &HeatBulkConfig::new(3),
    );
    close(out.stiff[0], 2.0 * 3.0 / 6.0);
    close(out.stiff[1], -2.0 / 6.0);
    close(out.stiff[5], 2.0 * 1.0 / 6.0);
    close(out.stiff[6], 0.0);
}

#[test]
fn axisymmetric_weights_by_radius() {
    // Triangle (1,0) (2,0) (1,1): centroid r = 4/3, area 1/2.
    let xyz = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
    let mut mat = scalar_k(1.0, 3);
    mat.density = Some(vec![1.0; 3]);
    mat.heat_capacity = Some(vec![3.0; 3]);
    let mut cfg = HeatBulkConfig::new(2);
    cfg.transient = true;
    cfg.axisymmetric = true;
    let out = bulk(TRI, &xyz, &mat, &cfg);
    let total: f64 = out.mass.as_ref().unwrap().iter().sum();
    close(total, 3.0 * 0.5 * 4.0 / 3.0);
}

#[test]
fn heat_transfer_multiplier_scales_the_whole_local_system() {
    let mut mat = scalar_k(1.0, 3);
    mat.volumetric_heat_source = Some(vec![1.0; 3]);
    mat.density = Some(vec![1.0; 3]);
    mat.heat_capacity = Some(vec![1.0; 3]);
    let mut cfg = HeatBulkConfig::new(2);
    cfg.transient = true;
    let base = bulk(TRI, &unit_tri(), &mat, &cfg);
    mat.heat_transfer_multiplier = Some(vec![2.5; 3]);
    let out = bulk(TRI, &unit_tri(), &mat, &cfg);
    for (a, b) in base.stiff.iter().zip(&out.stiff) {
        close(2.5 * a, *b);
    }
    for (a, b) in base.force.iter().zip(&out.force) {
        close(2.5 * a, *b);
    }
    for (a, b) in base
        .mass
        .as_ref()
        .unwrap()
        .iter()
        .zip(out.mass.as_ref().unwrap())
    {
        close(2.5 * a, *b);
    }
}

#[test]
fn galerkin_convection_reproduces_v_dot_grad_t_for_a_linear_field() {
    // T = x, v = (2, 0), rho cp = 3: (rho cp v.grad T, phi_p) = 6 * area/3.
    let mut mat = scalar_k(1.0, 3);
    mat.density = Some(vec![1.5; 3]);
    mat.heat_capacity = Some(vec![2.0; 3]);
    mat.convection_velocity = Some(vec![[2.0, 0.0, 0.0]; 3]);
    let out = bulk(TRI, &unit_tri(), &mat, &HeatBulkConfig::new(2));
    let mut diffusion = mat.clone();
    diffusion.convection_velocity = None;
    let diff = bulk(TRI, &unit_tri(), &diffusion, &HeatBulkConfig::new(2));
    let t = [0.0, 1.0, 0.0];
    for p in 0..3 {
        let r: f64 = (0..3)
            .map(|q| (out.stiff[p * 3 + q] - diff.stiff[p * 3 + q]) * t[q])
            .sum();
        close(r, 6.0 * 0.5 / 3.0);
    }
}

#[test]
fn supg_is_consistent_for_a_linear_solution() {
    // Source = rho cp v.grad T for T = a.x makes the convective residual of the
    // exact linear field vanish only if the SUPG force matches its stiffness.
    let xyz = vec![[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.2, 0.9, 0.0]];
    let a = [0.7, -0.3];
    let v = [1.3, 0.4];
    let (rho, cp) = (2.0, 3.0);
    let t: Vec<f64> = xyz.iter().map(|x| a[0] * x[0] + a[1] * x[1]).collect();
    let mut mat = scalar_k(0.05, 3);
    mat.density = Some(vec![rho; 3]);
    mat.heat_capacity = Some(vec![cp; 3]);
    mat.convection_velocity = Some(vec![[v[0], v[1], 0.0]; 3]);
    mat.volumetric_heat_source = Some(vec![rho * cp * (v[0] * a[0] + v[1] * a[1]); 3]);
    let mut cfg = HeatBulkConfig::new(2);
    cfg.stabilization = Stabilization::Supg {
        h_k: 1.0,
        m_k: 1.0 / 3.0,
    };
    let out = bulk(TRI, &xyz, &mat, &cfg);
    // Remove the diffusion part, which only carries boundary terms here.
    let mut diffusion = mat.clone();
    diffusion.convection_velocity = None;
    let diff = bulk(TRI, &xyz, &diffusion, &HeatBulkConfig::new(2));
    for p in 0..3 {
        let r: f64 = (0..3)
            .map(|q| (out.stiff[p * 3 + q] - diff.stiff[p * 3 + q]) * t[q])
            .sum::<f64>()
            - out.force[p];
        assert!(r.abs() < 1e-13, "row {p}: {r:e}");
    }
    // And the stabilisation actually changed the matrix.
    let plain = bulk(TRI, &xyz, &mat, &HeatBulkConfig::new(2));
    assert!((plain.stiff[0] - out.stiff[0]).abs() > 1e-6);
}

#[test]
fn nodal_length_mismatch_is_an_error() {
    let mut out = LocalSystem::zeros(0);
    let mat = scalar_k(1.0, 2);
    let e = heat_bulk_local(
        &TestBasis,
        &TestQuadrature,
        TRI,
        &unit_tri(),
        &mat,
        &HeatBulkConfig::new(2),
        &mut out,
    )
    .unwrap_err();
    assert_eq!(
        e,
        HeatError::NodalLength {
            what: "heat conductivity",
            expected: 3,
            got: 2
        }
    );
}

fn bnd(xyz: &[[f64; 3]], bc: &HeatBcAtNodes, cfg: &HeatBoundaryConfig) -> LocalSystem {
    let mut out = LocalSystem::zeros(0);
    heat_boundary_local(&TestBasis, &TestQuadrature, LINE, xyz, bc, cfg, &mut out).unwrap();
    out
}

#[test]
fn boundary_flux_and_robin_on_a_line() {
    let xyz = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let flux = HeatBcAtNodes {
        heat_flux: Some(vec![3.0; 2]),
        ..Default::default()
    };
    let out = bnd(&xyz, &flux, &HeatBoundaryConfig::default());
    assert_eq!(out.stiff, vec![0.0; 4]);
    close(out.force[0], 3.0);
    close(out.force[1], 3.0);

    let robin = HeatBcAtNodes {
        heat_transfer_coefficient: Some(vec![2.0; 2]),
        external_temperature: Some(vec![5.0; 2]),
        ..Default::default()
    };
    let out = bnd(&xyz, &robin, &HeatBoundaryConfig::default());
    // c L / 6 [[2,1],[1,2]] and c Text L / 2.
    close(out.stiff[0], 2.0 * 2.0 * 2.0 / 6.0);
    close(out.stiff[1], 2.0 * 2.0 / 6.0);
    close(out.force[0], 2.0 * 5.0);
    assert!(out.mass.is_none());
}

#[test]
fn boundary_coefficient_without_external_temperature_has_zero_source() {
    let xyz = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
    let bc = HeatBcAtNodes {
        heat_transfer_coefficient: Some(vec![1.0; 2]),
        ..Default::default()
    };
    let out = bnd(&xyz, &bc, &HeatBoundaryConfig::default());
    assert_eq!(out.force, vec![0.0; 2]);
    // Line of length sqrt(2): diagonal entry c L / 3.
    close(out.stiff[0], 2f64.sqrt() / 3.0);
}

#[test]
fn boundary_axisymmetric_flux_is_weighted_by_radius() {
    let xyz = vec![[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let bc = HeatBcAtNodes {
        heat_flux: Some(vec![1.0; 2]),
        ..Default::default()
    };
    let cfg = HeatBoundaryConfig {
        axisymmetric: true,
        ..Default::default()
    };
    let out = bnd(&xyz, &bc, &cfg);
    close(out.force[0], 5.0 / 3.0);
    close(out.force[1], 7.0 / 3.0);
}

/// A convected triangle with uniform properties and velocity `(2, 1)`.
fn convected(rho: f64, compressibility: Option<Compressibility>) -> HeatMaterialAtNodes {
    HeatMaterialAtNodes {
        conductivity: Some(Conductivity::Scalar(vec![0.5; 3])),
        density: Some(vec![rho; 3]),
        heat_capacity: Some(vec![4.0; 3]),
        heat_source: Some(vec![1.5; 3]),
        convection_velocity: Some(vec![[2.0, 1.0, 0.0]; 3]),
        compressibility,
        ..Default::default()
    }
}

fn supg_cfg(transient: bool) -> HeatBulkConfig {
    let mut cfg = HeatBulkConfig::new(2);
    cfg.transient = transient;
    cfg.stabilization = Stabilization::Supg {
        h_k: 0.7,
        m_k: 1.0 / 3.0,
    };
    cfg
}

fn same(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y);
    }
}

#[test]
fn thermal_compressibility_scales_the_density_everywhere_it_is_used() {
    // rho (1 - beta (T - Tref)) = 2 (1 - 0.25 (3 - 1)) = 1: every term,
    // SUPG and mass included, sees density 1. Without the flow pressure
    // there is no pressure work.
    let thermal = Compressibility {
        model: CompressModel::Thermal {
            heat_expansion: Some(vec![0.25; 3]),
            reference_temperature: Some(vec![1.0; 3]),
        },
        temperature: vec![3.0; 3],
        pressure: None,
        prev_pressure: None,
    };
    let got = bulk(
        TRI,
        &unit_tri(),
        &convected(2.0, Some(thermal)),
        &supg_cfg(true),
    );
    let want = bulk(TRI, &unit_tri(), &convected(1.0, None), &supg_cfg(true));
    same(&got.stiff, &want.stiff);
    same(&got.force, &want.force);
    same(got.mass.as_ref().unwrap(), want.mass.as_ref().unwrap());
}

#[test]
fn perfect_gas_density_and_pressure_work() {
    // p = 1 + x + 2y at the nodes, p0 = 3, gamma = 1.4, cp = 4, T = 2.
    let pressure = vec![1.0, 2.0, 3.0];
    let gas = |prev: Option<(Vec<f64>, f64)>| Compressibility {
        model: CompressModel::PerfectGas {
            specific_heat_ratio: Some(vec![1.4; 3]),
            reference_pressure: Some(vec![3.0; 3]),
        },
        temperature: vec![2.0; 3],
        pressure: Some(pressure.clone()),
        prev_pressure: prev,
    };
    let cfg = HeatBulkConfig::new(2);
    let got = bulk(TRI, &unit_tri(), &convected(9.0, Some(gas(None))), &cfg);

    // The density at each point is (p + p0) / (R T) with R = (gamma-1)/gamma
    // cp: linear in p, so the stiffness is that of the incompressible
    // element with the density interpolated from its nodal values.
    let r_t = (1.4 - 1.0) / 1.4 * 4.0 * 2.0;
    let rho_nodes: Vec<f64> = pressure.iter().map(|p| (p + 3.0) / r_t).collect();
    let mut want_mat = convected(0.0, None);
    want_mat.density = Some(rho_nodes.clone());
    let want = bulk(TRI, &unit_tri(), &want_mat, &cfg);
    same(&got.stiff, &want.stiff);

    // Force: the heat source times rho, plus Pcoeff = 1 times
    // v.grad p = 2*1 + 1*2 = 4 against each basis function (area/3 each).
    let mut want_force = want.force.clone();
    want_force.iter_mut().for_each(|f| *f += 4.0 / 6.0);
    same(&got.force, &want_force);

    // Transient: dp/dt = (p - p_prev)/dt = 0.5 everywhere adds 0.5/6 more.
    let prev: Vec<f64> = pressure.iter().map(|p| p - 1.0).collect();
    let mut tcfg = HeatBulkConfig::new(2);
    tcfg.transient = true;
    let got_t = bulk(
        TRI,
        &unit_tri(),
        &convected(9.0, Some(gas(Some((prev, 2.0))))),
        &tcfg,
    );
    let mut want_t = want_force.clone();
    want_t.iter_mut().for_each(|f| *f += 0.5 / 6.0);
    same(&got_t.force, &want_t);
}

#[test]
fn thermal_pressure_work_coefficient() {
    // Pcoeff = T beta / (1 - beta (T - Tref)) = 3 * 0.25 / 0.5 = 1.5 and
    // v.grad p = 4: the pressure work adds 6 / 6 to each node.
    let with = |pressure: Option<Vec<f64>>| Compressibility {
        model: CompressModel::Thermal {
            heat_expansion: Some(vec![0.25; 3]),
            reference_temperature: Some(vec![1.0; 3]),
        },
        temperature: vec![3.0; 3],
        pressure,
        prev_pressure: None,
    };
    let cfg = HeatBulkConfig::new(2);
    let got = bulk(
        TRI,
        &unit_tri(),
        &convected(2.0, Some(with(Some(vec![1.0, 2.0, 3.0])))),
        &cfg,
    );
    let base = bulk(TRI, &unit_tri(), &convected(2.0, Some(with(None))), &cfg);
    same(&got.stiff, &base.stiff);
    let want: Vec<f64> = base.force.iter().map(|f| f + 1.0).collect();
    same(&got.force, &want);
}

#[test]
fn perfect_gas_without_the_flow_pressure_is_an_error() {
    let gas = Compressibility {
        model: CompressModel::PerfectGas {
            specific_heat_ratio: None,
            reference_pressure: None,
        },
        temperature: vec![2.0; 3],
        pressure: None,
        prev_pressure: None,
    };
    let mut out = LocalSystem::zeros(0);
    let err = heat_bulk_local(
        &TestBasis,
        &TestQuadrature,
        TRI,
        &unit_tri(),
        &convected(1.0, Some(gas)),
        &HeatBulkConfig::new(2),
        &mut out,
    )
    .unwrap_err();
    assert!(err.to_string().contains("computed convection"), "{err}");
}
