//! Postprocessing kernels and integration-point inputs against closed forms.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_elastic::post::{
    principal_values, stress_strain_local, tresca, von_mises, MassSystem, StressPostInput,
};
use elmer_physics_elastic::{stress_bulk_local, CoordSystem, ElasticLaw, IpProps, StressBulkInput};
use support::{Gauss, P1Basis};

const QUAD: [[f64; 3]; 4] = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.2, 1.0, 0.0], [0.1, 1.1, 0.0]];

fn iso(e: f64, nu: f64) -> ElasticLaw {
    ElasticLaw::Isotropic { young: vec![e; 4], poisson: vec![nu; 4] }
}

fn post<'a>(law: &'a ElasticLaw, disp: &'a [[f64; 3]]) -> StressPostInput<'a> {
    StressPostInput {
        etype: ElementType(404),
        xyz: &QUAD,
        cdim: 2,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law,
        disp,
        temperature: None,
        heat_expansion: None,
        ip: None,
        want_strain: true,
        want_stress: true,
        np: None,
    }
}

#[test]
fn principal_values_of_a_known_tensor() {
    // Rotated diagonal diag(5, 2, -1): the eigenvalues survive the rotation.
    let (c, s) = (0.6_f64, 0.8_f64);
    let r = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
    let d = [5.0, 2.0, -1.0];
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = (0..3).map(|k| r[i][k] * d[k] * r[j][k]).sum();
        }
    }
    let w = principal_values(&t);
    for (a, b) in w.iter().zip(d) {
        assert!((a - b).abs() < 1e-12, "{w:?}");
    }
    assert!((tresca(&w) - 3.0).abs() < 1e-12);
    // Uniaxial stress: von Mises is its magnitude.
    assert!((von_mises(&[[7.0, 0.0, 0.0], [0.0; 3], [0.0; 3]]) - 7.0).abs() < 1e-12);
}

#[test]
fn rigid_translation_and_free_thermal_expansion_carry_no_stress() {
    let law = iso(200.0, 0.3);
    let disp = [[0.3, -0.2, 0.0]; 4];
    let out = stress_strain_local(&P1Basis, &Gauss, &post(&law, &disp)).unwrap();
    assert!(out.strain.iter().chain(&out.stress).all(|v| v.abs() < 1e-13));

    // u = alpha dT x: the total strain is exactly the thermal one.
    let (alpha, dt) = (1.2e-5, 40.0);
    let disp: Vec<[f64; 3]> = QUAD.iter().map(|p| [alpha * dt * p[0], alpha * dt * p[1], 0.0]).collect();
    let temp = [dt; 4];
    let beta = [[alpha; 3]; 4];
    let mut inp = post(&law, &disp);
    inp.temperature = Some(&temp);
    inp.heat_expansion = Some(&beta);
    let out = stress_strain_local(&P1Basis, &Gauss, &inp).unwrap();
    assert!(out.strain.iter().chain(&out.stress).all(|v| v.abs() < 1e-12), "{:?}", out.stress);
}

#[test]
fn plane_strain_uniaxial_stretch_gives_the_out_of_plane_stress() {
    let (e, nu) = (210.0, 0.25);
    let law = iso(e, nu);
    let eps = 1e-3;
    let disp: Vec<[f64; 3]> = QUAD.iter().map(|p| [eps * p[0], 0.0, 0.0]).collect();
    let out = stress_strain_local(&P1Basis, &Gauss, &post(&law, &disp)).unwrap();
    // Components (11, 22, 33, 12) per node; mass-weighted sums recover the nodal value.
    let m: f64 = out.mass.iter().sum();
    let comp = |k: usize| (0..4).map(|p| out.stress[4 * p + k]).sum::<f64>() / m;
    let lame1 = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    assert!((comp(0) - (lame1 + 2.0 * mu) * eps).abs() < 1e-12);
    assert!((comp(1) - lame1 * eps).abs() < 1e-12);
    assert!((comp(2) - lame1 * eps).abs() < 1e-12);
    assert!(comp(3).abs() < 1e-12);
}

#[test]
fn mass_system_solves_the_projection_exactly() {
    // Two 1D linear elements on slots 0-1-2 reproduce a linear field.
    let mut sys = MassSystem::new(3);
    let h = 0.5;
    let local = [h / 3.0, h / 6.0, h / 6.0, h / 3.0];
    sys.add(&[Some(0), Some(1)], &local);
    sys.add(&[Some(1), Some(2)], &local);
    let f = |x: f64| 2.0 + 3.0 * x;
    // int N_i f over each element for linear f.
    let rhs = |a: f64, b: f64| [h * (2.0 * f(a) + f(b)) / 6.0, h * (f(a) + 2.0 * f(b)) / 6.0];
    let (r1, r2) = (rhs(0.0, 0.5), rhs(0.5, 1.0));
    let b = [r1[0], r1[1] + r2[0], r2[1]];
    let mut x = [0.0; 3];
    sys.solve(&b, &mut x).unwrap();
    for (xi, want) in x.iter().zip([f(0.0), f(0.5), f(1.0)]) {
        assert!((xi - want).abs() < 1e-12, "{x:?}");
    }
}

#[test]
fn integration_point_data_equal_to_nodal_data_changes_nothing() {
    let zeros = vec![[0.0; 3]; 4];
    let dens = vec![1.0; 4];
    let law = iso(150.0, 0.2);
    let mut load = zeros.clone();
    load.iter_mut().for_each(|l| *l = [0.7, -0.4, 0.0]);
    let make = |ip: Option<&IpProps>, with_load: bool| {
        let inp = StressBulkInput {
            etype: ElementType(404),
            xyz: &QUAD,
            cdim: 2,
            coord: CoordSystem::Cartesian,
            plane_stress: false,
            law: &law,
            density: &dens,
            damping: None,
            load: if with_load { &load } else { &zeros },
            pressure_load: None,
            stress_load: None,
            strain_load: None,
            temperature: None,
            heat_expansion: None,
            want_mass: false,
            np: None,
            ip,
        };
        let mut out = LocalSystem::zeros(0);
        stress_bulk_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
        out
    };
    let nodal = make(None, true);
    let ip = IpProps {
        young: Some(vec![150.0; 4]),
        poisson: Some(vec![0.2; 4]),
        beta: None,
        load: Some(vec![[0.7, -0.4, 0.0]; 4]),
        pressure: Some(vec![0.0; 4]),
    };
    let at_ip = make(Some(&ip), false);
    for (a, b) in nodal.stiff.iter().zip(&at_ip.stiff).chain(nodal.force.iter().zip(&at_ip.force)) {
        assert!((a - b).abs() <= 1e-13 * a.abs().max(1.0), "{a} {b}");
    }
    let bad = IpProps { young: Some(vec![1.0; 3]), ..IpProps::default() };
    let inp = StressBulkInput {
        etype: ElementType(404),
        xyz: &QUAD,
        cdim: 2,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law: &law,
        density: &dens,
        damping: None,
        load: &zeros,
        pressure_load: None,
        stress_load: None,
        strain_load: None,
        temperature: None,
        heat_expansion: None,
        want_mass: false,
        np: None,
        ip: Some(&bad),
    };
    assert!(stress_bulk_local(&P1Basis, &Gauss, &inp, &mut LocalSystem::zeros(0)).is_err());
}
