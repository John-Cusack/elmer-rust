//! Self-consistency of the element kernels: the stiffness the Newton assembly
//! returns is the exact derivative of the internal force, for St. Venant-Kirchhoff
//! and neo-Hookean kinematics, and without `Large Deflection` one step reproduces
//! the linear small-strain stiffness of `elmer-physics-elastic`.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_elastic::{CoordSystem, ElasticLaw};
use elmer_physics_elasticsolve::{elastic_bulk_local, neo_hookean_local, BulkInput, NeoHookeInput};
use support::{Gauss, P1Basis};

const TRI: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.3, 0.9, 0.0]];
const QUAD: [[f64; 3]; 4] = [[0.5, 0.0, 0.0], [1.6, 0.1, 0.0], [1.5, 1.1, 0.0], [0.4, 0.9, 0.0]];
const HEX: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.1, 0.0, 0.1],
    [1.0, 1.2, 0.0],
    [0.0, 1.0, 0.1],
    [0.1, 0.0, 1.0],
    [1.0, 0.1, 1.1],
    [1.1, 1.0, 1.0],
    [0.0, 1.1, 1.2],
];

/// A displacement field that is neither small nor symmetric.
fn displacement(n: usize, cdim: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|k| std::array::from_fn(|a| if a < cdim { 0.08 * ((1 + k * 3 + a * 5) % 7) as f64 - 0.2 } else { 0.0 }))
        .collect()
}

fn iso(n: usize) -> ElasticLaw {
    ElasticLaw::Isotropic { young: vec![200.0; n], poisson: vec![0.3; n] }
}

fn bulk(etype: u16, xyz: &[[f64; 3]], cdim: usize, coord: CoordSystem, large: bool, disp: &[[f64; 3]]) -> LocalSystem {
    let n = xyz.len();
    let law = iso(n);
    let (density, load) = (vec![1.0; n], vec![[0.0; 3]; n]);
    let inp = BulkInput {
        etype: ElementType(etype),
        xyz,
        cdim,
        coord,
        plane_stress: false,
        large_deflection: large,
        law: &law,
        disp,
        density: &density,
        damping: None,
        load: &load,
        inertial: None,
        pressure_load: None,
        stress_load: None,
        strain_load: None,
        temperature: None,
        heat_expansion: None,
        want_mass: false,
        np: None,
        ip: None,
    };
    let mut out = LocalSystem::zeros(0);
    elastic_bulk_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
    out
}

fn neo(etype: u16, xyz: &[[f64; 3]], cdim: usize, coord: CoordSystem, disp: &[[f64; 3]]) -> LocalSystem {
    let n = xyz.len();
    let (young, poisson, density, load) = (vec![200.0; n], vec![0.3; n], vec![1.0; n], vec![[0.0; 3]; n]);
    let inp = NeoHookeInput {
        etype: ElementType(etype),
        xyz,
        cdim,
        coord,
        plane_stress: false,
        young: &young,
        poisson: &poisson,
        disp,
        density: &density,
        damping: None,
        load: &load,
        inertial: None,
        want_mass: false,
        np: None,
    };
    let mut out = LocalSystem::zeros(0);
    neo_hookean_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
    out
}

/// The internal force `K(u) u - F(u)` of the Newton form, for `u` flattened node-major.
fn internal(local: &LocalSystem, u: &[f64]) -> Vec<f64> {
    (0..local.n)
        .map(|i| (0..local.n).map(|j| local.stiff[i * local.n + j] * u[j]).sum::<f64>() - local.force[i])
        .collect()
}

/// `|K - dR/du|` by central differences, relative to `|K|`.
fn tangent_error(assemble: &dyn Fn(&[[f64; 3]]) -> LocalSystem, n: usize, cdim: usize) -> f64 {
    let disp = displacement(n, cdim);
    let flat = |d: &[[f64; 3]]| -> Vec<f64> { d.iter().flat_map(|v| v[..cdim].to_vec()).collect() };
    let unflat = |u: &[f64]| -> Vec<[f64; 3]> {
        u.chunks(cdim).map(|c| std::array::from_fn(|a| if a < cdim { c[a] } else { 0.0 })).collect()
    };
    let u0 = flat(&disp);
    let k = assemble(&disp).stiff;
    let h = 1e-6;
    let (mut err, mut scale) = (0.0f64, 0.0f64);
    for j in 0..u0.len() {
        let (mut up, mut um) = (u0.clone(), u0.clone());
        up[j] += h;
        um[j] -= h;
        let (rp, rm) = (internal(&assemble(&unflat(&up)), &up), internal(&assemble(&unflat(&um)), &um));
        for i in 0..u0.len() {
            let fd = (rp[i] - rm[i]) / (2.0 * h);
            err = err.max((fd - k[i * u0.len() + j]).abs());
            scale = scale.max(k[i * u0.len() + j].abs());
        }
    }
    err / scale
}

#[test]
fn st_venant_kirchhoff_tangent_is_consistent() {
    for (etype, xyz, cdim, coord) in [
        (303, TRI.to_vec(), 2, CoordSystem::Cartesian),
        (404, QUAD.to_vec(), 2, CoordSystem::Cartesian),
        (404, QUAD.to_vec(), 2, CoordSystem::AxiSymmetric),
        (808, HEX.to_vec(), 3, CoordSystem::Cartesian),
    ] {
        let n = xyz.len();
        let e = tangent_error(&|d| bulk(etype, &xyz, cdim, coord, true, d), n, cdim);
        assert!(e < 1e-6, "{etype} {coord:?}: tangent error {e:e}");
    }
}

#[test]
fn neo_hookean_tangent_is_consistent() {
    for (etype, xyz, cdim, coord) in [
        (303, TRI.to_vec(), 2, CoordSystem::Cartesian),
        (404, QUAD.to_vec(), 2, CoordSystem::Cartesian),
        (404, QUAD.to_vec(), 2, CoordSystem::AxiSymmetric),
        (808, HEX.to_vec(), 3, CoordSystem::Cartesian),
    ] {
        let n = xyz.len();
        let e = tangent_error(&|d| neo(etype, &xyz, cdim, coord, d), n, cdim);
        assert!(e < 1e-6, "{etype} {coord:?}: tangent error {e:e}");
    }
}

#[test]
fn small_strain_step_is_independent_of_the_iterate() {
    // Without Large Deflection the Newton form solves K u = f whatever u it starts from.
    let n = QUAD.len();
    let zero = vec![[0.0; 3]; n];
    let a = bulk(404, &QUAD, 2, CoordSystem::Cartesian, false, &zero);
    let b = bulk(404, &QUAD, 2, CoordSystem::Cartesian, false, &displacement(n, 2));
    let scale = a.stiff.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (x, y) in a.stiff.iter().zip(&b.stiff) {
        assert!((x - y).abs() < 1e-12 * scale);
    }
    for (x, y) in a.force.iter().zip(&b.force) {
        assert!((x - y).abs() < 1e-10 * scale);
    }
}

#[test]
fn small_strain_stiffness_is_the_linear_one() {
    use elmer_physics_elastic::{stress_bulk_local, StressBulkInput};
    for (etype, xyz, cdim) in [(303, TRI.to_vec(), 2), (404, QUAD.to_vec(), 2), (808, HEX.to_vec(), 3)] {
        let n = xyz.len();
        let law = iso(n);
        let (density, load) = (vec![1.0; n], vec![[0.0; 3]; n]);
        let inp = StressBulkInput {
            etype: ElementType(etype),
            xyz: &xyz,
            cdim,
            coord: CoordSystem::Cartesian,
            plane_stress: false,
            law: &law,
            density: &density,
            damping: None,
            load: &load,
            pressure_load: None,
            stress_load: None,
            strain_load: None,
            temperature: None,
            heat_expansion: None,
            want_mass: false,
            np: None,
            ip: None,
        };
        let mut lin = LocalSystem::zeros(0);
        stress_bulk_local(&P1Basis, &Gauss, &inp, &mut lin).unwrap();
        let mine = bulk(etype, &xyz, cdim, CoordSystem::Cartesian, false, &vec![[0.0; 3]; n]);
        let scale = lin.stiff.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for (x, y) in lin.stiff.iter().zip(&mine.stiff) {
            assert!((x - y).abs() < 1e-12 * scale, "{etype}: {x} vs {y}");
        }
    }
}
