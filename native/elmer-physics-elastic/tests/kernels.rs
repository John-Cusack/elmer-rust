//! Local-kernel physics checks on small meshes with a test-local P1 basis.

mod support;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_elastic::*;
use support::{elem_center, scatter, Gauss, P1Basis};

const TOL: f64 = 1e-12;

struct Mesh {
    nodes: Vec<[f64; 3]>,
    elems: Vec<Vec<usize>>,
    etype: u16,
}

/// A 2x2 block of unit squares, interior node (1,1) pushed off-centre, split
/// into triangles (or kept as quads).
fn patch2d(quads: bool) -> Mesh {
    let mut nodes = vec![];
    for j in 0..3 {
        for i in 0..3 {
            nodes.push([i as f64, j as f64, 0.0]);
        }
    }
    nodes[4] = [1.2, 0.85, 0.0];
    let id = |i: usize, j: usize| 3 * j + i;
    let mut elems = vec![];
    for j in 0..2 {
        for i in 0..2 {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            if quads {
                elems.push(vec![a, b, c, d]);
            } else {
                elems.push(vec![a, b, c]);
                elems.push(vec![a, c, d]);
            }
        }
    }
    Mesh { nodes, elems, etype: if quads { 404 } else { 303 } }
}

/// Big tetrahedron split into four around an interior point.
fn patch3d() -> Mesh {
    let nodes = vec![
        [0.0, 0.0, 0.0],
        [2.0, 0.1, 0.0],
        [0.2, 1.8, 0.0],
        [0.3, 0.2, 2.2],
        [0.55, 0.5, 0.45],
    ];
    let elems = vec![
        vec![0, 1, 2, 4],
        vec![0, 1, 4, 3],
        vec![0, 4, 2, 3],
        vec![4, 1, 2, 3],
    ];
    Mesh { nodes, elems, etype: 504 }
}

struct Setup {
    law: ElasticLaw,
    cdim: usize,
    coord: CoordSystem,
    plane_stress: bool,
}

fn iso(e: f64, nu: f64, n: usize) -> ElasticLaw {
    ElasticLaw::Isotropic { young: vec![e; n], poisson: vec![nu; n] }
}

#[derive(Default)]
struct Opt<'a> {
    load: Option<&'a [[f64; 3]]>,
    temperature: Option<&'a [f64]>,
    heat_expansion: Option<&'a [[f64; 3]]>,
    strain_load: Option<&'a [[f64; 6]]>,
    want_mass: bool,
}

fn bulk(m: &Mesh, e: usize, s: &Setup, o: &Opt<'_>) -> LocalSystem {
    let n = m.elems[e].len();
    let xyz: Vec<[f64; 3]> = m.elems[e].iter().map(|&k| m.nodes[k]).collect();
    let zeros3 = vec![[0.0; 3]; n];
    let dens = vec![1.0; n];
    let inp = StressBulkInput {
        etype: ElementType(m.etype),
        xyz: &xyz,
        cdim: s.cdim,
        coord: s.coord,
        plane_stress: s.plane_stress,
        law: &s.law,
        density: &dens,
        damping: None,
        load: o.load.unwrap_or(&zeros3),
        pressure_load: None,
        stress_load: None,
        strain_load: o.strain_load,
        temperature: o.temperature,
        heat_expansion: o.heat_expansion,
        want_mass: o.want_mass,
        np: None,
    };
    let mut out = LocalSystem::zeros(0);
    stress_bulk_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
    out
}

fn assemble(m: &Mesh, s: &Setup) -> (Vec<f64>, Vec<f64>) {
    let ng = m.nodes.len() * s.cdim;
    let (mut k, mut f) = (vec![0.0; ng * ng], vec![0.0; ng]);
    for e in 0..m.elems.len() {
        let loc = bulk(m, e, s, &Opt::default());
        scatter(&mut k, &mut f, ng, s.cdim, &m.elems[e], &loc);
    }
    (k, f)
}

fn matvec(k: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n).map(|i| (0..n).map(|j| k[i * n + j] * x[j]).sum()).collect()
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |a, x| a.max(x.abs()))
}

#[test]
fn constant_strain_patch_2d_plane_strain_and_stress() {
    for quads in [false, true] {
        for plane_stress in [false, true] {
            let m = patch2d(quads);
            let s = Setup { law: iso(210.0, 0.3, if quads { 4 } else { 3 }), cdim: 2, coord: CoordSystem::Cartesian, plane_stress };
            let (k, _) = assemble(&m, &s);
            // u = A x + b, a general constant-strain field with rotation.
            let a = [[1e-3, 4e-4], [-2e-4, 7e-4]];
            let u: Vec<f64> = m
                .nodes
                .iter()
                .flat_map(|x| (0..2).map(move |i| a[i][0] * x[0] + a[i][1] * x[1] + 0.01 * (i as f64 + 1.0)))
                .collect();
            let r = matvec(&k, &u);
            // Interior node 4: no internal force without loads.
            assert!(r[8].abs() < TOL && r[9].abs() < TOL, "quads={quads} ps={plane_stress} r={:?}", &r[8..10]);
        }
    }
}

#[test]
fn constant_strain_patch_3d() {
    let m = patch3d();
    let s = Setup { law: iso(70.0, 0.33, 4), cdim: 3, coord: CoordSystem::Cartesian, plane_stress: false };
    let (k, _) = assemble(&m, &s);
    let a = [[1e-3, 2e-4, -3e-4], [5e-4, -7e-4, 1e-4], [-2e-4, 3e-4, 9e-4]];
    let u: Vec<f64> = m
        .nodes
        .iter()
        .flat_map(|x| (0..3).map(move |i| (0..3).map(|j| a[i][j] * x[j]).sum::<f64>() + 0.02))
        .collect();
    let r = matvec(&k, &u);
    assert!(max_abs(&r[12..15]) < TOL, "{:?}", &r[12..15]);
}

#[test]
fn rigid_motion_has_no_energy_before_constraints() {
    let m = patch3d();
    let s = Setup { law: iso(70.0, 0.33, 4), cdim: 3, coord: CoordSystem::Cartesian, plane_stress: false };
    let (k, _) = assemble(&m, &s);
    let n = m.nodes.len();
    // Translation along each axis and infinitesimal rotation about z.
    for c in 0..3 {
        let u: Vec<f64> = (0..n * 3).map(|d| if d % 3 == c { 1.0 } else { 0.0 }).collect();
        assert!(max_abs(&matvec(&k, &u)) < TOL);
    }
    let rot: Vec<f64> = m.nodes.iter().flat_map(|x| [-x[1], x[0], 0.0]).collect();
    assert!(max_abs(&matvec(&k, &rot)) < TOL);
    // Symmetry.
    let ng = n * 3;
    for i in 0..ng {
        for j in 0..ng {
            assert!((k[i * ng + j] - k[j * ng + i]).abs() < TOL * 100.0);
        }
    }
}

#[test]
fn plane_stress_equals_plane_strain_with_modified_constants() {
    let m = patch2d(false);
    let (e, nu) = (210.0, 0.3);
    let ps = Setup { law: iso(e, nu, 3), cdim: 2, coord: CoordSystem::Cartesian, plane_stress: true };
    let pe = Setup {
        law: iso(e * (1.0 + 2.0 * nu) / ((1.0 + nu) * (1.0 + nu)), nu / (1.0 + nu), 3),
        cdim: 2,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
    };
    let same = Setup { law: iso(e, nu, 3), cdim: 2, coord: CoordSystem::Cartesian, plane_stress: false };
    let (kps, _) = assemble(&m, &ps);
    let (kpe, _) = assemble(&m, &pe);
    let (ksame, _) = assemble(&m, &same);
    let d1 = kps.iter().zip(&kpe).fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
    let d2 = kps.iter().zip(&ksame).fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
    assert!(d1 < 1e-11, "{d1}");
    assert!(d2 > 1.0, "plane stress and strain must differ for the same E, nu: {d2}");
}

#[test]
fn anisotropic_c_reproduces_isotropic_in_3d_and_2d() {
    let (e, nu) = (70.0, 0.33);
    let l1 = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        c[i][..3].fill(l1);
        c[i][i] = l1 + 2.0 * mu;
        c[i + 3][i + 3] = mu;
    }
    let ani = |n| ElasticLaw::Anisotropic { c: vec![c; n], rotation: None };
    let m3 = patch3d();
    for (m, cdim, plane_stress, n) in [(m3, 3, false, 4), (patch2d(false), 2, false, 3), (patch2d(false), 2, true, 3)] {
        let a = Setup { law: ani(n), cdim, coord: CoordSystem::Cartesian, plane_stress };
        let i = Setup { law: iso(e, nu, n), cdim, coord: CoordSystem::Cartesian, plane_stress };
        let (ka, _) = assemble(&m, &a);
        let (ki, _) = assemble(&m, &i);
        let d = ka.iter().zip(&ki).fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
        assert!(d < 1e-11, "cdim {cdim} ps {plane_stress}: {d}");
    }
}

#[test]
fn orthotropic_rotation_by_quarter_turn_swaps_axes() {
    // Diagonal-ish orthotropic C with distinct E_x, E_y; rotating the material
    // frame by 90 degrees about z equals swapping the in-plane moduli.
    let mut c = [[0.0; 6]; 6];
    let (cx, cy, cz, cxy) = (200.0, 100.0, 80.0, 30.0);
    c[0][0] = cx;
    c[1][1] = cy;
    c[2][2] = cz;
    c[0][1] = 20.0;
    c[1][0] = 20.0;
    c[0][2] = 15.0;
    c[2][0] = 15.0;
    c[1][2] = 12.0;
    c[2][1] = 12.0;
    c[3][3] = cxy;
    c[4][4] = 25.0;
    c[5][5] = 22.0;
    let mut swapped = c;
    swapped[0][0] = cy;
    swapped[1][1] = cx;
    swapped[0][2] = 12.0;
    swapped[2][0] = 12.0;
    swapped[1][2] = 15.0;
    swapped[2][1] = 15.0;
    swapped[4][4] = 22.0;
    swapped[5][5] = 25.0;
    let t = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let m = patch3d();
    let rot = Setup { law: ElasticLaw::Anisotropic { c: vec![c; 4], rotation: Some(t) }, cdim: 3, coord: CoordSystem::Cartesian, plane_stress: false };
    let sw = Setup { law: ElasticLaw::Anisotropic { c: vec![swapped; 4], rotation: None }, cdim: 3, coord: CoordSystem::Cartesian, plane_stress: false };
    let (kr, _) = assemble(&m, &rot);
    let (ks, _) = assemble(&m, &sw);
    let d = kr.iter().zip(&ks).fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
    assert!(d < 1e-11, "{d}");
}

#[test]
fn axisymmetric_patch_is_exact() {
    // u_r = a r, u_z = c z: constant (rr, zz, hoop) strain, no net force on the
    // interior node.
    let m = patch2d(true);
    let s = Setup { law: iso(100.0, 0.25, 4), cdim: 2, coord: CoordSystem::AxiSymmetric, plane_stress: false };
    // Shift the mesh off the axis so r > 0 everywhere.
    let m = Mesh { nodes: m.nodes.iter().map(|x| [x[0] + 1.0, x[1], 0.0]).collect(), ..m };
    let (k, _) = assemble(&m, &s);
    let (a, c) = (1e-3, -4e-4);
    let u: Vec<f64> = m.nodes.iter().flat_map(|x| [a * x[0], c * x[1]]).collect();
    let r = matvec(&k, &u);
    assert!(r[8].abs() < 1e-10 && r[9].abs() < 1e-10, "{:?}", &r[8..10]);
    // Compare with the closed form for the force on a boundary node region:
    // total radial reaction is not zero, so the matrix is not trivially null.
    assert!(max_abs(&r) > 1e-3);
}

#[test]
fn body_force_and_thermal_load_balance() {
    let m = patch3d();
    let s = Setup { law: iso(70.0, 0.33, 4), cdim: 3, coord: CoordSystem::Cartesian, plane_stress: false };
    // Body force: total load equals f * volume.
    let vol: f64 = (0..4)
        .map(|e| {
            let x: Vec<_> = m.elems[e].iter().map(|&k| m.nodes[k]).collect();
            let d = |a: usize| [x[a][0] - x[0][0], x[a][1] - x[0][1], x[a][2] - x[0][2]];
            let (a, b, c) = (d(1), d(2), d(3));
            ((a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0).abs()
        })
        .sum();
    let mut tot = [0.0; 3];
    for e in 0..4 {
        let l = bulk(&m, e, &s, &Opt { load: Some(&[[1.0, -2.0, 3.0]; 4]), ..Default::default() });
        for p in 0..4 {
            for c in 0..3 {
                tot[c] += l.force[p * 3 + c];
            }
        }
    }
    assert!((tot[0] - vol).abs() < TOL && (tot[1] + 2.0 * vol).abs() < TOL && (tot[2] - 3.0 * vol).abs() < TOL);

    // Free thermal expansion u = beta T x leaves the interior node in
    // equilibrium (stress-free), so K u = f there.
    let (beta, temp) = (1.2e-5, 100.0);
    let u: Vec<f64> = m.nodes.iter().flat_map(|x| x.map(|c| beta * temp * c)).collect();
    let ng = 15;
    let (mut k, mut f) = (vec![0.0; ng * ng], vec![0.0; ng]);
    for e in 0..4 {
        let l = bulk(
            &m,
            e,
            &s,
            &Opt { temperature: Some(&[temp; 4]), heat_expansion: Some(&[[beta; 3]; 4]), ..Default::default() },
        );
        scatter(&mut k, &mut f, ng, 3, &m.elems[e], &l);
    }
    let ku = matvec(&k, &u);
    for d in 12..15 {
        assert!((ku[d] - f[d]).abs() < 1e-9, "dof {d}: {} vs {}", ku[d], f[d]);
    }
    assert!(max_abs(&f) > 1e-3, "thermal load must be nonzero");

    // Strain load is the same thing written directly.
    let eps = beta * temp;
    let mut f2 = vec![0.0; ng];
    let mut k2 = vec![0.0; ng * ng];
    for e in 0..4 {
        let l = bulk(&m, e, &s, &Opt { strain_load: Some(&[[eps, eps, eps, 0.0, 0.0, 0.0]; 4]), ..Default::default() });
        scatter(&mut k2, &mut f2, ng, 3, &m.elems[e], &l);
    }
    for d in 0..ng {
        assert!((f[d] - f2[d]).abs() < 1e-11);
    }
}

#[test]
fn mass_matrix_totals_rho_volume() {
    let m = patch2d(true);
    let s = Setup { law: iso(1.0, 0.3, 4), cdim: 2, coord: CoordSystem::Cartesian, plane_stress: false };
    let mut tot = 0.0;
    for e in 0..4 {
        let l = bulk(&m, e, &s, &Opt { want_mass: true, ..Default::default() });
        tot += l.mass.unwrap().iter().sum::<f64>();
    }
    // rho = 1, area 4, two components.
    assert!((tot - 8.0).abs() < 1e-12, "{tot}");
}

#[derive(Default)]
struct EdgeOpt<'a> {
    traction: Option<&'a [[f64; 3]; 2]>,
    normal_traction: Option<&'a [f64; 2]>,
    spring: Option<Spring<'a>>,
}

fn edge_bc(p: [[f64; 3]; 2], center: [f64; 3], o: EdgeOpt<'_>) -> LocalSystem {
    let tr = [[0.0; 3]; 2];
    let nt = [0.0; 2];
    let inp = StressBoundaryInput {
        etype: ElementType(202),
        xyz: &p,
        cdim: 2,
        coord: CoordSystem::Cartesian,
        traction: o.traction.unwrap_or(&tr),
        normal_traction: o.normal_traction.unwrap_or(&nt),
        spring: o.spring.unwrap_or(Spring::None),
        parent_center: center,
        np: None,
    };
    let mut out = LocalSystem::zeros(0);
    stress_boundary_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
    out
}

#[test]
fn boundary_traction_pressure_and_orientation() {
    // Edge along +x at y=0 of a body above it: outward normal is -y.
    let p = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let above = [1.0, 1.0, 0.0];
    let below = [1.0, -1.0, 0.0];
    let l = edge_bc(p, above, EdgeOpt { normal_traction: Some(&[3.0, 3.0]), ..Default::default() });
    let fy: f64 = (0..2).map(|k| l.force[k * 2 + 1]).sum();
    assert!((fy + 6.0).abs() < TOL, "pressure 3 over length 2 pushes along -y (outward): {fy}");
    let l = edge_bc(p, below, EdgeOpt { normal_traction: Some(&[3.0, 3.0]), ..Default::default() });
    let fy: f64 = (0..2).map(|k| l.force[k * 2 + 1]).sum();
    assert!((fy - 6.0).abs() < TOL);

    let l = edge_bc(p, above, EdgeOpt { traction: Some(&[[1.0, 2.0, 0.0], [1.0, 2.0, 0.0]]), ..Default::default() });
    let (fx, fy): (f64, f64) = (0..2).fold((0.0, 0.0), |a, k| (a.0 + l.force[k * 2], a.1 + l.force[k * 2 + 1]));
    assert!((fx - 2.0).abs() < TOL && (fy - 4.0).abs() < TOL);
    // Consistent nodal load of a constant traction is half each.
    assert!((l.force[0] - 1.0).abs() < TOL && (l.force[2] - 1.0).abs() < TOL);
}

#[test]
fn boundary_normal_traction_on_closed_triangle_is_in_equilibrium() {
    let tri = [[0.0, 0.0, 0.0], [2.0, 0.1, 0.0], [0.3, 1.5, 0.0]];
    let center = elem_center(&tri);
    let mut tot = [0.0; 2];
    for (a, b) in [(0, 1), (1, 2), (2, 0)] {
        let l = edge_bc([tri[a], tri[b]], center, EdgeOpt { normal_traction: Some(&[5.0, 5.0]), ..Default::default() });
        for k in 0..2 {
            tot[0] += l.force[k * 2];
            tot[1] += l.force[k * 2 + 1];
        }
    }
    assert!(tot[0].abs() < TOL && tot[1].abs() < TOL, "{tot:?}");
}

#[test]
fn spring_boundary_matrices() {
    let p = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let center = [1.0, 1.0, 0.0];
    // Normal spring k: u = n (unit, outward) over the edge stores k * length.
    let l = edge_bc(p, center, EdgeOpt { spring: Some(Spring::Normal(&[4.0, 4.0])), ..Default::default() });
    let nvec = [0.0, -1.0];
    let u = [nvec[0], nvec[1], nvec[0], nvec[1]];
    let ku = (0..4).map(|i| (0..4).map(|j| l.stiff[i * 4 + j] * u[j]).sum::<f64>());
    let energy: f64 = ku.zip(u).map(|(a, b)| a * b).sum();
    assert!((energy - 4.0 * 2.0).abs() < TOL, "{energy}");
    // Component springs act per direction.
    let mut k = [[0.0; 3]; 3];
    k[0][0] = 2.0;
    k[1][1] = 6.0;
    let kt = [k; 2];
    let l = edge_bc(p, center, EdgeOpt { spring: Some(Spring::Tensor(&kt)), ..Default::default() });
    let ux = [1.0, 0.0, 1.0, 0.0];
    let ex: f64 = (0..4).map(|i| (0..4).map(|j| ux[i] * l.stiff[i * 4 + j] * ux[j]).sum::<f64>()).sum();
    assert!((ex - 2.0 * 2.0).abs() < TOL, "{ex}");
}

#[test]
fn unsupported_shapes_are_refused() {
    let m = patch2d(false);
    let s = Setup { law: iso(1.0, 0.3, 3), cdim: 2, coord: CoordSystem::Cartesian, plane_stress: false };
    let xyz: Vec<[f64; 3]> = m.elems[0].iter().map(|&k| m.nodes[k]).collect();
    let z = vec![[0.0; 3]; 3];
    let d = vec![1.0; 3];
    let inp = StressBulkInput {
        etype: ElementType(303),
        xyz: &xyz,
        cdim: 1,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law: &s.law,
        density: &d,
        damping: None,
        load: &z,
        pressure_load: None,
        stress_load: None,
        strain_load: None,
        temperature: None,
        heat_expansion: None,
        want_mass: false,
        np: None,
    };
    let mut out = LocalSystem::zeros(0);
    assert!(matches!(
        stress_bulk_local(&P1Basis, &Gauss, &inp, &mut out),
        Err(StressError::Unsupported(_))
    ));
}
