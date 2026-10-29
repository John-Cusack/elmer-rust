//! Globally assembled free-DOF matrix and right-hand side against the
//! `Linear System Save` dump of pristine ElmerSolver, at 1e-13 relative.
//!
//! Each `tests/data/<case>/` holds the mesh ElmerGrid wrote, the sif that was
//! run, and `linsys_{a,b,perm}.dat` from the Elmer run. The sif is mirrored by
//! the `Case` below. Elmer applies Dirichlet conditions before the dump (row
//! zeroed, diagonal kept, free rows keep their constrained-column couplings
//! and their load), so only the free rows are compared. The test also copes
//! with a dump that folds constrained values into the right-hand side.

mod support;

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use elmer_core::{ElementType, LocalSystem};
use elmer_physics_elastic::*;
use support::{elem_center, Gauss, P1Basis};

struct Mesh {
    nodes: Vec<[f64; 3]>,
    /// (etype, body, nodes 0-based)
    elems: Vec<(u16, Vec<usize>)>,
    /// (tag, parent, etype, nodes 0-based)
    bnd: Vec<(i32, usize, u16, Vec<usize>)>,
}

fn data(case: &str, file: &str) -> String {
    let p: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "data", case, file].iter().collect();
    fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn read_mesh(case: &str) -> Mesh {
    let nodes = data(case, "mesh.nodes")
        .lines()
        .map(|l| {
            let t: Vec<f64> = l.split_whitespace().skip(2).map(|x| x.parse().unwrap()).collect();
            [t[0], t[1], t[2]]
        })
        .collect();
    let elems = data(case, "mesh.elements")
        .lines()
        .map(|l| {
            let t: Vec<usize> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            (t[2] as u16, t[3..].iter().map(|n| n - 1).collect())
        })
        .collect();
    let bnd = data(case, "mesh.boundary")
        .lines()
        .map(|l| {
            let t: Vec<usize> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            // id tag left right type nodes...; Elmer takes Left, else Right.
            let parent = if t[2] != 0 { t[2] } else { t[3] } - 1;
            (t[1] as i32, parent, t[4] as u16, t[5..].iter().map(|n| n - 1).collect())
        })
        .collect();
    Mesh { nodes, elems, bnd }
}

#[derive(Default, Clone)]
struct Bc {
    traction: [f64; 3],
    normal_traction: f64,
    spring_normal: Option<f64>,
    spring_tensor: Option<[[f64; 3]; 3]>,
}

struct Case {
    name: &'static str,
    cdim: usize,
    coord: CoordSystem,
    plane_stress: bool,
    law: fn(usize) -> ElasticLaw,
    load: [f64; 3],
    pressure: f64,
    /// boundary tag -> prescribed displacement
    dirichlet: Vec<(i32, [f64; 3])>,
    bcs: Vec<(i32, Bc)>,
}

fn iso(e: f64, nu: f64) -> impl Fn(usize) -> ElasticLaw {
    move |n| ElasticLaw::Isotropic { young: vec![e; n], poisson: vec![nu; n] }
}

fn run(c: &Case) {
    let mesh = read_mesh(c.name);
    let nn = mesh.nodes.len();
    let ndof = nn * c.cdim;
    let mut k = vec![0.0; ndof * ndof];
    let mut f = vec![0.0; ndof];

    let scatter = |k: &mut Vec<f64>, f: &mut Vec<f64>, nodes: &[usize], l: &LocalSystem| {
        for (p, &gp) in nodes.iter().enumerate() {
            for i in 0..c.cdim {
                f[gp * c.cdim + i] += l.force[p * c.cdim + i];
                for (q, &gq) in nodes.iter().enumerate() {
                    for j in 0..c.cdim {
                        k[(gp * c.cdim + i) * ndof + gq * c.cdim + j] +=
                            l.stiff[(p * c.cdim + i) * l.n + q * c.cdim + j];
                    }
                }
            }
        }
    };

    for (et, nodes) in &mesh.elems {
        let n = nodes.len();
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&i| mesh.nodes[i]).collect();
        let law = (c.law)(n);
        let load = vec![c.load; n];
        let dens = vec![1.0; n];
        let pr = vec![c.pressure; n];
        let inp = StressBulkInput {
            etype: ElementType(*et),
            xyz: &xyz,
            cdim: c.cdim,
            coord: c.coord,
            plane_stress: c.plane_stress,
            law: &law,
            density: &dens,
            damping: None,
            load: &load,
            pressure_load: if c.pressure != 0.0 { Some(&pr) } else { None },
            stress_load: None,
            strain_load: None,
            temperature: None,
            heat_expansion: None,
            want_mass: false,
            np: None,
        };
        let mut out = LocalSystem::zeros(0);
        stress_bulk_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
        scatter(&mut k, &mut f, nodes, &out);
    }

    let mut constrained: HashMap<usize, f64> = HashMap::new();
    for (tag, bnodes, et, parent) in mesh.bnd.iter().map(|(t, p, e, n)| (*t, n, *e, *p)) {
        for (dt, val) in &c.dirichlet {
            if *dt == tag {
                for &nd in bnodes {
                    for i in 0..c.cdim {
                        constrained.insert(nd * c.cdim + i, val[i]);
                    }
                }
            }
        }
        let Some((_, bc)) = c.bcs.iter().find(|(t, _)| *t == tag) else { continue };
        let n = bnodes.len();
        let xyz: Vec<[f64; 3]> = bnodes.iter().map(|&i| mesh.nodes[i]).collect();
        let pxyz: Vec<[f64; 3]> = mesh.elems[parent].1.iter().map(|&i| mesh.nodes[i]).collect();
        let tr = vec![bc.traction; n];
        let nt = vec![bc.normal_traction; n];
        let sn = vec![bc.spring_normal.unwrap_or(0.0); n];
        let st = vec![bc.spring_tensor.unwrap_or_default(); n];
        let spring = if bc.spring_normal.is_some() {
            Spring::Normal(&sn)
        } else if bc.spring_tensor.is_some() {
            Spring::Tensor(&st)
        } else {
            Spring::None
        };
        let inp = StressBoundaryInput {
            etype: ElementType(et),
            xyz: &xyz,
            cdim: c.cdim,
            coord: c.coord,
            traction: &tr,
            normal_traction: &nt,
            spring,
            parent_center: elem_center(&pxyz),
            np: None,
        };
        let mut out = LocalSystem::zeros(0);
        stress_boundary_local(&P1Basis, &Gauss, &inp, &mut out).unwrap();
        scatter(&mut k, &mut f, bnodes, &out);
    }
    assert!(!constrained.is_empty());

    // Elmer's numbering.
    let perm: Vec<usize> = data(c.name, "linsys_perm.dat")
        .lines()
        .map(|l| l.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect();
    let row_of = |d: usize| c.cdim * (perm[d / c.cdim] - 1) + d % c.cdim;
    let mut dump: HashMap<(usize, usize), f64> = HashMap::new();
    for l in data(c.name, "linsys_a.dat").lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        dump.insert((t[0].parse::<usize>().unwrap() - 1, t[1].parse::<usize>().unwrap() - 1), t[2].parse().unwrap());
    }
    let mut rhs: Vec<f64> = vec![0.0; ndof];
    for l in data(c.name, "linsys_b.dat").lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        rhs[t[0].parse::<usize>().unwrap() - 1] = t[1].parse().unwrap();
    }

    let kmax = k.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    let free: Vec<usize> = (0..ndof).filter(|d| !constrained.contains_key(d)).collect();
    // Does the dump keep the free-row/constrained-column couplings? If so the
    // right-hand side is the plain load, otherwise the couplings are folded in.
    let keeps_columns = free.iter().any(|&i| {
        constrained.keys().any(|&j| dump.get(&(row_of(i), row_of(j))).is_some_and(|v| *v != 0.0))
    });
    let fmax = free.iter().fold(0.0f64, |a: f64, &d| a.max(f[d].abs()).max(rhs[row_of(d)].abs()));
    let (mut worst_k, mut worst_f) = (0.0f64, 0.0f64);
    let nfree = free.len();
    for &i in &free {
        let folded: f64 = if keeps_columns {
            0.0
        } else {
            constrained.iter().map(|(&j, v)| k[i * ndof + j] * v).sum()
        };
        worst_f = worst_f.max((rhs[row_of(i)] - (f[i] - folded)).abs() / fmax);
        let cols = if keeps_columns { (0..ndof).collect::<Vec<_>>() } else { free.clone() };
        for j in cols {
            let theirs = dump.get(&(row_of(i), row_of(j))).copied().unwrap_or(0.0);
            worst_k = worst_k.max((theirs - k[i * ndof + j]).abs() / kmax);
        }
    }
    println!("{}: keeps constrained columns: {keeps_columns}", c.name);
    println!("{}: {nfree} free dofs, matrix {worst_k:.2e}, rhs {worst_f:.2e}", c.name);
    assert!(worst_k < 1e-13, "{}: matrix mismatch {worst_k:e}", c.name);
    assert!(worst_f < 1e-13, "{}: rhs mismatch {worst_f:e}", c.name);
}

fn fixed() -> (i32, [f64; 3]) {
    (1, [0.0; 3])
}

fn beam_bcs() -> Vec<(i32, Bc)> {
    vec![(2, Bc { traction: [0.0, -5.0e7, 0.0], normal_traction: 1.0e6, ..Default::default() })]
}

#[test]
fn a2d_plane_strain_triangles_via_stresssolve_front() {
    run(&Case {
        name: "a2d",
        cdim: 2,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law: |n| iso(70e9, 0.3)(n),
        load: [1.0e3, -2.0e3, 0.0],
        pressure: 0.0,
        dirichlet: vec![fixed()],
        bcs: beam_bcs(),
    });
}

fn b2d_bcs() -> Vec<(i32, Bc)> {
    let mut b = beam_bcs();
    b.push((3, Bc { spring_normal: Some(2.0e9), normal_traction: 3.0e5, ..Default::default() }));
    let mut t = [[0.0; 3]; 3];
    t[0][0] = 1.0e9;
    t[1][1] = 4.0e9;
    t[0][1] = 5.0e8;
    b.push((4, Bc { spring_tensor: Some(t), traction: [2.0e5, 0.0, 0.0], ..Default::default() }));
    b
}

#[test]
fn b2d_plane_stress_quads_springs_pressure() {
    run(&Case {
        name: "b2d",
        cdim: 2,
        coord: CoordSystem::Cartesian,
        plane_stress: true,
        law: |n| iso(70e9, 0.3)(n),
        load: [1.0e3, -2.0e3, 0.0],
        pressure: 5.0e2,
        dirichlet: vec![fixed()],
        bcs: b2d_bcs(),
    });
}

#[test]
fn e2d_axisymmetric_quads() {
    run(&Case {
        name: "e2d",
        cdim: 2,
        coord: CoordSystem::AxiSymmetric,
        plane_stress: false,
        law: |n| iso(70e9, 0.3)(n),
        load: [1.0e3, -2.0e3, 0.0],
        pressure: 5.0e2,
        dirichlet: vec![fixed()],
        bcs: b2d_bcs(),
    });
}

#[test]
fn c3d_anisotropic_hexahedra() {
    run(&Case {
        name: "c3d",
        cdim: 3,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law: |n| ElasticLaw::Anisotropic {
            c: vec![
                [
                    [1346.154, 576.923, 576.923, 0.0, 0.0, 0.0],
                    [576.923, 1346.154, 576.923, 0.0, 0.0, 0.0],
                    [576.923, 576.923, 1346.154, 200.0, 0.0, 0.0],
                    [0.0, 0.0, 200.0, 384.615, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 384.615, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 384.615],
                ];
                n
            ],
            rotation: None,
        },
        load: [0.0; 3],
        pressure: 0.0,
        dirichlet: vec![(5, [0.0; 3]), (6, [1.0e-3, 2.0e-3, 0.5e-3])],
        bcs: vec![],
    });
}

#[test]
fn d3d_isotropic_hexahedra_loads_and_spring() {
    run(&Case {
        name: "d3d",
        cdim: 3,
        coord: CoordSystem::Cartesian,
        plane_stress: false,
        law: |n| iso(2.0e5, 0.33)(n),
        load: [10.0, -5.0, 2.0],
        pressure: 0.0,
        dirichlet: vec![(5, [0.0; 3])],
        bcs: vec![
            (6, Bc { traction: [100.0, 0.0, -40.0], normal_traction: 25.0, ..Default::default() }),
            (1, Bc { spring_normal: Some(300.0), ..Default::default() }),
        ],
    });
}
