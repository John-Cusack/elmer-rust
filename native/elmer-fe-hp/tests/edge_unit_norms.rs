//! Reproduces the `ElementUnitTest_edge`, `ElementUnitTest_edge_piola` and
//! `ElementUnitTest_edge_quadratic` quick tests of the Elmer test suite.
//!
//! Each test integrates the mass matrix (sum of all entries) and the
//! stiffness matrix (sum of absolute entries) of the edge basis over one
//! unit element, exactly as `SingleEdgeElements.F90` does, with the element's
//! default Gauss rule (elements.def `Gauss Points`: 3, 4, 4, 8, 8, 8 for
//! 303, 404, 504, 605, 706, 808). The geometry and node numbering are those of
//! the test's `unitelems` mesh (bulk elements first, then the 2-D boundary
//! elements 404 and 303 on the face z = 0 of the brick).
//!
//! The reference values are the 17-digit `Norm for element` lines printed by
//! the pristine ElmerSolver (the `.sif` files only carry ~8 digits). The Rust
//! results agree to at most ~2e-15 relative and many are bitwise equal; they
//! are not asserted bitwise because the Gauss points here are closed forms
//! (Elmer computes them by an eigenvalue iteration) and summation order inside
//! the intrinsics may differ. Only these aggregate norms were captured, not
//! per-function values, so pointwise agreement is not asserted.

// Reference literals are printed with 17 digits to match the Fortran output.
#![allow(clippy::excessive_precision)]

use elmer_core::mesh::ElementType;
use elmer_fe_hp::edge::{edge_basis, edge_element_info, get_edge_basis, EdgeEval};

type P3 = [f64; 3];

struct Unit {
    code: u16,
    nodes: Vec<P3>,
    gindex: Vec<i64>,
}

fn units() -> Vec<Unit> {
    let brick = vec![
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [1., 0., 1.],
        [1., 1., 1.],
        [0., 1., 1.],
    ];
    vec![
        Unit { code: 808, nodes: brick.clone(), gindex: (1..=8).collect() },
        Unit {
            code: 706,
            nodes: vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [1., 1., 0.],
                [0., 0., 1.],
                [1., 0., 1.],
                [1., 1., 1.],
            ],
            gindex: (9..=14).collect(),
        },
        Unit {
            code: 605,
            nodes: vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.]],
            gindex: (15..=19).collect(),
        },
        Unit {
            code: 504,
            nodes: vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 0., 1.]],
            gindex: (20..=23).collect(),
        },
        Unit { code: 404, nodes: brick[..4].to_vec(), gindex: (1..=4).collect() },
        Unit { code: 303, nodes: brick[..3].to_vec(), gindex: (1..=3).collect() },
    ]
}

/// `GaussPointsAdapt(Element)` for the default rules: `(u, v, w, weight)`.
fn gauss(code: u16) -> Vec<[f64; 4]> {
    let g = 1.0 / 3.0f64.sqrt();
    let pts = [-g, g];
    let wts = [1.0, 1.0];
    let mut out = Vec::new();
    match code {
        303 => {
            let u = [0.16666666666666667, 0.66666666666666667, 0.16666666666666667];
            let v = [0.16666666666666667, 0.16666666666666667, 0.66666666666666667];
            for t in 0..3 {
                out.push([u[t], v[t], 0.0, 0.33333333333333333 / 2.0]);
            }
        }
        404 => {
            for i in 0..2 {
                for j in 0..2 {
                    out.push([pts[j], pts[i], 0.0, wts[i] * wts[j]]);
                }
            }
        }
        504 => {
            let u = [0.1757281246520584, 0.2445310270213291, 0.5556470949048655, 0.0240937534217468];
            let v = [0.5656137776620919, 0.0501800797762026, 0.1487681308666864, 0.2354380116950194];
            let w = [0.2180665126782654, 0.5635595064952189, 0.0350112499848832, 0.1833627308416330];
            for t in 0..4 {
                out.push([u[t], v[t], w[t], 0.25 / 6.0]);
            }
        }
        605 => {
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let w = (pts[i] + 1.0) / 2.0;
                        out.push([
                            pts[k] * (1.0 - w),
                            pts[j] * (1.0 - w),
                            w,
                            wts[i] * wts[j] * wts[k] * (1.0 - w) * (1.0 - w) / 2.0,
                        ]);
                    }
                }
            }
        }
        706 => {
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let v = (pts[j] + 1.0) / 2.0;
                        let u = (pts[k] + 1.0) / 2.0 * (1.0 - v);
                        out.push([u, v, pts[i], wts[i] * wts[j] * wts[k] * (1.0 - v) / 4.0]);
                    }
                }
            }
        }
        808 => {
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        out.push([pts[k], pts[j], pts[i], wts[i] * wts[j] * wts[k]]);
                    }
                }
            }
        }
        _ => unreachable!(),
    }
    out
}

#[derive(Clone, Copy)]
enum Kind {
    /// ElementUnitTest_edge: `GetEdgeBasis`
    Whitney,
    /// ElementUnitTest_edge_piola / _quadratic: `EdgeElementInfo`
    Piola(u8),
}

fn dot(a: &P3, b: &P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `(mass norm, stiffness norm)` of one unit element.
fn norms(u: &Unit, kind: Kind) -> (f64, f64) {
    let et = ElementType(u.code);
    let nd = 64;
    let mut mass = vec![vec![0.0f64; nd]; nd];
    let mut stiff = vec![vec![0.0f64; nd]; nd];
    for ip in gauss(u.code) {
        let ev: EdgeEval = match kind {
            Kind::Whitney => get_edge_basis(et, &u.nodes, &u.gindex, [ip[0], ip[1], ip[2]], 3).unwrap(),
            Kind::Piola(d) => edge_element_info(et, &u.nodes, &u.gindex, [ip[0], ip[1], ip[2]], d, 3).unwrap(),
        };
        let weight = ip[3] * ev.det;
        for p in 0..ev.dofs {
            for q in 0..ev.dofs {
                stiff[p][q] += weight * dot(&ev.curl[p], &ev.curl[q]);
                mass[p][q] += weight * dot(&ev.edge[p], &ev.edge[q]);
            }
        }
    }
    let m: f64 = mass.iter().flatten().sum();
    let s: f64 = stiff.iter().flatten().map(|x| x.abs()).sum();
    (m, s)
}

fn check(kind: Kind, expected: [f64; 12]) {
    let us = units();
    for (idx, u) in us.iter().enumerate() {
        let (m, s) = norms(u, kind);
        for (what, got, want) in [("mass", m, expected[idx]), ("stiffness", s, expected[6 + idx])] {
            let rel = (got - want).abs() / want.abs();
            eprintln!("{} {what}: {got:.17} vs {want:.17} rel {rel:.2e}", u.code);
            assert!(rel < 1e-13, "{} {what}: {got} vs {want} (rel {rel:e})", u.code);
        }
    }
}

// Element order: 808, 706, 605, 504, 404, 303; six mass norms then six
// stiffness norms.
#[test]
fn edge_whitney_lowest_order() {
    check(
        Kind::Whitney,
        [
            2.3333333333333326, 0.99999999999999967, 0.89660493827160537, 0.48333333333333328,
            1.3333333333333335, 0.49999999999999994, 39.999999999999993, 31.333333333333343,
            26.666666666666650, 22.666666666666671, 16.000000000000000, 18.000000000000000,
        ],
    );
}

#[test]
fn edge_piola_lowest_order() {
    check(
        Kind::Piola(1),
        [
            2.7777777777777795, 0.32665885066844397, 0.99002266632306113, 0.10387605587898596,
            1.7777777777777781, 0.18232128337103762, 127.11111111111153, 19.762947354836061,
            36.908864303575939, 4.0069384267237691, 50.666666666666657, 5.1961524227066347,
        ],
    );
}

#[test]
fn edge_piola_second_order() {
    check(
        Kind::Piola(2),
        [
            87.037037037037152, 4.9299745098231327, 12.442365435833340, 0.46645189992006059,
            20.222222222222218, 0.94870061348874046, 1152.8888888888928, 207.85023449490080,
            328.93897739881140, 52.492523723186075, 146.66666666666666, 51.712819757647623,
        ],
    );
}

#[test]
fn edge_basis_reference_counts() {
    let mut out = [[0.0; 3]; 64];
    let mut curl = [[0.0; 3]; 64];
    let p = [0.2, 0.3, 0.25];
    for (code, lo, hi, w) in [
        (303u16, 3usize, 8usize, 3usize),
        (404, 6, 12, 4),
        (504, 6, 20, 6),
        (605, 10, 31, 8),
        (706, 15, 36, 9),
        (808, 27, 54, 12),
    ] {
        let et = ElementType(code);
        assert_eq!(edge_basis(et, 1, true, p, &mut out, &mut curl).unwrap(), lo);
        assert_eq!(edge_basis(et, 2, true, p, &mut out, &mut curl).unwrap(), hi);
        assert_eq!(edge_basis(et, 1, false, p, &mut out, &mut curl).unwrap(), w);
    }
    assert!(edge_basis(ElementType(808), 3, true, p, &mut out, &mut curl).is_err());
    assert!(edge_basis(ElementType(808), 2, false, p, &mut out, &mut curl).is_err());
    assert!(edge_basis(ElementType(202), 1, true, p, &mut out, &mut curl).is_err());
    assert!(edge_basis(ElementType(808), 1, true, p, &mut out[..3], &mut curl).is_err());
}
