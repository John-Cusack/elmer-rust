//! Every supported element type, one affinely mapped element at a time:
//! point location recovers the reference coordinates and interpolation is
//! exact for the polynomial degree the element carries.

mod support;

use elmer_core::NullLogger;
use elmer_interp::shape::Shape;
use elmer_interp::{interpolate, Extras, FieldData, InterpOptions, Locator};
use support::*;

const CODES: [u16; 20] = [
    202, 203, 204, 303, 306, 310, 404, 408, 409, 412, 416, 504, 510, 605, 706, 715, 718, 808, 820,
    827,
];

/// Polynomial degree reproduced exactly by each code.
fn degree(code: u16) -> i32 {
    match code {
        202 | 303 | 404 | 504 | 605 | 706 | 808 => 1,
        204 | 310 => 3,
        _ => 2,
    }
}

/// Deterministic points strictly inside the reference element.
fn interior_points(code: u16) -> Vec<[f64; 3]> {
    let fam = code / 100;
    let mut pts = Vec::new();
    let g = [0.15, 0.3, 0.5, 0.7];
    for &a in &g {
        for &b in &g {
            for &c in &g {
                let p = match fam {
                    2 => [2.0 * a - 1.0, 0.0, 0.0],
                    3 => [a * 0.5, b * 0.5, 0.0],
                    4 => [2.0 * a - 1.0, 2.0 * b - 1.0, 0.0],
                    5 => [a * 0.3, b * 0.3, c * 0.3],
                    6 => {
                        let w = c * 0.8;
                        [
                            (2.0 * a - 1.0) * (1.0 - w) * 0.9,
                            (2.0 * b - 1.0) * (1.0 - w) * 0.9,
                            w,
                        ]
                    }
                    7 => [a * 0.5, b * 0.4, 2.0 * c - 1.0],
                    _ => [2.0 * a - 1.0, 2.0 * b - 1.0, 2.0 * c - 1.0],
                };
                pts.push(p);
            }
        }
    }
    pts.dedup_by(|x, y| x == y);
    pts
}

fn affine(u: [f64; 3]) -> [f64; 3] {
    // Non-degenerate, non-orthogonal map with a shift.
    [
        1.5 * u[0] + 0.2 * u[1] - 0.1 * u[2] + 3.0,
        -0.3 * u[0] + 1.1 * u[1] + 0.25 * u[2] - 2.0,
        0.15 * u[0] - 0.2 * u[1] + 0.9 * u[2] + 0.5,
    ]
}

fn poly(deg: i32, x: [f64; 3]) -> f64 {
    let mut f = 1.0 + 0.7 * x[0] - 1.3 * x[1] + 0.4 * x[2];
    if deg >= 2 {
        f += 0.5 * x[0] * x[0] - 0.2 * x[1] * x[1] + 0.3 * x[0] * x[1] + 0.1 * x[2] * x[1]
            - 0.15 * x[2] * x[0];
    }
    if deg >= 3 {
        f += 0.02 * x[0] * x[0] * x[0] - 0.03 * x[0] * x[0] * x[1] + 0.01 * x[1] * x[1] * x[1];
    }
    f
}

#[test]
fn every_element_type_locates_and_reproduces_its_polynomials() {
    for code in CODES {
        let sh = Shape::get(code).unwrap_or_else(|| panic!("no shape {code}"));
        let dim = sh.dim();
        // 1-D and 2-D elements live in a coordinate subspace so the Jacobian is square.
        let map = |u: [f64; 3]| {
            let a = affine(u);
            let mut x = [0.0; 3];
            x[..dim].copy_from_slice(&a[..dim]);
            x
        };
        let coords: Vec<[f64; 3]> = (0..sh.n_nodes()).map(|i| map(sh.node(i))).collect();
        let old = TestMesh::new(
            dim,
            coords,
            vec![(code, 1, (0..sh.n_nodes() as u32).collect())],
        );
        let pts = interior_points(code);
        let new_coords: Vec<[f64; 3]> = pts.iter().map(|&u| map(u)).collect();
        let new = TestMesh::new(dim, new_coords.clone(), vec![]);

        let deg = degree(code);
        let of = vec![FieldData::nodal(
            "f",
            old.coords.iter().map(|&c| poly(deg, c)).collect(),
        )];
        let mut nf = vec![FieldData::nodal("f", vec![f64::NAN; new.coords.len()])];
        let loc = Locator::new(&old, InterpOptions::default()).unwrap();
        let rep = interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap();
        assert_eq!(rep.n_not_found, 0, "code {code}");
        for (i, &c) in new_coords.iter().enumerate() {
            let hit = loc.locate(c, None).unwrap();
            for d in 0..dim {
                assert_close(hit.uvw[d], pts[i][d], 1e-9, &format!("{code} uvw[{d}]"));
            }
            assert_close(
                nf[0].values[i],
                poly(deg, c),
                1e-9,
                &format!("{code} value"),
            );
        }
    }
}

#[test]
fn tilted_surface_element_in_three_space() {
    // A planar triangle in 3-D: GlobalToLocal uses the normal equations.
    let t = |u: f64, v: f64| [u + 0.3 * v, 0.5 * u + v, 0.4 * u - 0.6 * v];
    let coords = vec![t(0.0, 0.0), t(2.0, 0.0), t(0.0, 2.0)];
    let old = TestMesh::new(3, coords, vec![(303, 1, vec![0, 1, 2])]);
    let new = TestMesh::new(3, vec![t(0.5, 0.5), t(0.2, 1.1)], vec![]);
    let f = |p: [f64; 3]| 2.0 * p[0] - p[1] + 0.5 * p[2];
    let of = vec![FieldData::nodal(
        "f",
        old.coords.iter().map(|&c| f(c)).collect(),
    )];
    let mut nf = vec![FieldData::nodal("f", vec![0.0; 2])];
    let opts = InterpOptions {
        coordinate_system_dim: Some(3),
        ..Default::default()
    };
    let loc = Locator::new(&old, opts).unwrap();
    let rep = interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap();
    assert_eq!(rep.n_not_found, 0);
    for i in 0..2 {
        assert_close(nf[0].values[i], f(new.coords[i]), 1e-12, "tilted");
    }
}

#[test]
fn points_outside_are_rejected_for_every_family() {
    for code in CODES {
        let sh = Shape::get(code).unwrap();
        let dim = sh.dim();
        let coords: Vec<[f64; 3]> = (0..sh.n_nodes()).map(|i| sh.node(i)).collect();
        let old = TestMesh::new(
            dim,
            coords,
            vec![(code, 1, (0..sh.n_nodes() as u32).collect())],
        );
        let loc = Locator::new(&old, InterpOptions::default()).unwrap();
        let far = [3.7, 2.9, 4.1];
        let mut p = [0.0; 3];
        p[..dim].copy_from_slice(&far[..dim]);
        assert!(loc.locate(p, None).is_none(), "{code} accepted a far point");
        // Just outside the closed reference element along the first axis.
        let q = if code / 100 == 6 {
            [1.2, 0.0, 0.1]
        } else {
            [1.2, 0.0, 0.0]
        };
        assert!(
            loc.locate(q, None).is_none(),
            "{code} accepted an outside point"
        );
    }
}
