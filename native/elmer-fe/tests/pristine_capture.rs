//! Bitwise comparison against values captured from the pristine Elmer
//! binary (Debug build, `HAVE_QP`) with gdb: Gauss rules returned by the
//! Integration.F90 rule builders, `ElementInfo` outputs on curved and
//! distorted elements of every Lagrange type in 1D, 2D and 3D models, the
//! basis coefficients of each type, and `StabParam` / `ElementDiameter`.
//! Fixture format: tests/data/*.txt, reals as f64 bit patterns in hex.

use elmer_core::{BasisEval, ElementType};
use elmer_fe::quadrature::{
    gauss_points_1d, gauss_points_brick, gauss_points_pyramid, gauss_points_quad,
    gauss_points_tetra, gauss_points_triangle, gauss_points_wedge, gauss_points_wedge2,
    gauss_points_wedge_economic,
};
use elmer_fe::{
    element_info, element_type, stab_param, stab_param_mk_only, type_stabilization_mk, GaussRule,
    MetricStatus,
};

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap_or_else(|_| panic!("bad hex {s}")))
}

fn data(name: &str) -> String {
    let p = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{p}: {e}"))
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[test]
fn gauss_rules_bitwise() {
    let text = data("rules.txt");
    let mut lines = text.lines();
    let mut checked = 0;
    while let Some(head) = lines.next() {
        let t: Vec<&str> = head.split_whitespace().collect();
        assert_eq!(t[0], "rule");
        let kind = t[1];
        let args: Vec<usize> = t[2..t.len() - 1]
            .iter()
            .map(|x| x.parse().unwrap())
            .collect();
        let n: usize = t[t.len() - 1].parse().unwrap();
        let want: Vec<[f64; 4]> = (0..n)
            .map(|_| {
                let r: Vec<f64> = lines.next().unwrap().split_whitespace().map(hex).collect();
                [r[0], r[1], r[2], r[3]]
            })
            .collect();
        let got: GaussRule = match kind {
            "line" => gauss_points_1d(args[0]),
            "tri" => gauss_points_triangle(args[0]),
            "quad" => gauss_points_quad(args[0]),
            "tet" => gauss_points_tetra(args[0]),
            "pyr" => gauss_points_pyramid(args[0]),
            "wedge" => gauss_points_wedge(args[0]),
            "brick" => gauss_points_brick(args[0]),
            "wedge2" => gauss_points_wedge2(args[0], args[1]),
            "wecon" => gauss_points_wedge_economic(args[0]),
            k => panic!("unknown rule kind {k}"),
        }
        .unwrap_or_else(|e| panic!("{head}: {e}"));
        assert_eq!(got.len(), n, "{head}: point count");
        // GaussPointsQuad never writes w: Elmer's value is whatever the
        // previous rule left in the shared workspace.
        let check_w = kind != "quad";
        for (i, w) in want.iter().enumerate() {
            let g = [got.u[i], got.v[i], got.w[i], got.s[i]];
            for c in 0..4 {
                if c == 2 && !check_w {
                    continue;
                }
                assert!(
                    same(g[c], w[c]),
                    "{head}: point {i} component {c}: got {:e} want {:e}",
                    g[c],
                    w[c]
                );
            }
        }
        checked += 1;
    }
    assert!(checked > 80, "only {checked} rules");
}

struct Elem {
    case: String,
    code: u16,
    cdim: usize,
    xyz: Vec<[f64; 3]>,
    pts: Vec<Pt>,
}

/// One captured `ElementInfo` call.
struct Pt {
    /// Called from `StabParam` (on the reference element) rather than at an
    /// integration point of SingleElements.
    from_stab: bool,
    uvw: [f64; 3],
    detj: f64,
    basis: Vec<f64>,
    dbasis: Vec<[f64; 3]>,
}

fn read_info() -> Vec<Elem> {
    let text = data("info.txt");
    let mut out: Vec<Elem> = Vec::new();
    let mut lines = text.lines();
    while let Some(l) = lines.next() {
        let t: Vec<&str> = l.split_whitespace().collect();
        match t[0] {
            "elem" => {
                let n: usize = t[3].parse().unwrap();
                let x: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(hex)
                    .collect();
                out.push(Elem {
                    case: t[1].to_string(),
                    code: t[2].parse().unwrap(),
                    cdim: t[4].parse().unwrap(),
                    xyz: (0..n).map(|i| [x[i], x[n + i], x[2 * n + i]]).collect(),
                    pts: Vec::new(),
                });
            }
            "pt" => {
                let v: Vec<f64> = t[2..].iter().map(|s| hex(s)).collect();
                let b: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(hex)
                    .collect();
                let d: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(hex)
                    .collect();
                let d = d.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
                out.last_mut().unwrap().pts.push(Pt {
                    from_stab: t[1] == "stab",
                    uvw: [v[0], v[1], v[2]],
                    detj: v[3],
                    basis: b,
                    dbasis: d,
                });
            }
            k => panic!("unexpected line {k}"),
        }
    }
    out
}

#[test]
fn element_info_bitwise() {
    let elems = read_info();
    let mut types = std::collections::BTreeSet::new();
    let mut points = 0;
    for e in &elems {
        let def = element_type(e.code).unwrap();
        let mut status = MetricStatus::Unchecked;
        let mut ev = BasisEval::default();
        for (
            k,
            Pt {
                uvw,
                detj,
                basis,
                dbasis: db,
                ..
            },
        ) in e.pts.iter().enumerate()
        {
            element_info(def, &e.xyz, *uvw, e.cdim, &mut status, &mut ev)
                .unwrap_or_else(|err| panic!("{} {} pt {k}: {err}", e.case, e.code));
            let ctx = format!(
                "{} type {} cdim {} pt {k} uvw {uvw:?}",
                e.case, e.code, e.cdim
            );
            assert!(
                same(ev.det_j, *detj),
                "{ctx}: detJ {:e} vs {:e}",
                ev.det_j,
                detj
            );
            for i in 0..def.n_nodes {
                assert!(
                    same(ev.basis[i], basis[i]),
                    "{ctx}: basis {i} {:e} vs {:e}",
                    ev.basis[i],
                    basis[i]
                );
                for j in 0..3 {
                    assert!(
                        same(ev.dbasis_dx[i][j], db[i][j]),
                        "{ctx}: dBasisdx({i},{j}) {:e} vs {:e}",
                        ev.dbasis_dx[i][j],
                        db[i][j]
                    );
                }
            }
            points += 1;
        }
        assert_eq!(status, MetricStatus::Double, "{} {}", e.case, e.code);
        types.insert(e.code);
    }
    let want: Vec<u16> = vec![
        202, 203, 204, 303, 306, 310, 404, 408, 409, 412, 416, 504, 510, 605, 613, 706, 715, 718,
        808, 820, 827,
    ];
    assert_eq!(types.into_iter().collect::<Vec<_>>(), want);
    assert!(points > 700, "only {points} points");
}

#[test]
fn basis_coefficients_bitwise() {
    let text = data("coeff.txt");
    let mut lines = text.lines();
    let mut n_types = 0;
    while let Some(head) = lines.next() {
        let t: Vec<&str> = head.split_whitespace().collect();
        let code: u16 = t[1].parse().unwrap();
        let n: usize = t[2].parse().unwrap();
        let def = element_type(code).unwrap();
        for i in 0..n {
            let r: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
            let k: usize = r[2].parse().unwrap();
            assert_eq!(def.basis[i].len(), k, "type {code} basis {i}");
            for (j, term) in r[3..].iter().enumerate() {
                let f: Vec<&str> = term.split(',').collect();
                let bt = def.basis[i][j];
                assert_eq!(
                    (bt.p as u32, bt.q as u32, bt.r as u32),
                    (
                        f[0].parse().unwrap(),
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap()
                    ),
                    "type {code} basis {i} term {j} powers"
                );
                assert!(
                    same(bt.coeff, hex(f[3])),
                    "type {code} basis {i} term {j}: {:e} vs {:e}",
                    bt.coeff,
                    hex(f[3])
                );
            }
        }
        n_types += 1;
    }
    assert!(n_types >= 21, "only {n_types} types");
}

/// Units in the last place between two finite values of one sign.
fn ulps(a: f64, b: f64) -> u64 {
    (a.to_bits() as i64 - b.to_bits() as i64).unsigned_abs()
}

#[test]
fn stabilization_parameters() {
    let elems = read_info();
    let geom = |case: &str, code: u16| {
        elems
            .iter()
            .find(|e| e.case == case && e.code == code && e.pts.iter().any(|p| !p.from_stab))
            .unwrap_or_else(|| panic!("no element {case} {code}"))
    };
    let text = data("stab.txt");
    let mut n = 0;
    for l in text.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        let (case, kind, code) = (t[0], t[1], t[2].parse::<u16>().unwrap());
        let def = element_type(code).unwrap();
        // Linear types are closed-form and must match bit for bit; the
        // eigenvalue path of higher-order types goes through LAPACK DSYGV in
        // Elmer, which is not transliterated (yet), so those get 64 ulps.
        let tol = if def.basis_degree <= 1 { 0 } else { 64 };
        let e = geom(case, code);
        let check = |what: &str, got: f64, want: f64| {
            let d = ulps(got, want);
            assert!(d <= tol, "{l}: {what} got {got:e} want {want:e} ({d} ulps)");
        };
        match kind {
            "typemk" => check(
                "type mK",
                type_stabilization_mk(ElementType(code)).unwrap(),
                hex(t[3]),
            ),
            "stab" => {
                let long = t[3] == "1";
                let (mk, hk) = stab_param(ElementType(code), &e.xyz, e.cdim, long).unwrap();
                check("mK", mk, hex(t[4]));
                check("hK", hk, hex(t[5]));
            }
            "stabmk" => check(
                "mK (no hK)",
                stab_param_mk_only(ElementType(code), &e.xyz, e.cdim).unwrap(),
                hex(t[3]),
            ),
            k => panic!("unexpected {k}"),
        }
        n += 1;
    }
    assert!(n > 100, "only {n} stab lines");
}
