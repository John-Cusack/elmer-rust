//! Lowest-order p basis against the aggregates printed by the pristine
//! ElmerSolver for `ElementUnitTest_p1` (and `ElementUnitTest_h1`, whose
//! Lagrange elements span the same space):
//! `Norm = SUM(MASS)` (mass pass) and `SUM(ABS(STIFF))` (stiffness pass) of a
//! single unit element. Only these aggregates are printed, so this is not a
//! bitwise per-function comparison; tolerance is 1e-11 relative. Reference
//! values are the 17-digit `Norm for element` prints of the `p1` run.
//! The quadrature below is test-local (collapsed Gauss-Legendre on the p
//! reference element; Elmer's own tables are exact for these integrands, the
//! pyramid's rational gradients included).

use elmer_core::fe::{BasisEval, ElementBasis};
use elmer_core::mesh::ElementType;
use elmer_fe_hp::{PElementBasis, PElementSpec};

fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let pn = if n == 1 { x } else { p1 };
            let pn1 = if n == 1 { 1.0 } else { p0 };
            dp = n as f64 * (x * pn - pn1) / (x * x - 1.0);
            let dx = pn / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

/// (u,v,w,weight) on the p reference element of `ecode`.
fn rule(ecode: u16, n: usize) -> Vec<([f64; 3], f64)> {
    let g = gauss_legendre(n);
    // GL on [0,1]
    let g01: Vec<(f64, f64)> = g.iter().map(|&(x, w)| ((x + 1.0) / 2.0, w / 2.0)).collect();
    let s3 = 3f64.sqrt();
    let mut pts = Vec::new();
    let tri = |a: f64, b: f64| -> ([f64; 2], f64) {
        // collapsed: l2 = a, l3 = b(1-a), l1 = rest; the unit-simplex measure
        // integrates to 1/2, so scale by 2 * area (sqrt3)
        let l3 = b * (1.0 - a);
        let l2 = a;
        // vertices (-1,0),(1,0),(0,sqrt3): x = -l1 + l2, y = sqrt3 l3
        let l1 = 1.0 - l2 - l3;
        ([-l1 + l2, s3 * l3], 2.0 * (1.0 - a) * s3)
    };
    match ecode {
        303 => {
            for &(a, wa) in &g01 {
                for &(b, wb) in &g01 {
                    let (p, j) = tri(a, b);
                    pts.push(([p[0], p[1], 0.0], wa * wb * j));
                }
            }
        }
        404 => {
            for &(a, wa) in &g {
                for &(b, wb) in &g {
                    pts.push(([a, b, 0.0], wa * wb));
                }
            }
        }
        504 => {
            // vertices V1(-1,0,0) V2(1,0,0) V3(0,s3,0) V4(0,1/s3,sqrt(8/3)); volume 2 sqrt2 / 3
            let v4 = [0.0, 1.0 / s3, (8.0f64 / 3.0).sqrt()];
            let vol = 2.0 * 2f64.sqrt() / 3.0;
            for &(a, wa) in &g01 {
                for &(b, wb) in &g01 {
                    for &(c, wc) in &g01 {
                        // l4 = c, l3 = b(1-c), l2 = a(1-c)(1-b)
                        let l4 = c;
                        let l3 = b * (1.0 - c);
                        let l2 = a * (1.0 - c) * (1.0 - b);
                        let l1 = 1.0 - l2 - l3 - l4;
                        let x = -l1 + l2;
                        let y = s3 * l3 + v4[1] * l4;
                        let z = v4[2] * l4;
                        let jac = (1.0 - c) * (1.0 - c) * (1.0 - b) * 6.0 * vol;
                        pts.push(([x, y, z], wa * wb * wc * jac));
                    }
                }
            }
        }
        605 => {
            // base [-1,1]^2 at w=0, apex at w = sqrt2; volume 4 sqrt2 / 3
            for &(t, wt) in &g01 {
                for &(a, wa) in &g {
                    for &(b, wb) in &g {
                        let sc = 1.0 - t;
                        let jac = 2f64.sqrt() * sc * sc;
                        pts.push(([a * sc, b * sc, 2f64.sqrt() * t], wt * wa * wb * jac));
                    }
                }
            }
        }
        706 => {
            for &(a, wa) in &g01 {
                for &(b, wb) in &g01 {
                    let (p, j) = tri(a, b);
                    for &(c, wc) in &g {
                        pts.push(([p[0], p[1], c], wa * wb * wc * j));
                    }
                }
            }
        }
        808 => {
            for &(a, wa) in &g {
                for &(b, wb) in &g {
                    for &(c, wc) in &g {
                        pts.push(([a, b, c], wa * wb * wc));
                    }
                }
            }
        }
        _ => unreachable!(),
    }
    pts
}

fn geometry(ecode: u16) -> Vec<[f64; 3]> {
    match ecode {
        808 => vec![
            [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
            [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.],
        ],
        706 => vec![
            [0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
            [0., 0., 1.], [1., 0., 1.], [0., 1., 1.],
        ],
        605 => vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.]],
        504 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        404 => vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
        303 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
        _ => unreachable!(),
    }
}

fn norms(ecode: u16, n: usize) -> (f64, f64) {
    let basis = PElementBasis::new(PElementSpec::uniform(ecode, 1, None, None).unwrap());
    let xyz = geometry(ecode);
    let nfun = PElementSpec::uniform(ecode, 1, None, None).unwrap().n_basis().unwrap();
    let mut mass = 0.0;
    let mut stiff = vec![0.0; nfun * nfun];
    let mut out = BasisEval::default();
    for (uvw, wgt) in rule(ecode, n) {
        basis.eval(ElementType(ecode), &xyz, uvw, &mut out).unwrap();
        let w = wgt * out.det_j;
        let nb = out.basis.len();
        for p in 0..nb {
            for q in 0..nb {
                mass += w * out.basis[p] * out.basis[q];
                let d: f64 = (0..3).map(|c| out.dbasis_dx[p][c] * out.dbasis_dx[q][c]).sum();
                stiff[p * nb + q] += w * d;
            }
        }
    }
    // SUM(ABS(STIFF)) of the accumulated matrix, as SingleElements does.
    (mass, stiff.iter().map(|x| x.abs()).sum())
}

fn check(got: f64, want: f64, what: &str) {
    assert!(
        (got - want).abs() <= 1e-11 * want.abs(),
        "{what}: got {got:.17}, want {want:.17}"
    );
}

#[test]
fn lowest_order_unit_element_norms() {
    // (ecode, mass norm, stiffness norm) from ElementUnitTest_p1 run.log.
    let refs: [(u16, f64, f64); 6] = [
        (808, 1.0000000000000004, 5.3333333333333330),
        (706, 0.50000000000000033, 4.3333333333333348),
        (605, 0.33333333333333320, 3.4074074074074074),
        (504, 0.16666666666666666, 2.0000000000000004),
        (404, 1.0000000000000002, 5.3333333333333339),
        (303, 0.50000000000000011, 4.0000000000000000),
    ];
    for (code, m, s) in refs {
        let (mass, stiff) = norms(code, 10);
        check(mass, m, &format!("{code} mass"));
        check(stiff, s, &format!("{code} stiffness"));
    }
}
