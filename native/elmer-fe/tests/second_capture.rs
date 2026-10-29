//! `ElementInfo(..., ddBasisddx, SecondDerivatives=.TRUE.)` bit for bit
//! against the pristine binary. Captured with gdb after
//! `GlobalSecondDerivatives` while `AdvectionDiffusion` (`Stabilize = True`,
//! constant convection) assembled the same curved, distorted elements as
//! `pristine_capture.rs`, in 1D, 2D and 3D models; at most 8 points per
//! element are kept. Format: tests/data/second.txt, reals as f64 bits in hex.

use elmer_core::{BasisEval, ElementType};
use elmer_fe::{element_info_second, element_type, LagrangeBasis, MetricStatus};

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

#[test]
fn second_derivatives_bitwise() {
    let p = format!("{}/tests/data/second.txt", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&p).unwrap();
    let mut lines = text.lines();
    let (mut code, mut cdim, mut case) = (0u16, 0usize, String::new());
    let mut xyz: Vec<[f64; 3]> = Vec::new();
    let mut status = MetricStatus::Unchecked;
    let (mut ev, mut dd) = (BasisEval::default(), Vec::new());
    let mut types = std::collections::BTreeSet::new();
    let mut points = 0;
    while let Some(l) = lines.next() {
        let t: Vec<&str> = l.split_whitespace().collect();
        match t[0] {
            "elem" => {
                case = t[1].to_string();
                code = t[2].parse().unwrap();
                let n: usize = t[3].parse().unwrap();
                cdim = t[4].parse().unwrap();
                let x: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(hex)
                    .collect();
                xyz = (0..n).map(|i| [x[i], x[n + i], x[2 * n + i]]).collect();
                status = MetricStatus::Unchecked;
                types.insert(code);
            }
            "pt" => {
                let v: Vec<f64> = t[1..].iter().map(|s| hex(s)).collect();
                let want: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(hex)
                    .collect();
                let def = element_type(code).unwrap();
                let uvw = [v[0], v[1], v[2]];
                element_info_second(def, &xyz, uvw, cdim, &mut status, &mut ev, &mut dd).unwrap();
                let ctx = format!("{case} type {code} uvw {uvw:?}");
                assert_eq!(ev.det_j.to_bits(), v[3].to_bits(), "{ctx}: detJ");
                for q in 0..def.n_nodes {
                    for i in 0..3 {
                        for j in 0..3 {
                            // Fortran ddBasisddx(q,i,j), column-major.
                            let w = want[q * 9 + i * 3 + j];
                            assert!(
                                dd[q][i][j].to_bits() == w.to_bits(),
                                "{ctx}: ddBasisddx({q},{i},{j}) {:e} vs {w:e}",
                                dd[q][i][j]
                            );
                        }
                    }
                }
                points += 1;
            }
            k => panic!("unexpected {k}"),
        }
    }
    let want: Vec<u16> = vec![
        202, 203, 204, 303, 306, 310, 404, 408, 409, 412, 416, 504, 510, 605, 613, 706, 715, 718,
        808, 820, 827,
    ];
    assert_eq!(types.into_iter().collect::<Vec<_>>(), want);
    assert!(points > 150, "only {points} points");
}

#[test]
fn hessian_of_a_quadratic_field_is_exact_and_symmetric() {
    // sum_p f(x_p) d2N_p reproduces the Hessian of f = x^2 + 3xy - y^2 + 2yz
    // on a sheared 820 element: under an affine map f stays quadratic in the
    // reference coordinates, which the serendipity space contains.
    let def = element_type(820).unwrap();
    let xyz: Vec<[f64; 3]> = (0..def.n_nodes)
        .map(|k| {
            let (u, v, w) = (def.node_u[k], def.node_v[k], def.node_w[k]);
            [u + 0.1 * v + 0.3, 0.9 * v + 0.05 * w, w + 0.02 * u - 0.1]
        })
        .collect();
    let f = |p: &[f64; 3]| p[0] * p[0] + 3.0 * p[0] * p[1] - p[1] * p[1] + 2.0 * p[1] * p[2];
    let mut dd = Vec::new();
    LagrangeBasis::new(3)
        .eval_hessian(ElementType(820), &xyz, [0.3, -0.2, 0.1], &mut dd)
        .unwrap();
    let mut h = [[0.0; 3]; 3];
    for (p, d) in xyz.iter().zip(&dd) {
        for i in 0..3 {
            for j in 0..3 {
                h[i][j] += f(p) * d[i][j];
            }
        }
    }
    let want = [[2.0, 3.0, 0.0], [3.0, -2.0, 2.0], [0.0, 2.0, 0.0]];
    for i in 0..3 {
        for j in 0..3 {
            assert!((h[i][j] - want[i][j]).abs() < 1e-10, "{h:?}");
            // Elmer does not symmetrise: (i,j) and (j,i) agree to rounding.
            assert!((dd[0][i][j] - dd[0][j][i]).abs() < 1e-14);
        }
    }
}
