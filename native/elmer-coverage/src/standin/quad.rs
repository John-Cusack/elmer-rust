//! Gauss rules on the Elmer reference elements. Weights sum to the reference
//! element measure (line 2, triangle 1/2, quad 4, tet 1/6).

use elmer_core::{ElementType, QuadPoint, Quadrature};

pub struct TestQuadrature;

fn gauss_line(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => {
            let a = 1.0 / 3f64.sqrt();
            vec![(-a, 1.0), (a, 1.0)]
        }
        3 => {
            let a = (0.6f64).sqrt();
            vec![(-a, 5.0 / 9.0), (0.0, 8.0 / 9.0), (a, 5.0 / 9.0)]
        }
        4 => {
            let a = (3.0 / 7.0 - 2.0 / 7.0 * (1.2f64).sqrt()).sqrt();
            let b = (3.0 / 7.0 + 2.0 / 7.0 * (1.2f64).sqrt()).sqrt();
            let wa = (18.0 + 30f64.sqrt()) / 36.0;
            let wb = (18.0 - 30f64.sqrt()) / 36.0;
            vec![(-b, wb), (-a, wa), (a, wa), (b, wb)]
        }
        _ => panic!("line rule {n}"),
    }
}

impl Quadrature for TestQuadrature {
    fn gauss_points(&self, etype: ElementType, np: Option<usize>) -> Vec<QuadPoint> {
        let pt = |uvw: [f64; 3], weight: f64| QuadPoint { uvw, weight };
        match etype.0 {
            202 | 203 => {
                let n = np.unwrap_or(if etype.0 == 202 { 2 } else { 3 });
                gauss_line(n)
                    .into_iter()
                    .map(|(u, w)| pt([u, 0.0, 0.0], w))
                    .collect()
            }
            303 | 306 => match np.unwrap_or(if etype.0 == 303 { 3 } else { 6 }) {
                1 => vec![pt([1.0 / 3.0, 1.0 / 3.0, 0.0], 0.5)],
                3 => vec![
                    pt([1.0 / 6.0, 1.0 / 6.0, 0.0], 1.0 / 6.0),
                    pt([2.0 / 3.0, 1.0 / 6.0, 0.0], 1.0 / 6.0),
                    pt([1.0 / 6.0, 2.0 / 3.0, 0.0], 1.0 / 6.0),
                ],
                6 => {
                    let a = 0.445948490915965;
                    let b = 0.091576213509771;
                    let wa = 0.223381589678011 / 2.0;
                    let wb = 0.109951743655322 / 2.0;
                    vec![
                        pt([a, a, 0.0], wa),
                        pt([1.0 - 2.0 * a, a, 0.0], wa),
                        pt([a, 1.0 - 2.0 * a, 0.0], wa),
                        pt([b, b, 0.0], wb),
                        pt([1.0 - 2.0 * b, b, 0.0], wb),
                        pt([b, 1.0 - 2.0 * b, 0.0], wb),
                    ]
                }
                n => panic!("triangle rule {n}"),
            },
            404 | 408 => {
                let np = np.unwrap_or(if etype.0 == 404 { 4 } else { 9 });
                let n = (np as f64).sqrt().round() as usize;
                let g = gauss_line(n);
                let mut v = Vec::new();
                for (vv, wv) in &g {
                    for (uu, wu) in &g {
                        v.push(pt([*uu, *vv, 0.0], wu * wv));
                    }
                }
                v
            }
            504 => match np.unwrap_or(4) {
                1 => vec![pt([0.25, 0.25, 0.25], 1.0 / 6.0)],
                4 => {
                    let a = 0.1381966011250105;
                    let b = 0.5854101966249685;
                    vec![
                        pt([a, a, a], 1.0 / 24.0),
                        pt([b, a, a], 1.0 / 24.0),
                        pt([a, b, a], 1.0 / 24.0),
                        pt([a, a, b], 1.0 / 24.0),
                    ]
                }
                n => panic!("tet rule {n}"),
            },
            t => panic!("quadrature for element type {t}"),
        }
    }
}
