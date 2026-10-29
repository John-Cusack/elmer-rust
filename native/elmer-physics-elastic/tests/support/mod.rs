//! Test-local P1 Lagrange basis and Gauss rules, standing in for `elmer-fe`
//! until S1 lands. Test-only by design (swarm plan 4.4 rule 5).
#![allow(dead_code)]

use elmer_core::{BasisEval, ElementBasis, ElementType, QuadPoint, Quadrature};

/// Reference shape functions and local derivatives for 202, 303, 404, 504.
fn shape(etype: ElementType, uvw: [f64; 3]) -> (Vec<f64>, Vec<[f64; 3]>) {
    let [u, v, w] = uvw;
    match etype.0 {
        202 => (
            vec![(1.0 - u) / 2.0, (1.0 + u) / 2.0],
            vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]],
        ),
        303 => (
            vec![1.0 - u - v, u, v],
            vec![[-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ),
        404 => (
            vec![
                (1.0 - u) * (1.0 - v) / 4.0,
                (1.0 + u) * (1.0 - v) / 4.0,
                (1.0 + u) * (1.0 + v) / 4.0,
                (1.0 - u) * (1.0 + v) / 4.0,
            ],
            vec![
                [-(1.0 - v) / 4.0, -(1.0 - u) / 4.0, 0.0],
                [(1.0 - v) / 4.0, -(1.0 + u) / 4.0, 0.0],
                [(1.0 + v) / 4.0, (1.0 + u) / 4.0, 0.0],
                [-(1.0 + v) / 4.0, (1.0 - u) / 4.0, 0.0],
            ],
        ),
        504 => (
            vec![1.0 - u - v - w, u, v, w],
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ),
        808 => {
            const S: [[f64; 3]; 8] = [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ];
            let f = |s: f64, x: f64| 1.0 + s * x;
            let mut n = vec![];
            let mut d = vec![];
            for q in &S {
                n.push(f(q[0], u) * f(q[1], v) * f(q[2], w) / 8.0);
                d.push([
                    q[0] * f(q[1], v) * f(q[2], w) / 8.0,
                    f(q[0], u) * q[1] * f(q[2], w) / 8.0,
                    f(q[0], u) * f(q[1], v) * q[2] / 8.0,
                ]);
            }
            (n, d)
        }
        c => panic!("support basis: element {c}"),
    }
}

fn edim(etype: ElementType) -> usize {
    match etype.family() {
        2 => 1,
        3 | 4 => 2,
        5 | 8 => 3,
        f => panic!("family {f}"),
    }
}

#[allow(clippy::needless_range_loop)]
fn inv(m: &[[f64; 3]; 3], n: usize) -> [[f64; 3]; 3] {
    let mut o = [[0.0; 3]; 3];
    match n {
        1 => o[0][0] = 1.0 / m[0][0],
        2 => {
            let d = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            o[0][0] = m[1][1] / d;
            o[1][1] = m[0][0] / d;
            o[0][1] = -m[0][1] / d;
            o[1][0] = -m[1][0] / d;
        }
        _ => {
            let d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            for i in 0..3 {
                for j in 0..3 {
                    let (r1, r2) = ((j + 1) % 3, (j + 2) % 3);
                    let (c1, c2) = ((i + 1) % 3, (i + 2) % 3);
                    o[i][j] = (m[r1][c1] * m[r2][c2] - m[r1][c2] * m[r2][c1]) / d;
                }
            }
        }
    }
    o
}

/// ElementInfo-like evaluation: physical gradients, tangential for elements of
/// lower dimension than their embedding (`dN = dL (J^T J)^-1 J^T`).
pub struct P1Basis;

impl ElementBasis for P1Basis {
    fn eval(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let (n, dl) = shape(etype, uvw);
        let ed = edim(etype);
        // J[a][e] = d x_a / d xi_e
        let mut j = [[0.0; 3]; 3];
        for (k, d) in dl.iter().enumerate() {
            for a in 0..3 {
                for e in 0..ed {
                    j[a][e] += xyz[k][a] * d[e];
                }
            }
        }
        let mut g = [[0.0; 3]; 3];
        for e in 0..ed {
            for f in 0..ed {
                g[e][f] = (0..3).map(|a| j[a][e] * j[a][f]).sum();
            }
        }
        let gi = inv(&g, ed);
        let det_g = match ed {
            1 => g[0][0],
            2 => g[0][0] * g[1][1] - g[0][1] * g[1][0],
            _ => {
                g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                    - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
                    + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0])
            }
        };
        out.det_j = det_g.sqrt();
        out.basis = n;
        out.dbasis_dx = dl
            .iter()
            .map(|d| {
                let mut o = [0.0; 3];
                for a in 0..3 {
                    for e in 0..ed {
                        for f in 0..ed {
                            o[a] += d[e] * gi[e][f] * j[a][f];
                        }
                    }
                }
                o
            })
            .collect();
        Ok(())
    }
}

/// Gauss rules. `np = None` is the Elmer default for the linear element.
pub struct Gauss;

fn pt(u: f64, v: f64, w: f64, weight: f64) -> QuadPoint {
    QuadPoint { uvw: [u, v, w], weight }
}

impl Quadrature for Gauss {
    fn gauss_points(&self, etype: ElementType, np: Option<usize>) -> Vec<QuadPoint> {
        let g = 1.0 / 3f64.sqrt();
        match (etype.0, np) {
            (202, _) => vec![pt(-g, 0.0, 0.0, 1.0), pt(g, 0.0, 0.0, 1.0)],
            (303, None | Some(3)) => vec![
                pt(1.0 / 6.0, 1.0 / 6.0, 0.0, 1.0 / 6.0),
                pt(2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 6.0),
                pt(1.0 / 6.0, 2.0 / 3.0, 0.0, 1.0 / 6.0),
            ],
            (303, Some(1)) => vec![pt(1.0 / 3.0, 1.0 / 3.0, 0.0, 0.5)],
            (404, _) => vec![
                pt(-g, -g, 0.0, 1.0),
                pt(g, -g, 0.0, 1.0),
                pt(g, g, 0.0, 1.0),
                pt(-g, g, 0.0, 1.0),
            ],
            (504, _) => vec![pt(0.25, 0.25, 0.25, 1.0 / 6.0)],
            (808, _) => {
                let mut v = vec![];
                for k in [-g, g] {
                    for j in [-g, g] {
                        for i in [-g, g] {
                            v.push(pt(i, j, k, 1.0));
                        }
                    }
                }
                v
            }
            (c, n) => panic!("support quadrature: {c} np {n:?}"),
        }
    }
}

pub fn elem_center(xyz: &[[f64; 3]]) -> [f64; 3] {
    let n = xyz.len() as f64;
    std::array::from_fn(|a| xyz.iter().map(|p| p[a]).sum::<f64>() / n)
}

/// Scatter an element block into a dense global system (node-major, `cdim`
/// components), for patch tests.
pub fn scatter(
    k: &mut [f64],
    f: &mut [f64],
    ndof_global: usize,
    cdim: usize,
    nodes: &[usize],
    local: &elmer_core::LocalSystem,
) {
    for (p, &gp) in nodes.iter().enumerate() {
        for i in 0..cdim {
            f[gp * cdim + i] += local.force[p * cdim + i];
            for (q, &gq) in nodes.iter().enumerate() {
                for j in 0..cdim {
                    k[(gp * cdim + i) * ndof_global + gq * cdim + j] +=
                        local.stiff[(p * cdim + i) * local.n + q * cdim + j];
                }
            }
        }
    }
}

pub mod cases;
pub mod sim;
