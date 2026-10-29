//! Lagrange bases in Elmer node order: 202, 203, 204, 303, 306, 404, 408, 504.
//!
//! `det_j` is the surface/length measure `sqrt(det(J J^T))`, so boundary
//! elements embedded in a larger space integrate correctly. Gradients are
//! taken from the leading `d x d` block of `J` (`d` = reference dimension),
//! which is right for bulk elements in the xy / xyz coordinate axes and is not
//! meaningful for embedded boundary elements (the kernels never read it there).

use elmer_core::{BasisEval, ElementBasis, ElementType};

pub struct TestBasis;

/// Reference dimension, shape values and reference gradients at `uvw`.
fn shape(etype: u16, uvw: [f64; 3]) -> Option<(usize, Vec<f64>, Vec<[f64; 3]>)> {
    let [u, v, w] = uvw;
    Some(match etype {
        202 => (
            1,
            vec![(1.0 - u) / 2.0, (1.0 + u) / 2.0],
            vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]],
        ),
        203 => (
            1,
            vec![u * (u - 1.0) / 2.0, u * (u + 1.0) / 2.0, 1.0 - u * u],
            vec![
                [u - 0.5, 0.0, 0.0],
                [u + 0.5, 0.0, 0.0],
                [-2.0 * u, 0.0, 0.0],
            ],
        ),
        204 => {
            // Nodes at u = -1, 1, -1/3, 1/3 (elements.def).
            let (n, g) = lagrange_line(&[-1.0, 1.0, -1.0 / 3.0, 1.0 / 3.0], u);
            (1, n, g.into_iter().map(|d| [d, 0.0, 0.0]).collect())
        }
        303 => (
            2,
            vec![1.0 - u - v, u, v],
            vec![[-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ),
        306 => {
            let l = [1.0 - u - v, u, v];
            let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            let mut n = vec![0.0; 6];
            let mut g = vec![[0.0; 3]; 6];
            for i in 0..3 {
                n[i] = l[i] * (2.0 * l[i] - 1.0);
                for k in 0..2 {
                    g[i][k] = (4.0 * l[i] - 1.0) * dl[i][k];
                }
            }
            for (m, (a, b)) in [(0usize, 1usize), (1, 2), (2, 0)].iter().enumerate() {
                n[3 + m] = 4.0 * l[*a] * l[*b];
                for k in 0..2 {
                    g[3 + m][k] = 4.0 * (dl[*a][k] * l[*b] + l[*a] * dl[*b][k]);
                }
            }
            (2, n, g)
        }
        404 => {
            let su = [-1.0, 1.0, 1.0, -1.0];
            let sv = [-1.0, -1.0, 1.0, 1.0];
            let mut n = vec![0.0; 4];
            let mut g = vec![[0.0; 3]; 4];
            for i in 0..4 {
                n[i] = (1.0 + su[i] * u) * (1.0 + sv[i] * v) / 4.0;
                g[i][0] = su[i] * (1.0 + sv[i] * v) / 4.0;
                g[i][1] = sv[i] * (1.0 + su[i] * u) / 4.0;
            }
            (2, n, g)
        }
        408 => {
            // Serendipity: corners 1-4, mid-edge 5 (1-2), 6 (2-3), 7 (3-4), 8 (4-1).
            let su = [-1.0, 1.0, 1.0, -1.0];
            let sv = [-1.0, -1.0, 1.0, 1.0];
            let mut n = vec![0.0; 8];
            let mut g = vec![[0.0; 3]; 8];
            for i in 0..4 {
                let (a, b) = (1.0 + su[i] * u, 1.0 + sv[i] * v);
                let c = su[i] * u + sv[i] * v - 1.0;
                n[i] = a * b * c / 4.0;
                g[i][0] = su[i] * b * c / 4.0 + a * b * su[i] / 4.0;
                g[i][1] = a * sv[i] * c / 4.0 + a * b * sv[i] / 4.0;
            }
            n[4] = (1.0 - u * u) * (1.0 - v) / 2.0;
            g[4] = [-u * (1.0 - v), -(1.0 - u * u) / 2.0, 0.0];
            n[5] = (1.0 + u) * (1.0 - v * v) / 2.0;
            g[5] = [(1.0 - v * v) / 2.0, -(1.0 + u) * v, 0.0];
            n[6] = (1.0 - u * u) * (1.0 + v) / 2.0;
            g[6] = [-u * (1.0 + v), (1.0 - u * u) / 2.0, 0.0];
            n[7] = (1.0 - u) * (1.0 - v * v) / 2.0;
            g[7] = [-(1.0 - v * v) / 2.0, -(1.0 - u) * v, 0.0];
            (2, n, g)
        }
        504 => (
            3,
            vec![1.0 - u - v - w, u, v, w],
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ),
        _ => return None,
    })
}

/// One-dimensional Lagrange polynomials through `nodes` and their
/// derivatives at `u`.
fn lagrange_line(nodes: &[f64], u: f64) -> (Vec<f64>, Vec<f64>) {
    let mut n = Vec::with_capacity(nodes.len());
    let mut g = Vec::with_capacity(nodes.len());
    for (i, &xi) in nodes.iter().enumerate() {
        let others: Vec<f64> = nodes.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &x)| x).collect();
        let den: f64 = others.iter().map(|&xj| xi - xj).product();
        n.push(others.iter().map(|&xj| u - xj).product::<f64>() / den);
        let mut d = 0.0;
        for k in 0..others.len() {
            d += others.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &xj)| u - xj).product::<f64>();
        }
        g.push(d / den);
    }
    (n, g)
}

impl ElementBasis for TestBasis {
    fn eval(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let (d, n, g) = shape(etype.0, uvw).ok_or_else(|| format!("element type {}", etype.0))?;
        if xyz.len() != n.len() {
            return Err(format!("{} nodes for type {}", xyz.len(), etype.0));
        }
        // J[a][c] = sum_i dN_i/du_a x_i[c]
        let mut j = [[0.0; 3]; 3];
        for (i, x) in xyz.iter().enumerate() {
            for a in 0..d {
                for c in 0..3 {
                    j[a][c] += g[i][a] * x[c];
                }
            }
        }
        // Surface measure sqrt(det(J J^T)).
        let mut jjt = [[0.0; 3]; 3];
        for a in 0..d {
            for b in 0..d {
                jjt[a][b] = (0..3).map(|c| j[a][c] * j[b][c]).sum();
            }
        }
        let det_j = match d {
            1 => jjt[0][0].sqrt(),
            2 => (jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0]).sqrt(),
            _ => {
                let m = &jjt;
                (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
                    .sqrt()
            }
        };
        // Inverse of the leading d x d block of J.
        let inv = invert(&j, d)?;
        out.det_j = det_j;
        out.basis = n;
        out.dbasis_dx = g
            .iter()
            .map(|gi| {
                let mut dx = [0.0; 3];
                for c in 0..d {
                    for a in 0..d {
                        dx[c] += inv[c][a] * gi[a];
                    }
                }
                dx
            })
            .collect();
        Ok(())
    }
}

fn invert(j: &[[f64; 3]; 3], d: usize) -> Result<[[f64; 3]; 3], String> {
    let mut r = [[0.0; 3]; 3];
    match d {
        1 => r[0][0] = 1.0 / j[0][0],
        2 => {
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            if det == 0.0 {
                return Err("singular jacobian".into());
            }
            r[0][0] = j[1][1] / det;
            r[0][1] = -j[0][1] / det;
            r[1][0] = -j[1][0] / det;
            r[1][1] = j[0][0] / det;
        }
        _ => {
            let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
                - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
                + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
            if det == 0.0 {
                return Err("singular jacobian".into());
            }
            for a in 0..3 {
                for b in 0..3 {
                    let (a1, a2) = ((a + 1) % 3, (a + 2) % 3);
                    let (b1, b2) = ((b + 1) % 3, (b + 2) % 3);
                    // Cofactor transpose.
                    r[b][a] = (j[a1][b1] * j[a2][b2] - j[a1][b2] * j[a2][b1]) / det;
                }
            }
        }
    }
    Ok(r)
}
