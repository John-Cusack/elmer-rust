//! Nodal Lagrange shape functions on the reference element.
//!
//! `GlobalToLocal` and `InterpolateInElement` in the Fortran use the *nodal*
//! basis (`NodalBasisFunctions`, `NodalFirstDerivatives`), built at start-up
//! from the node coordinates and monomial lists of `elements.def` by inverting
//! a Vandermonde matrix. This module does the same from the generated
//! [`crate::element_defs`] table. The pyramid 605 is the closed rational form
//! of `InterpolateInElement3D`, not a polynomial.
//!
//! These tables exist because the Newton iteration needs the reference
//! gradient, which [`elmer_core::ElementBasis`] (physical gradients only)
//! cannot supply. Weights written into a projector still come from the
//! caller's `ElementBasis`.

use crate::element_defs::DEFS;
use std::sync::OnceLock;

/// Element dimension (`Element % TYPE % DIMENSION`) of an Elmer element code.
pub fn element_dim(code: u16) -> usize {
    match code / 100 {
        1 => 0,
        2 => 1,
        3 | 4 => 2,
        _ => 3,
    }
}

/// Reference-space node coordinates and nodal basis of one element type.
pub struct Shape {
    code: u16,
    dim: usize,
    nodes: Vec<[f64; 3]>,
    /// Per node: `(p, q, r, coeff)` terms of `coeff * u^p v^q w^r`.
    terms: Vec<Vec<(u8, u8, u8, f64)>>,
}

const MAX_DEG: usize = 4;

impl Shape {
    /// Shape table for an Elmer element code, if supported.
    pub fn get(code: u16) -> Option<&'static Shape> {
        static TABLE: OnceLock<Vec<Shape>> = OnceLock::new();
        let table = TABLE.get_or_init(|| DEFS.iter().map(Shape::build).collect());
        table.iter().find(|s| s.code == code)
    }

    pub fn code(&self) -> u16 {
        self.code
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Reference coordinates of local node `i`.
    pub fn node(&self, i: usize) -> [f64; 3] {
        self.nodes[i]
    }

    fn build(def: &crate::element_defs::Def) -> Shape {
        let n = def.basis.len();
        let dim = def.dim as usize;
        let nodes: Vec<[f64; 3]> = (0..n).map(|i| [def.u[i], def.v[i], def.w[i]]).collect();
        // Monomial powers of the j-th basis term (AddElementDescription).
        let pow = |k: usize| -> (u8, u8, u8) {
            let k = k - 1;
            match dim {
                1 => (k as u8, 0, 0),
                2 => ((k % MAX_DEG) as u8, (k / MAX_DEG) as u8, 0),
                _ => (
                    (k % MAX_DEG) as u8,
                    ((k / MAX_DEG) % MAX_DEG) as u8,
                    (k / (MAX_DEG * MAX_DEG)) as u8,
                ),
            }
        };
        let powers: Vec<(u8, u8, u8)> = def.basis.iter().map(|&b| pow(b as usize)).collect();
        // A[i][j] = monomial_j(node_i); coefficients are columns of A^-1.
        let mut a = vec![vec![0.0f64; n]; n];
        for (i, nd) in nodes.iter().enumerate() {
            for (j, &(p, q, r)) in powers.iter().enumerate() {
                a[i][j] = nd[0].powi(p as i32) * nd[1].powi(q as i32) * nd[2].powi(r as i32);
            }
        }
        let inv = invert(a);
        // Basis function of node i: sum_j inv[j][i] * monomial_j.
        let terms = (0..n)
            .map(|i| {
                powers
                    .iter()
                    .enumerate()
                    .map(|(j, &(p, q, r))| (p, q, r, inv[j][i]))
                    .collect()
            })
            .collect();
        Shape {
            code: def.code,
            dim,
            nodes,
            terms,
        }
    }

    /// Basis values `N_i(u,v,w)` into `out[..n_nodes]`.
    pub fn basis(&self, uvw: [f64; 3], out: &mut [f64]) {
        if self.code == 605 {
            return pyramid5(uvw, Some(out), None);
        }
        if linear_basis(self.code, uvw, Some(out), None) {
            return;
        }
        for (o, t) in out.iter_mut().zip(&self.terms) {
            *o = t
                .iter()
                .map(|&(p, q, r, c)| {
                    c * uvw[0].powi(p as i32) * uvw[1].powi(q as i32) * uvw[2].powi(r as i32)
                })
                .sum();
        }
    }

    /// Reference gradients `dN_i/d(u,v,w)` into `out[..n_nodes]`.
    pub fn dbasis(&self, uvw: [f64; 3], out: &mut [[f64; 3]]) {
        if self.code == 605 {
            return pyramid5(uvw, None, Some(out));
        }
        if linear_basis(self.code, uvw, None, Some(out)) {
            return;
        }
        let pw = |x: f64, p: u8| x.powi(p as i32);
        let dp = |x: f64, p: u8| {
            if p == 0 {
                0.0
            } else {
                p as f64 * x.powi(p as i32 - 1)
            }
        };
        for (o, t) in out.iter_mut().zip(&self.terms) {
            let mut g = [0.0; 3];
            for &(p, q, r, c) in t {
                g[0] += c * dp(uvw[0], p) * pw(uvw[1], q) * pw(uvw[2], r);
                g[1] += c * pw(uvw[0], p) * dp(uvw[1], q) * pw(uvw[2], r);
                g[2] += c * pw(uvw[0], p) * pw(uvw[1], q) * dp(uvw[2], r);
            }
            *o = g;
        }
    }
}

/// Closed forms the Fortran uses for the linear element types
/// (`*NodalLBasisAll`, `*NodalPBasisAll` in `PElementBase.F90`). They differ
/// from the table evaluation only in rounding, but `GlobalToLocal` stops on a
/// residual threshold at machine precision, so the same expressions keep the
/// iteration count, and hence the result, identical. Returns false for other
/// element types.
fn linear_basis(
    code: u16,
    uvw: [f64; 3],
    basis: Option<&mut [f64]>,
    grad: Option<&mut [[f64; 3]]>,
) -> bool {
    let [u, v, w] = uvw;
    const US: [f64; 8] = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
    const VS: [f64; 8] = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
    const WS: [f64; 8] = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
    match code {
        202 => {
            if let Some(b) = basis {
                for k in 0..2 {
                    b[k] = 0.5 * (1.0 + US[k] * u);
                }
            }
            if let Some(g) = grad {
                g[0] = [-0.5, 0.0, 0.0];
                g[1] = [0.5, 0.0, 0.0];
            }
        }
        303 => {
            if let Some(b) = basis {
                b[0] = 1.0 - u - v;
                b[1] = u;
                b[2] = v;
            }
            if let Some(g) = grad {
                g[0] = [-1.0, -1.0, 0.0];
                g[1] = [1.0, 0.0, 0.0];
                g[2] = [0.0, 1.0, 0.0];
            }
        }
        404 => {
            if let Some(b) = basis {
                for k in 0..4 {
                    b[k] = 0.25 * (1.0 + US[k] * u) * (1.0 + VS[k] * v);
                }
            }
            if let Some(g) = grad {
                for k in 0..4 {
                    g[k] = [
                        0.25 * US[k] * (1.0 + VS[k] * v),
                        0.25 * (1.0 + US[k] * u) * VS[k],
                        0.0,
                    ];
                }
            }
        }
        504 => {
            if let Some(b) = basis {
                b[0] = 1.0 - u - v - w;
                b[1] = u;
                b[2] = v;
                b[3] = w;
            }
            if let Some(g) = grad {
                g[0] = [-1.0, -1.0, -1.0];
                g[1] = [1.0, 0.0, 0.0];
                g[2] = [0.0, 1.0, 0.0];
                g[3] = [0.0, 0.0, 1.0];
            }
        }
        706 => {
            let tri = [1.0 - u - v, u, v];
            let line = [0.5 * (1.0 - w), 0.5 * (1.0 + w)];
            if let Some(b) = basis {
                for k in 0..3 {
                    b[k] = line[0] * tri[k];
                    b[3 + k] = line[1] * tri[k];
                }
            }
            if let Some(g) = grad {
                let gt = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
                for k in 0..3 {
                    g[k] = [gt[k][0] * line[0], gt[k][1] * line[0], tri[k] * -0.5];
                    g[3 + k] = [gt[k][0] * line[1], gt[k][1] * line[1], tri[k] * 0.5];
                }
            }
        }
        808 => {
            if let Some(b) = basis {
                for k in 0..8 {
                    b[k] = 0.125 * (1.0 + US[k] * u) * (1.0 + VS[k] * v) * (1.0 + WS[k] * w);
                }
            }
            if let Some(g) = grad {
                for k in 0..8 {
                    g[k] = [
                        0.125 * US[k] * (1.0 + VS[k] * v) * (1.0 + WS[k] * w),
                        0.125 * (1.0 + US[k] * u) * VS[k] * (1.0 + WS[k] * w),
                        0.125 * (1.0 + US[k] * u) * (1.0 + VS[k] * v) * WS[k],
                    ];
                }
            }
        }
        _ => return false,
    }
    true
}

/// Pyramid 605 nodal basis and gradient (`InterpolateInElement3D`,
/// `FirstDerivativeIn{U,V,W}3D`), with the Fortran's `w == 1` guard.
fn pyramid5(uvw: [f64; 3], basis: Option<&mut [f64]>, grad: Option<&mut [[f64; 3]]>) {
    let [u, v, mut w] = uvw;
    if w == 1.0 {
        w = 1.0 - 1.0e-12;
    }
    let s = 1.0 / (1.0 - w);
    // (sign of the uvw/(1-w) term, sign of u, sign of v) for nodes 1..4.
    const SG: [(f64, f64, f64); 4] = [
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, 1.0),
    ];
    if let Some(b) = basis {
        for (k, &(t, su, sv)) in SG.iter().enumerate() {
            b[k] = ((1.0 + su * u) * (1.0 + sv * v) - w + t * u * v * w * s) / 4.0;
        }
        b[4] = w;
    }
    if let Some(g) = grad {
        for (k, &(t, su, sv)) in SG.iter().enumerate() {
            g[k] = [
                (su * (1.0 + sv * v) + t * v * w * s) / 4.0,
                (sv * (1.0 + su * u) + t * u * w * s) / 4.0,
                (-1.0 + t * u * v * s * s) / 4.0,
            ];
        }
        g[4] = [0.0, 0.0, 1.0];
    }
}

/// Gauss-Jordan inverse with partial pivoting (the table is at most 27x27).
fn invert(mut a: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut inv: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for c in 0..n {
        let piv = (c..n)
            .max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))
            .unwrap_or(c);
        a.swap(c, piv);
        inv.swap(c, piv);
        let d = 1.0 / a[c][c];
        for j in 0..n {
            a[c][j] *= d;
            inv[c][j] *= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r][c];
                if f != 0.0 {
                    for j in 0..n {
                        a[r][j] -= f * a[c][j];
                        inv[r][j] -= f * inv[c][j];
                    }
                }
            }
        }
    }
    inv
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_of_unity_and_kronecker() {
        for def in DEFS {
            let sh = Shape::get(def.code).unwrap();
            let n = sh.n_nodes();
            let mut b = vec![0.0; n];
            for i in 0..n {
                sh.basis(sh.node(i), &mut b);
                for (j, bj) in b.iter().enumerate() {
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (bj - want).abs() < 1e-12,
                        "{} node {i} fn {j}: {bj}",
                        def.code
                    );
                }
            }
            // an interior-ish point
            let p = match sh.dim() {
                1 => [0.2, 0.0, 0.0],
                2 if def.code / 100 == 3 => [0.2, 0.3, 0.0],
                2 => [0.2, -0.3, 0.0],
                _ if def.code / 100 == 5 => [0.2, 0.3, 0.1],
                _ if def.code / 100 == 7 => [0.2, 0.3, 0.1],
                _ => [0.2, -0.3, 0.1],
            };
            sh.basis(p, &mut b);
            let s: f64 = b.iter().sum();
            assert!((s - 1.0).abs() < 1e-12, "{} sum {s}", def.code);
        }
    }

    #[test]
    fn gradient_matches_finite_difference() {
        for def in DEFS {
            let sh = Shape::get(def.code).unwrap();
            let n = sh.n_nodes();
            let p = [0.21, 0.17, 0.13];
            let mut g = vec![[0.0; 3]; n];
            sh.dbasis(p, &mut g);
            let h = 1e-6;
            for d in 0..sh.dim() {
                let (mut bp, mut bm) = (vec![0.0; n], vec![0.0; n]);
                let (mut pp, mut pm) = (p, p);
                pp[d] += h;
                pm[d] -= h;
                sh.basis(pp, &mut bp);
                sh.basis(pm, &mut bm);
                for k in 0..n {
                    let fd = (bp[k] - bm[k]) / (2.0 * h);
                    assert!(
                        (fd - g[k][d]).abs() < 1e-6,
                        "{} fn {k} d{d}: {fd} vs {}",
                        def.code,
                        g[k][d]
                    );
                }
            }
        }
    }
}
