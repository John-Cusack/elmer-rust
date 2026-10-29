//! Reference-element nodal basis functions and their local gradients:
//! `NodalBasisFunctions` and `NodalFirstDerivatives` (ElementBasis.F90:2148,
//! 2253) for non-p elements.
//!
//! The linear line, triangle, quad, tetra, wedge and brick take Elmer's
//! closed-form fast path (PElementBase.F90 `*NodalPBasisAll` /
//! `*NodalLBasisAll`); the other types sum the elements.def monomials, and the
//! 605 and 613 pyramids use their rational special forms. Expression order
//! follows the Fortran so the results agree bit for bit.

use crate::element_type::{powi, BasisTerm, ElementTypeDef};

/// Sign tables of PElementBase.F90 (`usgn`, `vsgn`, `wsgn`).
const QUAD_U: [f64; 4] = [-1.0, 1.0, 1.0, -1.0];
const QUAD_V: [f64; 4] = [-1.0, -1.0, 1.0, 1.0];
const BRICK_U: [f64; 8] = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
const BRICK_V: [f64; 8] = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
const BRICK_W: [f64; 8] = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];

/// Whether `NodalBasisFunctions` takes its closed-form fast path.
fn fast_path(def: &ElementTypeDef) -> bool {
    let family = def.code / 100;
    family != 6 && family >= def.code % 100
}

/// Elmer's pyramid apex guard: `IF ( w == 1 ) w = 1.0d0-1.0d-12`.
fn apex_guard(w: f64) -> f64 {
    if w == 1.0 {
        1.0 - 1.0e-12
    } else {
        w
    }
}

/// Basis values at `(u,v,w)` into `basis[0..n]`.
pub fn nodal_basis(def: &ElementTypeDef, u: f64, v: f64, w: f64, basis: &mut [f64]) {
    let n = def.n_nodes;
    if fast_path(def) {
        match def.family() {
            2 => {
                let c = 1.0 / 2.0;
                // c*(1+usgn*u) with usgn = -1, 1; 1 + (-1)*u rounds as 1 - u.
                basis[0] = c * (1.0 - u);
                basis[1] = c * (1.0 + u);
            }
            3 => {
                basis[0] = 1.0 - u - v;
                basis[1] = u;
                basis[2] = v;
            }
            4 => {
                let c = 1.0 / 4.0;
                for i in 0..4 {
                    basis[i] = c * (1.0 + QUAD_U[i] * u) * (1.0 + QUAD_V[i] * v);
                }
            }
            5 => {
                basis[0] = 1.0 - u - v - w;
                basis[1] = u;
                basis[2] = v;
                basis[3] = w;
            }
            7 => {
                let half = 1.0 / 2.0;
                let tri = [1.0 - u - v, u, v];
                let line = [half * (1.0 - w), half * (1.0 + w)];
                for i in 0..3 {
                    basis[i] = line[0] * tri[i];
                    basis[i + 3] = line[1] * tri[i];
                }
            }
            8 => {
                let c = 1.0 / 8.0;
                for i in 0..8 {
                    basis[i] = c
                        * (1.0 + BRICK_U[i] * u)
                        * (1.0 + BRICK_V[i] * v)
                        * (1.0 + BRICK_W[i] * w);
                }
            }
            _ => unreachable!("fast path family"),
        }
        return;
    }
    match def.dim {
        1 => {
            for (b, terms) in basis.iter_mut().zip(&def.basis) {
                let mut s = 0.0;
                for t in terms {
                    if t.p == 0 {
                        s += t.coeff;
                    } else {
                        s += t.coeff * powi(u, t.p as u32);
                    }
                }
                *b = s;
            }
        }
        2 => {
            let (ult, vlt) = (powers(u, def.basis_degree), powers(v, def.basis_degree));
            for (b, terms) in basis.iter_mut().zip(&def.basis) {
                let mut s = 0.0;
                for t in terms {
                    s += t.coeff * ult[t.p as usize] * vlt[t.q as usize];
                }
                *b = s;
            }
        }
        _ if def.family() == 6 => {
            let w = apex_guard(w);
            for (q, b) in basis.iter_mut().enumerate().take(n) {
                *b = pyramid_value(def, q, u, v, w);
            }
        }
        _ => {
            let d = def.basis_degree;
            let (ult, vlt, wlt) = (powers(u, d), powers(v, d), powers(w, d));
            for (b, terms) in basis.iter_mut().zip(&def.basis) {
                let mut s = 0.0;
                for t in terms {
                    s += t.coeff * ult[t.p as usize] * vlt[t.q as usize] * wlt[t.r as usize];
                }
                *b = s;
            }
        }
    }
}

/// Local gradients at `(u,v,w)` into `dbasis[0..n]`; only the first
/// `def.dim` components are written.
pub fn nodal_first_derivatives(
    def: &ElementTypeDef,
    u: f64,
    v: f64,
    w: f64,
    dbasis: &mut [[f64; 3]],
) {
    let n = def.n_nodes;
    if fast_path(def) {
        match def.family() {
            2 => {
                let c = 1.0 / 2.0;
                dbasis[0][0] = -c;
                dbasis[1][0] = c;
            }
            3 => {
                dbasis[0][0] = -1.0;
                dbasis[0][1] = -1.0;
                dbasis[1][0] = 1.0;
                dbasis[1][1] = 0.0;
                dbasis[2][0] = 0.0;
                dbasis[2][1] = 1.0;
            }
            4 => {
                let c = 1.0 / 4.0;
                for i in 0..4 {
                    dbasis[i][0] = c * QUAD_U[i] * (1.0 + QUAD_V[i] * v);
                    dbasis[i][1] = c * (1.0 + QUAD_U[i] * u) * QUAD_V[i];
                }
            }
            5 => {
                dbasis[0] = [-1.0, -1.0, -1.0];
                dbasis[1] = [1.0, 0.0, 0.0];
                dbasis[2] = [0.0, 1.0, 0.0];
                dbasis[3] = [0.0, 0.0, 1.0];
            }
            7 => {
                let half = 1.0 / 2.0;
                let tri = [1.0 - u - v, u, v];
                let line = [half * (1.0 - w), half * (1.0 + w)];
                let gradtri = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
                let gradline = [-half, half];
                for i in 0..3 {
                    for j in 0..2 {
                        dbasis[i][j] = gradtri[i][j] * line[0];
                        dbasis[i + 3][j] = gradtri[i][j] * line[1];
                    }
                    dbasis[i][2] = tri[i] * gradline[0];
                    dbasis[i + 3][2] = tri[i] * gradline[1];
                }
            }
            8 => {
                let c = 1.0 / 8.0;
                for i in 0..8 {
                    let (su, sv, sw) = (BRICK_U[i], BRICK_V[i], BRICK_W[i]);
                    dbasis[i][0] = c * su * (1.0 + sv * v) * (1.0 + sw * w);
                    dbasis[i][1] = c * (1.0 + su * u) * sv * (1.0 + sw * w);
                    dbasis[i][2] = c * (1.0 + su * u) * (1.0 + sv * v) * sw;
                }
            }
            _ => unreachable!("fast path family"),
        }
        return;
    }
    match def.dim {
        1 => {
            for (d, terms) in dbasis.iter_mut().zip(&def.basis) {
                let mut s = 0.0;
                for t in terms {
                    if t.p >= 1 {
                        s += t.p as f64 * t.coeff * powi(u, t.p as u32 - 1);
                    }
                }
                d[0] = s;
            }
        }
        2 => {
            let (ult, vlt) = (powers(u, def.basis_degree), powers(v, def.basis_degree));
            for (d, terms) in dbasis.iter_mut().zip(&def.basis) {
                let (mut s, mut t2) = (0.0, 0.0);
                for t in terms {
                    let (p, q) = (t.p as usize, t.q as usize);
                    if p >= 1 {
                        s += p as f64 * t.coeff * ult[p - 1] * vlt[q];
                    }
                    if q >= 1 {
                        t2 += q as f64 * t.coeff * ult[p] * vlt[q - 1];
                    }
                }
                d[0] = s;
                d[1] = t2;
            }
        }
        _ if def.family() == 6 => {
            let w = apex_guard(w);
            for (q, d) in dbasis.iter_mut().enumerate().take(n) {
                *d = pyramid_gradient(def, q, u, v, w);
            }
        }
        _ => {
            let deg = def.basis_degree;
            let (ult, vlt, wlt) = (powers(u, deg), powers(v, deg), powers(w, deg));
            for (d, terms) in dbasis.iter_mut().zip(&def.basis) {
                let (mut s, mut t2, mut z) = (0.0, 0.0, 0.0);
                for t in terms {
                    let (p, q, r) = (t.p as usize, t.q as usize, t.r as usize);
                    if p >= 1 {
                        s += p as f64 * t.coeff * ult[p - 1] * vlt[q] * wlt[r];
                    }
                    if q >= 1 {
                        t2 += q as f64 * t.coeff * ult[p] * vlt[q - 1] * wlt[r];
                    }
                    if r >= 1 {
                        z += r as f64 * t.coeff * ult[p] * vlt[q] * wlt[r - 1];
                    }
                }
                *d = [s, t2, z];
            }
        }
    }
}

/// `ult(0:6)`: 1, x, then `x**i` for `i = 2..=degree`.
fn powers(x: f64, degree: usize) -> [f64; 7] {
    let mut lt = [0.0; 7];
    lt[0] = 1.0;
    lt[1] = x;
    for (i, l) in lt.iter_mut().enumerate().take(degree + 1).skip(2) {
        *l = powi(x, i as u32);
    }
    lt
}

/// `InterpolateInElement3D` generic branch for one node (x(q) = 1).
fn interp_terms(terms: &[BasisTerm], u: f64, v: f64, w: f64) -> f64 {
    let mut s = 0.0;
    for t in terms {
        s += t.coeff * powi(u, t.p as u32) * powi(v, t.q as u32) * powi(w, t.r as u32);
    }
    0.0 + s * 1.0
}

/// `InterpolateInElement3D(element, e_q, u, v, w)` for a pyramid; `w`
/// has been through [`apex_guard`].
fn pyramid_value(def: &ElementTypeDef, q: usize, u: f64, v: f64, w: f64) -> f64 {
    let s = 1.0 / (1.0 - w);
    match def.code {
        605 => {
            let t = match q {
                0 => ((1.0 - u) * (1.0 - v) - w + u * v * w * s) / 4.0,
                1 => ((1.0 + u) * (1.0 - v) - w - u * v * w * s) / 4.0,
                2 => ((1.0 + u) * (1.0 + v) - w + u * v * w * s) / 4.0,
                3 => ((1.0 - u) * (1.0 + v) - w - u * v * w * s) / 4.0,
                _ => w,
            };
            0.0 + t
        }
        613 => {
            let t = match q {
                0 => (-u - v - 1.0) * ((1.0 - u) * (1.0 - v) - w + u * v * w * s) / 4.0,
                1 => (u - v - 1.0) * ((1.0 + u) * (1.0 - v) - w - u * v * w * s) / 4.0,
                2 => (u + v - 1.0) * ((1.0 + u) * (1.0 + v) - w + u * v * w * s) / 4.0,
                3 => (-u + v - 1.0) * ((1.0 - u) * (1.0 + v) - w - u * v * w * s) / 4.0,
                4 => w * (2.0 * w - 1.0),
                5 => (1.0 + u - w) * (1.0 - u - w) * (1.0 - v - w) * s / 2.0,
                6 => (1.0 + v - w) * (1.0 - v - w) * (1.0 + u - w) * s / 2.0,
                7 => (1.0 + u - w) * (1.0 - u - w) * (1.0 + v - w) * s / 2.0,
                8 => (1.0 + v - w) * (1.0 - v - w) * (1.0 - u - w) * s / 2.0,
                9 => w * (1.0 - u - w) * (1.0 - v - w) * s,
                10 => w * (1.0 + u - w) * (1.0 - v - w) * s,
                11 => w * (1.0 + u - w) * (1.0 + v - w) * s,
                _ => w * (1.0 - u - w) * (1.0 + v - w) * s,
            };
            0.0 + t
        }
        _ => interp_terms(&def.basis[q], u, v, w),
    }
}

/// `FirstDerivativeIn{U,V,W}3D(element, e_q, u, v, w)` for a pyramid.
fn pyramid_gradient(def: &ElementTypeDef, q: usize, u: f64, v: f64, w: f64) -> [f64; 3] {
    let s = 1.0 / (1.0 - w);
    match def.code {
        605 => match q {
            0 => [
                0.0 + (-(1.0 - v) + v * w * s) / 4.0,
                0.0 + (-(1.0 - u) + u * w * s) / 4.0,
                0.0 + (-1.0 + u * v * (s * s)) / 4.0,
            ],
            1 => [
                0.0 + ((1.0 - v) - v * w * s) / 4.0,
                0.0 + (-(1.0 + u) - u * w * s) / 4.0,
                0.0 + (-1.0 - u * v * (s * s)) / 4.0,
            ],
            2 => [
                0.0 + ((1.0 + v) + v * w * s) / 4.0,
                0.0 + ((1.0 + u) + u * w * s) / 4.0,
                0.0 + (-1.0 + u * v * (s * s)) / 4.0,
            ],
            3 => [
                0.0 + (-(1.0 + v) - v * w * s) / 4.0,
                0.0 + ((1.0 - u) - u * w * s) / 4.0,
                0.0 + (-1.0 - u * v * (s * s)) / 4.0,
            ],
            _ => [0.0, 0.0, 0.0 + 1.0],
        },
        613 => pyramid13_gradient(q, u, v, w, s),
        _ => {
            let terms = &def.basis[q];
            let (mut du, mut dv, mut dw) = (0.0, 0.0, 0.0);
            for t in terms {
                let (p, q, r) = (t.p as u32, t.q as u32, t.r as u32);
                if p >= 1 {
                    du += p as f64 * t.coeff * powi(u, p - 1) * powi(v, q) * powi(w, r);
                }
            }
            for t in terms {
                let (p, q, r) = (t.p as u32, t.q as u32, t.r as u32);
                if q >= 1 {
                    dv += q as f64 * t.coeff * powi(u, p) * powi(v, q - 1) * powi(w, r);
                }
            }
            for t in terms {
                let (p, q, r) = (t.p as u32, t.q as u32, t.r as u32);
                if r >= 1 {
                    dw += r as f64 * t.coeff * powi(u, p) * powi(v, q) * powi(w, r - 1);
                }
            }
            [0.0 + du * 1.0, 0.0 + dv * 1.0, 0.0 + dw * 1.0]
        }
    }
}

/// The 613 branches of `FirstDerivativeIn{U,V,W}3D`.
fn pyramid13_gradient(q: usize, u: f64, v: f64, w: f64, s: f64) -> [f64; 3] {
    let s2 = s * s;
    match q {
        0 => [
            0.0 + (-((1.0 - u) * (1.0 - v) - w + u * v * w * s)
                + (-u - v - 1.0) * (-(1.0 - v) + v * w * s))
                / 4.0,
            0.0 + (-((1.0 - u) * (1.0 - v) - w + u * v * w * s)
                + (-u - v - 1.0) * (-(1.0 - u) + u * w * s))
                / 4.0,
            0.0 + (-u - v - 1.0) * (-1.0 + u * v * s2) / 4.0,
        ],
        1 => [
            0.0 + (((1.0 + u) * (1.0 - v) - w - u * v * w * s)
                + (u - v - 1.0) * ((1.0 - v) - v * w * s))
                / 4.0,
            0.0 + (-((1.0 + u) * (1.0 - v) - w - u * v * w * s)
                + (u - v - 1.0) * (-(1.0 + u) - u * w * s))
                / 4.0,
            0.0 + (u - v - 1.0) * (-1.0 - u * v * s2) / 4.0,
        ],
        2 => [
            0.0 + (((1.0 + u) * (1.0 + v) - w + u * v * w * s)
                + (u + v - 1.0) * ((1.0 + v) + v * w * s))
                / 4.0,
            0.0 + (((1.0 + u) * (1.0 + v) - w + u * v * w * s)
                + (u + v - 1.0) * ((1.0 + u) + u * w * s))
                / 4.0,
            0.0 + (u + v - 1.0) * (-1.0 + u * v * s2) / 4.0,
        ],
        3 => [
            0.0 + (-((1.0 - u) * (1.0 + v) - w - u * v * w * s)
                + (-u + v - 1.0) * (-(1.0 + v) - v * w * s))
                / 4.0,
            0.0 + (((1.0 - u) * (1.0 + v) - w - u * v * w * s)
                + (-u + v - 1.0) * ((1.0 - u) - u * w * s))
                / 4.0,
            0.0 + (-u + v - 1.0) * (-1.0 - u * v * s2) / 4.0,
        ],
        4 => [0.0, 0.0, 0.0 + (4.0 * w - 1.0)],
        5 => [
            0.0 + ((1.0 - u - w) * (1.0 - v - w) - (1.0 + u - w) * (1.0 - v - w)) * s / 2.0,
            0.0 - (1.0 + u - w) * (1.0 - u - w) * s / 2.0,
            0.0 + ((-((1.0 - u - w) * (1.0 - v - w))
                - (1.0 + u - w) * (1.0 - v - w)
                - (1.0 + u - w) * (1.0 - u - w))
                * s
                + (1.0 + u - w) * (1.0 - u - w) * (1.0 - v - w) * s2)
                / 2.0,
        ],
        6 => [
            0.0 + ((1.0 + v - w) * (1.0 - v - w)) * s / 2.0,
            0.0 + ((1.0 - v - w) * (1.0 + u - w) - (1.0 + v - w) * (1.0 + u - w)) * s / 2.0,
            0.0 + ((-((1.0 - v - w) * (1.0 + u - w))
                - (1.0 + v - w) * (1.0 + u - w)
                - (1.0 + v - w) * (1.0 - v - w))
                * s
                + (1.0 + v - w) * (1.0 - v - w) * (1.0 + u - w) * s2)
                / 2.0,
        ],
        7 => [
            0.0 + ((1.0 - u - w) * (1.0 + v - w) - (1.0 + u - w) * (1.0 + v - w)) * s / 2.0,
            0.0 + (1.0 + u - w) * (1.0 - u - w) * s / 2.0,
            0.0 + ((-((1.0 - u - w) * (1.0 + v - w))
                - (1.0 + u - w) * (1.0 + v - w)
                - (1.0 + u - w) * (1.0 - u - w))
                * s
                + (1.0 + u - w) * (1.0 - u - w) * (1.0 + v - w) * s2)
                / 2.0,
        ],
        8 => [
            0.0 + (-((1.0 + v - w) * (1.0 - v - w))) * s / 2.0,
            0.0 + ((1.0 - v - w) * (1.0 - u - w) - (1.0 + v - w) * (1.0 - u - w)) * s / 2.0,
            0.0 + ((-((1.0 - v - w) * (1.0 - u - w))
                - (1.0 + v - w) * (1.0 - u - w)
                - (1.0 + v - w) * (1.0 - v - w))
                * s
                + (1.0 + v - w) * (1.0 - v - w) * (1.0 - u - w) * s2)
                / 2.0,
        ],
        9 => [
            0.0 - w * (1.0 - v - w) * s,
            0.0 - w * (1.0 - u - w) * s,
            0.0 + (((1.0 - u - w) * (1.0 - v - w) - w * (1.0 - v - w) - w * (1.0 - u - w)) * s
                + w * (1.0 - u - w) * (1.0 - v - w) * s2),
        ],
        10 => [
            0.0 + w * (1.0 - v - w) * s,
            0.0 - w * (1.0 + u - w) * s,
            0.0 + (((1.0 + u - w) * (1.0 - v - w) - w * (1.0 - v - w) - w * (1.0 + u - w)) * s
                + w * (1.0 + u - w) * (1.0 - v - w) * s2),
        ],
        11 => [
            0.0 + w * (1.0 + v - w) * s,
            0.0 + w * (1.0 + u - w) * s,
            0.0 + (((1.0 + u - w) * (1.0 + v - w) - w * (1.0 + v - w) - w * (1.0 + u - w)) * s
                + w * (1.0 + u - w) * (1.0 + v - w) * s2),
        ],
        _ => [
            0.0 - w * (1.0 + v - w) * s,
            0.0 + w * (1.0 - u - w) * s,
            0.0 + (((1.0 - u - w) * (1.0 + v - w) - w * (1.0 + v - w) - w * (1.0 - u - w)) * s
                + w * (1.0 - u - w) * (1.0 + v - w) * s2),
        ],
    }
}
