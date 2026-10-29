//! Second derivatives of the nodal basis: `EvalSecondDerivativesRef`
//! (ElemInfo.F90:528) with `SecondDerivatives1D/2D/3D` (ElementBasis.F90) on
//! the reference element, and `GlobalSecondDerivatives` (ElemInfo.F90:10663)
//! mapping them to physical coordinates.
//!
//! Like Elmer, every non-p type differentiates its elements.def monomials
//! (also the linear types whose values take the closed-form fast path); the
//! 605 and 613 pyramids use their rational forms.

use crate::element_type::{powi, ElementTypeDef};

/// `ddLBasisddx(q,i,j)` for every node `q` at `(u,v,w)`: `out[q][i][j]`,
/// zero outside `dim x dim`.
pub fn local_second_derivatives(def: &ElementTypeDef, uvw: [f64; 3]) -> Vec<[[f64; 3]; 3]> {
    let [u, v, w] = uvw;
    let n = def.n_nodes;
    let mut out = vec![[[0.0; 3]; 3]; n];
    for (q, dd) in out.iter_mut().enumerate() {
        match def.dim {
            1 => dd[0][0] = second_1d(def, q, u),
            2 => {
                let d = second_2d(def, q, u, v);
                for i in 0..2 {
                    dd[i][..2].copy_from_slice(&d[i]);
                }
            }
            _ => *dd = second_3d(def, q, u, v, w),
        }
    }
    out
}

/// `SecondDerivatives1D( element, e_q, u )`.
fn second_1d(def: &ElementTypeDef, q: usize, u: f64) -> f64 {
    let mut s = 0.0;
    for t in &def.basis[q] {
        if t.p >= 2 {
            let p = t.p as i32;
            s += (p * (p - 1)) as f64 * t.coeff * powi(u, t.p as u32 - 2);
        }
    }
    0.0 + s * 1.0
}

/// `SecondDerivatives2D( element, e_q, u, v )`.
fn second_2d(def: &ElementTypeDef, q: usize, u: f64, v: f64) -> [[f64; 2]; 2] {
    let terms = &def.basis[q];
    let (mut s11, mut s12, mut s22) = (0.0, 0.0, 0.0);
    for t in terms {
        let (p, r) = (t.p as u32, t.q as u32);
        if p >= 2 {
            s11 += ((p * (p - 1)) as f64) * t.coeff * powi(u, p - 2) * powi(v, r);
        }
    }
    for t in terms {
        let (p, r) = (t.p as u32, t.q as u32);
        if p >= 1 && r >= 1 {
            s12 += ((p * r) as f64) * t.coeff * powi(u, p - 1) * powi(v, r - 1);
        }
    }
    for t in terms {
        let (p, r) = (t.p as u32, t.q as u32);
        if r >= 2 {
            s22 += ((r * (r - 1)) as f64) * t.coeff * powi(u, p) * powi(v, r - 2);
        }
    }
    let d12 = 0.0 + s12 * 1.0;
    [[0.0 + s11 * 1.0, d12], [d12, 0.0 + s22 * 1.0]]
}

/// `SecondDerivatives3D( element, e_q, u, v, w )`.
fn second_3d(def: &ElementTypeDef, q: usize, u: f64, v: f64, w: f64) -> [[f64; 3]; 3] {
    match def.code {
        605 => {
            let w = if w == 1.0 { 1.0 - 1.0e-12 } else { w };
            let s = 1.0 / (1.0 - w);
            let mut x = [0.0f64; 5];
            x[q] = 1.0;
            let a = x[0] - x[1] + x[2] - x[3];
            let mut ddx = [[0.0f64; 3]; 3];
            ddx[0][1] = a * (1.0 + w * s);
            ddx[1][0] = ddx[0][1];
            ddx[0][2] = a * v * (s * s);
            ddx[2][0] = ddx[0][2];
            ddx[1][2] = a * u * (s * s);
            ddx[2][1] = ddx[1][2];
            ddx.map(|r| r.map(|d| d / 4.0))
        }
        613 => {
            let w = if w == 1.0 { 1.0 - 1.0e-12 } else { w };
            pyramid13_second(q, u, v, w, 1.0 / (1.0 - w))
        }
        _ => {
            let terms = &def.basis[q];
            let pw = |x: f64, k: u8, d: u8| powi(x, (k - d) as u32);
            let mut ddx = [[0.0f64; 3]; 3];
            let mut s = 0.0;
            for t in terms {
                if t.p >= 2 {
                    let c = (t.p as u32 * (t.p as u32 - 1)) as f64;
                    s += c * t.coeff * pw(u, t.p, 2) * pw(v, t.q, 0) * pw(w, t.r, 0);
                }
            }
            ddx[0][0] = 0.0 + s * 1.0;
            s = 0.0;
            for t in terms {
                if t.p >= 1 && t.q >= 1 {
                    let c = (t.p as u32 * t.q as u32) as f64;
                    s += c * t.coeff * pw(u, t.p, 1) * pw(v, t.q, 1) * pw(w, t.r, 0);
                }
            }
            ddx[0][1] = 0.0 + s * 1.0;
            // Elmer's loop for (1,3) runs `DO i = 2,k` and so skips term 1.
            s = 0.0;
            for t in terms.iter().take(def.n_nodes).skip(1) {
                if t.p >= 1 && t.r >= 1 {
                    let c = (t.p as u32 * t.r as u32) as f64;
                    s += c * t.coeff * pw(u, t.p, 1) * pw(v, t.q, 0) * pw(w, t.r, 1);
                }
            }
            ddx[0][2] = 0.0 + s * 1.0;
            s = 0.0;
            for t in terms {
                if t.q >= 2 {
                    let c = (t.q as u32 * (t.q as u32 - 1)) as f64;
                    s += c * t.coeff * pw(u, t.p, 0) * pw(v, t.q, 2) * pw(w, t.r, 0);
                }
            }
            ddx[1][1] = 0.0 + s * 1.0;
            s = 0.0;
            for t in terms {
                if t.q >= 1 && t.r >= 1 {
                    let c = (t.q as u32 * t.r as u32) as f64;
                    s += c * t.coeff * pw(u, t.p, 0) * pw(v, t.q, 1) * pw(w, t.r, 1);
                }
            }
            ddx[1][2] = 0.0 + s * 1.0;
            s = 0.0;
            for t in terms {
                if t.r >= 2 {
                    let c = (t.r as u32 * (t.r as u32 - 1)) as f64;
                    s += c * t.coeff * pw(u, t.p, 0) * pw(v, t.q, 0) * pw(w, t.r, 2);
                }
            }
            ddx[2][2] = 0.0 + s * 1.0;
            ddx[1][0] = ddx[0][1];
            ddx[2][0] = ddx[0][2];
            ddx[2][1] = ddx[1][2];
            ddx
        }
    }
}

/// `GlobalSecondDerivatives( elm, nodes, values, u, v, w, Metric, dLBasisdx,
/// ddLBasisddx, nd )`: `out[q][i][j]` for `i, j < cdim`.
#[allow(clippy::needless_range_loop)]
pub fn global_second_derivatives(
    dim: usize,
    cdim: usize,
    xyz: &[[f64; 3]],
    dl: &[[f64; 3]],
    ddl: &[[[f64; 3]; 3]],
    metric: &[[f64; 3]; 3],
    out: &mut [[[f64; 3]; 3]],
) {
    let nd = xyz.len();
    let mut dx = [[0.0f64; 3]; 3];
    for c in 0..cdim.min(3) {
        for i in 0..dim {
            let mut s = 0.0;
            for q in 0..nd {
                s += xyz[q][c] * dl[q][i];
            }
            dx[c][i] = s;
        }
    }
    let mut ddx = [[[0.0f64; 3]; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            for c in 0..3 {
                let mut s = 0.0;
                for q in 0..nd {
                    s += ddl[q][i][j] * xyz[q][c];
                }
                ddx[c][i][j] = s;
            }
        }
    }
    let mut c2 = [[[0.0f64; 3]; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let mut s = 0.0;
                for l in 0..cdim {
                    s += ddx[l][i][j] * dx[l][k];
                }
                c2[i][j][k] = s;
            }
        }
    }
    let mut c1 = [[[0.0f64; 3]; 3]; 3];
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let mut s = 0.0;
                for l in 0..dim {
                    s += metric[k][l] * c2[i][j][l];
                }
                c1[i][j][k] = s;
            }
        }
    }
    for q in 0..nd {
        let df = dl[q];
        let mut ddf = ddl[q];
        for i in 0..dim {
            for j in 0..dim {
                let mut s = 0.0;
                for k in 0..dim {
                    s -= c1[i][j][k] * df[k];
                }
                ddf[i][j] += s;
            }
        }
        let mut cddf = [[0.0f64; 3]; 3];
        for i in 0..dim {
            for j in 0..dim {
                let mut s = 0.0;
                for k in 0..dim {
                    for l in 0..dim {
                        s += metric[i][k] * metric[j][l] * ddf[k][l];
                    }
                }
                cddf[i][j] = s;
            }
        }
        out[q] = [[0.0; 3]; 3];
        for i in 0..cdim {
            for j in 0..cdim {
                let mut s = 0.0;
                for k in 0..dim {
                    for l in 0..dim {
                        s += dx[i][k] * dx[j][l] * cddf[k][l];
                    }
                }
                out[q][i][j] = s;
            }
        }
    }
}
/// The 613 branch of `SecondDerivatives3D` for `x = e_q` (generated from
/// ElementBasis.F90 by the lane tool gen_613_dd.py; each `x(N) *` factor is 1).
#[rustfmt::skip]
#[allow(clippy::all)]
fn pyramid13_second(q: usize, u: f64, v: f64, w: f64, s: f64) -> [[f64; 3]; 3] {
    let (s2, s3) = (s * s, powi(s, 3));
    let mut ddx = [[0.0f64; 3]; 3];
    let n = q + 1;
    let mut t = 0.0f64;
    match n {
        1 => { t = t - (-(1.0-v) + v*w*s)/2.0; }
        2 => { t = t + ( (1.0-v) - v*w*s)/2.0; }
        3 => { t = t + ( (1.0+v) + v*w*s)/2.0; }
        4 => { t = t - (-(1.0+v) - v*w*s)/2.0; }
        6 => { t = t - (1.0-v-w) * s; }
        8 => { t = t - (1.0+v-w) * s; }
        _ => {}
    }
    ddx[0][0] = 0.0 + t;
    let mut t = 0.0f64;
    match n {
        1 => { t = t + ((1.0-u) - u*w*s)/4.0; t = t + ((1.0-v) - v*w*s)/4.0; t = t + (-u-v-1.0)*(1.0+w*s)/4.0; }
        2 => { t = t + (-(1.0+u) - u*w*s)/4.0; t = t + ( -(1.0-v) + v*w*s)/4.0; t = t + ( u-v-1.0)*(-1.0-w*s)/4.0; }
        3 => { t = t + ( (1.0+u) + u*w*s)/4.0; t = t + ( (1.0+v) + v*w*s)/4.0; t = t + ( u+v-1.0)*(1.0+w*s)/4.0; }
        4 => { t = t + ( -(1.0-u) + u*w*s)/4.0; t = t + (-(1.0+v) - v*w*s)/4.0; t = t + (-u+v-1.0)*(-1.0-w*s)/4.0; }
        6 => { t = t - (1.0-u-w)*s/2.0; t = t + (1.0+u-w)*s/2.0; }
        7 => { t = t + (1.0-v-w)*s/2.0; t = t - (1.0+v-w)*s/2.0; }
        8 => { t = t + (1.0-u-w)*s/2.0; t = t - (1.0+u-w)*s/2.0; }
        9 => { t = t - (1.0-v-w)*s/2.0; t = t + (1.0+v-w)*s/2.0; }
        10 => { t = t + w*s; }
        11 => { t = t - w*s; }
        12 => { t = t + w*s; }
        13 => { t = t - w*s; }
        _ => {}
    }
    ddx[0][1] = 0.0 + t;
    let mut t = 0.0f64;
    match n {
        1 => { t = t - (-1.0 + u*v*s2) / 4.0; t = t + (-u-v-1.0) * (v*s2) / 4.0; }
        2 => { t = t + (-1.0 - u*v*s2) / 4.0; t = t + ( u-v-1.0) * (-v*s2) / 4.0; }
        3 => { t = t + (-1.0 + u*v*s2) / 4.0; t = t + ( u+v-1.0) * (v*s2) / 4.0; }
        4 => { t = t - (-1.0 - u*v*s2) / 4.0; t = t + (-u+v-1.0) * (-v*s2) / 4.0; }
        6 => { t = t - (1.0-v-w) * s / 2.0; t = t - (1.0-u-w) * s / 2.0; t = t + (1.0-u-w)*(1.0-v-w) * s2 / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0-v-w) * s2 / 2.0; }
        7 => { t = t - (1.0-v-w) * s / 2.0; t = t - (1.0+v-w) * s / 2.0; t = t + (1.0+v-w)*(1.0-v-w) * s2 / 2.0; }
        8 => { t = t - (1.0+v-w) * s / 2.0; t = t - (1.0-u-w) * s / 2.0; t = t + (1.0-u-w)*(1.0+v-w) * s2 / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0+v-w) * s2 / 2.0; }
        9 => { t = t + (1.0-v-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0-v-w) * s2 / 2.0; }
        10 => { t = t + w * s; t = t - (1.0-v-w) * s2; }
        11 => { t = t - w * s; t = t + (1.0-v-w) * s2; }
        12 => { t = t - w * s; t = t + (1.0+v-w) * s2; }
        13 => { t = t + w * s; t = t - (1.0+v-w) * s2; }
        _ => {}
    }
    ddx[0][2] = 0.0 + t;
    let mut t = 0.0f64;
    match n {
        1 => { t = t - (-(1.0-u) + u*w*s)/2.0; }
        2 => { t = t - (-(1.0+u) - u*w*s)/2.0; }
        3 => { t = t + ( (1.0+u) + u*w*s)/2.0; }
        4 => { t = t + ( (1.0-u) - u*w*s)/2.0; }
        7 => { t = t - (1.0+u-w)*s; }
        9 => { t = t - (1.0-u-w)*s; }
        _ => {}
    }
    ddx[1][1] = 0.0 + t;
    let mut t = 0.0f64;
    match n {
        1 => { t = t - (-1.0 + u*v*s2) / 4.0; t = t + (-u-v-1.0) * (u*s2) / 4.0; }
        2 => { t = t - (-1.0 - u*v*s2) / 4.0; t = t + ( u-v-1.0) * (-u*s2) / 4.0; }
        3 => { t = t + (-1.0 + u*v*s2) / 4.0; t = t + ( u+v-1.0) * (u*s2) / 4.0; }
        4 => { t = t + (-1.0 - u*v*s2) / 4.0; t = t + (-u+v-1.0) * (-u*s2) / 4.0; }
        6 => { t = t + (1.0-u-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0-u-w) * s2 / 2.0; }
        7 => { t = t - (1.0+u-w) * s / 2.0; t = t - (1.0-v-w) * s / 2.0; t = t + (1.0-v-w)*(1.0+u-w) * s2 / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0+u-w) * s2 / 2.0; }
        8 => { t = t - (1.0-u-w) * s / 2.0; t = t - (1.0+u-w) * s / 2.0; t = t + (1.0+u-w)*(1.0-u-w) * s2 / 2.0; }
        9 => { t = t - (1.0-u-w) * s / 2.0; t = t - (1.0-v-w) * s / 2.0; t = t + (1.0-v-w)*(1.0-u-w) * s2 / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0-u-w) * s2 / 2.0; }
        10 => { t = t + w * s; t = t - (1.0-u-w) * s2; }
        11 => { t = t + w * s; t = t - (1.0+u-w) * s2; }
        12 => { t = t - w * s; t = t + (1.0+u-w) * s2; }
        13 => { t = t - w * s; t = t + (1.0-u-w) * s2; }
        _ => {}
    }
    ddx[1][2] = 0.0 + t;
    let mut t = 0.0f64;
    match n {
        1 => { t = t + (-u-v-1.0) * ( u*v*2.0*s3) / 4.0; }
        2 => { t = t + ( u-v-1.0) * (-u*v*2.0*s3) / 4.0; }
        3 => { t = t + ( u+v-1.0) * ( u*v*2.0*s3) / 4.0; }
        4 => { t = t + (-u+v-1.0) * (-u*v*2.0*s3) / 4.0; }
        5 => { t = t + 4.0; }
        6 => { t = t + (1.0-v-w) * s / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t - (1.0-u-w)*(1.0-v-w) * s2 / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0-v-w) * s2 / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0-u-w) * s2 / 2.0; t = t - (1.0-u-w)*(1.0-v-w) * s2 / 2.0; t = t - (1.0+u-w)*(1.0-v-w) * s2 / 2.0; t = t - (1.0+u-w)*(1.0-u-w) * s2 / 2.0; t = t + (1.0+u-w)*(1.0-u-w)*(1.0-v-w) * 2.0*s3 / 2.0; }
        7 => { t = t + (1.0+u-w) * s / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t - (1.0-v-w)*(1.0+u-w) * s2 / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0+u-w) * s2 / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0-v-w) * s2 / 2.0; t = t - (1.0-v-w)*(1.0+u-w) * s2 / 2.0; t = t - (1.0+v-w)*(1.0+u-w) * s2 / 2.0; t = t - (1.0+v-w)*(1.0-v-w) * s2 / 2.0; t = t + (1.0+v-w)*(1.0-v-w)*(1.0+u-w) * 2.0*s3 / 2.0; }
        8 => { t = t + (1.0+v-w) * s / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t - (1.0-u-w)*(1.0+v-w) * s2 / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0+v-w) * s2 / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t + (1.0+u-w) * s / 2.0; t = t - (1.0+u-w)*(1.0-u-w) * s2 / 2.0; t = t - (1.0-u-w)*(1.0+v-w) * s2 / 2.0; t = t - (1.0+u-w)*(1.0+v-w) * s2 / 2.0; t = t - (1.0+u-w)*(1.0-u-w) * s2 / 2.0; t = t + (1.0+u-w)*(1.0-u-w)*(1.0+v-w) * 2.0*s3 / 2.0; }
        9 => { t = t + (1.0-u-w) * s / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t - (1.0-v-w)*(1.0-u-w) * s2 / 2.0; t = t + (1.0-u-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0-u-w) * s2 / 2.0; t = t + (1.0-v-w) * s / 2.0; t = t + (1.0+v-w) * s / 2.0; t = t - (1.0+v-w)*(1.0-v-w) * s2 / 2.0; t = t - (1.0-v-w)*(1.0-u-w) * s2 / 2.0; t = t - (1.0+v-w)*(1.0-u-w) * s2 / 2.0; t = t - (1.0+v-w)*(1.0-v-w) * s2 / 2.0; t = t + (1.0+v-w)*(1.0-v-w)*(1.0-u-w) * 2.0*s3 / 2.0; }
        10 => { t = t + w * s; t = t - (1.0-v-w) * s2; t = t + w * s; t = t - (1.0-u-w) * s2; t = t - (1.0-v-w) * s2; t = t - (1.0-u-w) * s2; t = t + (1.0-u-w) * (1.0-v-w) * 2.0*s3; }
        11 => { t = t + w * s; t = t - (1.0-v-w) * s2; t = t + w * s; t = t - (1.0+u-w) * s2; t = t - (1.0-v-w) * s2; t = t - (1.0+u-w) * s2; t = t + (1.0+u-w) * (1.0-v-w) * 2.0*s3; }
        12 => { t = t + w * s; t = t - (1.0+v-w) * s2; t = t + w * s; t = t - (1.0+u-w) * s2; t = t - (1.0+v-w) * s2; t = t - (1.0+u-w) * s2; t = t + (1.0+u-w) * (1.0+v-w) * 2.0*s3; }
        13 => { t = t + w*s; t = t - (1.0+v-w) * s2; t = t + w*s; t = t - (1.0-u-w) * s2; t = t - (1.0+v-w) * s2; t = t - (1.0-u-w) * s2; t = t + (1.0-u-w) * (1.0+v-w) * 2.0*s3; }
        _ => {}
    }
    ddx[2][2] = 0.0 + t;
    ddx[1][0] = ddx[0][1];
    ddx[2][0] = ddx[0][2];
    ddx[2][1] = ddx[1][2];
    ddx
}
