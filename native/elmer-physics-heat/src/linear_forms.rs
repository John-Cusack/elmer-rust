//! The bilinear and linear forms of Elmer's `LinearForms.F90`, in its exact
//! summation order, so the vectorised heat kernel rounds as
//! `LocalMatrixVec` does.
//!
//! Points are processed in blocks of [`BLOCK`] (`VECTOR_BLOCK_LENGTH`). A
//! block shorter than [`SMALL`] (`VECTOR_SMALL_THRESH`) uses the plain loops,
//! adding each point's term into the entry directly, `k` (the derivative
//! direction) outer and the point inner. A longer block goes through BLAS;
//! the order here is the reference (netlib) `DGEMM`/`DGEMV` one: a dot product
//! over the block's points from zero, then added to the entry.
//!
//! Layouts: point values `u[t*n + p]` (`Basis(t,p)`), gradients
//! `du[(k*ngp + t)*n + p]` (`dBasisdx(t,p,k)`), matrices row-major
//! `g[i*n + j]` (`G(i,j)`), and a per-point vector field `beta[t][k]`.
#![allow(clippy::needless_range_loop)]

/// `VECTOR_BLOCK_LENGTH` (Types.F90).
pub(crate) const BLOCK: usize = 128;
/// `VECTOR_SMALL_THRESH` (Types.F90).
pub(crate) const SMALL: usize = 9;

fn blocks(m: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..m).step_by(BLOCK).map(move |a| (a, (a + BLOCK).min(m)))
}

/// Reference `DGEMM('T', 'N', n, n, b-a, 1, A(a,1), .., wrk, .., 1, G, ..)`:
/// `G(i,j) = sum_l A(l,i)*wrk(l,j) + G(i,j)`, with `wrk` holding the block's
/// rows only.
fn dgemm_tn(n: usize, a: usize, b: usize, amat: &[f64], wrk: &[f64], g: &mut [f64]) {
    for j in 0..n {
        for i in 0..n {
            let mut temp = 0.0;
            for l in a..b {
                temp += amat[l * n + i] * wrk[(l - a) * n + j];
            }
            g[i * n + j] += temp;
        }
    }
}

/// `LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)`:
/// `G(i,j) += (alpha grad u_j, grad u_i)`.
pub(crate) fn grad_u_dot_grad_u(
    ngp: usize,
    n: usize,
    dim: usize,
    du: &[f64],
    weight: &[f64],
    g: &mut [f64],
    alpha: &[f64],
) {
    let d = |k: usize| &du[k * ngp * n..(k + 1) * ngp * n];
    for (a, b) in blocks(ngp) {
        if b - a < SMALL {
            for i in 0..n {
                for j in 0..n {
                    let mut s = g[i * n + j];
                    for k in 0..dim {
                        let dk = d(k);
                        for l in a..b {
                            s += dk[l * n + i] * dk[l * n + j] * weight[l] * alpha[l];
                        }
                    }
                    g[i * n + j] = s;
                }
            }
        } else {
            let mut wrk = vec![0.0; (b - a) * n];
            for k in 0..dim {
                let dk = d(k);
                for l in a..b {
                    for j in 0..n {
                        wrk[(l - a) * n + j] = weight[l] * alpha[l] * dk[l * n + j];
                    }
                }
                dgemm_tn(n, a, b, dk, &wrk, g);
            }
        }
    }
}

/// `LinearForms_GradUdotU(m, n, dim, GradU, U, weight, G, alpha, beta)`:
/// `G(i,j) += (alpha beta.grad u_j, u_i)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn grad_u_dot_u(
    ngp: usize,
    n: usize,
    dim: usize,
    du: &[f64],
    u: &[f64],
    weight: &[f64],
    g: &mut [f64],
    alpha: &[f64],
    beta: &[[f64; 3]],
) {
    let d = |k: usize| &du[k * ngp * n..(k + 1) * ngp * n];
    for (a, b) in blocks(ngp) {
        if b - a < SMALL {
            for i in 0..n {
                for j in 0..n {
                    let mut s = g[i * n + j];
                    for k in 0..dim {
                        let dk = d(k);
                        for l in a..b {
                            s += dk[l * n + j] * u[l * n + i] * weight[l] * alpha[l] * beta[l][k];
                        }
                    }
                    g[i * n + j] = s;
                }
            }
        } else {
            let mut wrk = vec![0.0; (b - a) * n];
            for k in 0..dim {
                let dk = d(k);
                for l in a..b {
                    for j in 0..n {
                        wrk[(l - a) * n + j] = weight[l] * alpha[l] * beta[l][k] * dk[l * n + j];
                    }
                }
                dgemm_tn(n, a, b, u, &wrk, g);
            }
        }
    }
}

/// `LinearForms_UdotV(m, n, dim, U, V, weight, G[, alpha])`, and with `v = u`
/// `LinearForms_UdotU`, whose loops are the same: `G(i,j) += (alpha v_j, u_i)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn u_dot_v(
    ngp: usize,
    n: usize,
    u: &[f64],
    v: &[f64],
    weight: &[f64],
    g: &mut [f64],
    alpha: Option<&[f64]>,
) {
    for (a, b) in blocks(ngp) {
        if b - a < SMALL {
            for i in 0..n {
                for j in 0..n {
                    let mut s = g[i * n + j];
                    for l in a..b {
                        s += match alpha {
                            None => u[l * n + i] * v[l * n + j] * weight[l],
                            Some(al) => u[l * n + i] * v[l * n + j] * weight[l] * al[l],
                        };
                    }
                    g[i * n + j] = s;
                }
            }
        } else {
            let mut wrk = vec![0.0; (b - a) * n];
            for l in a..b {
                for j in 0..n {
                    wrk[(l - a) * n + j] = match alpha {
                        None => weight[l] * v[l * n + j],
                        Some(al) => weight[l] * al[l] * v[l * n + j],
                    };
                }
            }
            dgemm_tn(n, a, b, u, &wrk, g);
        }
    }
}

/// `LinearForms_UdotF(m, n, U, weight, F, UdotF)`: `out(i) += (f, u_i)`.
pub(crate) fn u_dot_f(ngp: usize, n: usize, u: &[f64], weight: &[f64], f: &[f64], out: &mut [f64]) {
    for (a, b) in blocks(ngp) {
        if b - a < SMALL {
            for i in 0..n {
                let mut s = out[i];
                for l in a..b {
                    s += u[l * n + i] * f[l] * weight[l];
                }
                out[i] = s;
            }
        } else {
            // Reference DGEMV('T', ..., 1, U, .., wrk, 1, 1, UdotF, 1).
            let wrk: Vec<f64> = (a..b).map(|l| weight[l] * f[l]).collect();
            for i in 0..n {
                let mut temp = 0.0;
                for l in a..b {
                    temp += u[l * n + i] * wrk[l - a];
                }
                out[i] += temp;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Form<'a> = &'a dyn Fn(usize, &[usize], &mut [f64]);

    /// The BLAS path (a block of at least `SMALL` points) and the plain loops
    /// compute the same form, up to rounding.
    #[test]
    fn both_paths_agree() {
        let n = 3;
        let dim = 2;
        let big = 12;
        let u: Vec<f64> = (0..big * n)
            .map(|i| 0.1 + (i as f64 * 0.37).sin())
            .collect();
        let du: Vec<f64> = (0..dim * big * n)
            .map(|i| (i as f64 * 0.71).cos())
            .collect();
        let w: Vec<f64> = (0..big).map(|t| 0.5 + 0.01 * t as f64).collect();
        let al: Vec<f64> = (0..big).map(|t| 1.0 + 0.1 * t as f64).collect();
        let be: Vec<[f64; 3]> = (0..big).map(|t| [t as f64, 1.0 - t as f64, 0.0]).collect();

        // The same points split into blocks of 4 (each below SMALL) by
        // calling the form once per block on a sliced copy.
        let split = |f: Form<'_>, len: usize| {
            let mut g = vec![0.0; len];
            for c in [0..4, 4..8, 8..12] {
                f(c.len(), &c.collect::<Vec<_>>(), &mut g);
            }
            g
        };
        let pick = |v: &[f64], idx: &[usize], per: usize| -> Vec<f64> {
            idx.iter()
                .flat_map(|&t| v[t * per..(t + 1) * per].to_vec())
                .collect()
        };
        let pick_grad = |idx: &[usize]| -> Vec<f64> {
            (0..dim)
                .flat_map(|k| pick(&du[k * big * n..(k + 1) * big * n], idx, n))
                .collect()
        };

        let close = |a: &[f64], b: &[f64]| {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() <= 1e-13 * (1.0 + y.abs()), "{x} vs {y}");
            }
        };

        let mut g = vec![0.0; n * n];
        grad_u_dot_grad_u(big, n, dim, &du, &w, &mut g, &al);
        let s = split(
            &|m, idx, g| {
                grad_u_dot_grad_u(
                    m,
                    n,
                    dim,
                    &pick_grad(idx),
                    &pick(&w, idx, 1),
                    g,
                    &pick(&al, idx, 1),
                )
            },
            n * n,
        );
        close(&g, &s);

        let mut g = vec![0.0; n * n];
        grad_u_dot_u(big, n, dim, &du, &u, &w, &mut g, &al, &be);
        let s = split(
            &|m, idx, g| {
                let b: Vec<[f64; 3]> = idx.iter().map(|&t| be[t]).collect();
                grad_u_dot_u(
                    m,
                    n,
                    dim,
                    &pick_grad(idx),
                    &pick(&u, idx, n),
                    &pick(&w, idx, 1),
                    g,
                    &pick(&al, idx, 1),
                    &b,
                )
            },
            n * n,
        );
        close(&g, &s);

        let mut g = vec![0.0; n * n];
        u_dot_v(big, n, &u, &u, &w, &mut g, Some(&al));
        let s = split(
            &|m, idx, g| {
                let uu = pick(&u, idx, n);
                u_dot_v(
                    m,
                    n,
                    &uu,
                    &uu,
                    &pick(&w, idx, 1),
                    g,
                    Some(&pick(&al, idx, 1)),
                )
            },
            n * n,
        );
        close(&g, &s);

        let mut f = vec![0.0; n];
        u_dot_f(big, n, &u, &w, &al, &mut f);
        let s = split(
            &|m, idx, f| {
                u_dot_f(
                    m,
                    n,
                    &pick(&u, idx, n),
                    &pick(&w, idx, 1),
                    &pick(&al, idx, 1),
                    f,
                )
            },
            n,
        );
        close(&f, &s);
    }
}
