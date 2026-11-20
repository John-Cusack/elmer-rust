//! Static condensation of the local bubble functions of a HeatSolveVec
//! element: `CondensateP` (`CondensatePR`, MatrixAssembly.F90) for a steady
//! run, and `CondensatePTransient` (the BDF(1) form `Default1stOrderTime`
//! takes with `Nb > 0`) for a transient one, in the Fortran's arithmetic
//! order. One scalar dof per basis function; the bubbles are the trailing
//! `nb` of the `n + nb` local functions.

// Index loops spell out the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::LocalSystem;

use crate::error::HeatError;

/// `LUDecomp` (LinearAlgebra.F90:56): in-place LU with column pivoting, the
/// diagonal replaced by its reciprocal. Column-major `a[i + j*n]`. A copy of
/// elmer-fe's private port (`lu.rs`).
fn lu_decomp(a: &mut [f64], n: usize, pivot: &mut [usize]) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[ix(i, k)].abs() > a[ix(i, j)].abs() {
                j = k;
            }
        }
        if a[ix(i, j)].abs() == 0.0 {
            return None;
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a.swap(ix(k, j), ix(k, i));
            }
        }
        for k in i + 1..n {
            a[ix(i, k)] /= a[ix(i, i)];
        }
        for k in i + 1..n {
            if j != i {
                a.swap(ix(k, i), ix(k, j));
            }
            for l in i + 1..n {
                a[ix(k, l)] -= a[ix(k, i)] * a[ix(i, l)];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[ix(i, i)].abs() == 0.0 {
            return None;
        }
        a[ix(i, i)] = 1.0 / a[ix(i, i)];
    }
    Some(())
}

/// `InvertMatrix` (LinearAlgebra.F90:191): `a` (column-major) by its inverse.
fn invert_matrix(a: &mut [f64], n: usize) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    let mut pivot = vec![0usize; n];
    lu_decomp(a, n, &mut pivot)?;
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = -a[ix(i, j)];
            for k in i + 1..j {
                s -= a[ix(i, k)] * a[ix(k, j)];
            }
            a[ix(i, j)] = s;
        }
    }
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s -= a[ix(j, k)] * a[ix(k, i)];
            }
            a[ix(j, i)] = a[ix(i, i)] * s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                if k != i {
                    s += a[ix(i, k)] * a[ix(k, j)];
                } else {
                    s += a[ix(k, j)];
                }
            }
            a[ix(i, j)] = s;
        }
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            for j in 0..n {
                a.swap(ix(i, j), ix(pivot[i], j));
            }
        }
    }
    Some(())
}

/// The blocks `CondensatePR` forms: `Kbb^-1` (column-major `nb x nb`) and
/// `Fb`.
struct Bubble {
    inv: Vec<f64>,
    fb: Vec<f64>,
}

fn bubble_block(local: &LocalSystem, n: usize, nb: usize) -> Result<Bubble, HeatError> {
    let m = n + nb;
    let mut inv = vec![0.0; nb * nb];
    for p in 0..nb {
        for q in 0..nb {
            inv[p + q * nb] = local.stiff[(n + p) * m + n + q];
        }
    }
    invert_matrix(&mut inv, nb)
        .ok_or_else(|| HeatError::Invalid("InvertMatrix: the bubble block is singular".into()))?;
    let fb = (0..nb).map(|p| local.force[n + p]).collect();
    Ok(Bubble { inv, fb })
}

/// `F(1:n) -= MATMUL(Klb, MATMUL(Kbb, Fb))` and `K(1:n,1:n) -= MATMUL(Klb,
/// MATMUL(Kbb, Kbl))`, then `local` keeps the `n` retained functions.
fn eliminate(local: &mut LocalSystem, n: usize, nb: usize, b: &Bubble) {
    let m = n + nb;
    let k = |i: usize, j: usize| local.stiff[i * m + j];
    // MATMUL(Kbb, Fb)
    let t: Vec<f64> = (0..nb)
        .map(|p| {
            let mut s = 0.0;
            for q in 0..nb {
                s += b.inv[p + q * nb] * b.fb[q];
            }
            s
        })
        .collect();
    // MATMUL(Kbb, Kbl), nb x n
    let mut w = vec![0.0; nb * n];
    for p in 0..nb {
        for c in 0..n {
            let mut s = 0.0;
            for q in 0..nb {
                s += b.inv[p + q * nb] * k(n + q, c);
            }
            w[p * n + c] = s;
        }
    }
    let mut stiff = vec![0.0; n * n];
    let mut force = vec![0.0; n];
    for r in 0..n {
        let mut s = 0.0;
        for p in 0..nb {
            s += k(r, n + p) * t[p];
        }
        force[r] = local.force[r] - s;
        for c in 0..n {
            let mut s = 0.0;
            for p in 0..nb {
                s += k(r, n + p) * w[p * n + c];
            }
            stiff[r * n + c] = k(r, c) - s;
        }
    }
    local.n = n;
    local.stiff = stiff;
    local.force = force;
    local.mass = None;
    local.damp = None;
}

/// `CondensateP(n, nb, STIFF, FORCE)`: the trailing `nb` functions of the
/// `n + nb` in `local` are eliminated; `local` keeps the first `n`.
pub(crate) fn condensate_p(local: &mut LocalSystem, n: usize, nb: usize) -> Result<(), HeatError> {
    if nb == 0 {
        return Ok(());
    }
    check(local, n, nb)?;
    let b = bubble_block(local, n, nb)?;
    eliminate(local, n, nb, &b);
    Ok(())
}

/// `CondensatePTransient(n, nb, 1, dt, M, K, F, xprev, x, BubblePrev,
/// BubbleCur)`: implicit Euler on the whole local system (`K += M/dt`,
/// `F += M xlprev/dt`, `xlprev` being `xprev` at the retained functions and
/// `bubble_prev` at the bubbles), then [`condensate_p`]. Returns the current
/// bubble values `Kbb^-1 (Fb - Kbl x)` for the history, `x` being the
/// current iterate at the retained functions.
pub(crate) fn condensate_p_transient(
    local: &mut LocalSystem,
    n: usize,
    nb: usize,
    dt: f64,
    xprev: &[f64],
    x: &[f64],
    bubble_prev: &[f64],
) -> Result<Vec<f64>, HeatError> {
    check(local, n, nb)?;
    if xprev.len() != n || x.len() != n || bubble_prev.len() != nb {
        return Err(HeatError::Invalid(
            "CondensatePTransient: inconsistent local sizes".into(),
        ));
    }
    let m = n + nb;
    let Some(mass) = local.mass.take() else {
        return Err(HeatError::Invalid(
            "CondensatePTransient: a transient element without a mass matrix".into(),
        ));
    };
    let xl: Vec<f64> = xprev.iter().chain(bubble_prev).copied().collect();
    for i in 0..m {
        for j in 0..m {
            local.stiff[i * m + j] += mass[i * m + j] / dt;
        }
    }
    for i in 0..m {
        let mut s = 0.0;
        for j in 0..m {
            s += mass[i * m + j] * xl[j];
        }
        local.force[i] += s / dt;
    }
    let b = bubble_block(local, n, nb)?;
    // MATMUL(Kbb, Fb - MATMUL(Kbl, x))
    let r: Vec<f64> = (0..nb)
        .map(|p| {
            let mut s = 0.0;
            for c in 0..n {
                s += local.stiff[(n + p) * m + c] * x[c];
            }
            b.fb[p] - s
        })
        .collect();
    let cur = (0..nb)
        .map(|p| {
            let mut s = 0.0;
            for q in 0..nb {
                s += b.inv[p + q * nb] * r[q];
            }
            s
        })
        .collect();
    eliminate(local, n, nb, &b);
    Ok(cur)
}

fn check(local: &LocalSystem, n: usize, nb: usize) -> Result<(), HeatError> {
    let m = n + nb;
    if local.n != m || local.stiff.len() != m * m || local.force.len() != m {
        return Err(HeatError::Invalid(format!(
            "CondensateP: a local system of {} functions for {n} + {nb}",
            local.n
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system() -> LocalSystem {
        // A symmetric positive definite 4x4, the last two functions bubbles.
        let stiff = vec![
            4.0, 1.0, 0.5, 0.2, //
            1.0, 3.0, 0.3, 0.1, //
            0.5, 0.3, 2.0, 0.4, //
            0.2, 0.1, 0.4, 5.0,
        ];
        LocalSystem {
            n: 4,
            stiff,
            force: vec![1.0, 2.0, 3.0, 4.0],
            mass: None,
            damp: None,
        }
    }

    /// The condensed system's solution is the full system's at the retained
    /// functions.
    #[test]
    fn condensation_is_the_schur_complement() {
        let full = system();
        let mut c = full.clone();
        condensate_p(&mut c, 2, 2).unwrap();
        // Solve the 2x2 condensed system.
        let (a, b, cc, d) = (c.stiff[0], c.stiff[1], c.stiff[2], c.stiff[3]);
        let det = a * d - b * cc;
        let x = [
            (d * c.force[0] - b * c.force[1]) / det,
            (a * c.force[1] - cc * c.force[0]) / det,
        ];
        // The full system's residual at the retained rows vanishes with the
        // bubbles recovered as Kbb^-1 (Fb - Kbl x).
        let k = |i: usize, j: usize| full.stiff[i * 4 + j];
        let rb = [
            full.force[2] - k(2, 0) * x[0] - k(2, 1) * x[1],
            full.force[3] - k(3, 0) * x[0] - k(3, 1) * x[1],
        ];
        let (p, q, r, s) = (k(2, 2), k(2, 3), k(3, 2), k(3, 3));
        let dd = p * s - q * r;
        let xb = [(s * rb[0] - q * rb[1]) / dd, (p * rb[1] - r * rb[0]) / dd];
        let xs = [x[0], x[1], xb[0], xb[1]];
        for i in 0..4 {
            let res: f64 = full.force[i] - (0..4).map(|j| k(i, j) * xs[j]).sum::<f64>();
            assert!(res.abs() < 1e-13, "row {i}: {res}");
        }
    }

    /// The transient form is the steady one of `K + M/dt`, `F + M xlprev/dt`.
    #[test]
    fn transient_is_implicit_euler_then_condensation() {
        let mut t = system();
        let mass: Vec<f64> = (0..16)
            .map(|i| if i % 5 == 0 { 0.5 } else { 0.1 })
            .collect();
        t.mass = Some(mass.clone());
        let (dt, xprev, bprev, x) = (0.25, [1.0, -1.0], [0.5, 0.25], [2.0, 0.5]);
        let cur = condensate_p_transient(&mut t, 2, 2, dt, &xprev, &x, &bprev).unwrap();

        let mut s = system();
        let xl = [xprev[0], xprev[1], bprev[0], bprev[1]];
        for i in 0..4 {
            for j in 0..4 {
                s.stiff[i * 4 + j] += mass[i * 4 + j] / dt;
                s.force[i] += mass[i * 4 + j] * xl[j] / dt;
            }
        }
        let k = s.clone();
        condensate_p(&mut s, 2, 2).unwrap();
        for (a, b) in t
            .stiff
            .iter()
            .zip(&s.stiff)
            .chain(t.force.iter().zip(&s.force))
        {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
        // The recovered bubbles solve the bubble rows for the given x.
        for p in 0..2 {
            let row = 2 + p;
            let res = k.force[row]
                - (0..2).map(|c| k.stiff[row * 4 + c] * x[c]).sum::<f64>()
                - (0..2)
                    .map(|q| k.stiff[row * 4 + 2 + q] * cur[q])
                    .sum::<f64>();
            assert!(res.abs() < 1e-12, "bubble row {p}: {res}");
        }
    }
}
