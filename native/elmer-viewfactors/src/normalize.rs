//! `NormalizeFactors` and the file mask of `WriteOutputFile`
//! (`fem/src/ViewFactors.F90`).

/// `EPSILON(1._dp)`.
const EPS: f64 = f64::EPSILON;

/// Normalises the `n * n` row-major `factors` of surfaces with the given
/// `areas` so that every surface sees exactly the half space, and leaves
/// them divided by the area of the viewing surface.
///
/// The areas times the factors are first forced symmetric (a non-finite
/// entry is replaced by its partner). With `open` surfaces that is all: the
/// factors need not sum to one. Otherwise the scaling `D F D = A` is solved
/// by Newton iteration (at most 20 steps, stopping at `sum(r^2/A)/n <=
/// EPSILON`), each step by conjugate gradients with the diagonal as the
/// preconditioner, started from the current scaling, to a residual of 1e-6
/// of the right-hand side in at most 100 iterations, as `IterSolv` does.
pub fn normalize_factors(factors: &mut [f64], areas: &[f64], open: bool) {
    let n = areas.len();
    assert_eq!(factors.len(), n * n);

    for i in 0..n {
        for j in i..n {
            let si = areas[i] * factors[i * n + j];
            let sj = areas[j] * factors[j * n + i];
            let li = si.abs() < f64::MAX;
            let lj = sj.abs() < f64::MAX;
            let s = match (li, lj) {
                (true, true) => (si + sj) / 2.0,
                (true, false) => si,
                (false, true) => sj,
                (false, false) => 0.0,
            };
            factors[i * n + j] = s;
            factors[j * n + i] = s;
        }
    }

    if open {
        for i in 0..n {
            for j in 0..n {
                factors[i * n + j] /= areas[i];
            }
        }
        return;
    }

    let mut sol = vec![1.0; n];
    let mut rhs = vec![0.0; n];
    let mut jac = vec![0.0; n * n];
    let mut jdiag = vec![0.0; n];
    for _ in 0..20 {
        for i in 0..n {
            let mut cum = 0.0;
            for j in 0..n {
                cum += factors[i * n + j] * sol[j];
            }
            cum *= sol[i];
            rhs[i] = areas[i] - cum;
        }
        let cum: f64 = rhs.iter().zip(areas).map(|(r, a)| r * r / a).sum::<f64>() / n as f64;
        if cum <= EPS {
            break;
        }

        for i in 0..n {
            for j in 0..n {
                jac[i * n + j] = factors[i * n + j] * sol[i];
            }
            for j in 0..n {
                jac[i * n + i] += factors[i * n + j] * sol[j];
            }
            jdiag[i] = 1.0 / jac[i * n + i];
        }

        let prev = sol.clone();
        cg(n, &jac, &jdiag, &mut sol, &rhs);
        for (s, p) in sol.iter_mut().zip(&prev) {
            *s += p;
        }
    }

    for i in 0..n {
        for j in 0..n {
            factors[i * n + j] = factors[i * n + j] * sol[i] * sol[j] / areas[i];
        }
    }
}

/// `HUTI_D_CG` with `fm_DiagPrec`: `x` holds the initial guess and receives
/// the solution of `a x = b`.
fn cg(n: usize, a: &[f64], dinv: &[f64], x: &mut [f64], b: &[f64]) {
    const TOL: f64 = 1.0e-6;
    const MAXIT: usize = 100;
    let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let matvec = |v: &[f64], out: &mut [f64]| {
        for i in 0..n {
            out[i] = a[i * n..(i + 1) * n]
                .iter()
                .zip(v)
                .map(|(a, v)| a * v)
                .sum();
        }
    };
    let bnorm = norm(b);
    let mut r = vec![0.0; n];
    matvec(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    if norm(&r) <= TOL * bnorm {
        return;
    }
    let mut z: Vec<f64> = r.iter().zip(dinv).map(|(r, d)| r * d).collect();
    let mut p = z.clone();
    let mut rho: f64 = r.iter().zip(&z).map(|(r, z)| r * z).sum();
    let mut q = vec![0.0; n];
    for _ in 0..MAXIT {
        matvec(&p, &mut q);
        let pq: f64 = p.iter().zip(&q).map(|(p, q)| p * q).sum();
        let alpha = rho / pq;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        if norm(&r) <= TOL * bnorm {
            break;
        }
        for i in 0..n {
            z[i] = r[i] * dinv[i];
        }
        let rho_new: f64 = r.iter().zip(&z).map(|(r, z)| r * z).sum();
        let beta = rho_new / rho;
        rho = rho_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
}

/// The rows `WriteOutputFile` saves: of the `n_rows` rows of the `n`-column
/// `factors`, the columns whose factor exceeds `min_factor / 10`
/// (`Minimum View Factor`, 1e-20 by default, loosened by a factor of ten
/// "as the errors cannot be renormalised any more"). With `single_precision`
/// the values are rounded to `f32`, as the binary single precision file
/// stores them.
pub fn saved_rows(
    factors: &[f64],
    n: usize,
    n_rows: usize,
    min_factor: f64,
    single_precision: bool,
) -> Vec<Vec<(usize, f64)>> {
    let limit = min_factor / 10.0;
    (0..n_rows)
        .map(|i| {
            (0..n)
                .filter_map(|j| {
                    let v = factors[i * n + j];
                    (v > limit).then_some((j, if single_precision { v as f32 as f64 } else { v }))
                })
                .collect()
        })
        .collect()
}
