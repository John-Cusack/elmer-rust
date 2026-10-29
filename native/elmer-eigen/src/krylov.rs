//! The restarted Krylov iteration that stands in for ARPACK's `DSAUPD`/`DSEUPD`
//! (symmetric) and `DNAUPD`/`DNEUPD` (general).
//!
//! Both keep a Krylov decomposition `OP V_k = V_k S + v_k b^T` with a basis
//! that is orthonormal in the `B` inner product (`BMAT = 'G'`), extend it to
//! `ncv` vectors with Arnoldi steps (classical Gram-Schmidt with the DGKS
//! correction, as `dsaitr`/`dnaitr`), and restart by keeping the wanted Ritz
//! vectors:
//!
//! * symmetric: thick-restart Lanczos (Wu and Simon, 2000), the symmetric
//!   Krylov-Schur method; the kept block is `diag(theta)`.
//! * general: Krylov-Schur style thick restart on the real invariant subspace
//!   of the wanted Ritz values (real and imaginary parts of the Ritz vectors,
//!   re-orthonormalised; Morgan 1996, Stewart 2001), so a conjugate pair is
//!   never split.
//!
//! Convergence is ARPACK's (`dsconv`/`dnconv`): a wanted Ritz value `theta`
//! is converged when its Ritz estimate is at most `tol * max(eps^(2/3), |theta|)`.

use elmer_core::LinalgError;
use elmer_direct::dense::{geev, syev};
use num_complex::Complex64;

use crate::config::ArpackWhich;

/// The operator as ARPACK's reverse communication sees it.
pub trait EigenOperator {
    fn dim(&self) -> usize;
    /// `BMAT = 'G'`: the basis is orthonormal in the `B` inner product.
    fn generalized(&self) -> bool;
    /// `y = B x`; only called when [`generalized`](Self::generalized).
    fn apply_b(&self, x: &[f64], y: &mut [f64]);
    /// `y = OP x`. `bx` is `B x` when generalized (ARPACK's `IDO = 1` hands it
    /// over so shift-invert need not recompute it), else `x`.
    fn apply_op(&mut self, x: &[f64], bx: &[f64], y: &mut [f64]) -> Result<(), LinalgError>;
}

#[derive(Clone, Debug)]
pub struct KrylovParams {
    pub nev: usize,
    pub ncv: usize,
    pub which: ArpackWhich,
    /// `<= 0` means machine precision.
    pub tol: f64,
    /// ARPACK `IPARAM(3)`: restarts allowed beyond the first factorisation.
    pub max_iter: usize,
    pub symmetric: bool,
    /// Replace each Ritz vector `y` by `OP y / theta` (the mode 3 purification
    /// of `dseupd`/`dneupd`), which removes components in the null space of a
    /// singular `B`.
    pub purify: bool,
    /// Seed of the start vector.
    pub seed: u64,
}

/// Wanted Ritz pairs of `OP`, in `which` order (conjugate pairs adjacent,
/// positive imaginary part first).
#[derive(Clone, Debug)]
pub struct RitzPairs {
    pub theta: Vec<Complex64>,
    /// Unit `B`-norm before purification.
    pub vectors: Vec<Vec<Complex64>>,
    /// ARPACK `IPARAM(3)` on output: factorisation rounds used.
    pub iterations: usize,
    /// Number of `OP` applications.
    pub op_applications: usize,
    /// ARPACK `IPARAM(5)`: wanted Ritz values meeting the tolerance.
    pub converged: usize,
}

/// `dlamch('E')`: relative machine precision with rounding.
const EPS: f64 = f64::EPSILON / 2.0;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

/// A basis vector and `B` times it.
type Direction = (Vec<f64>, Vec<f64>);

/// SplitMix64, uniform on (-1, 1) like `dlarnv(2, ...)` in `dgetv0`.
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }
}

/// Ordering key of `which`: entries that sort first are wanted. Ties are
/// broken so that a conjugate pair stays adjacent, positive imaginary first.
fn which_order(which: ArpackWhich, theta: &[Complex64]) -> Vec<usize> {
    let key = |z: &Complex64| -> f64 {
        match which {
            ArpackWhich::LM => -z.norm(),
            ArpackWhich::SM => z.norm(),
            ArpackWhich::LR => -z.re,
            ArpackWhich::SR => z.re,
            ArpackWhich::LI => -z.im.abs(),
            ArpackWhich::SI => z.im.abs(),
        }
    };
    let mut order: Vec<usize> = (0..theta.len()).collect();
    order.sort_by(|&p, &q| {
        let (a, b) = (&theta[p], &theta[q]);
        key(a)
            .total_cmp(&key(b))
            .then(a.re.total_cmp(&b.re))
            .then(b.im.total_cmp(&a.im))
    });
    order
}

struct Iteration<'a> {
    op: &'a mut dyn EigenOperator,
    n: usize,
    generalized: bool,
    ops: usize,
    rng: Rng,
}

impl Iteration<'_> {
    fn b(&self, x: &[f64]) -> Vec<f64> {
        if self.generalized {
            let mut y = vec![0.0; self.n];
            self.op.apply_b(x, &mut y);
            y
        } else {
            x.to_vec()
        }
    }

    fn op(&mut self, x: &[f64], bx: &[f64]) -> Result<Vec<f64>, LinalgError> {
        let mut y = vec![0.0; self.n];
        self.op.apply_op(x, bx, &mut y)?;
        self.ops += 1;
        Ok(y)
    }

    /// `dgetv0`: a random vector, forced into the range of `OP` for the
    /// generalized problem (B possibly singular), `B`-orthogonal to `basis`
    /// and of unit `B`-norm. `None` when every attempt is (numerically) in
    /// the span of `basis`.
    fn start_vector(&mut self, basis: &[Vec<f64>]) -> Result<Option<Direction>, LinalgError> {
        for _ in 0..5 {
            let mut x: Vec<f64> = (0..self.n).map(|_| self.rng.uniform()).collect();
            if self.generalized {
                let bx = self.b(&x);
                x = self.op(&x, &bx)?;
            }
            let mut bx = self.b(&x);
            let norm0 = dot(&x, &bx).max(0.0).sqrt();
            if norm0 == 0.0 {
                continue;
            }
            let mut norm = norm0;
            for _ in 0..2 {
                if basis.is_empty() {
                    break;
                }
                let c: Vec<f64> = basis.iter().map(|v| dot(v, &bx)).collect();
                for (v, ci) in basis.iter().zip(&c) {
                    axpy(-ci, v, &mut x);
                }
                bx = self.b(&x);
                norm = dot(&x, &bx).max(0.0).sqrt();
                if norm > 0.717 * norm0 {
                    break;
                }
            }
            if norm > 1e3 * EPS * norm0 {
                let s = 1.0 / norm;
                x.iter_mut().for_each(|t| *t *= s);
                bx.iter_mut().for_each(|t| *t *= s);
                return Ok(Some((x, bx)));
            }
        }
        Ok(None)
    }
}

/// Dense `(m+1) x m` projected matrix, column-major.
struct Hbar {
    m: usize,
    a: Vec<f64>,
}

impl Hbar {
    fn new(m: usize) -> Self {
        Hbar {
            m,
            a: vec![0.0; (m + 1) * m],
        }
    }
    fn get(&self, i: usize, j: usize) -> f64 {
        self.a[i + j * (self.m + 1)]
    }
    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i + j * (self.m + 1)] = v;
    }
}

/// Run the iteration. `p.ncv` vectors of length `op.dim()` are stored.
pub fn krylov_schur(
    op: &mut dyn EigenOperator,
    p: &KrylovParams,
) -> Result<RitzPairs, LinalgError> {
    let n = op.dim();
    let (nev, m) = (p.nev, p.ncv);
    // DSAUPD/DNAUPD argument checks (INFO = -1 .. -3).
    if nev == 0 || n == 0 {
        return Err(LinalgError::Dim("eigen: NEV and N must be positive".into()));
    }
    let gap = if p.symmetric { 1 } else { 2 };
    if m < nev + gap || m > n {
        return Err(LinalgError::Dim(format!(
            "eigen: need NEV + {gap} <= NCV <= N, got NEV = {nev}, NCV = {m}, N = {n} \
             (Eigen System Lanczos Vectors)"
        )));
    }
    if p.symmetric && !matches!(p.which, ArpackWhich::LM | ArpackWhich::SM) {
        // DSAUPD INFO = -5: Elmer maps only LM/SM for the symmetric solver.
        return Err(LinalgError::Unsupported(format!(
            "eigen: WHICH = {:?} with a symmetric matrix (DSAUPD info = -5)",
            p.which
        )));
    }
    let tol = if p.tol > 0.0 { p.tol } else { EPS };
    let eps23 = EPS.powf(2.0 / 3.0);

    let generalized = op.generalized();
    let mut it = Iteration {
        op,
        n,
        generalized,
        ops: 0,
        rng: Rng(p.seed),
    };

    let mut v: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
    let mut bv: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
    let (x0, bx0) = it
        .start_vector(&[])?
        .ok_or_else(|| LinalgError::Breakdown {
            method: "eigen",
            detail: "zero start vector".into(),
        })?;
    v.push(x0);
    bv.push(bx0);

    let mut h = Hbar::new(m);
    let mut kept = 0usize;
    let mut iter = 1usize;
    loop {
        // Extend the factorisation from `kept` to `m` columns.
        for j in kept..m {
            let mut w = it.op(&v[j], &bv[j])?;
            let mut bw = it.b(&w);
            let wnorm = dot(&w, &bw).max(0.0).sqrt();
            let mut hj: Vec<f64> = v[..=j].iter().map(|vi| dot(vi, &bw)).collect();
            for (vi, c) in v[..=j].iter().zip(&hj) {
                axpy(-c, vi, &mut w);
            }
            bw = it.b(&w);
            let mut beta = dot(&w, &bw).max(0.0).sqrt();
            let mut prev = wnorm;
            let mut rounds = 0;
            while beta <= 0.717 * prev && rounds < 2 {
                let c: Vec<f64> = v[..=j].iter().map(|vi| dot(vi, &bw)).collect();
                for ((vi, ci), hi) in v[..=j].iter().zip(&c).zip(hj.iter_mut()) {
                    axpy(-ci, vi, &mut w);
                    *hi += ci;
                }
                bw = it.b(&w);
                prev = beta;
                beta = dot(&w, &bw).max(0.0).sqrt();
                rounds += 1;
            }
            for (i, hi) in hj.iter().enumerate() {
                h.set(i, j, *hi);
            }
            let invariant = beta <= 0.717 * prev || beta <= EPS * wnorm;
            if invariant && j + 1 == m {
                // Exact last step: the residual is zero and no new direction is
                // needed (the restart supplies one if it must).
                h.set(j + 1, j, 0.0);
                v.truncate(j + 1);
                bv.truncate(j + 1);
                v.push(vec![0.0; n]);
                bv.push(vec![0.0; n]);
            } else if invariant {
                // Invariant subspace found: continue from a fresh direction.
                h.set(j + 1, j, 0.0);
                let fresh = it
                    .start_vector(&v[..=j])?
                    .ok_or_else(|| LinalgError::Breakdown {
                        method: "eigen",
                        detail: format!("Krylov space exhausted at dimension {}", j + 1),
                    })?;
                v.truncate(j + 1);
                bv.truncate(j + 1);
                v.push(fresh.0);
                bv.push(fresh.1);
            } else {
                h.set(j + 1, j, beta);
                let s = 1.0 / beta;
                w.iter_mut().for_each(|t| *t *= s);
                bw.iter_mut().for_each(|t| *t *= s);
                v.truncate(j + 1);
                bv.truncate(j + 1);
                v.push(w);
                bv.push(bw);
            }
        }

        // Ritz pairs of the m x m block.
        let beta_m = h.get(m, m - 1);
        let (theta, y) = ritz(&h, p.symmetric)?;
        let order = which_order(p.which, &theta);
        let est = |k: usize| beta_m.abs() * y[k][m - 1].norm();
        let nconv = order[..nev]
            .iter()
            .filter(|&&k| est(k) <= tol * eps23.max(theta[k].norm()))
            .count();

        if nconv >= nev {
            let mut out_theta = Vec::with_capacity(nev);
            let mut out_vec = Vec::with_capacity(nev);
            for &k in &order[..nev] {
                let mut x = vec![Complex64::new(0.0, 0.0); n];
                for (j, vj) in v[..m].iter().enumerate() {
                    let c = y[k][j];
                    for (xi, &vij) in x.iter_mut().zip(vj) {
                        *xi += c * vij;
                    }
                }
                if p.purify && theta[k].norm() > 0.0 {
                    let c = Complex64::new(beta_m, 0.0) * y[k][m - 1] / theta[k];
                    for (xi, &r) in x.iter_mut().zip(&v[m]) {
                        *xi += c * r;
                    }
                }
                out_theta.push(theta[k]);
                out_vec.push(x);
            }
            return Ok(RitzPairs {
                theta: out_theta,
                vectors: out_vec,
                iterations: iter,
                op_applications: it.ops,
                converged: nconv,
            });
        }
        if iter > p.max_iter {
            return Err(LinalgError::Breakdown {
                method: "eigen",
                detail: format!(
                    "Maximum number of iterations reached ({} of {nev} Ritz values converged)",
                    nconv
                ),
            });
        }

        // Restart size: dsaup2/dnaup2 keep NEV + min(NCONV, NP/2).
        let mut keep = nev + nconv.min((m - nev) / 2);
        if keep == 1 && m >= 6 {
            keep = m / 2;
        } else if keep == 1 && m > 2 {
            keep = 2;
        }
        keep = keep.min(m - gap);
        if !p.symmetric && theta[order[keep - 1]].im > 0.0 {
            // Do not split a conjugate pair: its partner is next in `order`.
            if keep + 1 < m {
                keep += 1;
            } else {
                keep -= 1;
            }
        }
        let q = restart_basis(&order[..keep], &theta, &y, m, p.symmetric);
        let keep = q.len();

        // V_new = V_m Q, v_keep = v_m; BV likewise.
        let comb = |src: &[Vec<f64>]| -> Vec<Vec<f64>> {
            q.iter()
                .map(|col| {
                    let mut out = vec![0.0; n];
                    for (s, &c) in src.iter().zip(col) {
                        axpy(c, s, &mut out);
                    }
                    out
                })
                .collect()
        };
        let mut nv = comb(&v[..m]);
        let mut nbv = comb(&bv[..m]);
        if beta_m == 0.0 {
            // Zero residual: b = 0 and the next direction is a fresh one.
            let (x, bx) = it
                .start_vector(&nv)?
                .ok_or_else(|| LinalgError::Breakdown {
                    method: "eigen",
                    detail: format!("Krylov space exhausted at dimension {keep}"),
                })?;
            nv.push(x);
            nbv.push(bx);
        } else {
            nv.push(v[m].clone());
            nbv.push(bv[m].clone());
        }

        // S = Q^T H_m Q, b = beta_m * Q[m-1, :].
        let mut nh = Hbar::new(m);
        for (b, qb) in q.iter().enumerate() {
            if p.symmetric {
                // Q holds eigenvectors of the symmetric block: S = diag(theta).
                nh.set(b, b, theta[order[b]].re);
            } else {
                let hq: Vec<f64> = (0..m)
                    .map(|i| (0..m).map(|l| h.get(i, l) * qb[l]).sum())
                    .collect();
                for (a, qa) in q.iter().enumerate() {
                    nh.set(a, b, dot(qa, &hq));
                }
            }
            nh.set(keep, b, beta_m * qb[m - 1]);
        }
        v = nv;
        bv = nbv;
        h = nh;
        kept = keep;
        iter += 1;
    }
}

/// Eigenpairs of the leading `m x m` block. Vectors have unit 2-norm.
fn ritz(h: &Hbar, symmetric: bool) -> Result<(Vec<Complex64>, Vec<Vec<Complex64>>), LinalgError> {
    let m = h.m;
    let mut a = vec![0.0; m * m];
    for j in 0..m {
        for i in 0..m {
            a[i + j * m] = if symmetric && i < j {
                h.get(j, i)
            } else {
                h.get(i, j)
            };
        }
    }
    if symmetric {
        let e = syev(m, &a)?;
        let theta = e.values.iter().map(|&t| Complex64::new(t, 0.0)).collect();
        let y = (0..m)
            .map(|k| {
                e.vectors[k * m..(k + 1) * m]
                    .iter()
                    .map(|&t| Complex64::new(t, 0.0))
                    .collect()
            })
            .collect();
        Ok((theta, y))
    } else {
        let e = geev(m, &a)?;
        let y = (0..m)
            .map(|k| e.vectors[k * m..(k + 1) * m].to_vec())
            .collect();
        Ok((e.values, y))
    }
}

/// Orthonormal real basis (columns of length `m`) of the span of the kept
/// Ritz vectors. A conjugate pair contributes its real and imaginary parts.
fn restart_basis(
    kept: &[usize],
    theta: &[Complex64],
    y: &[Vec<Complex64>],
    m: usize,
    symmetric: bool,
) -> Vec<Vec<f64>> {
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(kept.len());
    for &k in kept {
        if symmetric || theta[k].im == 0.0 {
            cols.push(y[k].iter().map(|z| z.re).collect());
        } else if theta[k].im > 0.0 {
            cols.push(y[k].iter().map(|z| z.re).collect());
            cols.push(y[k].iter().map(|z| z.im).collect());
        }
        // The negative-imaginary partner is covered by its pair.
    }
    if symmetric {
        return cols;
    }
    // Modified Gram-Schmidt, twice; drop (numerically) dependent columns.
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(cols.len());
    for mut c in cols {
        let n0 = dot(&c, &c).sqrt();
        for _ in 0..2 {
            for qi in &q {
                let d = dot(qi, &c);
                axpy(-d, qi, &mut c);
            }
        }
        let nrm = dot(&c, &c).sqrt();
        if nrm > 1e-10 * n0 {
            c.iter_mut().for_each(|t| *t /= nrm);
            q.push(c);
        }
    }
    debug_assert!(q.iter().all(|c| c.len() == m));
    q
}
