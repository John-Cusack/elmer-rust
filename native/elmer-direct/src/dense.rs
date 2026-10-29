//! Dense LAPACK replacements on faer (the `d/zgetrf`, `d/zgetrs`, `dpotrf`,
//! `dpotrs`, `dsyev`, `dsygv`, `dgeev`, `dgesvd` import surface of
//! `libelmersolver.so`, strategy report 2.2).
//!
//! Matrices are column-major contiguous slices, as LAPACK sees them from
//! Fortran (`lda == m`). Pivot indices are 0-based and in LAPACK's sequential
//! interchange form: row `i` was swapped with row `ipiv[i]`.

use elmer_core::LinalgError;
use faer::linalg::triangular_solve::{
    solve_lower_triangular_in_place, solve_unit_lower_triangular_in_place,
    solve_upper_triangular_in_place,
};
use faer::{Mat, MatMut, MatRef, Par, Side};
use num_complex::Complex64;

use crate::scalar::DirectScalar;

fn view<T: DirectScalar>(a: &[T], m: usize, n: usize) -> Result<MatRef<'_, T>, LinalgError> {
    if a.len() != m * n {
        return Err(LinalgError::Dim(format!("{m}x{n} matrix needs {} entries, got {}", m * n, a.len())));
    }
    Ok(MatRef::from_column_major_slice(a, m, n))
}

fn copy_out<T: DirectScalar>(src: MatRef<'_, T>, dst: &mut [T]) {
    let (m, n) = (src.nrows(), src.ncols());
    for j in 0..n {
        for i in 0..m {
            dst[j * m + i] = src[(i, j)];
        }
    }
}

fn to_vec<T: DirectScalar>(src: MatRef<'_, T>) -> Vec<T> {
    let mut v = vec![T::zero(); src.nrows() * src.ncols()];
    copy_out(src, &mut v);
    v
}

/// LU with partial pivoting of an `m x n` matrix (`DGETRF`/`ZGETRF`). On return
/// `a` holds `L` (unit diagonal, below) and `U` (on and above) and the result is
/// the interchange vector (length `min(m, n)`). A zero `U[j,j]` is reported as
/// `Err(Singular(j))` after the factors have been written, as LAPACK's `INFO > 0`.
pub fn getrf<T: DirectScalar>(m: usize, n: usize, a: &mut [T]) -> Result<Vec<usize>, LinalgError> {
    let k = m.min(n);
    if k == 0 {
        return Ok(Vec::new());
    }
    let lu = view(a, m, n)?.partial_piv_lu();
    let (l, u) = (lu.L(), lu.U());
    for j in 0..n {
        for i in 0..m {
            a[j * m + i] = if i > j { l[(i, j)] } else if i < k { u[(i, j)] } else { T::zero() };
        }
    }
    // faer returns the final row order; replay it as sequential interchanges.
    let (fwd, _) = lu.P().arrays();
    let mut pos: Vec<usize> = (0..m).collect(); // pos[r] = current slot of original row r
    let mut at: Vec<usize> = (0..m).collect(); // at[s] = original row currently in slot s
    let mut ipiv = Vec::with_capacity(k);
    for i in 0..k {
        let want = fwd[i];
        let s = pos[want];
        ipiv.push(s);
        let displaced = at[i];
        at.swap(i, s);
        pos[want] = i;
        pos[displaced] = s;
    }
    for j in 0..k {
        if a[j * m + j] == T::zero() {
            return Err(LinalgError::Singular(j));
        }
    }
    Ok(ipiv)
}

/// Solve `A X = B` from [`getrf`] factors of a square `A` (`DGETRS('N')`).
/// `b` is `n x nrhs`, column-major.
pub fn getrs<T: DirectScalar>(n: usize, nrhs: usize, lu: &[T], ipiv: &[usize], b: &mut [T]) -> Result<(), LinalgError> {
    let lu = view(lu, n, n)?;
    if ipiv.len() != n || b.len() != n * nrhs {
        return Err(LinalgError::Dim("getrs: ipiv/rhs size mismatch".into()));
    }
    for col in b.chunks_exact_mut(n).take(nrhs) {
        for (i, &p) in ipiv.iter().enumerate() {
            col.swap(i, p);
        }
    }
    let mut rhs = MatMut::from_column_major_slice_mut(b, n, nrhs);
    solve_unit_lower_triangular_in_place(lu, rhs.as_mut(), Par::Seq);
    solve_upper_triangular_in_place(lu, rhs.as_mut(), Par::Seq);
    Ok(())
}

/// Cholesky `A = L L^H` of the lower triangle (`DPOTRF('L')`/`ZPOTRF('L')`).
/// Overwrites the lower triangle of `a`; the strict upper triangle is untouched.
/// `Err(Singular(j))` when the leading minor of order `j+1` is not positive definite.
pub fn potrf<T: DirectScalar>(n: usize, a: &mut [T]) -> Result<(), LinalgError> {
    if n == 0 {
        return Ok(());
    }
    let llt = view(a, n, n)?.llt(Side::Lower).map_err(|e| match e {
        faer::linalg::solvers::LltError::NonPositivePivot { index } => LinalgError::Singular(index),
    })?;
    let l = llt.L();
    for j in 0..n {
        for i in j..n {
            a[j * n + i] = l[(i, j)];
        }
    }
    Ok(())
}

/// Solve `A X = B` from [`potrf`] factors (`DPOTRS('L')`).
pub fn potrs<T: DirectScalar>(n: usize, nrhs: usize, l: &[T], b: &mut [T]) -> Result<(), LinalgError> {
    let l = view(l, n, n)?;
    if b.len() != n * nrhs {
        return Err(LinalgError::Dim("potrs: rhs size mismatch".into()));
    }
    let mut rhs = MatMut::from_column_major_slice_mut(b, n, nrhs);
    solve_lower_triangular_in_place(l, rhs.as_mut(), Par::Seq);
    solve_upper_triangular_in_place(l.adjoint(), rhs.as_mut(), Par::Seq);
    Ok(())
}

/// Eigen-decomposition result: `values[k]` pairs with column `k` of `vectors`
/// (column-major `n x n`).
#[derive(Clone, Debug)]
pub struct Eigen<V> {
    pub values: Vec<V>,
    pub vectors: Vec<V>,
}

fn evd_err(what: &str, e: impl core::fmt::Debug) -> LinalgError {
    LinalgError::Breakdown { method: "dense eigen", detail: format!("{what}: {e:?}") }
}

/// Symmetric eigenproblem `A v = w v` from the lower triangle (`DSYEV('V','L')`).
/// Eigenvalues ascending, eigenvectors orthonormal.
pub fn syev(n: usize, a: &[f64]) -> Result<Eigen<f64>, LinalgError> {
    if n == 0 {
        return Ok(Eigen { values: vec![], vectors: vec![] });
    }
    let evd = view(a, n, n)?.self_adjoint_eigen(Side::Lower).map_err(|e| evd_err("syev", e))?;
    Ok(Eigen {
        values: (0..n).map(|i| evd.S().column_vector()[i]).collect(),
        vectors: to_vec(evd.U()),
    })
}

/// Generalized symmetric-definite eigenproblem `A x = w B x`, `B` SPD, lower
/// triangles read (`DSYGV(1,'V','L')`). Eigenvalues ascending; eigenvectors are
/// `B`-orthonormal (`X^T B X = I`). `Err(Singular(j))` if `B` is not positive definite.
pub fn sygv(n: usize, a: &[f64], b: &[f64]) -> Result<Eigen<f64>, LinalgError> {
    if n == 0 {
        return Ok(Eigen { values: vec![], vectors: vec![] });
    }
    let (a, b) = (view(a, n, n)?, view(b, n, n)?);
    let l = b
        .llt(Side::Lower)
        .map_err(|e| match e {
            faer::linalg::solvers::LltError::NonPositivePivot { index } => LinalgError::Singular(index),
        })?;
    let l = l.L();
    // Symmetrise the lower triangle of A, then C = L^-1 A L^-T.
    let mut c = Mat::<f64>::from_fn(n, n, |i, j| if i >= j { a[(i, j)] } else { a[(j, i)] });
    solve_lower_triangular_in_place(l, c.as_mut(), Par::Seq); // L^-1 A
    let mut ct = c.transpose().to_owned(); // A L^-T transposed = (L^-1 A)^T
    solve_lower_triangular_in_place(l, ct.as_mut(), Par::Seq); // L^-1 (L^-1 A)^T = C^T
    let c = ct; // C symmetric
    let evd = c.as_ref().self_adjoint_eigen(Side::Lower).map_err(|e| evd_err("sygv", e))?;
    let mut x = evd.U().to_owned();
    solve_upper_triangular_in_place(l.transpose(), x.as_mut(), Par::Seq); // x = L^-T y
    Ok(Eigen {
        values: (0..n).map(|i| evd.S().column_vector()[i]).collect(),
        vectors: to_vec(x.as_ref()),
    })
}

/// General real eigenproblem `A v = w v` (`DGEEV('N','V')`), returning complex
/// eigenvalues and right eigenvectors. Each vector has unit Euclidean norm and
/// its largest component real, the `DGEEV` normalisation. Conjugate pairs are
/// adjacent, positive imaginary part first.
pub fn geev(n: usize, a: &[f64]) -> Result<Eigen<Complex64>, LinalgError> {
    if n == 0 {
        return Ok(Eigen { values: vec![], vectors: vec![] });
    }
    let evd = view(a, n, n)?.eigen().map_err(|e| evd_err("geev", e))?;
    let mut order: Vec<usize> = (0..n).collect();
    let w: Vec<Complex64> = (0..n).map(|i| evd.S().column_vector()[i]).collect();
    // Pair conjugates: stable sort keeping input order, then +imag before -imag.
    order.sort_by(|&p, &q| {
        let (a, b) = (w[p], w[q]);
        if a.re == b.re && a.im.abs() == b.im.abs() && a.im != b.im {
            b.im.partial_cmp(&a.im).unwrap_or(core::cmp::Ordering::Equal)
        } else {
            p.cmp(&q)
        }
    });
    let u = evd.U();
    let mut values = Vec::with_capacity(n);
    let mut vectors = Vec::with_capacity(n * n);
    for &k in &order {
        values.push(w[k]);
        let col: Vec<Complex64> = (0..n).map(|i| u[(i, k)]).collect();
        let norm = col.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        let big = col.iter().copied().fold(Complex64::new(0.0, 0.0), |m, z| if z.norm_sqr() > m.norm_sqr() { z } else { m });
        let phase = if big.norm() > 0.0 { big.conj() / big.norm() } else { Complex64::new(1.0, 0.0) };
        vectors.extend(col.iter().map(|&z| z * phase / norm));
    }
    Ok(Eigen { values, vectors })
}

/// Singular value decomposition `A = U diag(s) V^T` (`DGESVD('A','A')`).
#[derive(Clone, Debug)]
pub struct Svd {
    /// Singular values, non-increasing, length `min(m, n)`.
    pub s: Vec<f64>,
    /// `m x m` left singular vectors, column-major.
    pub u: Vec<f64>,
    /// `n x n` transposed right singular vectors `V^T`, column-major.
    pub vt: Vec<f64>,
}

pub fn gesvd(m: usize, n: usize, a: &[f64]) -> Result<Svd, LinalgError> {
    if m == 0 || n == 0 {
        return Ok(Svd { s: vec![], u: vec![], vt: vec![] });
    }
    let svd = view(a, m, n)?
        .svd()
        .map_err(|e| LinalgError::Breakdown { method: "gesvd", detail: format!("{e:?}") })?;
    let k = m.min(n);
    Ok(Svd {
        s: (0..k).map(|i| svd.S().column_vector()[i]).collect(),
        u: to_vec(svd.U()),
        vt: to_vec(svd.V().transpose()),
    })
}

/// Solve a dense square system `A x = b` with LU (convenience over [`getrf`] + [`getrs`]).
pub fn gesv<T: DirectScalar>(n: usize, nrhs: usize, a: &mut [T], b: &mut [T]) -> Result<(), LinalgError> {
    let ipiv = getrf(n, n, a)?;
    getrs(n, nrhs, a, &ipiv, b)
}
