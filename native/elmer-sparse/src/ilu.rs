//! Incomplete factorizations: `CRS_IncompleteLU`, `CRS_ComplexIncompleteLU`,
//! `CRS_ILUT`, `CRS_ComplexILUT`, the incomplete Cholesky branch and the
//! diagonal preconditioner, with the `CRS_LUSolve` triangular solves.
//!
//! Loop bodies are transliterations of CRSMatrix.F90 (lines 4065-5570) with
//! the same operation order, so factors and solves match the serial Fortran
//! build bit for bit. Factors are stored as Elmer stores them: strictly lower
//! multipliers, the upper part, and the *reciprocal* of the pivot on the
//! diagonal slot.

use crate::scalar::{diag_index, norm2, ElmerScalar, AEPS};
use elmer_core::{CsrMatrix, LinalgError, Preconditioner, Scalar};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    /// ILU(n) or ILUT factors.
    Lu,
    /// Incomplete Cholesky (`A % Cholesky`): lower triangle only.
    Cholesky,
    /// `CRS_DiagPrecondition`: `vals` holds the matrix diagonal.
    Jacobi,
}

/// An incomplete-factorization (or diagonal) preconditioner.
#[derive(Clone, Debug)]
pub struct Ilu<T: Scalar = f64> {
    n: usize,
    kind: Kind,
    rows: Vec<usize>,
    cols: Vec<u32>,
    diag: Vec<usize>,
    vals: Vec<T>,
}

/// Structure of the matrix in a form the factorization loops use.
struct Pattern<'a> {
    n: usize,
    rows: &'a [usize],
    cols: &'a [u32],
    diag: Vec<usize>,
}

fn pattern<'a, T: Scalar>(
    a: &'a CsrMatrix<T>,
    need_diag: bool,
) -> Result<Pattern<'a>, LinalgError> {
    if a.nrows != a.ncols {
        return Err(LinalgError::Dim(format!(
            "incomplete factorization needs a square matrix, got {}x{}",
            a.nrows, a.ncols
        )));
    }
    if a.row_ptr.len() != a.nrows + 1 || a.col_idx.len() != a.values.len() {
        return Err(LinalgError::Dim("inconsistent CSR arrays".into()));
    }
    let mut diag = Vec::with_capacity(a.nrows);
    for i in 0..a.nrows {
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        if lo == hi {
            return Err(LinalgError::Singular(i));
        }
        if !a.col_idx[lo..hi].windows(2).all(|w| w[0] < w[1]) {
            return Err(LinalgError::Dim(format!(
                "row {i} columns are not strictly ascending (call sort_rows)"
            )));
        }
        match diag_index(a, i) {
            Some(d) => diag.push(d),
            None if need_diag => return Err(LinalgError::Singular(i)),
            None => diag.push(usize::MAX),
        }
    }
    Ok(Pattern {
        n: a.nrows,
        rows: &a.row_ptr,
        cols: &a.col_idx,
        diag,
    })
}

/// One round of `InitializeILU1`: add the fill from eliminating the strictly
/// lower entries of each row against the rows above.
fn ilu1_pattern(
    n: usize,
    rows: &[usize],
    cols: &[u32],
    diag: &[usize],
) -> (Vec<usize>, Vec<u32>, Vec<usize>) {
    let mut c = vec![0u8; n];
    let mut out_rows = vec![0usize; n + 1];
    let mut out_cols: Vec<u32> = Vec::with_capacity(cols.len() * 2);
    let mut out_diag = vec![0usize; n];
    for i in 0..n {
        for k in rows[i]..rows[i + 1] {
            c[cols[k] as usize] = 1;
        }
        let row_min = cols[rows[i]] as usize;
        let mut row_max = cols[rows[i + 1] - 1] as usize;
        for k in row_min..i {
            if c[k] == 1 {
                for l in diag[k] + 1..rows[k + 1] {
                    let j = cols[l] as usize;
                    if c[j] == 0 {
                        c[j] = 2;
                        row_max = row_max.max(j);
                    }
                }
            }
        }
        for k in row_min..=row_max {
            if c[k] > 0 {
                c[k] = 0;
                if k == i {
                    out_diag[i] = out_cols.len();
                }
                out_cols.push(k as u32);
            }
        }
        out_rows[i + 1] = out_cols.len();
    }
    (out_rows, out_cols, out_diag)
}

/// Pattern of ILU(`level`): the matrix pattern for 0, else `level` rounds of fill.
fn ilu_n_pattern(p: &Pattern<'_>, level: usize) -> (Vec<usize>, Vec<u32>, Vec<usize>) {
    if level == 0 {
        return (p.rows.to_vec(), p.cols.to_vec(), p.diag.clone());
    }
    let mut cur = ilu1_pattern(p.n, p.rows, p.cols, &p.diag);
    for _ in 1..level {
        cur = ilu1_pattern(p.n, &cur.0, &cur.1, &cur.2);
    }
    cur
}

impl<T: ElmerScalar> Ilu<T> {
    /// ILU(`level`) of `a` (`CRS_IncompleteLU` / `CRS_ComplexIncompleteLU`).
    /// `a` must be square with sorted rows that all store their diagonal.
    pub fn ilu_n(a: &CsrMatrix<T>, level: usize) -> Result<Self, LinalgError> {
        Self::ilu_n_with_values(a, &a.values, level)
    }

    /// As [`Ilu::ilu_n`], factorizing `values` (which share `a`'s structure)
    /// instead of `a.values` (Elmer's `PrecValues`).
    pub fn ilu_n_with_values(
        a: &CsrMatrix<T>,
        values: &[T],
        level: usize,
    ) -> Result<Self, LinalgError> {
        if values.len() != a.values.len() {
            return Err(LinalgError::Dim("values length differs from nnz".into()));
        }
        let p = pattern(a, true)?;
        let n = p.n;
        let (lr, lc, ld) = ilu_n_pattern(&p, level);
        let mut lv = vec![T::zero(); lc.len()];
        let mut s = vec![T::zero(); n];
        let mut c = vec![false; n];
        for i in 0..n {
            for k in p.rows[i]..p.rows[i + 1] {
                s[p.cols[k] as usize] = values[k];
            }
            for k in lr[i]..lr[i + 1] {
                c[lc[k] as usize] = true;
            }
            for m in lr[i]..ld[i] {
                let k = lc[m] as usize;
                if s[k] == T::zero() {
                    continue;
                }
                let dk = lv[ld[k]];
                if dk.abs() > AEPS {
                    s[k] = s[k].elmer_div(dk);
                }
                for l in ld[k] + 1..lr[k + 1] {
                    let j = lc[l] as usize;
                    if c[j] {
                        s[j] = s[j] - s[k] * lv[l];
                    }
                }
            }
            for k in lr[i]..lr[i + 1] {
                let col = lc[k] as usize;
                lv[k] = s[col];
                s[col] = T::zero();
                c[col] = false;
            }
        }
        invert_pivots(&ld, &mut lv);
        Ok(Ilu {
            n,
            kind: Kind::Lu,
            rows: lr,
            cols: lc,
            diag: ld,
            vals: lv,
        })
    }

    /// ILUT with drop tolerance `tol` relative to the row norm of `a`
    /// (`CRS_ILUT` / `CRS_ComplexILUT`).
    pub fn ilut(a: &CsrMatrix<T>, tol: f64) -> Result<Self, LinalgError> {
        Self::ilut_with_values(a, &a.values, tol)
    }

    /// As [`Ilu::ilut`], factorizing `values` instead of `a.values`.
    pub fn ilut_with_values(a: &CsrMatrix<T>, values: &[T], tol: f64) -> Result<Self, LinalgError> {
        if values.len() != a.values.len() {
            return Err(LinalgError::Dim("values length differs from nnz".into()));
        }
        let p = pattern(a, false)?;
        let n = p.n;
        // The Fortran norm of a complex row runs over the interleaved (re, -im) reals.
        let row_norms: Vec<f64> = (0..n)
            .map(|i| {
                let mut comps = Vec::with_capacity(T::COMPONENTS * (p.rows[i + 1] - p.rows[i]));
                for k in p.rows[i]..p.rows[i + 1] {
                    values[k].for_each_component(&mut |x| comps.push(x));
                }
                norm2(comps.into_iter())
            })
            .collect();

        let mut lr = vec![0usize; n + 1];
        let mut lc: Vec<u32> = Vec::new();
        let mut lv: Vec<T> = Vec::new();
        let mut ld = vec![usize::MAX; n];
        let mut s = vec![T::zero(); n];
        let mut c = vec![false; n];
        for i in 0..n {
            for k in p.rows[i]..p.rows[i + 1] {
                let col = p.cols[k] as usize;
                c[col] = true;
                s[col] = values[k];
            }
            let row_min = p.cols[p.rows[i]] as usize;
            let mut row_max = p.cols[p.rows[i + 1] - 1] as usize;
            for k in row_min..i {
                if c[k] {
                    let dk = lv[ld[k]];
                    if dk.abs() > AEPS {
                        s[k] = s[k].elmer_div(dk);
                    }
                    for l in ld[k] + 1..lr[k + 1] {
                        let j = lc[l] as usize;
                        if !c[j] {
                            c[j] = true;
                            row_max = row_max.max(j);
                        }
                        s[j] = s[j] - s[k] * lv[l];
                    }
                }
            }
            let thresh = tol * row_norms[i];
            for k in row_min..=row_max {
                if c[k] {
                    let m = s[k].abs();
                    let keep = if T::DROP_INCLUSIVE {
                        m >= thresh
                    } else {
                        m > thresh
                    };
                    if keep || k == i {
                        if k == i {
                            ld[i] = lc.len();
                        }
                        lc.push(k as u32);
                        lv.push(s[k]);
                    }
                    s[k] = T::zero();
                    c[k] = false;
                }
            }
            if ld[i] == usize::MAX {
                return Err(LinalgError::Singular(i));
            }
            lr[i + 1] = lc.len();
        }
        invert_pivots(&ld, &mut lv);
        Ok(Ilu {
            n,
            kind: Kind::Lu,
            rows: lr,
            cols: lc,
            diag: ld,
            vals: lv,
        })
    }

    /// `CRS_DiagPrecondition` / `CRS_ComplexDiagPrecondition`: `z = r / diag(A)`.
    /// A tiny real pivot (`<= AEPS`) leaves that entry of `r` unchanged.
    pub fn diagonal(a: &CsrMatrix<T>) -> Self {
        let vals = (0..a.nrows)
            .map(|i| diag_index(a, i).map_or(T::zero(), |j| a.values[j]))
            .collect();
        Ilu {
            n: a.nrows,
            kind: Kind::Jacobi,
            rows: Vec::new(),
            cols: Vec::new(),
            diag: Vec::new(),
            vals,
        }
    }

    /// Order of the system.
    pub fn n(&self) -> usize {
        self.n
    }
    /// Stored factor entries (zero for the diagonal preconditioner).
    pub fn nnz(&self) -> usize {
        self.cols.len()
    }
    /// Row pointers of the factor pattern (`ILURows`, 0-based).
    pub fn row_ptr(&self) -> &[usize] {
        &self.rows
    }
    /// Column indices of the factor pattern (`ILUCols`, 0-based).
    pub fn col_idx(&self) -> &[u32] {
        &self.cols
    }
    /// Diagonal slots of the factor (`ILUDiag`, 0-based).
    pub fn diag_idx(&self) -> &[usize] {
        &self.diag
    }
    /// Factor values (`ILUValues`): reciprocal pivots on the diagonal slots.
    pub fn values(&self) -> &[T] {
        &self.vals
    }

    /// `b <- M^-1 b` (`CRS_LUSolve` / `CRS_ComplexLUSolve`).
    pub fn solve_in_place(&self, b: &mut [T]) {
        assert!(b.len() >= self.n, "rhs shorter than the system");
        match self.kind {
            Kind::Jacobi => {
                for (bi, &d) in b.iter_mut().zip(&self.vals) {
                    if !T::JACOBI_GUARD || d.abs() > AEPS {
                        *bi = bi.elmer_div(d);
                    }
                }
            }
            Kind::Lu => lu_solve(self, b),
            Kind::Cholesky => cholesky_solve(self, b),
        }
    }

    /// `b <- M^-T b` (plain transpose, no conjugation).
    pub fn solve_transpose_in_place(&self, b: &mut [T]) {
        assert!(b.len() >= self.n, "rhs shorter than the system");
        match self.kind {
            Kind::Jacobi | Kind::Cholesky => self.solve_in_place(b),
            Kind::Lu => {
                for i in 0..self.n {
                    let wi = b[i] * self.vals[self.diag[i]];
                    b[i] = wi;
                    for j in self.diag[i] + 1..self.rows[i + 1] {
                        let c = self.cols[j] as usize;
                        b[c] = b[c] - self.vals[j] * wi;
                    }
                }
                for i in (0..self.n).rev() {
                    let zi = b[i];
                    for j in self.rows[i]..self.diag[i] {
                        let c = self.cols[j] as usize;
                        b[c] = b[c] - self.vals[j] * zi;
                    }
                }
            }
        }
    }
}

/// "Prescale the diagonal for the LU solve": pivots below `AEPS` become 1.
fn invert_pivots<T: ElmerScalar>(ld: &[usize], lv: &mut [T]) {
    for &d in ld {
        lv[d] = if lv[d].abs() < AEPS {
            T::one()
        } else {
            T::one().elmer_div(lv[d])
        };
    }
}

fn lu_solve<T: ElmerScalar>(f: &Ilu<T>, b: &mut [T]) {
    let (rows, cols, diag, vals) = (&f.rows, &f.cols, &f.diag, &f.vals);
    for i in 0..f.n {
        let (mut s1, mut s2, mut s3, mut s4) = (b[i], T::zero(), T::zero(), T::zero());
        let k = diag[i] - rows[i];
        let end4 = rows[i] + 4 * (k / 4);
        let mut j = rows[i];
        while j < end4 {
            s1 = s1 - vals[j] * b[cols[j] as usize];
            s2 = s2 - vals[j + 1] * b[cols[j + 1] as usize];
            s3 = s3 - vals[j + 2] * b[cols[j + 2] as usize];
            s4 = s4 - vals[j + 3] * b[cols[j + 3] as usize];
            j += 4;
        }
        for j in end4..diag[i] {
            s1 = s1 - vals[j] * b[cols[j] as usize];
        }
        b[i] = s1 + s2 + s3 + s4;
    }
    for i in (0..f.n).rev() {
        let (mut s1, mut s2, mut s3, mut s4) = (b[i], T::zero(), T::zero(), T::zero());
        let k = rows[i + 1] - diag[i] - 1;
        let start = diag[i] + 1;
        let end4 = start + 4 * (k / 4);
        let mut j = start;
        while j < end4 {
            s1 = s1 - vals[j] * b[cols[j] as usize];
            s2 = s2 - vals[j + 1] * b[cols[j + 1] as usize];
            s3 = s3 - vals[j + 2] * b[cols[j + 2] as usize];
            s4 = s4 - vals[j + 3] * b[cols[j + 3] as usize];
            j += 4;
        }
        for j in end4..rows[i + 1] {
            s1 = s1 - vals[j] * b[cols[j] as usize];
        }
        b[i] = vals[diag[i]] * (s1 + s2 + s3 + s4);
    }
}

fn cholesky_solve<T: ElmerScalar>(f: &Ilu<T>, b: &mut [T]) {
    let (rows, cols, diag, vals) = (&f.rows, &f.cols, &f.diag, &f.vals);
    for i in 0..f.n {
        let (mut s1, mut s2, mut s3, mut s4) = (b[i], T::zero(), T::zero(), T::zero());
        let k = diag[i] - rows[i];
        let end4 = rows[i] + 4 * (k / 4);
        let mut j = rows[i];
        while j < end4 {
            s1 = s1 - vals[j] * b[cols[j] as usize];
            s2 = s2 - vals[j + 1] * b[cols[j + 1] as usize];
            s3 = s3 - vals[j + 2] * b[cols[j + 2] as usize];
            s4 = s4 - vals[j + 3] * b[cols[j + 3] as usize];
            j += 4;
        }
        for j in end4..diag[i] {
            s1 = s1 - vals[j] * b[cols[j] as usize];
        }
        b[i] = (s1 + s2 + s3 + s4) * vals[diag[i]];
    }
    for i in (0..f.n).rev() {
        b[i] = b[i] * vals[diag[i]];
        for j in rows[i]..diag[i] {
            let c = cols[j] as usize;
            b[c] = b[c] - vals[j] * b[i];
        }
    }
}

impl Ilu<f64> {
    /// Incomplete Cholesky on the ILU(`level`) pattern (`CRS_IncompleteLU` with
    /// `A % Cholesky`), for symmetric positive definite matrices. Only the lower
    /// triangle of `a` is read. `retry` is `Linear System Cholesky Retry` with
    /// its `Linear System Cholesky Factor`: on a nonpositive pivot the diagonal
    /// is scaled up by `1 + factor` and the factorization restarts.
    pub fn ic_n(a: &CsrMatrix<f64>, level: usize, retry: Option<f64>) -> Result<Self, LinalgError> {
        let p = pattern(a, true)?;
        let n = p.n;
        let (lr, lc, ld) = ilu_n_pattern(&p, level);
        let mut lv = vec![0.0; lc.len()];
        let mut s = vec![0.0; n];
        let mut t = vec![0.0; n];
        let c_factor = retry.map(|f| f + 1.0);
        let mut scl = 1.0_f64;
        'attempt: loop {
            s.fill(0.0);
            t.fill(0.0);
            for i in 0..n {
                for k in p.rows[i]..=p.diag[i] {
                    t[p.cols[k] as usize] = a.values[k];
                }
                for k in lr[i]..=ld[i] {
                    s[lc[k] as usize] = lv[k];
                }
                s[i] = scl * t[i];
                for m in lr[i]..ld[i] {
                    let j = lc[m] as usize;
                    let (mut r1, mut r2, mut r3, mut r4) = (t[j], 0.0, 0.0, 0.0);
                    let pn = ld[j] - lr[j];
                    let end4 = lr[j] + 4 * (pn / 4);
                    let mut l = lr[j];
                    while l < end4 {
                        r1 -= s[lc[l] as usize] * lv[l];
                        r2 -= s[lc[l + 1] as usize] * lv[l + 1];
                        r3 -= s[lc[l + 2] as usize] * lv[l + 2];
                        r4 -= s[lc[l + 3] as usize] * lv[l + 3];
                        l += 4;
                    }
                    for l in end4..ld[j] {
                        r1 -= s[lc[l] as usize] * lv[l];
                    }
                    s[j] = (r1 + r2 + r3 + r4) * lv[ld[j]];
                    s[i] -= s[j] * s[j];
                }
                if s[i] <= AEPS {
                    s[i] = 1.0;
                    if let Some(cf) = c_factor {
                        scl *= cf;
                        continue 'attempt;
                    }
                } else {
                    s[i] = 1.0 / s[i].sqrt();
                }
                for k in p.rows[i]..=p.diag[i] {
                    t[p.cols[k] as usize] = 0.0;
                }
                for k in lr[i]..=ld[i] {
                    let j = lc[k] as usize;
                    lv[k] = s[j];
                    s[j] = 0.0;
                }
            }
            break;
        }
        Ok(Ilu {
            n,
            kind: Kind::Cholesky,
            rows: lr,
            cols: lc,
            diag: ld,
            vals: lv,
        })
    }
}

impl<T: ElmerScalar> Preconditioner<T> for Ilu<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z[..self.n].copy_from_slice(&r[..self.n]);
        self.solve_in_place(z);
    }

    fn apply_transpose(&self, r: &[T], z: &mut [T]) -> Result<(), LinalgError> {
        z[..self.n].copy_from_slice(&r[..self.n]);
        self.solve_transpose_in_place(z);
        Ok(())
    }
}
