//! CRS operations from `CRSMatrix.F90`, in 0-based form.
//!
//! The sums keep the Fortran association order (four running accumulators for
//! the general `spmv`), so results agree bit for bit with the serial,
//! non-MKL, non-OpenMP build rather than only up to rounding.

use elmer_core::{CsrMatrix, LinalgError, LinearOperator, Scalar};

/// Marker in `CsrMatrix::diag` for a row that stores no diagonal entry
/// (Fortran's `Diag(i) = 0`).
pub const NO_DIAG: usize = usize::MAX;

/// Skip marker for `glue_local` indices (Fortran's `Indeces(i) <= 0`).
pub const NO_INDEX: u32 = u32::MAX;

fn check_spmv<T: Scalar>(a: &CsrMatrix<T>, values: &[T], x: &[T], y: &[T]) {
    assert_eq!(values.len(), a.col_idx.len(), "values/col_idx length");
    assert!(x.len() >= a.ncols, "x shorter than ncols");
    assert!(y.len() >= a.nrows, "y shorter than nrows");
}

/// `y = A x` in the order of `CRS_MatrixVectorMultiply`'s general path
/// (`ndeg = 1`): groups of four entries into four accumulators, the tail into
/// the first, combined as `((r1 + r2) + r3) + r4`.
pub fn spmv<T: Scalar>(a: &CsrMatrix<T>, x: &[T], y: &mut [T]) {
    spmv_with_values(a, &a.values, x, y);
}

/// [`spmv`] against another coefficient array sharing `a`'s structure
/// (the `UseValues` argument: mass, damping, preconditioner values).
pub fn spmv_with_values<T: Scalar>(a: &CsrMatrix<T>, values: &[T], x: &[T], y: &mut [T]) {
    check_spmv(a, values, x, y);
    for i in 0..a.nrows {
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        let (mut r1, mut r2, mut r3, mut r4) = (T::zero(), T::zero(), T::zero(), T::zero());
        let k = hi - lo;
        let end4 = lo + 4 * (k / 4);
        let mut j = lo;
        while j < end4 {
            r1 = r1 + x[a.col_idx[j] as usize] * values[j];
            r2 = r2 + x[a.col_idx[j + 1] as usize] * values[j + 1];
            r3 = r3 + x[a.col_idx[j + 2] as usize] * values[j + 2];
            r4 = r4 + x[a.col_idx[j + 3] as usize] * values[j + 3];
            j += 4;
        }
        for j in end4..hi {
            r1 = r1 + x[a.col_idx[j] as usize] * values[j];
        }
        y[i] = r1 + r2 + r3 + r4;
    }
}

/// `y = A x` for matrices stored as dense `ndeg`-blocks per row, in the order
/// of the `SELECT CASE( A % ndeg )` branches (`ndeg` in 2, 3, 4, 5, 6, 8, 10).
/// Each row must hold a multiple of `ndeg` entries and the columns of each
/// group of `w` entries must be consecutive. `ndeg = 1` is plain [`spmv`].
pub fn spmv_blocked<T: Scalar>(
    a: &CsrMatrix<T>,
    ndeg: usize,
    x: &[T],
    y: &mut [T],
) -> Result<(), LinalgError> {
    check_spmv(a, &a.values, x, y);
    let width = match ndeg {
        5 | 10 => 5,
        4 | 8 => 4,
        3 | 6 => 3,
        2 => 2,
        _ => {
            spmv(a, x, y);
            return Ok(());
        }
    };
    for i in 0..a.nrows {
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        if (hi - lo) % ndeg != 0 {
            return Err(LinalgError::Dim(format!(
                "row {i} holds {} entries, which is not a multiple of ndeg = {ndeg}",
                hi - lo
            )));
        }
        let mut r = [T::zero(); 5];
        let mut j = lo;
        while j < hi {
            let l = a.col_idx[j] as usize;
            for (t, rt) in r.iter_mut().enumerate().take(width) {
                *rt = *rt + x[l + t] * a.values[j + t];
            }
            j += width;
        }
        let mut s = r[0];
        for rt in &r[1..width] {
            s = s + *rt;
        }
        y[i] = s;
    }
    Ok(())
}

/// `y = A^T x` (`CRS_TransposeMatrixVectorMultiply`, no conjugation).
pub fn spmv_transpose<T: Scalar>(a: &CsrMatrix<T>, x: &[T], y: &mut [T]) {
    assert!(x.len() >= a.nrows && y.len() >= a.ncols);
    y[..a.ncols].fill(T::zero());
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            let c = a.col_idx[j] as usize;
            y[c] = y[c] + a.values[j] * x[i];
        }
    }
}

/// [`LinearOperator`] over a borrowed matrix using the Elmer-ordered [`spmv`]
/// (core's own impl on `CsrMatrix` is a plain sequential sum).
#[derive(Clone, Copy, Debug)]
pub struct CsrOp<'a, T: Scalar = f64>(pub &'a CsrMatrix<T>);

impl<T: Scalar> LinearOperator<T> for CsrOp<'_, T> {
    fn nrows(&self) -> usize {
        self.0.nrows
    }
    fn ncols(&self) -> usize {
        self.0.ncols
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        spmv(self.0, x, y);
    }
    fn apply_transpose(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() < self.0.nrows || y.len() < self.0.ncols {
            return Err(LinalgError::Dim("transpose vector too short".into()));
        }
        spmv_transpose(self.0, x, y);
        Ok(())
    }
}

/// `CRS_SortMatrix`: sort the columns of every row ascending (values follow)
/// and recompute `diag` when present.
pub fn sort_rows<T: Scalar>(a: &mut CsrMatrix<T>) {
    for i in 0..a.nrows {
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        if a.col_idx[lo..hi].windows(2).all(|w| w[0] <= w[1]) {
            continue;
        }
        let mut pairs: Vec<(u32, T)> = (lo..hi).map(|k| (a.col_idx[k], a.values[k])).collect();
        pairs.sort_by_key(|p| p.0);
        for (k, (c, v)) in (lo..hi).zip(pairs) {
            a.col_idx[k] = c;
            a.values[k] = v;
        }
    }
    if a.diag.is_some() {
        a.diag = Some(find_diag(a));
    }
}

/// Diagonal positions of a sorted matrix; [`NO_DIAG`] where a row has none.
pub fn find_diag<T: Scalar>(a: &CsrMatrix<T>) -> Vec<usize> {
    (0..a.nrows)
        .map(|i| {
            let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
            a.col_idx[lo..hi]
                .iter()
                .position(|&c| c as usize == i)
                .map_or(NO_DIAG, |p| lo + p)
        })
        .collect()
}

/// Build a sorted CSR matrix from `(row, col, value)` triplets. Duplicates are
/// summed in input order, as repeated `List_AddToMatrixElement` calls would.
/// `diag` is filled in.
pub fn from_triplets<T: Scalar>(
    nrows: usize,
    ncols: usize,
    triplets: impl IntoIterator<Item = (usize, u32, T)>,
) -> CsrMatrix<T> {
    let mut t: Vec<(usize, u32, T)> = triplets.into_iter().collect();
    for &(r, c, _) in &t {
        assert!(r < nrows && (c as usize) < ncols, "triplet out of range");
    }
    t.sort_by_key(|e| (e.0, e.1)); // stable: duplicates keep input order
    let mut row_ptr = vec![0usize; nrows + 1];
    let mut col_idx: Vec<u32> = Vec::with_capacity(t.len());
    let mut values: Vec<T> = Vec::with_capacity(t.len());
    let mut prev: Option<(usize, u32)> = None;
    for (r, c, v) in t {
        if prev == Some((r, c)) {
            let last = values.len() - 1;
            values[last] = values[last] + v;
        } else {
            col_idx.push(c);
            values.push(v);
            row_ptr[r + 1] += 1;
            prev = Some((r, c));
        }
    }
    for i in 0..nrows {
        row_ptr[i + 1] += row_ptr[i];
    }
    let mut a = CsrMatrix {
        nrows,
        ncols,
        row_ptr,
        col_idx,
        values,
        diag: None,
    };
    a.diag = Some(find_diag(&a));
    a
}

/// `CRS_Search`: position of column `col` in row `i`, if stored.
pub fn find_entry<T: Scalar>(a: &CsrMatrix<T>, i: usize, col: u32) -> Option<usize> {
    let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
    a.col_idx[lo..hi].binary_search(&col).ok().map(|p| lo + p)
}

/// `CRS_CheckMatrixElement`.
pub fn contains<T: Scalar>(a: &CsrMatrix<T>, i: usize, j: u32) -> bool {
    find_entry(a, i, j).is_some()
}

/// `CRS_GetMatrixElement`: the stored value, zero when absent.
pub fn get_element<T: Scalar>(a: &CsrMatrix<T>, i: usize, j: u32) -> T {
    find_entry(a, i, j).map_or(T::zero(), |k| a.values[k])
}

/// `CRS_AddToMatrixElement`. Returns `false` when `(i, j)` is not in the
/// structure and `val` is nonzero (Fortran then converts the matrix to list
/// form; callers decide what to do). Adding zero to a missing entry is a no-op.
pub fn add_to_element<T: Scalar>(a: &mut CsrMatrix<T>, i: usize, j: u32, val: T) -> bool {
    match find_entry(a, i, j) {
        Some(k) => {
            a.values[k] = a.values[k] + val;
            true
        }
        None => val == T::zero(),
    }
}

/// `CRS_SetMatrixElement`. Returns `false` when `(i, j)` is not stored.
pub fn set_element<T: Scalar>(a: &mut CsrMatrix<T>, i: usize, j: u32, val: T) -> bool {
    match find_entry(a, i, j) {
        Some(k) => {
            a.values[k] = val;
            true
        }
        None => false,
    }
}

/// `CRS_ZeroMatrix`: zero every coefficient, keep the structure.
pub fn zero_matrix<T: Scalar>(a: &mut CsrMatrix<T>) {
    a.values.fill(T::zero());
}

/// `CRS_ZeroRow`: zero the coefficients of row `i`.
pub fn zero_row<T: Scalar>(a: &mut CsrMatrix<T>, i: usize) {
    let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
    a.values[lo..hi].fill(T::zero());
}

/// `CRS_RowSum`.
pub fn row_sum<T: Scalar>(a: &CsrMatrix<T>, i: usize) -> T {
    let mut s = T::zero();
    for k in a.row_ptr[i]..a.row_ptr[i + 1] {
        s = s + a.values[k];
    }
    s
}

/// `CRS_RowSumAbs`.
pub fn row_sum_abs<T: Scalar>(a: &CsrMatrix<T>, i: usize) -> f64 {
    let mut s = 0.0;
    for k in a.row_ptr[i]..a.row_ptr[i + 1] {
        s += a.values[k].abs();
    }
    s
}

/// `CRS_GlueLocalMatrix`: add the dense `(n*dofs) x (n*dofs)` row-major
/// `local` matrix into `a`. `idx[i]` is the global node of local node `i`
/// ([`NO_INDEX`] skips it); global dof `dofs*node + c` pairs with local dof
/// `dofs*i + c`. Rows without a stored diagonal are skipped and entries that are
/// not in the structure are dropped, as in Fortran.
pub fn glue_local(a: &mut CsrMatrix<f64>, dofs: usize, idx: &[u32], local: &[f64]) {
    let mut values = std::mem::take(&mut a.values);
    glue_local_into(a, &mut values, dofs, idx, local);
    a.values = values;
}

/// [`glue_local`] into a separate coefficient array sharing `a`'s structure
/// (the `GlobalValues` argument).
pub fn glue_local_into(
    a: &CsrMatrix<f64>,
    values: &mut [f64],
    dofs: usize,
    idx: &[u32],
    local: &[f64],
) {
    let n = idx.len();
    let nl = n * dofs;
    assert_eq!(local.len(), nl * nl, "local matrix must be (n*dofs)^2");
    assert_eq!(values.len(), a.col_idx.len(), "values/col_idx length");
    for i in 0..n {
        if idx[i] == NO_INDEX {
            continue;
        }
        // Fortran loops k = 0..dofs-1 over `Dofs*Indeces(i) - k`.
        for kk in (0..dofs).rev() {
            let row = dofs * idx[i] as usize + kk;
            let lrow = dofs * i + kk;
            let (lo, hi) = (a.row_ptr[row], a.row_ptr[row + 1]);
            if let Some(d) = &a.diag {
                if d[row] == NO_DIAG {
                    continue;
                }
            }
            for j in 0..n {
                if idx[j] == NO_INDEX {
                    continue;
                }
                for ll in (0..dofs).rev() {
                    let col = (dofs * idx[j] as usize + ll) as u32;
                    if let Ok(p) = a.col_idx[lo..hi].binary_search(&col) {
                        values[lo + p] += local[lrow * nl + dofs * j + ll];
                    }
                }
            }
        }
    }
}
