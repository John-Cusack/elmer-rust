//! Band storage and band solvers (`BandMatrix.F90`, `SolveBand.F90`,
//! `SolveSBand.F90`).
//!
//! The storage is the LAPACK band layout Elmer uses. For a general band matrix
//! with `kl` sub- and `ku` super-diagonals, column `j` occupies
//! `ldab = 2*kl + ku + 1` slots (`kl` spare rows for LU fill-in) and
//! `A[i,j]` lives at `ab[j*ldab + kl + ku + i - j]` (0-based). Elmer's
//! `BAND_INDEX(i,j)` is this with `kl = ku = Subband`. The symmetric variant
//! keeps the lower triangle only, `SBAND_INDEX(i,j) = j*(kd+1) + i - j`.
//!
//! [`BandLu`] is `DGBTRF`/`DGBTRS` (`ZGBTRF`/`ZGBTRS` for complex): the
//! blocked factorisation of LAPACK 3.12 when its block size applies (`NB =
//! 32` from `ILAENV` for `KU > 64`, used when `NB <= KL`), with its reference
//! BLAS kernels, else `DGBTF2`; and
//! [`BandCholesky`] is `DPBTF2`/`DPBTRS` (`ZPBTF2`/`ZPBTRS`), written in the
//! reference operation order (reciprocal pivot scaling, rank-1 updates column by
//! column) so results match the netlib reference routines bit for bit.

// `Scalar` has no `*Assign` bounds, so the reference-order `x = x + y` updates stay;
// index loops mirror the Fortran/LAPACK loops they transliterate.
#![allow(clippy::assign_op_pattern, clippy::needless_range_loop)]

use elmer_core::{CsrMatrix, DirectSolver, LinalgError, Scalar};

use crate::scalar::DirectScalar;

/// General band matrix in LAPACK `GB` storage with room for LU fill-in.
#[derive(Clone, Debug)]
pub struct BandMatrix<T: DirectScalar = f64> {
    n: usize,
    kl: usize,
    ku: usize,
    ab: Vec<T>,
}

impl<T: DirectScalar> BandMatrix<T> {
    /// Zero `n x n` band matrix with `kl` sub- and `ku` super-diagonals.
    pub fn new(n: usize, kl: usize, ku: usize) -> Self {
        let ldab = 2 * kl + ku + 1;
        Self { n, kl, ku, ab: vec![T::zero(); ldab * n] }
    }

    /// Elmer's `Band_CreateMatrix` shape: `kl = ku = subband`.
    pub fn new_subband(n: usize, subband: usize) -> Self {
        Self::new(n, subband, subband)
    }

    /// Band matrix holding the entries of `a`; `kl`/`ku` are the actual lower
    /// and upper bandwidths of its pattern. Duplicate entries are summed.
    pub fn from_csr(a: &CsrMatrix<T>) -> Result<Self, LinalgError> {
        if a.nrows != a.ncols {
            return Err(LinalgError::Dim(format!("band solve needs a square matrix, got {}x{}", a.nrows, a.ncols)));
        }
        check_csr(a)?;
        let (mut kl, mut ku) = (0, 0);
        for i in 0..a.nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k] as usize;
                if j >= a.ncols {
                    return Err(LinalgError::Dim(format!("column index {j} out of range")));
                }
                kl = kl.max(i.saturating_sub(j));
                ku = ku.max(j.saturating_sub(i));
            }
        }
        let mut m = Self::new(a.nrows, kl, ku);
        for i in 0..a.nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                m.add(i, a.col_idx[k] as usize, a.values[k]);
            }
        }
        Ok(m)
    }

    pub fn n(&self) -> usize {
        self.n
    }
    /// [`Self::from_csr`] in the shape Elmer's `BandSolver` gives `DGBTRF`:
    /// `kl = ku = Subband`, the largest `|i - j|` of the pattern.
    pub fn from_csr_subband(a: &CsrMatrix<T>) -> Result<Self, LinalgError> {
        let b = Self::from_csr(a)?;
        if b.kl == b.ku {
            return Ok(b);
        }
        let s = b.kl.max(b.ku);
        let mut m = Self::new(a.nrows, s, s);
        for i in 0..a.nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                m.add(i, a.col_idx[k] as usize, a.values[k]);
            }
        }
        Ok(m)
    }

    pub fn kl(&self) -> usize {
        self.kl
    }
    pub fn ku(&self) -> usize {
        self.ku
    }
    /// Leading dimension of the band array.
    pub fn ldab(&self) -> usize {
        2 * self.kl + self.ku + 1
    }
    /// Raw band array (column-major, `ldab` rows).
    pub fn data(&self) -> &[T] {
        &self.ab
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.ldab() + self.kl + self.ku + i - j
    }

    #[inline]
    fn in_band(&self, i: usize, j: usize) -> bool {
        i < self.n && j < self.n && i <= j + self.kl && j <= i + self.ku
    }

    /// `A[i,j]` (zero outside the band).
    pub fn get(&self, i: usize, j: usize) -> T {
        if self.in_band(i, j) {
            self.ab[self.idx(i, j)]
        } else {
            T::zero()
        }
    }

    /// Raw entry, including the LU fill-in rows above the original upper band
    /// (zero outside the `kl + ku` upper bandwidth).
    pub fn factor_entry(&self, i: usize, j: usize) -> T {
        if i < self.n && j < self.n && i <= j + self.kl && j <= i + self.ku + self.kl {
            self.ab[self.idx(i, j)]
        } else {
            T::zero()
        }
    }

    /// `A[i,j] = v`. Panics outside the band.
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        assert!(self.in_band(i, j), "BandMatrix::set({i},{j}) outside the band");
        let k = self.idx(i, j);
        self.ab[k] = v;
    }

    /// `A[i,j] += v` (`Band_AddToMatrixElement`). Panics outside the band.
    pub fn add(&mut self, i: usize, j: usize, v: T) {
        assert!(self.in_band(i, j), "BandMatrix::add({i},{j}) outside the band");
        let k = self.idx(i, j);
        self.ab[k] = self.ab[k] + v;
    }

    /// `Band_ZeroMatrix`.
    pub fn zero(&mut self) {
        self.ab.fill(T::zero());
    }

    /// `Band_ZeroRow`: zero row `n` within the band.
    pub fn zero_row(&mut self, n: usize) {
        for j in n.saturating_sub(self.kl)..=(n + self.ku).min(self.n - 1) {
            let k = self.idx(n, j);
            self.ab[k] = T::zero();
        }
    }

    /// `Band_GlueLocalMatrix`: add a local `(nodes*dofs)^2` row-major matrix.
    /// `indices` are 0-based node numbers; row `dofs*node + (dofs-1-k)` pairs with
    /// local row `dofs*i + (dofs-1-k)` exactly as the 1-based Fortran loop does.
    pub fn glue_local(&mut self, dofs: usize, indices: &[usize], local: &[T]) {
        let nl = indices.len() * dofs;
        assert_eq!(local.len(), nl * nl, "local matrix must be (nodes*dofs)^2");
        for (i, &ni) in indices.iter().enumerate() {
            for k in 0..dofs {
                let row = dofs * ni + k;
                for (j, &nj) in indices.iter().enumerate() {
                    for l in 0..dofs {
                        let col = dofs * nj + l;
                        self.add(row, col, local[(dofs * i + k) * nl + dofs * j + l]);
                    }
                }
            }
        }
    }

    /// `v = A u` (`Band_MatrixVectorMultiply`).
    pub fn matvec(&self, u: &[T], v: &mut [T]) {
        for i in 0..self.n {
            let mut s = T::zero();
            for j in i.saturating_sub(self.kl)..=(i + self.ku).min(self.n - 1) {
                s = s + u[j] * self.ab[self.idx(i, j)];
            }
            v[i] = s;
        }
    }

    /// `Band_MatrixVectorProd` with the transposed operator: `v = A^T u`.
    pub fn matvec_transpose(&self, u: &[T], v: &mut [T]) {
        v[..self.n].fill(T::zero());
        for i in 0..self.n {
            let s = u[i];
            for j in i.saturating_sub(self.kl)..=(i + self.ku).min(self.n - 1) {
                v[j] = v[j] + s * self.ab[self.idx(i, j)];
            }
        }
    }

    /// `Band_DiagPrecondition`: `u = D^-1 v`, leaving entries with
    /// `|D| <= aeps` unscaled.
    pub fn diag_precondition(&self, v: &[T], u: &mut [T], aeps: f64) {
        for i in 0..self.n {
            let d = self.ab[self.idx(i, i)];
            u[i] = if Scalar::abs(d) > aeps { v[i] / d } else { v[i] };
        }
    }
}

fn check_csr<T: Scalar>(a: &CsrMatrix<T>) -> Result<(), LinalgError> {
    if a.row_ptr.len() != a.nrows + 1
        || a.row_ptr.last().copied() != Some(a.col_idx.len())
        || a.col_idx.len() != a.values.len()
    {
        return Err(LinalgError::Dim("inconsistent CSR arrays".into()));
    }
    Ok(())
}

/// Band LU with partial pivoting: `DGBTRF`/`DGBTRS` (`ZGBTRF`/`ZGBTRS`), the
/// factorisation behind `Linear System Direct Method = Banded`.
#[derive(Clone, Debug)]
pub struct BandLu<T: DirectScalar = f64> {
    lu: Option<BandMatrix<T>>,
    ipiv: Vec<usize>,
}

impl<T: DirectScalar> Default for BandLu<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: DirectScalar> BandLu<T> {
    pub fn new() -> Self {
        Self { lu: None, ipiv: Vec::new() }
    }

    /// Factor a band matrix in place (consumes it; the L and U factors replace
    /// the band entries, exactly as `DGBTRF` overwrites `AB`).
    pub fn factor_band(&mut self, mut m: BandMatrix<T>) -> Result<(), LinalgError> {
        self.lu = None;
        let ipiv = gbtrf(&mut m)?;
        self.ipiv = ipiv;
        self.lu = Some(m);
        Ok(())
    }

    /// Row interchanges of the last factorisation: row `j` was swapped with `ipiv()[j]`.
    pub fn ipiv(&self) -> &[usize] {
        &self.ipiv
    }

    /// The factored band array (`L` multipliers below, `U` on and above the diagonal).
    pub fn factors(&self) -> Option<&BandMatrix<T>> {
        self.lu.as_ref()
    }

    /// Solve for `nrhs` right-hand sides stored column after column in `b`.
    pub fn solve_many(&self, b: &mut [T], nrhs: usize) -> Result<(), LinalgError> {
        let lu = self
            .lu
            .as_ref()
            .ok_or_else(|| LinalgError::Unsupported("solve before factor".into()))?;
        if b.len() != lu.n * nrhs {
            return Err(LinalgError::Dim(format!("rhs length {} != n*nrhs = {}", b.len(), lu.n * nrhs)));
        }
        gbtrs(lu, &self.ipiv, b, nrhs);
        Ok(())
    }
}

impl<T: DirectScalar> DirectSolver<T> for BandLu<T> {
    fn factor(&mut self, a: &CsrMatrix<T>) -> Result<(), LinalgError> {
        // BandSolver (DirectSolve.F90): DGBTRF( N, N, Subband, Subband, ... ).
        self.factor_band(BandMatrix::from_csr_subband(a)?)
    }

    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError> {
        self.solve_many(b, 1)
    }
}

/// `DGBTRF(M = N = n, KL, KU, AB, LDAB, IPIV, INFO)` (LAPACK 3.12): the
/// blocked right-looking band LU, with `DGBTF2` below the block size.
/// `NB = ILAENV(1, 'DGBTRF', ...)` is 32 when `KU > 64` and 1 otherwise;
/// the blocked code runs when `1 < NB <= KL`. Fortran's aliased views of `AB`
/// (`LDAB - 1` strides along a row of the matrix) are kept as linear indices
/// into one array, in the reference kernels' read and write order: `IDAMAX`,
/// `DSWAP`, `DSCAL`, `DGER`, `DCOPY`, `DLASWP`, `DTRSM('L','L','N','U')` and
/// `DGEMM('N','N')`. Returns the pivot rows (0-based).
fn gbtrf<T: DirectScalar>(m: &mut BandMatrix<T>) -> Result<Vec<usize>, LinalgError> {
    const NBMAX: usize = 64;
    const LDWORK: usize = NBMAX + 1;
    let (n, kl, ku) = (m.n, m.kl, m.ku);
    if n == 0 {
        return Ok(Vec::new());
    }
    let nb = if ku <= 64 { 1 } else { 32 }.min(NBMAX);
    if nb <= 1 || nb > kl {
        return gbtf2(m);
    }
    let mm = n; // M
    let kv = ku + kl;
    let ldab = m.ldab();
    let ld = ldab - 1; // the stride along a row of the matrix
    // 1-based (row, column) of AB, and of the work arrays.
    let at = |r: usize, c: usize| (c - 1) * ldab + (r - 1);
    let wt = |r: usize, c: usize| (c - 1) * LDWORK + (r - 1);
    let ab = &mut m.ab;
    let mut work13 = vec![T::zero(); LDWORK * NBMAX];
    let mut work31 = vec![T::zero(); LDWORK * NBMAX];
    let mut ipiv = vec![0usize; n]; // 1-based, as the Fortran
    let mut info = 0usize;
    // Zero the superdiagonal elements of the work array WORK13 and the
    // subdiagonal elements of WORK31 (they are already zero here).
    // Gaussian elimination with partial pivoting: set fill-in elements in
    // columns KU+2 to KV to zero.
    for j in ku + 2..=kv.min(n) {
        for i in kv + 2 - j..=kl {
            ab[at(i, j)] = T::zero();
        }
    }
    // JU is the index of the last column affected by the current stage of
    // the factorization.
    let mut ju = 1usize;
    let mn = mm.min(n);
    let mut j = 1usize;
    while j <= mn {
        let jb = nb.min(mn - j + 1);
        // The active part of the matrix is partitioned
        //    A11   A12   A13
        //    A21   A22   A23
        //    A31   A32   A33
        // Here A11, A21 and A31 denote the current block of JB columns
        // which is about to be factorized. The number of rows in the
        // partitioning are JB, I2, I3 respectively, and the numbers of
        // columns are JB, J2, J3. The superdiagonal elements of A13 and the
        // subdiagonal elements of A31 lie outside the band.
        let i2 = (kl as isize - jb as isize).min(mm as isize - j as isize - jb as isize + 1);
        let i3 = (jb as isize).min(mm as isize - j as isize - kl as isize + 1);
        // J2 and J3 are computed after JU has been updated.
        // Factorize the current block of JB columns.
        for jj in j..j + jb {
            // Set fill-in elements in column JJ+KV to zero.
            if jj + kv <= n {
                for i in 1..=kl {
                    ab[at(i, jj + kv)] = T::zero();
                }
            }
            // Find pivot and test for singularity. KM is the number of
            // subdiagonal elements in the current column.
            let km = kl.min(mm - jj);
            // IDAMAX(KM+1, AB(KV+1,JJ), 1)
            let mut jp = 1usize;
            let mut dmax = ab[at(kv + 1, jj)].abs1();
            for t in 2..=km + 1 {
                let v = ab[at(kv + t, jj)].abs1();
                if v > dmax {
                    jp = t;
                    dmax = v;
                }
            }
            ipiv[jj - 1] = jp + jj - j;
            if ab[at(kv + jp, jj)] != T::zero() {
                ju = ju.max((jj + ku + jp - 1).min(n));
                // Apply interchange to columns J to J+JB-1.
                if jp != 1 {
                    if jp + jj - 1 < j + kl {
                        let (x0, y0) = (at(kv + 1 + jj - j, j), at(kv + jp + jj - j, j));
                        for t in 0..jb {
                            ab.swap(x0 + t * ld, y0 + t * ld);
                        }
                    } else {
                        // The interchange affects columns J to JJ-1 of A31
                        // which are stored in the work array WORK31.
                        let x0 = at(kv + 1 + jj - j, j);
                        let w0 = wt(jp + jj - j - kl, 1);
                        for t in 0..jj - j {
                            std::mem::swap(&mut ab[x0 + t * ld], &mut work31[w0 + t * LDWORK]);
                        }
                        let (x0, y0) = (at(kv + 1, jj), at(kv + jp, jj));
                        for t in 0..j + jb - jj {
                            ab.swap(x0 + t * ld, y0 + t * ld);
                        }
                    }
                }
                // Compute multipliers: DSCAL(KM, ONE/AB(KV+1,JJ), AB(KV+2,JJ), 1).
                let r = T::one() / ab[at(kv + 1, jj)];
                let x0 = at(kv + 2, jj);
                for t in 0..km {
                    ab[x0 + t] = r * ab[x0 + t];
                }
                // Update trailing submatrix within the band and within the
                // current block. JM is the index of the last column which
                // needs to be updated.
                let jm = ju.min(j + jb - 1);
                if jm > jj {
                    // DGER(KM, JM-JJ, -ONE, AB(KV+2,JJ), 1, AB(KV,JJ+1),
                    //      LDAB-1, AB(KV+1,JJ+1), LDAB-1)
                    let (x0, y0, a0) = (at(kv + 2, jj), at(kv, jj + 1), at(kv + 1, jj + 1));
                    for c in 0..jm - jj {
                        let y = ab[y0 + c * ld];
                        if y != T::zero() {
                            let temp = -y;
                            for i in 0..km {
                                let k = a0 + c * ld + i;
                                ab[k] = ab[k] + ab[x0 + i] * temp;
                            }
                        }
                    }
                }
            } else if info == 0 {
                // If pivot is zero, set INFO to the index of the pivot unless
                // a zero pivot has already been found.
                info = jj;
            }
            // Copy current column of A31 into the work array WORK31.
            let nw = (jj as isize - j as isize + 1).min(i3);
            if nw > 0 {
                let x0 = at(kv + kl + 1 + j - jj, jj);
                let w0 = wt(1, jj - j + 1);
                let nw = nw as usize;
                work31[w0..w0 + nw].copy_from_slice(&ab[x0..x0 + nw]);
            }
        }
        if j + jb <= n {
            // Apply the row interchanges to the other blocks.
            // J2 may be negative, as in the Fortran (DLASWP and the
            // updates then do nothing).
            let j2 = (ju as isize - j as isize + 1).min(kv as isize) - jb as isize;
            let j2u = j2.max(0) as usize;
            let j3 = (ju as isize - j as isize - kv as isize + 1).max(0) as usize;
            // Use DLASWP to apply the row interchanges to A12, A22, and A32.
            // DLASWP(J2, AB(KV+1-JB,J+JB), LDAB-1, 1, JB, IPIV(J), 1)
            {
                let a0 = at(kv + 1 - jb, j + jb);
                for i in 1..=jb {
                    let ip = ipiv[j + i - 2];
                    if ip != i {
                        for c in 0..j2u {
                            ab.swap(a0 + c * ld + (i - 1), a0 + c * ld + (ip - 1));
                        }
                    }
                }
            }
            // Adjust the pivot indices.
            for i in j..j + jb {
                ipiv[i - 1] += j - 1;
            }
            // Apply the row interchanges to A13, A23, and A33 columnwise.
            let k2 = (j as isize - 1 + jb as isize + j2) as usize;
            for i in 1..=j3 {
                let jj = k2 + i;
                for ii in j + i - 1..j + jb {
                    let ip = ipiv[ii - 1];
                    if ip != ii {
                        ab.swap(at(kv + 1 + ii - jj, jj), at(kv + 1 + ip - jj, jj));
                    }
                }
            }
            // Update the relevant part of the trailing submatrix.
            if j2 > 0 {
                // Update A12: DTRSM('Left','Lower','No transpose','Unit', JB,
                // J2, ONE, AB(KV+1,J), LDAB-1, AB(KV+1-JB,J+JB), LDAB-1)
                let (a0, b0) = (at(kv + 1, j), at(kv + 1 - jb, j + jb));
                trsm_llnu(ab, a0, ld, None, b0, ld, jb, j2u);
                if i2 > 0 {
                    // Update A22: DGEMM(I2, J2, JB, -ONE, AB(KV+1+JB,J), LDAB-1,
                    // AB(KV+1-JB,J+JB), LDAB-1, ONE, AB(KV+1,J+JB), LDAB-1)
                    let (a0, b0, c0) = (at(kv + 1 + jb, j), at(kv + 1 - jb, j + jb), at(kv + 1, j + jb));
                    gemm_nn_minus(ab, None, a0, ld, b0, ld, c0, ld, i2 as usize, j2u, jb);
                }
                if i3 > 0 {
                    // Update A32: DGEMM(I3, J2, JB, -ONE, WORK31, LDWORK,
                    // AB(KV+1-JB,J+JB), LDAB-1, ONE, AB(KV+KL+1-JB,J+JB), LDAB-1)
                    let (b0, c0) = (at(kv + 1 - jb, j + jb), at(kv + kl + 1 - jb, j + jb));
                    gemm_nn_minus(ab, Some(&work31), 0, LDWORK, b0, ld, c0, ld, i3 as usize, j2u, jb);
                }
            }
            if j3 > 0 {
                // Copy the lower triangle of A13 into the work array WORK13.
                for jj in 1..=j3 {
                    for ii in jj..=jb {
                        work13[wt(ii, jj)] = ab[at(ii - jj + 1, jj + j + kv - 1)];
                    }
                }
                // Update A13 in the work array: DTRSM(..., JB, J3, ONE,
                // AB(KV+1,J), LDAB-1, WORK13, LDWORK)
                let a0 = at(kv + 1, j);
                trsm_llnu(ab, a0, ld, Some(&mut work13), 0, LDWORK, jb, j3);
                if i2 > 0 {
                    // Update A23: DGEMM(I2, J3, JB, -ONE, AB(KV+1+JB,J), LDAB-1,
                    // WORK13, LDWORK, ONE, AB(1+JB,J+KV), LDAB-1)
                    let (a0, c0) = (at(kv + 1 + jb, j), at(1 + jb, j + kv));
                    gemm_nn_b(ab, a0, ld, &work13, LDWORK, c0, ld, i2 as usize, j3, jb, None);
                }
                if i3 > 0 {
                    // Update A33: DGEMM(I3, J3, JB, -ONE, WORK31, LDWORK, WORK13,
                    // LDWORK, ONE, AB(1+KL,J+KV), LDAB-1)
                    let c0 = at(1 + kl, j + kv);
                    gemm_nn_b(ab, 0, LDWORK, &work13, LDWORK, c0, ld, i3 as usize, j3, jb, Some(&work31));
                }
                // Copy the lower triangle of A13 back into place.
                for jj in 1..=j3 {
                    for ii in jj..=jb {
                        ab[at(ii - jj + 1, jj + j + kv - 1)] = work13[wt(ii, jj)];
                    }
                }
            }
        } else {
            // Adjust the pivot indices.
            for i in j..j + jb {
                ipiv[i - 1] += j - 1;
            }
        }
        // Partially undo the interchanges in the current block to restore
        // the upper triangular form of A31 and copy the upper triangle of
        // A31 back into place.
        for jj in (j..j + jb).rev() {
            let jp = ipiv[jj - 1] - jj + 1;
            if jp != 1 {
                // Apply interchange to columns J to JJ-1.
                if jp + jj - 1 < j + kl {
                    // The interchange does not affect A31.
                    let (x0, y0) = (at(kv + 1 + jj - j, j), at(kv + jp + jj - j, j));
                    for t in 0..jj - j {
                        ab.swap(x0 + t * ld, y0 + t * ld);
                    }
                } else {
                    // The interchange does affect A31.
                    let x0 = at(kv + 1 + jj - j, j);
                    let w0 = wt(jp + jj - j - kl, 1);
                    for t in 0..jj - j {
                        std::mem::swap(&mut ab[x0 + t * ld], &mut work31[w0 + t * LDWORK]);
                    }
                }
            }
            // Copy the current column of A31 back into place.
            let nw = i3.min(jj as isize - j as isize + 1);
            if nw > 0 {
                let w0 = wt(1, jj - j + 1);
                let x0 = at(kv + kl + 1 + j - jj, jj);
                let nw = nw as usize;
                ab[x0..x0 + nw].copy_from_slice(&work31[w0..w0 + nw]);
            }
        }
        j += nb;
    }
    let ipiv: Vec<usize> = ipiv.into_iter().map(|p| p - 1).collect();
    if info != 0 {
        return Err(LinalgError::Singular(info - 1));
    }
    Ok(ipiv)
}

/// Reference `DTRSM('Left', 'Lower', 'No transpose', 'Unit', M, N, ONE, A,
/// LDA, B, LDB)`: `B := inv(A) B`. `A` lives in `ab` at `a0`; `B` in `ab` at
/// `b0` (`work = None`, possibly overlapping, as Fortran's aliased call) or
/// in `work` at `b0`.
#[allow(clippy::too_many_arguments)]
fn trsm_llnu<T: DirectScalar>(ab: &mut [T], a0: usize, lda: usize, mut work: Option<&mut Vec<T>>, b0: usize, ldb: usize, m: usize, n: usize) {
    for jc in 0..n {
        for k in 0..m {
            let bk = match work.as_deref() {
                Some(w) => w[b0 + jc * ldb + k],
                None => ab[b0 + jc * ldb + k],
            };
            if bk != T::zero() {
                for i in k + 1..m {
                    let aik = ab[a0 + k * lda + i];
                    match work.as_deref_mut() {
                        Some(w) => {
                            let p = b0 + jc * ldb + i;
                            w[p] = w[p] - bk * aik;
                        }
                        None => {
                            let p = b0 + jc * ldb + i;
                            ab[p] = ab[p] - bk * aik;
                        }
                    }
                }
            }
        }
    }
}

/// Reference `DGEMM('N', 'N', M, N, K, -ONE, A, LDA, B, LDB, ONE, C, LDC)`
/// with `B` and `C` in `ab`, `A` in `ab` (`aw = None`) or in `aw`:
/// `C(:,j) += TEMP*A(:,l)` with `TEMP = -B(l,j)`, for `l` then `i`.
#[allow(clippy::too_many_arguments)]
fn gemm_nn_minus<T: DirectScalar>(
    ab: &mut [T],
    aw: Option<&[T]>,
    a0: usize,
    lda: usize,
    b0: usize,
    ldb: usize,
    c0: usize,
    ldc: usize,
    m: usize,
    n: usize,
    k: usize,
) {
    for jc in 0..n {
        for l in 0..k {
            let temp = -ab[b0 + jc * ldb + l];
            for i in 0..m {
                let a = match aw {
                    Some(w) => w[a0 + l * lda + i],
                    None => ab[a0 + l * lda + i],
                };
                let p = c0 + jc * ldc + i;
                ab[p] = ab[p] + temp * a;
            }
        }
    }
}

/// [`gemm_nn_minus`] with `B` in `bw` (WORK13) and `A` in `ab` (`aw = None`)
/// or in `aw` (WORK31).
#[allow(clippy::too_many_arguments)]
fn gemm_nn_b<T: DirectScalar>(
    ab: &mut [T],
    a0: usize,
    lda: usize,
    bw: &[T],
    ldb: usize,
    c0: usize,
    ldc: usize,
    m: usize,
    n: usize,
    k: usize,
    aw: Option<&[T]>,
) {
    for jc in 0..n {
        for l in 0..k {
            let temp = -bw[jc * ldb + l];
            for i in 0..m {
                let a = match aw {
                    Some(w) => w[a0 + l * lda + i],
                    None => ab[a0 + l * lda + i],
                };
                let p = c0 + jc * ldc + i;
                ab[p] = ab[p] + temp * a;
            }
        }
    }
}

/// `DGBTF2`. Returns the pivot rows; `Err(Singular(j))` when `U[j,j]` is exactly
/// zero (LAPACK `INFO = j+1`, after completing the factorisation).
fn gbtf2<T: DirectScalar>(m: &mut BandMatrix<T>) -> Result<Vec<usize>, LinalgError> {
    let (n, kl, ku) = (m.n, m.kl, m.ku);
    let mut ipiv = vec![0usize; n];
    let mut info = None;
    let mut ju = 0usize;
    for j in 0..n {
        let km = kl.min(n - 1 - j);
        let mut jp = 0;
        let mut best = m.get(j, j).abs1();
        for t in 1..=km {
            let v = m.get(j + t, j).abs1();
            if v > best {
                best = v;
                jp = t;
            }
        }
        ipiv[j] = j + jp;
        if m.get(j + jp, j) != T::zero() {
            ju = ju.max((j + ku + jp).min(n - 1));
            if jp != 0 {
                for c in j..=ju {
                    let (p, q) = (m.idx(j + jp, c), m.idx(j, c));
                    m.ab.swap(p, q);
                }
            }
            if km > 0 {
                let r = T::one() / m.get(j, j);
                for t in 1..=km {
                    let k = m.idx(j + t, j);
                    m.ab[k] = m.ab[k] * r;
                }
                for c in j + 1..=ju {
                    let temp = -m.ab[m.idx(j, c)]; // may lie in the fill-in region, outside `get`'s band
                    for t in 1..=km {
                        let k = m.idx(j + t, c);
                        let l = m.get(j + t, j);
                        m.ab[k] = m.ab[k] + l * temp;
                    }
                }
            }
        } else if info.is_none() {
            info = Some(j);
        }
    }
    match info {
        Some(j) => Err(LinalgError::Singular(j)),
        None => Ok(ipiv),
    }
}

/// `DGBTRS('N')`.
fn gbtrs<T: DirectScalar>(lu: &BandMatrix<T>, ipiv: &[usize], b: &mut [T], nrhs: usize) {
    let (n, kl, ku) = (lu.n, lu.kl, lu.ku);
    if n == 0 {
        return;
    }
    for rhs in b.chunks_exact_mut(n).take(nrhs) {
        if kl > 0 {
            for j in 0..n - 1 {
                let lm = kl.min(n - 1 - j);
                let l = ipiv[j];
                if l != j {
                    rhs.swap(l, j);
                }
                let temp = -rhs[j];
                for t in 1..=lm {
                    rhs[j + t] = rhs[j + t] + lu.get(j + t, j) * temp;
                }
            }
        }
        // DTBSV('U','N','N') with k = kl + ku.
        let k = kl + ku;
        for j in (0..n).rev() {
            if rhs[j] != T::zero() {
                rhs[j] = rhs[j] / lu.ab[lu.idx(j, j)];
                let temp = rhs[j];
                for i in (j.saturating_sub(k)..j).rev() {
                    rhs[i] = rhs[i] - temp * lu.ab[lu.idx(i, j)];
                }
            }
        }
    }
}

/// Symmetric (Hermitian) band matrix, lower triangle only (`SBAND_INDEX`).
#[derive(Clone, Debug)]
pub struct SymBandMatrix<T: DirectScalar = f64> {
    n: usize,
    kd: usize,
    ab: Vec<T>,
}

impl<T: DirectScalar> SymBandMatrix<T> {
    /// Zero `n x n` matrix with `kd` sub-diagonals (Elmer's `Subband`).
    pub fn new(n: usize, kd: usize) -> Self {
        Self { n, kd, ab: vec![T::zero(); (kd + 1) * n] }
    }

    /// Lower triangle of `a` (entries with `col <= row`); the upper triangle is
    /// ignored, as Elmer only stores the lower half of a symmetric band matrix.
    pub fn from_csr(a: &CsrMatrix<T>) -> Result<Self, LinalgError> {
        if a.nrows != a.ncols {
            return Err(LinalgError::Dim(format!("band solve needs a square matrix, got {}x{}", a.nrows, a.ncols)));
        }
        check_csr(a)?;
        let mut kd = 0;
        for i in 0..a.nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k] as usize;
                if j >= a.ncols {
                    return Err(LinalgError::Dim(format!("column index {j} out of range")));
                }
                kd = kd.max(i.saturating_sub(j));
            }
        }
        let mut m = Self::new(a.nrows, kd);
        for i in 0..a.nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k] as usize;
                if j <= i {
                    m.add(i, j, a.values[k]);
                }
            }
        }
        Ok(m)
    }

    pub fn n(&self) -> usize {
        self.n
    }
    pub fn kd(&self) -> usize {
        self.kd
    }
    pub fn data(&self) -> &[T] {
        &self.ab
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * (self.kd + 1) + i - j
    }

    /// `A[i,j]` for `j <= i`, mirrored for `j > i` (zero outside the band).
    pub fn get(&self, i: usize, j: usize) -> T {
        let (r, c, mirrored) = if j <= i { (i, j, false) } else { (j, i, true) };
        if r >= self.n || r - c > self.kd {
            return T::zero();
        }
        let v = self.ab[self.idx(r, c)];
        if mirrored {
            v.conj()
        } else {
            v
        }
    }

    /// `A[i,j] += v`; like `Band_AddToMatrixElement`, entries above the diagonal are ignored.
    pub fn add(&mut self, i: usize, j: usize, v: T) {
        if j <= i {
            assert!(i < self.n && i - j <= self.kd, "SymBandMatrix::add({i},{j}) outside the band");
            let k = self.idx(i, j);
            self.ab[k] = self.ab[k] + v;
        }
    }

    /// `Band_SetMatrixElement`; entries above the diagonal are ignored.
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        if j <= i {
            assert!(i < self.n && i - j <= self.kd, "SymBandMatrix::set({i},{j}) outside the band");
            let k = self.idx(i, j);
            self.ab[k] = v;
        }
    }

    /// `Band_ZeroRow` for the symmetric format (lower part of the row only).
    pub fn zero_row(&mut self, n: usize) {
        for j in n.saturating_sub(self.kd)..=n {
            let k = self.idx(n, j);
            self.ab[k] = T::zero();
        }
    }

    /// `v = A u` (`Band_MatrixVectorMultiply`, symmetric branch; the upper
    /// triangle is the mirror, without conjugation as in Elmer).
    pub fn matvec(&self, u: &[T], v: &mut [T]) {
        for i in 0..self.n {
            let mut s = T::zero();
            for j in i.saturating_sub(self.kd)..=i {
                s = s + u[j] * self.ab[self.idx(i, j)];
            }
            for j in i + 1..=(i + self.kd).min(self.n - 1) {
                s = s + u[j] * self.ab[self.idx(j, i)];
            }
            v[i] = s;
        }
    }

    /// `SBand_SetDirichlet`: impose `x_n = value`, moving the column to the
    /// right-hand side and keeping symmetry.
    pub fn set_dirichlet(&mut self, b: &mut [T], n: usize, value: T) {
        for j in n.saturating_sub(self.kd)..n {
            let k = self.idx(n, j);
            b[j] = b[j] - value * self.ab[k];
            self.ab[k] = T::zero();
        }
        for j in n + 1..=(n + self.kd).min(self.n - 1) {
            let k = self.idx(j, n);
            b[j] = b[j] - value * self.ab[k];
            self.ab[k] = T::zero();
        }
        b[n] = value;
        let k = self.idx(n, n);
        self.ab[k] = T::one();
    }
}

/// Band Cholesky `A = L L^H`: `DPBTRF('L')`/`DPBTRS` (`ZPBTRF`/`ZPBTRS`), the
/// `symmetric banded` direct method.
#[derive(Clone, Debug)]
pub struct BandCholesky<T: DirectScalar = f64> {
    l: Option<SymBandMatrix<T>>,
}

impl<T: DirectScalar> Default for BandCholesky<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: DirectScalar> BandCholesky<T> {
    pub fn new() -> Self {
        Self { l: None }
    }

    /// Factor in place (`DPBTF2`). `Err(Singular(j))` when the leading minor of
    /// order `j+1` is not positive definite.
    pub fn factor_band(&mut self, mut m: SymBandMatrix<T>) -> Result<(), LinalgError> {
        self.l = None;
        pbtf2(&mut m)?;
        self.l = Some(m);
        Ok(())
    }

    pub fn factors(&self) -> Option<&SymBandMatrix<T>> {
        self.l.as_ref()
    }

    pub fn solve_many(&self, b: &mut [T], nrhs: usize) -> Result<(), LinalgError> {
        let l = self
            .l
            .as_ref()
            .ok_or_else(|| LinalgError::Unsupported("solve before factor".into()))?;
        if b.len() != l.n * nrhs {
            return Err(LinalgError::Dim(format!("rhs length {} != n*nrhs = {}", b.len(), l.n * nrhs)));
        }
        pbtrs(l, b, nrhs);
        Ok(())
    }
}

impl<T: DirectScalar> DirectSolver<T> for BandCholesky<T> {
    fn factor(&mut self, a: &CsrMatrix<T>) -> Result<(), LinalgError> {
        self.factor_band(SymBandMatrix::from_csr(a)?)
    }

    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError> {
        self.solve_many(b, 1)
    }
}

/// `DPBTF2('L')` / `ZPBTF2('L')`.
fn pbtf2<T: DirectScalar>(m: &mut SymBandMatrix<T>) -> Result<(), LinalgError> {
    let (n, kd) = (m.n, m.kd);
    for j in 0..n {
        let d = m.ab[m.idx(j, j)].re();
        if d <= 0.0 || d.is_nan() {
            return Err(LinalgError::Singular(j));
        }
        let ajj = d.sqrt();
        let kj = m.idx(j, j);
        m.ab[kj] = T::from_f64(ajj);
        let kn = kd.min(n - 1 - j);
        if kn > 0 {
            let r = 1.0 / ajj;
            for t in 1..=kn {
                let k = m.idx(j + t, j);
                m.ab[k] = m.ab[k] * T::from_f64(r);
            }
            // ZHER/DSYR('L', kn, -1, x, A(j+1.., j+1..)).
            for c in 0..kn {
                let xc = m.ab[m.idx(j + 1 + c, j)];
                if xc != T::zero() {
                    let temp = -xc.conj();
                    for r in c..kn {
                        let xr = m.ab[m.idx(j + 1 + r, j)];
                        let k = m.idx(j + 1 + r, j + 1 + c);
                        if r == c {
                            m.ab[k] = T::from_f64(m.ab[k].re() + (xr * temp).re());
                        } else {
                            m.ab[k] = m.ab[k] + xr * temp;
                        }
                    }
                } else {
                    let k = m.idx(j + 1 + c, j + 1 + c);
                    m.ab[k] = T::from_f64(m.ab[k].re());
                }
            }
        }
    }
    Ok(())
}

/// `DPBTRS('L')` / `ZPBTRS('L')`: forward `L y = b`, then backward `L^H x = y`.
fn pbtrs<T: DirectScalar>(l: &SymBandMatrix<T>, b: &mut [T], nrhs: usize) {
    let (n, kd) = (l.n, l.kd);
    if n == 0 {
        return;
    }
    for x in b.chunks_exact_mut(n).take(nrhs) {
        for j in 0..n {
            if x[j] != T::zero() {
                x[j] = x[j] / l.ab[l.idx(j, j)];
                let temp = x[j];
                for i in j + 1..=(j + kd).min(n - 1) {
                    x[i] = x[i] - temp * l.ab[l.idx(i, j)];
                }
            }
        }
        for j in (0..n).rev() {
            let mut temp = x[j];
            for i in (j + 1..=(j + kd).min(n - 1)).rev() {
                temp = temp - l.ab[l.idx(i, j)].conj() * x[i];
            }
            x[j] = temp / l.ab[l.idx(j, j)].conj();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// xorshift64, the generator `tests/data/gbtrf_lapack.txt` was made with.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 * 2f64.powi(-53) * 2.0 - 1.0
        }
    }

    fn fnv(v: &[f64]) -> u64 {
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        for x in v {
            for b in x.to_bits().to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x100_0000_01b3);
            }
        }
        h
    }

    /// `DGBTRF` and `DGBTRS` against the system reference LAPACK 3.12 the
    /// pristine Elmer links, bit for bit: the factor array (as an FNV-1a
    /// digest of its bits), the pivots, `INFO` and one solve. Each line of
    /// the data file is `seed n kl ku zero_column info digest ipiv [x]` for a
    /// random band matrix (60% fill, diagonal scaled by 10, column
    /// `zero_column` all zero when >= 0) and right-hand side drawn in that
    /// order. The cases cover the blocked path (`KU > 64`, `NB = 32 <= KL`),
    /// `KL != KU`, `N < KL`, a zero pivot and `DGBTF2` (`KU = 64`).
    #[test]
    fn gbtrf_matches_the_system_lapack_bit_for_bit() {
        let data = include_str!("../tests/data/gbtrf_lapack.txt");
        let mut blocked = 0;
        for line in data.lines() {
            let t: Vec<&str> = line.split_whitespace().collect();
            let p = |k: usize| t[k].parse::<i64>().unwrap();
            let (seed, n, kl, ku, zc) = (
                p(0) as u64,
                p(1) as usize,
                p(2) as usize,
                p(3) as usize,
                p(4),
            );
            let mut r = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut m = BandMatrix::<f64>::new(n, kl, ku);
            for j in 0..n {
                for i in j.saturating_sub(ku)..n.min(j + kl + 1) {
                    let keep = r.next() % 10 < 6;
                    let mut v = r.unit();
                    if i == j {
                        v *= 10.0;
                    }
                    if (keep || i == j) && j as i64 != zc {
                        m.set(i, j, v);
                    }
                }
            }
            let mut b: Vec<f64> = (0..n).map(|_| r.unit()).collect();
            blocked += usize::from(ku > 64 && kl >= 32);
            let res = gbtrf(&mut m);
            let ipiv: Vec<String> = match &res {
                Ok(ip) => ip.iter().map(|q| (q + 1).to_string()).collect(),
                Err(_) => Vec::new(),
            };
            assert_eq!(
                format!("{:016x}", fnv(m.data())),
                t[6],
                "factor array, case {seed}"
            );
            let info = p(5);
            if info != 0 {
                assert!(
                    matches!(res, Err(LinalgError::Singular(j)) if j as i64 == info - 1),
                    "INFO, case {seed}"
                );
                continue;
            }
            assert_eq!(ipiv.join(","), t[7], "pivots, case {seed}");
            let ip = res.unwrap();
            gbtrs(&m, &ip, &mut b, 1);
            let x: Vec<String> = b.iter().map(|x| format!("{:016x}", x.to_bits())).collect();
            assert_eq!(x.join(","), t[8], "solution, case {seed}");
        }
        assert!(blocked >= 5);
    }
}
