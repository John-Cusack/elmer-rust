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
//! [`BandLu`] is `DGBTF2`/`DGBTRS` (`ZGBTF2`/`ZGBTRS` for complex) and
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
        let ipiv = gbtf2(&mut m)?;
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
        self.factor_band(BandMatrix::from_csr(a)?)
    }

    fn solve_in_place(&self, b: &mut [T]) -> Result<(), LinalgError> {
        self.solve_many(b, 1)
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
