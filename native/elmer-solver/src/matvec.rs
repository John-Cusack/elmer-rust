//! The matrix-vector product the Krylov methods see: `CRS_MatrixVectorProd`
//! (CRSMatrix.F90), i.e. `CRS_MatrixVectorMultiply` (CRSMatrix.F90:1648) in
//! its non-MKL form.
//!
//! Elmer does not sum a row left to right. A matrix from `CRS_CreateMatrix`
//! carries `A % ndeg`, the dofs of its variable, and the product steps
//! through each row a block at a time with one accumulator per block
//! column (`ndeg` 2 to 6, 8 and 10; 8 and 10 step by 4 and 5); any other
//! matrix (`ndeg` -1: EliminateDirichlet's reduced system, the constraint
//! collection matrix, list matrices) has its rows unrolled by four. The
//! partial sums are added at the end, so the rounding differs from a plain
//! row sum, and a Krylov method run to a loose tolerance carries that
//! difference into the solution.
//!
//! `elmer-sparse` has both loops ([`elmer_sparse::spmv`],
//! [`elmer_sparse::spmv_blocked`]); this picks the one the matrix's `ndeg`
//! selects.

use elmer_core::{CsrMatrix, LinalgError, LinearOperator};

/// `CRS_MatrixVectorMultiply( A, u, v )`: `v = A u` with Elmer's product
/// for a matrix of `A % ndeg = ndeg`.
pub fn crs_matrix_vector_multiply(a: &CsrMatrix<f64>, ndeg: usize, u: &[f64], v: &mut [f64]) {
    ElmerMatvec::new(a, ndeg).apply(u, v);
}

/// `A` with Elmer's product. `ndeg` is the matrix's `A % ndeg` (0 when it
/// has none).
pub(crate) struct ElmerMatvec<'a> {
    a: &'a CsrMatrix<f64>,
    /// The `ndeg` the product uses: 1 (the unrolled-by-four default) unless
    /// the matrix holds whole blocks.
    ndeg: usize,
}

impl<'a> ElmerMatvec<'a> {
    pub(crate) fn new(a: &'a CsrMatrix<f64>, ndeg: usize) -> Self {
        let step = match ndeg {
            5 | 10 => 5,
            4 | 8 => 4,
            3 | 6 => 3,
            2 => 2,
            _ => 0,
        };
        // The blocked loops take the columns of a step to be consecutive
        // (`u(l+1)` for `Cols(j+1)`), which holds for every matrix
        // InitializeMatrix builds; Elmer only checks the row lengths
        // (CRSMatrix.F90:1702) and stops on a partial block. A matrix
        // here without whole blocks gets the default loop.
        let ndeg = if step > 0 && whole_blocks(a, ndeg, step) {
            ndeg
        } else {
            1
        };
        ElmerMatvec { a, ndeg }
    }
}

fn whole_blocks(a: &CsrMatrix<f64>, ndeg: usize, step: usize) -> bool {
    (0..a.nrows).all(|i| {
        let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
        (hi - lo) % ndeg == 0
            && (lo..hi).step_by(step).all(|j| {
                let c = a.col_idx[j];
                (1..step).all(|k| a.col_idx[j + k] == c + k as u32)
            })
    })
}

impl LinearOperator<f64> for ElmerMatvec<'_> {
    fn nrows(&self) -> usize {
        self.a.nrows
    }

    fn ncols(&self) -> usize {
        self.a.ncols
    }

    fn apply(&self, u: &[f64], v: &mut [f64]) {
        if self.ndeg > 1 {
            elmer_sparse::spmv_blocked(self.a, self.ndeg, u, v)
                .expect("whole blocks were checked in ElmerMatvec::new");
        } else {
            elmer_sparse::spmv(self.a, u, v);
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) -> Result<(), LinalgError> {
        self.a.apply_transpose(x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn csr(rows: &[&[(u32, f64)]]) -> CsrMatrix<f64> {
        let mut a = CsrMatrix {
            nrows: rows.len(),
            ncols: 0,
            row_ptr: vec![0],
            col_idx: Vec::new(),
            values: Vec::new(),
            diag: None,
        };
        for r in rows {
            for &(c, v) in *r {
                a.col_idx.push(c);
                a.values.push(v);
                a.ncols = a.ncols.max(c as usize + 1);
            }
            a.row_ptr.push(a.values.len());
        }
        a
    }

    /// Five entries chosen so that the left-to-right row sum, the
    /// four-accumulator sum and the block sums all round differently.
    const ROW: [(u32, f64); 6] = [
        (0, 1.0),
        (1, 1.0e16),
        (2, 1.0),
        (3, -1.0e16),
        (4, 1.0),
        (5, 1.0),
    ];

    #[test]
    fn default_loop_unrolls_by_four() {
        let a = csr(&[&ROW]);
        let u = vec![1.0; 6];
        let mut v = [0.0];
        ElmerMatvec::new(&a, 0).apply(&u, &mut v);
        // r1 = 1 + 1 + 1 (the remainder joins r1) = 3, r2 = 1e16, r3 = 1,
        // r4 = -1e16: ((3 + 1e16) + 1) - 1e16 = 4 in binary64.
        let (r1, r2, r3, r4) = (3.0f64, 1.0e16f64, 1.0f64, -1.0e16f64);
        assert_eq!(v[0], r1 + r2 + r3 + r4);
        // The plain row sum rounds differently, so the test bites.
        let mut plain = [0.0];
        a.apply(&u, &mut plain);
        assert_ne!(v[0], plain[0]);
    }

    #[test]
    fn blocked_loop_keeps_one_sum_per_block_column() {
        let a = csr(&[&ROW]);
        let u = vec![1.0; 6];
        let mut v = [0.0];
        ElmerMatvec::new(&a, 3).apply(&u, &mut v);
        // Blocks (0,1,2) and (3,4,5): r1 = 1 - 1e16, r2 = 1e16 + 1, r3 = 2.
        let (r1, r2, r3) = (1.0f64 - 1.0e16, 1.0e16f64 + 1.0, 2.0f64);
        assert_eq!(v[0], r1 + r2 + r3);
        let mut four = [0.0];
        ElmerMatvec::new(&a, 0).apply(&u, &mut four);
        assert_ne!(v[0], four[0]);
    }

    #[test]
    fn partial_blocks_fall_back_to_the_default_loop() {
        // Five entries with ndeg = 2, and a gap inside a block with ndeg = 3.
        let short = csr(&[&ROW[..5]]);
        let gap = csr(&[&[
            (0, 1.0),
            (2, 1.0e16),
            (3, 1.0),
            (4, -1.0e16),
            (5, 1.0),
            (6, 1.0),
        ]]);
        let u = vec![1.0; 7];
        for (a, ndeg) in [(&short, 2), (&gap, 3)] {
            let (mut v, mut d) = ([0.0], [0.0]);
            ElmerMatvec::new(a, ndeg).apply(&u, &mut v);
            ElmerMatvec::new(a, 0).apply(&u, &mut d);
            assert_eq!(v, d);
        }
    }

    #[test]
    fn products_agree_with_the_matrix() {
        let a = csr(&[
            &[(0, 2.0), (1, -1.0)],
            &[(0, -1.0), (1, 2.0), (2, -1.0)],
            &[(1, -1.0), (2, 2.0)],
        ]);
        let u = [1.0, 2.0, 3.0];
        let (mut v, mut w) = ([0.0; 3], [0.0; 3]);
        ElmerMatvec::new(&a, 1).apply(&u, &mut v);
        a.apply(&u, &mut w);
        assert_eq!(v, w);
    }
}
