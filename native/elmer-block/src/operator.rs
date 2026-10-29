//! `BlockMatrixVectorProd` (BlockSolve.F90:2979).

use elmer_core::LinearOperator;

use crate::BlockSystem;

/// `v = A u` over the block-ordered vector: row block `i` sums
/// `A(i,j) u(j)` over `j` in order, skipping absent submatrices and the rows
/// of a `Dummy block i`.
pub struct BlockOperator<'a> {
    sys: &'a BlockSystem,
    dummy: Vec<bool>,
}

impl<'a> BlockOperator<'a> {
    /// `dummy[i]` is `Dummy block i+1` as the outer iteration reads it.
    pub fn new(sys: &'a BlockSystem, dummy: Vec<bool>) -> Self {
        BlockOperator { sys, dummy }
    }
}

impl LinearOperator<f64> for BlockOperator<'_> {
    fn nrows(&self) -> usize {
        self.sys.size()
    }

    fn ncols(&self) -> usize {
        self.sys.size()
    }

    fn apply(&self, u: &[f64], v: &mut [f64]) {
        let sys = self.sys;
        v[..sys.size()].fill(0.0);
        let mut s = Vec::new();
        for (i, &dummy) in self.dummy.iter().enumerate() {
            if dummy {
                continue;
            }
            let ri = sys.range(i);
            for j in 0..sys.nblocks() {
                let Some(a) = sys.block(i, j) else { continue };
                s.clear();
                s.resize(a.nrows, 0.0);
                elmer_sparse::spmv(a, &u[sys.range(j)], &mut s);
                for (vk, sk) in v[ri.clone()].iter_mut().zip(&s) {
                    *vk += sk;
                }
            }
        }
    }
}
