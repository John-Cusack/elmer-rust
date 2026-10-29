//! Splitting a monolithic system into blocks: `BlockPickMatrix`
//! (BlockSolve.F90:719) with `CRS_BlockMatrixPick` and `CRS_BlockMatrixPick2`
//! (CRSMatrix.F90:3409, :3682), and the block/monolithic vector permutation
//! (`BlockPerm`).
//!
//! The monolithic unknowns are node-major, component-minor (`comps`
//! components per node). Each block collects some components; within a block
//! the unknowns stay node-major in the order the components appear. Without
//! `Block Structure` every component is its own block.

use elmer_core::{CsrMatrix, Params};

use crate::BlockError;

/// Which block each component of the monolithic system belongs to.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BlockLayout {
    /// Components per node (`SIZE(BlockStruct)`, else the variable dofs).
    pub comps: usize,
    /// Block of each component, 0-based (`BlockStruct - 1`).
    pub block_of: Vec<usize>,
    /// Number of blocks (`NoVar`).
    pub nblocks: usize,
    /// Given by `Block Structure` (`GotBlockStruct`).
    pub given: bool,
}

impl BlockLayout {
    /// One block per component (`BlockInitMatrix` without `Block Structure`).
    pub fn by_component(comps: usize) -> Self {
        BlockLayout {
            comps,
            block_of: (0..comps).collect(),
            nblocks: comps,
            given: false,
        }
    }

    /// `Block Structure` (BlockSolve.F90:332-358): one 1-based block number
    /// per component, each in `1..=comps`.
    pub fn from_structure(structure: &[i32]) -> Result<Self, BlockError> {
        let comps = structure.len();
        if structure.iter().any(|&b| b < 1 || b as usize > comps) {
            return Err(BlockError::fatal(
                "BlockInitMatrix",
                "Incompatible values in > Block Structure < given!",
            ));
        }
        let nblocks = structure.iter().copied().max().unwrap_or(0) as usize;
        let block_of: Vec<usize> = structure.iter().map(|&b| b as usize - 1).collect();
        if (0..nblocks).any(|k| !block_of.contains(&k)) {
            return Err(BlockError::fatal(
                "BlockInitMatrix",
                "Invalid > Block Structure < given!",
            ));
        }
        Ok(BlockLayout {
            comps,
            block_of,
            nblocks,
            given: true,
        })
    }

    /// The layout a solver section asks for: `Block Structure` when given
    /// (its size must equal `dofs`), else one block per component.
    pub fn from_params(p: &dyn Params, dofs: usize) -> Result<Self, BlockError> {
        match p.get("block structure") {
            Some(elmer_core::Value::Integer(s)) => {
                if s.len() != dofs {
                    return Err(BlockError::fatal(
                        "BlockInitMatrix",
                        "Incompatible size of > Block Structure < given!",
                    ));
                }
                Self::from_structure(s)
            }
            Some(_) => Err(BlockError::fatal(
                "BlockInitMatrix",
                "Block Structure is not an integer array",
            )),
            None => Ok(Self::by_component(dofs)),
        }
    }

    /// Components of block `k`, in increasing order (`Irow`/`Icol`).
    fn members(&self, k: usize) -> Vec<usize> {
        (0..self.comps).filter(|&c| self.block_of[c] == k).collect()
    }
}

/// A block system: `nblocks x nblocks` submatrices, split from a monolithic
/// system ([`split`](Self::split)) or assembled per block
/// ([`from_blocks`](Self::from_blocks)).
#[derive(Clone, Debug)]
pub struct BlockSystem {
    /// The components behind the blocks of a split system; `None` for a
    /// system assembled per block.
    pub layout: Option<BlockLayout>,
    nblocks: usize,
    /// Start of each block in the block-ordered vector; `nblocks + 1` entries.
    pub offsets: Vec<usize>,
    /// Row-major submatrices; `None` where the pick found no entries.
    pub blocks: Vec<Option<CsrMatrix<f64>>>,
    /// `perm[k]` is the monolithic index of block-ordered index `k`
    /// (`BlockPerm`, BlockSolve.F90:788-829, 0-based).
    pub perm: Vec<usize>,
}

impl BlockSystem {
    /// Picks every submatrix of `a`. With `eliminate_zero` (`Block Eliminate
    /// Zero Submatrices`) a submatrix whose absolute sum is below
    /// `sqrt(tiny)` is dropped.
    pub fn split(a: &CsrMatrix<f64>, layout: BlockLayout, eliminate_zero: bool) -> Result<Self, BlockError> {
        let n = a.nrows;
        if layout.nblocks < 2 {
            return Err(BlockError::fatal(
                "CRS_BlockMatrixPick",
                "No applicable to just one block!",
            ));
        }
        if a.ncols != n || !n.is_multiple_of(layout.comps) {
            return Err(BlockError::fatal(
                "BlockPickMatrix",
                format!(
                    "matrix of size {n}x{} does not split into {} components",
                    a.ncols, layout.comps
                ),
            ));
        }
        let nv = n / layout.comps;
        let members: Vec<Vec<usize>> = (0..layout.nblocks).map(|k| layout.members(k)).collect();
        let mut offsets = vec![0];
        for m in &members {
            offsets.push(offsets.last().unwrap() + m.len() * nv);
        }
        let mut blocks = Vec::with_capacity(layout.nblocks * layout.nblocks);
        for r in 0..layout.nblocks {
            for c in 0..layout.nblocks {
                let b = pick(a, &layout, &members[r], &members[c], r == c);
                let b = b.filter(|b| {
                    !(eliminate_zero && b.values.iter().map(|v| v.abs()).sum::<f64>() < f64::MIN_POSITIVE.sqrt())
                });
                blocks.push(b);
            }
        }
        // BlockPerm: block k holds its components node by node.
        let mut perm = vec![0; n];
        for (k, m) in members.iter().enumerate() {
            for (j, &comp) in m.iter().enumerate() {
                for i in 0..nv {
                    perm[offsets[k] + m.len() * i + j] = layout.comps * i + comp;
                }
            }
        }
        Ok(BlockSystem {
            nblocks: layout.nblocks,
            layout: Some(layout),
            offsets,
            blocks,
            perm,
        })
    }

    /// A system assembled block by block (`BlockSolver` without a
    /// `Procedure`): `blocks` row-major, `nblocks x nblocks`, every diagonal
    /// block present and square, the off-diagonal ones sized to fit. The
    /// block-ordered vector is the concatenation of the block variables.
    pub fn from_blocks(nblocks: usize, blocks: Vec<Option<CsrMatrix<f64>>>) -> Result<Self, BlockError> {
        if blocks.len() != nblocks * nblocks {
            return Err(BlockError::fatal(
                "BlockInitMatrix",
                format!("{} submatrices for {nblocks} blocks", blocks.len()),
            ));
        }
        let mut offsets = vec![0];
        for i in 0..nblocks {
            let Some(d) = &blocks[i * nblocks + i] else {
                return Err(BlockError::fatal(
                    "BlockKrylovIter",
                    format!("Submatrix {} not associated!", 11 * (i + 1)),
                ));
            };
            if d.nrows != d.ncols {
                return Err(BlockError::fatal(
                    "BlockInitMatrix",
                    format!("block {} is not square", 11 * (i + 1)),
                ));
            }
            offsets.push(offsets[i] + d.nrows);
        }
        for i in 0..nblocks {
            for j in 0..nblocks {
                if let Some(b) = &blocks[i * nblocks + j] {
                    if b.nrows != offsets[i + 1] - offsets[i] || b.ncols != offsets[j + 1] - offsets[j] {
                        return Err(BlockError::fatal(
                            "BlockInitMatrix",
                            format!("block {}{} does not fit its row and column blocks", i + 1, j + 1),
                        ));
                    }
                }
            }
        }
        let n = offsets[nblocks];
        Ok(BlockSystem {
            layout: None,
            nblocks,
            offsets,
            blocks,
            perm: (0..n).collect(),
        })
    }

    pub fn nblocks(&self) -> usize {
        self.nblocks
    }

    /// Size of the block-ordered vector (`TotSize`).
    pub fn size(&self) -> usize {
        *self.offsets.last().unwrap()
    }

    /// Submatrix `(i, j)`, 0-based.
    pub fn block(&self, i: usize, j: usize) -> Option<&CsrMatrix<f64>> {
        self.blocks[i * self.nblocks() + j].as_ref()
    }

    /// The rows of block `i` in a block-ordered vector.
    pub fn range(&self, i: usize) -> std::ops::Range<usize> {
        self.offsets[i]..self.offsets[i + 1]
    }

    /// Monolithic vector to block order.
    pub fn gather(&self, mono: &[f64]) -> Vec<f64> {
        self.perm.iter().map(|&m| mono[m]).collect()
    }

    /// Block-ordered vector back into the monolithic one.
    pub fn scatter(&self, blocked: &[f64], mono: &mut [f64]) {
        for (k, &m) in self.perm.iter().enumerate() {
            mono[m] = blocked[k];
        }
    }
}

/// `CRS_BlockMatrixPick2` for the components `rows` x `cols` (the plain
/// `CRS_BlockMatrixPick` is the one-component case). Rows and columns keep
/// node-major order; the diagonal index is set for a diagonal block whose
/// every row stores its diagonal.
fn pick(
    a: &CsrMatrix<f64>,
    layout: &BlockLayout,
    rows: &[usize],
    cols: &[usize],
    diagonal: bool,
) -> Option<CsrMatrix<f64>> {
    let comps = layout.comps;
    let nv = a.nrows / comps;
    let mut col_slot = vec![usize::MAX; comps];
    for (mc, &c) in cols.iter().enumerate() {
        col_slot[c] = mc;
    }
    let mut b = CsrMatrix {
        nrows: rows.len() * nv,
        ncols: cols.len() * nv,
        row_ptr: Vec::with_capacity(rows.len() * nv + 1),
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::new();
    let mut all_diag = diagonal;
    b.row_ptr.push(0);
    for isub in 0..nv {
        for (mr, &r) in rows.iter().enumerate() {
            let i = comps * isub + r;
            let row = rows.len() * isub + mr;
            let mut d = None;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let l = a.col_idx[k] as usize;
                let mc = col_slot[l % comps];
                if mc == usize::MAX {
                    continue;
                }
                let lsub = cols.len() * (l / comps) + mc;
                if diagonal && lsub == row {
                    d = Some(b.values.len());
                }
                b.col_idx.push(lsub as u32);
                b.values.push(a.values[k]);
            }
            b.row_ptr.push(b.values.len());
            match d {
                Some(d) => diag.push(d),
                None => all_diag = false,
            }
        }
    }
    if b.values.is_empty() {
        return None;
    }
    if all_diag {
        b.diag = Some(diag);
    }
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 3 nodes x 2 components with entry value 10*row + col (1-based).
    fn mono() -> CsrMatrix<f64> {
        let n = 6;
        let mut a = CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr: vec![0],
            col_idx: Vec::new(),
            values: Vec::new(),
            diag: None,
        };
        for i in 0..n {
            for j in 0..n {
                if (i / 2).abs_diff(j / 2) <= 1 {
                    a.col_idx.push(j as u32);
                    a.values.push((10 * (i + 1) + j + 1) as f64);
                }
            }
            a.row_ptr.push(a.values.len());
        }
        a
    }

    #[test]
    fn component_blocks_pick_strided_entries() {
        let a = mono();
        let s = BlockSystem::split(&a, BlockLayout::by_component(2), false).unwrap();
        let b11 = s.block(0, 0).unwrap();
        assert_eq!(b11.nrows, 3);
        assert_eq!(b11.row_ptr, vec![0, 2, 5, 7]);
        assert_eq!(b11.col_idx, vec![0, 1, 0, 1, 2, 1, 2]);
        assert_eq!(b11.values, vec![11.0, 13.0, 31.0, 33.0, 35.0, 53.0, 55.0]);
        assert_eq!(b11.diag, Some(vec![0, 3, 6]));
        let b21 = s.block(1, 0).unwrap();
        assert_eq!(b21.values, vec![21.0, 23.0, 41.0, 43.0, 45.0, 63.0, 65.0]);
        assert_eq!(b21.diag, None);
        assert_eq!(s.perm, vec![0, 2, 4, 1, 3, 5]);
        let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let g = s.gather(&x);
        assert_eq!(g, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
        let mut back = [0.0; 6];
        s.scatter(&g, &mut back);
        assert_eq!(back, x);
    }

    #[test]
    fn block_structure_groups_components() {
        // 2 nodes x 3 components, Block Structure = 1 1 2.
        let n = 6;
        let mut a = CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr: vec![0],
            col_idx: Vec::new(),
            values: Vec::new(),
            diag: None,
        };
        for i in 0..n {
            for j in 0..n {
                a.col_idx.push(j as u32);
                a.values.push((10 * (i + 1) + j + 1) as f64);
            }
            a.row_ptr.push(a.values.len());
        }
        let layout = BlockLayout::from_structure(&[1, 1, 2]).unwrap();
        let s = BlockSystem::split(&a, layout, false).unwrap();
        assert_eq!(s.offsets, vec![0, 4, 6]);
        assert_eq!(s.perm, vec![0, 1, 3, 4, 2, 5]);
        let b11 = s.block(0, 0).unwrap();
        // Row 0 is monolithic row 1 (node 1, comp 1): columns 1 2 4 5.
        assert_eq!(&b11.values[..4], &[11.0, 12.0, 14.0, 15.0]);
        assert_eq!(&b11.col_idx[..4], &[0, 1, 2, 3]);
        assert_eq!(b11.diag, Some(vec![0, 5, 10, 15]));
        let b12 = s.block(0, 1).unwrap();
        assert_eq!(b12.ncols, 2);
        assert_eq!(&b12.values[..2], &[13.0, 16.0]);
        // A dense pick of every block reproduces A x in block order.
        let x: Vec<f64> = (1..=6).map(f64::from).collect();
        let mut y = vec![0.0; 6];
        elmer_core::LinearOperator::apply(&a, &x, &mut y);
        let xb = s.gather(&x);
        let mut yb = vec![0.0; 6];
        for i in 0..2 {
            for j in 0..2 {
                let m = s.block(i, j).unwrap();
                let mut t = vec![0.0; m.nrows];
                elmer_core::LinearOperator::apply(m, &xb[s.range(j)], &mut t);
                for (k, v) in s.range(i).zip(t) {
                    yb[k] += v;
                }
            }
        }
        assert_eq!(s.gather(&y), yb);
    }

    #[test]
    fn bad_structures_are_fatal() {
        assert!(BlockLayout::from_structure(&[1, 3]).is_err());
        assert!(BlockLayout::from_structure(&[2, 2]).is_err());
    }
}
