//! Test-only operators that reproduce Elmer's floating-point order, and the
//! fixture loader. The real CSR kernels belong to elmer-sparse (lane L01).
#![allow(dead_code)]

use std::path::PathBuf;

use elmer_core::{LinalgError, LinearOperator, Preconditioner, Scalar};

/// CSR with `CRS_MatrixVectorMultiply`'s default (`ndeg = 1`) kernel:
/// four interleaved row accumulators summed as `r1 + r2 + r3 + r4`
/// (CRSMatrix.F90:1781-1797).
pub struct ElmerCsr<T> {
    pub n: usize,
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub vals: Vec<T>,
}

impl<T: Scalar> ElmerCsr<T> {
    /// From 0-based `(row, col, value)` triplets already in row-major storage order.
    pub fn from_triplets(n: usize, t: &[(usize, usize, T)]) -> Self {
        let mut rows = vec![0; n + 1];
        for &(i, _, _) in t {
            rows[i + 1] += 1;
        }
        for i in 0..n {
            rows[i + 1] += rows[i];
        }
        ElmerCsr {
            n,
            rows,
            cols: t.iter().map(|e| e.1).collect(),
            vals: t.iter().map(|e| e.2).collect(),
        }
    }

    pub fn diag(&self) -> Vec<T> {
        (0..self.n)
            .map(|i| {
                (self.rows[i]..self.rows[i + 1])
                    .find(|&j| self.cols[j] == i)
                    .map_or(T::zero(), |j| self.vals[j])
            })
            .collect()
    }
}

impl<T: Scalar> LinearOperator<T> for ElmerCsr<T> {
    fn nrows(&self) -> usize {
        self.n
    }
    fn ncols(&self) -> usize {
        self.n
    }
    fn apply(&self, u: &[T], v: &mut [T]) {
        for i in 0..self.n {
            let (lo, hi) = (self.rows[i], self.rows[i + 1]);
            let k = hi - lo;
            let (mut r1, mut r2, mut r3, mut r4) = (T::zero(), T::zero(), T::zero(), T::zero());
            let mut j = lo;
            while j < lo + 4 * (k / 4) {
                r1 = r1 + u[self.cols[j]] * self.vals[j];
                r2 = r2 + u[self.cols[j + 1]] * self.vals[j + 1];
                r3 = r3 + u[self.cols[j + 2]] * self.vals[j + 2];
                r4 = r4 + u[self.cols[j + 3]] * self.vals[j + 3];
                j += 4;
            }
            while j < hi {
                r1 = r1 + u[self.cols[j]] * self.vals[j];
                j += 1;
            }
            v[i] = r1 + r2 + r3 + r4;
        }
    }
    /// The `HUTI_MAT_TRPSED` branch of `CRS_MatrixVectorProd` (CRSMatrix.F90).
    fn apply_transpose(&self, u: &[T], v: &mut [T]) -> Result<(), LinalgError> {
        v.fill(T::zero());
        for i in 0..self.n {
            let s = u[i];
            for j in self.rows[i]..self.rows[i + 1] {
                v[self.cols[j]] = v[self.cols[j]] + s * self.vals[j];
            }
        }
        Ok(())
    }
}

/// `CRS_DiagPrecondition`: `u = v / a_ii` where `|a_ii| > AEPS`, else `u = v`.
pub struct DiagPrecond<T>(pub Vec<T>);

impl<T: Scalar> Preconditioner<T> for DiagPrecond<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        // AEPS = 10 * EPSILON(1.0_dp) (Types.F90).
        let aeps = 10.0 * f64::EPSILON;
        for i in 0..r.len() {
            z[i] = if self.0[i].abs() > aeps {
                r[i] / self.0[i]
            } else {
                r[i]
            };
        }
    }
}

pub fn fixture_dir(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(case)
}

/// A `Linear System Save` dump: 1-based `i j value` lines and one rhs value per row.
pub fn load_system(case: &str) -> (ElmerCsr<f64>, Vec<f64>) {
    let dir = fixture_dir(case);
    let b: Vec<f64> = std::fs::read_to_string(dir.join("b.txt"))
        .unwrap()
        .lines()
        .map(|l| l.trim().parse().unwrap())
        .collect();
    let t: Vec<(usize, usize, f64)> = std::fs::read_to_string(dir.join("a.txt"))
        .unwrap()
        .lines()
        .map(|l| {
            let mut it = l.split_whitespace();
            let i: usize = it.next().unwrap().parse().unwrap();
            let j: usize = it.next().unwrap().parse().unwrap();
            let v: f64 = it.next().unwrap().parse().unwrap();
            (i - 1, j - 1, v)
        })
        .collect();
    (ElmerCsr::from_triplets(b.len(), &t), b)
}

pub fn load_lines(case: &str, run: &str) -> Vec<String> {
    std::fs::read_to_string(fixture_dir(case).join(format!("{run}.out")))
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}
