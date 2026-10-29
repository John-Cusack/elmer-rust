#![allow(dead_code)]
use elmer_core::CsrMatrix;
use num_complex::Complex64;

/// Deterministic generator (no rand dependency).
pub struct Lcg(pub u64);
impl Lcg {
    pub fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    }
}

/// 2-D 5-point Laplacian on an `m x m` grid (SPD, sparse).
pub fn laplacian(m: usize) -> CsrMatrix<f64> {
    let n = m * m;
    let (mut row_ptr, mut col_idx, mut values) = (vec![0usize], vec![], vec![]);
    for r in 0..m {
        for c in 0..m {
            let i = r * m + c;
            let mut push = |j: usize, v: f64| {
                col_idx.push(j as u32);
                values.push(v);
            };
            if r > 0 {
                push(i - m, -1.0);
            }
            if c > 0 {
                push(i - 1, -1.0);
            }
            push(i, 4.0);
            if c + 1 < m {
                push(i + 1, -1.0);
            }
            if r + 1 < m {
                push(i + m, -1.0);
            }
            row_ptr.push(col_idx.len());
        }
    }
    CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values, diag: None }
}

/// Banded unsymmetric matrix with the given bandwidths and a weak diagonal (forces pivoting).
pub fn random_band(n: usize, kl: usize, ku: usize, seed: u64) -> CsrMatrix<f64> {
    let mut g = Lcg(seed);
    let (mut row_ptr, mut col_idx, mut values) = (vec![0usize], vec![], vec![]);
    for i in 0..n {
        for j in i.saturating_sub(kl)..=(i + ku).min(n - 1) {
            col_idx.push(j as u32);
            values.push(if i == j { 0.1 * g.next_f64() } else { g.next_f64() });
        }
        row_ptr.push(col_idx.len());
    }
    CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values, diag: None }
}

pub fn to_complex(a: &CsrMatrix<f64>, seed: u64) -> CsrMatrix<Complex64> {
    let mut g = Lcg(seed);
    CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: a.values.iter().map(|&v| Complex64::new(v, 0.3 * g.next_f64())).collect(),
        diag: None,
    }
}

pub fn dense_of(a: &CsrMatrix<f64>) -> Vec<f64> {
    let mut d = vec![0.0; a.nrows * a.ncols];
    for i in 0..a.nrows {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            d[a.col_idx[k] as usize * a.nrows + i] += a.values[k];
        }
    }
    d
}

pub fn spmv<T: elmer_core::Scalar>(a: &CsrMatrix<T>, x: &[T]) -> Vec<T> {
    use elmer_core::LinearOperator;
    let mut y = vec![T::zero(); a.nrows];
    a.apply(x, &mut y);
    y
}

pub fn rel_residual<T: elmer_core::Scalar>(a: &CsrMatrix<T>, x: &[T], b: &[T]) -> f64 {
    let ax = spmv(a, x);
    let num: f64 = ax.iter().zip(b).map(|(p, q)| { let d = *p - *q; elmer_core::Scalar::abs(d).powi(2) }).sum();
    let den: f64 = b.iter().map(|q| elmer_core::Scalar::abs(*q).powi(2)).sum();
    (num / den).sqrt()
}

/// A system dumped by Elmer's `Linear System Save` (`linsys_{a,b,sol}.dat`).
pub struct Dump {
    pub a: CsrMatrix<f64>,
    pub b: Vec<f64>,
    pub sol: Option<Vec<f64>>,
}

/// Read a `Linear System Save` (`linsys_{a,b,sol}.dat`) from `dir` with
/// `elmer-sparse`'s reader.
pub fn read_dump(dir: &std::path::Path) -> Dump {
    let s = elmer_sparse::save::read_linear_system_save(dir, "linsys", elmer_sparse::save::Columns::default())
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    Dump { a: s.stiff.expect("linsys_a.dat"), b: s.rhs.expect("linsys_b.dat"), sol: s.solution }
}

pub fn fixture(name: &str) -> Dump {
    read_dump(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name))
}

pub fn is_symmetric(a: &CsrMatrix<f64>) -> bool {
    let d = dense_of(a);
    let n = a.nrows;
    (0..n).all(|i| (0..n).all(|j| d[j * n + i] == d[i * n + j]))
}
