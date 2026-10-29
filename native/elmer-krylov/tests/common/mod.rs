//! Fixture loading for the tests. Operators and preconditioners come from
//! elmer-sparse (lane L01), whose kernels keep Elmer's floating-point order.
#![allow(dead_code)]

use std::path::PathBuf;

use elmer_core::{CsrMatrix, Scalar};

pub fn fixture_dir(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(case)
}

/// 0-based `(row, col, value)` triplets into a square CSR matrix.
pub fn csr<T: Scalar>(n: usize, t: &[(usize, u32, T)]) -> CsrMatrix<T> {
    elmer_sparse::ops::from_triplets(n, n, t.iter().copied())
}

/// A `Linear System Save` dump: 1-based `i j value` lines and one rhs value per row.
pub fn load_system(case: &str) -> (CsrMatrix<f64>, Vec<f64>) {
    let dir = fixture_dir(case);
    let b: Vec<f64> = std::fs::read_to_string(dir.join("b.txt"))
        .unwrap()
        .lines()
        .map(|l| l.trim().parse().unwrap())
        .collect();
    let t: Vec<(usize, u32, f64)> = std::fs::read_to_string(dir.join("a.txt"))
        .unwrap()
        .lines()
        .map(|l| {
            let mut it = l.split_whitespace();
            let i: usize = it.next().unwrap().parse().unwrap();
            let j: u32 = it.next().unwrap().parse().unwrap();
            let v: f64 = it.next().unwrap().parse().unwrap();
            (i - 1, j - 1, v)
        })
        .collect();
    (csr(b.len(), &t), b)
}

pub fn load_lines(case: &str, run: &str) -> Vec<String> {
    std::fs::read_to_string(fixture_dir(case).join(format!("{run}.out")))
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}
