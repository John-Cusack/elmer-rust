#![allow(dead_code)]

use elmer_core::CsrMatrix;
use num_complex::Complex64;
use std::path::PathBuf;

pub struct Fixture {
    pub k: CsrMatrix<f64>,
    pub m: CsrMatrix<f64>,
    pub tight: Vec<f64>,
    pub asis: Vec<f64>,
    pub modes: Vec<Vec<f64>>,
}

fn dir(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn data_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
}

fn parse_f64(t: &str) -> f64 {
    t.replace(['D', 'd'], "E")
        .parse()
        .unwrap_or_else(|_| panic!("bad number {t}"))
}

fn csr(n: usize, row_ptr: Vec<usize>, col_idx: Vec<u32>, values: Vec<f64>) -> CsrMatrix<f64> {
    let diag = (0..n)
        .map(|i| {
            (row_ptr[i]..row_ptr[i + 1])
                .find(|&p| col_idx[p] as usize == i)
                .expect("diagonal")
        })
        .collect();
    CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr,
        col_idx,
        values,
        diag: Some(diag),
    }
}

pub fn fixture(name: &str) -> Fixture {
    let text = std::fs::read_to_string(dir(name).join("system.txt")).unwrap();
    let mut lines = data_lines(&text);
    let head: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect();
    let (n, nnz) = (head[0], head[1]);
    let mut row_ptr = vec![0usize; n + 1];
    let (mut col_idx, mut kv, mut mv) = (
        Vec::with_capacity(nnz),
        Vec::with_capacity(nnz),
        Vec::with_capacity(nnz),
    );
    let mut prev_row = 0;
    for line in lines {
        let t: Vec<&str> = line.split_whitespace().collect();
        let row: usize = t[0].parse().unwrap();
        assert!(row > prev_row, "rows ascending");
        for r in prev_row..row {
            row_ptr[r + 1] = col_idx.len();
        }
        for e in t[1..].chunks(3) {
            col_idx.push(e[0].parse::<u32>().unwrap() - 1);
            kv.push(parse_f64(e[1]));
            mv.push(parse_f64(e[2]));
        }
        row_ptr[row] = col_idx.len();
        prev_row = row;
    }
    for r in prev_row..n {
        row_ptr[r + 1] = col_idx.len();
    }
    assert_eq!(col_idx.len(), nnz);
    let k = csr(n, row_ptr.clone(), col_idx.clone(), kv);
    let m = csr(n, row_ptr, col_idx, mv);

    let ev = std::fs::read_to_string(dir(name).join("eigenvalues.txt")).unwrap();
    let mut tight = Vec::new();
    let mut asis = Vec::new();
    for line in data_lines(&ev) {
        let mut t = line.split_whitespace();
        let tag = t.next().unwrap();
        let vals: Vec<f64> = t.map(parse_f64).collect();
        match tag {
            "tight" => tight = vals,
            "asis" => asis = vals,
            _ => panic!("tag {tag}"),
        }
    }

    let mt = std::fs::read_to_string(dir(name).join("modes.txt")).unwrap();
    let mut modes: Vec<Vec<f64>> = Vec::new();
    for line in data_lines(&mt).skip(1) {
        if line.starts_with("mode") {
            modes.push(Vec::with_capacity(n));
        } else {
            modes.last_mut().unwrap().push(parse_f64(line.trim()));
        }
    }
    assert!(modes.iter().all(|x| x.len() == n));
    Fixture {
        k,
        m,
        tight,
        asis,
        modes,
    }
}

pub fn re(x: &[Complex64]) -> Vec<f64> {
    x.iter().map(|z| z.re).collect()
}

/// Max relative difference after choosing the sign of `x` that best matches `y`.
pub fn sign_normalised_diff(x: &[f64], y: &[f64]) -> f64 {
    let scale = y.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let d = |s: f64| {
        x.iter()
            .zip(y)
            .map(|(a, b)| (s * a - b).abs())
            .fold(0.0, f64::max)
            / scale
    };
    d(1.0).min(d(-1.0))
}

pub struct Lcg(pub u64);

impl Lcg {
    pub fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }
}

/// Tridiagonal CSR with constant bands `(sub, diag, sup)`.
pub fn tridiag(n: usize, sub: f64, diag: f64, sup: f64) -> CsrMatrix<f64> {
    let (mut rp, mut ci, mut v) = (vec![0], Vec::new(), Vec::new());
    for i in 0..n {
        if i > 0 {
            ci.push((i - 1) as u32);
            v.push(sub);
        }
        ci.push(i as u32);
        v.push(diag);
        if i + 1 < n {
            ci.push((i + 1) as u32);
            v.push(sup);
        }
        rp.push(ci.len());
    }
    csr(n, rp, ci, v)
}
