mod common;
use common::*;

use elmer_core::{CsrMatrix, DirectSolver, LinalgError, ParamMap, Value};
use elmer_direct::{
    kind_for_method, select_direct, try_select_direct, BandCholesky, BandLu, DirectKind,
    FaerCholesky, FaerLu,
};
use num_complex::Complex64;

fn rhs(n: usize, seed: u64) -> Vec<f64> {
    let mut g = Lcg(seed);
    (0..n).map(|_| g.next_f64()).collect()
}

fn solve_with(s: &mut dyn DirectSolver<f64>, a: &CsrMatrix<f64>, b: &[f64]) -> Vec<f64> {
    s.factor(a).unwrap();
    let mut x = b.to_vec();
    s.solve_in_place(&mut x).unwrap();
    x
}

#[test]
fn sparse_lu_and_cholesky_solve_laplacian() {
    elmer_direct::set_single_thread();
    let a = laplacian(12);
    let b = rhs(a.nrows, 1);
    let x = solve_with(&mut FaerLu::new(), &a, &b);
    assert!(rel_residual(&a, &x, &b) < 1e-13);
    let y = solve_with(&mut FaerCholesky::new(), &a, &b);
    assert!(rel_residual(&a, &y, &b) < 1e-13);
    for (p, q) in x.iter().zip(&y) {
        assert!((p - q).abs() < 1e-11);
    }
}

#[test]
fn sparse_lu_unsymmetric_and_refactor_same_pattern() {
    let a = random_band(60, 3, 5, 7);
    let b = rhs(60, 2);
    let mut lu = FaerLu::new();
    let x = solve_with(&mut lu, &a, &b);
    assert!(rel_residual(&a, &x, &b) < 1e-9);
    // New values, same pattern: reuses the symbolic analysis.
    let mut a2 = a.clone();
    for v in &mut a2.values {
        *v *= 2.0;
    }
    lu.factor(&a2).unwrap();
    let mut x2 = b.clone();
    lu.solve_in_place(&mut x2).unwrap();
    for (p, q) in x.iter().zip(&x2) {
        assert!((p - 2.0 * q).abs() < 1e-8 * (1.0 + p.abs()));
    }
}

#[test]
fn sparse_lu_complex() {
    let a = to_complex(&laplacian(8), 3);
    let mut g = Lcg(5);
    let b: Vec<Complex64> = (0..a.nrows).map(|_| Complex64::new(g.next_f64(), g.next_f64())).collect();
    let mut lu = FaerLu::<Complex64>::new();
    lu.factor(&a).unwrap();
    let mut x = b.clone();
    lu.solve_in_place(&mut x).unwrap();
    assert!(rel_residual(&a, &x, &b) < 1e-12);
}

#[test]
fn cholesky_rejects_indefinite() {
    let mut a = laplacian(4);
    a.values[a.row_ptr[3] + 1] = -10.0; // diagonal of row 3 goes negative
    let r = FaerCholesky::<f64>::new().factor(&a);
    assert!(matches!(r, Err(LinalgError::Breakdown { .. })), "{r:?}");
}

#[test]
fn solve_before_factor_is_an_error() {
    let mut b = [1.0];
    assert!(FaerLu::<f64>::new().solve_in_place(&mut b).is_err());
    assert!(BandLu::<f64>::new().solve_in_place(&mut b).is_err());
}

#[test]
fn band_lu_matches_dense_pivoted_reference_bitwise() {
    // Reference: dense LU with the same pivot rule and operation order, run on
    // the full matrix (no band bookkeeping).
    let n = 25;
    let a = random_band(n, 3, 2, 11);
    let mut d = dense_of(&a); // column-major
    let mut ipiv = vec![0; n];
    for j in 0..n {
        let mut p = j;
        for i in j + 1..n {
            if d[j * n + i].abs() > d[j * n + p].abs() {
                p = i;
            }
        }
        ipiv[j] = p;
        if p != j {
            for c in 0..n {
                d.swap(c * n + j, c * n + p);
            }
        }
        let r = 1.0 / d[j * n + j];
        for i in j + 1..n {
            d[j * n + i] *= r;
        }
        for c in j + 1..n {
            let t = -d[c * n + j];
            for i in j + 1..n {
                d[c * n + i] += d[j * n + i] * t;
            }
        }
    }
    let mut lu = BandLu::<f64>::new();
    lu.factor(&a).unwrap();
    assert_eq!(lu.ipiv(), &ipiv[..]);
    let f = lu.factors().unwrap();
    for j in 0..n {
        // U only: LAPACK band storage keeps each column's multipliers in the
        // order they were computed (unswapped), unlike the dense layout.
        for i in j.saturating_sub(f.ku() + f.kl())..=j {
            assert_eq!(f.factor_entry(i, j).to_bits(), d[j * n + i].to_bits(), "entry ({i},{j})");
        }
    }
    // And the solve agrees with the residual.
    let b = rhs(n, 4);
    let mut x = b.clone();
    lu.solve_in_place(&mut x).unwrap();
    assert!(rel_residual(&a, &x, &b) < 1e-10);
}

#[test]
fn band_lu_complex_and_singular() {
    let a = to_complex(&random_band(30, 2, 2, 13), 9);
    let mut g = Lcg(8);
    let b: Vec<Complex64> = (0..30).map(|_| Complex64::new(g.next_f64(), g.next_f64())).collect();
    let mut lu = BandLu::<Complex64>::new();
    lu.factor(&a).unwrap();
    let mut x = b.clone();
    lu.solve_in_place(&mut x).unwrap();
    assert!(rel_residual(&a, &x, &b) < 1e-9);

    let mut s = random_band(6, 1, 1, 3);
    for k in s.row_ptr[2]..s.row_ptr[3] {
        s.values[k] = 0.0; // zero row => singular
    }
    assert!(matches!(BandLu::<f64>::new().factor(&s), Err(LinalgError::Singular(_))));
}

#[test]
fn band_cholesky_real_and_complex_hermitian() {
    let a = laplacian(6);
    let b = rhs(a.nrows, 6);
    let x = solve_with(&mut BandCholesky::new(), &a, &b);
    assert!(rel_residual(&a, &x, &b) < 1e-13);
    let y = solve_with(&mut BandLu::new(), &a, &b);
    for (p, q) in x.iter().zip(&y) {
        assert!((p - q).abs() < 1e-11);
    }

    // Hermitian positive definite: diagonally dominant with complex off-diagonals.
    let n: usize = 12;
    let (mut row_ptr, mut col_idx, mut values) = (vec![0usize], vec![], vec![]);
    for i in 0..n {
        for j in i.saturating_sub(2)..=(i + 2).min(n - 1) {
            col_idx.push(j as u32);
            values.push(if i == j {
                Complex64::new(10.0, 0.0)
            } else {
                let t = (i.min(j) * 7 + i.max(j)) as f64;
                let z = Complex64::new(t.sin(), t.cos());
                if j < i { z } else { z.conj() }
            });
        }
        row_ptr.push(col_idx.len());
    }
    let h = CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values, diag: None };
    let bb: Vec<Complex64> = (0..n).map(|i| Complex64::new(i as f64, 1.0)).collect();
    let mut ch = BandCholesky::<Complex64>::new();
    ch.factor(&h).unwrap();
    let mut z = bb.clone();
    ch.solve_in_place(&mut z).unwrap();
    assert!(rel_residual(&h, &z, &bb) < 1e-13);
}

#[test]
fn band_cholesky_rejects_indefinite() {
    let mut a = laplacian(3);
    a.values[a.row_ptr[4] + 2] = -1.0;
    assert!(matches!(BandCholesky::<f64>::new().factor(&a), Err(LinalgError::Singular(_))));
}

#[test]
fn band_matrix_fortran_ops() {
    use elmer_direct::{BandMatrix, SymBandMatrix};
    // BAND_INDEX(i,j) with Subband=2, 1-based: (j-1)*7 + i-j+5 -> 0-based j*7 + i-j+4.
    let mut m = BandMatrix::<f64>::new_subband(5, 2);
    assert_eq!(m.ldab(), 7);
    m.set(3, 1, 4.5);
    assert_eq!(m.data()[7 + 3 - 1 + 4], 4.5);
    m.add(3, 1, 0.5);
    assert_eq!(m.get(3, 1), 5.0);
    m.zero_row(3);
    assert_eq!(m.get(3, 1), 0.0);

    // glue_local with dofs=2 equals scalar adds.
    let mut g = BandMatrix::<f64>::new_subband(6, 5);
    let local: Vec<f64> = (0..16).map(|v| v as f64).collect();
    g.glue_local(2, &[0, 2], &local);
    assert_eq!(g.get(0, 0), local[0]);
    assert_eq!(g.get(1, 4), local[4 + 2]);
    assert_eq!(g.get(5, 0), local[3 * 4]);

    // SBand Dirichlet keeps symmetry and matches the dense computation.
    let mut s = SymBandMatrix::<f64>::new(4, 1);
    for i in 0..4 {
        s.set(i, i, 2.0);
        if i > 0 {
            s.set(i, i - 1, -1.0);
        }
    }
    let mut b = vec![1.0; 4];
    s.set_dirichlet(&mut b, 1, 3.0);
    assert_eq!(b, vec![4.0, 3.0, 4.0, 1.0]);
    assert_eq!(s.get(1, 0), 0.0);
    assert_eq!(s.get(1, 1), 1.0);
    let mut v = vec![0.0; 4];
    s.matvec(&[1.0, 2.0, 3.0, 4.0], &mut v);
    assert_eq!(v, vec![2.0, 2.0, 2.0 * 3.0 - 4.0, -3.0 + 8.0]);
}

#[test]
fn select_direct_keywords() {
    assert_eq!(kind_for_method("Banded").unwrap(), DirectKind::BandLu);
    assert_eq!(kind_for_method("symmetric banded").unwrap(), DirectKind::BandCholesky);
    for m in ["umfpack", "Big UMFPACK", "mumps", "superlu", "pardiso"] {
        assert_eq!(kind_for_method(m).unwrap(), DirectKind::SparseLu, "{m}");
    }
    assert_eq!(kind_for_method("cholmod").unwrap(), DirectKind::SparseCholesky);
    assert!(kind_for_method("bogus").is_err());

    let a = laplacian(5);
    let b = rhs(a.nrows, 3);
    for method in [None, Some("umfpack"), Some("cholmod"), Some("symmetric banded")] {
        let mut p = ParamMap::new();
        if let Some(m) = method {
            p.insert("Linear System Direct Method", Value::Str(m.into()));
        }
        let mut s = select_direct(&p);
        let x = solve_with(s.as_mut(), &a, &b);
        assert!(rel_residual(&a, &x, &b) < 1e-12, "{method:?}");
    }
    let mut p = ParamMap::new();
    p.insert("linear  system direct method", Value::Str("nope".into()));
    assert!(try_select_direct(&p).is_err());
}
