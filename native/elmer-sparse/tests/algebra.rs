//! Self-contained checks of the numerical kernels against dense arithmetic.

#![allow(clippy::needless_range_loop)]

use elmer_core::{CsrMatrix, LinearOperator, Preconditioner};
use elmer_sparse::scalar::ElmerScalar;
use elmer_sparse::*;
use num_complex::Complex64;

/// 2D five-point operator on an `m x m` grid with convection (nonsymmetric).
fn grid(m: usize, conv: f64) -> CsrMatrix<f64> {
    let mut t = Vec::new();
    for y in 0..m {
        for x in 0..m {
            let i = y * m + x;
            t.push((i, i as u32, 4.0 + 0.1 * (i % 3) as f64));
            if x > 0 {
                t.push((i, (i - 1) as u32, -1.0 - conv));
            }
            if x + 1 < m {
                t.push((i, (i + 1) as u32, -1.0 + conv));
            }
            if y > 0 {
                t.push((i, (i - m) as u32, -1.0));
            }
            if y + 1 < m {
                t.push((i, (i + m) as u32, -1.0));
            }
        }
    }
    from_triplets(m * m, m * m, t)
}

fn tridiag(n: usize) -> CsrMatrix<f64> {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i as u32, 2.0 + 0.3 * i as f64));
        if i > 0 {
            t.push((i, (i - 1) as u32, -1.0));
        }
        if i + 1 < n {
            t.push((i, (i + 1) as u32, -0.7));
        }
    }
    from_triplets(n, n, t)
}

fn vec_n(n: usize) -> Vec<f64> {
    (0..n).map(|i| ((i * 7 + 3) % 11) as f64 - 4.5).collect()
}

fn residual<T: ElmerScalar>(a: &CsrMatrix<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::zero(); a.nrows];
    spmv(a, x, &mut ax);
    ax.iter()
        .zip(b)
        .map(|(p, q)| (*p - *q).abs())
        .fold(0.0, f64::max)
}

#[test]
fn ilu0_of_tridiagonal_is_exact() {
    let a = tridiag(20);
    let b = vec_n(20);
    let mut x = b.clone();
    Ilu::ilu_n(&a, 0).unwrap().solve_in_place(&mut x);
    assert!(residual(&a, &x, &b) < 1e-12);
}

#[test]
fn ilut_with_zero_tolerance_is_complete_lu() {
    let a = grid(5, 0.3);
    let b = vec_n(25);
    let mut x = b.clone();
    Ilu::ilut(&a, 0.0).unwrap().solve_in_place(&mut x);
    assert!(residual(&a, &x, &b) < 1e-11);
}

#[test]
fn ilu_high_level_on_small_grid_is_complete_lu() {
    let a = grid(4, 0.2);
    let b = vec_n(16);
    let mut x = b.clone();
    Ilu::ilu_n(&a, 9).unwrap().solve_in_place(&mut x);
    assert!(residual(&a, &x, &b) < 1e-11);
}

fn arrow(n: usize) -> CsrMatrix<f64> {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i as u32, 5.0 + i as f64));
        if i > 0 {
            t.push((0, i as u32, 1.0));
            t.push((i, 0, 2.0));
        }
    }
    from_triplets(n, n, t)
}

#[test]
fn known_fill_pattern() {
    let a = arrow(4);
    let ilu0 = Ilu::ilu_n(&a, 0).unwrap();
    assert_eq!(ilu0.row_ptr(), &a.row_ptr[..]);
    assert_eq!(ilu0.col_idx(), &a.col_idx[..]);
    let ilu1 = Ilu::ilu_n(&a, 1).unwrap();
    // eliminating column 0 fills every row of the arrow: the factor is dense
    assert_eq!(ilu1.row_ptr(), &[0, 4, 8, 12, 16]);
    assert_eq!(ilu1.diag_idx(), &[0, 5, 10, 15]);
    let b = vec_n(4);
    let mut x = b.clone();
    ilu1.solve_in_place(&mut x);
    assert!(residual(&a, &x, &b) < 1e-12);
    // level 1 is a no-op beyond that
    assert_eq!(Ilu::ilu_n(&a, 2).unwrap().nnz(), 16);
}

#[test]
fn ilut_drops_small_entries() {
    let a = grid(6, 0.0);
    let loose = Ilu::ilut(&a, 0.0).unwrap().nnz();
    let tight = Ilu::ilut(&a, 0.2).unwrap().nnz();
    assert!(tight < loose);
}

#[test]
fn missing_diagonal_is_an_error() {
    let a = from_triplets(2, 2, [(0, 0u32, 1.0), (1, 0, 1.0)]);
    assert!(Ilu::ilu_n(&a, 0).is_err());
    assert!(Ilu::ilut(&a, 0.1).is_err());
}

#[test]
fn transpose_solve_is_adjoint() {
    let a = grid(5, 0.4);
    for ilu in [
        Ilu::ilu_n(&a, 0).unwrap(),
        Ilu::ilu_n(&a, 1).unwrap(),
        Ilu::ilut(&a, 0.05).unwrap(),
    ] {
        let (x, y) = (
            vec_n(25),
            (0..25).map(|i| (i as f64).sin()).collect::<Vec<_>>(),
        );
        let (mut mx, mut mty) = (vec![0.0; 25], vec![0.0; 25]);
        ilu.apply(&x, &mut mx);
        ilu.apply_transpose(&y, &mut mty).unwrap();
        let l: f64 = mx.iter().zip(&y).map(|(p, q)| p * q).sum();
        let r: f64 = x.iter().zip(&mty).map(|(p, q)| p * q).sum();
        assert!((l - r).abs() < 1e-12 * l.abs().max(1.0));
    }
}

#[test]
fn incomplete_cholesky_exact_on_tridiagonal() {
    let n = 12;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i as u32, 3.0));
        if i > 0 {
            t.push((i, (i - 1) as u32, -1.0));
        }
        if i + 1 < n {
            t.push((i, (i + 1) as u32, -1.0));
        }
    }
    let a = from_triplets(n, n, t);
    let b = vec_n(n);
    let mut x = b.clone();
    Ilu::ic_n(&a, 0, None).unwrap().solve_in_place(&mut x);
    assert!(residual(&a, &x, &b) < 1e-12);
}

#[test]
fn diagonal_preconditioner_guards_tiny_pivots() {
    let a = from_triplets(2, 2, [(0, 0u32, 2.0), (1, 1, 1e-300)]);
    let d = Ilu::<f64>::diagonal(&a);
    let mut z = [0.0; 2];
    d.apply(&[4.0, 3.0], &mut z);
    assert_eq!(z, [2.0, 3.0]);
}

fn cgrid(m: usize) -> CsrMatrix<Complex64> {
    let g = grid(m, 0.25);
    let values = g
        .values
        .iter()
        .enumerate()
        .map(|(k, &v)| {
            Complex64::new(
                v,
                if g.col_idx[k] as usize == 0 {
                    0.3
                } else {
                    0.1 * v
                },
            )
        })
        .collect();
    CsrMatrix {
        nrows: g.nrows,
        ncols: g.ncols,
        row_ptr: g.row_ptr,
        col_idx: g.col_idx,
        values,
        diag: g.diag,
    }
}

#[test]
fn complex_factorizations() {
    let a = cgrid(4);
    let b: Vec<Complex64> = (0..16)
        .map(|i| Complex64::new(i as f64, 1.0 - i as f64))
        .collect();
    for ilu in [Ilu::ilu_n(&a, 9).unwrap(), Ilu::ilut(&a, 0.0).unwrap()] {
        let mut x = b.clone();
        ilu.solve_in_place(&mut x);
        assert!(residual(&a, &x, &b) < 1e-11);
    }
    let ilu0 = Ilu::ilu_n(&a, 0).unwrap();
    let (x, y) = (b.clone(), b.iter().rev().cloned().collect::<Vec<_>>());
    let (mut mx, mut mty) = (
        vec![Complex64::default(); 16],
        vec![Complex64::default(); 16],
    );
    ilu0.apply(&x, &mut mx);
    ilu0.apply_transpose(&y, &mut mty).unwrap();
    let l: Complex64 = mx.iter().zip(&y).map(|(p, q)| *p * *q).sum();
    let r: Complex64 = x.iter().zip(&mty).map(|(p, q)| *p * *q).sum();
    assert!((l - r).norm() < 1e-11 * l.norm());
}

#[test]
fn complex_division_matches_reference_values() {
    let q = Complex64::new(1.0, 2.0).elmer_div(Complex64::new(3.0, -4.0));
    assert!((q.re - -0.2).abs() < 1e-16 && (q.im - 0.4).abs() < 1e-16);
    let q = Complex64::new(1.0, 0.0).elmer_div(Complex64::new(0.0, 2.0));
    assert!(q.re.abs() < 1e-16 && (q.im + 0.5).abs() < 1e-16);
}

fn roundtrip<T: ElmerScalar + std::fmt::Debug>(
    a0: &CsrMatrix<T>,
    b0: &[T],
    x0: &[T],
    method: ScalingMethod,
) {
    let (mut a, mut b, mut x) = (a0.clone(), b0.to_vec(), x0.to_vec());
    let opts = ScaleOptions {
        method,
        ..Default::default()
    };
    let rs = scale_system(&mut a, Some(&mut b), Some(&mut x), &opts);
    // vectors transform like the system does
    let (mut bv, mut xv) = (b0.to_vec(), x0.to_vec());
    rs.scale_vectors(&mut bv, Some(&mut xv), false);
    for i in 0..b.len() {
        assert!((b[i] - bv[i]).abs() < 1e-12, "{method:?} b[{i}]");
        assert!((x[i] - xv[i]).abs() < 1e-12, "{method:?} x[{i}]");
    }
    rs.unscale(&mut a, Some(&mut b), Some(&mut x));
    for (p, q) in a.values.iter().zip(&a0.values) {
        assert!((*p - *q).abs() <= 1e-13 * q.abs().max(1.0), "{method:?}");
    }
    for i in 0..b.len() {
        assert!((b[i] - b0[i]).abs() <= 1e-12 * b0[i].abs().max(1.0));
        assert!((x[i] - x0[i]).abs() <= 1e-12 * x0[i].abs().max(1.0));
    }
}

#[test]
fn scaling_round_trips() {
    let a = grid(4, 0.3);
    let (b, x) = (
        vec_n(16),
        vec_n(16).iter().map(|v| v + 0.5).collect::<Vec<_>>(),
    );
    for m in [
        ScalingMethod::Diagonal,
        ScalingMethod::RowEquilibration,
        ScalingMethod::Constant,
        ScalingMethod::None,
    ] {
        roundtrip(&a, &b, &x, m);
    }
    let ca = cgrid(4);
    let cb: Vec<Complex64> = b.iter().map(|&v| Complex64::new(v, 1.0)).collect();
    let cx: Vec<Complex64> = x.iter().map(|&v| Complex64::new(1.0, v)).collect();
    for m in [
        ScalingMethod::Diagonal,
        ScalingMethod::RowEquilibration,
        ScalingMethod::Constant,
    ] {
        roundtrip(&ca, &cb, &cx, m);
    }
}

#[test]
fn scaling_properties() {
    let mut a = grid(4, 0.3);
    scale_system::<f64>(&mut a, None, None, &ScaleOptions::default());
    for i in 0..a.nrows {
        assert!((get_element(&a, i, i as u32) - 1.0).abs() < 1e-14);
    }
    let mut a = grid(4, 0.3);
    scale_system::<f64>(
        &mut a,
        None,
        None,
        &ScaleOptions {
            method: ScalingMethod::RowEquilibration,
            ..Default::default()
        },
    );
    for i in 0..a.nrows {
        assert!((row_sum_abs(&a, i) - 1.0).abs() < 1e-14);
    }
    assert_eq!(
        ScalingMethod::from_keyword("RowSum").unwrap(),
        ScalingMethod::RowEquilibration
    );
    assert!(ScalingMethod::from_keyword("bogus").is_err());
}

#[test]
fn scaling_preserves_the_solution() {
    let a0 = grid(4, 0.3);
    let b0 = vec_n(16);
    let (mut a, mut b) = (a0.clone(), b0.clone());
    let rs = scale_system(&mut a, Some(&mut b), None, &ScaleOptions::default());
    let mut y = b.clone();
    Ilu::ilut(&a, 0.0).unwrap().solve_in_place(&mut y);
    rs.unscale(&mut a, None, Some(&mut y));
    assert!(residual(&a0, &y, &b0) < 1e-11);
}

#[test]
fn spmv_variants_agree() {
    let a = grid(5, 0.3);
    let x = vec_n(25);
    let (mut y1, mut y2) = (vec![0.0; 25], vec![0.0; 25]);
    spmv(&a, &x, &mut y1);
    CsrOp(&a).apply(&x, &mut y2);
    let mut y3 = vec![0.0; 25];
    elmer_core::LinearOperator::apply(&a, &x, &mut y3);
    for i in 0..25 {
        assert!((y1[i] - y2[i]).abs() == 0.0 && (y1[i] - y3[i]).abs() < 1e-13);
    }
    let mut yt = vec![0.0; 25];
    spmv_transpose(&a, &x, &mut yt);
    let mut yt2 = vec![0.0; 25];
    elmer_core::LinearOperator::apply_transpose(&a, &x, &mut yt2).unwrap();
    assert_eq!(yt, yt2);
}

#[test]
fn spmv_blocked_matches_plain() {
    // 3 dofs per node, dense 3x3 blocks on a path of 4 nodes
    let nodes: usize = 4;
    let mut t = Vec::new();
    for nd in 0..nodes {
        for nb in [nd.wrapping_sub(1), nd, nd + 1] {
            if nb >= nodes {
                continue;
            }
            for r in 0..3 {
                for c in 0..3 {
                    t.push((
                        3 * nd + r,
                        (3 * nb + c) as u32,
                        1.0 + (3 * nd + r) as f64 * 0.1 - c as f64 * 0.07,
                    ));
                }
            }
        }
    }
    let a = from_triplets(12, 12, t);
    let x = vec_n(12);
    let (mut y, mut yb) = (vec![0.0; 12], vec![0.0; 12]);
    spmv(&a, &x, &mut y);
    spmv_blocked(&a, 3, &x, &mut yb).unwrap();
    for i in 0..12 {
        assert!((y[i] - yb[i]).abs() < 1e-13);
    }
    // a row with 4 entries is not a multiple of 3 in the end rows
    assert!(spmv_blocked(&tridiag(5), 2, &vec_n(5), &mut [0.0; 5]).is_err());
}

#[test]
fn triplets_sort_and_element_access() {
    let a = from_triplets(
        3,
        3,
        [
            (1, 2u32, 1.0),
            (1, 0, 2.0),
            (1, 2, 3.0),
            (0, 0, 5.0),
            (2, 2, 7.0),
        ],
    );
    assert_eq!(a.row_ptr, vec![0, 1, 3, 4]);
    assert_eq!(a.col_idx, vec![0, 0, 2, 2]);
    assert_eq!(a.values, vec![5.0, 2.0, 4.0, 7.0]);
    assert_eq!(a.diag.as_ref().unwrap(), &vec![0, NO_DIAG, 3]);
    let mut b = a.clone();
    assert!(add_to_element(&mut b, 1, 2, 1.0));
    assert!(!add_to_element(&mut b, 1, 1, 1.0));
    assert!(add_to_element(&mut b, 1, 1, 0.0));
    assert!(set_element(&mut b, 0, 0, 9.0));
    assert!(contains(&b, 1, 0) && !contains(&b, 0, 2));
    assert_eq!(get_element(&b, 1, 2), 5.0);
    assert_eq!(row_sum(&b, 1), 7.0);
    zero_row(&mut b, 1);
    assert_eq!(row_sum_abs(&b, 1), 0.0);
    zero_matrix(&mut b);
    assert!(b.values.iter().all(|&v| v == 0.0));
    // unsorted input matrix gets sorted
    let mut u = CsrMatrix {
        nrows: 2,
        ncols: 2,
        row_ptr: vec![0, 2, 3],
        col_idx: vec![1, 0, 1],
        values: vec![1.0, 2.0, 3.0],
        diag: Some(vec![0, 2]),
    };
    sort_rows(&mut u);
    assert_eq!(u.col_idx, vec![0, 1, 1]);
    assert_eq!(u.values, vec![2.0, 1.0, 3.0]);
    assert_eq!(u.diag.unwrap(), vec![0, 2]);
}

#[test]
fn glue_local_matches_dense_scatter() {
    // two triangles of a unit square, 2 dofs per node
    let elems: [[u32; 3]; 2] = [[0, 1, 2], [1, 3, 2]];
    let dofs = 2;
    let mut t = Vec::new();
    for e in &elems {
        for &i in e {
            for &j in e {
                for r in 0..dofs {
                    for c in 0..dofs {
                        t.push((dofs * i as usize + r, (dofs * j as usize + c) as u32, 0.0));
                    }
                }
            }
        }
    }
    let mut a = from_triplets(8, 8, t);
    let mut dense = vec![0.0; 64];
    for (k, e) in elems.iter().enumerate() {
        let local: Vec<f64> = (0..36)
            .map(|q| (q + 1) as f64 * (k + 1) as f64 * 0.5)
            .collect();
        glue_local(&mut a, dofs, e, &local);
        for (li, &i) in e.iter().enumerate() {
            for (lj, &j) in e.iter().enumerate() {
                for r in 0..dofs {
                    for c in 0..dofs {
                        dense[(dofs * i as usize + r) * 8 + dofs * j as usize + c] +=
                            local[(dofs * li + r) * 6 + dofs * lj + c];
                    }
                }
            }
        }
    }
    for i in 0..8 {
        for j in 0..8 {
            assert_eq!(get_element(&a, i, j as u32), dense[i * 8 + j]);
        }
    }
    // skipped nodes and absent entries are ignored
    let before = a.values.clone();
    glue_local(&mut a, 1, &[NO_INDEX, 0], &[1.0, 1.0, 1.0, 1.0]);
    assert_eq!(a.values[0], before[0] + 1.0);
}

#[test]
fn list_directed_formatting_matches_gfortran() {
    use elmer_sparse::save::list_f64;
    // captured from gfortran `WRITE(*,*)`
    let cases = [
        (0.0, "   0.0000000000000000     "),
        (1.0, "   1.0000000000000000     "),
        (12.5, "   12.500000000000000     "),
        (123456.789, "   123456.78900000000     "),
        (0.1, "  0.10000000000000001     "),
        (0.09999999999999999, "   9.9999999999999992E-002"),
        (1e15, "   1000000000000000.0     "),
        (1e16, "   10000000000000000.     "),
        (1e17, "   1.0000000000000000E+017"),
        (-2.5e-5, "  -2.5000000000000001E-005"),
        (9.9999999999999999, "   10.000000000000000     "),
        (-1e100, "  -1.0000000000000000E+100"),
    ];
    for (v, s) in cases {
        assert_eq!(list_f64(v), s, "{v:e}");
    }
}

#[test]
fn save_with_damp_mass_and_sum_round_trips() {
    let a = grid(3, 0.2);
    let damp: Vec<f64> = (0..a.nnz()).map(|k| 0.5 * k as f64).collect();
    let mass: Vec<f64> = (0..a.nnz())
        .map(|k| if k % 4 == 0 { 0.0 } else { 1e-7 * k as f64 })
        .collect();
    let sys = SavedSystem {
        stiff: Some(a.clone()),
        damp: Some(damp.clone()),
        mass: Some(mass.clone()),
        rhs: Some(vec_n(9)),
        perm: Some(
            (0..9)
                .map(|i| if i == 4 { elmer_core::NO_DOF } else { i })
                .collect(),
        ),
        solution: None,
    };
    let dir = std::env::temp_dir().join(format!("elmer-sparse-t2-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let opts = SaveOptions {
        save_damp: true,
        save_mass: true,
        offset_ind: Some(7),
        ..Default::default()
    };
    write_linear_system_save(&dir, "s", &sys, &opts).unwrap();
    let cols = Columns {
        stiff: true,
        damp: true,
        mass: true,
    };
    let back = read_linear_system_save(&dir, "s", cols).unwrap();
    let m = back.stiff.unwrap();
    assert_eq!(m.row_ptr, a.row_ptr);
    assert_eq!(m.col_idx, a.col_idx);
    assert_eq!(m.values, a.values);
    assert_eq!(back.damp.unwrap(), damp);
    assert_eq!(back.mass.unwrap(), mass);
    assert_eq!(back.rhs.unwrap(), sys.rhs.clone().unwrap());
    assert_eq!(back.perm, sys.perm);
    let sizes = std::fs::read_to_string(dir.join("s_sizes.dat")).unwrap();
    assert_eq!(sizes.lines().count(), 6);
    // skip zeros + sum
    let opts = SaveOptions {
        skip_zeros: true,
        save_sum: true,
        save_mass: true,
        no_stiff: true,
        ..Default::default()
    };
    write_linear_system_save(&dir, "z", &sys, &opts).unwrap();
    let az = std::fs::read_to_string(dir.join("z_a.dat")).unwrap();
    assert_eq!(
        az.lines().count(),
        a.nnz() - (0..a.nnz()).filter(|k| k % 4 == 0).count()
    );
    let bz = std::fs::read_to_string(dir.join("z_b.dat")).unwrap();
    assert_eq!(bz.lines().next().unwrap().split_whitespace().count(), 4);
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(
        elmer_sparse::save::parse_fortran_f64("1.5D+02"),
        Some(150.0)
    );
    assert_eq!(
        elmer_sparse::save::parse_fortran_f64("1.5-100"),
        Some(1.5e-100)
    );
}
