//! Checks against values captured from the pristine ElmerSolver (see
//! `tests/fixtures/README.md`). The factors themselves are not dumped by
//! Elmer, so the preconditioners are pinned through one CG iteration from
//! zero, `x1 = alpha * M^-1 b`, printed to 8 digits in the VTU.

#![allow(clippy::needless_range_loop)]

use elmer_core::{CsrMatrix, Preconditioner};
use elmer_sparse::{
    read_linear_system_save, scale_system, spmv, write_linear_system_save, Columns, Ilu,
    SaveOptions, ScaleOptions,
};
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn load() -> (CsrMatrix<f64>, Vec<f64>, Vec<u32>) {
    let s = read_linear_system_save(&fixtures(), "sq", Columns::default()).unwrap();
    (s.stiff.unwrap(), s.rhs.unwrap(), s.perm.unwrap())
}

fn x1(name: &str) -> Vec<f64> {
    std::fs::read_to_string(fixtures().join(format!("x1_{name}.txt")))
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn check_cg1(a: &CsrMatrix<f64>, b: &[f64], perm: &[u32], m: &dyn Preconditioner, name: &str) {
    let n = a.nrows;
    let mut z = vec![0.0; n];
    m.apply(b, &mut z);
    let mut az = vec![0.0; n];
    spmv(a, &z, &mut az);
    let alpha = dot(b, &z) / dot(&z, &az);
    let want = x1(name);
    let scale = want.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for (node, &w) in want.iter().enumerate() {
        let got = alpha * z[perm[node] as usize];
        assert!(
            (got - w).abs() <= 2e-7 * scale,
            "{name}: node {node}: got {got:e}, Fortran {w:e}"
        );
    }
}

#[test]
fn saved_system_shape() {
    let (a, b, perm) = load();
    assert_eq!((a.nrows, a.nnz(), b.len(), perm.len()), (49, 289, 49, 49));
}

#[test]
fn cg_first_iteration_matches_fortran() {
    let (a, b, perm) = load();
    check_cg1(&a, &b, &perm, &Ilu::<f64>::diagonal(&a), "diag");
    for (level, name) in [(0, "ilu0"), (1, "ilu1"), (2, "ilu2")] {
        let ilu = Ilu::ilu_n(&a, level).unwrap();
        check_cg1(&a, &b, &perm, &ilu, name);
    }
}

/// Elmer scales the system (`ScaleLinearSystem`, diagonal) before factorizing
/// and iterating, so the ILUT drop decisions see the scaled rows.
#[test]
fn ilut_on_the_scaled_system_matches_fortran() {
    let (a0, b0, perm) = load();
    let (mut a, mut b) = (a0, b0);
    let rs = scale_system(&mut a, Some(&mut b), None, &ScaleOptions::default());
    let ilut = Ilu::ilut(&a, 0.05).unwrap();
    // "ILU(T) (Real), NOF nonzeros: 237" in the Fortran log (ILURows(n+1), 1-based)
    assert_eq!(ilut.nnz(), 236);
    let n = a.nrows;
    let mut z = vec![0.0; n];
    ilut.apply(&b, &mut z);
    let mut az = vec![0.0; n];
    spmv(&a, &z, &mut az);
    let alpha = dot(&b, &z) / dot(&z, &az);
    let mut x: Vec<f64> = z.iter().map(|v| alpha * v).collect();
    rs.scale_vectors(&mut vec![0.0; n], Some(&mut x), true);
    let want = x1("ilut");
    let scale = want.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for (node, &w) in want.iter().enumerate() {
        let got = x[perm[node] as usize];
        assert!(
            (got - w).abs() <= 2e-7 * scale,
            "node {node}: got {got:e}, Fortran {w:e}"
        );
    }
}

#[test]
fn ilu_levels_differ_from_each_other() {
    let (a, _, _) = load();
    let n0 = Ilu::ilu_n(&a, 0).unwrap().nnz();
    let n1 = Ilu::ilu_n(&a, 1).unwrap().nnz();
    let n2 = Ilu::ilu_n(&a, 2).unwrap().nnz();
    assert_eq!(n0, a.nnz());
    assert!(n0 < n1 && n1 < n2, "{n0} {n1} {n2}");
}

#[test]
fn save_files_round_trip_byte_for_byte() {
    let s = read_linear_system_save(&fixtures(), "sq", Columns::default()).unwrap();
    let dir = std::env::temp_dir().join(format!("elmer-sparse-rt-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    write_linear_system_save(&dir, "sq", &s, &SaveOptions::default()).unwrap();
    for f in ["sq_a.dat", "sq_b.dat", "sq_perm.dat", "sq_sizes.dat"] {
        let want = std::fs::read_to_string(fixtures().join(f)).unwrap();
        let got = std::fs::read_to_string(dir.join(f)).unwrap();
        assert_eq!(got, want, "{f}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
