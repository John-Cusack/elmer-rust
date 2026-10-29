//! Eigenvalues and modes of the two pristine `ElmerSolver` runs under
//! `tests/fixtures` (see the README there).
//!
//! The configuration is the one Elmer used: `Eigen System Values = 10`,
//! tolerance 1e-6, smallest magnitude, consistent mass, no shift, and the
//! non-symmetric iteration (`Linear System Symmetric` unset, so `DNAUPD`).
//! Eigenvalues must match the tight reference run to 1e-10 relative, or to
//! the precision double arithmetic determines them to when that is coarser:
//! `eps * ||K||_inf * ||x||^2 / |lambda|` for a unit `M`-norm mode `x` (first
//! order perturbation of `K`). Only the softest geomstiff mode needs the
//! latter (about 6e-9). The first three modes are compared sign-normalised.
mod common;
use common::*;

use elmer_core::Logger;
use elmer_direct::{set_single_thread, FaerLu};
use elmer_eigen::{eigs, log_summary, EigenConfig, Mass};
use std::sync::Mutex;

fn run(name: &str, symmetric: bool) {
    set_single_thread();
    let f = fixture(name);
    let mut cfg = EigenConfig::new(10, symmetric);
    cfg.tol = 1e-6;
    let res = eigs(&f.k, Mass::Consistent(&f.m), &mut FaerLu::new(), &cfg).unwrap();
    assert_eq!(res.values.len(), 10);

    let k_inf = (0..f.k.nrows)
        .map(|i| {
            f.k.values[f.k.row_ptr[i]..f.k.row_ptr[i + 1]]
                .iter()
                .map(|v| v.abs())
                .sum::<f64>()
        })
        .fold(0.0, f64::max);
    let mut worst_tight = 0.0_f64;
    let mut worst_asis = 0.0_f64;
    let mut coarse = Vec::new();
    for (i, v) in res.values.iter().enumerate() {
        let (t, a) = (f.tight[i], f.asis[i]);
        let rt = (v.re - t).abs() / t.abs();
        worst_tight = worst_tight.max(rt);
        worst_asis = worst_asis.max((v.re - a).abs() / a.abs());
        let x2: f64 = res.vectors[i].iter().map(|z| z.norm_sqr()).sum();
        let determined = f64::EPSILON * k_inf * x2 / v.re.abs();
        let tol = determined.max(1e-10);
        if tol > 1e-10 {
            coarse.push(format!("mode {} diff {rt:.2e} <= {tol:.2e}", i + 1));
        }
        assert!(v.im.abs() <= 1e-10 * v.re.abs(), "{name} {i}: {v}");
        assert!(
            rt <= tol,
            "{name} mode {}: {:.13e} vs Fortran {t:.13e} ({rt:.2e} > {tol:.2e})",
            i + 1,
            v.re
        );
    }
    let mut worst_vec = 0.0_f64;
    for (i, x) in f.modes.iter().enumerate() {
        let d = sign_normalised_diff(&re(&res.vectors[i]), x);
        worst_vec = worst_vec.max(d);
        assert!(d < 1e-6, "{name} mode {}: vector differs by {d:.2e}", i + 1);
    }
    println!(
        "{name} (symmetric = {symmetric}): n = {}, {} update iterations, {} OP applications; \
         max rel eigenvalue diff vs tight {worst_tight:.2e}, vs asis {worst_asis:.2e}; \
         max sign-normalised mode diff {worst_vec:.2e}; coarser than 1e-10: {coarse:?}",
        f.k.nrows, res.iterations, res.op_applications
    );
}

#[test]
fn plates_eigen_matches_fortran() {
    run("plates_eigen", false);
}

#[test]
fn geomstiff_matches_fortran() {
    run("geomstiff", false);
}

#[test]
fn symmetric_iteration_matches_fortran_too() {
    run("plates_eigen", true);
    run("geomstiff", true);
}

struct Lines(Mutex<Vec<String>>);

impl Logger for Lines {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        // Max Output Level = 5 in the test; Elmer trims trailing blanks.
        if level <= 5 {
            self.0
                .lock()
                .unwrap()
                .push(format!("{caller}: {msg}").trim_end().to_string());
        }
    }
    fn warn(&self, _caller: &str, _msg: &str) {}
}

/// `ElmerSolver` stdout of the geomstiff run (asis.out.log).
const GEOMSTIFF_STDOUT: &str = "\
EigenSolve:
EigenSolve: Eigen system solution complete:
EigenSolve:
EigenSolve: Number of eigensystem iterations is: 33
EigenSolve: Number of converged Ritz values is: 10
EigenSolve: Number of update iterations taken: 1
EigenSolve: Computed 10 Eigen Values
EigenSolve: 1:    1.979368E+02   0.000000E+00
EigenSolve: 2:    4.853725E+03   0.000000E+00
EigenSolve: 3:    3.334773E+04   0.000000E+00
EigenSolve: 4:    1.058255E+05   0.000000E+00
EigenSolve: 5:    1.188349E+05   0.000000E+00
EigenSolve: 6:    3.023169E+05   0.000000E+00
EigenSolve: 7:    6.258318E+05   0.000000E+00
EigenSolve: 8:    9.517211E+05   0.000000E+00
EigenSolve: 9:    1.128225E+06   0.000000E+00
EigenSolve: 10:    1.842892E+06   0.000000E+00
EigenSolve: --------------------------------";

#[test]
fn geomstiff_summary_matches_elmer_stdout() {
    set_single_thread();
    let f = fixture("geomstiff");
    let mut cfg = EigenConfig::new(10, false);
    cfg.tol = 1e-6;
    let res = eigs(&f.k, Mass::Consistent(&f.m), &mut FaerLu::new(), &cfg).unwrap();
    let log = Lines(Mutex::new(Vec::new()));
    log_summary(&res, cfg.tol, &log);
    let got = log.0.into_inner().unwrap();
    let want: Vec<&str> = GEOMSTIFF_STDOUT.lines().collect();
    assert_eq!(got, want);
}

/// WaveEquEigen (unit square, structured mesh) has double eigenvalues. The
/// tight ARPACK run returns only one copy of 130.286 and of 173.064 and fills
/// up with 180.96 and 202.98. The test's own as-shipped run (iterative inner
/// solves) finds both copies but only to about 1e-9. Here all ten smallest
/// come back, each matching the tight run's value to 1e-10.
#[test]
fn wave_equ_eigen_keeps_multiplicities() {
    set_single_thread();
    let f = fixture("wave_equ_eigen");
    let t = &f.tight;
    let want = [t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[6], t[7], t[7]];
    for symmetric in [false, true] {
        let mut cfg = EigenConfig::new(10, symmetric);
        cfg.tol = 1e-6;
        let res = eigs(&f.k, Mass::Consistent(&f.m), &mut FaerLu::new(), &cfg).unwrap();
        for (i, (v, w)) in res.values.iter().zip(&want).enumerate() {
            let r = (v.re - w).abs() / w;
            assert!(
                r <= 1e-10,
                "symmetric {symmetric} mode {}: {} vs {w} ({r:.2e})",
                i + 1,
                v.re
            );
        }
        let n = f.k.nrows;
        let apply = |a: &elmer_core::CsrMatrix<f64>, x: &[f64]| {
            let mut y = vec![0.0; n];
            elmer_core::LinearOperator::apply(a, x, &mut y);
            y
        };
        for (v, x) in res.values.iter().zip(&res.vectors) {
            let x = re(x);
            let (kx, mx) = (apply(&f.k, &x), apply(&f.m, &x));
            let r: f64 = kx
                .iter()
                .zip(&mx)
                .map(|(a, b)| (a - v.re * b).powi(2))
                .sum::<f64>()
                .sqrt();
            let k: f64 = kx.iter().map(|a| a * a).sum::<f64>().sqrt();
            assert!(
                r / k < 1e-8,
                "symmetric {symmetric} lambda {}: residual {:.2e}",
                v.re,
                r / k
            );
        }
        // The two modes of a double eigenvalue span its eigenspace: M-orthogonal
        // from the symmetric iteration, merely independent from the general one
        // (as from DNEUPD).
        for (a, b) in [(1, 2), (4, 5), (6, 7), (8, 9)] {
            let (xa, xb) = (re(&res.vectors[a]), re(&res.vectors[b]));
            let dot = |p: &[f64], q: &[f64]| -> f64 {
                p.iter().zip(&apply(&f.m, q)).map(|(s, t)| s * t).sum()
            };
            let cos = dot(&xa, &xb) / (dot(&xa, &xa) * dot(&xb, &xb)).sqrt();
            let limit = if symmetric { 1e-8 } else { 0.5 };
            assert!(
                cos.abs() < limit,
                "symmetric {symmetric} modes {a}/{b}: M-cosine {cos:e}"
            );
        }
        let d = sign_normalised_diff(&re(&res.vectors[0]), &f.modes[0]);
        assert!(d < 1e-8, "mode 1 vector differs by {d:.2e}");
    }
}
