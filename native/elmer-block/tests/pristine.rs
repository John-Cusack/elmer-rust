//! Block solves of monolithic systems the pristine ElmerSolver saved
//! (`Linear System Save`, full-precision list-directed output) for two quick
//! tests, against the iteration counts and norms pristine printed.
//!
//! - `blocklinelast2b`: BlockLinElast2b (StressSolve, 2 dofs, `Additive
//!   Namespaces`): outer GCR over block Jacobi with BiCGStab and ILU0
//!   blocks, no scaling.
//! - `step_stokes`: Step_stokes_vec_blockprec `stabilized.sif`
//!   (IncompressibleNS, 3 dofs): outer GCR over block Gauss-Seidel with
//!   IDR(s) and ILU0 blocks, the pressure block preconditioned by the Schur
//!   matrix, inside the monolithic diagonal scaling.

use std::path::Path;

use elmer_block::{block_solve, comp_norm, BlockOptions, OuterScaling};
use elmer_core::{CsrMatrix, NullLogger, ParamMap, SerialComm, Value};
use elmer_sparse::{read_linear_system_save, Columns};

fn load(case: &str) -> (CsrMatrix<f64>, Vec<f64>) {
    let s = read_linear_system_save(&fixture(case), "linsys", Columns::default()).unwrap();
    (s.stiff.unwrap(), s.rhs.unwrap())
}

fn fixture(case: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(case)
}

fn section(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn s(v: &str) -> Value {
    Value::Str(v.into())
}

fn i(v: i32) -> Value {
    Value::Integer(vec![v])
}

fn l(v: bool) -> Value {
    Value::Logical(v)
}

fn r(v: f64) -> Value {
    Value::real(v)
}

#[test]
fn blocklinelast2b_matches_pristine() {
    let (a, b) = load("blocklinelast2b");
    let p = section(&[
        ("Linear System Solver", s("block")),
        ("Linear System Scaling", l(false)),
        ("Block: Linear System Solver", s("iterative")),
        ("Block: Linear System Refactorize", l(false)),
        ("Block: Linear System Convergence Tolerance", r(0.1e-2)),
        ("Block: Linear System Iterative Method", s("bicgstab")),
        ("Block: Linear System Max Iterations", i(500)),
        ("Block: Linear System Symmetric", l(true)),
        ("Block: Linear System Residual Output", i(20)),
        ("Block: Linear System Preconditioning", s("ILU0")),
        ("Linear System Block Method", s("other")),
        ("Outer: Linear System Solver", s("Iterative")),
        ("Outer: Linear System Convergence Tolerance", r(1e-7)),
        ("Outer: Linear System Iterative Method", s("GCR")),
        ("Outer: Linear System Symmetric", l(true)),
        ("Outer: Linear System GCR Restart", i(25)),
        ("Outer: Linear System Residual Output", i(1)),
        ("Outer: Linear System Max Iterations", i(500)),
    ]);
    let opts = BlockOptions {
        additive_namespaces: true,
        ..Default::default()
    };
    let mut x = vec![0.0; a.nrows];
    let rep = block_solve(&p, 2, &a, &b, &mut x, &opts, &NullLogger, &SerialComm).unwrap();

    // Pristine stdout: 21 GCR iterations; per sweep a BiCGStab solve of each
    // block, the very first one "Solution trivially zero!".
    let outer = rep.outer.unwrap();
    assert_eq!(rep.nblocks, 2);
    assert_eq!(outer.iterations, 21);
    assert_eq!(rep.sweeps, 21);
    let inner: Vec<usize> = rep.inner.iter().map(|r| r.iterations).collect();
    #[rustfmt::skip]
    let pristine = [
        0, 7, 6, 7, 7, 7, 6, 7, 8, 8, 5, 8, 5, 7, 5, 8, 9, 8, 5, 11, 5,
        5, 5, 6, 8, 10, 6, 11, 9, 11, 8, 13, 6, 6, 9, 6, 10, 6, 10, 11, 10, 11,
    ];
    assert_eq!(inner, pristine);
    assert!((outer.residual - 0.7621e-07).abs() < 0.0001e-07, "{}", outer.residual);
    // ComputeChange NRM of Disp: 6.09028287E+00.
    let nrm = comp_norm(&x);
    assert!((nrm - 6.090_282_87).abs() < 1e-8, "{nrm}");
}

#[test]
fn step_stokes_matches_pristine() {
    let (mut a, b) = load("step_stokes");
    let mut p = section(&[
        ("Linear System Solver", s("Block")),
        ("Block Gauss-Seidel", l(true)),
        ("Block Matrix Reuse", l(false)),
        ("Block Scaling", l(false)),
        ("Block Preconditioner", l(true)),
        // IncompressibleNS adds this with Block Preconditioner.
        ("Block Matrix Schur Variable", s("schur")),
        ("Outer: Linear System Solver", s("Iterative")),
        ("Outer: Linear System Iterative Method", s("GCR")),
        ("Outer: Linear System GCR Restart", i(250)),
        ("Outer: Linear System Residual Output", i(1)),
        ("Outer: Linear System Max Iterations", i(500)),
        ("Outer: Linear System Abort Not Converged", l(false)),
        ("Outer: Linear System Convergence Tolerance", r(1e-6)),
    ]);
    for blk in [11, 22, 33] {
        for (k, v) in [
            ("Linear System Convergence Tolerance", r(0.01)),
            ("Linear System Solver", s("iterative")),
            ("Linear System Scaling", l(false)),
            ("Linear System Preconditioning", s("ilu")),
            ("Linear System Residual Output", i(100)),
            ("Linear System Max Iterations", i(1000)),
            ("Linear System Iterative Method", s("idrs")),
            ("Calculate Loads", l(false)),
        ] {
            p.insert(&format!("block {blk}: {k}"), v);
        }
    }
    // The initial guess holds the Dirichlet values (SolveLinearSystem's x):
    // rows with only a diagonal entry.
    let mut x: Vec<f64> = (0..a.nrows)
        .map(|row| {
            let (lo, hi) = (a.row_ptr[row], a.row_ptr[row + 1]);
            let only_diag = (lo..hi).all(|k| a.col_idx[k] as usize == row || a.values[k] == 0.0);
            let d = (lo..hi).find(|&k| a.col_idx[k] as usize == row).map(|k| a.values[k]);
            match d {
                Some(d) if only_diag => b[row] / d,
                _ => 0.0,
            }
        })
        .collect();
    // SolveLinearSystem's Linear System Scaling (default on) around the solve.
    let mut bs = b.clone();
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let scaling = elmer_sparse::scale_system(
        &mut a,
        Some(&mut bs[..]),
        Some(&mut x[..]),
        &elmer_sparse::ScaleOptions {
            rhs_scaling: bnorm != 0.0,
            ..Default::default()
        },
    );
    // The pressure Schur approximation pristine assembled (`prec 33: Linear
    // System Save`): -1/mu times the pressure mass matrix, unscaled.
    let schur = read_linear_system_save(&fixture("step_stokes"), "prec33", Columns::default())
        .unwrap()
        .stiff
        .unwrap();
    let opts = BlockOptions {
        outer_scaling: Some(OuterScaling {
            diag: scaling.diag(),
            rhs_scaling: scaling.rhs_scaling(),
        }),
        schur: Some(&schur),
        ..Default::default()
    };
    let rep = block_solve(&p, 3, &a, &bs, &mut x, &opts, &NullLogger, &SerialComm).unwrap();
    scaling.unscale(&mut a, Some(&mut bs[..]), Some(&mut x[..]));

    // BlockKrylovIter at output level 10: Rhs and Solution norm at start
    // 2.36691E-02, Solution norm 2.69841E-01, Relative change 1.67743E+00.
    assert!((rep.rhs_norm - 2.36691e-2).abs() < 1e-7, "{}", rep.rhs_norm);
    assert!((rep.start_norm - 2.36691e-2).abs() < 1e-7, "{}", rep.start_norm);
    assert!((rep.solution_norm - 2.69841e-1).abs() < 1e-6, "{}", rep.solution_norm);
    assert!((rep.max_change - 1.67743).abs() < 1e-5, "{}", rep.max_change);
    let outer = rep.outer.unwrap();
    assert_eq!(rep.nblocks, 3);
    assert_eq!(outer.iterations, 36);
    let inner: Vec<usize> = rep.inner.iter().map(|r| r.iterations).collect();
    #[rustfmt::skip]
    let pristine = [
        3, 3, 1, 4, 3, 1, 4, 3, 1, 3, 4, 1, 5, 2, 1, 3, 5, 1, 5, 3, 1, 4, 5, 1,
        5, 6, 1, 6, 5, 1, 6, 4, 1, 5, 4, 1, 7, 4, 1, 5, 3, 1, 6, 4, 1, 5, 6, 1,
        5, 4, 1, 7, 6, 1, 6, 4, 1, 6, 4, 1, 6, 6, 1, 6, 4, 1, 5, 4, 1, 7, 5, 1,
        5, 6, 1, 5, 4, 1, 6, 5, 1, 6, 5, 1, 5, 5, 1, 6, 7, 1, 6, 5, 1, 6, 7, 1,
        6, 5, 1, 6, 5, 1, 5, 4, 1, 5, 6, 1,
    ];
    // Pristine's IDR(s) draws its shadow space from gfortran's global
    // RANDOM_NUMBER stream; elmer-itermethods uses its own generator, so the
    // velocity blocks' counts vary by an iteration or two. The block order
    // and the Schur block (one iteration each) match exactly.
    assert_eq!(inner.len(), pristine.len());
    for (k, (got, want)) in inner.iter().zip(pristine).enumerate() {
        assert_eq!(rep.inner[k].block, k % 3);
        if k % 3 == 2 {
            assert_eq!(*got, want, "Schur block solve {k}");
        } else {
            assert!(got.abs_diff(want) <= 2, "solve {k}: {got} vs {want}");
        }
    }
    assert!((outer.residual - 0.5148e-06).abs() < 0.05e-06, "{}", outer.residual);
    // ComputeChange NRM of the flow solution, 7.55401834E-02: the RMS of the
    // velocity components (`Nonlinear System Norm Dofs = dim`).
    let vel: Vec<f64> = x
        .iter()
        .enumerate()
        .filter(|(k, _)| k % 3 != 2)
        .map(|(_, v)| *v)
        .collect();
    let nrm = comp_norm(&vel);
    assert!((nrm - 7.554_018_34e-2).abs() < 1e-9, "{nrm}");
}

/// BlockLinElast2b with `Block Preconditioner = False` (`BlockStandardIter`)
/// and the outer tolerance at 1e-5, run on pristine for this test.
fn standard_iteration(gauss_seidel: bool) -> (elmer_block::BlockReport, Vec<f64>) {
    let (a, b) = load("blocklinelast2b");
    let p = section(&[
        ("Linear System Solver", s("block")),
        ("Linear System Scaling", l(false)),
        ("Block Preconditioner", l(false)),
        ("Block Gauss-Seidel", l(gauss_seidel)),
        ("Block: Linear System Solver", s("iterative")),
        ("Block: Linear System Refactorize", l(false)),
        ("Block: Linear System Convergence Tolerance", r(0.1e-2)),
        ("Block: Linear System Iterative Method", s("bicgstab")),
        ("Block: Linear System Max Iterations", i(500)),
        ("Block: Linear System Preconditioning", s("ILU0")),
        ("Outer: Linear System Solver", s("Iterative")),
        ("Outer: Linear System Convergence Tolerance", r(1e-5)),
        ("Outer: Linear System Iterative Method", s("GCR")),
        ("Outer: Linear System Max Iterations", i(500)),
    ]);
    let opts = BlockOptions {
        additive_namespaces: true,
        ..Default::default()
    };
    let mut x = vec![0.0; a.nrows];
    let rep = block_solve(&p, 2, &a, &b, &mut x, &opts, &NullLogger, &SerialComm).unwrap();
    (rep, x)
}

#[test]
fn blocklinelast2b_gauss_seidel_iteration_matches_pristine() {
    let (rep, x) = standard_iteration(true);
    // "Converged after iterations: 28"; each iteration solves block 1 then 2.
    assert_eq!(rep.outer.as_ref().unwrap().iterations, 28);
    let inner: Vec<usize> = rep.inner.iter().map(|r| r.iterations).collect();
    let mut pristine = vec![0, 7, 19, 7, 17, 7, 18, 7, 18, 7, 18, 7, 15, 7, 14, 7, 16, 7];
    pristine.extend([17, 7].repeat(19));
    assert_eq!(inner, pristine);
    // Last "Block 1 norms: 1.181E+00 8.231E-06", block 2 change 4.149E-06.
    assert!((rep.max_change - 8.231e-6).abs() < 0.001e-6, "{}", rep.max_change);
    // ComputeChange NRM 6.0902323.
    assert!((comp_norm(&x) - 6.090_232_3).abs() < 1e-7, "{}", comp_norm(&x));
}

#[test]
fn blocklinelast2b_jacobi_iteration_matches_pristine() {
    let (rep, x) = standard_iteration(false);
    assert_eq!(rep.outer.as_ref().unwrap().iterations, 51);
    let inner: Vec<usize> = rep.inner.iter().map(|r| r.iterations).collect();
    #[rustfmt::skip]
    let pristine = [
        0, 7, 19, 8, 16, 7, 17, 7, 20, 7, 18, 8, 19, 7, 18, 8, 18, 7, 18, 8, 19, 7, 14, 7, 19, 7,
        15, 8, 18, 7, 17, 8, 18, 7, 17, 8, 18, 7, 17, 8, 18, 7, 17, 8, 18, 7, 17, 8, 18, 7, 17, 8,
        18, 7, 17, 8, 18, 7, 17, 8, 19, 7, 17, 8, 19, 7, 17, 8, 19, 7, 17, 8, 19, 7, 17, 8, 20, 7,
        17, 8, 20, 7, 17, 9, 21, 7, 17, 9, 21, 7, 17, 9, 21, 7, 17, 9, 24, 7, 17, 9, 27, 7,
    ];
    assert_eq!(inner, pristine);
    // Last iteration: block 1 change 1.359E-06, block 2 8.298E-06.
    assert!((rep.max_change - 8.298e-6).abs() < 0.001e-6, "{}", rep.max_change);
    // ComputeChange NRM 6.0901689.
    assert!((comp_norm(&x) - 6.090_168_9).abs() < 1e-7, "{}", comp_norm(&x));
}
