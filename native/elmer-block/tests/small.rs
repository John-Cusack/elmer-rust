//! Block solves of small coupled systems: correctness of the solution
//! (residual of the monolithic system), the Jacobi/Gauss-Seidel difference,
//! `Block Structure`, `Block Order`, the block Jacobi/Gauss-Seidel outer
//! iteration, a system assembled per block, and the refusals.

use elmer_block::{block_solve, block_solve_system, BlockError, BlockOptions, BlockSystem, ElmerInner};
use elmer_core::{CsrMatrix, LinearOperator, NullLogger, ParamMap, SerialComm, Value};

/// A 1D chain of `nodes` nodes with `dofs` components: a Laplacian per
/// component plus a weaker coupling between components at each node.
fn coupled(nodes: usize, dofs: usize) -> CsrMatrix<f64> {
    let n = nodes * dofs;
    let mut a = CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::new();
    for i in 0..n {
        let (node, comp) = (i / dofs, i % dofs);
        for j in 0..n {
            let (nj, cj) = (j / dofs, j % dofs);
            let v = if nj == node && cj == comp {
                4.0 + comp as f64
            } else if nj == node {
                0.3
            } else if nj.abs_diff(node) == 1 && cj == comp {
                -1.0
            } else {
                continue;
            };
            if i == j {
                diag.push(a.values.len());
            }
            a.col_idx.push(j as u32);
            a.values.push(v);
        }
        a.row_ptr.push(a.values.len());
    }
    a.diag = Some(diag);
    a
}

fn params(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn base(inner: &str) -> ParamMap {
    let mut p = params(&[
        ("Linear System Solver", Value::Str("block".into())),
        ("Outer: Linear System Iterative Method", Value::Str("gcr".into())),
        ("Outer: Linear System Max Iterations", Value::Integer(vec![200])),
        ("Outer: Linear System Convergence Tolerance", Value::real(1e-10)),
        ("Linear System Scaling", Value::Logical(false)),
    ]);
    for blk in [11, 22, 33] {
        p.insert(&format!("block {blk}: Linear System Solver"), Value::Str(inner.into()));
        p.insert(
            &format!("block {blk}: Linear System Direct Method"),
            Value::Str("umfpack".into()),
        );
        p.insert(
            &format!("block {blk}: Linear System Iterative Method"),
            Value::Str("bicgstab".into()),
        );
        p.insert(
            &format!("block {blk}: Linear System Preconditioning"),
            Value::Str("ilu0".into()),
        );
        p.insert(
            &format!("block {blk}: Linear System Max Iterations"),
            Value::Integer(vec![500]),
        );
        p.insert(
            &format!("block {blk}: Linear System Convergence Tolerance"),
            Value::real(1e-12),
        );
    }
    p
}

fn rhs(n: usize) -> Vec<f64> {
    (0..n).map(|i| ((i * 7 % 5) as f64) - 1.5).collect()
}

fn residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> f64 {
    let mut r = vec![0.0; a.nrows];
    a.apply(x, &mut r);
    let rn: f64 = r.iter().zip(b).map(|(r, b)| (r - b).powi(2)).sum::<f64>().sqrt();
    rn / b.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn direct_blocks_solve_the_monolithic_system() {
    let a = coupled(12, 2);
    let b = rhs(a.nrows);
    for gs in [false, true] {
        let mut p = base("direct");
        p.insert("Block Gauss-Seidel", Value::Logical(gs));
        let mut x = vec![0.0; a.nrows];
        let rep = block_solve(
            &p,
            2,
            &a,
            &b,
            &mut x,
            &BlockOptions::default(),
            &NullLogger,
            &SerialComm,
        )
        .unwrap();
        assert!(residual(&a, &x, &b) < 1e-9, "gs={gs}: {}", residual(&a, &x, &b));
        assert_eq!(rep.nblocks, 2);
        assert_eq!(rep.inner.len(), 2 * rep.sweeps);
        assert!(rep.inner.iter().all(|r| r.iterations == 0));
    }
}

#[test]
fn gauss_seidel_needs_fewer_sweeps_than_jacobi() {
    let a = coupled(20, 3);
    let b = rhs(a.nrows);
    let sweeps = |gs: bool| {
        let mut p = base("iterative");
        p.insert("Block Gauss-Seidel", Value::Logical(gs));
        let mut x = vec![0.0; a.nrows];
        let rep = block_solve(
            &p,
            3,
            &a,
            &b,
            &mut x,
            &BlockOptions::default(),
            &NullLogger,
            &SerialComm,
        )
        .unwrap();
        assert!(residual(&a, &x, &b) < 1e-8);
        rep.sweeps
    };
    let (jacobi, gs) = (sweeps(false), sweeps(true));
    assert!(gs < jacobi, "Gauss-Seidel {gs} vs Jacobi {jacobi}");
}

#[test]
fn block_structure_and_order() {
    let a = coupled(10, 3);
    let b = rhs(a.nrows);
    let mut p = base("direct");
    p.insert("Block Structure", Value::Integer(vec![1, 1, 2]));
    p.insert("Block Gauss-Seidel", Value::Logical(true));
    p.insert("Block Order", Value::Integer(vec![2, 1]));
    let mut x = vec![0.0; a.nrows];
    let rep = block_solve(
        &p,
        3,
        &a,
        &b,
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    assert_eq!(rep.nblocks, 2);
    assert_eq!(rep.inner[0].block, 1);
    assert_eq!(rep.inner[1].block, 0);
    assert!(residual(&a, &x, &b) < 1e-9);
}

#[test]
fn dummy_block_passes_its_rhs_through() {
    // `block:` is the namespace of the sweep only: block 2 is preconditioned
    // by the identity while the product keeps it, and GCR still converges.
    let a = coupled(8, 2);
    let b = rhs(a.nrows);
    let mut p = base("direct");
    p.insert("Block: Dummy Block 2", Value::Logical(true));
    let mut x = vec![0.0; a.nrows];
    let rep = block_solve(
        &p,
        2,
        &a,
        &b,
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    assert!(rep.inner.iter().all(|r| r.block == 0));
    assert!(residual(&a, &x, &b) < 1e-9);
}

#[test]
fn zero_rhs_skips_the_solver() {
    let a = coupled(4, 2);
    let mut x = vec![0.0; a.nrows];
    let rep = block_solve(
        &base("direct"),
        2,
        &a,
        &vec![0.0; a.nrows],
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    assert!(rep.outer.is_none());
    assert_eq!(rep.sweeps, 0);
}

#[test]
fn unported_options_are_refused() {
    let a = coupled(4, 2);
    let b = rhs(a.nrows);
    let solve = |p: &ParamMap, dofs: usize| {
        let mut x = vec![0.0; a.nrows];
        block_solve(
            p,
            dofs,
            &a,
            &b,
            &mut x,
            &BlockOptions::default(),
            &NullLogger,
            &SerialComm,
        )
    };
    for key in [
        "Block Scaling",
        "Block Monolithic",
        "Block A-V System",
        "Linear System Complex",
    ] {
        let mut p = base("direct");
        p.insert(key, Value::Logical(true));
        assert!(matches!(solve(&p, 2), Err(BlockError::Unsupported(_))), "{key}");
    }
    let mut p = base("direct");
    p.insert("Schur Operator 1", Value::Integer(vec![2]));
    assert!(matches!(solve(&p, 2), Err(BlockError::Unsupported(_))));
    // One block: BlockSolveInt's plain DefaultSolve route.
    assert!(matches!(solve(&base("direct"), 1), Err(BlockError::Unsupported(_))));
    // The Schur variable needs its matrix.
    let mut p = base("direct");
    p.insert("Block Matrix Schur Variable", Value::Str("schur".into()));
    assert!(matches!(solve(&p, 2), Err(BlockError::Fatal { .. })));
}

#[test]
fn inner_iteration_limit_is_a_numerical_error() {
    let a = coupled(30, 2);
    let b = rhs(a.nrows);
    let mut p = base("iterative");
    p.insert("block 11: Linear System Max Iterations", Value::Integer(vec![1]));
    p.insert("block 11: Linear System Preconditioning", Value::Str("none".into()));
    let mut x = vec![0.0; a.nrows];
    let r = block_solve(
        &p,
        2,
        &a,
        &b,
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    );
    assert!(matches!(r, Err(BlockError::Numerical { .. })), "{r:?}");
    // Linear System Abort Not Converged = False lets it continue.
    p.insert("block 11: Linear System Abort Not Converged", Value::Logical(false));
    let mut x = vec![0.0; a.nrows];
    block_solve(
        &p,
        2,
        &a,
        &b,
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    assert!(residual(&a, &x, &b) < 1e-8);
}

#[test]
fn standard_iteration_converges() {
    let a = coupled(10, 2);
    let b = rhs(a.nrows);
    for gs in [false, true] {
        let mut p = base("direct");
        p.insert("Block Preconditioner", Value::Logical(false));
        p.insert("Block Gauss-Seidel", Value::Logical(gs));
        p.insert("Outer: Linear System Convergence Tolerance", Value::real(1e-12));
        let mut x = vec![0.0; a.nrows];
        let rep = block_solve(
            &p,
            2,
            &a,
            &b,
            &mut x,
            &BlockOptions::default(),
            &NullLogger,
            &SerialComm,
        )
        .unwrap();
        let out = rep.outer.unwrap();
        assert_eq!(out.status, elmer_core::SolveStatus::Converged, "gs={gs}");
        assert_eq!(rep.inner.len(), 2 * out.iterations);
        assert!(residual(&a, &x, &b) < 1e-9, "gs={gs}: {}", residual(&a, &x, &b));
    }
}

#[test]
fn assembled_blocks_match_the_split_solve() {
    // The same system given block by block (BlockSolver without a
    // Procedure) solves to the same answer as the split monolithic one.
    let a = coupled(9, 2);
    let b = rhs(a.nrows);
    let split = BlockSystem::split(&a, elmer_block::BlockLayout::by_component(2), false).unwrap();
    let blocks = (0..4).map(|k| split.block(k / 2, k % 2).cloned()).collect();
    let sys = BlockSystem::from_blocks(2, blocks).unwrap();
    let p = base("direct");
    let bb = split.gather(&b);
    let mut xb = vec![0.0; a.nrows];
    let inner = Box::new(ElmerInner::new(&NullLogger, &SerialComm, true));
    let rep = block_solve_system(
        &p,
        &sys,
        &bb,
        &mut xb,
        &BlockOptions::default(),
        inner,
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    let mut x = vec![0.0; a.nrows];
    let rep2 = block_solve(
        &p,
        2,
        &a,
        &b,
        &mut x,
        &BlockOptions::default(),
        &NullLogger,
        &SerialComm,
    )
    .unwrap();
    assert_eq!(rep.outer.unwrap().iterations, rep2.outer.unwrap().iterations);
    assert_eq!(xb, split.gather(&x));
    // A missing diagonal block is fatal.
    assert!(matches!(
        BlockSystem::from_blocks(2, vec![None, None, None, None]),
        Err(BlockError::Fatal { .. })
    ));
}
