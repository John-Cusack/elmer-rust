//! Iteration-by-iteration comparison with the pristine ElmerSolver.
//!
//! Each fixture system was dumped with `Linear System Save` from a quick test
//! and solved by the reference binary with the method edited; the expected
//! files are its stdout iteration lines verbatim. See tests/fixtures/README.md.

mod common;

use std::sync::{Arc, Mutex};

use common::{load_lines, load_system, DiagPrecond};
use elmer_core::{IdentityPrecond, IterativeSolver, ParamMap, Preconditioner, SerialComm, Value};
use elmer_krylov::{Huti, HutiConfig};

/// The keyword set the capture script wrote into the solver section.
fn params(method: &str) -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Linear System Iterative Method", Value::Str(method.into()));
    p.insert("Linear System Residual Output", Value::Integer(vec![1]));
    p.insert("Linear System Max Iterations", Value::Integer(vec![500]));
    p.insert("Linear System Convergence Tolerance", Value::real(1.0e-8));
    p.insert("Linear System Abort Not Converged", Value::Logical(false));
    p
}

/// Runs one captured configuration and returns the printed lines.
fn run(case: &str, method: &str, diag: bool, max_iter: i32) -> Vec<String> {
    let (a, b) = load_system(case);
    let mut p = params(method);
    p.insert("Linear System Max Iterations", Value::Integer(vec![max_iter]));
    let mut solver: Huti = Huti::new(HutiConfig::from_params(&p).unwrap());
    let lines = Arc::new(Mutex::new(Vec::new()));
    let sink = lines.clone();
    solver.set_output(Box::new(move |l| sink.lock().unwrap().push(l.to_owned())));
    let mut x = vec![0.0; b.len()];
    let dp;
    let m: &dyn Preconditioner<f64> = if diag {
        dp = DiagPrecond(a.diag());
        &dp
    } else {
        &IdentityPrecond
    };
    solver.solve(&a, m, &b, &mut x, &SerialComm).unwrap();
    let out = lines.lock().unwrap().clone();
    out
}

fn check(case: &str, method: &str, prec: &str, max_iter: i32) {
    let want = load_lines(case, &format!("{method}_{prec}"));
    let got = run(case, method, prec == "diagonal", max_iter);
    for (k, (g, w)) in got.iter().zip(&want).enumerate() {
        assert_eq!(g, w, "{case} {method}/{prec}: line {} differs", k + 1);
    }
    assert_eq!(
        got.len(),
        want.len(),
        "{case} {method}/{prec}: line count (last ours {:?}, Fortran {:?})",
        got.last(),
        want.last()
    );
}

macro_rules! parity {
    ($name:ident, $case:literal, $method:literal, $prec:literal, $max:literal) => {
        #[test]
        fn $name() {
            check($case, $method, $prec, $max);
        }
    };
}

// HeatAniso (symmetric, 1001 dofs): CG/BiCGStab/BiCGStab(2)/GMRES converge,
// CGS overflows to Infinity (divergence), TFQMR stagnates to the limit.
parity!(heataniso_cg, "heataniso", "cg", "none", 2000);
parity!(heataniso_cgs, "heataniso", "cgs", "none", 2000);
parity!(heataniso_bicgstab, "heataniso", "bicgstab", "none", 2000);
parity!(heataniso_bicgstab2, "heataniso", "bicgstab2", "none", 2000);
parity!(heataniso_gmres, "heataniso", "gmres", "none", 2000);
parity!(heataniso_tfqmr, "heataniso", "tfqmr", "none", 500);

// HeatAniso with constant convection (nonsymmetric), no preconditioner.
parity!(conv_cg, "heataniso_conv", "cg", "none", 500);
parity!(conv_cgs, "heataniso_conv", "cgs", "none", 500);
parity!(conv_bicgstab, "heataniso_conv", "bicgstab", "none", 500);
parity!(conv_bicgstab2, "heataniso_conv", "bicgstab2", "none", 500);
parity!(conv_gmres, "heataniso_conv", "gmres", "none", 500);
parity!(conv_tfqmr, "heataniso_conv", "tfqmr", "none", 500);

// Same system with the diagonal preconditioner: right-applied for CG, CGS and
// BiCGStab, left-applied for BiCGStab(2), GMRES and TFQMR.
parity!(conv_diag_cg, "heataniso_conv", "cg", "diagonal", 500);
parity!(conv_diag_cgs, "heataniso_conv", "cgs", "diagonal", 500);
parity!(conv_diag_bicgstab, "heataniso_conv", "bicgstab", "diagonal", 500);
parity!(conv_diag_bicgstab2, "heataniso_conv", "bicgstab2", "diagonal", 500);
parity!(conv_diag_gmres, "heataniso_conv", "gmres", "diagonal", 500);
parity!(conv_diag_tfqmr, "heataniso_conv", "tfqmr", "diagonal", 500);
