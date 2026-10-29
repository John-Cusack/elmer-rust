//! Solve a system written by Elmer's `Linear System Save` with a HUTI method.
//!
//! ```text
//! cargo run -p elmer-krylov --example solve_csr -- [DIR] [key=value ...]
//! ```
//!
//! `DIR` holds `linsys_a.dat` / `linsys_b.dat` (default: the bundled 1D Poisson
//! system in `examples/data/poisson1d`). Keys: `method` (IterSolve names),
//! `prec` (none, diagonal, ilu0, ilu1, ...), `tol`, `max_iter`, `restart`,
//! `output` (residual output interval). The HUTI iteration lines go to stdout,
//! followed by a summary. The exit status is 0 only when the solve converged.

use std::path::PathBuf;
use std::process::ExitCode;

use elmer_core::{IdentityPrecond, IterativeSolver, ParamMap, Preconditioner, SerialComm, Value};
use elmer_krylov::{Huti, HutiConfig};
use elmer_sparse::{read_linear_system_save, Columns, CsrOp, Ilu};

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("solve_csr: {e}");
            ExitCode::from(2)
        }
    }
}

fn run() -> Result<bool, String> {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/data/poisson1d");
    let mut p = ParamMap::new();
    p.insert("Linear System Iterative Method", Value::Str("cg".into()));
    p.insert("Linear System Max Iterations", Value::Integer(vec![1000]));
    p.insert("Linear System Convergence Tolerance", Value::real(1e-10));
    let mut prec = "none".to_string();
    for arg in std::env::args().skip(1) {
        let Some((k, v)) = arg.split_once('=') else {
            dir = PathBuf::from(arg);
            continue;
        };
        let int = || v.parse::<i32>().map_err(|e| format!("{k}: {e}"));
        match k {
            "method" => {
                p.insert("Linear System Iterative Method", Value::Str(v.into()));
            }
            "prec" => prec = v.to_ascii_lowercase(),
            "tol" => {
                let t = v.parse::<f64>().map_err(|e| format!("tol: {e}"))?;
                p.insert("Linear System Convergence Tolerance", Value::real(t));
            }
            "max_iter" => {
                p.insert("Linear System Max Iterations", Value::Integer(vec![int()?]));
            }
            "restart" => {
                p.insert("Linear System GMRES Restart", Value::Integer(vec![int()?]));
            }
            "output" => {
                p.insert(
                    "Linear System Residual Output",
                    Value::Integer(vec![int()?]),
                );
            }
            _ => return Err(format!("unknown key {k}")),
        }
    }

    let sys = read_linear_system_save(&dir, "linsys", Columns::default())
        .map_err(|e| format!("{}: {e}", dir.display()))?;
    let a = sys.stiff.ok_or("no matrix in linsys_a.dat")?;
    let b = sys.rhs.ok_or("no linsys_b.dat")?;
    if b.len() != a.nrows {
        return Err(format!("rhs has {} rows, matrix {}", b.len(), a.nrows));
    }

    let ilu;
    let m: &dyn Preconditioner<f64> = match prec.as_str() {
        "none" => &IdentityPrecond,
        "diagonal" => {
            ilu = Ilu::diagonal(&a);
            &ilu
        }
        s => match s.strip_prefix("ilu").and_then(|l| l.parse().ok()) {
            Some(level) => {
                ilu = Ilu::ilu_n(&a, level).map_err(|e| e.to_string())?;
                &ilu
            }
            None => return Err(format!("unknown preconditioner {s}")),
        },
    };

    let cfg = HutiConfig::from_params(&p).map_err(|e| e.to_string())?;
    let method = cfg.method;
    let mut solver: Huti = Huti::new(cfg);
    let mut x = vec![0.0; b.len()];
    let stats = solver
        .solve(&CsrOp(&a), m, &b, &mut x, &SerialComm)
        .map_err(|e| e.to_string())?;

    // Independent check of the answer: ||b - A x|| / ||b||.
    let mut ax = vec![0.0; b.len()];
    elmer_sparse::ops::spmv(&a, &x, &mut ax);
    let rr: f64 = ax.iter().zip(&b).map(|(p, q)| (q - p) * (q - p)).sum();
    let bb: f64 = b.iter().map(|q| q * q).sum();

    println!(
        "system      {} ({} rows, {} nnz)",
        dir.display(),
        a.nrows,
        a.values.len()
    );
    println!("method      {} / preconditioning {prec}", method.name());
    println!("status      {:?}", stats.status);
    println!("iterations  {}", stats.iterations);
    println!("monitored   {:e}", stats.residual);
    println!("true resid  {:e}", (rr / bb).sqrt());
    Ok(stats.status == elmer_core::SolveStatus::Converged)
}
