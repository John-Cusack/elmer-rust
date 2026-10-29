//! The distributed solve under real MPI (`harness = false`). Run by cargo, this
//! binary relaunches itself with `mpirun -np 1` and `-np 2` (never more,
//! always `HWLOC_COMPONENTS=-gl`) and compares the ranks' results:
//!
//! * the partitioned Poisson fixture: np=1 and np=2 HUTI iteration counts and
//!   solutions agree;
//! * with `ELMER_DIST_CAPTURE=<dir>`, the per-rank `linsys_{a,b}.dat.R` files
//!   that a pristine `ElmerSolver_mpi` run at np=2 saved there: CG converges to
//!   the norm in that run's `np2.log` (the AngleMetisLayer quick test, see the
//!   crate README for how the capture is made).
//!
//! Skipped with a message when `mpirun` is not installed.

mod support;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use elmer_core::{Communicator, ParamMap, Value};
use elmer_dist::{
    dist_solve, owners_by_lowest_rank, read_rank_system, DistPrecond, DistSolveOptions, DofLayout,
    MpiExchange,
};
use elmer_krylov::{HutiConfig, HutiMethod};
use elmer_mpi::MpiRuntime;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("worker") {
        worker(&args[2..]);
        return;
    }
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("elmer-dist-mpi");
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).expect("tmp dir");

    let mesh = support::poisson::data("sq");
    for (method, precond) in [("cg", "none"), ("cg", "diagonal"), ("bicgstab", "none")] {
        let mut runs = Vec::new();
        for np in [1, 2] {
            let dir = out.join(format!("poisson-{method}-{precond}-np{np}"));
            let w = ["poisson", mesh.to_str().unwrap(), method, precond];
            let Some(r) = mpirun(np, &dir, &w) else {
                return;
            };
            runs.push(r);
        }
        let (one, two) = (&runs[0], &runs[1]);
        let it = |r: &[HashMap<String, String>]| r[0]["iterations"].clone();
        assert_eq!(
            it(one),
            it(two),
            "{method}/{precond}: np=1 vs np=2 iterations"
        );
        let x1 = solution(one);
        let x2 = solution(two);
        assert_eq!(x1.len(), 81);
        assert_eq!(x2.len(), 81);
        let worst = x1
            .iter()
            .map(|(g, v)| (v - x2[g]).abs())
            .fold(0.0, f64::max);
        assert!(worst < 1e-10, "{method}/{precond}: max |x1 - x2| = {worst}");
        println!(
            "poisson {method}/{precond}: np=1 and np=2 both {} iterations, max |x1 - x2| = {worst:.2e}",
            it(one)
        );
    }

    if let Some(cap) = std::env::var_os("ELMER_DIST_CAPTURE") {
        let cap = PathBuf::from(cap);
        let w = ["linsys", cap.to_str().unwrap(), "linsys"];
        let Some(r) = mpirun(2, &out.join("capture"), &w) else {
            return;
        };
        let norm: f64 = r[0]["norm"].parse().unwrap();
        let iters: usize = r[0]["iterations"].parse().unwrap();
        let (e_iters, e_norm) = elmer_log(&cap.join("np2.log"));
        println!(
            "capture {}: {iters} iterations, norm {norm:.10e}; ElmerSolver_mpi np=2: {e_iters} iterations, norm {e_norm:.8e}",
            cap.display()
        );
        assert!(
            (norm - e_norm).abs() <= 1e-7 * e_norm,
            "norm {norm} vs Elmer {e_norm}"
        );
        assert!(
            iters.abs_diff(e_iters) <= 2,
            "{iters} iterations vs Elmer {e_iters}"
        );
    } else {
        println!("ELMER_DIST_CAPTURE not set; skipping the captured np=2 system");
    }
}

/// `mpirun -np <np> <this binary> worker <w...> <dir>`; each rank's
/// `key value` lines from `<dir>/rank.R`. `None` when mpirun is missing.
fn mpirun(np: usize, dir: &Path, w: &[&str]) -> Option<Vec<HashMap<String, String>>> {
    fs::create_dir_all(dir).expect("run dir");
    let exe = std::env::current_exe().expect("own path");
    let mut cmd = Command::new("mpirun");
    cmd.args(["--oversubscribe", "-np", &np.to_string()])
        .arg(&exe)
        .arg("worker")
        .args(w)
        .arg(dir)
        .env("HWLOC_COMPONENTS", "-gl")
        .env("OMP_NUM_THREADS", "1");
    let o = match cmd.output() {
        Ok(o) => o,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            println!("mpirun not found; skipping the MPI tests");
            return None;
        }
        Err(e) => panic!("mpirun: {e}"),
    };
    let stdout = String::from_utf8_lossy(&o.stdout);
    assert!(
        o.status.success(),
        "mpirun -np {np} {w:?} failed:\n{stdout}\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
    fs::write(dir.join("stdout"), stdout.as_bytes()).expect("save stdout");
    Some(
        (0..np)
            .map(|r| {
                fs::read_to_string(dir.join(format!("rank.{r}")))
                    .expect("rank result")
                    .lines()
                    .filter_map(|l| l.split_once(' '))
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect()
            })
            .collect(),
    )
}

/// `x` entries (`gid value`, owned dofs) from every rank.
fn solution(ranks: &[HashMap<String, String>]) -> HashMap<u64, f64> {
    let mut m = HashMap::new();
    for r in ranks {
        for v in r["x"].split(';').filter(|s| !s.is_empty()) {
            let (g, x) = v.split_once('=').unwrap();
            assert!(
                m.insert(g.parse().unwrap(), x.parse().unwrap()).is_none(),
                "dof {g} owned twice"
            );
        }
    }
    m
}

/// Last HUTI iteration line and the solver's `NRM` from an ElmerSolver log.
fn elmer_log(path: &Path) -> (usize, f64) {
    let log = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut iters = None;
    let mut norm = None;
    for l in log.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() == 2 && t[1].contains('E') {
            if let (Ok(i), Ok(_)) = (t[0].parse::<usize>(), t[1].parse::<f64>()) {
                iters = Some(i);
            }
        }
        if l.starts_with("ComputeChange: NS") && norm.is_none() {
            norm = l
                .split("( ")
                .nth(1)
                .and_then(|s| s.split_whitespace().next()?.parse().ok());
        }
    }
    (
        iters.expect("HUTI lines in log"),
        norm.expect("ComputeChange line in log"),
    )
}

fn worker(a: &[String]) {
    let rt = MpiRuntime::init().expect("MPI init");
    let world = rt.world();
    let comm = MpiExchange::new(&world);
    let rank = comm.rank();
    let dir = PathBuf::from(a.last().expect("output dir"));
    let mut lines = Vec::new();
    match a[0].as_str() {
        "poisson" => {
            let p = support::poisson::assemble(&comm, Path::new(&a[1]));
            let method = match a[2].as_str() {
                "cg" => HutiMethod::Cg,
                _ => HutiMethod::BiCGStab,
            };
            let mut huti = HutiConfig::new(method);
            huti.tol = 1e-10;
            huti.max_iter = 500;
            let precond = if a[3] == "diagonal" {
                DistPrecond::Diagonal
            } else {
                DistPrecond::None
            };
            let o = DistSolveOptions {
                huti,
                precond,
                scaling: true,
            };
            let mut x = vec![0.0; p.layout.n_local()];
            let s = dist_solve(&comm, &p.layout, &p.a, &p.b, &mut x, &o, None).expect("solve");
            lines.push(format!("iterations {}", s.iterations));
            lines.push(format!("status {:?}", s.status));
            let g = p.layout.gids();
            let xs: Vec<String> = p
                .layout
                .owned()
                .iter()
                .map(|&l| format!("{}={:e}", g[l as usize], x[l as usize]))
                .collect();
            lines.push(format!("x {}", xs.join(";")));
        }
        "linsys" => {
            let sys = read_rank_system(Path::new(&a[1]), &a[2], rank).expect("read capture");
            let owners = owners_by_lowest_rank(&comm, &sys.gids);
            let layout = DofLayout::new(&comm, sys.gids.clone(), owners).expect("layout");
            // The AngleMetisLayer solver section.
            let mut prm = ParamMap::new();
            prm.insert("Linear System Iterative Method", Value::Str("cg".into()));
            prm.insert("Linear System Max Iterations", Value::Integer(vec![1000]));
            prm.insert("Linear System Convergence Tolerance", Value::real(1.0e-8));
            prm.insert("Linear System Preconditioning", Value::Str("none".into()));
            prm.insert("Linear System Residual Output", Value::Integer(vec![20]));
            let o = DistSolveOptions::from_params(&prm).expect("options");
            let mut x = vec![0.0; layout.n_local()];
            let s =
                dist_solve(&comm, &layout, &sys.matrix, &sys.rhs, &mut x, &o, None).expect("solve");
            // Nonlinear System Consistent Norm: RMS over owned dofs.
            let mut acc = [0.0, layout.n_owned() as f64];
            acc[0] = layout.owned().iter().map(|&l| x[l as usize].powi(2)).sum();
            comm.sum_f64(&mut acc);
            lines.push(format!("iterations {}", s.iterations));
            lines.push(format!("status {:?}", s.status));
            lines.push(format!("norm {:.17e}", (acc[0] / acc[1]).sqrt()));
            lines.push(format!("owned {}", layout.n_owned()));
            lines.push(format!("dofs {}", acc[1]));
        }
        w => panic!("unknown worker {w}"),
    }
    fs::write(dir.join(format!("rank.{rank}")), lines.join("\n") + "\n").expect("write result");
}
