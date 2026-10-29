//! `elmer-guicompat`: ElmerGUI process-boundary checker.
//!
//! ```text
//! elmer-guicompat fixtures --edf-dir DIR --out DIR
//! elmer-guicompat check --solver /abs/ElmerSolver [--solver-mpi /abs/ElmerSolver_mpi]
//!                       --edf-dir DIR --work DIR [--timeout SECS] [--only ID]... [--native-only]
//!                       [--strict] [--json]
//! elmer-guicompat parse-log FILE
//! ```
//!
//! Exit: 0 all good, 1 a check failed, 2 usage error. Unsupported (solver
//! exit 78) is reported but only fails the run under `--strict`.
//! `--native-only` strips the Fortran-delegation variables from the solver's
//! environment; a solve of unsupported physics is then a failure.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use elmer_guicompat::check::{run_all, Config};
use elmer_guicompat::fixtures::{build_all, EdfSet};
use elmer_guicompat::logparse::parse_log;
use elmer_guicompat::run::require_absolute;

const USAGE: &str = "usage: elmer-guicompat fixtures --edf-dir DIR --out DIR\n\
       elmer-guicompat check --solver ABS [--solver-mpi ABS] --edf-dir DIR --work DIR \
[--timeout SECS] [--only ID]... [--native-only] [--strict] [--json]\n\
       elmer-guicompat parse-log FILE";

fn usage(msg: &str) -> ExitCode {
    eprintln!("{msg}\n{USAGE}");
    ExitCode::from(2)
}

#[derive(Default)]
struct Opts {
    solver: Option<String>,
    solver_mpi: Option<String>,
    edf_dir: Option<String>,
    work: Option<String>,
    out: Option<String>,
    timeout: Option<u64>,
    only: Vec<String>,
    strict: bool,
    native_only: bool,
    json: bool,
    positional: Vec<String>,
}

fn parse(args: impl Iterator<Item = String>) -> Result<Opts, String> {
    let mut o = Opts::default();
    let mut it = args;
    while let Some(a) = it.next() {
        let mut val = |name: &str| it.next().ok_or_else(|| format!("{name} needs a value"));
        match a.as_str() {
            "--solver" => o.solver = Some(val("--solver")?),
            "--solver-mpi" => o.solver_mpi = Some(val("--solver-mpi")?),
            "--edf-dir" => o.edf_dir = Some(val("--edf-dir")?),
            "--work" => o.work = Some(val("--work")?),
            "--out" => o.out = Some(val("--out")?),
            "--only" => o.only.push(val("--only")?),
            "--timeout" => {
                o.timeout = Some(
                    val("--timeout")?
                        .parse()
                        .map_err(|_| "bad --timeout".to_string())?,
                )
            }
            "--strict" => o.strict = true,
            "--native-only" => o.native_only = true,
            "--json" => o.json = true,
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            _ => o.positional.push(a),
        }
    }
    Ok(o)
}

fn main() -> ExitCode {
    let mut args = std::env::args().skip(1);
    let Some(cmd) = args.next() else {
        return usage("missing subcommand");
    };
    let o = match parse(args) {
        Ok(o) => o,
        Err(e) => return usage(&e),
    };
    match cmd.as_str() {
        "fixtures" => {
            let (Some(edf), Some(out)) = (&o.edf_dir, &o.out) else {
                return usage("fixtures needs --edf-dir and --out");
            };
            let set = match EdfSet::load(Path::new(edf)) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("{e}");
                    return ExitCode::from(2);
                }
            };
            for fx in build_all(&set) {
                if let Err(e) = fx.write_to(&Path::new(out).join(&fx.id)) {
                    eprintln!("{}: {e}", fx.id);
                    return ExitCode::from(2);
                }
                println!("{} ({})", fx.id, fx.class.as_str());
            }
            ExitCode::SUCCESS
        }
        "check" => {
            let (Some(solver), Some(edf), Some(work)) = (&o.solver, &o.edf_dir, &o.work) else {
                return usage("check needs --solver, --edf-dir and --work");
            };
            let solver = match require_absolute(solver) {
                Ok(p) => p,
                Err(e) => return usage(&e),
            };
            let solver_mpi = match o.solver_mpi.as_deref().map(require_absolute).transpose() {
                Ok(p) => p,
                Err(e) => return usage(&e),
            };
            if let Err(e) = std::fs::create_dir_all(work) {
                eprintln!("{work}: {e}");
                return ExitCode::from(2);
            }
            let cfg = Config {
                solver,
                solver_mpi,
                work: PathBuf::from(work),
                edf_dir: PathBuf::from(edf),
                timeout: Duration::from_secs(o.timeout.unwrap_or(120)),
                only: o.only.clone(),
                native_only: o.native_only,
            };
            let rep = run_all(&cfg);
            if o.json {
                print!("{}", rep.to_json());
            } else {
                print!("{}", rep.to_text());
            }
            if rep.ok(o.strict) {
                ExitCode::SUCCESS
            } else {
                ExitCode::from(1)
            }
        }
        "parse-log" => {
            let Some(f) = o.positional.first() else {
                return usage("parse-log needs a file");
            };
            let text = match std::fs::read_to_string(f) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("{f}: {e}");
                    return ExitCode::from(2);
                }
            };
            let log = parse_log(&text);
            println!("finished: {}", log.finished);
            println!("times: {:?}", log.times);
            for c in &log.changes {
                println!(
                    "change: {} norm={:e} relc={:e} {:?}",
                    c.name, c.norm, c.relative_change, c.kind
                );
            }
            ExitCode::SUCCESS
        }
        other => usage(&format!("unknown subcommand {other}")),
    }
}
