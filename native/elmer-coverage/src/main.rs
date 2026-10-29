//! `elmer-coverage run | report | one`. See the crate docs.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use elmer_coverage::{child, discover, record, report, runner};

const DEFAULT_BUILD: &str = "/home/john/runs/elmer-fork/orig_src/build";
const DEFAULT_SCRATCH: &str = "/home/john/firstmate/data/elmer-swarm-w4coverage/scratch";

fn usage() -> ExitCode {
    eprintln!(
        "usage:\n  elmer-coverage run [--build-dir D] [--scratch D] [--out D] [--only a,b] [--timeout S] [--commit SHA] [--keep]\n  \
         elmer-coverage report --results results.tsv [--out COVERAGE.md] [--commit SHA]\n  \
         elmer-coverage one --dir D --sif S [--modules a,b] [--static]   (internal)"
    );
    ExitCode::from(2)
}

struct Args(Vec<String>);

impl Args {
    fn value(&self, flag: &str) -> Option<&str> {
        let p = self.0.iter().position(|a| a == flag)?;
        self.0.get(p + 1).map(String::as_str)
    }
    fn has(&self, flag: &str) -> bool {
        self.0.iter().any(|a| a == flag)
    }
}

fn today() -> String {
    std::process::Command::new("date").arg("+%Y-%m-%d").output().ok().map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()).unwrap_or_default()
}

fn main() -> ExitCode {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = raw.first().cloned() else { return usage() };
    let args = Args(raw[1..].to_vec());
    match cmd.as_str() {
        "one" => {
            let (Some(dir), Some(sif)) = (args.value("--dir"), args.value("--sif")) else { return usage() };
            let modules: Vec<String> = args.value("--modules").unwrap_or("").split(',').filter(|m| !m.is_empty()).map(str::to_string).collect();
            let report = child::run_one(Path::new(dir), sif, &modules, args.has("--static"));
            print!("{}", report.render());
            ExitCode::SUCCESS
        }
        "run" => {
            let build = PathBuf::from(args.value("--build-dir").unwrap_or(DEFAULT_BUILD));
            let out = PathBuf::from(args.value("--out").unwrap_or("."));
            let cfg = runner::Config {
                exe: std::env::current_exe().expect("current exe"),
                elmergrid: build.join("elmergrid/src/ElmerGrid"),
                scratch: PathBuf::from(args.value("--scratch").unwrap_or(DEFAULT_SCRATCH)),
                timeout: Duration::from_secs(args.value("--timeout").and_then(|t| t.parse().ok()).unwrap_or(120)),
                clean: !args.has("--keep"),
            };
            let only: Vec<&str> = args.value("--only").map(|o| o.split(',').collect()).unwrap_or_default();
            let tests = match discover::list(&build, "quick") {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("elmer-coverage: {e}");
                    return ExitCode::FAILURE;
                }
            };
            let mut records = Vec::new();
            for (i, t) in tests.iter().enumerate().filter(|(_, t)| only.is_empty() || only.contains(&t.name.as_str())) {
                let r = runner::run_test(&cfg, t);
                eprintln!("[{}/{}] {} {} {}", i + 1, tests.len(), t.name, r.verdict.name(), r.first().unwrap_or(&r.detail));
                records.push(r);
            }
            records.sort_by(|a, b| a.name.cmp(&b.name));
            let meta = report::Meta { commit: args.value("--commit").unwrap_or("unknown").to_string(), date: today() };
            if let Err(e) = write_outputs(&out, &records, &meta) {
                eprintln!("elmer-coverage: {e}");
                return ExitCode::FAILURE;
            }
            ExitCode::SUCCESS
        }
        "report" => {
            let Some(results) = args.value("--results") else { return usage() };
            let records = match std::fs::read_to_string(results).map_err(|e| e.to_string()).and_then(|t| record::from_tsv(&t)) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("elmer-coverage: {e}");
                    return ExitCode::FAILURE;
                }
            };
            let meta = report::Meta { commit: args.value("--commit").unwrap_or("unknown").to_string(), date: today() };
            let md = report::render(&records, &meta);
            match args.value("--out") {
                Some(p) => std::fs::write(p, md).map_err(|e| eprintln!("elmer-coverage: {e}")).map(|_| ExitCode::SUCCESS).unwrap_or(ExitCode::FAILURE),
                None => {
                    print!("{md}");
                    ExitCode::SUCCESS
                }
            }
        }
        _ => usage(),
    }
}

fn write_outputs(out: &Path, records: &[record::Record], meta: &report::Meta) -> std::io::Result<()> {
    std::fs::create_dir_all(out)?;
    std::fs::write(out.join("results.tsv"), record::to_tsv(records))?;
    std::fs::write(out.join("COVERAGE.md"), report::render(records, meta))
}
