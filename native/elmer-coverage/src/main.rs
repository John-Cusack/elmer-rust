//! `elmer-coverage run | report | one`. See the crate docs.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use elmer_coverage::{child, discover, record, report, runner};

const DEFAULT_BUILD: &str = "/home/john/runs/elmer-fork/orig_src/build";

const USAGE: &str = "usage:\n  \
     elmer-coverage run [--build-dir D] [--scratch D] [--out D] [--ctest-cache F] [--only a,b] [--timeout S] [--commit SHA] [--keep] [--mpi-run BIN]\n                       \
     [--golden ROOT --tolerances JSON --oracle BIN]   (field check of matches)\n  \
     elmer-coverage report --results results.tsv [--out COVERAGE.md] [--commit SHA]\n  \
     elmer-coverage one --dir D --sif S [--modules a,b] [--static]   (internal)\n\n\
     --scratch defaults to a new directory under the system temp dir. A scratch\n\
     directory must be new, empty or one elmer-coverage created (it holds a\n\
     .elmer-coverage-scratch marker): nothing is deleted anywhere else.\n";

fn usage() -> ExitCode {
    eprint!("{USAGE}");
    ExitCode::from(2)
}

/// The flags of each subcommand: those that take a value, and switches.
fn flags(cmd: &str) -> Option<(&'static [&'static str], &'static [&'static str])> {
    Some(match cmd {
        "run" => (
            &[
                "--build-dir",
                "--scratch",
                "--out",
                "--ctest-cache",
                "--only",
                "--timeout",
                "--commit",
                "--mpi-run",
                "--golden",
                "--tolerances",
                "--oracle",
            ],
            &["--keep"],
        ),
        "report" => (&["--results", "--out", "--commit"], &[]),
        "one" => (&["--dir", "--sif", "--modules"], &["--static"]),
        _ => return None,
    })
}

/// What the command line asks for, before anything runs.
#[derive(Debug, PartialEq)]
enum Parsed {
    /// `-h` / `--help` anywhere, or no command: print the usage, nothing else.
    Help,
    /// A known command whose flags all check out.
    Run(String, Vec<String>),
    /// An unknown command or flag, a flag without its value, or a stray word.
    Bad(String),
}

fn parse(raw: &[String]) -> Parsed {
    if raw.is_empty() || raw.iter().any(|a| a == "-h" || a == "--help") {
        return Parsed::Help;
    }
    let cmd = raw[0].as_str();
    let Some((with_value, switches)) = flags(cmd) else {
        return Parsed::Bad(format!("unknown command {cmd:?}"));
    };
    let rest = &raw[1..];
    let mut i = 0;
    while i < rest.len() {
        let a = rest[i].as_str();
        if with_value.contains(&a) {
            if rest.get(i + 1).is_none_or(|v| v.starts_with("--")) {
                return Parsed::Bad(format!("{a} needs a value"));
            }
            i += 2;
        } else if switches.contains(&a) {
            i += 1;
        } else {
            return Parsed::Bad(format!("unknown argument {a:?} for {cmd}"));
        }
    }
    Parsed::Run(cmd.to_string(), rest.to_vec())
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
    let (cmd, args) = match parse(&raw) {
        Parsed::Help => {
            print!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Parsed::Bad(why) => {
            eprintln!("elmer-coverage: {why}");
            return usage();
        }
        Parsed::Run(cmd, rest) => (cmd, Args(rest)),
    };
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
                scratch: args
                    .value("--scratch")
                    .map(PathBuf::from)
                    .unwrap_or_else(|| std::env::temp_dir().join(format!("elmer-coverage-{}", std::process::id()))),
                timeout: Duration::from_secs(args.value("--timeout").and_then(|t| t.parse().ok()).unwrap_or(120)),
                clean: !args.has("--keep"),
                mpi_run: args.value("--mpi-run").map(PathBuf::from).or_else(|| {
                    let sibling = std::env::current_exe().ok()?.with_file_name("elmer-mpi-run");
                    sibling.is_file().then_some(sibling)
                }),
                fields: match (args.value("--golden"), args.value("--tolerances"), args.value("--oracle")) {
                    (Some(g), Some(t), Some(o)) => Some(elmer_coverage::fields::FieldCheck {
                        golden: PathBuf::from(g),
                        tolerances: PathBuf::from(t),
                        oracle: PathBuf::from(o),
                    }),
                    (None, None, None) => None,
                    _ => return usage(),
                },
            };
            if let Err(e) = elmer_coverage::stage::prepare_scratch(&cfg.scratch) {
                eprintln!("elmer-coverage: {e}");
                return ExitCode::FAILURE;
            }
            let only: Vec<&str> = args.value("--only").map(|o| o.split(',').collect()).unwrap_or_default();
            let cache = args.value("--ctest-cache").map(PathBuf::from).unwrap_or_else(|| out.join("ctest-quick.txt"));
            let tests = match discover::list_cached(&build, "quick", &cache) {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    /// `--help` / `-h` anywhere and an empty command line only print the
    /// usage; unknown commands or flags and missing values are refused
    /// before anything runs.
    #[test]
    fn help_and_unknown_flags_have_no_side_effects() {
        assert_eq!(parse(&v(&["run", "--help"])), Parsed::Help);
        assert_eq!(parse(&v(&["-h"])), Parsed::Help);
        assert_eq!(parse(&v(&[])), Parsed::Help);
        assert_eq!(parse(&v(&["run", "--only", "a", "-h"])), Parsed::Help);
        assert!(matches!(parse(&v(&["run", "--frobnicate"])), Parsed::Bad(_)));
        assert!(matches!(parse(&v(&["run", "--scratch"])), Parsed::Bad(_)));
        assert!(matches!(parse(&v(&["run", "--scratch", "--keep"])), Parsed::Bad(_)));
        assert!(matches!(parse(&v(&["run", "stray"])), Parsed::Bad(_)));
        assert!(matches!(parse(&v(&["delete-everything"])), Parsed::Bad(_)));
        assert!(matches!(parse(&v(&["report", "--keep"])), Parsed::Bad(_)));
        assert_eq!(
            parse(&v(&["run", "--scratch", "/tmp/x", "--keep", "--only", "a,b"])),
            Parsed::Run("run".into(), v(&["--scratch", "/tmp/x", "--keep", "--only", "a,b"]))
        );
        assert_eq!(parse(&v(&["one", "--dir", "d", "--sif", "s", "--static"])).clone_name(), Some("one"));
    }

    impl Parsed {
        fn clone_name(&self) -> Option<&str> {
            match self {
                Parsed::Run(c, _) => Some(c),
                _ => None,
            }
        }
    }
}
