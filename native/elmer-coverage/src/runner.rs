//! Runs the quick suite: stage each test, run its solver steps in child
//! processes, and fold the child reports into a [`Record`].

use std::fs;
use std::path::{Path, PathBuf};
use std::os::unix::process::CommandExt;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crate::analyze::{self, sort_blockers};
use crate::child::{failure_lines, fatal_message, ChildReport, Verdict};
use crate::discover::QuickTest;
use crate::record::Record;
use crate::stage::{self, Recipe, Step};

pub struct Config {
    /// The harness binary, re-run as the per-solve child.
    pub exe: PathBuf,
    pub elmergrid: PathBuf,
    pub scratch: PathBuf,
    /// Per solver run.
    pub timeout: Duration,
    /// Delete a staged directory once its test is classified.
    pub clean: bool,
    /// `elmer-mpi-run`, which runs the tests with more than one process
    /// under `mpirun`; without it they are refused as `mpi: np=N`.
    pub mpi_run: Option<PathBuf>,
}

/// `ELMERSOLVER_STARTINFO`'s first line: the case the test runs.
fn startinfo(dir: &Path) -> Option<String> {
    let text = fs::read_to_string(dir.join("ELMERSOLVER_STARTINFO")).ok()?;
    text.lines().next().map(|l| l.trim().to_string()).filter(|l| !l.is_empty())
}

/// `src/sif.rs:33` out of `thread 'main' panicked at src/sif.rs:33:12:`.
fn panic_location(stderr: &str) -> Option<String> {
    let rest = stderr.lines().find_map(|l| l.split_once("panicked at ").map(|(_, r)| r))?;
    let loc = rest.trim_end_matches(':');
    let mut parts = loc.rsplitn(3, ':');
    let (_col, line, file) = (parts.next()?, parts.next()?, parts.next()?);
    Some(format!("{file}:{line}"))
}

/// The panic line and the message after it.
fn panic_summary(stderr: &str) -> String {
    let mut lines = stderr.lines().skip_while(|l| !l.contains("panicked at "));
    let at = lines.next().unwrap_or("").trim();
    let msg = lines.next().unwrap_or("").trim();
    if at.is_empty() { stderr.lines().next().unwrap_or("").to_string() } else { format!("{at} {msg}") }
}

/// Spawns `exe one ...` and waits up to `limit`.
fn run_child(cfg: &Config, dir: &Path, sif: &str, recipe: &Recipe, static_only: bool) -> ChildReport {
    let out_path = dir.join("child.out");
    let err_path = dir.join("child.err");
    let fail = |detail: String| ChildReport { verdict: Verdict::Error, detail, blockers: vec![], libs: vec![] };
    let (Ok(out), Ok(err)) = (fs::File::create(&out_path), fs::File::create(&err_path)) else {
        return fail("cannot create child output files".into());
    };
    let mut cmd = Command::new("nice");
    cmd.args(["-n", "10"])
        .arg(&cfg.exe)
        .arg("one")
        .arg("--dir")
        .arg(dir)
        .arg("--sif")
        .arg(sif)
        .arg("--modules")
        .arg(recipe.local_modules.join(","))
        .env("OMP_NUM_THREADS", "1")
        .env("RAYON_NUM_THREADS", "1")
        .env("RUST_BACKTRACE", "0")
        .stdin(Stdio::null())
        .stdout(out)
        .stderr(err);
    if static_only {
        cmd.arg("--static");
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) => return fail(format!("spawn: {e}")),
    };
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(st)) => break Some(st),
            Ok(None) if start.elapsed() > cfg.timeout => {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(10)),
            Err(e) => return fail(format!("wait: {e}")),
        }
    };
    let stdout = fs::read_to_string(&out_path).unwrap_or_default();
    let stderr = fs::read_to_string(&err_path).unwrap_or_default();
    let mut report = ChildReport::parse(&stdout);
    match (status, &mut report) {
        (None, r) => {
            let mut r = r.take().unwrap_or_else(|| fail(String::new()));
            r.verdict = Verdict::Error;
            r.detail = format!("timeout after {}s", cfg.timeout.as_secs());
            r.blockers.push("crash: timeout".to_string());
            r
        }
        (Some(_), Some(r)) => r.clone(),
        (Some(st), None) => {
            let mut r = fail(format!("child {st}: {}", panic_summary(&stderr)));
            if let Some(loc) = panic_location(&stderr) {
                r.blockers.push(format!("crash: panic at {loc}"));
            }
            r
        }
    }
}

/// `mpirun -np N elmer-mpi-run sif` in `dir`, waiting up to the timeout.
/// Rank 0 writes `TEST.PASSED_N`; the exit codes are `elmer-cli`'s (0
/// solved, 78 refused with `elmer-mpi-run: refused: ...` lines on stderr).
fn run_mpi(cfg: &Config, mpi_run: &Path, dir: &Path, sif: &str, np: u32) -> ChildReport {
    let fail = |detail: String| ChildReport { verdict: Verdict::Error, detail, blockers: vec![], libs: vec![] };
    let (out_path, err_path) = (dir.join("mpi.out"), dir.join("mpi.err"));
    let (Ok(out), Ok(err)) = (fs::File::create(&out_path), fs::File::create(&err_path)) else {
        return fail("cannot create mpirun output files".into());
    };
    let passed = dir.join(format!("TEST.PASSED_{np}"));
    let _ = fs::remove_file(&passed);
    let mut cmd = Command::new("nice");
    cmd.args(["-n", "10", "mpirun", "--oversubscribe", "-np", &np.to_string()])
        .arg(mpi_run)
        .arg(sif)
        .current_dir(dir)
        .env("HWLOC_COMPONENTS", "-gl")
        .env("OMP_NUM_THREADS", "1")
        .env("RAYON_NUM_THREADS", "1")
        .env("RUST_BACKTRACE", "0")
        .stdin(Stdio::null())
        .stdout(out)
        .stderr(err)
        // Its own process group, so a timeout kills the ranks with mpirun.
        .process_group(0);
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) => return fail(format!("spawn mpirun: {e}")),
    };
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(st)) => break st,
            Ok(None) if start.elapsed() > cfg.timeout => {
                let _ = Command::new("kill").args(["-KILL", "--", &format!("-{}", child.id())]).status();
                let _ = child.wait();
                let mut r = fail(format!("timeout after {}s", cfg.timeout.as_secs()));
                r.blockers.push("crash: timeout".to_string());
                return r;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(10)),
            Err(e) => return fail(format!("wait: {e}")),
        }
    };
    let stdout = fs::read_to_string(&out_path).unwrap_or_default();
    let stderr = fs::read_to_string(&err_path).unwrap_or_default();
    let refusals: Vec<&str> = stderr.lines().filter_map(|l| l.split_once("refused: ").map(|(_, r)| r.trim())).collect();
    let passed = fs::read_to_string(&passed).ok().map(|s| s.trim().to_string());
    let mut r = ChildReport { verdict: Verdict::Error, detail: String::new(), blockers: vec![], libs: vec![] };
    match (status.code(), passed.as_deref()) {
        (Some(78), _) if !refusals.is_empty() => {
            r.blockers = refusals.iter().map(|m| analyze::from_refusal(m)).collect();
            r.verdict = Verdict::Refused;
            r.detail = refusals.join("; ");
        }
        (Some(0), Some("1")) => r.verdict = Verdict::Match,
        (Some(0), Some(_)) => {
            r.verdict = Verdict::Mismatch;
            r.detail = failure_lines(&stdout);
        }
        (Some(0), None) => r.verdict = Verdict::NoOracle,
        (code, _) => match fatal_message(&stdout) {
            Some(m) => match analyze::capability_gap(&m) {
                Some(gap) => {
                    r.blockers.push(gap);
                    r.verdict = Verdict::Refused;
                    r.detail = m;
                }
                None => {
                    r.blockers.push(format!("fatal: {}", analyze::generic(&m)));
                    r.detail = format!("fatal: {m}");
                }
            },
            None => {
                let last = stderr.lines().rev().find(|l| !l.trim().is_empty()).unwrap_or("");
                r.detail = format!("mpirun exit {code:?}: {last}");
            }
        },
    }
    r
}

/// Runs one test and classifies it.
pub fn run_test(cfg: &Config, test: &QuickTest) -> Record {
    let recipe = stage::recipe(test);
    let mut blockers = Vec::new();
    let mpi_run = cfg.mpi_run.as_deref().filter(|_| test.nprocs > 1);
    if test.nprocs > 1 && mpi_run.is_none() {
        blockers.push(format!("mpi: np={}", test.nprocs));
    }
    for s in &recipe.steps {
        if let Step::Other(tool) = s {
            blockers.push(format!("prestep: {tool}"));
        }
    }
    // Only a test that can be set up runs; the others get the static screen.
    let static_only = !blockers.is_empty();
    let mut record = Record {
        name: test.name.clone(),
        nprocs: test.nprocs,
        verdict: Verdict::Refused,
        secs: 0.0,
        libs: vec![],
        blockers: vec![],
        detail: String::new(),
    };
    let mesh_steps: Vec<Step> = if static_only { vec![] } else { recipe.steps.clone() };
    let dir = match stage::stage(test, &mesh_steps, &cfg.scratch, &cfg.elmergrid) {
        Ok(d) => d,
        Err(e) => {
            record.verdict = Verdict::Error;
            record.detail = format!("stage: {e}");
            record.blockers = blockers;
            return record;
        }
    };
    let mut sifs: Vec<String> = recipe
        .steps
        .iter()
        .filter_map(|s| match s {
            Step::Solve(Some(sif)) => Some(Some(sif.clone())),
            Step::Solve(None) => Some(startinfo(&dir)),
            _ => None,
        })
        .map(|s| s.unwrap_or_default())
        .collect();
    if sifs.is_empty() {
        sifs.push(startinfo(&dir).unwrap_or_default());
    }
    if sifs.iter().any(String::is_empty) {
        record.verdict = Verdict::Error;
        record.detail = "no ELMERSOLVER_STARTINFO".into();
        record.blockers = blockers;
        return record;
    }

    let start = Instant::now();
    let mut verdict: Option<Verdict> = None;
    let mut libs = Vec::new();
    for sif in &sifs {
        let r = match mpi_run {
            // The static screen first: a case it refuses is not run.
            Some(bin) if !static_only => {
                let screen = run_child(cfg, &dir, sif, &recipe, true);
                if screen.verdict == Verdict::Refused && screen.blockers.is_empty() {
                    ChildReport { libs: screen.libs, ..run_mpi(cfg, bin, &dir, sif, test.nprocs) }
                } else {
                    screen
                }
            }
            _ => run_child(cfg, &dir, sif, &recipe, static_only),
        };
        verdict = Some(match verdict {
            None => r.verdict.clone(),
            Some(v) => v.worst(r.verdict.clone()),
        });
        if !r.detail.is_empty() && record.detail.is_empty() {
            record.detail = if sifs.len() > 1 { format!("{sif}: {}", r.detail) } else { r.detail.clone() };
        }
        blockers.extend(r.blockers);
        libs.extend(r.libs);
    }
    record.secs = start.elapsed().as_secs_f64();
    record.verdict = verdict.unwrap_or(Verdict::Error);
    if static_only && record.verdict != Verdict::Error {
        record.verdict = Verdict::Refused;
    }
    libs.sort();
    libs.dedup();
    record.libs = libs;
    // A refusal's runtime message restates a solver blocker the screen
    // already named; keep the screen's.
    let has_solver = blockers.iter().any(|b| b.starts_with("solver:") && !b.contains("(no Procedure)"));
    blockers.retain(|b| !(has_solver && b.contains("(no Procedure)")));
    sort_blockers(&mut blockers);
    record.blockers = blockers;
    if cfg.clean {
        let _ = fs::remove_dir_all(&dir);
    }
    record
}
