//! One solver run, in its own process. Elmer gives every run a fresh MATC
//! environment and the native parser links the same C library, so the
//! harness never parses or solves two cases in one process.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use elmer_cli::args;
use elmer_cli::driver::{Engine, EngineOutcome, Environment, RunRequest};
use elmer_cli::logger::OutputConfig;
use elmer_cli::{ElmerLogger, ExitCode};
use elmer_run::NativeEngine;

use crate::analyze;
use crate::frontend;

/// What one run came to.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// Ran; every reference check passed (`TEST.PASSED` = 1).
    Match,
    /// Ran; a reference check failed (`TEST.PASSED` = 0) or the solver
    /// stopped on a fatal error.
    Mismatch,
    /// Ran to the end but the case has no reference to compare with.
    NoOracle,
    /// The engine declined the case.
    Refused,
    /// The process crashed, timed out or could not be set up.
    Error,
}

impl Verdict {
    pub fn name(&self) -> &'static str {
        match self {
            Verdict::Match => "match",
            Verdict::Mismatch => "mismatch",
            Verdict::NoOracle => "no-oracle",
            Verdict::Refused => "refused",
            Verdict::Error => "error",
        }
    }

    pub fn parse(s: &str) -> Option<Verdict> {
        [Verdict::Match, Verdict::Mismatch, Verdict::NoOracle, Verdict::Refused, Verdict::Error]
            .into_iter()
            .find(|v| v.name() == s)
    }

    /// Worst of two runs of one test (a test passes only if every run does).
    pub fn worst(self, other: Verdict) -> Verdict {
        let rank = |v: &Verdict| match v {
            Verdict::Match => 0,
            Verdict::NoOracle => 1,
            Verdict::Mismatch => 2,
            Verdict::Refused => 3,
            Verdict::Error => 4,
        };
        if rank(&other) > rank(&self) {
            other
        } else {
            self
        }
    }
}

/// The child's report, printed as tab-separated lines on stdout.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChildReport {
    pub verdict: Verdict,
    pub detail: String,
    pub blockers: Vec<String>,
    /// Solver library stems of the case's `Solver` sections.
    pub libs: Vec<String>,
}

impl ChildReport {
    pub fn render(&self) -> String {
        let mut s = String::new();
        for b in &self.blockers {
            s.push_str(&format!("BLOCKER\t{}\n", one_line(b)));
        }
        s.push_str(&format!("LIBS\t{}\n", self.libs.join(",")));
        s.push_str(&format!("RESULT\t{}\t{}\n", self.verdict.name(), one_line(&self.detail)));
        s
    }

    pub fn parse(stdout: &str) -> Option<ChildReport> {
        let mut blockers = Vec::new();
        let mut libs = Vec::new();
        let mut result = None;
        for line in stdout.lines() {
            if let Some(l) = line.strip_prefix("LIBS\t") {
                libs = l.split(',').filter(|s| !s.is_empty()).map(str::to_string).collect();
            } else if let Some(b) = line.strip_prefix("BLOCKER\t") {
                blockers.push(b.to_string());
            } else if let Some(r) = line.strip_prefix("RESULT\t") {
                let (v, d) = r.split_once('\t').unwrap_or((r, ""));
                result = Some((Verdict::parse(v)?, d.to_string()));
            }
        }
        let (verdict, detail) = result?;
        Some(ChildReport { verdict, detail, blockers, libs })
    }
}

fn one_line(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

#[derive(Clone, Default)]
struct Buf(Arc<Mutex<Vec<u8>>>);

impl Write for Buf {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Buf {
    fn text(&self) -> String {
        String::from_utf8_lossy(&self.0.lock().unwrap()).into_owned()
    }
}

/// The message of the log's first `ERROR::` line (`ERROR:: Caller: msg`).
fn fatal_message(log: &str) -> Option<String> {
    let line = log.lines().find(|l| l.contains("ERROR::"))?;
    Some(line.split_once("ERROR::").map(|(_, m)| m.trim().to_string()).unwrap_or_default())
}

/// Lines of the log that carry the verdict: failed comparisons and the
/// fatal error.
fn failure_lines(log: &str) -> String {
    let lines: Vec<&str> = log
        .lines()
        .filter(|l| l.contains("FAILED") || l.contains("ERROR::"))
        .map(str::trim)
        .take(3)
        .collect();
    lines.join(" | ")
}

/// Parses `sif` (relative to `dir`), screens it, and unless `static_only`
/// runs it through `NativeEngine` with the stand-in front end. The log goes
/// to `dir/native.log`.
pub fn run_one(dir: &Path, sif: &str, local_modules: &[String], static_only: bool) -> ChildReport {
    if let Err(e) = std::env::set_current_dir(dir) {
        return ChildReport { verdict: Verdict::Error, detail: format!("cd {}: {e}", dir.display()), blockers: vec![], libs: vec![] };
    }
    let sif_path = PathBuf::from(sif);
    // Staged under the test's name, which is the `fem/tests` case name.
    let case = dir.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let registry = elmer_run::solver_registry_for(case);
    let (mut blockers, libs) = match frontend::read_model(&sif_path) {
        Ok((model, refusal)) => {
            let mut b = analyze::blockers(&model, &registry, local_modules, case);
            b.extend(refusal.map(|m| analyze::from_refusal(&m)));
            (b, analyze::libs(&model))
        }
        Err(e) => {
            let blocker = analyze::from_refusal(&e);
            return ChildReport { verdict: Verdict::Refused, detail: e, blockers: vec![blocker], libs: vec![] };
        }
    };
    if static_only {
        return ChildReport { verdict: Verdict::Refused, detail: "static screen only".into(), blockers, libs };
    }

    let engine = NativeEngine::new(Box::new(frontend::frontend()));
    let buf = Buf::default();
    let log = ElmerLogger::new(OutputConfig::default(), Box::new(buf.clone()));
    let cli = args::parse(&[sif.to_string()]).expect("a sif path is a valid command line");
    let env = Environment::new(dir.to_path_buf());
    let _ = std::fs::remove_file("TEST.PASSED");
    let outcome = engine.run(&RunRequest { sif: &sif_path, args: &cli, env: &env }, &log);
    let text = buf.text();
    let _ = std::fs::write("native.log", &text);

    let (verdict, detail) = match outcome {
        EngineOutcome::Refused(reasons) => {
            for r in &reasons {
                blockers.push(analyze::from_refusal(r));
            }
            (Verdict::Refused, reasons.join("; "))
        }
        EngineOutcome::Done(code) => {
            let passed = std::fs::read_to_string("TEST.PASSED").ok().map(|s| s.trim().to_string());
            match (code, passed.as_deref()) {
                (ExitCode::Pass, Some("1")) => (Verdict::Match, String::new()),
                (ExitCode::Pass, Some(_)) => (Verdict::Mismatch, failure_lines(&text)),
                (ExitCode::Pass, None) => (Verdict::NoOracle, String::new()),
                // The solver stopped on a fatal error. One that says a piece
                // is missing is a refusal; the rest is a failed run.
                (_, _) => match fatal_message(&text) {
                    Some(m) => match analyze::capability_gap(&m) {
                        Some(gap) => {
                            blockers.push(gap);
                            (Verdict::Refused, m)
                        }
                        None => {
                            blockers.push(format!("fatal: {}", analyze::generic(&m)));
                            (Verdict::Error, format!("fatal: {m}"))
                        }
                    },
                    None => (Verdict::Error, "fatal exit without message".to_string()),
                },
            }
        }
    };
    analyze::sort_blockers(&mut blockers);
    ChildReport { verdict, detail, blockers, libs }
}
