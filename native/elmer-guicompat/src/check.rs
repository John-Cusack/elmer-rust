//! The GUI contract run: probes, typing, one run per fixture, log scraping
//! and output loading, against a solver executable given by absolute path.

use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::fixtures::{build_all, Class, EdfSet, Fixture};
use crate::logparse::parse_log;
use crate::run::run;
use crate::typecheck::check_fixture;
use crate::vtu::{load_ep, load_vtu};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Pass,
    /// The solver refused the case (exit 78). Expected for physics that has
    /// no native implementation; a gap for the ones that should have one.
    Unsupported,
    Fail,
    Skipped,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Pass => "pass",
            Status::Unsupported => "unsupported",
            Status::Fail => "FAIL",
            Status::Skipped => "skipped",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Step {
    pub name: String,
    pub status: Status,
    pub detail: String,
}

#[derive(Clone, Debug, Default)]
pub struct Report {
    pub steps: Vec<Step>,
}

impl Report {
    fn add(&mut self, name: impl Into<String>, status: Status, detail: impl Into<String>) {
        self.steps.push(Step {
            name: name.into(),
            status,
            detail: detail.into(),
        });
    }

    /// No failures; with `strict`, no refusal of a native-target fixture
    /// either (heat, elasticity, output).
    pub fn ok(&self, strict: bool) -> bool {
        self.steps.iter().all(|s| match s.status {
            Status::Fail => false,
            Status::Unsupported => {
                !(strict && s.name.starts_with("run:") && !s.name.contains("unsupported-physics"))
            }
            _ => true,
        })
    }

    pub fn to_text(&self) -> String {
        let mut o = String::new();
        for s in &self.steps {
            o.push_str(&format!(
                "[{:>11}] {} - {}\n",
                s.status.as_str(),
                s.name,
                s.detail
            ));
        }
        o
    }

    pub fn to_json(&self) -> String {
        let esc = |s: &str| {
            let mut o = String::from("\"");
            for c in s.chars() {
                match c {
                    '"' => o.push_str("\\\""),
                    '\\' => o.push_str("\\\\"),
                    '\n' => o.push_str("\\n"),
                    c if c < ' ' => o.push_str(&format!("\\u{:04x}", c as u32)),
                    c => o.push(c),
                }
            }
            o.push('"');
            o
        };
        let rows: Vec<String> = self
            .steps
            .iter()
            .map(|s| {
                format!(
                    "  {{\"name\": {}, \"status\": {}, \"detail\": {}}}",
                    esc(&s.name),
                    esc(s.status.as_str()),
                    esc(&s.detail)
                )
            })
            .collect();
        format!("[\n{}\n]\n", rows.join(",\n"))
    }
}

pub struct Config {
    pub solver: PathBuf,
    pub solver_mpi: Option<PathBuf>,
    pub work: PathBuf,
    pub edf_dir: PathBuf,
    pub timeout: Duration,
    /// Run only these fixture ids (all when empty).
    pub only: Vec<String>,
    /// Remove the Fortran-delegation variables from the solver's
    /// environment, so a solve of unsupported physics can only be a fake.
    pub native_only: bool,
}

fn first_line(log: &str) -> String {
    log.lines()
        .find(|l| !l.trim().is_empty())
        .unwrap_or("")
        .trim()
        .to_string()
}

/// Text after `Backend:` in the `-v` output, if the solver prints one.
fn backend_of(log: &str) -> Option<String> {
    log.lines()
        .find_map(|l| l.split_once("Backend:").map(|(_, b)| b.trim().to_string()))
}

pub fn run_all(cfg: &Config) -> Report {
    let mut rep = Report::default();
    if let Err(e) = std::fs::create_dir_all(&cfg.work) {
        rep.add("work", Status::Fail, format!("{}: {e}", cfg.work.display()));
        return rep;
    }

    // ElmerGUI::mainwindow.cpp:1236-1281. Note the accumulating `args`: the
    // `_mpi` probe is started as `-v -v -v`.
    match run(
        &cfg.solver,
        &["-v"],
        &cfg.work,
        Duration::from_secs(20),
        cfg.native_only,
    ) {
        Ok(o) if o.exit == Some(0) => {
            let backend = backend_of(&o.log);
            rep.add(
                "probe:ElmerSolver -v",
                Status::Pass,
                format!(
                    "{} (backend: {})",
                    first_line(&o.log),
                    backend.as_deref().unwrap_or("not reported")
                ),
            );
        }
        Ok(o) => rep.add(
            "probe:ElmerSolver -v",
            Status::Fail,
            format!("exit {:?}, timed out: {}", o.exit, o.timed_out),
        ),
        Err(e) => rep.add(
            "probe:ElmerSolver -v",
            Status::Fail,
            format!("not started: {e}"),
        ),
    }
    match &cfg.solver_mpi {
        Some(mpi) => match run(
            mpi,
            &["-v", "-v", "-v"],
            &cfg.work,
            Duration::from_secs(20),
            cfg.native_only,
        ) {
            Ok(o) if o.exit == Some(0) => rep.add(
                "probe:ElmerSolver_mpi -v -v -v",
                Status::Pass,
                first_line(&o.log),
            ),
            Ok(o) => rep.add(
                "probe:ElmerSolver_mpi -v -v -v",
                Status::Fail,
                format!("exit {:?}", o.exit),
            ),
            Err(e) => rep.add(
                "probe:ElmerSolver_mpi -v -v -v",
                Status::Fail,
                format!("not started: {e}"),
            ),
        },
        None => rep.add(
            "probe:ElmerSolver_mpi -v -v -v",
            Status::Skipped,
            "no --solver-mpi given",
        ),
    }

    let set = match EdfSet::load(&cfg.edf_dir) {
        Ok(s) => s,
        Err(e) => {
            rep.add("fixtures", Status::Fail, e);
            return rep;
        }
    };
    for fx in build_all(&set) {
        if !cfg.only.is_empty() && !cfg.only.contains(&fx.id) {
            continue;
        }
        let findings = check_fixture(
            &fx,
            set.get(if fx.id == "resultoutput" {
                "heatequation"
            } else {
                &fx.id
            }),
        );
        let fmt =
            |f: &crate::typecheck::Finding| format!("{} / {}: {}", f.section, f.key, f.problem);
        let errors: Vec<String> = findings.iter().filter(|f| !f.note).map(fmt).collect();
        let notes: Vec<String> = findings.iter().filter(|f| f.note).map(fmt).collect();
        let name = format!("typing:{}", fx.id);
        if errors.is_empty() {
            let mut d = String::from("every keyword resolves");
            if !notes.is_empty() {
                d.push_str(&format!(" (note: {})", notes.join("; ")));
            }
            rep.add(name, Status::Pass, d);
        } else {
            rep.add(name, Status::Fail, errors.join("; "));
        }
        run_fixture(cfg, &fx, &mut rep);
    }
    rep
}

fn run_fixture(cfg: &Config, fx: &Fixture, rep: &mut Report) {
    let label = if fx.class == Class::Unsupported {
        format!("run:{} (unsupported-physics)", fx.id)
    } else {
        format!("run:{}", fx.id)
    };
    let dir = cfg.work.join(&fx.id);
    let _ = std::fs::remove_dir_all(&dir);
    if let Err(e) = fx.write_to(&dir) {
        rep.add(label, Status::Fail, format!("cannot write fixture: {e}"));
        return;
    }
    // The GUI starts `ElmerSolver` with no arguments; the SIF is found via
    // ELMERSOLVER_STARTINFO in the working directory.
    let out = match run(&cfg.solver, &[], &dir, cfg.timeout, cfg.native_only) {
        Ok(o) => o,
        Err(e) => {
            rep.add(label, Status::Fail, format!("not started: {e}"));
            return;
        }
    };
    if out.timed_out {
        rep.add(
            label,
            Status::Fail,
            format!("timed out after {:?}", cfg.timeout),
        );
        return;
    }
    match out.exit {
        Some(78) => {
            rep.add(label, Status::Unsupported, first_line(&out.log));
            return;
        }
        Some(0) => {}
        other => {
            let tail: Vec<&str> = out.log.lines().rev().take(3).collect();
            rep.add(
                label,
                Status::Fail,
                format!("exit {other:?}; tail: {}", tail.join(" | ")),
            );
            return;
        }
    }
    if fx.class == Class::Unsupported && cfg.native_only {
        rep.add(
            label,
            Status::Fail,
            "exit 0 without Fortran delegation for physics with no native implementation",
        );
        return;
    }
    let log = parse_log(&out.log);
    let mut problems = Vec::new();
    if !log.finished {
        problems.push("no end-of-run line the GUI recognises".to_string());
    }
    if log.changes.is_empty() {
        problems.push("no ComputeChange line".to_string());
    }
    if log
        .changes
        .iter()
        .any(|c| !c.norm.is_finite() || !c.relative_change.is_finite())
    {
        problems.push("non-finite ComputeChange value".to_string());
    }
    let mut notes = vec![format!("{} ComputeChange lines", log.changes.len())];
    for rel in &fx.expected_outputs {
        match check_output(&dir, rel) {
            Ok(n) => notes.push(n),
            Err(e) => problems.push(format!("{rel}: {e}")),
        }
    }
    if problems.is_empty() {
        rep.add(label, Status::Pass, notes.join("; "));
    } else {
        rep.add(label, Status::Fail, problems.join("; "));
    }
}

fn check_output(dir: &Path, rel: &str) -> Result<String, String> {
    let path = dir.join(rel);
    if rel.ends_with(".vtu") {
        let bytes = std::fs::read(&path).map_err(|e| format!("missing: {e}"))?;
        let v = load_vtu(&bytes).map_err(|e| e.to_string())?;
        let t = v
            .point_array("temperature")
            .ok_or("no temperature point array")?;
        if t.values.iter().any(|x| !x.is_finite()) {
            return Err("non-finite temperature".into());
        }
        Ok(format!(
            "{rel}: {} points, {} cells, temperature loaded",
            v.points, v.cells
        ))
    } else if rel.ends_with(".ep") {
        let text = std::fs::read_to_string(&path).map_err(|e| format!("missing: {e}"))?;
        let ep = load_ep(&text).map_err(|e| e.to_string())?;
        if ep.fields.is_empty() {
            return Err("no field declared".into());
        }
        Ok(format!(
            "{rel}: {} nodes, {} elements, {} fields",
            ep.nodes,
            ep.elements,
            ep.fields.len()
        ))
    } else {
        std::fs::metadata(&path)
            .map(|_| rel.to_string())
            .map_err(|e| e.to_string())
    }
}
