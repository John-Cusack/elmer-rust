//! Test-only stand-ins (ADR 003, rule 6) for S1's `elmer-mesh` (mesh-file
//! reader), `elmer-fe` (basis, quadrature; copied from `elmer-run`'s tests,
//! plus the cubic line element) and a logger that keeps what Elmer would
//! print. Replace them with the real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop)]

pub mod basis;
pub mod fixtures;
pub mod mesh;
pub mod poisson;
pub mod quad;

use std::sync::Mutex;

use elmer_core::{ParamMap, Value};
use elmer_solver::Model;

/// Collects every message as `caller: msg`, warnings with Elmer's prefix.
#[derive(Default)]
pub struct Capture(Mutex<Vec<String>>);

impl elmer_core::Logger for Capture {
    fn info(&self, caller: &str, msg: &str, _level: u8) {
        self.0.lock().unwrap().push(format!("{caller}: {msg}"));
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.lock().unwrap().push(format!("WARNING:: {caller}: {msg}"));
    }
}

impl Capture {
    /// Lines starting with any of `prefixes`, in order.
    pub fn lines(&self, prefixes: &[&str]) -> Vec<String> {
        let v = self.0.lock().unwrap();
        v.iter().filter(|l| prefixes.iter().any(|p| l.starts_with(p))).cloned().collect()
    }
}

/// A section from `(key, value)` pairs.
pub fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

pub fn s(v: &str) -> Value {
    Value::Str(v.into())
}

pub fn i(v: &[i32]) -> Value {
    Value::Integer(v.to_vec())
}

pub fn r(v: f64) -> Value {
    Value::real(v)
}

/// Splits a `ComputeChange` line into its fixed fields and the RELC value,
/// which an iterative solve reproduces only to a few digits when it is tiny.
pub fn split_change(l: &str) -> (String, String, f64, String) {
    let (head, rest) = l.split_once("(NRM,RELC): (").unwrap();
    let nrm = rest[..15].to_string();
    let relc: f64 = rest[15..30].trim().parse().unwrap();
    (head.to_string(), nrm, relc, rest[30..].to_string())
}

/// The sections of a `.sif` as the solver framework reads them, through
/// `elmer-sif`'s parser (as `elmer-savedata`'s tests do). Elmer's
/// `ListAddString` lowercases string values; `elmer-sif` keeps the case of
/// quoted text, so this lowercases them.
pub fn sif_model(path: &std::path::Path) -> Model {
    let parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
    to_model(&parsed.sif.merged())
}

/// One `Model` per `RUN` block, each holding only that block's entries.
/// A block with no entries after the last `RUN` is dropped: Elmer's
/// `ReloadInputFile` finds nothing and the session ends.
pub fn sif_runs(path: &std::path::Path) -> Vec<Model> {
    let parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
    let mut runs: Vec<Model> = parsed.sif.runs.iter().map(to_model).collect();
    let empty = |m: &Model| {
        [&m.simulation, &m.constants].iter().all(|p| p.0.is_empty())
            && [
                &m.solvers,
                &m.equations,
                &m.bodies,
                &m.materials,
                &m.body_forces,
                &m.initial_conditions,
                &m.boundary_conditions,
            ]
            .iter()
            .all(|v| v.iter().all(|p| p.0.is_empty()))
    };
    while runs.len() > 1 && runs.last().is_some_and(empty) {
        runs.pop();
    }
    runs
}

fn to_model(m: &elmer_sif::Model) -> Model {
    let convert = |s: &elmer_sif::Section| -> ParamMap {
        let mut p = elmer_sif::section_to_param_map(s);
        for v in p.0.values_mut() {
            if let Value::Str(t) = v {
                *t = t.to_lowercase();
            }
        }
        p
    };
    let numbered = |v: &[elmer_sif::Section]| -> Vec<ParamMap> {
        let n = v.iter().map(|s| s.id as usize).max().unwrap_or(0);
        let mut out = vec![ParamMap::new(); n];
        for s in v {
            out[s.id as usize - 1] = convert(s);
        }
        out
    };
    let one = |s: &Option<elmer_sif::Section>| s.as_ref().map(convert).unwrap_or_default();
    Model {
        header: ParamMap::new(),
        simulation: one(&m.simulation),
        constants: one(&m.constants),
        solvers: numbered(&m.solvers),
        equations: numbered(&m.equations),
        bodies: numbered(&m.bodies),
        materials: numbered(&m.materials),
        body_forces: numbered(&m.body_forces),
        initial_conditions: numbered(&m.ics),
        boundary_conditions: numbered(&m.bcs),
    }
}
