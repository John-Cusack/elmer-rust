//! Test-only stand-ins (ADR 003, rule 6) for S1's `elmer-mesh` (mesh-file
//! reader), `elmer-fe` (basis, quadrature; copied from `elmer-run`'s tests,
//! plus the cubic line element) and a logger that keeps what Elmer would
//! print. Replace them with the real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop)]

pub mod basis;
pub mod mesh;
pub mod poisson;
pub mod quad;

use std::sync::Mutex;

use elmer_core::{ParamMap, Value};

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
