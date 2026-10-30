//! σ_F measurement (report §4.1 O3): Fortran's own sensitivity to
//! legitimate reordering, per (test, result file, field, saved step).
//!
//! The golden arm's captured fields are compared against each sensitivity
//! arm's capture of the same test; σ_F is the largest relative-Linf
//! distance in node space. A field whose structure differs between arms
//! (missing variable or step, different length) has no finite σ_F: it is
//! recorded as `structure-differs` and cannot be graded at Tier A.

use crate::grade::{relative_linf_floor, to_node_space, variable_scale};
use crate::result::{read_result_file, ResultFile};
use crate::tolerance::ToleranceTable;
use crate::Result;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// One measured field.
#[derive(Debug, Clone, PartialEq)]
pub struct SigmaEntry {
    pub test: String,
    /// Result file path relative to the run dir (`square/case.result`).
    pub file: String,
    pub field: String,
    /// Which copy of `field` in the step (0 unless the name repeats).
    pub occurrence: usize,
    pub step: usize,
    /// `None` when the arms disagree structurally.
    pub sigma_f: Option<f64>,
    /// Arm that set σ_F (largest distance), or the arm that broke structure.
    pub worst_arm: String,
}

impl SigmaEntry {
    /// Tier-A bar and the rule that produced it.
    pub fn bar(&self) -> Option<(f64, &'static str)> {
        self.sigma_f.map(ToleranceTable::tier_a_bar)
    }
}

/// Whether `bytes` start like an Elmer result file of either form.
pub fn looks_like_result(bytes: &[u8]) -> bool {
    let head = &bytes[..bytes.len().min(16)];
    let head = String::from_utf8_lossy(head);
    let head = head.trim_start();
    head.starts_with("BINARY 3.") || head.starts_with("ASCII 3")
}

/// Result files a captured run wrote: every file under `run_dir` that
/// parses as a result header and is not a staged input (a path that also
/// exists in the test's source dir, such as a restart file).
pub fn result_files(run_dir: &Path, source_dir: &Path) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let mut stack = vec![run_dir.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir).map_err(|e| crate::io_with(&dir, e))?;
        for entry in entries {
            let entry = entry.map_err(|e| crate::io_with(&dir, e))?;
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            let Ok(rel) = path.strip_prefix(run_dir) else {
                continue;
            };
            let rel = rel.display().to_string();
            if rel.ends_with(".pos") || source_dir.join(&rel).exists() {
                continue;
            }
            let mut head = [0u8; 16];
            let n = std::fs::File::open(&path)
                .and_then(|mut f| std::io::Read::read(&mut f, &mut head))
                .unwrap_or(0);
            if looks_like_result(&head[..n]) {
                out.push(rel);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Compare one golden result file against the same file from each arm.
pub fn measure_file(
    test: &str,
    file: &str,
    golden: &ResultFile,
    arms: &[(String, Option<ResultFile>)],
) -> Vec<SigmaEntry> {
    let mut out = Vec::new();
    for (step, g_step) in golden.steps.iter().enumerate() {
        for (k, g_var) in g_step.vars.iter().enumerate() {
            let occurrence = g_step.vars[..k].iter().filter(|v| v.name == g_var.name).count();
            let mut entry = SigmaEntry {
                test: test.to_string(),
                file: file.to_string(),
                field: g_var.name.clone(),
                occurrence,
                step,
                sigma_f: Some(0.0),
                worst_arm: String::new(),
            };
            let g_nodes = to_node_space(g_var).ok();
            let floor = variable_scale(golden, step, &g_var.name);
            for (arm, arm_file) in arms {
                let a_var = arm_file
                    .as_ref()
                    .and_then(|f| f.steps.get(step))
                    .and_then(|s| s.vars.iter().filter(|v| v.name == g_var.name).nth(occurrence));
                let a_nodes = a_var.and_then(|v| to_node_space(v).ok());
                let d = match (&g_nodes, &a_nodes) {
                    (Some(g), Some(a)) => relative_linf_floor(g, a, floor),
                    _ => None,
                };
                match d {
                    None => {
                        entry.sigma_f = None;
                        entry.worst_arm = arm.clone();
                        break;
                    }
                    Some(d) => {
                        if entry.worst_arm.is_empty() || entry.sigma_f.is_some_and(|s| d > s) {
                            entry.sigma_f = Some(d);
                            entry.worst_arm = arm.clone();
                        }
                    }
                }
            }
            out.push(entry);
        }
    }
    // A step the arm saved beyond the golden's count is structure too.
    for (arm, arm_file) in arms {
        if let Some(f) = arm_file {
            if f.steps.len() != golden.steps.len() {
                out.push(SigmaEntry {
                    test: test.to_string(),
                    file: file.to_string(),
                    field: String::from("<steps>"),
                    occurrence: 0,
                    step: golden.steps.len().min(f.steps.len()),
                    sigma_f: None,
                    worst_arm: arm.clone(),
                });
            }
        }
    }
    out
}

/// Measure every result file of one captured test. `golden_run` and each
/// arm's run dir hold the same test captured on different arms.
pub fn measure_test(
    test: &str,
    golden_run: &Path,
    source_dir: &Path,
    arm_runs: &[(String, PathBuf)],
) -> Result<Vec<SigmaEntry>> {
    let mut out = Vec::new();
    for file in result_files(golden_run, source_dir)? {
        let golden = match read_result_file(&golden_run.join(&file)) {
            Ok(g) => g,
            Err(e) => {
                return Err(crate::OracleError::Parse(format!("{test}/{file}: {e}")));
            }
        };
        let arms: Vec<(String, Option<ResultFile>)> = arm_runs
            .iter()
            .map(|(name, run)| (name.clone(), read_result_file(&run.join(&file)).ok()))
            .collect();
        out.extend(measure_file(test, &file, &golden, &arms));
    }
    Ok(out)
}

/// Per-test roll-up: the worst σ_F over every field and step (infinite
/// when any field is structure-differs), with the field count.
pub fn per_test_worst(entries: &[SigmaEntry]) -> BTreeMap<String, (f64, usize)> {
    let mut out: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    for e in entries {
        let s = e.sigma_f.unwrap_or(f64::INFINITY);
        let slot = out.entry(e.test.clone()).or_insert((0.0, 0));
        slot.0 = slot.0.max(s);
        slot.1 += 1;
    }
    out
}

/// Serialize entries as `tolerances.json` (version 1). Structure-differs
/// fields are listed with a null bar so the grader refuses them instead of
/// falling back silently.
pub fn tolerances_json(
    entries: &[SigmaEntry],
    arms: &serde_json::Value,
    source_pin: &str,
    formula_inputs: &str,
) -> String {
    let mut lines = Vec::new();
    for e in entries {
        let (sigma, bar, rule) = match e.bar() {
            Some((bar, rule)) => (
                serde_json::json!(e.sigma_f.unwrap_or(0.0)),
                serde_json::json!(bar),
                rule,
            ),
            None => (serde_json::Value::Null, serde_json::Value::Null, "structure-differs"),
        };
        let mut obj = serde_json::json!({
            "test": e.test,
            "file": e.file,
            "field": e.field,
            "step": e.step,
            "sigma_f": sigma,
            "bar": bar,
            "rule": rule,
            "worst_arm": e.worst_arm,
        });
        if e.occurrence > 0 {
            obj["occurrence"] = serde_json::json!(e.occurrence);
        }
        lines.push(format!("    {obj}"));
    }
    let head = serde_json::json!({
        "version": 1,
        "note": "Tier-A bars from measured Fortran sensitivity to legitimate reordering. bar = max(10.sigma_f, cap); cap = 1e-10 applies only where every arm agrees bitwise; rule structure-differs means no Tier-A bar exists for that field. Changes only by ADR with captain sign-off.",
        "arms": arms,
        "source_pin": source_pin,
        "formula_inputs": formula_inputs,
    });
    let head = serde_json::to_string_pretty(&head).unwrap_or_default();
    let head = head.trim_end().trim_end_matches('}').trim_end();
    format!("{head},\n  \"entries\": [\n{}\n  ]\n}}\n", lines.join(",\n"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result::{TimeStep, VariableData};

    fn file(values: Vec<f64>) -> ResultFile {
        ResultFile {
            groups: Vec::new(),
            steps: vec![TimeStep {
                saves_done: 1,
                ntime: 1,
                sim_time: 0.0,
                vars: vec![VariableData {
                    name: "temperature".to_string(),
                    perm: None,
                    perm_size: 0,
                    values,
                }],
            }],
        }
    }

    #[test]
    fn sigma_is_worst_arm_distance() {
        let g = file(vec![1.0, 2.0, 4.0]);
        let arms = vec![
            ("fma".to_string(), Some(file(vec![1.0, 2.0, 4.0 + 4e-12]))),
            ("openblas".to_string(), Some(file(vec![1.0 + 4e-13, 2.0, 4.0]))),
        ];
        let e = measure_file("t", "case.result", &g, &arms);
        assert_eq!(e.len(), 1);
        assert_eq!(e[0].worst_arm, "fma");
        assert!((e[0].sigma_f.unwrap() - 1e-12).abs() < 1e-15);
        let (bar, rule) = e[0].bar().unwrap();
        assert_eq!(rule, "max10sigma");
        assert!((bar - 1e-11).abs() < 1e-14);
    }

    #[test]
    fn bitwise_agreement_takes_the_cap() {
        let g = file(vec![1.0, 2.0]);
        let arms = vec![("fma".to_string(), Some(file(vec![1.0, 2.0])))];
        let e = measure_file("t", "case.result", &g, &arms);
        assert_eq!(e[0].sigma_f, Some(0.0));
        assert_eq!(e[0].bar(), Some((1e-10, "cap-agree")));
    }

    #[test]
    fn missing_arm_file_is_structure_differs() {
        let g = file(vec![1.0, 2.0]);
        let arms = vec![("fma".to_string(), None)];
        let e = measure_file("t", "case.result", &g, &arms);
        assert_eq!(e[0].sigma_f, None);
        assert_eq!(e[0].bar(), None);
        let json = tolerances_json(&e, &serde_json::json!({}), "pin", "x");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["entries"][0]["rule"], "structure-differs");
        assert!(v["entries"][0]["bar"].is_null());
    }

    #[test]
    fn result_header_sniffing() {
        assert!(looks_like_result(b" BINARY 3.L\n"));
        assert!(looks_like_result(b" ASCII 3\n"));
        assert!(!looks_like_result(b"<?xml"));
    }
}
