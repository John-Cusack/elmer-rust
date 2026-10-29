//! Tier bars derived from Fortran's measured sensitivity to legitimate
//! reordering (σ_F table, report §4.2).
//!
//! The bar for a field is `max(10·σ_F, cap)` where σ_F is the largest
//! relative-Linf distance between the golden Release arm and any
//! sensitivity arm (FMA contraction, OpenBLAS) on the same case, and
//! cap = 1e-10 applies only where every arm agrees bitwise. The table is
//! committed oracle data; it changes only by ADR with captain sign-off.

use crate::{OracleError, Result};

/// One measured field bar.
#[derive(Debug, Clone, PartialEq)]
pub struct TierEntry {
    pub test: String,
    pub field: String,
    pub step: usize,
    pub sigma_f: f64,
    pub bar: f64,
    pub rule: String,
}

/// Committed tolerance table.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ToleranceTable {
    entries: Vec<TierEntry>,
}

impl ToleranceTable {
    /// Parse the committed `tolerances.json` format.
    pub fn parse(text: &str) -> Result<ToleranceTable> {
        let v: serde_json::Value =
            serde_json::from_str(text).map_err(|e| OracleError::Parse(format!("tolerances: {e}")))?;
        let version = v
            .get("version")
            .and_then(|x| x.as_u64())
            .ok_or_else(|| OracleError::Parse(String::from("tolerances: missing version")))?;
        if version != 1 {
            return Err(OracleError::Parse(format!(
                "tolerances: unsupported version {version}"
            )));
        }
        let mut entries = Vec::new();
        let list = v
            .get("entries")
            .and_then(|x| x.as_array())
            .ok_or_else(|| OracleError::Parse(String::from("tolerances: missing entries")))?;
        for e in list {
            let get_str = |k: &str| {
                e.get(k)
                    .and_then(|x| x.as_str())
                    .ok_or_else(|| OracleError::Parse(format!("tolerances: entry missing {k}")))
            };
            let get_f64 = |k: &str| {
                e.get(k)
                    .and_then(|x| x.as_f64())
                    .ok_or_else(|| OracleError::Parse(format!("tolerances: entry missing {k}")))
            };
            entries.push(TierEntry {
                test: get_str("test")?.to_string(),
                field: get_str("field")?.to_string(),
                step: e
                    .get("step")
                    .and_then(|x| x.as_u64())
                    .ok_or_else(|| OracleError::Parse(String::from("tolerances: entry missing step")))?
                    as usize,
                sigma_f: get_f64("sigma_f")?,
                bar: get_f64("bar")?,
                rule: get_str("rule")?.to_string(),
            });
        }
        Ok(ToleranceTable { entries })
    }

    /// Tier bar for one field, if the table covers it.
    pub fn bar_for(&self, test: &str, field: &str, step: usize) -> Option<f64> {
        self.entries
            .iter()
            .find(|e| e.test == test && e.field == field && e.step == step)
            .map(|e| e.bar)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Recompute the Tier-A bar from a measured σ_F.
    pub fn tier_a_bar(sigma_f: f64) -> (f64, &'static str) {
        if sigma_f == 0.0 {
            (1e-10, "cap-agree")
        } else {
            (10.0 * sigma_f, "max10sigma")
        }
    }
}

/// Grade with per-field bars from the table, falling back to `fallback`
/// for uncovered fields (reported, never silent).
pub fn grade_with_table(
    golden: &crate::result::ResultFile,
    candidate: &crate::result::ResultFile,
    test: &str,
    table: &ToleranceTable,
    fallback: f64,
) -> (Vec<crate::grade::Mismatch>, Vec<String>) {
    use std::collections::HashMap;
    let mut out = Vec::new();
    let mut warnings = Vec::new();
    for (step, (g_step, c_step)) in golden.steps.iter().zip(candidate.steps.iter()).enumerate() {
        let c_vars: HashMap<&str, &crate::result::VariableData> =
            c_step.vars.iter().map(|v| (v.name.as_str(), v)).collect();
        for g_var in &g_step.vars {
            let bar = match table.bar_for(test, &g_var.name, step) {
                Some(b) => b,
                None => {
                    warnings.push(format!(
                        "no Tier bar for {test}/{} step {step}; using fallback {fallback:e}",
                        g_var.name
                    ));
                    fallback
                }
            };
            out.extend(crate::grade::grade_variable(
                &g_var.name,
                step,
                g_var,
                c_vars.get(g_var.name.as_str()).copied(),
                bar,
            ));
        }
    }
    if golden.steps.len() != candidate.steps.len() {
        out.push(crate::grade::Mismatch {
            variable: String::from("<steps>"),
            step: golden.steps.len().min(candidate.steps.len()),
            kind: crate::grade::MismatchKind::LengthDiffers {
                golden: golden.steps.len(),
                candidate: candidate.steps.len(),
            },
        });
    }
    (out, warnings)
}


#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"{
      "version": 1,
      "arms": {"release": "abc"},
      "entries": [
        {"test": "t", "field": "temperature", "step": 0,
         "sigma_f": 1.5e-8, "bar": 1.5e-7, "rule": "max10sigma"}
      ]
    }"#;

    #[test]
    fn parses_table_and_applies_bars() {
        let table = ToleranceTable::parse(SAMPLE).unwrap();
        assert_eq!(table.len(), 1);
        assert_eq!(table.bar_for("t", "temperature", 0), Some(1.5e-7));
        assert_eq!(table.bar_for("t", "other", 0), None);
        assert_eq!(table.bar_for("t", "temperature", 1), None);
    }

    #[test]
    fn rejects_bad_version_and_shape() {
        assert!(ToleranceTable::parse(r#"{"version": 2, "entries": []}"#).is_err());
        assert!(ToleranceTable::parse(r#"{"version": 1}"#).is_err());
        assert!(ToleranceTable::parse("not json").is_err());
    }

    #[test]
    fn tier_a_bar_rules() {
        assert_eq!(ToleranceTable::tier_a_bar(0.0), (1e-10, "cap-agree"));
        assert_eq!(ToleranceTable::tier_a_bar(1.5e-8), (1.5e-7, "max10sigma"));
    }
}
