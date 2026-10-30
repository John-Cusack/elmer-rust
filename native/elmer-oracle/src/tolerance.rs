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
    /// Result file relative to the run dir; `case.result` when the table
    /// predates per-file keys.
    pub file: String,
    pub field: String,
    /// Which copy of `field` in the step (0 unless the name repeats).
    pub occurrence: usize,
    pub step: usize,
    /// `None` (with `bar` `None`) when the arms disagree structurally.
    pub sigma_f: Option<f64>,
    /// `None` means no Tier-A bar exists: the field is graded by O1/O6.
    pub bar: Option<f64>,
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
            // Null is meaningful (structure-differs); absence is an error.
            let get_opt_f64 = |k: &str| match e.get(k) {
                Some(serde_json::Value::Null) => Ok(None),
                Some(x) => x
                    .as_f64()
                    .map(Some)
                    .ok_or_else(|| OracleError::Parse(format!("tolerances: {k} not a number"))),
                None => Err(OracleError::Parse(format!("tolerances: entry missing {k}"))),
            };
            entries.push(TierEntry {
                test: get_str("test")?.to_string(),
                file: e
                    .get("file")
                    .and_then(|x| x.as_str())
                    .unwrap_or("case.result")
                    .to_string(),
                field: get_str("field")?.to_string(),
                occurrence: e.get("occurrence").and_then(|x| x.as_u64()).unwrap_or(0) as usize,
                step: e
                    .get("step")
                    .and_then(|x| x.as_u64())
                    .ok_or_else(|| OracleError::Parse(String::from("tolerances: entry missing step")))?
                    as usize,
                sigma_f: get_opt_f64("sigma_f")?,
                bar: get_opt_f64("bar")?,
                rule: get_str("rule")?.to_string(),
            });
        }
        Ok(ToleranceTable { entries })
    }

    /// Table entry for one field of one result file, if covered.
    pub fn entry_for(
        &self,
        test: &str,
        file: &str,
        field: &str,
        occurrence: usize,
        step: usize,
    ) -> Option<&TierEntry> {
        self.entries.iter().find(|e| {
            e.test == test
                && e.file == file
                && e.field == field
                && e.occurrence == occurrence
                && e.step == step
        })
    }

    /// Tier bar for one field of `case.result`, if the table covers it.
    pub fn bar_for(&self, test: &str, field: &str, step: usize) -> Option<f64> {
        self.entry_for(test, "case.result", field, 0, step).and_then(|e| e.bar)
    }

    pub fn entries(&self) -> &[TierEntry] {
        &self.entries
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
    file: &str,
    table: &ToleranceTable,
    fallback: f64,
) -> (Vec<crate::grade::Mismatch>, Vec<String>) {
    let mut out = Vec::new();
    let mut warnings = Vec::new();
    for (step, (g_step, c_step)) in golden.steps.iter().zip(candidate.steps.iter()).enumerate() {
        for (k, (g_var, c_var)) in
            crate::grade::pair_vars(&g_step.vars, &c_step.vars).into_iter().enumerate()
        {
            let occurrence = g_step.vars[..k].iter().filter(|v| v.name == g_var.name).count();
            let bar = match table.entry_for(test, file, &g_var.name, occurrence, step) {
                Some(TierEntry { bar: Some(b), .. }) => *b,
                Some(TierEntry { rule, .. }) => {
                    warnings.push(format!(
                        "{test}/{file}/{} step {step}: no Tier-A bar ({rule}); graded by O1/O6, not here",
                        g_var.name
                    ));
                    continue;
                }
                None => {
                    warnings.push(format!(
                        "no Tier bar for {test}/{} step {step}; using fallback {fallback:e}",
                        g_var.name
                    ));
                    fallback
                }
            };
            let floor = crate::grade::variable_scale(golden, step, &g_var.name);
            out.extend(crate::grade::grade_variable_floor(
                &g_var.name,
                step,
                g_var,
                c_var,
                bar,
                floor,
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
        assert_eq!(table.entries()[0].file, "case.result");
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
