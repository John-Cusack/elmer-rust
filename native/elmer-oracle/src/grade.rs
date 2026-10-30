//! External grading of candidate output against frozen golden fields.
//!
//! The grader compares full-precision captured fields, never the
//! self-graded `TEST.PASSED`. Errors are measured in node space: values
//! stored in Perm order are mapped back through the Perm pairs before
//! comparison, and the per-field bar is `max(10·σ_F, cap)` from the
//! tolerance table (Tier A). This module owns the comparison mechanics;
//! tolerance values themselves are oracle data, loaded from the committed
//! table.

use crate::result::VariableData;

/// One field-level mismatch found by the grader.
#[derive(Debug, Clone, PartialEq)]
pub struct Mismatch {
    pub variable: String,
    pub step: usize,
    pub kind: MismatchKind,
}

/// How a field comparison failed.
#[derive(Debug, Clone, PartialEq)]
pub enum MismatchKind {
    MissingVariable,
    MissingStep,
    LengthDiffers { golden: usize, candidate: usize },
    /// Relative Linf error above the allowed bar.
    ExceedsTolerance { error: f64, allowed: f64 },
    /// A non-finite value where the golden field is finite.
    NonFinite { index: usize },
    /// The file itself is inconsistent (e.g. one dof, two values).
    CorruptFile { detail: String },
}

/// Map Perm-ordered values back to dof space. `perm` holds 1-based
/// `(node, dof)` pairs in stored order: the writer emits `Values(Perm(i))`
/// for each positive entry, so the dof index is the pair's second
/// element. Duplicate nodes share one dof and must agree bitwise;
/// disagreement means a corrupt file. Proved against the pristine
/// Release arm: this mapping plus sequential RMS reproduces the
/// solver-internal norm exactly (Tolerance-0 rerun, rel err 0).
pub fn to_node_space(var: &VariableData) -> crate::Result<Vec<f64>> {
    match &var.perm {
        None => Ok(var.values.clone()),
        Some(pairs) => {
            let max_dof = pairs.iter().map(|&(_, k)| k).max().unwrap_or(0).max(0) as usize;
            let mut dof = vec![f64::NAN; max_dof];
            for (i, v) in var.values.iter().enumerate() {
                let Some(&(_, k)) = pairs.get(i) else {
                    continue;
                };
                if k <= 0 || (k as usize) > max_dof {
                    continue;
                }
                let slot = &mut dof[k as usize - 1];
                if slot.is_nan() {
                    *slot = *v;
                } else if *slot != *v {
                    return Err(crate::OracleError::GradeFailed(format!(
                        "variable {}: dof {k} has conflicting values",
                        var.name
                    )));
                }
            }
            Ok(dof)
        }
    }
}

/// Relative Linf error of `candidate` against `golden`, normalized by
/// ‖golden‖∞. Returns `None` when the golden field is all zeros (exact
/// zero then requires exact zero; any nonzero is reported as infinite).
pub fn relative_linf(golden: &[f64], candidate: &[f64]) -> Option<f64> {
    if golden.len() != candidate.len() {
        return None;
    }
    let norm = golden.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    if norm == 0.0 {
        if candidate.iter().all(|&v| v == 0.0) {
            return Some(0.0);
        }
        return None;
    }
    let mut worst = 0.0;
    for (g, c) in golden.iter().zip(candidate.iter()) {
        if !c.is_finite() && g.is_finite() {
            return None;
        }
        let e = (c - g).abs() / norm;
        if e > worst {
            worst = e;
        }
    }
    Some(worst)
}

/// Grade one variable at one step. `allowed` is the Tier bar for this
/// field and step.
pub fn grade_variable(
    variable: &str,
    step: usize,
    golden: &VariableData,
    candidate: Option<&VariableData>,
    allowed: f64,
) -> Vec<Mismatch> {
    let mk = |kind| vec![Mismatch {
        variable: variable.to_string(),
        step,
        kind,
    }];
    let cand = match candidate {
        None => return mk(MismatchKind::MissingVariable),
        Some(c) => c,
    };
    let g = match to_node_space(golden) {
        Ok(g) => g,
        Err(e) => return mk(MismatchKind::CorruptFile { detail: e.to_string() }),
    };
    let c = match to_node_space(cand) {
        Ok(c) => c,
        Err(e) => return mk(MismatchKind::CorruptFile { detail: e.to_string() }),
    };
    if g.len() != c.len() {
        return mk(MismatchKind::LengthDiffers {
            golden: g.len(),
            candidate: c.len(),
        });
    }
    for (i, (&gv, &cv)) in g.iter().zip(c.iter()).enumerate() {
        if gv.is_finite() && !cv.is_finite() {
            return mk(MismatchKind::NonFinite { index: i });
        }
    }
    match relative_linf(&g, &c) {
        None => mk(MismatchKind::ExceedsTolerance {
            error: f64::INFINITY,
            allowed,
        }),
        Some(e) if e > allowed => mk(MismatchKind::ExceedsTolerance { error: e, allowed }),
        _ => Vec::new(),
    }
}

/// Pair each golden variable of one step with the candidate variable of
/// the same name and the same occurrence: a step can carry one name twice
/// (a solver-specific mesh writing its own copy), and the k-th golden
/// copy grades against the k-th candidate copy, never the first.
pub fn pair_vars<'a>(
    golden: &'a [VariableData],
    candidate: &'a [VariableData],
) -> Vec<(&'a VariableData, Option<&'a VariableData>)> {
    golden
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let k = golden[..i].iter().filter(|v| v.name == g.name).count();
            (g, candidate.iter().filter(|c| c.name == g.name).nth(k))
        })
        .collect()
}

/// Grade every saved step of a candidate `.result` file against a golden
/// one: variables matched by exact name, steps by index. `allowed` is
/// the Tier bar applied uniformly; per-field tables refine it later.
pub fn grade_results(
    golden: &crate::result::ResultFile,
    candidate: &crate::result::ResultFile,
    allowed: f64,
) -> Vec<Mismatch> {
    let mut out = Vec::new();
    if golden.steps.len() != candidate.steps.len() {
        out.push(Mismatch {
            variable: String::from("<steps>"),
            step: golden.steps.len().min(candidate.steps.len()),
            kind: MismatchKind::LengthDiffers {
                golden: golden.steps.len(),
                candidate: candidate.steps.len(),
            },
        });
    }
    for (step, (g_step, c_step)) in golden.steps.iter().zip(candidate.steps.iter()).enumerate() {
        for (g_var, c_var) in pair_vars(&g_step.vars, &c_step.vars) {
            out.extend(grade_variable(&g_var.name, step, g_var, c_var, allowed));
        }
    }
    out
}
/// Perturbations for the mutation suite (§4.3 item 3): every mutated golden
/// copy must fail grading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mutation {
    /// Add 1e-6 relative noise to every value.
    Noise,
    /// Reverse the value order (node permutation).
    Reverse,
    /// Change one boundary value by 10%.
    OneValue,
    /// Drop the last value (length change).
    DropLast,
}

/// Apply a mutation to node-space values.
pub fn mutate(values: &[f64], mutation: Mutation) -> Vec<f64> {
    match mutation {
        Mutation::Noise => values.iter().map(|&v| v * (1.0 + 1e-6)).collect(),
        Mutation::Reverse => values.iter().rev().copied().collect(),
        Mutation::OneValue => {
            let mut out = values.to_vec();
            if let Some(first) = out.first_mut() {
                *first *= 1.1;
            }
            out
        }
        Mutation::DropLast => values[..values.len().saturating_sub(1)].to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn var(values: Vec<f64>) -> VariableData {
        VariableData {
            name: "t".to_string(),
            perm: None,
            perm_size: values.len() as i32,
            values,
        }
    }

    #[test]
    fn identical_fields_pass() {
        let g = var(vec![1.0, 2.0, 3.0]);
        let c = var(vec![1.0, 2.0, 3.0]);
        assert!(grade_variable("t", 0, &g, Some(&c), 1e-10).is_empty());
    }

    #[test]
    fn each_mutation_is_caught() {
        let golden = var(vec![300.0, 310.0, 295.0, 305.0]);
        for m in [
            Mutation::Noise,
            Mutation::Reverse,
            Mutation::OneValue,
            Mutation::DropLast,
        ] {
            let mutated = var(mutate(&golden.values, m));
            let mismatches = grade_variable("t", 0, &golden, Some(&mutated), 1e-10);
            assert!(!mismatches.is_empty(), "mutation {m:?} escaped grading");
        }
    }

    #[test]
    fn small_noise_within_loose_bar_passes() {
        let golden = var(vec![300.0, 310.0]);
        let close = var(vec![300.0 + 1e-9, 310.0]);
        assert!(grade_variable("t", 0, &golden, Some(&close), 1e-6).is_empty());
    }

    #[test]
    fn missing_and_short_candidates_fail() {
        let golden = var(vec![1.0, 2.0]);
        assert_eq!(
            grade_variable("t", 1, &golden, None, 1e-10)[0].kind,
            MismatchKind::MissingVariable
        );
        let short = var(vec![1.0]);
        assert!(matches!(
            grade_variable("t", 1, &golden, Some(&short), 1e-10)[0].kind,
            MismatchKind::LengthDiffers { .. }
        ));
    }

    #[test]
    fn perm_order_maps_to_dof_space() {
        // Stored order follows positive Perm entries; the dof index is
        // the pair's second element: Values(Perm(i)).
        let v = VariableData {
            name: "t".to_string(),
            perm: Some(vec![(3, 2), (1, 1)]),
            perm_size: 3,
            values: vec![30.0, 10.0],
        };
        assert_eq!(to_node_space(&v).unwrap(), vec![10.0, 30.0]);
    }

    #[test]
    fn conflicting_duplicate_dofs_are_corrupt() {
        let v = VariableData {
            name: "t".to_string(),
            perm: Some(vec![(1, 1), (2, 1)]),
            perm_size: 2,
            values: vec![10.0, 11.0],
        };
        assert!(to_node_space(&v).is_err());
        assert!(matches!(
            grade_variable("t", 0, &v, Some(&v), 1e-10)[0].kind,
            MismatchKind::CorruptFile { .. }
        ));
    }

    fn result_file(vars: Vec<VariableData>) -> crate::result::ResultFile {
        crate::result::ResultFile {
            steps: vec![crate::result::TimeStep {
                saves_done: 1,
                ntime: 1,
                sim_time: 1.0,
                vars,
            }],
        }
    }

    #[test]
    fn grade_results_matches_by_name_and_step() {
        let golden = result_file(vec![var(vec![1.0, 2.0])]);
        let same = result_file(vec![var(vec![1.0, 2.0])]);
        assert!(grade_results(&golden, &same, 1e-10).is_empty());
        let mut drifted = var(vec![1.0, 2.0 + 1e-3]);
        drifted.name = "t".to_string();
        let bad = result_file(vec![drifted]);
        let mm = grade_results(&golden, &bad, 1e-10);
        assert_eq!(mm.len(), 1);
        assert!(matches!(
            mm[0].kind,
            MismatchKind::ExceedsTolerance { .. }
        ));
    }

    #[test]
    fn grade_results_flags_missing_variables_and_steps() {
        let golden = result_file(vec![var(vec![1.0])]);
        let empty = crate::result::ResultFile { steps: vec![] };
        let mm = grade_results(&golden, &empty, 1e-10);
        assert!(!mm.is_empty());
        let other = result_file(vec![VariableData {
            name: "u".to_string(),
            perm: None,
            perm_size: 1,
            values: vec![1.0],
        }]);
        let mm = grade_results(&golden, &other, 1e-10);
        assert!(matches!(mm[0].kind, MismatchKind::MissingVariable));
    }
}
