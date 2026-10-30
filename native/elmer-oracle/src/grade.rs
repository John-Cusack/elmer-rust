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
    relative_linf_floor(golden, candidate, 0.0)
}

/// [`relative_linf`] normalised by `max(‖golden‖∞, floor)`. The floor is
/// the scale of the solver variable the field belongs to
/// ([`variable_scale`]): a component that is round-off around zero is
/// judged against its whole variable, as Elmer's own norm and Reference
/// Norm judge it, instead of against its own noise.
pub fn relative_linf_floor(golden: &[f64], candidate: &[f64], floor: f64) -> Option<f64> {
    if golden.len() != candidate.len() {
        return None;
    }
    let own = golden
        .iter()
        .filter(|v| v.is_finite())
        .fold(0.0f64, |m, &v| m.max(v.abs()));
    if own == 0.0 {
        // An identically zero golden field keeps the exact-zero rule: any
        // value beyond one ulp of its variable's scale is a real error, so
        // the scaled floor never applies here.
        let round_off = f64::EPSILON * floor;
        let ok = golden
            .iter()
            .zip(candidate.iter())
            .filter(|(g, _)| g.is_finite())
            .all(|(_, c)| c.abs() <= round_off);
        return ok.then_some(0.0);
    }
    let norm = own.max(floor);
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

/// Scale of the solver variable `field` belongs to at `step`: the largest
/// `‖·‖∞` over its component fields (from the file's preamble groups).
/// Zero for a field that is a variable of its own.
pub fn variable_scale(file: &crate::result::ResultFile, step: usize, field: &str) -> f64 {
    let Some(group) = file.groups.iter().find(|g| g.iter().any(|f| f == field)) else {
        return 0.0;
    };
    let Some(st) = file.steps.get(step) else {
        return 0.0;
    };
    let mut scale = 0.0f64;
    for member in group {
        if let Some(v) = st.vars.iter().find(|v| &v.name == member) {
            if let Ok(nodes) = to_node_space(v) {
                for x in nodes.iter().filter(|x| x.is_finite()) {
                    scale = scale.max(x.abs());
                }
            }
        }
    }
    scale
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
    grade_variable_floor(variable, step, golden, candidate, allowed, 0.0)
}

/// [`grade_variable`] with the near-zero floor of [`relative_linf_floor`].
pub fn grade_variable_floor(
    variable: &str,
    step: usize,
    golden: &VariableData,
    candidate: Option<&VariableData>,
    allowed: f64,
    floor: f64,
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
    match relative_linf_floor(&g, &c, floor) {
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
            let floor = variable_scale(golden, step, &g_var.name);
            out.extend(grade_variable_floor(&g_var.name, step, g_var, c_var, allowed, floor));
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
    fn near_zero_component_is_judged_against_its_variable() {
        // velocity 2 is round-off around zero; velocity 1 carries the flow.
        let mk = |name: &str, values: Vec<f64>| VariableData {
            name: name.to_string(),
            perm: None,
            perm_size: values.len() as i32,
            values,
        };
        let file = |v2: Vec<f64>| crate::result::ResultFile {
            steps: vec![crate::result::TimeStep {
                saves_done: 1,
                ntime: 1,
                sim_time: 0.0,
                vars: vec![mk("velocity 1", vec![1.0, 2.0]), mk("velocity 2", v2)],
            }],
            groups: vec![vec!["velocity 1".into(), "velocity 2".into()]],
        };
        let golden = file(vec![5e-16, -3e-16]);
        let noisy = file(vec![4e-16, -2e-16]);
        assert_eq!(variable_scale(&golden, 0, "velocity 2"), 2.0);
        assert!(grade_results(&golden, &noisy, 1e-10).is_empty());
        // Without the group the same noise is an O(1) relative error.
        let mut lone = golden.clone();
        lone.groups.clear();
        let mut lone_noisy = noisy.clone();
        lone_noisy.groups.clear();
        assert!(!grade_results(&lone, &lone_noisy, 1e-10).is_empty());
        // A real error in the small component is still caught.
        let wrong = file(vec![1e-3, -3e-16]);
        assert!(!grade_results(&golden, &wrong, 1e-10).is_empty());
        // An identically zero component keeps the exact-zero rule: one ulp
        // of the variable's scale passes, 1e-6 does not.
        let zero = file(vec![0.0, 0.0]);
        assert!(grade_results(&zero, &file(vec![2.0 * f64::EPSILON, 0.0]), 1e-10).is_empty());
        assert!(!grade_results(&zero, &file(vec![1e-6, 0.0]), 1e-10).is_empty());
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
            groups: Vec::new(),
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
        let empty = crate::result::ResultFile { steps: vec![], ..Default::default() };
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
