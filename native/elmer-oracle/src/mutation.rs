//! Mutation suite over golden result files (report §4.3 item 3).
//!
//! Each mutation perturbs a copy of a golden file the way a wrong
//! candidate would, and the grader, with the committed Tier bars, must
//! reject every copy it can see. Perturbations act in node space and are
//! written back through the Perm, so a mutated file is always internally
//! consistent and is caught for the right reason.
//!
//! Field mutations (noise, permutation, one value) are applied to one
//! field of one step at a time. A mutation is *detectable* when its size
//! in the grader's own metric (relative Linf against the golden field)
//! exceeds that field's Tier bar; detectable mutations must be caught, and
//! the rest are reported as below the floor, which is the measured price
//! of a field whose Fortran reorder noise (10·σ_F) is already larger than
//! the perturbation. A mutation that cannot change the file at all (a
//! permutation of a constant field) is vacuous. File mutations (a dropped
//! step, swapped variables) change structure and are always detectable.

use crate::grade::{relative_linf_floor, to_node_space, variable_scale};
use crate::result::{ResultFile, VariableData};
use crate::tolerance::{grade_with_table, TierEntry, ToleranceTable};

/// The §4.3 perturbations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mutation {
    /// Every value of the field moved by 1e-6·‖u‖∞, alternating sign.
    Noise,
    /// The field's largest and smallest nodal values swap places.
    Permute,
    /// One nodal value (the field's largest) moved by 1e-6·‖u‖∞.
    OneValue,
    /// The last saved step removed.
    DropStep,
    /// Two fields of one step exchange values.
    SwapVariables,
}

impl Mutation {
    pub fn name(self) -> &'static str {
        match self {
            Mutation::Noise => "noise-1e-6",
            Mutation::Permute => "node-permutation",
            Mutation::OneValue => "one-value",
            Mutation::DropStep => "dropped-step",
            Mutation::SwapVariables => "swapped-variable",
        }
    }
}

/// Outcome of one mutation.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Caught,
    /// Detectable, yet the grader accepted it: a grader bug.
    Escaped,
    /// Its size in the grader metric is within the field's bar.
    BelowFloor,
    /// The mutation cannot change this file.
    Vacuous(&'static str),
}

/// One mutation of one golden file.
#[derive(Debug, Clone, PartialEq)]
pub struct Trial {
    pub mutation: Mutation,
    /// `field@step` for field mutations, empty for file mutations.
    pub target: String,
    pub outcome: Outcome,
}

/// Write node-space values back into Perm order.
fn from_node_space(var: &VariableData, nodes: &[f64]) -> VariableData {
    let mut out = var.clone();
    match &var.perm {
        None => out.values = nodes.to_vec(),
        Some(pairs) => {
            for (i, v) in out.values.iter_mut().enumerate() {
                if let Some(&(_, k)) = pairs.get(i) {
                    if k > 0 && (k as usize) <= nodes.len() {
                        *v = nodes[k as usize - 1];
                    }
                }
            }
        }
    }
    out
}

fn scale(nodes: &[f64]) -> f64 {
    let m = nodes.iter().filter(|v| v.is_finite()).fold(0.0f64, |m, v| m.max(v.abs()));
    if m > 0.0 {
        m
    } else {
        1.0
    }
}

/// Mutate one field's node-space values; `Err` names why it cannot change.
pub fn mutate_nodes(nodes: &[f64], m: Mutation) -> Result<Vec<f64>, &'static str> {
    let s = scale(nodes);
    let finite: Vec<usize> = (0..nodes.len()).filter(|&i| nodes[i].is_finite()).collect();
    let mut out = nodes.to_vec();
    match m {
        Mutation::Noise => {
            for (i, v) in out.iter_mut().enumerate() {
                *v += if i % 2 == 0 { 1e-6 } else { -1e-6 } * s;
            }
        }
        Mutation::Permute => {
            let lo = finite.iter().copied().min_by(|&a, &b| nodes[a].total_cmp(&nodes[b]));
            let hi = finite.iter().copied().max_by(|&a, &b| nodes[a].total_cmp(&nodes[b]));
            match (lo, hi) {
                (Some(a), Some(b)) if nodes[a] != nodes[b] => out.swap(a, b),
                _ => return Err("constant field"),
            }
        }
        Mutation::OneValue => {
            let at = finite
                .iter()
                .copied()
                .max_by(|&a, &b| nodes[a].abs().total_cmp(&nodes[b].abs()))
                .ok_or("empty field")?;
            out[at] += 1e-6 * s;
        }
        Mutation::DropStep | Mutation::SwapVariables => return Err("not a field mutation"),
    }
    Ok(out)
}

fn bar_of(table: &ToleranceTable, test: &str, file: &str, field: &str, occ: usize, step: usize, fallback: f64) -> Option<f64> {
    match table.entry_for(test, file, field, occ, step) {
        Some(TierEntry { bar, .. }) => *bar,
        None => Some(fallback),
    }
}

fn grade(golden: &ResultFile, mutated: &ResultFile, test: &str, file: &str, table: &ToleranceTable, fallback: f64) -> Outcome {
    let (mismatches, _) = grade_with_table(golden, mutated, test, file, table, fallback);
    if mismatches.is_empty() {
        Outcome::Escaped
    } else {
        Outcome::Caught
    }
}

/// Run every mutation against one golden file with the committed bars.
pub fn run_suite(
    golden: &ResultFile,
    test: &str,
    file: &str,
    table: &ToleranceTable,
    fallback: f64,
) -> Vec<Trial> {
    let mut trials = Vec::new();
    for (si, step) in golden.steps.iter().enumerate() {
        for (vi, var) in step.vars.iter().enumerate() {
            let occ = step.vars[..vi].iter().filter(|v| v.name == var.name).count();
            let target = format!("{}@{si}", var.name);
            let nodes = to_node_space(var).ok();
            let bar = bar_of(table, test, file, &var.name, occ, si, fallback);
            for m in [Mutation::Noise, Mutation::Permute, Mutation::OneValue] {
                let outcome = match (&nodes, bar) {
                    (None, _) => Outcome::Vacuous("unmappable field"),
                    (_, None) => Outcome::Vacuous("no Tier-A bar (structure-differs)"),
                    (Some(nodes), Some(bar)) => match mutate_nodes(nodes, m) {
                        Err(why) => Outcome::Vacuous(why),
                        Ok(mutated_nodes) => {
                            let floor = variable_scale(golden, si, &var.name);
                            let size = relative_linf_floor(nodes, &mutated_nodes, floor)
                                .unwrap_or(f64::INFINITY);
                            if size <= bar {
                                Outcome::BelowFloor
                            } else {
                                let mut mutated = golden.clone();
                                mutated.steps[si].vars[vi] = from_node_space(var, &mutated_nodes);
                                grade(golden, &mutated, test, file, table, fallback)
                            }
                        }
                    },
                };
                trials.push(Trial { mutation: m, target: target.clone(), outcome });
            }
        }
    }
    // Structure mutations.
    let mut dropped = golden.clone();
    dropped.steps.pop();
    trials.push(Trial {
        mutation: Mutation::DropStep,
        target: String::new(),
        outcome: grade(golden, &dropped, test, file, table, fallback),
    });
    let outcome = match swap_first_distinct(golden) {
        Err(why) => Outcome::Vacuous(why),
        Ok((s, si, a, b)) => {
            // Detectable when either field, in the other's place, moves
            // further than its own bar (or the lengths differ).
            let (va, vb) = (&golden.steps[si].vars[a], &golden.steps[si].vars[b]);
            let size = match (to_node_space(va), to_node_space(vb)) {
                (Ok(x), Ok(y)) => {
                    let fa = variable_scale(golden, si, &va.name);
                    let fb = variable_scale(golden, si, &vb.name);
                    relative_linf_floor(&x, &y, fa)
                        .unwrap_or(f64::INFINITY)
                        .max(relative_linf_floor(&y, &x, fb).unwrap_or(f64::INFINITY))
                }
                _ => f64::INFINITY,
            };
            let occ = |vi: usize| {
                let n = &golden.steps[si].vars[vi].name;
                golden.steps[si].vars[..vi].iter().filter(|v| &v.name == n).count()
            };
            let bar_a = bar_of(table, test, file, &va.name, occ(a), si, fallback);
            let bar_b = bar_of(table, test, file, &vb.name, occ(b), si, fallback);
            let bar = bar_a.unwrap_or(f64::INFINITY).min(bar_b.unwrap_or(f64::INFINITY));
            if size <= bar {
                Outcome::BelowFloor
            } else {
                grade(golden, &s, test, file, table, fallback)
            }
        }
    };
    trials.push(Trial { mutation: Mutation::SwapVariables, target: String::new(), outcome });
    trials
}

/// Node-space equality where unmapped dofs (NaN) match each other.
fn same_field(a: &VariableData, b: &VariableData) -> bool {
    match (to_node_space(a), to_node_space(b)) {
        (Ok(x), Ok(y)) => {
            x.len() == y.len()
                && x.iter().zip(&y).all(|(p, q)| p == q || (p.is_nan() && q.is_nan()))
        }
        _ => false,
    }
}

/// Exchange the values of the first two fields of a step that differ.
/// Returns the mutated file plus (step, field a, field b).
fn swap_first_distinct(
    golden: &ResultFile,
) -> Result<(ResultFile, usize, usize, usize), &'static str> {
    let mut out = golden.clone();
    for (si, step) in out.steps.iter_mut().enumerate() {
        let n = step.vars.len();
        for a in 0..n {
            for b in a + 1..n {
                if !same_field(&step.vars[a], &step.vars[b]) {
                    let (va, vb) = (step.vars[a].clone(), step.vars[b].clone());
                    step.vars[a] = VariableData { name: va.name, ..vb.clone() };
                    step.vars[b] = VariableData { name: vb.name, ..va };
                    return Ok((out, si, a, b));
                }
            }
        }
    }
    Err("fewer than two distinct fields")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result::TimeStep;

    fn var(name: &str, values: Vec<f64>) -> VariableData {
        let n = values.len() as i32;
        VariableData {
            name: name.to_string(),
            // Stored order reversed relative to dof order.
            perm: Some((1..=n).map(|i| (i, n + 1 - i)).collect()),
            perm_size: n,
            values,
        }
    }

    fn golden() -> ResultFile {
        let step = |t: f64| TimeStep {
            saves_done: 1,
            ntime: 1,
            sim_time: t,
            vars: vec![
                var("temperature", vec![300.0, 310.0, 295.0]),
                var("flux", vec![1.0, 2.0, 3.0]),
            ],
        };
        ResultFile { steps: vec![step(0.0), step(1.0)], ..Default::default() }
    }

    #[test]
    fn every_mutation_is_caught_at_a_tight_bar() {
        let table = ToleranceTable::default();
        let trials = run_suite(&golden(), "t", "case.result", &table, 1e-10);
        assert_eq!(trials.len(), 2 * 2 * 3 + 2);
        for t in trials {
            assert_eq!(t.outcome, Outcome::Caught, "{t:?}");
        }
    }

    #[test]
    fn mutations_within_a_loose_bar_are_below_floor_not_escapes() {
        let table = ToleranceTable::default();
        for t in run_suite(&golden(), "t", "case.result", &table, 1e-5) {
            match t.mutation {
                Mutation::Noise | Mutation::OneValue => assert_eq!(t.outcome, Outcome::BelowFloor),
                _ => assert_eq!(t.outcome, Outcome::Caught, "{t:?}"),
            }
        }
    }

    #[test]
    fn mutated_fields_stay_perm_consistent() {
        let g = golden();
        let v = &g.steps[0].vars[0];
        let nodes = to_node_space(v).unwrap();
        for m in [Mutation::Noise, Mutation::Permute, Mutation::OneValue] {
            let back = from_node_space(v, &mutate_nodes(&nodes, m).unwrap());
            assert_eq!(to_node_space(&back).unwrap(), mutate_nodes(&nodes, m).unwrap());
        }
    }

    #[test]
    fn constant_fields_make_permutation_vacuous() {
        assert!(mutate_nodes(&[5.0, 5.0, 5.0], Mutation::Permute).is_err());
        let mut g = golden();
        for s in &mut g.steps {
            s.vars = vec![var("c", vec![5.0, 5.0, 5.0])];
        }
        assert!(swap_first_distinct(&g).is_err());
        // Unmapped dofs are NaN in node space; identical fields with the
        // same gaps are still the same field.
        let gap = |name: &str| VariableData {
            name: name.to_string(),
            perm: Some(vec![(1, 1), (2, 3)]),
            perm_size: 2,
            values: vec![1.0, 2.0],
        };
        assert!(same_field(&gap("a"), &gap("b")));
    }
}
