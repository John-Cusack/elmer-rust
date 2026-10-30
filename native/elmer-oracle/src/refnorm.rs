//! Bitwise grader validation (report §4.3 item 2, §5.2 item 4).
//!
//! For each solver a sif checks with `Reference Norm`, the grader
//! recomputes the norm from the captured golden field and writes it back
//! as the reference with `Reference Norm Tolerance = 0`. Elmer's own
//! `CompareToReferenceSolution` then passes only on exact equality, which
//! proves both the capture (the `.result` holds the solver's exact values,
//! in an order the grader can map) and the grader's `ComputeNorm`.
//!
//! The eligible class is deliberately narrow and is reported as such:
//! solvers with an explicit scalar `Variable`, default norm options
//! (no `Nonlinear System Norm *`, no `Norm Permutation`), set inside a
//! `Solver N` section of a single-run sif.

use crate::norm::{compute_norm, NormOpts};
use crate::result::ResultFile;
use crate::runner::{keyword_of, opens_table, strip_sif_comment};

/// One `Reference Norm` solver of a sif.
#[derive(Debug, Clone, PartialEq)]
pub struct RefSolver {
    pub number: u32,
    /// Lowercased variable name, when the solver names a scalar one.
    pub variable: Result<String, String>,
}

/// Keywords that change `ComputeNorm`; their presence takes a solver out
/// of the validated class.
const NORM_KEYS: [&str; 7] = [
    // SaveScalars and friends report a chosen scalar as the solver norm.
    "show norm index",
    "show norm",
    "nonlinear system norm degree",
    "nonlinear system norm dofs",
    "nonlinear system norm components",
    "nonlinear system nodal norm",
    "norm permutation",
];

/// The `Reference Norm` solvers of one sif, or why the sif is out of the
/// validated class as a whole. Keys come from `Solver N` sections and from
/// `Solver N :: key = value` one-liners (the usual home of `Reference
/// Norm`), merged in file order.
pub fn reference_solvers(text: &str) -> Result<Vec<RefSolver>, String> {
    let mut keys: std::collections::BTreeMap<u32, Vec<(String, String)>> = Default::default();
    let solver_number = |sec: &str| -> Option<u32> { sec.trim().strip_prefix("solver")?.trim().parse().ok() };
    let mut runs = 0;
    let mut seen = std::collections::BTreeSet::new();
    let mut section: Option<String> = None;
    let mut block: Vec<(String, String)> = Vec::new();
    let mut in_table = false;
    let mut in_comment = false;
    for raw in text.lines() {
        let code = strip_sif_comment(raw).trim();
        let lower = code.to_lowercase();
        if in_comment {
            in_comment = lower != "*/";
            continue;
        }
        if in_table {
            in_table = !lower.starts_with("end");
            continue;
        }
        if let Some(sec) = &section {
            if lower == "end" {
                if let Some(n) = solver_number(sec) {
                    if !seen.insert(n) {
                        return Err(format!("two Solver {n} sections (merged by Elmer)"));
                    }
                    keys.entry(n).or_default().append(&mut block);
                }
                block.clear();
                section = None;
            } else {
                if let (Some(k), Some((_, v))) = (keyword_of(code), code.split_once('=')) {
                    block.push((k, v.trim().to_string()));
                }
                in_table = opens_table(code);
            }
            continue;
        }
        if lower.is_empty() || lower.starts_with('$') {
            continue;
        }
        if let Some((sec, rest)) = code.split_once("::") {
            if let (Some(n), Some(k), Some((_, v))) =
                (solver_number(&sec.to_lowercase()), keyword_of(rest), rest.split_once('='))
            {
                keys.entry(n).or_default().push((k, v.trim().to_string()));
            }
            continue;
        }
        if lower == "run" {
            runs += 1;
            continue;
        }
        if lower.starts_with("include") {
            return Err(String::from("sif uses Include"));
        }
        if lower.starts_with("echo ") || lower.starts_with("check ") || lower == "stop" {
            continue;
        }
        if lower == "/*" {
            in_comment = true;
            continue;
        }
        section = Some(lower);
    }
    if runs > 1 {
        return Err(format!("{runs} RUN blocks"));
    }
    Ok(keys
        .into_iter()
        .filter(|(_, k)| k.iter().any(|(k, _)| k == "reference norm"))
        .map(|(number, k)| RefSolver { number, variable: solver_variable(&k) })
        .collect())
}

fn solver_variable(keys: &[(String, String)]) -> Result<String, String> {
    if let Some((k, _)) = keys.iter().find(|(k, _)| NORM_KEYS.contains(&k.as_str())) {
        return Err(format!("norm option {k}"));
    }
    let Some((_, value)) = keys.iter().rev().find(|(k, _)| k == "variable") else {
        return Err(String::from("no explicit Variable (legacy or default name)"));
    };
    let v = value.trim().trim_start_matches("-global").trim();
    let v = v
        .strip_prefix("String")
        .or_else(|| v.strip_prefix("string"))
        .unwrap_or(v)
        .trim()
        .trim_matches('"')
        .trim();
    if v.contains('[') {
        return Err(String::from("vector variable"));
    }
    let lower = v.to_lowercase();
    let mut parts = lower.split_whitespace().peekable();
    let mut name = Vec::new();
    while let Some(p) = parts.next() {
        if p == "-dofs" {
            match parts.next().and_then(|n| n.parse::<u32>().ok()) {
                Some(1) => continue,
                _ => return Err(String::from("multi-dof variable")),
            }
        }
        if p.starts_with('-') {
            return Err(format!("variable option {p}"));
        }
        name.push(p);
    }
    if name.is_empty() {
        return Err(String::from("empty Variable"));
    }
    Ok(name.join(" "))
}

/// Whether `variable`'s Perm reaches past the mesh nodes (edge, face or
/// DG unknowns), where the grader's nodal reconstruction may not follow
/// the solver's own norm path.
pub fn non_nodal(file: &ResultFile, nodes: Option<usize>, variable: &str) -> bool {
    file.steps.last().is_some_and(|s| {
        s.vars.iter().any(|v| {
            v.name.to_lowercase() == variable && nodes.is_none_or(|n| v.perm_size as usize > n)
        })
    })
}

/// Norm of `variable` at the last saved step of `file`, computed as the
/// serial `ComputeNorm` does from the rebuilt `Values` and `Perm` arrays.
pub fn golden_norm(file: &ResultFile, sizes: &[(String, usize)], variable: &str) -> Result<f64, String> {
    let step = file.steps.last().ok_or("no saved step")?;
    let mut matches = step.vars.iter().filter(|v| v.name.to_lowercase() == variable);
    let var = matches.next().ok_or_else(|| format!("{variable} not in the result file"))?;
    if matches.next().is_some() {
        return Err(format!("{variable} appears twice"));
    }
    let pairs = var.perm.as_ref().ok_or("NULL Perm")?;

    let n_values = pairs.iter().map(|&(_, k)| k).max().unwrap_or(0).max(0) as usize;
    let mut values = vec![f64::NAN; n_values];
    let mut perm = vec![0i32; var.perm_size.max(0) as usize];
    for (i, &(node, k)) in pairs.iter().enumerate() {
        let (Some(slot), Some(v)) = (perm.get_mut(node as usize - 1), var.values.get(i)) else {
            return Err(String::from("Perm entry out of range"));
        };
        *slot = k;
        values[k as usize - 1] = *v;
    }
    if values.iter().any(|v| v.is_nan()) {
        return Err(String::from("values not covered by the Perm (not all saved)"));
    }
    // Edge, face, bubble or DG unknowns beyond the nodal Perm are not in
    // the file, yet ComputeNorm sums them: out of class.
    match sizes.iter().find(|(n, _)| n.to_lowercase() == variable) {
        Some(&(_, size)) if size == n_values => {}
        Some(&(_, size)) => {
            return Err(format!("{size} values, {n_values} saved (non-nodal unknowns)"));
        }
        None => return Err(format!("{variable} not in the preamble")),
    }
    compute_norm(&values, Some(&perm), 1, perm.len(), &NormOpts::default())
        .map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result::{TimeStep, VariableData};

    #[test]
    fn finds_reference_solvers_and_their_variables() {
        let sif = "Header\n Mesh DB \".\" \"m\"\nEnd\nSolver 1\n  Equation = Heat\n  Variable = Temperature\n  Reference Norm = 1.0\nEnd\nSolver 2\n  Variable = -dofs 1 Potential\n  Reference Norm = 2.0\n  Nonlinear System Norm Degree = 1\nEnd\nSolver 3\n  Equation = Heat Equation\n  Reference Norm = 3.0\nEnd\nSolver 4\n  Variable = Flow[Velocity:2 Pressure:1]\nEnd\n";
        let s = reference_solvers(sif).unwrap();
        assert_eq!(s.len(), 3);
        assert_eq!(s[0], RefSolver { number: 1, variable: Ok("temperature".into()) });
        assert!(s[1].variable.as_ref().unwrap_err().contains("norm option"));
        assert!(s[2].variable.is_err());
    }

    #[test]
    fn one_liners_merge_and_multiple_runs_are_out_of_class() {
        let s = reference_solvers("Solver 1\n Variable = T\nEnd\nSolver 1 :: Reference Norm = 1\n").unwrap();
        assert_eq!(s, vec![RefSolver { number: 1, variable: Ok("t".into()) }]);
        assert!(reference_solvers("Solver 1\n Variable = T\n Reference Norm = 1\nEnd\nRUN\nRUN\n").is_err());
        assert!(reference_solvers("Solver 3\n Variable = T\nEnd\nSolver 3\n Show Norm Index = 2\nEnd\n").is_err());
    }

    #[test]
    fn golden_norm_rebuilds_values_in_dof_order() {
        // Nodes 1..3 stored with dofs 3,1,2: Values = [b, c, a].
        let file = ResultFile {
            steps: vec![TimeStep {
                saves_done: 1,
                ntime: 1,
                sim_time: 0.0,
                vars: vec![VariableData {
                    name: "Temperature".into(),
                    perm: Some(vec![(1, 3), (2, 1), (3, 2)]),
                    perm_size: 3,
                    values: vec![3.0, 1.0, 2.0],
                }],
            }],
        };
        let n = golden_norm(&file, &[("temperature".into(), 3)], "temperature").unwrap();
        assert!(golden_norm(&file, &[("temperature".into(), 4)], "temperature").is_err());
        assert!(!non_nodal(&file, Some(3), "temperature"));
        assert!(non_nodal(&file, Some(2), "temperature"));
        let expect = compute_norm(&[1.0, 2.0, 3.0], Some(&[3, 1, 2]), 1, 3, &NormOpts::default()).unwrap();
        assert_eq!(n, expect);
    }
}
