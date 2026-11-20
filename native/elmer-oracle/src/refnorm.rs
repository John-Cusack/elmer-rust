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
//! solvers with default norm options (no `Nonlinear System Norm *`, no
//! `Norm Permutation`), set inside a `Solver N` section of a single-run
//! sif. The variable is the solver's explicit `Variable` (scalar,
//! `-dofs N` or a bracketed vector, rebuilt from its saved one-dof
//! components), or, for a solver without one, the single variable the
//! result preamble assigns to the solver's `Equation`.

use crate::norm::{compute_norm, NormOpts};
use crate::result::{PreambleVar, ResultFile};
use crate::runner::{keyword_of, opens_table, strip_sif_comment};

/// One `Reference Norm` solver of a sif.
#[derive(Debug, Clone, PartialEq)]
pub struct RefSolver {
    pub number: u32,
    /// Lowercased variable name without its bracketed component list
    /// (`flow` for `Flow[Velocity:2 Pressure:1]`), or why the solver is out
    /// of class. [`NO_VARIABLE`] marks a solver that names none.
    pub variable: Result<String, String>,
    /// Lowercased `Equation` value, for a default variable.
    pub equation: Option<String>,
    /// Lowercased basename of the `Procedure` file, or the legacy
    /// solver's when the equation name selects one.
    pub procedure: Option<String>,
    /// Lowercased entry subroutine of the `Procedure` (its second string).
    pub entry: Option<String>,
    /// `Coupled Solver = True`: the norm runs over the coupled unknowns.
    pub coupled: bool,
    /// Lowercased solver-level `Mesh` (the solver runs on another mesh).
    pub mesh: Option<String>,
}

/// Solver entries that assign `Solver % Variable % Norm` themselves
/// instead of leaving it to `ComputeNorm` (pristine `fem/src/modules`,
/// grep `Variable % Norm =`, keyed by the enclosing solver subroutine).
/// A mismatch in one of these is a grader gap, not a capture failure.
const SELF_NORM: [&str; 34] = [
    "acousticssolver",
    "compressibilityscale",
    "solidphasecons",
    "craigbamptonsolver",
    "distancesolver1",
    "emportsolver",
    "fixtangentvelo",
    "fluxsolver",
    "freesurfacereduced",
    "helmholtzbemsolver",
    "isosurfacesolver",
    "levelsetdistance",
    "levelsetintegrate",
    "bsolver",
    "magnetodynamicscalcfields",
    "whitneyavsolver",
    "advdiffsolver",
    "particleadvector",
    "particledynamics",
    "phasechangesolve",
    "poissonbemsolver",
    "poissonboltzmannsolve",
    "parallelprojecttoplane",
    "projecttoplane",
    "elmerpostoutputsolver",
    "richardspostprocess",
    "savedependence",
    "saveline",
    "savescalars",
    "statcurrentsolverlegacy",
    "steadyphasechange",
    "transientcost",
    "vectorhelmholtznodal",
    "vorticitysolver",
];

/// Solver entries that write into another solver's variable after that
/// solver computed its norm (`FixTangentVelo.F90:208` adds to the flow
/// solution's values), so the saved field is not the one Elmer normed.
pub const MODIFIES_OTHER: [&str; 1] = ["fixtangentvelo"];

impl RefSolver {
    /// Whether the solver's entry sets its own norm (see [`SELF_NORM`]).
    pub fn self_norm(&self) -> bool {
        self.entry.as_deref().is_some_and(|e| SELF_NORM.contains(&e))
    }

    /// Whether the solver's entry rewrites another solver's variable
    /// (see [`MODIFIES_OTHER`]).
    pub fn modifies_other(&self) -> bool {
        self.entry.as_deref().is_some_and(|e| MODIFIES_OTHER.contains(&e))
    }
}

/// Elmer's legacy equation names (`MainUtils.F90`, `AddEquationBasics`):
/// (equation, default variable, procedure file).
const LEGACY: [(&str, &str, &str); 5] = [
    ("navier-stokes", "flow solution", "flowsolve"),
    ("magnetic induction", "magnetic field", "magneticsolve"),
    ("stress analysis", "displacement", "stresssolve"),
    ("mesh update", "mesh update", "meshsolve"),
    ("heat equation", "temperature", "heatsolve"),
];

/// Procedures that add `Nonlinear System Norm DOFs = DOFs - 1` (velocity
/// only, pressure last) unless the sif sets it: `FlowSolve.F90`,
/// `IncompressibleNS.F90`, `FilmFlowSolver.F90`.
const VELOCITY_NORM: [&str; 3] = ["flowsolve", "incompressiblens", "filmflowsolver"];

impl RefSolver {
    /// The norm options the solver runs with: defaults, except the
    /// velocity-only norm the flow solvers set internally.
    pub fn norm_opts(&self, dofs: usize) -> NormOpts {
        let velocity = self.procedure.as_deref().is_some_and(|p| VELOCITY_NORM.contains(&p));
        NormOpts { norm_dofs: (velocity && dofs > 1).then(|| dofs - 1), ..NormOpts::default() }
    }
}

/// Reason given for a solver without an explicit `Variable`.
pub const NO_VARIABLE: &str = "no explicit Variable (legacy or default name)";

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

/// Whether a solver asks for the parallel consistent norm (owned rows
/// only). Serially Elmer ignores it.
pub fn consistent(text: &str, number: u32) -> bool {
    let mut section: Option<u32> = None;
    for raw in text.lines() {
        let code = strip_sif_comment(raw).trim();
        let lower = code.to_lowercase();
        let (owner, key) = match code.split_once("::") {
            Some((sec, rest)) => (sec.trim().to_lowercase().strip_prefix("solver").and_then(|n| n.trim().parse().ok()), rest),
            None => (section, code),
        };
        if code.split_once("::").is_none() {
            if lower == "end" {
                section = None;
            } else if let Some(n) = lower.strip_prefix("solver").and_then(|r| r.trim().parse().ok()) {
                section = Some(n);
            }
        }
        if owner == Some(number)
            && keyword_of(key).as_deref() == Some("nonlinear system consistent norm")
            && key.to_lowercase().contains("true")
        {
            return true;
        }
    }
    false
}

/// Owned flag per local node of rank `rank`, from ElmerGrid's
/// `partitioning.<np>/part.<rank+1>.{nodes,shared}`: a shared node belongs
/// to the first partition its `.shared` line lists (Elmer's `Neighbours(1)`,
/// `ReadSharedFile`), any other node to the rank itself.
pub fn owned_nodes(partitioning: &std::path::Path, rank: usize) -> Result<Vec<bool>, String> {
    let read = |name: &str| {
        std::fs::read_to_string(partitioning.join(name)).map_err(|e| format!("{name}: {e}"))
    };
    let nodes = read(&format!("part.{}.nodes", rank + 1))?;
    let shared = read(&format!("part.{}.shared", rank + 1))?;
    let mut owner = std::collections::HashMap::new();
    for l in shared.lines() {
        let v: Vec<i64> = l.split_whitespace().filter_map(|x| x.parse().ok()).collect();
        if v.len() >= 3 {
            owner.insert(v[0], v[2] - 1);
        }
    }
    nodes
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let g: i64 = l.split_whitespace().next().and_then(|x| x.parse().ok()).ok_or_else(|| format!("bad node line {l:?}"))?;
            Ok(owner.get(&g).is_none_or(|&o| o == rank as i64))
        })
        .collect()
}

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
        .map(|(number, k)| {
            let equation = equation(&k);
            RefSolver {
                number,
                variable: solver_variable(&k),
                procedure: procedure(&k).or_else(|| {
                    let eq = equation.as_deref()?;
                    LEGACY.iter().find(|l| l.0 == eq).map(|l| l.2.to_string())
                }),
                entry: keys_entry(&k),
                mesh: k.iter().rev().find(|(k, _)| k == "mesh").and_then(|(_, v)| {
                    v.split_whitespace().next().map(|m| m.trim_matches('"').to_lowercase())
                }),
                coupled: k.iter().any(|(k, v)| k == "coupled solver" && v.to_lowercase().contains("true")),
                equation,
            }
        })
        .collect())
}

/// Expand top-level `include <file>` lines, as Elmer's sif reader does
/// textually, when the included file only defines MATC variables (`$`
/// lines), comments and blank lines: `params.sif` constants, or the
/// `mesh.names` that ElmerGrid writes. Each name is tried against `bases`
/// in order (the golden run dir holds generated meshes, the source dir the
/// rest). An include that defines sections, or that cannot be read, is an
/// error, because the override goes into the main sif.
pub fn expand_includes(text: &str, bases: &[&std::path::Path]) -> Result<String, String> {
    let mut out = String::with_capacity(text.len());
    let mut in_section = false;
    for raw in text.lines() {
        let code = strip_sif_comment(raw).trim();
        let lower = code.to_lowercase();
        if !in_section && lower.starts_with("include") && !lower.starts_with("include path") {
            let name = code[7..].trim().trim_matches('"').trim();
            let body = bases
                .iter()
                .find_map(|b| std::fs::read(b.join(name)).ok())
                .ok_or_else(|| format!("sif uses Include {name} (not found)"))?;
            let body: String = body.iter().map(|&b| char::from(b)).collect();
            for line in body.lines() {
                let c = strip_sif_comment(line).trim();
                if !(c.is_empty() || c.starts_with('$')) {
                    return Err(format!("sif uses Include {name} (defines more than $ variables)"));
                }
            }
            out.push_str(&body);
            if !body.ends_with('\n') {
                out.push('\n');
            }
            continue;
        }
        if in_section {
            in_section = lower != "end";
        } else if !(lower.is_empty()
            || lower.starts_with('$')
            || lower.contains("::")
            || lower == "run"
            || lower == "stop"
            || lower.starts_with("echo ")
            || lower.starts_with("check "))
        {
            // One-line directives (`Check Keywords "Warn"`) open no section.
            in_section = true;
        }
        out.push_str(raw);
        out.push('\n');
    }
    Ok(out)
}

/// For each top-level `RUN` of a sif, the text Elmer has read when it
/// compares after that RUN: everything before it, with the earlier `RUN`
/// lines dropped (keys merge in file order, later ones win).
pub fn run_prefixes(text: &str) -> Vec<String> {
    let lines: Vec<&str> = text.lines().collect();
    let runs = crate::runner::run_lines(text);
    runs.iter()
        .map(|&r| {
            lines[..r]
                .iter()
                .enumerate()
                .filter(|(i, _)| !runs.contains(i))
                .map(|(_, l)| format!("{l}\n"))
                .collect()
        })
        .collect()
}

/// The `Mesh DB` mesh name of a sif (its Header, or a `Header ::`
/// one-liner), lowercased.
pub fn primary_mesh(text: &str) -> Option<String> {
    text.lines().find_map(|l| {
        let c = strip_sif_comment(l);
        // `Mesh DB "dir" "name"` carries no `=`; the name is its second
        // (quoted or bare) argument.
        let tail = c.split("::").last()?.trim().to_lowercase().strip_prefix("mesh db")?.to_string();
        let quoted: Vec<&str> = tail.split('"').skip(1).step_by(2).collect();
        let name = if quoted.len() >= 2 { quoted[1].to_string() } else { tail.split_whitespace().nth(1)?.to_string() };
        Some(name.trim().to_string())
    })
}

/// Which of several result files holding a solver's variable is the
/// solver's own: the one in its `Mesh`'s directory, after adaptive
/// refinement the last refined mesh (`RefinedMesh<N>`, highest N), else the
/// one in the `Mesh DB` mesh's directory. `None` when that leaves no single
/// file.
pub fn own_file<'a>(paths: &[&'a str], solver_mesh: Option<&str>, primary: Option<&str>, adaptive: bool) -> Option<&'a str> {
    let dir = |p: &str| p.split('/').next().unwrap_or("").to_lowercase();
    let pick = |m: &str| {
        let hits: Vec<&'a str> = paths.iter().copied().filter(|p| dir(p) == m).collect();
        (hits.len() == 1).then(|| hits[0])
    };
    if let Some(m) = solver_mesh {
        return pick(m);
    }
    if adaptive {
        let refined = paths
            .iter()
            .copied()
            .filter_map(|p| Some((dir(p).strip_prefix("refinedmesh")?.parse::<u32>().ok()?, p)))
            .max_by_key(|(n, _)| *n);
        if let Some((_, p)) = refined {
            return Some(p);
        }
    }
    pick(primary?)
}

/// Every solver's lowercased `Equation`, by solver number.
pub fn solver_equations(text: &str) -> Vec<(u32, String)> {
    let mut out = Vec::new();
    let mut section: Option<u32> = None;
    for raw in text.lines() {
        let code = strip_sif_comment(raw).trim();
        let lower = code.to_lowercase();
        if let Some(n) = section {
            if lower == "end" {
                section = None;
            } else if keyword_of(code).as_deref() == Some("equation") {
                if let Some(eq) = equation(&[(String::from("equation"), code.split_once('=').map_or("", |(_, v)| v).to_string())]) {
                    out.push((n, eq));
                }
            }
            continue;
        }
        if let Some(n) = lower.strip_prefix("solver").and_then(|r| r.trim().parse().ok()) {
            section = Some(n);
        }
    }
    out
}

/// Entry subroutine of the `Procedure` (its second quoted string).
fn keys_entry(keys: &[(String, String)]) -> Option<String> {
    let (_, v) = keys.iter().rev().find(|(k, _)| k == "procedure")?;
    Some(v.split('"').nth(3)?.trim().to_lowercase())
}

/// Basename of the `Procedure` file, without directory or extension.
fn procedure(keys: &[(String, String)]) -> Option<String> {
    let (_, v) = keys.iter().rev().find(|(k, _)| k == "procedure")?;
    let first = v.split('"').nth(1)?.trim().to_lowercase();
    let base = first.rsplit('/').next()?.to_string();
    Some(base.split('.').next().unwrap_or(&base).to_string())
}

fn equation(keys: &[(String, String)]) -> Option<String> {
    let (_, v) = keys.iter().rev().find(|(k, _)| k == "equation")?;
    let v = v.trim();
    let v = v.strip_prefix("String").or_else(|| v.strip_prefix("string")).unwrap_or(v);
    Some(v.trim().trim_matches('"').trim().to_lowercase())
}

/// The default variable of a solver without `Variable`: the legacy
/// solver's when the equation name selects one, else the one preamble
/// root variable owned by its `Equation`, when exactly one exists and no
/// other solver of the sif shares that equation name.
pub fn default_variable(
    vars: &[PreambleVar],
    solver: &RefSolver,
    equations: &[(u32, String)],
) -> Result<String, String> {
    if solver.coupled {
        return Err(String::from("no explicit Variable; coupled solver"));
    }
    let eq = solver.equation.as_deref().ok_or("no explicit Variable and no Equation")?;
    if let Some(l) = LEGACY.iter().find(|l| l.0 == eq) {
        return Ok(l.1.to_string());
    }
    if equations.iter().any(|(n, e)| *n != solver.number && e == eq) {
        return Err(format!("no explicit Variable; equation {eq:?} shared by several solvers"));
    }
    let component_names: std::collections::BTreeSet<String> = vars
        .iter()
        .filter(|v| v.dofs > 1)
        .flat_map(|v| v.components().unwrap_or_default())
        .collect();
    let roots: Vec<&PreambleVar> = vars
        .iter()
        .filter(|v| v.equation == eq && !component_names.contains(&v.name.to_lowercase()))
        .collect();
    match roots.as_slice() {
        [v] => Ok(v.base_name()),
        [] => Err(format!("no explicit Variable; no saved variable of equation {eq:?}")),
        _ => Err(format!("no explicit Variable; {} variables of equation {eq:?}", roots.len())),
    }
}

fn solver_variable(keys: &[(String, String)]) -> Result<String, String> {
    if let Some((k, _)) = keys.iter().find(|(k, _)| NORM_KEYS.contains(&k.as_str())) {
        return Err(format!("norm option {k}"));
    }
    let Some((_, value)) = keys.iter().rev().find(|(k, _)| k == "variable") else {
        return Err(String::from(NO_VARIABLE));
    };
    let v = value.trim().trim_start_matches("-global").trim();
    let v = v
        .strip_prefix("String")
        .or_else(|| v.strip_prefix("string"))
        .unwrap_or(v)
        .trim()
        .trim_matches('"')
        .trim();
    // A bracketed vector is named by its part before the bracket; quotes
    // may also follow `-dofs N` (`-dofs 3 "Displacement"`).
    let lower = v.replace('"', " ").split('[').next().unwrap_or("").to_lowercase();
    let mut parts = lower.split_whitespace().peekable();
    let mut name = Vec::new();
    while let Some(p) = parts.next() {
        if p == "-dofs" {
            match parts.next().and_then(|n| n.parse::<u32>().ok()) {
                Some(_) => continue,
                None => return Err(String::from("bad -dofs")),
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

/// The preamble entry a solver variable name refers to.
fn preamble_entry<'a>(vars: &'a [PreambleVar], variable: &str) -> Result<&'a PreambleVar, String> {
    let mut m = vars.iter().filter(|v| v.base_name() == variable);
    let v = m.next().ok_or_else(|| format!("{variable} not in the preamble"))?;
    if m.next().is_some() {
        return Err(format!("{variable} twice in the preamble"));
    }
    Ok(v)
}

/// Whether `variable`'s Perm reaches past the mesh nodes (edge, face or
/// DG unknowns), where the grader's nodal reconstruction may not follow
/// the solver's own norm path.
pub fn non_nodal(vars: &[PreambleVar], nodes: Option<usize>, variable: &str) -> bool {
    preamble_entry(vars, variable).is_ok_and(|v| nodes.is_some_and(|n| v.perm_size > n))
}

/// Norm of `variable` at the last saved step of `file`, computed as the
/// serial `ComputeNorm` does from the rebuilt `Values` and `Perm` arrays.
/// Only one-dof variables are saved, so a vector's `Values` are rebuilt
/// from its components: component `j` at Perm index `k` is
/// `Values(dofs*(k-1)+j)`.
pub fn golden_norm(
    file: &ResultFile,
    vars: &[PreambleVar],
    variable: &str,
    opts: impl Fn(usize) -> NormOpts,
) -> Result<f64, String> {
    let (values, perm, dofs) = golden_values(file, vars, variable)?;
    compute_norm(&values, Some(&perm), dofs, perm.len(), &opts(dofs)).map_err(|e| e.to_string())
}

/// The solver's `Values` (dof order), `Perm` and DOFs of `variable`,
/// rebuilt from the last saved step of `file`.
pub fn golden_values(file: &ResultFile, vars: &[PreambleVar], variable: &str) -> Result<(Vec<f64>, Vec<i32>, usize), String> {
    let entry = preamble_entry(vars, variable)?;
    let comps = entry
        .components()
        .ok_or_else(|| format!("{variable}: components of {:?} not listable", entry.name))?;
    let dofs = entry.dofs;
    let step = file.steps.last().ok_or("no saved step")?;
    let mut values: Vec<f64> = Vec::new();
    let mut perm: Option<Vec<i32>> = None;
    let mut first_pairs: Option<&Vec<(i32, i32)>> = None;
    for (j, comp) in comps.iter().enumerate() {
        let mut matches = step.vars.iter().filter(|v| v.name.to_lowercase() == *comp);
        let var = matches.next().ok_or_else(|| format!("{comp} not in the last saved step"))?;
        if matches.next().is_some() {
            return Err(format!("{comp} appears twice"));
        }
        let pairs = var.perm.as_ref().ok_or("NULL Perm")?;
        match first_pairs {
            None => first_pairs = Some(pairs),
            Some(p) if p == pairs => {}
            Some(_) => return Err(format!("{comp}: Perm differs from the first component")),
        }
        if perm.is_none() {
            let n_values = pairs.iter().map(|&(_, k)| k).max().unwrap_or(0).max(0) as usize;
            values = vec![f64::NAN; n_values * dofs];
            let mut p = vec![0i32; var.perm_size.max(0) as usize];
            for &(node, k) in pairs {
                let slot = p.get_mut((node as usize).wrapping_sub(1)).ok_or("Perm entry out of range")?;
                *slot = k;
            }
            perm = Some(p);
        }
        if var.values.len() != pairs.len() {
            return Err(format!("{comp}: {} values for {} Perm entries", var.values.len(), pairs.len()));
        }
        for (&(_, k), &v) in pairs.iter().zip(&var.values) {
            let at = dofs * (k as usize - 1) + j;
            *values.get_mut(at).ok_or("Perm entry out of range")? = v;
        }
    }
    if values.iter().any(|v| v.is_nan()) {
        return Err(String::from("values not covered by the Perm (not all saved)"));
    }
    // Edge, face, bubble or DG unknowns beyond the nodal Perm are not in
    // the file, yet ComputeNorm sums them: out of class.
    if entry.size != values.len() {
        return Err(format!("{} values, {} saved (non-nodal unknowns)", entry.size, values.len()));
    }
    Ok((values, perm.unwrap_or_default(), dofs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result::{TimeStep, VariableData};

    fn pv(name: &str, size: usize, perm_size: usize, dofs: usize, eq: &str) -> PreambleVar {
        PreambleVar { name: name.into(), size, perm_size, dofs, equation: eq.into() }
    }

    fn rs(number: u32, variable: &str, equation: Option<&str>) -> RefSolver {
        RefSolver {
            number,
            variable: Ok(variable.into()),
            equation: equation.map(Into::into),
            procedure: None,
            entry: None,
            coupled: false,
            mesh: None,
        }
    }

    #[test]
    fn finds_reference_solvers_and_their_variables() {
        let sif = "Header\n Mesh DB \".\" \"m\"\nEnd\nSolver 1\n  Equation = Heat\n  Variable = Temperature\n  Reference Norm = 1.0\nEnd\nSolver 2\n  Variable = -dofs 1 Potential\n  Reference Norm = 2.0\n  Nonlinear System Norm Degree = 1\nEnd\nSolver 3\n  Equation = \"Heat Equation\"\n  Reference Norm = 3.0\nEnd\nSolver 4\n  Variable = Flow[Velocity:2 Pressure:1]\n  Reference Norm = 4.0\nEnd\nSolver 5\n  Variable = -dofs 3 Displacement\n  Reference Norm = 5.0\nEnd\n";
        let s = reference_solvers(sif).unwrap();
        assert_eq!(s.len(), 5);
        assert_eq!(s[0], rs(1, "temperature", Some("heat")));
        assert_eq!(s[2].procedure.as_deref(), Some("heatsolve"));
        let flux = reference_solvers("Solver 2\n Procedure = \"FluxSolver\" \"FluxSolver\"\n Reference Norm = 1\nEnd\n").unwrap();
        assert_eq!(flux[0].entry.as_deref(), Some("fluxsolver"));
        assert!(flux[0].self_norm() && !s[0].self_norm());
        assert!(s[1].variable.as_ref().unwrap_err().contains("norm option"));
        assert_eq!(s[2].variable.as_ref().unwrap_err(), NO_VARIABLE);
        assert_eq!(s[2].equation.as_deref(), Some("heat equation"));
        assert_eq!(s[3].variable.as_deref(), Ok("flow"));
        assert_eq!(s[4].variable.as_deref(), Ok("displacement"));
        let quoted = reference_solvers("Solver 1\n Variable = -dofs 3 \"Displacement\"\n Reference Norm = 1\nEnd\n").unwrap();
        assert_eq!(quoted[0].variable.as_deref(), Ok("displacement"));
        assert_eq!(solver_equations(sif), vec![(1, "heat".into()), (3, "heat equation".into())]);
    }

    #[test]
    fn includes_of_matc_variables_expand() {
        let dir = std::env::temp_dir().join(format!("refnorm-include-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("m")).unwrap();
        std::fs::write(dir.join("params.sif"), "! constants\n$ h = 0.1\n\n").unwrap();
        std::fs::write(dir.join("m/mesh.names"), "$ wall = 3\n").unwrap();
        std::fs::write(dir.join("bad.sif"), "Solver 1\n Variable = T\nEnd\n").unwrap();
        let sif = "Check Keywords \"Warn\"\ninclude params.sif\nINCLUDE m/mesh.names\nHeader\n  Include Path \"\"\nEnd\nSolver 1\n Variable = T\n Reference Norm = 1\nEnd\n";
        let e = expand_includes(sif, &[&dir]).unwrap();
        assert!(e.contains("$ h = 0.1") && e.contains("$ wall = 3") && e.contains("Include Path"));
        assert_eq!(reference_solvers(&e).unwrap().len(), 1);
        assert!(expand_includes("include bad.sif\n", &[&dir]).unwrap_err().contains("defines more"));
        assert!(expand_includes("include none.sif\n", &[&dir]).unwrap_err().contains("not found"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn run_prefixes_accumulate_keys_up_to_each_run() {
        let sif = "Solver 1\n Variable = T\nEnd\nSolver 1 :: Reference Norm = 1\nRUN\nSolver 2 :: Reference Norm = 2\nRUN\n";
        let p = run_prefixes(sif);
        assert_eq!(p.len(), 2);
        assert_eq!(reference_solvers(&p[0]).unwrap().len(), 1);
        let second = reference_solvers(&p[1]).unwrap();
        assert_eq!(second.iter().map(|s| s.number).collect::<Vec<_>>(), [1, 2]);
    }

    #[test]
    fn own_file_follows_the_solver_mesh() {
        let p = ["square-tri/case.result", "square-tri-isomesh/case.result", "RefinedMesh2/case.result", "RefinedMesh10/case.result"];
        assert_eq!(own_file(&p, None, Some("square-tri"), false), Some("square-tri/case.result"));
        assert_eq!(own_file(&p, Some("square-tri-isomesh"), Some("square-tri"), false), Some("square-tri-isomesh/case.result"));
        assert_eq!(own_file(&p, None, Some("square-tri"), true), Some("RefinedMesh10/case.result"));
        assert_eq!(own_file(&p, Some("global"), None, false), None);
        assert_eq!(primary_mesh("Header\n  Mesh DB \".\" \"Square-Tri\"\n!  Mesh DB \".\" \"q\"\nEnd\n").as_deref(), Some("square-tri"));
        assert_eq!(primary_mesh("Header :: Mesh DB \".\" \"strip\"\n").as_deref(), Some("strip"));
    }

    #[test]
    fn one_liners_merge_and_multiple_runs_are_out_of_class() {
        let s = reference_solvers("Solver 1\n Variable = T\nEnd\nSolver 1 :: Reference Norm = 1\n").unwrap();
        assert_eq!(s, vec![rs(1, "t", None)]);
        assert!(reference_solvers("Solver 1\n Variable = T\n Reference Norm = 1\nEnd\nRUN\nRUN\n").is_err());
        assert!(reference_solvers("Solver 3\n Variable = T\nEnd\nSolver 3\n Show Norm Index = 2\nEnd\n").is_err());
    }

    #[test]
    fn default_variable_is_the_single_root_of_the_equation() {
        let vars = [
            pv("flow solution[velocity:2 pressure:1]", 6, 2, 3, "navier-stokes"),
            pv("velocity 1", 2, 2, 1, "navier-stokes"),
            pv("velocity 2", 2, 2, 1, "navier-stokes"),
            pv("pressure", 2, 2, 1, "navier-stokes"),
            pv("temperature", 2, 2, 1, "heat equation"),
            pv("heat flux", 2, 2, 1, "heat equation"),
        ];
        let mut s = rs(1, "", Some("navier-stokes"));
        s.variable = Err(NO_VARIABLE.into());
        let eqs = [(1, "navier-stokes".to_string()), (2, "heat solver".to_string())];
        assert_eq!(default_variable(&vars, &s, &eqs).as_deref(), Ok("flow solution"));
        s.coupled = true;
        assert!(default_variable(&vars, &s, &eqs).unwrap_err().contains("coupled"));
        s.coupled = false;
        // Non-legacy: two roots of one equation are ambiguous; so is a
        // shared equation.
        let vars = vars.map(|mut v| {
            v.equation = v.equation.replace("heat equation", "heat solver");
            v
        });
        s.number = 2;
        s.equation = Some("heat solver".into());
        assert!(default_variable(&vars, &s, &eqs).unwrap_err().contains("2 variables"));
        let shared = [(2, "heat solver".to_string()), (3, "heat solver".to_string())];
        assert!(default_variable(&vars, &s, &shared).unwrap_err().contains("shared"));
    }

    fn step(vars: Vec<VariableData>) -> ResultFile {
        ResultFile {
            groups: Vec::new(),
            steps: vec![TimeStep { saves_done: 1, ntime: 1, sim_time: 0.0, vars }],
        }
    }

    fn var(name: &str, perm: &[(i32, i32)], values: &[f64]) -> VariableData {
        VariableData { name: name.into(), perm: Some(perm.to_vec()), perm_size: perm.len() as i32, values: values.to_vec() }
    }

    #[test]
    fn golden_norm_rebuilds_values_in_dof_order() {
        // Nodes 1..3 stored with dofs 3,1,2: Values = [b, c, a].
        let file = step(vec![var("Temperature", &[(1, 3), (2, 1), (3, 2)], &[3.0, 1.0, 2.0])]);
        let pre = |size| [pv("temperature", size, 3, 1, "heat equation")];
        let n = golden_norm(&file, &pre(3), "temperature", |_| NormOpts::default()).unwrap();
        assert!(golden_norm(&file, &pre(4), "temperature", |_| NormOpts::default()).is_err());
        assert!(!non_nodal(&pre(3), Some(3), "temperature"));
        assert!(non_nodal(&pre(3), Some(2), "temperature"));
        let expect = compute_norm(&[1.0, 2.0, 3.0], Some(&[3, 1, 2]), 1, 3, &NormOpts::default()).unwrap();
        assert_eq!(n, expect);
    }

    #[test]
    fn golden_norm_interleaves_vector_components() {
        // Two nodes with Perm (1->2, 2->1), dofs 3 = velocity 1, velocity 2,
        // pressure: Values = [u1(n2) u2(n2) p(n2) u1(n1) u2(n1) p(n1)].
        let perm = [(1, 2), (2, 1)];
        let file = step(vec![
            var("velocity 1", &perm, &[1.0, 4.0]),
            var("velocity 2", &perm, &[2.0, 5.0]),
            var("pressure", &perm, &[3.0, 6.0]),
        ]);
        let pre = [pv("flow solution[velocity:2 pressure:1]", 6, 2, 3, "navier-stokes")];
        let n = golden_norm(&file, &pre, "flow solution", |_| NormOpts::default()).unwrap();
        let expect =
            compute_norm(&[4.0, 5.0, 6.0, 1.0, 2.0, 3.0], Some(&[2, 1]), 3, 2, &NormOpts::default()).unwrap();
        assert_eq!(n, expect);
        // FlowSolve's velocity-only norm: Norm DOFs = 2 of 3.
        let mut flow = rs(1, "flow solution", Some("navier-stokes"));
        flow.procedure = Some("flowsolve".into());
        let v = golden_norm(&file, &pre, "flow solution", |d| flow.norm_opts(d)).unwrap();
        let opts = NormOpts { norm_dofs: Some(2), ..NormOpts::default() };
        assert_eq!(v, compute_norm(&[4.0, 5.0, 6.0, 1.0, 2.0, 3.0], None, 3, 2, &opts).unwrap());
        assert_ne!(v, n);
        let complex = pv("potential[potential re:1 potential im:1]", 4, 2, 2, "mgharm");
        assert_eq!(complex.components().unwrap(), ["potential re", "potential im"]);
        assert_eq!(complex.base_name(), "potential");
        let plain = [pv("displacement", 6, 2, 3, "stress analysis")];
        assert_eq!(plain[0].components().unwrap(), ["displacement 1", "displacement 2", "displacement 3"]);
        // A missing component is out of class, not a partial norm.
        let short = step(vec![var("velocity 1", &perm, &[1.0, 4.0])]);
        assert!(golden_norm(&short, &pre, "flow solution", |_| NormOpts::default()).is_err());
    }
}
