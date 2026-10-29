//! Static screen of a parsed case: every feature the native stack cannot
//! run, not only the first. The run itself (`child`) remains the authority
//! on whether a case matches; the screen supplies the full blocker set that
//! the ranking needs, because a refusal reports only the first problem.
//!
//! A blocker is a string `class: name`. Classes, in the order the report
//! names the "first missing" one:
//!
//! | class | meaning |
//! |---|---|
//! | `mpi` | the test runs on more than one rank |
//! | `prestep` | `runtest.cmake` runs a tool other than ElmerGrid first |
//! | `plugin` | a Fortran module the test compiles itself |
//! | `solver` | a `Procedure` the registry does not have |
//! | `userfunc` | a user-function `Procedure` property value |
//! | `feature` | a modelling feature (eigen, p/edge elements, mortar, ...) |
//! | `stand-in` | a limit of the test-local mesh reader or basis |
//! | `runtime` | a refusal the engine reported that no rule above predicts |
//! | `fatal` | the solver stopped on a fatal error (a bug or an unlisted gap) |
//! | `crash` | the process panicked or timed out, with the panic site |

use elmer_core::{normalize_key, DependentKind, ParamMap, Value};
use elmer_solver::solver::procedure_of;
use elmer_solver::Model;

/// Order of the classes; the first present decides "first missing".
pub const CLASS_ORDER: [&str; 10] =
    ["mpi", "prestep", "plugin", "solver", "userfunc", "feature", "stand-in", "runtime", "fatal", "crash"];

pub fn class_of(blocker: &str) -> &str {
    blocker.split(':').next().unwrap_or("")
}

/// Sorts by class order, then name, and drops duplicates.
pub fn sort_blockers(blockers: &mut Vec<String>) {
    blockers.sort_by_key(|b| (CLASS_ORDER.iter().position(|c| *c == class_of(b)).unwrap_or(CLASS_ORDER.len()), b.clone()));
    blockers.dedup();
}

/// `./Poisson`, `Poisson.so`, `poisson` -> `poisson`.
fn module_stem(lib: &str) -> String {
    let s = normalize_key(lib);
    let s = s.rsplit('/').next().unwrap_or(&s);
    s.strip_suffix(".so").unwrap_or(s).to_string()
}

fn text(v: &Value) -> Option<String> {
    match v {
        Value::Str(s) | Value::File(s) => Some(normalize_key(s)),
        _ => None,
    }
}

/// `"Lib" "Entry"` or `Lib Entry`.
fn two_names(t: &str) -> Option<(String, String)> {
    let parts: Vec<&str> = t.split('"').map(str::trim).filter(|s| !s.is_empty()).flat_map(str::split_whitespace).collect();
    match parts[..] {
        [lib, entry] => Some((lib.to_string(), entry.to_string())),
        _ => None,
    }
}

fn is_true(v: &Value) -> bool {
    matches!(v, Value::Logical(true))
}

/// Features named by one section's keywords.
fn section_features(kind: &str, p: &ParamMap, out: &mut Vec<String>) {
    let mut feature = |name: &str| out.push(format!("feature: {name}"));
    for (key, v) in &p.0 {
        let key = normalize_key(key);
        let t = text(v).unwrap_or_default();
        match key.as_str() {
            "simulation type" if t == "scanning" => feature("scanning simulation"),
            "eigen analysis" if is_true(v) => feature("eigenvalue analysis"),
            "harmonic analysis" if is_true(v) => feature("harmonic analysis"),
            "linear system complex" if is_true(v) => feature("complex linear system"),
            "discontinuous galerkin" if is_true(v) => feature("discontinuous Galerkin"),
            "mesh levels" | "extruded mesh levels" => feature("mesh levels / extrusion"),
            "element" if t.contains("p:") => feature("p-hierarchic elements"),
            "element" if t.contains("e:") => feature("edge elements"),
            "linear system solver" if t == "multigrid" => feature("multigrid linear solver"),
            "linear system solver" if t == "block" => feature("block linear system"),
            _ => {}
        }
        if kind == "boundary condition" && (key.contains("mortar") || key.contains("periodic") || key.contains("contact")) {
            feature("mortar / periodic / contact BC");
        }
    }
}

/// All blockers of `model` visible without running it. `local_modules` are
/// the Fortran modules the test compiles itself.
pub fn blockers(model: &Model, registry: &elmer_solver::SolverRegistry, local_modules: &[String]) -> Vec<String> {
    let local: Vec<String> = local_modules.iter().map(|m| module_stem(m)).collect();
    let mut out = Vec::new();

    // Equations the Fortran injects a procedure for (MainUtils.F90:1225-1295);
    // `elmer-run` runs none of them without an explicit `Procedure`.
    const LEGACY: [&str; 5] = ["heat equation", "navier-stokes", "stress analysis", "magnetic induction", "mesh update"];
    for p in &model.solvers {
        // A section that never runs needs no implementation; VTU output
        // sections are served by elmer-run's output solver (refusals come
        // from its own check).
        if p.0.is_empty()
            || p.0.get("exec solver").and_then(text).as_deref() == Some("never")
            || elmer_run::output::is_result_output(p)
        {
            continue;
        }
        match procedure_of(p) {
            Some((lib, entry)) => {
                if local.contains(&module_stem(&lib)) || lib.starts_with("./") {
                    out.push(format!("plugin: {}", module_stem(&lib)));
                } else if registry.get(&lib, &entry).is_none() {
                    out.push(format!("solver: {} {}", module_stem(&lib), normalize_key(&entry)));
                }
            }
            None => {
                let eq = p.0.get("equation").and_then(text).unwrap_or_default();
                let completed = elmer_run::case::legacy_procedure(&eq).is_some_and(|(l, e)| registry.get(l, e).is_some());
                if completed {
                    // elmer-run supplies the procedure (AddEquationBasics).
                } else if LEGACY.contains(&eq.as_str()) {
                    out.push(format!("solver: implicit legacy solver `{eq}`"));
                } else {
                    out.push("solver: Coupled / Assembly solver without Procedure".to_string());
                }
            }
        }
    }

    let sections: [(&str, Vec<&ParamMap>); 8] = [
        ("simulation", vec![&model.simulation]),
        ("constants", vec![&model.constants]),
        ("solver", model.solvers.iter().collect()),
        ("equation", model.equations.iter().collect()),
        ("material", model.materials.iter().collect()),
        ("body force", model.body_forces.iter().collect()),
        ("initial condition", model.initial_conditions.iter().collect()),
        ("boundary condition", model.boundary_conditions.iter().collect()),
    ];
    for (kind, maps) in sections {
        for p in maps {
            section_features(kind, p, &mut out);
            for (key, v) in &p.0 {
                if kind == "solver" && key == "procedure" {
                    continue; // judged above
                }
                let named = match v {
                    Value::Dependent(d) => match &d.kind {
                        DependentKind::Procedure { library, entry } => Some((library.clone(), entry.clone())),
                        _ => None,
                    },
                    // `Bulk Assembly Procedure = File "./Lib" "Entry"`
                    Value::Str(t) | Value::File(t) if key.contains("procedure") => two_names(t),
                    _ => None,
                };
                if let Some((library, entry)) = named {
                    let stem = module_stem(&library);
                    if local.contains(&stem) || library.starts_with("./") {
                        out.push(format!("plugin: {stem}"));
                    } else {
                        out.push(format!("userfunc: {stem} {}", normalize_key(&entry)));
                    }
                }
            }
        }
    }
    sort_blockers(&mut out);
    out
}

/// Which crate on main provides the blocker, when one does and `elmer-run`
/// has not wired it: procedures registered by the other physics crates, and
/// the legacy equations whose library exists.
pub fn provider_of(blocker: &str) -> Option<&'static str> {
    let rest = blocker.strip_prefix("solver: ")?;
    if let Some(eq) = rest.strip_prefix("implicit legacy solver `") {
        return match eq.trim_end_matches('`') {
            "heat equation" => Some("elmer-physics-heat"),
            "stress analysis" => Some("elmer-physics-elastic"),
            _ => None,
        };
    }
    let (lib, entry) = rest.split_once(' ')?;
    type Register = fn(&mut elmer_solver::SolverRegistry);
    let providers: [(&str, Register); 3] = [
        ("elmer-physics-elastic", elmer_physics_elastic::register),
        ("elmer-physics-flux", elmer_physics_flux::register),
        ("elmer-savedata", elmer_savedata::register),
    ];
    providers.into_iter().find_map(|(name, register)| {
        let mut reg = elmer_solver::SolverRegistry::new();
        register(&mut reg);
        reg.get(lib, entry).is_some().then_some(name)
    })
}

/// `Solver 3: ...` -> `Solver N: ...`, so one gap in several sections is one
/// blocker.
pub fn generic(s: &str) -> String {
    let mut out = String::new();
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if c.is_ascii_digit() && out.ends_with(' ') && (out.contains("Solver") || out.contains("Condition")) {
            while chars.peek().is_some_and(|d| d.is_ascii_digit()) {
                chars.next();
            }
            out.push('N');
        } else {
            out.push(c);
        }
    }
    out
}

/// Solver library stems of the `Solver` sections (legacy `Equation` names
/// mapped to the library MainUtils injects), sorted.
pub fn libs(model: &Model) -> Vec<String> {
    let mut libs: Vec<String> = model
        .solvers
        .iter()
        .filter(|p| !p.0.is_empty() && !elmer_run::output::is_result_output(p))
        .filter_map(elmer_cli::capability::solver_lib_of)
        .map(|l| module_stem(&l))
        .collect();
    libs.sort();
    libs.dedup();
    libs
}

/// A fatal solver message that says a piece is not installed or not ported,
/// as the blocker it names. Such a run is a refusal in all but the exit path.
pub fn capability_gap(msg: &str) -> Option<String> {
    let m = msg.trim();
    if m.contains("no MATC engine installed") {
        return Some("runtime: MATC engine not installed in elmer-run (MATC value evaluated)".into());
    }
    if m.contains("not wired") || m.contains("NotWired") {
        return Some(format!("runtime: {}", generic(m)));
    }
    if m.contains("not ported") || m.contains("not implemented") || m.contains("unsupported") || m.contains("Unsupported") {
        return Some(format!("runtime: {}", generic(m)));
    }
    None
}

/// Maps an `elmer-run` refusal message onto the blocker the screen would
/// have named, so the two do not count one gap twice.
pub fn from_refusal(msg: &str) -> String {
    if let Some(rest) = msg.split("Procedure \"").nth(1) {
        let parts: Vec<&str> = rest.split('"').collect();
        if let [lib, _, entry, ..] = parts[..] {
            return format!("solver: {} {}", module_stem(lib), normalize_key(entry));
        }
    }
    if msg.contains("no Procedure") {
        return "solver: implicit legacy solver `(no Procedure)`".into();
    }
    if msg.starts_with("stand-in:") {
        return msg.to_string();
    }
    // elmer-sif returns `Mesh DB "d" "n"` as (d, "d/n") and elmer-run joins
    // both onto the case directory, so any `d` other than "." doubles.
    if msg.contains("mesh.nodes: No such file") {
        let parts: Vec<&str> = msg.split('/').collect();
        if parts.windows(2).any(|w| w[0] == w[1] && w[0] != ".") {
            return "runtime: Mesh DB directory other than `.` (elmer-run joins it twice)".into();
        }
    }
    // A refusal worded `unsupported: ...` or `not wired: ...` by SolverError.
    let msg = msg.strip_prefix("unsupported: ").or_else(|| msg.strip_prefix("not wired: ")).unwrap_or(msg);
    if msg.ends_with(" sections") || msg == "Control section" {
        return format!("feature: {msg}");
    }
    if let Some(gap) = capability_gap(msg) {
        return gap;
    }
    format!("runtime: {}", generic(msg))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pm(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    fn st(s: &str) -> Value {
        Value::Str(s.into())
    }

    #[test]
    fn solver_screen() {
        let model = Model {
            solvers: vec![
                pm(&[("procedure", st("\"HeatSolve\" \"HeatSolver\""))]),
                pm(&[("procedure", st("\"SaveData\" \"SaveScalars\""))]),
                pm(&[("procedure", st("\"Poisson\" \"PoissonSolver\""))]),
                pm(&[("equation", st("heat equation"))]),
                pm(&[("equation", st("magnetic induction"))]),
            ],
            ..Model::default()
        };
        let reg = elmer_run::solver_registry();
        let b = blockers(&model, &reg, &["Poisson".to_string()]);
        assert_eq!(
            b,
            // SaveScalars is registered and the legacy heat equation completed
            // by elmer-run, so neither is a blocker.
            vec!["plugin: poisson", "solver: implicit legacy solver `magnetic induction`"]
        );
    }

    #[test]
    fn features() {
        let model = Model {
            simulation: pm(&[("simulation type", st("scanning"))]),
            solvers: vec![pm(&[("procedure", st("\"HeatSolve\" \"HeatSolver\"")), ("eigen analysis", Value::Logical(true))])],
            boundary_conditions: vec![pm(&[("mortar bc", Value::Integer(vec![1]))])],
            ..Model::default()
        };
        let b = blockers(&model, &elmer_run::solver_registry(), &[]);
        assert_eq!(
            b,
            vec!["feature: eigenvalue analysis", "feature: mortar / periodic / contact BC", "feature: scanning simulation"]
        );
    }

    #[test]
    fn refusal_mapping() {
        assert_eq!(from_refusal("Solver 2: Procedure \"SaveData\" \"SaveScalars\""), "solver: savedata savescalars");
        assert_eq!(from_refusal("Solver 3: no Procedure"), "solver: implicit legacy solver `(no Procedure)`");
        assert_eq!(from_refusal("Boundary Condition 2: default target"), "runtime: Boundary Condition N: default target");
    }
}
