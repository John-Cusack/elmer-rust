//! Refusal of model features no provider implements. `Simulation` and the
//! physics kernels never read these keywords, so a case naming one runs to
//! a wrong answer instead of failing; each is refused before setup with the
//! feature it needs. Every rule here was found by a quick test that ran
//! and missed its reference norm, and confirmed by running the pristine
//! ElmerSolver without the keyword, which then reproduced the native norm.

use elmer_core::{normalize_key, Params, ParamMap, Value};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Model, SolverError};

enum Test {
    Present,
    True,
}

struct Rule {
    key: &'static str,
    test: Test,
    reason: &'static str,
}

const fn present(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::Present, reason }
}

const fn flag(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::True, reason }
}

const SIMULATION: &[Rule] = &[
    present("Mesh Levels", "uniform mesh refinement (SplitMeshEqual) is not implemented"),
    present("Extruded Mesh Levels", "mesh extrusion (MeshExtrude) is not implemented"),
    flag("Increase Element Order", "quadratic element generation is not implemented"),
];

/// Read by the solver's own procedure; skipped for `Exec Solver = Never`.
const SOLVER: &[Rule] = &[
    present("Mesh", "solver-specific meshes are not implemented"),
    present("Slave Solvers", "slave solvers are not implemented"),
    flag("Monolithic Slave", "monolithic slave solvers are not implemented"),
    flag("Apply Integral BCs", "integral constraints are not implemented"),
    flag("Apply Source Control", "source control (Control Target Value) is not implemented"),
    flag("Apply Conforming BCs", "conforming BCs are not implemented"),
    flag("Apply Mortar BCs", "mortar BCs are not implemented"),
];

const BODY_FORCE: &[Rule] = &[flag("Integral BC", "integral constraints are not implemented")];

const BOUNDARY: &[Rule] = &[
    // MapBodiesAndBCs fallbacks `Model::bc_of_tag` does not apply.
    flag("Default Target", "default BC targets are not implemented"),
    flag("Default BC", "default BC targets are not implemented"),
    flag("Default Internal BC", "default BC targets are not implemented"),
    flag("Default External BC", "default BC targets are not implemented"),
    present("Intersection BC", "intersection BCs are not implemented"),
    present("Intersection Body", "intersection BCs are not implemented"),
    flag("Discontinuous Boundary", "discontinuous boundaries (mesh splitting) are not implemented"),
    present("Periodic BC", "periodic BCs are not implemented"),
    present("Mortar BC", "mortar BCs are not implemented"),
    present("Conforming BC", "conforming BCs are not implemented"),
    flag("Integral BC", "integral constraints are not implemented"),
];

fn scan(section: &str, i: usize, params: &ParamMap, rules: &[Rule]) -> Result<(), SolverError> {
    for r in rules {
        let hit = match r.test {
            Test::Present => params.get(r.key).is_some(),
            Test::True => params.logical(r.key) == Some(true),
        };
        if hit {
            return Err(SolverError::Unsupported(format!("{section} {i}: {}: {}", r.key, r.reason)));
        }
    }
    Ok(())
}

fn runs(solver: &ParamMap) -> bool {
    !solver.0.is_empty() && !solver.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
}

/// The name a `Variable` / `Exported Variable` value declares.
fn declared(spec: &str) -> Option<String> {
    let spec = parse_variable_spec(&spec.replace('"', " "), None).ok()?;
    Some(spec.name).filter(|n| !n.is_empty())
}

/// Solver variables and their components, and the exported variables no
/// solver section declares as its own.
fn variables(model: &Model) -> (Vec<String>, Vec<(usize, String)>) {
    let mut primary = Vec::new();
    let mut exported = Vec::new();
    for (i, s) in model.solvers.iter().enumerate() {
        if let Some(v) = s.string("variable") {
            if let Ok(spec) = parse_variable_spec(&v.replace('"', " "), s.integer("variable dofs")) {
                primary.push(spec.name);
                primary.extend(spec.components);
            }
        }
        for (key, value) in &s.0 {
            let numbered = key.strip_prefix("exported variable ").is_some_and(|n| n.bytes().all(|b| b.is_ascii_digit()));
            if let (true, Value::Str(v)) = (numbered, value) {
                exported.extend(declared(v).map(|n| (i + 1, n)));
            }
        }
    }
    exported.retain(|(_, n)| !primary.contains(n));
    (primary, exported)
}

/// Refuses the first unimplemented feature `model` names.
pub fn check(model: &Model) -> Result<(), SolverError> {
    scan("Simulation", 1, &model.simulation, SIMULATION)?;
    for (i, s) in model.solvers.iter().enumerate().filter(|(_, s)| runs(s)) {
        scan("Solver", i + 1, s, SOLVER)?;
    }
    for (i, b) in model.body_forces.iter().enumerate() {
        scan("Body Force", i + 1, b, BODY_FORCE)?;
    }
    for (i, b) in model.boundary_conditions.iter().enumerate() {
        scan("Boundary Condition", i + 1, b, BOUNDARY)?;
    }

    let (primary, exported) = variables(model);
    // SetNodalLoads reads `<Variable> Load` from BCs; CheckPassiveElement
    // reads `<Variable> Passive` from body forces.
    for var in &primary {
        let load = format!("{var} load");
        if let Some(i) = model.boundary_conditions.iter().position(|b| b.0.contains_key(&load)) {
            return Err(SolverError::Unsupported(format!(
                "Boundary Condition {}: {load}: nodal loads (<Variable> Load) are not implemented",
                i + 1
            )));
        }
        let passive = format!("{var} passive");
        if let Some(i) = model.body_forces.iter().position(|b| b.0.contains_key(&passive)) {
            return Err(SolverError::Unsupported(format!(
                "Body Force {}: {passive}: passive elements (<Variable> Passive) are not implemented",
                i + 1
            )));
        }
    }

    // A dependency on an exported variable reads a variable that
    // `Simulation` never creates.
    let sections: [(&str, &[ParamMap]); 4] = [
        ("Material", &model.materials),
        ("Body Force", &model.body_forces),
        ("Initial Condition", &model.initial_conditions),
        ("Boundary Condition", &model.boundary_conditions),
    ];
    for (section, maps) in sections {
        for (i, p) in maps.iter().enumerate() {
            for (key, value) in &p.0 {
                let Value::Dependent(d) = value else { continue };
                let args = d.variables.iter().flat_map(|v| v.split(',')).map(normalize_key);
                for arg in args {
                    if let Some((solver, name)) = exported.iter().find(|(_, n)| *n == arg) {
                        return Err(SolverError::Unsupported(format!(
                            "{section} {}: {key} depends on Exported Variable {name} of Solver {solver}: \
                             exported variables are not implemented",
                            i + 1
                        )));
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{Dependent, DependentKind};

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

    fn heat(extra: &[(&str, Value)]) -> ParamMap {
        let mut p = pm(&[("procedure", st("\"HeatSolve\" \"HeatSolver\"")), ("variable", st("Temperature"))]);
        for (k, v) in extra {
            p.insert(k, v.clone());
        }
        p
    }

    fn refusal(model: &Model) -> String {
        match check(model) {
            Err(SolverError::Unsupported(m)) => m,
            other => panic!("{other:?}"),
        }
    }

    fn equals(var: &str) -> Value {
        Value::Dependent(Dependent {
            variables: vec![var.into()],
            kind: DependentKind::Table { x: vec![1.0], y: vec![1.0], cubic: false },
            rows: 1,
            cols: 1,
        })
    }

    #[test]
    fn plain_heat_case_passes() {
        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("temperature", Value::Real { rows: 1, cols: 1, data: vec![0.0] })])],
            ..Model::default()
        };
        assert!(check(&model).is_ok());
    }

    #[test]
    fn keyword_rules() {
        let model = Model {
            simulation: pm(&[("increase element order", Value::Logical(true))]),
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Simulation 1: Increase Element Order: quadratic element generation is not implemented"
        );

        let model = Model { solvers: vec![heat(&[("slave solvers", Value::Integer(vec![1]))])], ..Model::default() };
        assert_eq!(refusal(&model), "Solver 1: Slave Solvers: slave solvers are not implemented");

        let model = Model {
            solvers: vec![heat(&[]), heat(&[("exec solver", st("never")), ("mesh", st("local"))])],
            boundary_conditions: vec![pm(&[]), pm(&[("periodic bc", Value::Integer(vec![2]))])],
            ..Model::default()
        };
        assert_eq!(refusal(&model), "Boundary Condition 2: Periodic BC: periodic BCs are not implemented");

        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("default internal bc", Value::Logical(false))])],
            ..Model::default()
        };
        assert!(check(&model).is_ok());
    }

    #[test]
    fn variable_named_keywords() {
        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("temperature load", equals("tempb loads"))])],
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Boundary Condition 1: temperature load: nodal loads (<Variable> Load) are not implemented"
        );

        let model = Model {
            solvers: vec![heat(&[])],
            body_forces: vec![pm(&[]), pm(&[("temperature passive", Value::Logical(true))])],
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Body Force 2: temperature passive: passive elements (<Variable> Passive) are not implemented"
        );
    }

    #[test]
    fn dependency_on_exported_variable() {
        let model = Model {
            solvers: vec![heat(&[("exported variable 1", st("-ip \"ip x\"")), ("exported variable 2", st("temperature"))])],
            boundary_conditions: vec![
                pm(&[("heat flux", equals("temperature"))]),
                pm(&[("heat flux", equals("IP X"))]),
            ],
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Boundary Condition 2: heat flux depends on Exported Variable ip x of Solver 1: \
             exported variables are not implemented"
        );
    }
}
